import os
import numpy as np
import xarray as xr
import matplotlib 
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from obspy.core.event import read_events
import yaml
import pandas as pd
import xarray as xr
import obspy 
from mayavi import mlab

class element_output:
    def __init__(self, path:str, grid_format:list) -> None:
        """Initializes the element object for the given path

        Args:
            path (str): a string that contains the path to the 
            directory just above "input" and "output"
            grid_format (list): see the AxiSEM3D input templates
            for the element output. The options are:
            [2], [0,2,4],['all']
        """        
        self.path = path
        inparam_output_path = path + '/input/inparam.output.yaml'
        inparam_source_path = path + '/input/inparam.source.yaml'
        # Get element group names and channels from the input file
        with open(inparam_output_path, 'r') as file:
            output_yaml = yaml.load(file, Loader=yaml.FullLoader)
            for element_group in output_yaml['list_of_element_groups']:
                # this only works when only one element group has been used
                self.element_group_name = list(element_group.keys())[0]
                self.channel_type  = element_group[self.element_group_name]['wavefields']['coordinate_frame']
        # get lat lon of the event
        with open(inparam_source_path, 'r') as file:
            source_yaml = yaml.load(file, Loader=yaml.FullLoader)
            source_name = list(source_yaml['list_of_sources'][0].keys())[0]
            # assume a single point source
            source = source_yaml['list_of_sources'][0][source_name]
            [self.source_lat, self.source_lon]=  source['location']['latitude_longitude']
            source_yaml = yaml.load(file, Loader=yaml.FullLoader)
            
        self.path_to_elements_output = path + '/output/elements/' + self.element_group_name
        self.grid_format = grid_format
        
        # Load the data, without the dict_data_wave because that one is too big. So we only
        # load the data that allows us to locate the wavefield data so we can later only load 
        # the file containing the specific data point we are interested in
        self.na_grid, self.data_time, self.list_element_na, self.list_element_coords, self.\
        dict_list_element, self.dict_data_wave = self._read_element_output(path, load_wave_data=False)
        self.rotation_matrix = self._compute_rotation_matrix()

    def _compute_rotation_matrix(self):
        """Computes the rotation matrix that aligns the z axis with the source axis

        Returns:
            np.ndarray: 3D rotation matrix
        """
        # get real earth coordinates of the sources
        colatitude = np.pi/2 - np.deg2rad(self.source_lat)
        longitude = np.deg2rad(self.source_lon)
        # rotation matrix into the source frame
        return np.asarray([[np.cos(colatitude) * np.cos(longitude), -np.sin(longitude), np.sin(colatitude) * np.cos(longitude)],
                           [np.cos(colatitude) * np.sin(longitude), np.cos(longitude), np.sin(colatitude) * np.sin(longitude)],
                           [-np.sin(colatitude), 0, np.cos(colatitude)]])
        
        
    def _geo_to_cyl(self, point: list) -> list:
        """Converts geographical coordinates to cylindrical 
        coordinates in the seismic frame

        Args:
            point (list): [radial position in m, latitude in deg, longitude in deg]

        Returns:
            list: [radial position in m, azimuth from source in rad]
        """        
        s_earth = point[0]
        theta_earth = np.deg2rad(point[1])
        phi_earth = np.deg2rad(point[2])

        # Transform spherical coordinates in earth frame of the point to 
        # cartesian coordinates in the earth frame
        # station location in real earth (spherical coords)
        [x_earth, y_earth, z_earth] = [s_earth * np.cos(theta_earth) * np.cos(phi_earth), 
                                       s_earth * np.cos(theta_earth) * np.sin(phi_earth), 
                                       s_earth * np.sin(theta_earth)]
        # rotate coordinates of the point to source frame
        [x, y, z] = np.matmul(self.rotation_matrix.transpose(), np.asarray([x_earth, y_earth, z_earth]))
        # transform to cylindrical coords in the source frame
        # z is already fine 
        s = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        
        return [s, z, phi]

        
    def load_data_at_point(self, point: list) -> np.ndarray:
        """Expands an inplane point into the longitudinal direction
        using the fourier expansion

        Args:
            point (list): [radial position m, latitude deg, longitude deg] in geographical coords

        Returns:
            np.ndarray: _description_
        """
        
        # Transform geographical to cylindrical coords in source frame
        s, z, phi = self._geo_to_cyl(point)
        # Interpolate the data inplane
        interpolated_data = self.inplane_interpolation(point)
        
        # I don't fully understand how this fourier reconstruction works ...
        # set complex type
        complex_type = np.complex32 if interpolated_data.dtype == np.complex64 else np.complex128

        # find max fourier order
        max_Fourier_order = len(interpolated_data[:,0,0]) // 2

        # initialize result with 0th order 
        result = interpolated_data[0].copy()
        # add higher orders
        for order in np.arange(1, max_Fourier_order + 1):
            coeff = np.zeros(result.shape, dtype=complex_type)
            # real part
            coeff.real = interpolated_data[order * 2 - 1]
            # complex part of Fourier coefficients
            if order * 2 < len(interpolated_data): # check for Nyquist
                coeff.imag += interpolated_data[order * 2]
            result += (2. * np.exp(1j * order * phi) * coeff).real

        return result
    

    def inplane_interpolation(self, point: list, plot: bool = False)-> np.ndarray:
        """Takes in a point in spherical coordinates in the real earth 
        frame and outputs the displacement data in time for all the 
        available channels in the form of a numpy array. 

        Args:
            point (list): _description_

        Returns:
            np.ndarray: _description_
        """        
        
        # Transform geographical to cylindrical coords in source frame
        s, z, phi = self._geo_to_cyl(point)
        
        if self.grid_format == [0,2,4]:
            # The of the element are positioned like this (GLL point)
            # The numbers inbetween the points are the sub-element numbers
            # ^z
            # | (2)-(5)-(8)
            # |  | 1 | 3 |
            # | (1)-(4)-(7)
            # |  | 0 | 2 |
            # | (0)-(3)-(6)
            #  ____________>s
            
            # find the difference vector between our chosen point
            # and the center GLL point of every element
            difference_vectors = self.list_element_coords[:,4,0:2] - [s, z]
            # find the index of the central GLL point that is closest to our point 
            element_index = np.argmin((difference_vectors*difference_vectors).sum(axis=1))
            # grab the information about the element whose center we just found
            element_na = self.list_element_na[element_index]
            # Find the sub-element where our point resides
            corner_points = np.asarray([self.list_element_coords[element_index, 0, 0:2], 
                                        self.list_element_coords[element_index, 2, 0:2], 
                                        self.list_element_coords[element_index, 6, 0:2], 
                                        self.list_element_coords[element_index, 8, 0:2]])
            difference_vectors_to_corners = corner_points - [s, z]
            sub_element_index = np.argmin((difference_vectors_to_corners**2).sum(axis=1))
            
            # Get the subelement indices 
            if sub_element_index == 0:
                sub_element_indices = [0,1,3,4]
            elif sub_element_index == 1:
                sub_element_indices = [1,2,4,5]
            elif sub_element_index == 2:
                sub_element_indices = [3,4,6,7]
            elif sub_element_index == 3:
                sub_element_indices = [4,5,7,8]
                
            # Get the locations of the subelement points
            sub_element_points = self.list_element_coords[element_index,sub_element_indices,:]
            
            # Interpolate the function at the point using the 4 subelement GLL points
            # We compute weights based on distances between our point and the various element points
            # (I should change this linear intewrpolation scheme with a better one)
            sub_elements_distances = np.sqrt(((sub_element_points - [s,z])**2).sum(axis=1))
            total_distances = sub_elements_distances.sum()
            weights = (total_distances - sub_elements_distances) / (3 * total_distances)
            # From the element_na we find the global tag of our element and search in the dict_list_element[nag]
            # at which index we have this tag. This index will be the index_in_data_wave
            index_in_data_wave = np.nonzero(self.dict_list_element[element_na[2]] == element_na[0])[0][0]
            data_wave = self.dict_data_wave[element_na[2]][index_in_data_wave]
            # finally we interpolate at our point
            interpolated_data = np.zeros(data_wave[:,0,:,:].shape)
            for i in range(4):
                interpolated_data += weights[i] * data_wave[:,sub_element_indices[i],:,:]
            
            
            if plot is True:
                mlab.figure(bgcolor=(0,0,0))
                i = -1
                """ for element in self.list_element_coords:
                    i += 1
                    element_na_i = self.list_element_na[i]
                    index_in_data_wave_i = np.nonzero(self.dict_list_element[element_na_i[2]] == element_na_i[0])[0][0]
                    data_wave_i = self.dict_data_wave[element_na_i[2]][index_in_data_wave_i]
                    if np.max(np.abs(np.asarray(data_wave_i))) > 1e-7:
                        for point in element:
                            mlab.points3d(point[0], 0, point[1], color=(1,1,1), opacity=1, scale_factor=5000) """
                for point in sub_element_points:
                    mlab.points3d(point[0], 0, point[1], color=(1,0,0), opacity=1, scale_factor=5000)
                mlab.points3d(sub_element_points[4][0], 0, sub_element_points[4][1], color=(0,0,1), opacity=1, scale_factor=5000)
                mlab.points3d(s,0,z,color=(0,1,0), opacity=1, scale_factor=5000)
                mlab.show()
                
                
        elif self.grid_format == [2]:
            # The of the element are positioned like this (GLL point)
            # Each element has just one point (the central point)
            # ^z
            # | 
            # |  
            # |     (0)
            # |  
            # | 
            #  ____________>s
            
            # find the difference vector between our chosen point
            # and the center GLL point of every element
            inter_elements_distances = np.asarray([0,0,0,0])
            difference_vectors = self.list_element_coords[:,0,0:2] - [s, z]
            distances = (difference_vectors*difference_vectors).sum(axis=1)
            # find the index of the central GLL point that is closest to our point 
            element_index = np.argmin(distances)
            inter_elements_distances[0] = distances[element_index]
            distances[element_index] += 1e20
            # find second closest element 
            second_element_index = np.argmin(distances)
            inter_elements_distances[1] = distances[second_element_index]
            distances[second_element_index] += 1e20
            # find third closest element 
            third_element_index = np.argmin(distances)
            inter_elements_distances[2] = distances[third_element_index]
            distances[third_element_index] += 1e20
            # find second closest element 
            fourth_element_index = np.argmin(distances)
            inter_elements_distances[3] = distances[fourth_element_index]
            # grab the information about the elements whose centers we just found
            element_na = self.list_element_na[[element_index, second_element_index, third_element_index, fourth_element_index]]            
             
            # Interpolate the function at the point using the 4 intra element GLL points
            # We compute weights based on distances between our point and the various element points
            total_distances = inter_elements_distances.sum()
            weights = (total_distances - inter_elements_distances) / (3 * total_distances)
            # From the element_na we find the global tag of our elements and search in the dict_list_element[nag]
            # at which index we have this tag. This index will be the index_in_data_wave
            data_wave = []
            for i in range(4):
                index_in_data_wave = np.nonzero(self.dict_list_element[element_na[i][2]] == element_na[i][0])[0][0]
                data_wave.append(self.dict_data_wave[element_na[i][2]][index_in_data_wave])
            # finally we interpolate at our point
            interpolated_data = np.zeros(data_wave[0][:,0,:,:].shape)
            for i in range(4):
                interpolated_data += weights[i] * data_wave[i][:,0,:,:]
            
        return interpolated_data


    def stream(self, path_to_station_file: str) -> obspy.Stream:
        """Takes in the path to a station file used for axisem3d
        and returns a stram with the wavefields computed at all stations

        Args:
            path_to_station_file (str): path to station.txt file

        Returns:
            obspy.stream: stream
        """        
        
        # Open station file
        stations = (pd.read_csv(path_to_station_file, 
                    delim_whitespace=True, 
                    header=0, 
                    names=["name","network","latitude","longitude","useless","depth"]))
        # initiate stream that will hold data 
        stream = obspy.Stream()
        for _, station in stations.iterrows():
            stalat = station['latitude']
            stalon = station['longitude']
            stadepth = station['depth']
            starad = 6371e3 - stadepth
            # get the data at this station (assuming RTZ components)
            wave_data = self.load_data_at_point([starad, stalat, stalon])
            
            delta = self.data_time[1] - self.data_time[0]
            npts = len(self.data_time)
            network = station['network']
            station_name = station['name']
            print(station_name)
            for chn_index, chn in enumerate(['LXR', 'LXT', 'LXZ']):
                # form the traces at the channel level
                trace = obspy.Trace(wave_data[chn_index])
                trace.stats.delta = delta
                trace.stats.ntps = npts
                trace.stats.network = network
                trace.stats.station = station_name
                trace.stats.location = ''
                trace.stats.channel = chn
                trace.stats.starttime = obspy.UTCDateTime("1970-01-01T00:00:00.0Z") + self.data_time[0]
                stream.append(trace)

        return stream

    def _read_element_output(self, load_wave_data=False):
        """Reads a folder that contains the element output files form
        Axisem3D and outputs a more readable format of the data as numpy 
        arrays.

        Args:
            load_wave_data (bool, optional): _description_. Defaults to True.

        Returns:
            na_grid (numpy array): a 1D array that contains all "Nr"s used in the 
                                    fourier expansions in the D domain. 
            data_time (np array): global time steps of the simulation
            list_element_na (np array): For each element it gives a 1D array that
                                        contains:
                                        1. element tag in the mesh
                                        2. actual "Nr"
                                        3. stored "Nr" (in case you didn't want to store
                                        all the Nrs)
                                        4. element index in the data (local)
                                        5. element index in the data (global)
            list_element_coords (np array): For each element, for each grid point,
                                            gives the coordinates in the D domain as
                                            (s, z)
            dict_list_element (dict): Lists of element tags? arranged by Nr in the dict
            dict_data_wave (dict): For each number of Nr, for each element, for each Nr layer,
                                    for each gll point, for each channel, the time wave data
        """        
        ################ open files ################
        # filenames
        nc_fnames = [f for f in os.listdir(self.path_to_elements_output) if 'axisem3d_synthetics.nc' in f]

        # open files
        nc_files = []
        for nc_fname in nc_fnames:
            nc_files.append(xr.open_dataset(self.path_to_elements_output + '/' + nc_fname))
        
        ################ variables that are the same in the datasets ################
        # read Na grid (all azimuthal dimensions)
        na_grid = nc_files[0].data_vars['list_na_grid'].values.astype(int)

        # read time
        data_time = nc_files[0].data_vars['data_time'].values
        
        
        ################ variables to be concatenated over the datasets ################
        # define empty lists of xarray.DataArray objects
        xda_list_element_na = []
        xda_list_element_coords = []
        dict_xda_list_element = {}
        dict_xda_data_wave = {}
        for nag in na_grid:
            dict_xda_list_element[nag] = []
            dict_xda_data_wave[nag] = []

        # loop over nc files
        for nc_file in nc_files:
            # append DataArrays
            xda_list_element_na.append(nc_file.data_vars['list_element_na'])
            xda_list_element_coords.append(nc_file.data_vars['list_element_coords'])
            for nag in na_grid:
                dict_xda_list_element[nag].append(nc_file.data_vars['list_element__NaG=%d' % nag])
                dict_xda_data_wave[nag].append(nc_file.data_vars['data_wave__NaG=%d' % nag])
                
        # concat xarray.DataArray
        xda_list_element_na = xr.concat(xda_list_element_na, dim='dim_element')
        xda_list_element_coords = xr.concat(xda_list_element_coords, dim='dim_element')
        for nag in na_grid:
            dict_xda_list_element[nag] = xr.concat(dict_xda_list_element[nag], dim='dim_element__NaG=%d' % nag)
            dict_xda_data_wave[nag] = xr.concat(dict_xda_data_wave[nag], dim='dim_element__NaG=%d' % nag)
            
        # read data to numpy.ndarray
        list_element_na = xda_list_element_na.values.astype(int)
        list_element_coords = xda_list_element_coords.values
        dict_list_element = {}
        dict_data_wave = {} 
        for nag in na_grid:
            dict_list_element[nag] = dict_xda_list_element[nag].values.astype(int)
            if load_wave_data:
                dict_data_wave[nag] = dict_xda_data_wave[nag].values

        ############### return ################
        if load_wave_data:
            return na_grid, data_time, list_element_na, list_element_coords, dict_list_element, dict_data_wave
        else:
            return na_grid, data_time, list_element_na, list_element_coords, dict_list_element, dict_xda_data_wave
        
