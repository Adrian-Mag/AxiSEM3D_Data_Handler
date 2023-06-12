import os
import numpy as np
import xarray as xr
import matplotlib 
matplotlib.use('tkagg')
import yaml
import pandas as pd
import xarray as xr
import obspy 
from obspy.core.event import Catalog, Event, Origin, FocalMechanism, MomentTensor, Tensor
from obspy import UTCDateTime
from obspy.geodetics import FlinnEngdahl


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
        # Basic properties of the element object
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
            
        self.path_to_elements_output = path + '/output/elements/' + self.element_group_name
        self.grid_format = grid_format
        self.na_grid, self.data_time, self.list_element_na, self.list_element_coords, self.\
        dict_list_element, self.files, self.dict_file_no_na_grid = self._read_element_metadata()
        self.rotation_matrix = self._compute_rotation_matrix()


    def create_catalogue(self):
        # Create a catalogue
        source_path = self.path + '/input/inparam.source.yaml'
        with open(source_path, 'r') as file:
                source_yaml = yaml.load(file, Loader=yaml.FullLoader)
                cat = Catalog()
                for source in source_yaml['list_of_sources']:
                    for items in source.items():
                        event = Event()
                        origin = Origin()
                        
                        origin.time = UTCDateTime("1970-01-01T00:00:00.0Z") # default in obspy
                        origin.latitude = items[1]['location']['latitude_longitude'][0]
                        origin.longitude = items[1]['location']['latitude_longitude'][1]
                        origin.depth = items[1]['location']['depth']
                        origin.depth_type = "operator assigned"
                        origin.evaluation_mode = "manual"
                        origin.evaluation_status = "preliminary"
                        origin.region = FlinnEngdahl().get_region(origin.longitude, 
                                                                  origin.latitude)
                        if items[1]['mechanism']['type'] == 'FORCE_VECTOR':
                            m_rr = items[1]['mechanism']['data'][0]
                            m_tt = items[1]['mechanism']['data'][1]
                            m_pp = items[1]['mechanism']['data'][2]
                            m_rt = 0
                            m_rp = 0
                            m_tp = 0
                        else: 
                            m_rr = items[1]['mechanism']['data'][0]
                            m_tt = items[1]['mechanism']['data'][1]
                            m_pp = items[1]['mechanism']['data'][2]
                            m_rt = items[1]['mechanism']['data'][3]
                            m_rp = items[1]['mechanism']['data'][4]
                            m_tp = items[1]['mechanism']['data'][5]
                        
                        focal_mechanisms = FocalMechanism()
                        tensor = Tensor()
                        moment_tensor = MomentTensor()
                        tensor.m_rr = m_rr
                        tensor.m_tt = m_tt
                        tensor.m_pp = m_pp
                        tensor.m_rt = m_rt
                        tensor.m_rp = m_rp
                        tensor.m_tp = m_tp
                        moment_tensor.tensor = tensor
                        focal_mechanisms.moment_tensor = moment_tensor
                                            
                        # make associations, put everything together
                        cat.append(event)
                        event.origins = [origin]
                        event.focal_mechanisms = [focal_mechanisms]
        return cat


    def stream_STA(self, path_to_station_file: str) -> obspy.Stream:
        """Takes in the path to a station file used for axisem3d
        and returns a stream with the wavefields computed at all stations

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


    def stream(self, stadepth: float, stalat: float, stalon: float) -> obspy.Stream:
        """Takes in the location of a station in meters and degrees
        and returns a stream with the wavefields computed at all stations

        Args:
            stadepth: station depth in m
            stalat: station latitude in deg
            stalon: station longitude in deg

        Returns:
            obspy.stream: stream
        """        
        
        # initiate stream that will hold data 
        stream = obspy.Stream()
        starad = 6371e3 - stadepth
        # get the data at this station (assuming RTZ components)
        wave_data = self.load_data_at_point([starad, stalat, stalon])
        
        # Construct metadata 
        delta = self.data_time[1] - self.data_time[0]
        npts = len(self.data_time)
        network = str(np.random.randint(0, 100))
        station_name = str(np.random.randint(0, 100))
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


    def inplane_interpolation(self, point: list)-> np.ndarray:
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
        # spherical coordinates will be used for the GLL interpolation
        [r, theta] = self._cart_to_polar(s,z)
        
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
            
            # get the element points
            element_points = self.list_element_coords[element_index]
            radial_element_GLL_points = self._cart_to_polar(element_points[[0,1,2]][:,0], 
                                                            element_points[[0,1,2]][:,1])[0]
            theta_element_GLL_points = self._cart_to_polar(element_points[[0,3,6]][:,0], 
                                                           element_points[[0,3,6]][:,1])[1]

            # Now we get the data
            data_wave = self._read_element_data(element_na)
            # finally we interpolate at our point
            interpolated_data = np.zeros(data_wave[:,0,:,:].shape)
            # Now we interpolate using GLL
            for i in range(3):
                for j in range(3):
                    interpolated_data += self._lagrange(r, radial_element_GLL_points[j], radial_element_GLL_points) * \
                    self._lagrange(theta, theta_element_GLL_points[i], theta_element_GLL_points) * data_wave[:,3*i+j,:,:]

            
        return interpolated_data
    
    
    def _read_element_metadata(self):
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
        # filenames (sorted correctly)
        nc_fnames = sorted([f for f in os.listdir(self.path_to_elements_output) if 'axisem3d_synthetics.nc' in f])
        
        # open files
        nc_files = []
        for nc_fname in nc_fnames:
            nc_files.append(xr.open_dataset(self.path_to_elements_output + '/' + nc_fname))
        
        ################ variables that are the same in the datasets ################
        # read Na grid (all azimuthal dimensions)
        na_grid = nc_files[0].data_vars['list_na_grid'].values.astype(int)

        # read time
        data_time = nc_files[0].data_vars['data_time'].values
        
        ################ variables to be concatenated over the datasets minud the data itself################
        # define empty lists of xarray.DataArray objects
        xda_list_element_na = []
        xda_list_element_coords = []
        dict_xda_list_element = {}
        dict_file_no_na_grid = []
        ######dict_xda_data_wave = {}
        for nag in na_grid:
            dict_xda_list_element[nag] = []
        
        # loop over nc files
        for i, nc_file in enumerate(nc_files):
            # append DataArrays
            xda_list_element_na.append(nc_file.data_vars['list_element_na'])
            xda_list_element_coords.append(nc_file.data_vars['list_element_coords'])
            for nag in na_grid:
                dict_xda_list_element[nag].append(nc_file.data_vars['list_element__NaG=%d' % nag])
            dict_file_no_na_grid.append((nc_file.list_element_na.data[0][-1], nc_file.list_element_na.data[-1][-1]))
                
        # concat xarray.DataArray
        xda_list_element_na = xr.concat(xda_list_element_na, dim='dim_element')
        xda_list_element_coords = xr.concat(xda_list_element_coords, dim='dim_element')
        for nag in na_grid:
            dict_xda_list_element[nag] = xr.concat(dict_xda_list_element[nag], dim='dim_element__NaG=%d' % nag)
            
        # read data to numpy.ndarray
        list_element_na = xda_list_element_na.values.astype(int)
        list_element_coords = xda_list_element_coords.values
        dict_list_element = {}
        for nag in na_grid:
            dict_list_element[nag] = dict_xda_list_element[nag].values.astype(int)

        ############### return ################
        # Here we return the files only because in this format they are not being loaded into RAM
        # Since these files are huge we prefer to load into ram only the file where the data that we 
        # want is located and then close the file. 
        return na_grid, data_time, list_element_na, list_element_coords, dict_list_element, nc_files, dict_file_no_na_grid
        
    
    def _compute_rotation_matrix(self):
        """Computes the rotation matrix that aligns the z axis with the source axis

        Returns:
            np.ndarray: 3D rotation matrix
        """
        # get real earth coordinates of the sources
        colatitude = np.pi/2 - np.deg2rad(self.source_lat)
        longitude = np.deg2rad(self.source_lon)

        # rotation matrix into the source frame (based on Tarje's PhD)
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


    def _cart_to_polar(self, s: float, z: float) -> list:
        """Transforms inplane cylindrical coords (cartesian)
        to polar coords

        Args:
            s (float): distance from cylindarical axis
            z (float): distance along cylindrical axis

        Returns:
            list: [radius, theta]
        """       
        r = np.sqrt(s**2 + z**2)
        theta = np.arctan(z/s)
        
        return [r, theta]


    def _lagrange(self, evaluation_point, evaluation_GLL_point, GLL_points):
        """ Lagrange functino implementation
        """
        value = 1
        for point in GLL_points:
            if evaluation_GLL_point != point:
                value *= (evaluation_point - point) / (evaluation_GLL_point - point)
        return value


    def _read_element_data(self, element_na):
        # First we find in which file the data must be searched
        data_index = element_na[4]
        file_index = 0
        FOUND_FILE = False
        file_index = 0
        while FOUND_FILE is False:
            if (data_index >= self.dict_file_no_na_grid[file_index][0] 
                and data_index <= self.dict_file_no_na_grid[file_index][1]):
                FOUND_FILE = True
            else:
                file_index += 1
        # and extract the data from that file
        wave_data = self.files[file_index]['data_wave__NaG=%d' % element_na[2]][element_na[3]]

        return wave_data.values
    