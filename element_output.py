import matplotlib 
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import numpy as np
import xarray as xr
import matplotlib 
matplotlib.use('tkagg')
import yaml
import pandas as pd
import xarray as xr
import obspy 
from obspy.core.inventory import Inventory, Network, Station, Channel
from tqdm import tqdm
import concurrent.futures
import time

from .axisem3d_output import AxiSEM3DOutput
from AxiSEM3D_Kernels import sph2cart, cart2sph

class ElementOutput(AxiSEM3DOutput):
    def __init__(self, path_to_element_output:str) -> None:
        """Initializes the ElementOutput object for the given path to the element output directory.

        Args:
            path_to_element_output (str): Path to the element output directory.
            element_group (str, optional): Name of the element group. If None, the first element group found will be used.

        Attributes:
            path_to_elements_output (str): Path to the element output directory.
            na_grid (numpy.ndarray): NA grid information.
            data_time (numpy.ndarray): Data time information.
            list_element_na (numpy.ndarray): List of element NA values.
            list_element_coords (numpy.ndarray): List of element coordinates.
            dict_list_element (dict): Dictionary of list of elements.
            files (list): List of element output files.
            elements_index_limits (list): List containing element index limits.
            rotation_matrix (numpy.ndarray): Rotation matrix.
            coordinate_frame (str): Coordinate frame of the wavefields.
            channels (list): List of wavefield channels.
            detailed_channels (list): List of channels component by component
            grid_format (str): Grid format for the in-plane coordinates.
            source_lat (float): Latitude of the event located on the axis.
            source_lon (float): Longitude of the event located on the axis.
        """
        path_to_simulation = self._find_simulation_path(path_to_element_output)
        super().__init__(path_to_simulation)

        self.element_group_name = os.path.basename(path_to_element_output)
        # Get element group names, channels, and grid format from the input file
        with open(self.inparam_output, 'r') as file:
            output_yaml = yaml.load(file, Loader=yaml.FullLoader)
            element_group = output_yaml['list_of_element_groups'][0][self.element_group_name]
            self.coordinate_frame  = element_group['wavefields']['coordinate_frame']
            self.channels = element_group['wavefields']['channels']
            self.grid_format = element_group['inplane']['GLL_points_one_edge']

        # get lat lon of the event located on the axis
        with open(self.inparam_source, 'r') as file:
            source_yaml = yaml.load(file, Loader=yaml.FullLoader)
            source_name = list(source_yaml['list_of_sources'][0].keys())[0]
            # assume a single point source
            source = source_yaml['list_of_sources'][0][source_name]
            [self.source_lat, self.source_lon]=  source['location']['latitude_longitude']
            self.source_depth = float(source['location']['depth'])

        self.path_to_elements_output = path_to_element_output
        # Get metadata 
        self.na_grid, self.data_time, self.list_element_na, self.list_element_coords, self.\
        dict_list_element, self.files, self.elements_index_limits, self.detailed_channels = self._read_element_metadata()
        # Replace the numerical indicators of coordinates with letters based on the cordinate system
        self.detailed_channels = [element.replace('1', self.coordinate_frame[0]).\
                                  replace('2', self.coordinate_frame[1]).\
                                  replace('3', self.coordinate_frame[2]) \
                                  for element in self.detailed_channels]
        
        self.rotation_matrix = self._compute_rotation_matrix()


    def obspyfy(self, path_to_station_file: str):
        # Create obspyfy folder if not existent already
        obspyfy_path = self.path_to_elements_output + '/obspyfied'
        if not os.path.exists(obspyfy_path):
            os.mkdir(obspyfy_path) 
        cat = self.catalogue
        cat.write(obspyfy_path + '/cat.xml', format='QUAKEML')

        stations_file_name = os.path.basename(path_to_station_file).split('.')[0]
        inv = self.create_inventory(path_to_station_file)
        inv.write(obspyfy_path + '/' + stations_file_name + '_inv.xml', format="stationxml")

        stream = self.stream_STA(path_to_station_file)
        stream.write(obspyfy_path + '/' + self.element_group_name + '.mseed', format="MSEED") 


    def create_inventory(self, path_to_station_file: str):
        ##################
        # Create Inventory
        ##################

        networks = []
        station_names = []
        locations = []
        channels_list = []

        # Get path to where the new inventory will be saved, and coordinates

        # Create new empty inventory
        inv = Inventory(
            networks=[],
            source="Inventory from axisem STATIONS file")

        # Open station file
        stations = (pd.read_csv(path_to_station_file, 
                    delim_whitespace=True, 
                    header=0, 
                    names=["name","network","latitude","longitude","useless","depth"]))

        # Iterate over all stations in the stations file
        for _, station in stations.iterrows():
            # Create network if not already existent
            net_exists = False
            for network in inv:
                if network.code == station['network']:
                    net_exists = True
                    net = network
            if net_exists == False:
                net = Network(
                code=station['network'],
                stations=[])
                # add new network to inventory
                inv.networks.append(net)
            
            # Create station (should be unique!)
            sta = Station(
            code=station['name'],
            latitude=station['latitude'],
            longitude=station['longitude'],
            elevation=-station['depth'])
            net.stations.append(sta)
            
            # Create the channels
            for channel in self.detailed_channels:
                cha = Channel(
                code=channel,
                location_code="",
                latitude=station['latitude'],
                longitude=station['longitude'],
                elevation=-station['depth'],
                depth=station['depth'],
                azimuth=None,
                dip=None,
                sample_rate=None)
                sta.channels.append(cha)
            
            # Form the lists that will be used as inputs with read_netcdf
            # to get the stream of the wavefield data
            networks.append(station['network'])
            station_names.append(station['name'])
            locations.append('') # Axisem does not use locations
            channels_list.append(self.detailed_channels)

        return inv
    

    def _find_simulation_path(self, path: str):
        """Takes in the path to a station file used for axisem3d
        and returns a stream with the wavefields computed at all stations

        Args:
            path_to_station_file (str): path to station.txt file

        Returns:
            parent_directory
        """
        current_directory = os.path.abspath(path)
        while True:
            parent_directory = os.path.dirname(current_directory)
            if os.path.basename(current_directory) == 'output':
                return parent_directory
            elif current_directory == parent_directory:
                # Reached the root directory, "output" directory not found
                return None
            current_directory = parent_directory


    def stream_STA(self, path_to_station_file: str, 
                   channels: list = None,
                   time_limits: list = None) -> obspy.Stream:
        """Takes in the path to a station file used for axisem3d
        and returns a stream with the wavefields computed at all stations

        Args:
            path_to_station_file (str): path to station.txt file

        Returns:
            obspy.stream: stream
        """        
        # Get time slices from time limits
        time_slices = np.where((self.data_time >= time_limits[0]) & (self.data_time <= time_limits[1]))

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
            starad = self.Earth_Radius - stadepth
            # get the data at this station (assuming RTZ components)
            wave_data = self.load_data_at_point(point=[starad, stalat, stalon],
                                                coord_in_deg=True,
                                                channels=channels, 
                                                time_slices=time_slices)
            # COnstruct metadata
            delta = self.data_time[1] - self.data_time[0]
            npts = len(self.data_time)
            network = station['network']
            station_name = station['name']
            if channels is not None:
                selected_detailed_channels = [element for element in self.detailed_channels \
                                            if any(element.startswith(prefix) for prefix in channels)]
            else:
                selected_detailed_channels = self.detailed_channels
            for chn_index, chn in enumerate(selected_detailed_channels):
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


    def stream(self, point: list, coord_in_deg: bool = True, 
               channels: list = None,
               time_limits: list = None) -> obspy.Stream:
        """
        Generate a stream with the wavefields computed at all stations given the location.

        Args:
            point (list): The location of the station in meters and radians/degrees.
                        It should be a list with the following elements:
                        - radial position in meters (float)
                        - latitude in radians/degrees (float)
                        - longitude in radians/degrees (float)
            coord_in_deg (bool, optional): Specify whether the coordinates are in degrees.
                                        Defaults to True.
            channels (list, optional): List of wavefield channels to include.
                                    Defaults to None, which includes all channels.
            time_limits (list, optional): Time limits for the data.
                                        It should be a list with two elements:
                                        - start time in seconds (float)
                                        - end time in seconds (float)
                                        Defaults to None, which includes all times.
            fourier_order (int, optional): Fourier order. Defaults to None.

        Returns:
            obspy.Stream: A stream containing the wavefields computed at all stations.
        """      
        # Get time slices from time limits
        if time_limits is not None:
            time_slices = np.where((self.data_time >= time_limits[0]) & (self.data_time <= time_limits[1]))
        else:
            time_slices = None

        if coord_in_deg:
            point[1] = np.deg2rad(point[1])
            point[2] = np.deg2rad(point[2])

        # initiate stream that will hold data 
        stream = obspy.Stream()
        # get the data at this station (assuming RTZ components)
        wave_data = self.load_data_at_point(point=point,
                                            channels=channels, 
                                            time_slices=time_slices, 
                                            coord_in_deg=False)
        # Construct metadata 
        delta = self.data_time[1] - self.data_time[0]
        npts = len(self.data_time)
        network = str(np.random.randint(0, 100))
        station_name = str(np.random.randint(0, 100))

        if channels is not None:
                selected_detailed_channels = [element for element in self.detailed_channels \
                                            if any(element.startswith(prefix) for prefix in channels)]
        else:
            selected_detailed_channels = self.detailed_channels
        for chn_index, chn in enumerate(selected_detailed_channels):
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
    

    def load_data_at_point(self, point: list, coord_in_deg: bool = False, 
                        channels: list = None,
                        time_slices: list = None) -> np.ndarray:
        """
        Expand an in-plane point into the longitudinal direction using the Fourier expansion.

        Args:
            point (list): A list representing the point in geographical coordinates.
                        It should contain the following elements:
                        - radial position in meters (float)
                        - latitude in degrees/radians (float)
                        - longitude in degrees/radians (float)
            coord_in_deg (bool, optional): Specify whether the coordinates are in degrees.
                                        Defaults to True.
            channels (list, optional): List of channels to include.
                                    Defaults to None, which includes all channels.
            time_slices (list, optional): Points on the time axis where data should be loaded
                                        Defaults to None, which includes all times.
            fourier_order (int, optional): Maximum Fourier order.
                                        Defaults to None.

        Returns:
            np.ndarray: The result of the Fourier expansion, represented as a NumPy array.
        """
        # Transform geographical to cylindrical coords in source frame
        if coord_in_deg:
            point[1] = np.deg2rad(point[1])
            point[2] = np.deg2rad(point[2])

        _, _, phi = self._geo_to_cyl(point)

        # Get channel slices from channels
        if channels is not None:
            if(self._check_elements(channels, self.channels)):
                # Filter by channel chosen
                channel_slices = []
                for channel in channels:
                    channel_slices += [index for index, element in enumerate(self.detailed_channels) if element.startswith(channel)]
            else:
                raise Exception('Only the following channels exist: ' + ', '.join(self.channels))
        else:
            channel_slices = None

        # Interpolate the data in-plane
        interpolated_data = self._inplane_interpolation(point=point, channel_slices=channel_slices, 
                                                        time_slices=time_slices)

        # Set complex type
        complex_type = interpolated_data.dtype if np.iscomplexobj(interpolated_data) else np.complex128

        # Find max Fourier order
        max_fourier_order = len(interpolated_data[:, 0, 0]) // 2

        # Initialize result with 0th order
        result = interpolated_data[0].copy()

        # Add higher orders
        for order in range(1, max_fourier_order + 1):
            coeff = np.zeros(result.shape, dtype=complex_type)
            # Real part
            coeff.real = interpolated_data[order * 2 - 1]
            # Complex part of Fourier coefficients
            if order * 2 < len(interpolated_data):  # Check for Nyquist
                coeff.imag += interpolated_data[order * 2]
            result += (2.0 * np.exp(1j * order * phi) * coeff).real

        return result


    def _inplane_interpolation(self, point: list, channel_slices: list = None, 
                              time_slices: list = None)-> np.ndarray:
        """Takes in a point in spherical coordinates in the real earth frame
        and outputs the displacement data in time for all the available channels
        in the form of a NumPy array.

        Args:
            point (list): A list representing the point in spherical coordinates.
                        It should contain the following elements:
                        - radial position in meters (float)
                        - latitude in degrees (float)
                        - longitude in degrees (float)
            channels (list, optional): List of channels to include. Defaults to None, which includes all channels.
            time_limits (list, optional): Time limits for the data. It should be a list with two elements:
                                        - start time in seconds (float)
                                        - end time in seconds (float)
                                        Defaults to None, which includes all times.

        Returns:
            np.ndarray: The interpolated displacement data in time for all available channels,
                        represented as a NumPy array.
        """      
        
        # Transform geographical to cylindrical coords in source frame
        s, z, _ = self._geo_to_cyl(point)
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
            for i in range(len(self.elements_index_limits) - 1):
                if self.elements_index_limits[i] <= element_index < self.elements_index_limits[i+1]:
                    file_index = i
                    break
            # get the element points
            element_points = self.list_element_coords[element_index]
            radial_element_GLL_points = self._cart_to_polar(element_points[[0,1,2]][:,0], 
                                                            element_points[[0,1,2]][:,1])[0]
            theta_element_GLL_points = self._cart_to_polar(element_points[[0,3,6]][:,0], 
                                                           element_points[[0,3,6]][:,1])[1]

            # Now we get the data
            data_wave = self._read_element_data(element_na=element_na, file_index=file_index,
                                                channel_slices=channel_slices, time_slices=time_slices)
            # finally we interpolate at our point
            interpolated_data = np.zeros(data_wave[:,0,:,:].shape)
            # Now we interpolate using GLLelement_na, file_index, channels, time_limits
            for i in range(3):
                for j in range(3):
                    interpolated_data += self._lagrange(r, radial_element_GLL_points[j], radial_element_GLL_points) * \
                    self._lagrange(theta, theta_element_GLL_points[i], theta_element_GLL_points) * data_wave[:,3*i+j,:,:]

        return interpolated_data


    def _read_element_metadata(self):
        """Reads a folder that contains the element output files from Axisem3D
        and outputs the metadata needed to access any data point from the mesh.

        Returns:
            na_grid (numpy array): A 1D array that contains all "Nr"s used in the 
                                Fourier expansions in the D domain.
            data_time (np array): Global time steps of the simulation.
            list_element_na (np array): For each element, it gives a 1D array that
                                        contains:
                                        1. Element tag in the mesh
                                        2. Actual "Nr"
                                        3. Stored "Nr" (in case you didn't want to store
                                        all the Nrs)
                                        4. Element index in the data (local)
                                        5. Element index in the data (global)
            list_element_coords (np array): For each element, for each grid point,
                                            gives the coordinates in the D domain as
                                            (s, z).
            dict_list_element (dict): Lists of element tags arranged by Nr in the dict.
            nc_files (list): List of opened netCDF files containing the element output data.
            elements_index_limits (list): List of element index limits for each netCDF file.
            detailed_channels (list): List of detailed channel names.

        Note:
            This method assumes that the element output files are stored in the
            `path_to_elements_output` directory.
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
        detailed_channels = [str_byte.decode('utf-8') for str_byte in nc_files[0].list_channel.data]
        updated_array = []
        for element in detailed_channels:
            letter = element[0]
            digits = ''.join(sorted(element[1:]))
            updated_array.append(letter + digits)
        detailed_channels = updated_array
        elements_index_limits = [0]
        index_limit = 0
        ######dict_xda_data_wave = {}
        for nag in na_grid:
            dict_xda_list_element[nag] = []
        
        # loop over nc files
        for i, nc_file in enumerate(nc_files):
            # append DataArrays
            index_limit += nc_file.sizes['dim_element']
            elements_index_limits.append(index_limit)
            xda_list_element_na.append(nc_file.data_vars['list_element_na'])
            xda_list_element_coords.append(nc_file.data_vars['list_element_coords'])
            for nag in na_grid:
                dict_xda_list_element[nag].append(nc_file.data_vars['list_element__NaG=%d' % nag])

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
        return na_grid, data_time, list_element_na, \
               list_element_coords, dict_list_element, \
               nc_files, elements_index_limits, \
               detailed_channels
        

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
        """
        Convert geographical coordinates to cylindrical coordinates in the
        seismic frame.

        Args:
            point (list): [radial position in m, latitude in rad, longitude in rad]

        Returns:
            list: [radial position in m, vertical position in m, azimuth from
            source in rad]
        """
        radial_pos = point[0]
        latitude = point[1]
        longitude = point[2]

        cos_lat = np.cos(latitude)
        sin_lat = np.sin(latitude)
        cos_lon = np.cos(longitude)
        sin_lon = np.sin(longitude)

        # Transform geographical coordinates to Cartesian coordinates in the
        # Earth frame
        x_earth = radial_pos * cos_lat * cos_lon
        y_earth = radial_pos * cos_lat * sin_lon
        z_earth = radial_pos * sin_lat

        # Rotate coordinates from the Earth frame to the source frame
        rotated_coords = np.matmul(self.rotation_matrix.transpose(), np.asarray([x_earth, y_earth, z_earth]))
        x, y, z = rotated_coords

        # Convert to cylindrical coordinates in the source frame
        s = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)

        return [s, z, phi]



    def _cart_to_polar(self, s: float, z: float) -> list:
        """Transforms inplane cylindrical coords (cartesian)
        to polar coords

        Args:
            s (float): distance from cylindarical axis in m
            z (float): distance along cylindrical axis in m

        Returns:
            list: [radius, theta]
        """       
        r = np.sqrt(s**2 + z**2)
        theta = np.arctan(z/s)
        
        return [r, theta]


    def _lagrange(self, evaluation_point, evaluation_GLL_point, GLL_points):
        """ Lagrange function implementation
        """
        value = 1
        for point in GLL_points:
            if evaluation_GLL_point != point:
                value *= (evaluation_point - point) / (evaluation_GLL_point - point)
        return value


    def _read_element_data(self, element_na, file_index: int, 
                           channel_slices: list=None, time_slices: list=None):
        """Reads the element data from the specified file and returns the wave data.

        Args:
            element_na (tuple): Element information containing:
                - Element tag in the mesh
                - Actual "Nr"
                - Stored "Nr" (in case you didn't want to store all the Nrs)
                - Element index in the data (local)
                - Element index in the data (global)
            file_index (int): Index of the file containing the element data.
            channels (list, optional): List of channels to include. Defaults to None.
            time_limits (list, optional): List of time limits [t_min, t_max]. Defaults to None.

        Returns:
            np.ndarray: The wave data.

        Raises:
            Exception: If the specified channels or time limits are not available.

        Note:
            - If `channels` and `time_limits` are both None, the entire wave data is returned.
            - If only `time_limits` is provided, the wave data is filtered by the specified time range.
            - If only `channels` are provided, the wave data is filtered by the specified channels.
            - If both `channels` and `time_limits` are provided, the wave data is filtered by both.

        Note that the wave data is assumed to be stored in the `files` attribute, which is a list of opened netCDF files.
        """
        if channel_slices is None and time_slices is None:
            wave_data = self.files[file_index]['data_wave__NaG=%d' % element_na[2]][element_na[3]].values
        elif channel_slices is not None and time_slices is None:
            wave_data = self.files[file_index]['data_wave__NaG=%d' % element_na[2]][element_na[3]][:, :, channel_slices, :].values
        elif channel_slices is None and time_slices is not None:
            wave_data = self.files[file_index]['data_wave__NaG=%d' % element_na[2]][element_na[3]][:, :, :, time_slices].values
        elif channel_slices is not None and time_slices is not None:
            wave_data = self.files[file_index]['data_wave__NaG=%d' % element_na[2]][element_na[3]][:, :, channel_slices, time_slices].values
        
        return wave_data


    def _check_elements(self, list1, list2):
        """Checks if all elements in list1 can be found in list2.

        Args:
            list1 (list): The first list.
            list2 (list): The second list.

        Returns:
            bool: True if all elements in list1 are found in list2, False otherwise.
            list: List of elements from list1 that are not found in list2.
        """
        missing_elements = [element for element in list1 if element not in list2]
        if len(missing_elements) == 0:
            return True
        else:
            return False
    

    def animation(self, source_location: list, station_location: list, channels: list=['U'],
                          name: str='video', video_duration: int=20, frame_rate: int=10,
                          resolution: int=50, R_min: float=0, R_max: float=6371000,
                          lower_range: float=0.1, upper_range: float=0.9,
                          paralel_processing: bool=True, timeit: bool=False):
        """
        Generate an animation representing seismic data on a slice frame.

        Args:
            source_location (list): The coordinates [depth, lat, lon] of the seismic source in the Earth frame.
            station_location (list): The coordinates [depth, lat, lon] of the station location in the Earth frame.
            name (str, optional): The name of the output video file. Defaults to 'video'.
            video_duration (int, optional): The duration of the video in seconds. Defaults to 20.
            frame_rate (int, optional): The number of frames per second in the video. Defaults to 10.
            resolution (int, optional): The resolution of the slice mesh. Defaults to 100.
            R_min (float, optional): The minimum radius for data inclusion. Defaults to 0.
            R_max (float, optional): The maximum radius for data inclusion. Defaults to 6371000.
            lower_range (float, optional): The lower percentile range for the colorbar intensity. Defaults to 0.5.

        Returns:
            None
        """
        # Get time slices from frame rate and video_duration assuming that the
        # video will include the entire time axis 
        no_frames = frame_rate*video_duration
        time_slices = np.round(np.linspace(0, len(self.data_time) - 1, no_frames)).astype(int)

        if timeit is True:
            start_time = time.time()
        print('Loading data')
        if paralel_processing is True:
            inplane_field, point1, point2, \
            base1, base2, inplane_DIM1, \
            inplane_DIM2 = self.load_data_on_slice_parallel(source_location, station_location, 
                                                        R_max, R_min, resolution, channels, 
                                                        time_slices=time_slices, return_slice=True)
        else:
            inplane_field, point1, point2, \
            base1, base2, inplane_DIM1, \
            inplane_DIM2 = self.load_data_on_slice_serial(source_location, station_location, 
                                                        R_max, R_min, resolution, channels, 
                                                        time_slices=time_slices, return_slice=True)
        print('Create animation')
        # Create a figure and axis
        cbar_min, cbar_max = self._find_range(np.log10(np.abs(inplane_field)), lower_range, upper_range)

        fig, ax = plt.subplots()
        contour = ax.contourf(inplane_DIM1, inplane_DIM2, 
                              np.nan_to_num(np.log10(np.abs(inplane_field[:, :, 0, 0]))), 
                              levels=np.linspace(cbar_min, cbar_max, 100), cmap='RdBu_r', 
                              extend='both')

        def update(frame):
            ax.cla()
            contour = ax.contourf(inplane_DIM1, inplane_DIM2, 
                                  np.nan_to_num(np.log10(np.abs(inplane_field[:, :, 0, frame]))), 
                                  levels=np.linspace(cbar_min, cbar_max, 100), cmap='RdBu_r', extend='both')
            plt.scatter(np.dot(point1, base1), np.dot(point1, base2))
            plt.scatter(np.dot(point2, base1), np.dot(point2, base2))
            print(100 * frame / no_frames, '%')
            return contour


        cbar = plt.colorbar(contour)

        cbar_ticks = np.linspace(int(cbar_min), int(cbar_max), 5) # Example tick values
        cbar_ticklabels = [str(cbar_tick) for cbar_tick in cbar_ticks] # Example tick labels
        cbar.set_ticks(cbar_ticks)
        cbar.set_ticklabels(cbar_ticklabels)
        cbar.set_label('Intensity')


        frame_step = int(len(self.data_time) / ( video_duration * frame_rate ))
        ani = animation.FuncAnimation(fig, update, frames=video_duration * frame_rate, interval=1e3 / frame_rate)
        ani.save(self.path_to_elements_output + '/' + name + '_animation.mp4', writer='ffmpeg')
        if timeit is True:
            end_time = time.time()
            print(end_time - start_time)


    def _load_data_at_point_parallel_wrapper(self, point, channels, time_slices, indices):
        return (self.load_data_at_point(point=point, channels=channels, time_slices=time_slices, 
                                        coord_in_deg=False), indices)

    
    def load_data_on_slice_parallel(self, source_location, station_location,
                                R_max, R_min, resolution, channels, time_slices,
                                return_slice: bool=False):
        """
        Load data on a slice of points within a specified radius range and
        resolution using parallel processing and batching.

        Args:
            source_location (list): The source location [depth, latitude,
            longitude] in degrees. station_location (list): The station location
            [depth, latitude, longitude] in degrees. R_max (float): The maximum
            radius for the slice in Earth radii. R_min (float): The minimum
            radius for the slice in Earth radii. resolution (int): The
            resolution of the slice (number of points along each dimension).
            channels (list): The channels of data to load. return_slice (bool,
            optional): Whether to return additional slice information. 
                                        Defaults to False.

        Returns:
            numpy.ndarray or list: An ndarray containing the loaded data on the
            slice,
                                and optionally, additional slice information.

        """
        if return_slice is False:
            filtered_indices, filtered_slice_points = self._create_slice(source_location, station_location, 
                                                                        R_max, R_min, resolution)
        else:
            filtered_indices, filtered_slice_points, \
            point1, point2, base1, base2, \
            inplane_DIM1, inplane_DIM2 = self._create_slice(source_location, station_location, 
                                                            R_max, R_min, resolution, return_slice=True)
        
        inplane_field = np.zeros((resolution, resolution, 3, len(time_slices)))

    
        pbar = tqdm(total=len(filtered_indices))

        with concurrent.futures.ProcessPoolExecutor() as executor:
            batch_size = 1000  # Adjust the batch size as needed
            num_batches = len(filtered_indices) // batch_size
            remaining_tasks = len(filtered_indices) % batch_size
            
            for i in range(num_batches):
                batch_indices = filtered_indices[i*batch_size:(i+1)*batch_size]
                batch_points = filtered_slice_points[i*batch_size:(i+1)*batch_size]
                
                futures = [executor.submit(self._load_data_at_point_parallel_wrapper, point, channels, time_slices, indices) \
                        for point, indices in zip(batch_points, batch_indices)]
                
                for future in concurrent.futures.as_completed(futures):
                    [index1, index2] = future.result()[1]
                    inplane_field[int(index1), int(index2), :, :] = future.result()[0]
                    
                
                    pbar.update(1)
            
            # Process the remaining tasks (if any)
            if remaining_tasks > 0:
                batch_indices = filtered_indices[-remaining_tasks:]
                batch_points = filtered_slice_points[-remaining_tasks:]
                
                futures = [executor.submit(self._load_data_at_point_parallel_wrapper, point, channels, time_slices, indices) \
                        for point, indices in zip(batch_points, batch_indices)]
                
                for future in concurrent.futures.as_completed(futures):
                    [index1, index2] = future.result()[1]
                    inplane_field[int(index1), int(index2), :, :] = future.result()[0]
                    
                
                    pbar.update(1)

    
        pbar.close()

        if return_slice is False:
            return inplane_field
        else:
            return [inplane_field, point1, point2, 
                    base1, base2, 
                    inplane_DIM1, inplane_DIM2]


    def load_data_on_slice_serial(self, source_location: list, station_location: list, 
                                  R_max: float, R_min: float, resolution: int, 
                                  channels: list, time_slices: list, return_slice: bool=False):
        """
        Load data on a slice of points within a specified radius range and resolution.
        Not using multi-processing!

        Args:
            source_location (list): The source location [depth, latitude, longitude] in degrees.
            station_location (list): The station location [depth, latitude, longitude] in degrees.
            R_max (float): The maximum radius for the slice in Earth radii.
            R_min (float): The minimum radius for the slice in Earth radii.
            resolution (int): The resolution of the slice (number of points along each dimension).
            channels (list): The channels of data to load.
            return_slice (bool, optional): Whether to return additional slice information. 
                                        Defaults to False.

        Returns:
            numpy.ndarray or list: An ndarray containing the loaded data on the slice,
                                and optionally, additional slice information.

        """
        if return_slice is False:
            filtered_indices, filtered_slice_points = self._create_slice(source_location, station_location, 
                                                                        R_max, R_min, resolution)
        else:
            filtered_indices, filtered_slice_points, \
            point1, point2, base1, base2, \
            inplane_DIM1, inplane_DIM2 = self._create_slice(source_location, station_location, 
                                                            R_max, R_min, resolution, return_slice=True)
        inplane_field = np.zeros((resolution, resolution, 3, len(time_slices)))
        pbar = tqdm(total=len(filtered_indices))
        for [index1, index2], point in zip(filtered_indices, filtered_slice_points):
            inplane_field[index1, index2, :, :] = self.load_data_at_point(point=point, channels=channels, coord_in_deg=False, time_slices=time_slices)
            pbar.update(1)
        pbar.close()

        if return_slice is False:
            return inplane_field
        else:
            return [inplane_field, point1, point2, 
                    base1, base2, 
                    inplane_DIM1, inplane_DIM2]


    def _create_slice(self, source_location: list, station_location: list, 
                      R_max: float, R_min: float, resolution: int,
                      return_slice: bool=False):
        """
        Create a mesh for a slice of Earth within a specified radius range and resolution.

        Args:
            source_location (list): The source location [depth, latitude, longitude] in degrees.
            station_location (list): The station location [depth, latitude, longitude] in degrees.
            R_max (float): The maximum radius for the slice in Earth radii.
            R_min (float): The minimum radius for the slice in Earth radii.
            resolution (int): The resolution of the slice (number of points along each dimension).
            return_slice (bool, optional): Whether to return additional slice information. 
                                        Defaults to False.

        Returns:
            list: A list containing filtered indices and slice points, and optionally, additional slice information.

        """
            # Transform from depth lat lon in deg to rad lat lon in rad
        source_location[0] = self.Earth_Radius - source_location[0]
        source_location[1] = np.deg2rad(source_location[1])
        source_location[2] = np.deg2rad(source_location[2])
        station_location[0] = self.Earth_Radius - station_location[0]
        station_location[1] = np.deg2rad(station_location[1])
        station_location[2] = np.deg2rad(station_location[2])
        # Form vectors for the two points (Earth frame)
        point1 = sph2cart(source_location[0], source_location[1], source_location[2])
        point2 = sph2cart(station_location[0], station_location[1], station_location[2])

        # Do Gram-Schmidt orthogonalization to form slice basis (Earth frame)
        base1 = point1 / np.linalg.norm(point1)
        base2 = point2 - np.dot(point2, base1) * base1
        base2 /= np.linalg.norm(base2)

        # Generate index mesh
        indices_dim1 = np.arange(resolution)
        indices_dim2 = np.arange(resolution)
        # Generate in-plane mesh
        inplane_dim1 = np.linspace(-R_max, R_max, resolution)
        inplane_dim2 = np.linspace(-R_max, R_max, resolution)
        inplane_DIM1, inplane_DIM2 = np.meshgrid(inplane_dim1, inplane_dim2, indexing='ij')
        radii = np.sqrt(inplane_DIM1*inplane_DIM1 + inplane_DIM2*inplane_DIM2)

        filtered_indices = []
        filtered_slice_points = []
        # Generate slice mesh points
        for index1 in indices_dim1:
            for index2 in indices_dim2:
                if radii[index1, index2] < R_max and radii[index1, index2] > R_min:
                    [x, y, z] = inplane_dim1[index1] * base1 + inplane_dim2[index2] * base2  # Slice frame -> Earth frame
                    filtered_slice_points.append(cart2sph(x, y, z))
                    filtered_indices.append([index1, index2])
        if return_slice is False:
            return [filtered_indices, filtered_slice_points]
        else:
            return [filtered_indices, filtered_slice_points,
                    point1, point2, base1, base2, 
                    inplane_DIM1, inplane_DIM2]


    def _find_range(self, arr, percentage_min, percentage_max):
        """
        Find the smallest value in the array based on the given percentage.

        Args:
            arr (ndarray): The input array.
            percentage (float): The percentage of values to consider.

        Returns:
            smallest_value (float or None): The smallest value based on the given percentage,
                                        or None if the array is empty or contains no finite values.
        """
        # Flatten the array to a 1D array
        flattened = arr[np.isfinite(arr)].flatten()
        
        if len(flattened) == 0:
            return None

        # Sort the flattened array in ascending order
        sorted_arr = np.sort(flattened)
        
        # Compute the index that corresponds to percentage of the values
        percentile_index_min = int(len(sorted_arr) * percentage_min)        
        percentile_index_max= int(len(sorted_arr) * percentage_max)
        
        # Get the value at the computed index
        smallest_value = sorted_arr[percentile_index_min]
        biggest_value = sorted_arr[percentile_index_max]
        
        return [smallest_value, biggest_value] 

