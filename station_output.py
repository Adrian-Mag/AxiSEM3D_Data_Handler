""" 
Class for handling axisem station output in netcdf format.
"""
import numpy as np
import xarray as xr
import obspy 
from obspy import UTCDateTime
import pandas as pd
import netCDF4 as nc
import os 


class station_output:
    def __init__(self, path: str):
        """Class for handling axisem station output in netcdf format. 

        Args:
            path (str): Path to station output directory of simultaion
        """
        self.data = xr.open_mfdataset(path+"/axisem3d_synthetics.nc.*", 
                                      engine="netcdf4", 
                                      data_vars="different", 
                                      concat_dim="dim_station", 
                                      combine="nested")
        self.station_names_decoded = np.array([ls.decode("utf-8") for ls in self.data["list_station"].values])


    def station_to_stream(self, networks: list, station_names: list, locations: list, channels_list: list, data_type='RTZ') -> obspy.Stream:
        """Returns stream with data from the chosen stations and channels

        Args:
            networks (list): List of strings with networks
            station_names (list): List of strings with stations
            locations (list): List of strings with locations
            channels (str): String of channels
            data_type (str, optional): Type of data output (eg RTZ). Defaults to 'RTZ'.

        Returns:
            obspy.Stream: Obspy Stream of data 
        """
        channel_indices = []
        for channels in channels_list:
            channel_indices.append([data_type.find(x) for x in channels])
        # Initialize stream that will hold all the traces
        stream = obspy.Stream()
        # Compute sampling delta (it is the same for all traces) and no of points per trace
        delta = self.time[1] - self.time[0]
        npts = len(self.time)
        # Form full station name from components given (not including channel in the name)
        for i in range(len(networks)):
            network = networks[i]
            station_name = station_names[i]
            location = locations[i]
            
            for ichn, chn in enumerate(channel_indices[i]):
                if location == '':
                    full_name = str(network) + '.' + str(station_name)
                    index = self._find_station_index(full_name)
                    # extract data
                    displacement_data = self.data["data_wave"].values[index, chn, :]
                    # form obspy trace
                    trace = obspy.Trace(displacement_data)
                    trace.stats.delta = delta
                    trace.stats.ntps = npts
                    trace.stats.network = network
                    trace.stats.station = station_name
                    trace.stats.location = location
                    trace.stats.channel = 'LX' + channels[ichn]
                    trace.stats.starttime = self.starttime
                    stream.append(trace)
                else:
                    full_name = str(network) + '.' + str(station_name) + '_' + str(location)
                    index = self._find_station_index(full_name)
                    # extract data
                    displacement_data = self.data["data_wave"].values[index, chn, :]
                    # form obspy trace
                    trace = obspy.Trace(displacement_data)
                    trace.stats.delta = delta
                    trace.stats.ntps = npts
                    trace.stats.network = network
                    trace.stats.station = station_name
                    trace.stats.location = location
                    trace.stats.channel = channels[ichn]
                    trace.stats.starttime = self.starttime
                    stream.append(trace)

        return stream


    @property
    def time(self) -> np.ndarray:
        """Time axis of waveform

        Returns:
            np.ndarray: time
        """        
        return self.data['data_time'].values
    
    @property
    def starttime(self) -> UTCDateTime:
        default_event_time = UTCDateTime("1970-01-01T00:00:00.0Z")
        return default_event_time + self.time[0]      


    def _find_station_index(self, name: str) -> int:
        """Findx index of station in the data structure

        Args:
            name (str): Full name of station (including net, loc)

        Returns:
            int: Index of station
        """        
         # find indices of stations in the database 
        index = (name==self.station_names_decoded).argmax()
        if (index == 0) and (name != self.station_names_decoded[0]):
            raise Exception('Station could not be found')
        
        return index


def parse_station_output(path: str) -> obspy.Stream:
    """Function that takes the path to stations output, gathers the nc files 
    with the output from there and parses them into an obspy stream object.
    This could also be done with the station_output class, but less memory intensive
    because it tackles one station output file at a time. It also bypasses a possible
    problem with xarray concat function which seems to explode the RAM usage in certain
    cases.

    Args:
        path (str): path to the dir containing station output files

    Returns:
        obspy.Stream: just what it says
    """    
    # find all files with data in the given simulation dir
    pattern = 'axisem3d_synthetics.nc.rank' 
    matching_files = [f for f in os.listdir(path) if pattern in f]

    # read the rank list file
    rank_list = pd.read_csv(path + '/rank_station.info', sep=' ', header=None, names=['MPI_RANK', 'STATION_KEY', 'STATION_INDEX_IN_RANK'])

    # initiate stream object
    stream = obspy.Stream()
    
    # go file by file
    for file in matching_files:
        # get rank number of this file
        rank_no = file.split('.')[-1].strip('rank')
        
        # get wave data for this rank
        data = nc.Dataset(path + '/' + file)
        
        # get time for this rank (should be the same for all ranks
        # but I don't want to load a dataset outside the loop just to get the time)
        time = [float(t) for t in data['data_time'][:].data]
        delta = time[1] - time[0]
        npts = len(time)
        
        # get the list of statinos with this rank
        current_rank_list = rank_list[rank_list['MPI_RANK'].str.contains(rank_no)]
        # go station by station
        for index, trace_of_station in enumerate(data['data_wave']):
            station_key = current_rank_list.loc[current_rank_list['STATION_INDEX_IN_RANK'] == str(index), 'STATION_KEY'].iloc[0]
            network = station_key.split('.')[0]
            station_name = station_key.split('.')[-1]
            print(station_name)
            for chn_index, chn in enumerate(['LXR', 'LXT', 'LXZ']):
                # form the traces at the channel level
                trace = obspy.Trace(trace_of_station[chn_index].data)
                trace.stats.delta = delta
                trace.stats.ntps = npts
                trace.stats.network = network
                trace.stats.station = station_name
                trace.stats.location = ''
                trace.stats.channel = chn
                trace.stats.starttime = UTCDateTime("1970-01-01T00:00:00.0Z") + time[0]
                stream.append(trace)

    path_to_main_simulation_dir = path.partition("output")[0]
    name = path_to_main_simulation_dir.split('/')[-2]
    stream.write(path_to_main_simulation_dir + 'output/obspyfied/' + name + '.mseed', format="MSEED") 
