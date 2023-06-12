""" 
Code that parses the "stations" or "elements" netcdf output of AxiSEM3D
into obspy compatible files (catalogueXML, stationXML, mseed)
"""
import yaml
from obspy import UTCDateTime
from obspy.core.event import Catalog, Event, Origin, FocalMechanism, MomentTensor, Tensor
from obspy.geodetics import FlinnEngdahl
import os
from obspy.core.inventory import Inventory, Network, Station, Channel
import pandas as pd
import glob
from .station_output import parse_station_output
from .element_output import element_output


def obspyfy(path, output_type, stations_paths = None, grid_format = [0,2,4]):
    name = path.split('/')[-1]

    ##################
    # Create catalogue
    ##################

    source_path = path + '/input/inparam.source.yaml'
    cat = Catalog()

    with open(source_path, 'r') as file:
                source_yaml = yaml.load(file, Loader=yaml.FullLoader)
                
                for source in source_yaml['list_of_sources']:
                    for items in source.items():

                        event = Event()
                        
                        origin = Origin()
                        
                        origin.time = UTCDateTime("1970-01-01T00:00:00.0Z") # set as the default in obspy
                        origin.latitude = items[1]['location']['latitude_longitude'][0]
                        origin.longitude = items[1]['location']['latitude_longitude'][1]
                        origin.depth = items[1]['location']['depth']
                        origin.depth_type = "operator assigned"
                        origin.evaluation_mode = "manual"
                        origin.evaluation_status = "preliminary"
                        origin.region = FlinnEngdahl().get_region(origin.longitude, origin.latitude)
                        
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

    if not os.path.exists(path + '/output/obspyfied'):
        os.mkdir(path + '/output/obspyfied')         
    cat.write(path + '/output/obspyfied/cat.xml', format='QUAKEML')

    if output_type == 'stations':
        ########################################### STATIONS OUTPUTS ###########################################

        ##################
        # Create Inventory
        ##################

        stations_paths = glob.glob(path + '/input/STA*.txt')
        inparam_output_path = path + '/input/inparam.output.yaml'
        networks = []
        station_names = []
        locations = []
        channels_list = []

        for stations_path in stations_paths:

            # Get path to where the new inventory will be saved, and coordinates
            # of the stations channels (RTZ, etc)
            with open(inparam_output_path, 'r') as file:
                        output_yaml = yaml.load(file, Loader=yaml.FullLoader)
                        for station_grid in output_yaml['list_of_station_groups']:
                            station_grid_name = list(station_grid.keys())[0]
                            channel_type  = station_grid[station_grid_name]['wavefields']['coordinate_frame']

            # Create new empty inventory
            inv = Inventory(
                networks=[],
                source="Inventory from axisem STATIONS file")

            # Open station file
            stations = (pd.read_csv(stations_path, 
                        delim_whitespace=True, 
                        header=0, 
                        names=["name","network","latitude","longitude","useless","depth"]))

            # Iterate over all stations in the stations file
            for index, station in stations.iterrows():
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
                for channel in channel_type:
                    cha = Channel(
                    code='LX' + channel,
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
                channels_list.append(channel_type)

            # Save the newly formed inventory
            inv_name = stations_path.split('/')[-1].split('.')[0]
            inv.write(path + '/output/obspyfied/' + inv_name + '_inv.xml', format="stationxml")

            ###################
            # Create mseed file
            ###################

            # Get path to netcdf files containing the "traces" 
            streams_path = path + '/output/stations/' + station_grid_name + '/'
            parse_station_output(streams_path)
    else:
        
        ########################################### ELEMENT OUTPUTS ###########################################
        """ 
        The elements output is handled using the elements_output class. For this, 
        the grid format needs to be known and passed. It also requires a list of 
        paths to station.txt files. The code will go through each station.txt file,
        get the location of each station and get the interpolated wavefield there, 
        putting it in a stram object, which is finally saved. 
        """
        # create element object
        obj = element_output(path, grid_format=grid_format)
        ##################
        # Create Inventory
        ##################

        stations_path = stations_paths
        networks = []
        station_names = []
        locations = []
        channels_list = []
        
        channel_type = obj.channel_type
        
        # Create new empty inventory
        inv = Inventory(
            networks=[],
            source="Inventory from axisem STATIONS file")

        # Open station file
        stations = (pd.read_csv(stations_path, 
                    delim_whitespace=True, 
                    header=0, 
                    names=["name","network","latitude","longitude","useless","depth"]))

        # Iterate over all stations in the stations file
        for index, station in stations.iterrows():
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
            for channel in channel_type:
                cha = Channel(
                code='LX' + channel,
                location_code="",
                latitude=station['latitude'],
                longitude=station['longitude'],
                elevation=-station['depth'],
                depth=station['depth'],
                azimuth=None,
                dip=None,
                sample_rate=None)
                sta.channels.append(cha)
            
            # Form the lists that will be used as inputs with element_output
            # to get the stream of the wavefield data
            networks.append(station['network'])
            station_names.append(station['name'])
            locations.append('') # Axisem does not use locations
            channels_list.append(channel_type)

            # Save the newly formed inventory
            inv_name = stations_path.split('/')[-1].split('.')[0]
            inv.write(path + '/output/obspyfied/' + inv_name + '_inv.xml', format="stationxml")

        ###################
        # Create mseed file
        ###################

        
        stream = obj.stream_STA(path_to_station_file=stations_path)
        name = obj.element_group_name
        stream.write(path + '/output/obspyfied/' + name + '.mseed', format="MSEED") 
