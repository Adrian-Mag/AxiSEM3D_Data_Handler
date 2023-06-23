from .axisem3d_output import AxiSEM3DOutput
from .element_output import ElementOutput
from .station_output import StationOutput

path_to_simulation_output = '/disks/data/PhD/CMB/simu1D/FORWARD_CHECK'
path_to_element_output = '/disks/data/PhD/CMB/simu1D_element/FORWARD_DATA/output/elements/entire_earth'
point = [6371000, 0, 20]
path_to_station_output = '/disks/data/PhD/CMB/simu1D/FORWARD_CHECK/output/stations/Station_grid'
channels = ['U']

simulation = ElementOutput(path_to_element_output)
data = simulation.stream(point, channels=channels)
simulation.obspyfy('/disks/data/PhD/CMB/stations/STA_0_90_180_STATIONS.txt')

""" simulation = StationOutput(path_to_station_output)
data = simulation.parse_to_mseed()
simulation.obspyfy() """
print('A')