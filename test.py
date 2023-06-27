from .element_output import ElementOutput
from .station_output import StationOutput


path_to_element_output = "/disks/data/PhD/CMB/simu3D_CMB_element/test/output/elements/mantle"

point = [6371000, 0, 30]
#channels = ['E']
#fourier_order = 3
element_obj = ElementOutput(path_to_element_output)
#element_obj.stream(point=point, channels=channels).plot()
#element_obj.obspyfy('/disks/data/PhD/CMB/input3D_CMB_element/STA_10DEG_CROSS.txt')
element_obj.animation([6371000,0,0], [6371000,0,30], R_min=3480000, resolution=200, name='range_80', lower_range=0.7)


#path_to_station_output = '/disks/data/PhD/CMB/simu3D_CMB_element/test/output/stations/Station_grid'

#station_obj = StationOutput(path_to_station_output)
#print(station_obj.stream(networks='A', station_names='22'))
#station_obj.obspyfy()