from .element_output import ElementOutput
from .station_output import StationOutput
from line_profiler import LineProfiler

import time 

def profiled_code():
    path_to_element_output = "/disks/data/PhD/CMB/simu1D_element/BACKWARD_DATA/output/elements/entire_earth"
    element_obj = ElementOutput(path_to_element_output=path_to_element_output)
    element_obj.animation([0, 0, 0], [0, 0, 30], R_min=3400000, resolution=50, 
                          timeit=True, paralel_processing=True)
    #point = [6371000, 0, 45]
    #element_obj.load_data_at_point(point, coord_in_deg=True, channels=['U'])

lp = LineProfiler()
lp.add_function(ElementOutput.animation)
lp.add_function(ElementOutput.load_data_on_slice_parallel)

lp.wrapper = lp(profiled_code)
lp.wrapper()

lp.print_stats()