from .element_output import ElementOutput
from .station_output import StationOutput
from line_profiler import LineProfiler


path_to_element_output = "/disks/data/PhD/CMB/simu1D_element/BACKWARD_DATA/output/elements/entire_earth"

def profile_load_data_at_point():
    element_obj = ElementOutput(path_to_element_output=path_to_element_output)
    element_obj.load_data_at_point([6371000, 0, 0])
    # element_obj.animation([6371000, 0, 0], [6371000, 0, 30], R_min=3400000, resolution=50)    

lp = LineProfiler()
lp.add_function(ElementOutput.load_data_at_point)
lp.add_function(ElementOutput.inplane_interpolation)
lp.add_function(ElementOutput._read_element_data)
# Add other internal methods to profile

lp_wrapper = lp(profile_load_data_at_point)
lp_wrapper()

lp.print_stats()