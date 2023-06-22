from .axisem3d_output import AxiSEM3DOutput
from .element_output import ElementOutput


path_to_simulation_output = '/disks/data/PhD/CMB/simu1D/FORWARD_CHECK'
path_to_element_output = '/disks/data/PhD/CMB/simu3D_CMB_element/TEST/output/elements/mantle'
point = [6371000, 0, 20]
channels = ['S', 'U']
simulation = ElementOutput(path_to_element_output)
data = simulation.stream(point, channels=channels)
print('A')