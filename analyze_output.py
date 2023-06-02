import sys
sys.path.append('/home/adrian/PhD/AxiSEM3D/Output_Handlers')
from element_output import element_output

print(os.getpid())

# Get the forward and backward data 
path_to_backward = '/home/adrian/PhD/AxiSEM3D/CMB/simu1D_element/ADJOINT'

# Create element objects
backward_data = element_output(path_to_backward, [0,2,4])

# import inversion mesh
points = pd.read_csv('/home/adrian/PhD/AxiSEM3D/CMB/stations/3D_UNIFORM_STA.txt', sep=" ")