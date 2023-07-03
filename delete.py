import numpy as np

# Example code snippet to generate inplane_DIM1 and inplane_DIM2
R_max = 1.0
resolution = 3

inplane_dim1 = np.linspace(-R_max, R_max, resolution)
inplane_dim2 = np.linspace(-R_max, R_max, resolution)
inplane_DIM1, inplane_DIM2 = np.meshgrid(inplane_dim1, inplane_dim2, indexing='xy')

# Stack inplane_DIM1 and inplane_DIM2 along the third axis
stacked_array = np.dstack((inplane_DIM1, inplane_DIM2))
print(inplane_DIM1)
print(inplane_DIM2)
# Print the shape of the stacked array
print(stacked_array)
stacked_array_reshaped = stacked_array.reshape(-1, 2)
print(stacked_array_reshaped)