
'''
Section 1: Import Libraries
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

''''
Section 2: Make the grid
'''

#First the resolution
x_grid = 1e3 #x-direction
y_grid = 1e3 #y-direction
#Then the actual grid sizing
x_size = 20e3 #length in x-direction
y_size = 10e3 #lenth in y-direction

#Discretize in the x and y directions
x = np.arange(0,x_grid + x_size, x_grid)
y = np.arange(0, y_grid + y_size, y_grid)

print(x)
print(y)

#Well Info - km Easting, km Northing, m^2/s^2 Hydropotential
well_info = [
    [0, 0, 967],
    [0, 3, 985],
    [0, 6, 1000],
    [0, 8, 990],
    [20, 1, 965],
    [20, 4, 995],
    [20, 9, 970],
]

#Add values to the array
western_values = [np.array([0,3,6,8])*1e3, np.array([967,985,1000, 990])]
eastern_values = [np.array([1,4,9])*1e3, np.array([965, 995, 970])]

#Interpolation
potential_left = interp1d(western_values[0], western_values[1], kind = 'linear', fill_value="extrapolate")
potential_right = interp1d(eastern_values[0], eastern_values[1], kind = 'linear', fill_value = "extrapolate")

# Added by Yinjiao: evaluate the existing west/east interpolation
# functions at the grid points before interpolating north/south.

west_interp = potential_left(y)
east_interp = potential_right(y)

# Interpolate along the south (y = 0) and north (y = y_size) boundaries.
potential_bot = interp1d(
    [x[0], x[-1]],
    [west_interp[0], east_interp[0]],
    kind='linear'
)(x)

potential_top = interp1d(
    [x[0], x[-1]],
    [west_interp[-1], east_interp[-1]],
    kind='linear'
)(x)

print("Bottom boundary:", potential_bot)
print("Top boundary:", potential_top)
