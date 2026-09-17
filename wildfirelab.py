#Keeping track of what needs to happen in the lab
'''
What is actually happening?
    step 1: need to create a 3x3 forest with a
    buffer around the outside
    step 2: need to set the center of forest to be burning
    step 3: need to roll a dice to see where the fire will
    spread over a series of specific time steps, but it will
    not spread diagonally
    step 4: check how the model runs through visualizing
'''

#Section 1: Libraries

#The start of the script with a description of the purpose and
# a main place where libraries can be imported.

#!/usr/bin/env python3

'''
This file contains tools and scripts for completing Lab 1 for CLaSP410.
To reproduce the plots shown in the lab report, follow the steps below!

All instances of '# %%' are merely to divide the code into cells for each part.
The goal was to create ease when debugging individual sections of the cocde.
'''
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.gridspec as grid
import scipy as sp
from matplotlib.colors import ListedColormap #creating a custom color map

#Section 2: Defining Variables

#Define constants for the different states of the forest
burnt_bare = 1 #depicted as a tan color in the figures to show bare ground
alive_forest = 2 #depicted as a green color in the figures to show trees
on_fire = 3 #depicted as a red color in the figures to show flames

p_spread = 1.0 # Chance to spread to adjacent cells.
p_bare = 0.0 # Chance of cell to start as a bare patch.
p_ignite = 0.0 # Chance of cell to start on fire.

#Implement time as a varaible for the model
time_step = 1

#Section 3: Original Forest

#Create a 3x3 forest 

nx, ny = 5, 5 # Number of cells in X and Y direction, including buffer.

#creating the main forest array that the rest of the model will build from
ori_forest = np.zeros([ny, nx]) + 2

#setting the outside grid to be 1, so the fire will not spread into the ghost cells
for a in range(nx):
    ori_forest[0, a] = 1

for b in range(nx):
    ori_forest[b, 0] = 1

for c in range(ny):
    ori_forest[4, c] = 1

for d in range(ny):
    ori_forest[d, 4] = 1

print(ori_forest)

#Section 4: Begin the Burning

#Randomly generate forest with bare and forested cells
for i in range(1,nx-1): #Set range to avoid ghost nodes
    for j in range(1,ny-1):
        # Roll our "dice" to see if we get a bare spot:
        if np.random.rand() < p_bare:
            ori_forest[j, i] = burnt_bare # 1 is a bare spot randomly generated
        elif np.random.rand() < p_ignite:
            ori_forest[j,i] = on_fire

#Create a copy of the original forest to begin burning

ori_forest[2,2] = on_fire

burn_forest = np.copy(ori_forest)

print(burn_forest)

#Section 5: Spread the fire

'''
Create loop to begin rolling the dice to see where the randomly generated forest
will catch on fire and begin burning
'''

for k in range(0, time_step):
    for i in range(1,nx-1): #Set range to avoid ghost nodes
        for j in range(1,ny-1):
            print(f"Current Position: {j, i}")
            # Roll the dice to see if the fire will spread
            if ori_forest[j,i] == on_fire:
                if np.random.rand() < p_spread:
                    if ori_forest[j,i+1] == alive_forest:
                        burn_forest[j,i+1] = on_fire
                        print("Spreading!!")
                if np.random.rand() < p_spread:
                    if ori_forest[j+1,i] == alive_forest:
                        burn_forest[j+1, i] = on_fire
                        print("Spreading!!")
                if np.random.rand() < p_spread:
                    if ori_forest[j,i-1] == alive_forest:
                        burn_forest[j,i-1] = on_fire
                        print("Spreading!!")
                if np.random.rand() < p_spread:
                    if ori_forest[j-1,i] == alive_forest:
                        burn_forest[j-1,i] = on_fire
                        print("Spreading!!")

                burn_forest[j,i] = 1
                        
    ori_forest = np.copy(burn_forest)
   
    print(ori_forest)


#Section 6: Start visualizing the model

# Generate our custom segmented color map for this project.
# We can specify colors by names and then
# create a colormap that only uses those names. We have 3 fundamental
# states, so we want only 3 colors.
# Color info:
#https://matplotlib.org/stable/gallery/color/named_colors.html
forest_cmap = ListedColormap(['tan', 'darkgreen', 'firebrick'])

# Create figure and set of axes:
fig, ax = plt.subplots(1,1)

# Given our "forest" object, a 2D array that contains
# numbers 1, 2, or 3,
# Plot this using the "pcolor" method. Be sure to use our color map and
# set both *vmin* and *vmax*:
ax.pcolor(ori_forest, cmap=forest_cmap, vmin=1, vmax=3)


plt.show()

