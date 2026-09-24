#!/usr/bin/env python3

# Zombie Virus Model

'''
This part of the lab will modify the original wildfire model
to explore how we can use a similar model for tracking how 
diseases, such as the zombie virus, spreads.
'''

#Section 1: Libraries

#The start of the script with a description of the purpose and
# a main place where libraries can be imported.

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
import matplotlib.colors as mcolors


#Section 2: Defining Variables

#Define constants for the different states of the virus spread
deceased = 0 #depicted as a black color in the figures to show those who perished
immune = 1 #depicted as a tan color in the figures to show those who begin immune
healthy = 2 #depicted as a green color in the figures to show healthy individuals
infected = 3 #depicted as a red color in the figures to show where the virus has spread

p_spread = 0.8 # Chance to spread to adjacent cells
p_immune = 0.5 # Chance of person to start as immune
p_zombie = 0.05 # Chance of infection to begin in a cell
p_fatal = 0.01 # Chance of infected to perish instead of recover

#Implement time as a varaible for the model
time_step = 2

#Create an empty list so that the values from each time step can be saved
outbreaks = []

#Section 3: Original Population

#Create a population

nx, ny = 20, 20 # Number of cells in X and Y direction 
# **needs to be one step bigger in each direction to account for ghost cells

#creating the main array that the rest of the model will build from
no_outbreak = np.zeros([ny, nx]) + 2

#Setting the outside grid to be 1, so the virus will not spread into the ghost cells

no_outbreak[[ny-1], :] = 1 # Selecting the bottom row of the grid
no_outbreak[:, [nx-1]] = 1 # Selecting the last column of the grid
no_outbreak[0, :] = 1 # Selecting the top row of the grid
no_outbreak[:, 0] = 1 # Selecting the first column of the grid

print(no_outbreak)


#Section 4: Begin the outbreak

#Randomly generate population with immune and healthy people
for i in range(1,nx-1): #Set range to avoid ghost nodes
    for j in range(1,ny-1):
        # Roll our "dice" to see if we get an immune person:
        if np.random.rand() < p_immune:
            no_outbreak[j, i] = immune # Immune is a randomly generated immune person
        # Roll the dice to see if anyone will become infected:
        elif np.random.rand() < p_zombie:
            no_outbreak[j,i] = infected # Infected is a randomly generated carrier of the virus

#Create a copy of the original population to begin spreading

infected_pop = np.copy(no_outbreak)

print(infected_pop)

#Section 5: Spread the disease

'''
Create loop to begin rolling the dice to see where the randomly generated population
will have an infected person appear and begin spreading the virus
'''

# Create a funciton to roll the dice for percentage of infected who perish or become immune

#Create a function to roll specifically to see if an infected person will live or die
#This is separate from the main for loop as it will be applied within the loop later
def chance_perish():
    for i in range(1,nx-1): #Set range to avoid ghost nodes
        for j in range(1,ny-1):
            #Roll the dice to see if the infected person will die
            if np.random.rand() < p_fatal:
                infected_pop[j, i] = deceased
                print("You Died :(")
            #If they don't die, they become immune
            elif infected_pop[j, i] == immune:
                print("You Survived :)")

#Function for rolling the dice over all possible neighbors
for k in range(time_step):
    for i in range(1,nx-1): #Set range to avoid ghost nodes
        for j in range(1,ny-1):
            print(f"Current Position: {j, i}")
            # Roll the dice to see if the disease will spread
            if no_outbreak[j,i] == infected:
                if np.random.rand() < p_spread:
                    if no_outbreak[j,i+1] == healthy:
                        infected_pop[j,i+1] = infected
                        print("Spreading!!")
                if np.random.rand() < p_spread:
                    if no_outbreak[j+1,i] == healthy:
                        infected_pop[j+1, i] = infected
                        print("Spreading!!")
                if np.random.rand() < p_spread:
                    if no_outbreak[j,i-1] == healthy:
                        infected_pop[j,i-1] = infected
                        print("Spreading!!")
                if np.random.rand() < p_spread:
                    if no_outbreak[j-1,i] == healthy:
                        infected_pop[j-1,i] = infected
                        print("Spreading!!")
                # Roll the dice to see if any of the infected individuals perish or become immune
                infected_pop[j, i] = chance_perish() 
                #Any of the infected population that does survive becomes immune
                infected_pop[j, i] = immune
    
    no_outbreak = np.copy(infected_pop)
    outbreaks.append(no_outbreak)
   
    print(no_outbreak)
    print(outbreaks)

    
#Section 6: Start visualizing the model

# Generate our custom segmented color map for this project.
# We can specify colors by names and then
# create a colormap that only uses those names. We have 4 fundamental
# states, so we want only 4 colors.
# Color info:
#https://matplotlib.org/stable/gallery/color/named_colors.html
forest_cmap = ListedColormap(['black', 'tan', 'darkgreen', 'firebrick'])

# Create figure and set of axes:
fig, ax = plt.subplots(1,1, figsize=(5,5))
plt.axhline(0, color='k', linestyle='--', linewidth=0.7)
plt.title(f"Zombie Outbreak Spread Model {time_step}")

# Given our object, a 2D array that contains
# numbers 1, 2, or 3,
# Plot this using the "pcolor" method. Be sure to use our color map and
# set both *vmin* and *vmax*:
ax.pcolor(no_outbreak, cmap=forest_cmap, vmin=0.0, vmax=3.0)

plt.show()

