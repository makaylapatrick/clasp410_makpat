#!/usr/bin/env python3

'''
This file contains the tools and scripts needed to complete Lab 2
in CLIMATE410: Climate Modeling. To reproduce all plots from the
report follow the steps below!!

The script will be broken up into a series of sections 

'''

#Section 1: Import Libraries

import numpy as np
import matplotlib.pyplot as plt 
import scipy as sp

#Section 2: Create Function for Competition Equations

def lotka_comp(t, N, a, b, c, d):

    #Define the Lotka-Voltera Competition Eqns
    dN1dt = a * N[0] * (1 - N[0]) - b * N[0] * N[1]
    dN2dt = c * N[1] * (1 - N[1]) - d * N[0] * N[1]

    return (np.array([dN1dt, dN2dt]))


#Section 3: Introduce Initial Conditions

a = 1 # Reproduction rate of species 1
b = 2 # Varaible to scale the impact of species 2 on species 1
c = 1 # Reproduction rate of species 2
d = 3 # Varaible to scale the impact of species 1 on species 2

delta_t = 1 # Time step in years

#Section 4: Start trying to calculate the competition

time = np.arange(0, 100, delta_t) # Use np.arange to create an array of steps over a time series

N_init = np.array([0.3, 0.6]) # Create an array of the inital population sizes for species 1 and 2

#Create Euler Time Step Here
#Just did not get that far yet oops

competition = lotka_comp(time, N_init, a, b, c, d) # Run the function to calculate how species 1 and 2 compete

print(competition) # A check to make sure reasonable values are being prodcued

#species1 = competition[0] # Pull out the values for species 1 from the function
#species2 = competition[1] # Pull out the values for species 2 from the function

#plt.plot(time, species1, color = "blue", label = "Competitor 1") # Plot the curve of the 1st species
#plt.plot(time, species2, color = "red", label = "Competitor 2") # Plot the curve of the 2nd species






#Section 4
def lotka_pp(t, N, a=1, b=2, c=1, d=3):

    #Define the Lotka-Voltera Predator-Prey Eqns.
    dN1dt = a * N[0] - b * N[0] * N[1]




