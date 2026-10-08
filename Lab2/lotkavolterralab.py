#!/usr/bin/env python3

'''
This file contains the tools and scripts needed to complete Lab 2
in CLIMATE410: Climate Modeling. To reproduce all plots from the
report follow the steps below!!

The script will be broken up into a series of sections that compiles
each segment of the analysis into one identifiable area.
'''

#Section 1: Import Libraries

import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import ode

#Section 2: Introduce Initial Conditions

a = 1 # Reproduction rate of species 1
b = 2 # Varaible to scale the impact of species 2 on species 1
c = 1 # Reproduction rate of species 2
d = 3 # Varaible to scale the impact of species 1 on species 2

time_start = 0 # Initial time step
delta_t = 1 # Time step in years
time_total = 10 # Amount of time the model will iterate over

# Initial population values
species1_init = 0.3  
species2_init = 0.6
# Create an array of the inital population sizes for species 1 and 2
pop_init = np.array([species1_init, species2_init]) 

# Create empty lists to store model output for future graphing purposes
pop_history = [np.copy(pop_init)]
time_history = [0.0]

#Section 3: Create Function for Competition Equations

def lotka_comp(N, a, b, c, d):

    '''
    This function is designed to caluclate the Lotka-Volterra Competition
    equations for two species. 
    There are five required inputs for the equations:
    N: defines the inital population values for species 1 and species 2
    a, c: define the theorhetical repreoduction rates of species 1 and 2 respectively
    b, d: define the impact of each species on each other

    The return will provide an array of the calculated change in population values where 
    dN1dt is the change in population for species 1
    dN2dt is the change in population for species 2
    '''
    #Define the Lotka-Voltera Competition Eqns
    dN1dt = a * N[0] * (1 - N[0]) - b * N[0] * N[1]
    dN2dt= c * N[1] * (1 - N[1]) - d * N[0] * N[1]

    return(np.array(dN1dt, dN2dt))

first_step = lotka_comp(pop_init, a, b, c, d) # Run the function to see the first time step

print(first_step) # A check to make sure reasonable values are being prodcued

#Section 4: Calculating competition with Euler method

# A copy of the initial values must be created initally so that they 
# can be appeneded to thorughout the euler loop
pop_change = np.copy(pop_init) # Generate a copy of the inital population change value 
time_change = np.copy(time_start) # Generate a copy of the inital time step value 

# This for loop represents the euler time series calculation over the total time series
for t in range(time_total):
    # First we run the function we created earlier to calculate 
    # population change with our desired inputs
    competition = lotka_comp(pop_change, a, b, c, d)
    # Calculate how values for the next time steps
    pop_change = pop_change + competition
    time_change = time_change + delta_t
    # Append these new values to the lists of pop_change and time_change
    pop_history.append(pop_change)
    time_history.append(time_change)
    # Print the values calculated to ensure the calculations are running appropriately
    print(pop_history)
    print(time_history)

# Create an array of the population history calculated in the euler loop for easier plotting
pop_history = np.array(pop_history)

# We need to separate the data generated for each species for ease of plotting
species1 = pop_history[:,0] # Pull out the values for species 1 from the function
species2 = pop_history[:,1] # Pull out the values for species 2 from the function

# Section 5: Calculating competition with rk45 solver
'''
def rk45_comp(N, a, b, c, d, dt, time):
    r = ode(lotka_comp).set_integrator('dopri5')

    r.set_initial_value(N, 0).set_f_params(a, b, c, d)

    pop1 = []
    pop2 = []

    while r.successful() and r.t < time:
        r.integrate(r.t + dt)
        pop1.append(r.y[0])
        pop2.append(r.y[1])

    #print(pop1)
    #print(pop2)

    return[pop1, pop2]

rk45_pop_comp = rk45_comp(pop_init, a, b, c, d, delta_t, time_total)
'''
'''
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,8))

plt.style.use('tableau-colorblind10')

ax1.plot(time_history, species1, label = "Species 1")
ax1.plot(time_history, species2, label = "Species 2")

#ax2.plot(rk45_pop_comp[0], label = "Species 1")
#ax2.plot(rk45_pop_comp[1], label = "Species 2")

plt.show()
'''

# Generate a figure for the change in population over time we just calculated!!
plt.figure(figsize=(8,8))
plt.style.use('tableau-colorblind10')
plt.plot(time_history, species1, label = "Competitor 1") # Plot the curve of the 1st species
plt.plot(time_history, species2, label = "Competitor 2") # Plot the curve of the 2nd species

plt.title('Euler Change in Population')
plt.show()




#Section 4
def lotka_pp(t, N, a=1, b=2, c=1, d=3):

    #Define the Lotka-Voltera Predator-Prey Eqns.
    dN1dt = a * N[0] - b * N[0] * N[1]




