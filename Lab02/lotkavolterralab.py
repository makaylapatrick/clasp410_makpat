#!/usr/bin/env python3

'''
This file contains the tools and scripts needed to complete Lab 2
in CLIMATE410: Climate Modeling. To reproduce all plots from the
report follow the steps below!!

The script will be broken up into a series of sections that compiles
each segment of the analysis into one identifiable area.
'''


'''
Section 1: Import Libraries
'''

import numpy as np
import matplotlib.pyplot as plt ;#plt.ion()
from scipy.integrate import ode

'''
Section 2: Introduce Initial Conditions
'''

'''
 To replicate results of competion models set these values to be:
    a = 1, b = 2, c = 1, d = 3

 To replicate the results of predator-prey models set the values to be:
    Replicating species 1 exctinction
    a =
''' 

a = 1 # Reproduction rate of species 1
b = 2 # Varaible to scale the impact of species 2 on species 1
c = 1 # Reproduction rate of species 2
d = 3 # Varaible to scale the impact of species 1 on species 2

time_start = 0 # Initial time step
time_total = 100 # Amount of time the model will iterate over
delta_t_comp = 1 # Time step in years for the competition model
delta_t_pp = 0.001 # Time step in years for the predator-prey model

# Initial population values
species1_init = 0.3  
species2_init = 0.6
# Create an array of the inital population sizes for species 1 and 2
pop_init = np.array([species1_init, species2_init])

# Create empty lists to store model output for future plotting
# Competition model outputs
pop_history_comp = [np.copy(pop_init)]
time_history_comp = [0.0]
# Predator-Prey model outputs
pop_history_pp = [np.copy(pop_init)]
time_history_pp = [0.0]

# A copy of the initial values must be created so that they 
# can be appeneded to thorughout the euler loop without changing 
# the inital values that are applied elsewhere
# Competition input 
pop_change_comp = np.copy(pop_init) # Generate a copy of the inital population change value 
time_change_comp = np.copy(time_start) # Generate a copy of the inital time step value 
# Predator-Prey input
pop_change_pp = np.copy(pop_init) # Generate a copy of the inital population change value 
time_change_pp = np.copy(time_start)# Generate a copy of the inital time step value 

'''
Section 3: Create function for competition equations
'''

def lotka_comp(t, N, a, b, c, d):
    '''
    This function is designed to caluclate the Lotka-Volterra Competition
    equations for two species. 
    
    There are six required inputs for the equations:
        t: time_step necessary for argument consistency but is not used in equations
        N: defines the inital population values for species 1 and species 2
        a, c: define the theorhetical reproduction rates of species 1 and 2 respectively
        b, d: define the impact of each species on each other

    The return will provide an array of the calculated change in population values where 
        dN1dt is the change in population for species 1
        dN2dt is the change in population for species 2
    '''
    #Define the Lotka-Voltera Competition Eqns
    dN1dt = a * N[0] * (1 - N[0]) - b * N[0] * N[1]
    dN2dt= c * N[1] * (1 - N[1]) - d * N[0] * N[1]

    return(np.array([dN1dt, dN2dt]))

# Run the function to see the first time step
first_step_comp = lotka_comp(delta_t_comp, pop_init, a, b, c, d) 
# A check to make sure reasonable values are being prodcued
print(f"This is the population cahnge at the first time step of the competition model: {first_step_comp}") 

'''
Section 4: Calculating competition with Euler method
'''

# This for loop represents the euler time series calculation over the total time series
for t in range(time_total):
    # First we run the function we created earlier to calculate 
    # population change with our desired inputs
    competition = lotka_comp(delta_t_comp, pop_change_comp, a, b, c, d)
    # Calculate values for the next time steps
    pop_change_comp = pop_change_comp + competition * delta_t_comp
    time_change_comp = time_change_comp + delta_t_comp
    # Append these new values to the lists of pop_change and time_change
    pop_history_comp.append(pop_change_comp)
    time_history_comp.append(time_change_comp)
    # Uncomment to print the values calculated to ensure the calculations are running appropriately
    #print(pop_history_comp)
    #print(time_history_comp)

'''
Section 5: Calculating competition with rk45 solver
'''

def rk45_comp(N, a, b, c, d, dt, time):
    r = ode(lotka_comp).set_integrator('RK45', atol = 1.0e-14, rtol = 1.0e-12)

    r.set_initial_value(N, 0).set_f_params(a, b, c, d)

    # Need to create empty lists the function can append to
    pop1 = []
    pop2 = []
    # Need to append the inital conditons into the list as the first step
    pop1.append(pop_init[0]) 
    pop2.append(pop_init[1])

    # 
    while r.successful() and r.t < time:
        r.integrate(r.t + dt)
        pop1.append(r.y[0])
        pop2.append(r.y[1])
        
    return[pop1, pop2]

rk45_pop_comp = rk45_comp(pop_init, a, b, c, d, delta_t_comp, time_total)

'''
Section 6: Plotting the competition models
'''

# Create an array of the population history calculated in the euler loop for easier plotting
pop_history_comp = np.array(pop_history_comp)

# We need to separate the data generated for each species for ease of plotting
species1 = pop_history_comp[:,0] # Pull out the values for species 1 from the function
species2 = pop_history_comp[:,1] # Pull out the values for species 2 from the function

# Generate the figure size and to display both the euler and rk45 outputs
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,8))
# Utilize a color blind friendly color scheme for accessibility of figures
plt.style.use('tableau-colorblind10')
# Plot the euler method data on the left hand figure
ax1.plot(time_history_comp, species1, label = "Species 1")
ax1.plot(time_history_comp, species2, label = "Species 2")
ax1.set_title("Euler method")
# Plot the rk45 method data on the right hand figure
ax2.plot(rk45_pop_comp[0], label = "Species 1")
ax2.plot(rk45_pop_comp[1], label = "Species 2")
ax2.set_title("RK45 Solver Method")

plt.show()

'''
Section 7: Create function for Predator Prey Equations
'''

def lotka_pp(t, N, a, b, c, d):
    '''
    This function is designed to caluclate the Lotka-Volterra
    Predator-Prey equations for two species, where one predates the other. 

    There are six required inputs for the equations:
        t: time_step necessary for argument consistency but is not used in equations
        N: defines the inital population values for species 1 and species 2
        a, c: define the theorhetical reproduction rates of species 1 and 2 respectively
        b, d: define the impact of each species on each other
    
    The return will provide an array of the calculated change in population values where 
        dN1dt is the change in population for species 1
        dN2dt is the change in population for species 2
    '''
    #Define the Lotka-Voltera Predator-Prey Eqns.
    dN1dt = a * N[0] - b * N[0] * N[1]
    dN2dt = - c * N[1] + d * N[0] * N[1]

    return(np.array([dN1dt, dN2dt]))

# Run the function to see the first time step
first_step_pp = lotka_pp(delta_t_pp, pop_init, a, b, c, d) 
# A check to make sure reasonable values are being prodcued
print(f"This is the population change at the first time step of the predator-prey model: {first_step_pp}") 

'''
Section 8: Calculating predation with Euler method
'''

for t in range(time_total*1000):
    # First we run the function we created earlier to calculate 
    # population change with our desired inputs
    pred_prey = lotka_pp(delta_t_pp, pop_change_pp, a, b, c, d)
    # Calculate values for the next time steps
    pop_change_pp = pop_change_pp + pred_prey*delta_t_pp
    time_change_pp = time_change_pp + delta_t_pp
    # Append these new values to the lists of pop_change and time_change
    pop_history_pp.append(pop_change_pp)
    time_history_pp.append(time_change_pp)
    if time_change_pp>100:
        break
    # Uncomment to print the values calculated to ensure the calculations are running appropriately
    #print(pop_history_pp)
    #print(time_history_pp)

'''
Section 9: Calculating predation with rk45 method
'''

def rk45_pp(N, a, b, c, d, dt, time):
    r = ode(lotka_pp).set_integrator('RK45')

    r.set_initial_value(N, 0).set_f_params(a, b, c, d)

    # Need to create empty lists the function can append to
    prey = []
    predator= []
    # Need to append the inital conditons into the list as the first step
    prey.append(pop_init[0]) 
    predator.append(pop_init[1])

    # 
    while r.successful() and r.t < time:
        r.integrate(r.t + dt)
        prey.append(r.y[0])
        predator.append(r.y[1])
        
    return[prey, predator]

rk45_pop_pp = rk45_pp(pop_init, a, b, c, d, delta_t_pp, time_total)

# Create an array of the population history calculated in the euler loop for easier plotting
pop_history_pp = np.array(pop_history_pp)

# We need to separate the data generated for each species for ease of plotting
prey = pop_history_pp[:,0] # Pull out the values for species 1 from the function
predator = pop_history_pp[:,1] # Pull out the values for species 2 from the function

# Generate the figure size and to display both the euler and rk45 outputs
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,8))
# Utilize a color blind friendly color scheme for accessibility of figures
plt.style.use('tableau-colorblind10')
# Plot the euler method data on the left hand figure
ax1.plot(time_history_pp, prey, label = "Species 1")
ax1.plot(time_history_pp, predator, label = "Species 2")
# Plot the rk45 method data on the right hand figure
ax2.plot(rk45_pop_pp[0], label = "Species 1")
ax2.plot(rk45_pop_pp[1], label = "Species 2")

plt.show()

















