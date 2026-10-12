#!/usr/bin/env python3

'''
This file contains the tools and scripts needed to complete Lab 2
in CLIMATE410: Climate Modeling. To reproduce all plots from the
report follow the steps below!!

The script will be broken up into a series of sections that compiles
each segment of the analysis into one identifiable area.
'''

'''
All necessary inputs, additional libraries, and functions for solving Lorenz
equations can be found in the last two code sections of this file
'''

'''
Section 1: Import Libraries
'''

import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import ode

'''
Section 2: Introduce Initial Conditions
'''

'''
 To replicate results of competion models (Fig 1) set these values to be:
    a = 1, b = 2, c = 1, d = 3
    delta_t = 1
    N1 = 0.3, N2 = 0.6

To replicate species 1 exctinction (Fig 2) set these values to be:
    a = 0, b = 3, c = 0, d = 2
    delta_t = 1
    N1 = 0.3, N2 = 0.6
    
To replicate species 2 exctinction (Fig 2) set these values to be:
    a = 0, b = 1, c = 0, d = 3
    delta_t = 1
    N1 = 0.3, N2 = 0.6

To replicate both species exctintion (Fig 2) set these values to be:
    a = -1, b = 1, c = -1, d = 1
    delta_t = 1
    N1 = 0.3, N2 = 0.6

To replicate species uneasy truce (Fig 2) set these values to be:
     a = 1, b = 1, c = 1, d = 1
    delta_t = 1
    N1 = 0.3, N2 = 0.6

 To replicate the results of predator-prey models (Fig 1) set the values to be:
    a = 1, b = 2, c = 1, d = 3
    delta_t_pp = 0.001
    N1 = 0.3, N2 = 0.6


To replicate the results of predator-prey model phase diagrams (Fig 3a) set the values to be:
    a = 1, b = 1, c =1, d = 2
    delta_t_pp = 0.001
    N1 = 0.3, N2 = 0.6

To replicate the results of predator-prey model phase diagrams (Fig 3b) set the values to be:
    a = 1, b = 2, c =1, d = 1
    delta_t_pp = 0.001
    N1 = 0.3, N2 = 0.6

To replicate the results of predator-prey model phase diagrams (Fig 3c) set the values to be:
    a = 1, b = 2, c =1, d = 2
    delta_t_pp = 0.001
    N1 = 0.3, N2 = 0.6

To replicate the results of predator-prey model phase diagrams (Fig 3d) set the values to be:
    a = 1, b = 1, c =2, d = 1
    delta_t_pp = 0.001
    N1 = 0.3, N2 = 0.6

To replicate the results of predator-prey model phase diagrams (Fig 3e) set the values to be:
    a = 1, b = 1, c =3, d = 1
    delta_t_pp = 0.001
    N1 = 0.3, N2 = 0.6

To replicate the results of predator-prey model phase diagrams (Fig 3f) set the values to be:
    a = 2, b = 1, c =3, d = 1
    delta_t_pp = 0.001
    N1 = 0.3, N2 = 0.6
''' 

a = 1 # Reproduction rate of species 1
b = 2 # Varaible to scale the impact of species 2 on species 1
c = 1 # Reproduction rate of species 2
d = 2 # Varaible to scale the impact of species 1 on species 2

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
    '''
    This function is designed to calculate the Lotka-Volterra competition models, 
    but using the scipy integrater ode package.

    Instead of using the 'dopri5' method as shown in the lab documentation
    I used the 'RK45' method because I found this to provide more reasonable results.
    '''
    # Utilize the ode package to calculate the competition equations
    r = ode(lotka_comp).set_integrator('RK45')
    # Define the inital parameters for use in the ode solver
    r.set_initial_value(N, 0).set_f_params(a, b, c, d)
    # Need to create empty lists the function can append to over time
    pop1 = []
    pop2 = []
    # Need to append the inital conditons into the list as the first step
    pop1.append(pop_init[0]) 
    pop2.append(pop_init[1])
    # Run a loop to calculate competition over time
    while r.successful() and r.t < time:
        r.integrate(r.t + dt)
        # Append the calculated values to the empty lists created earlier
        pop1.append(r.y[0])
        pop2.append(r.y[1])
        
    return[pop1, pop2]

# Run the rk45 competition solver function!
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
fig, (ax1_comp, ax2_comp) = plt.subplots(1, 2, figsize=(12,8))
# Utilize a color blind friendly color scheme for accessibility of figures
plt.style.use('tableau-colorblind10')
# Plot the euler method data on the left hand figure
ax1_comp.plot(time_history_comp, species1, label = "Species 1")
ax1_comp.plot(time_history_comp, species2, label = "Species 2")
ax1_comp.set_xlim(0.0, 100)
ax1_comp.set_ylim(0,1)
ax1_comp.set_xlabel("Time (years)")
ax1_comp.set_ylabel("Population Density (N value)")
ax1_comp.set_title("Euler method")
# Plot the rk45 method data on the right hand figure
ax2_comp.plot(rk45_pop_comp[0], label = "Species 1")
ax2_comp.plot(rk45_pop_comp[1], label = "Species 2")
ax2_comp.set_xlim(0.0, 100)
ax2_comp.set_ylim(0,1)
ax2_comp.set_xlabel("Time (years)")
ax2_comp.set_ylabel("Population Density (N value)")
ax2_comp.set_title("RK45 Solver Method")

fig.suptitle("Lotka-Volterra Competition Model Outputs", fontsize = 20)
ax1_comp.legend()
ax2_comp.legend()
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
# To get an accurate representation of the time in 100 years over 
# very small time steps we need to multiply total_time by 1000.
predation_time = time_total * 1000

for t in range(predation_time):
    # First we run the function we created earlier to calculate 
    # population change with our desired inputs
    pred_prey = lotka_pp(delta_t_pp, pop_change_pp, a, b, c, d)
    # Calculate values for the next time steps
    pop_change_pp = pop_change_pp + pred_prey * delta_t_pp
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
    '''
    This function is designed to calculate the Lotka-Volterra predator-prey models, 
    but using the scipy integrater ode package.
    
    Instead of using the 'dopri5' method as shown in the lab documentation
    I used the 'RK45' method because I found this to provide more reasonable results.
    '''
    # Utilize the ode package to calculate the predator-prey models
    r = ode(lotka_pp).set_integrator('RK45')
    # Define the inital parameters for use in the ode solver
    r.set_initial_value(N, 0).set_f_params(a, b, c, d)
    # Need to create empty lists the function can append to
    prey = []
    predator= []
    # Need to append the inital conditons into the list as the first step
    prey.append(pop_init[0]) 
    predator.append(pop_init[1])

    # Run a loop to calculate competition over time
    while r.successful() and r.t < time:
        r.integrate(r.t + dt)
        # Append the calculated values to the empty lists created earlier
        prey.append(r.y[0])
        predator.append(r.y[1])
        
    return[prey, predator]
    
# Run the rk45 predator-prey solver function!
rk45_pop_pp = rk45_pp(pop_init, a, b, c, d, delta_t_pp, time_total)

'''
Section 10: Plotting the Predtor Prey Models
'''

# Create an array of the population history calculated in the euler loop for easier plotting
pop_history_pp = np.array(pop_history_pp)

# We need to separate the data generated for each species for ease of plotting
prey = pop_history_pp[:,0] # Pull out the values for species 1 from the function
predator = pop_history_pp[:,1] # Pull out the values for species 2 from the function

# Generate the figure size and to display both the euler and rk45 outputs
fig, (ax1_pp, ax2_pp) = plt.subplots(1, 2, figsize=(12,8))
# Utilize a color blind friendly color scheme for accessibility of figures
plt.style.use('tableau-colorblind10')
# Plot the euler method data on the left hand figure
ax1_pp.plot(time_history_pp, prey, label = "Prey")
ax1_pp.plot(time_history_pp, predator, label = "Predator")
ax1_pp.set_ylim(0.2,0.7)
ax1_pp.set_xlabel(f"Time (step of {delta_t_pp} over {time_total} years)")
ax1_pp.set_ylabel("Population Density (N value)")
ax1_pp.set_title("Euler method")
# Plot the rk45 method data on the right hand figure
ax2_pp.plot(rk45_pop_pp[0], label = "Prey")
ax2_pp.plot(rk45_pop_pp[1], label = "Predator")
ax2_pp.set_ylim(0.2,0.7)
ax2_pp.set_xlabel(f"Time (step of {delta_t_pp} over {time_total} years)")
ax2_pp.set_ylabel("Population Density (N value)")
ax2_pp.set_title("RK45 Solver Method")

fig.suptitle("Lotka-Volterra Predator-Prey Model Outputs", fontsize = 20)
ax1_pp.legend()
ax2_pp.legend()
plt.show()

'''
Section 11: Plotting Phase Diagrams
'''
# Generate the figure size and to display both the euler and rk45 outputs
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,8))
# Utilize a color blind friendly color scheme for accessibility of figures
plt.style.use('tableau-colorblind10')
# Plot the euler method data on the left hand figure
ax1.plot(prey, predator)
ax1.set_xlabel("Prey Population Density (N value)")
ax1.set_ylabel("Predator Population Density (N value)")
ax1.set_title("Euler method")
# Plot the rk45 method data on the right hand figure
ax2.set_xlabel("Prey Population Density (N value)")
ax2.set_ylabel("Predator Population Density (N value)")
ax2.plot(rk45_pop_pp[0], rk45_pop_pp[1])
ax2.set_title("RK45 Solver Method")

fig.suptitle("Predator-Prey Phase Diagrams", fontsize = 20)
plt.show()

'''
Section 12: Create a function to define the Lorenz equations
'''

def lorenz(init, t):
    '''
    This function is designed to solve the Lorenz equations for 
    weather prediction. 

    There are seven required inputs for solving these equations:
        t: time_step necessary for argument consistency but is not used in equations
        x: fluid convection rate
        y: horizontal temperature variation
        z: vertical temperature variation
        sigma, beta, and r: defined constants

    The function will return an array of the output of all three equations!

    '''
    # Define the inital condition variables
    x, y, z = init

    # Define the Lorenz equations
    dxdt = - sigma * x + sigma * y
    dydt = - x * z + r *x - y
    dzdt = x * y - beta * z 

    return [dxdt, dydt, dzdt]

'''
Section 13: Using the rk45 solver to solve Lorenz 
'''

'''
After running into many errors with the scipy ode package, I decided to follow a process
as seen in this post for solving Lorenze equations:
https://medium.com/quantaphy/the-butterfly-effect-meet-the-lorenz-attractor-381ed89b98ab

To replicate the results from the Lorenz model (Fig 4a) set the values to be:
    x = 1
    y = 1
    z = 1
    r =23

To replicate the results from the Lorenz model (Fig 4a) set the values to be:
    x = 1.0e-5
    y = 1.0e-5
    z = 1.0e-5
    r =23

To replicate the results from the Lorenz model (Fig 4a) set the values to be:
    x = 1
    y = 1
    z = 1
    r =25 

To replicate the results from the Lorenz model (Fig 4a) set the values to be:
    x = 75
    y = 400000
    z = -0.0000007
    r =25 
'''
# Import necessary aditional packages
from scipy.integrate import odeint

# Define the inital values for the equations
# Inital values for the Lorenz equations
x = 75 # Fluid convection rate
y = 400000 # Horizontal temperature variation
z = -0.0000007 # Vertical temperature variation
inital_conditions = [x, y, z]

r = 25 # Defined constant
sigma = 10 # Defnined constant
beta = 8 / 3 # Defined constant

# Define the time to integrate Lorenze equations over
time_start_lorenz = 0
time_final_lorenz = 100
time_steps_lorenz = np.linspace(time_start_lorenz, time_final_lorenz, time_final_lorenz * 100)

# Use odeint to integrate the Lorenz equations
# This package needs inputs of an inital solver function, the intial condtions in vector form,
# and the time_steps that will be integrated over for solving.
solve_lorenz = odeint(lorenz, inital_conditions, time_steps_lorenz)

# Extract the outputs for each indivual direction
dx = solve_lorenz[:, 0]
dy = solve_lorenz[:, 1]
dz = solve_lorenz[:, 2]

'''
Section 14: Plot the Lorenz Output
'''

# Plot the results in a 3D space because of three directional variables produced
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim((-30, 30))
ax.set_ylim((-30, 30))
ax.set_zlim((0, 50))

ax.plot(dx, dy, dz, color = 'teal')

ax.set_xlabel("Fluid Convection Rate (x-axis)")
ax.set_ylabel("Horizontal Temperature Variation (y-axis)")
ax.set_zlabel("Vertical Temperature Variation (z-axis)")
ax.set_title("Lorenz Equations Solution")

ax.grid(True)
plt.show()















