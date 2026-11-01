#!/usr/bin/env python3

'''
This file contains the tools and scripts needed to complete Lab 3
in CLIMATE410: Climate Modeling. To reproduce all plots from the
report follow the steps below!!

The script will be broken up into a series of sections that group similar
actions into one identifiable area
'''

''''
--------------------------------
Code Block 1: Building the model
--------------------------------
'''

'''
Section 1: Import Libraries
'''

import numpy as np
import matplotlib.pyplot as plt 

'''
Section 2: Define Constants for Code Block 1
'''
#Create a debugging feature as recommended to prevent influence of typos
debug = False

#Define the layers of the atmosphere here
nlayers = 1

#Define the emissivity of the Earth here
#The value of emissivity has to be a percentage between 0 and 1
emis = 1

#Define other important constants here
solar_rad = 1350   #The Solar Constant in W/m^2 
albedo = 0.30       #The Earth's Albedo
stef_bol = 5.67e-8  #The Stefan-Boltzmans constant

'''
Section 3: Part 1 - Derive an Energy Balance Model
'''
#Create a function for the matrix of coefficients (A)

def energy_balance(nlayers, epsilon):
    if debug:
        print("This is the start of the energy balance function")
    '''
    This function is designed to create the diagonal matrix of coefficients (A)
    from Ax = b model equation.

    nlayers: 
        the amount of atmosphere layers above the surface layer that will be calculated.
    example:
        for nlayers = 3 the expected output is a 4x4 matrix 
        with the surface layer as the first row.
    '''
    #Create the inital empty matrix
    matrix = np.zeros([nlayers+1, nlayers+1])
    #An equation to make calculating the impact of emissivity cleaner
    y = (1 - epsilon)
    #A loop to build each component of the atmosphere matrix
    for i in range(nlayers+1):
        for j in range(nlayers+1):
            if i == j:
                if i == 0:
                    matrix[i, j] = -1 #Set the first value to be -1
                else:
                    matrix[i, j] = -2 #Set the center diagonal to be -2
            elif i == 0:
                if j == 0:
                    matrix[i, j] = 1 #Set the second value on the surface layer to be 1
                else:
                    ex = np.abs(i - j) - 1
                    if debug:
                        print(f"this is the exponent value: {ex} at point {i}, {j}")
                    matrix[i, j] = y**ex #Set all other values of the surface layer
            elif i == j - 1:
                matrix[i, j] = epsilon #Set the values of the diagonal above the center line
            elif i == j + 1:
                matrix[i, j] = epsilon #Set the values of the diagonal below the center line
            else:
                ex = np.abs(i - j) - 1
                matrix[i, j] = epsilon * (y ** ex) #Set all other values of the model 

    return matrix

#Run a check to make sure the model is outputting the expected matrix
if debug:
    check = energy_balance(nlayers, emis)
    print("Making sure the matrix develops correctly:")
    print(check)

'''
Section 4: Part 2 - Write a function for the atmosphere model
'''

#Create a function for the layered atmosphere

def layered_atmos(nlayers, solar_radiation, albedo, emissivity, earth_surface = 0):
    if debug:
        print("This is the start of the layered atmos function")

    '''
    This function will solve the Stefan-Boltzman law to return 
    temperature values for each atmosphere layer.

    inputs:

    nlayers:
        The quantity of layers in the atmosphere model
    solar_radiation:
        The incoming solar radiation from the Sun. Defined above.
    albedo:
        The reflectivity of the Earth's surface. Defined above.
    emissivity:
        The measure of radiation absorbed and emitted by each atmosphere layer.
    A:
        The coefficient matrix generated in the energy_balance function above
    b:
        A zeros coeffecint matrix with the first row defined as the radiation absorption at the surface
    fluxes:
        The result of solving the linear equation for A and b
    
    returns:

    temp:
        The value of temperatures after converting fluxes to temperatures with the Stefan-Boltzman law
    '''

    #Utilize function for A
    A = energy_balance(nlayers, emissivity)

    #A check to make sure the matrix develops correctly in the function
    if debug:
        print("This is the matrix created for the flux calc: ")
        print(A)
    
    #For Ax = b we defined A in section 3, and now need to define b
    b = np.zeros((nlayers+1)) #An nlayer column full of zeros

    #Need to explictly define the surface layer value in b
    b[earth_surface] = -((solar_radiation*(1-albedo))/4) #Where Fsolar = -S0/4*(1-a)

    #Now solve for the flux values
    fluxes = np.linalg.solve(A, b) #An easier method of using linear solvers

    #A debug statment to make sure the calculation was successful
    if debug:
        print("This is the value of 'fluxes': ")
        print(fluxes)

    #Create an empty array for calculating temperature
    temp = np.zeros(len(fluxes))

    #Convert fluxes to temperature values using the Stefan-Boltzman Law
    #The law states F = e*sigma*T^4
    #Need to solve for temperature here
    for x in range(len(fluxes)):
        #Consulted AI to figure out why I was getting surface values of 400 K here
        #I needed to include a calculation to make sure the surface was calculated with the right emissivity
        if x == 0:
            temp_calc = (fluxes[x]/stef_bol) ** (1/4) #Emissivity not included because it should =1 at surface
        else:
            temp_calc = ((fluxes[x])/(stef_bol*emissivity)) ** (1/4)
        temp[x] = temp_calc

    #Have a debug statement in case the values of the fluxes or temperatures seem innacurate
    if debug:
        print("These are the flux values: ", fluxes)
        print("These are the temperatures: ", temp)
        
    return temp 

#A debug statment to ensure the layered_atmos function ran correctly
if debug:
    check2 = layered_atmos(nlayers, solar_rad, albedo, emis)
    print(check2)

'''
Section 5: A validation test for the model
'''
def validation():
    if debug:
        print("This is the start of the validation function")
    #Create a plot to show the matrix of 3 layer atmosphere with an emissivity of 1
    #I used the process outlined here to plot the matrix:
    #https://matplotlib.org/stable/gallery/images_contours_and_fields/matshow.html

    #And the process here to make a better figure for the matrix:
    #https://www.geeksforgeeks.org/python/matplotlib-pyplot-matshow-in-python/

    #Create the matrix
    val_model = energy_balance(nlayers, emis)

    #Plot the matrix
    figure1 = plt.figure()
    ax1 = figure1.add_subplot(111)

    plt_fig1 = ax1.matshow(val_model)
    figure1.colorbar(plt_fig1)

    plt.title("Validation Matrix", size = 15)
    plt.xlabel("Fluxes")
    plt.ylabel("Model Layers")

    plt.show()

if debug:
    validation()

'''
-------------------------
Code Block 2: Experiments
-------------------------
'''

'''
Section 6: Experiment 1 - Single Layer Atmosphere
'''

def ex1():
    if debug:
        print("This is the start of the ex1 function")
    '''
    This function is designed to run experiment one of the lab assignment:
    A single layer atmosphere surface temperature over a range of emissivity values

    Want to see what emissivity the model calculates at a temperature of 288 K
    '''
    #Following code ideas from classmate Jackson

    #Create list of emis values over a range from 0 - 1 with 10 steps
    emis = np.arange(0,1.1,0.1)

    #Create empty list for the surface layer
    surf_temperatures = []
    #Create empty list for the atmosphere layer
    atmos_temperatures = []

    #Run the model over the various emis ranges
    for e in emis:
        print("this is the e value")
        print(e)
        temp = layered_atmos(1, solar_rad, albedo, e, earth_surface=0)
        surf_temperatures.append(temp[0])
        atmos_temperatures.append(temp[1])

    print("These are the surface temp values: ")
    print(surf_temperatures)

    #Make the figure!
    fig1 = plt.figure(1, figsize = (5, 5))
    ax2 = fig1.add_subplot(111)
    #The data is strange so plot both surface and atmos layer for lab purposes
    ax2.plot(emis, surf_temperatures)
    ax2.hlines(288, 0, 1, color = 'red', label = 'expected 288 K')
    ax2.set_xlabel("Emissivities")
    ax2.set_ylabel("Temperature [K]")
    ax2.set_title("Temperature Change with Emissivity Change")
    plt.grid()
    plt.legend(loc = 'lower right')
    plt.show()

#Run experiment 1
ex1()

'''
Section 7: Experiment 2 - How many layers does the Earth need?
'''

def ex2():
    '''
    This function is designed to run experiment two of the lab:
    How many layers do you need to have a surface temperature of 288 K?
    '''
    #Following a different approach to try to avoid any problems that may have been in ex1

    #Define the emissivity from the lab document
    ex2_emis = 0.255

    #Assume a step range of 0 to 10
    steps = 10

    #Create empty lists to append the temp values and N values
    ex2_temps = []
    ex2_layers = []

    #Loop over the steps to determine layers needed
    for t in range(steps):
        N = t + 1 #Increase layer amount during each step
        ex2_calc = layered_atmos(N, solar_rad, albedo, ex2_emis, earth_surface=0)
        #Add values to the empty lists to save for plotting
        ex2_temps.append(ex2_calc[0])
        ex2_layers.append(N)

    if debug:
        print("These are the temps: ")
        print(ex2_temps)
        print("These are the layers: ")
        print(ex2_layers)

    #Plot the data!!
    fig1 = plt.figure(1, figsize = (7, 5))
    ax2 = fig1.add_subplot(111)
    ax2.plot(ex2_layers, ex2_temps)
    ax2.hlines(288, 0, 10, color = 'red', label = 'expected 288 K')
    ax2.set_xlabel("N Layers")
    ax2.set_ylabel("Temperature [K]")
    ax2.set_title("N Layers to Have Surface Temp of 288 K")
    plt.grid()
    plt.legend(loc = 'lower right')
    plt.show()

#Run experiment 2!
ex2()

def plot_alt():
    '''
    This function is designed to plot the altitude for ex2

    It has the exact same format but is using an expected N layer value of 5
    '''

    #Define the emissivity from the lab document
    ex2_emis = 0.255

    #Define the range of layers based on ex2 results
    layers = 5

    #Run the model for a five layer atmos with same emissivity 
    temps = layered_atmos(layers, solar_rad, albedo, ex2_emis, earth_surface=0)

    if debug:
        print("These are the temperatures: ")
        print(temps)

    #Plot the data!!
    fig1 = plt.figure(1, figsize = (7, 5))
    ax2 = fig1.add_subplot(111)
    ax2.plot(temps, np.arange(0, 6)) #Use np.arange here to make the layer list
    ax2.set_xlabel("Temperature [K]")
    ax2.set_ylabel("Atmospheric Layer")
    ax2.set_title("Altitude vs. Temperature")
    plt.grid()
    plt.show()

#Run the altitude plotting code! 
plot_alt()

'''
Section 8: Experiment 3 - Finding the Size of Venus Atmosphere
'''

def ex3():
    '''
    This function is used to find the amount of atmospheric layers on Venus
    '''
    #First define the new constants
    v_solar = 2600 #In W/m^2
    v_emis = 1
    steps = 110
    v_albedo = 0.77 #From NASA

    #Create empty lists to save values to
    venus_temp = []
    venus_layers = []

    #Run through enough steps to find 700 k surface in the dense venus atmosphere
    for l in range(steps):
        v_layers = l + 1
        venus = layered_atmos(v_layers, v_solar, v_albedo, v_emis, earth_surface=0)
        venus_temp.append(venus[0])
        venus_layers.append(v_layers)

    #Plot the data!!
    fig1 = plt.figure(1, figsize = (7, 5))
    ax2 = fig1.add_subplot(111)
    ax2.plot(venus_layers, venus_temp)
    ax2.hlines(700, 0, 110, color = 'red', label = 'expected 700 K')
    ax2.set_xlabel("Temperature [K]")
    ax2.set_ylabel("Atmospheric Layer")
    ax2.set_title("Layers of Venus Atmosphere")
    plt.grid()
    plt.legend(loc = 'lower right')
    plt.show()

#Run the function!
ex3()

'''
Section 9: Experiment 4 - Nuclear War
'''

def ex4():
    '''
    This function is designed to solve the nulcear winter problem
    '''
    #Define new variables for the scenario
    nuc_N = 5
    nuc_solar = 1350 #W/m^2
    nuc_albedo = 0 #No solar radiation will be influenced by the surface so this can be 0
    nuc_emis = 0.5

    #Run the model for a nuclear winter with the layer absorbing all the radiation being the top
    nuclear = layered_atmos(nuc_N, nuc_solar, nuc_albedo, nuc_emis, earth_surface=5)

    #Plot the data!!
    fig1 = plt.figure(1, figsize = (7, 5))
    ax2 = fig1.add_subplot(111)
    ax2.plot(nuclear, np.arange(0, nuc_N +1))
    ax2.set_xlabel("Temperature [K]")
    ax2.set_ylabel("Atmospheric Layer")
    ax2.set_title("Nuclear Winter")
    plt.grid()
    plt.show()

ex4()









    

