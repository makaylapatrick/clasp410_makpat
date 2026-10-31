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
Section 2: Constants
'''
#Define the layers of the atmosphere here
nlayers = 3

#Define the emissivity of the Earth here
#The value of emissivity has to be a percentage between 0 and 1
emis = 1     

#Define other important constants here
solar_rad = 1350   #The Solar Constant in W/m^2 
albedo = 0.30       #The Earth's Albedo
stef_bol = 5.67e-8  #The Stefan-Boltzmans constant

#Create a debugging feature as recommended to prevent influence of typos
debug = True

'''
Section 3: Part 1 - Derive an Energy Balance Model
'''
#Create a function for the matrix of coefficients (A)

def energy_balance(nlayers, epsilon):

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

def layered_atmos(nlayers, solar_radiation, albedo, emisivity, earth_surface = 0):

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
    emisivity:
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
    A = energy_balance(nlayers, emis)

    #A check to make sure the matrix develops correctly in the function
    if debug:
        print("This is the matrix created for the flux calc: ")
        print(A)
    
    #For Ax = b we defined A in section 3, and now need to define b
    b = np.zeros((nlayers+1)) #An nlayer column full of zeros

    #Need to explictly define the surface layer value in b
    b[earth_surface] = -(solar_radiation/4)*(1-albedo) #Where Fsolar = -S0/4*(1-a)

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
        temp_calc = ((fluxes[x])/(stef_bol*emisivity)) ** (1/4)
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


