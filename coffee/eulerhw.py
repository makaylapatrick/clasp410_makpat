#!/usr/bin/env python3

'''
For for this homework assignment we are trying to figure out how to calculate
how long it takes for a coffee cup to cool to the temperature we desire.
We use Euler's Method to run the caluclation within this code.

The follow is my pseudo-code of my plan:
1. Define what we know

temp_init = initial temperature (145)

temp_env =environmental temperature (18)

temp_final = final drinking temperature (80)

cons_k = value of cooling constant k (0.1)

delta_t = estimated time step (10)

tolerance = |temp_curr - temp_env|

2. Define what we don't know
the value of t at future time step when coffee temperature is equal to the final temperature
3. What can we do?
There would be an if/else loop through time steps while the temperature is greater than the final. 
It ends when final coffee temp is reached.
But the temperature will never reach exactly the final so you could define a tolerance
within an absolute value or a range of absolute values the final temperature can be accepted at
Need to create a copy of initial copy cup labeled like: new_cup
the loop will append to the copy until the final temperature is reached
the loop will have a statement that tests to see if the desired final temp range has been reached.
if it is not the loop will go through an additional statement that calculates:

temp_new = temp_curr - delta_t*cons_k*(temp_curr - temp_env)

4. Visualize the data
The final step will be to plot the data as a curve so we can see how the coffee cools over time
'''

#Section 1: Import Libraries

import numpy as np
import matplotlib.pyplot as plt

#Section 2: Define Variables

temp_init = 145     #The initial temperature value for the coffee
temp_env = 22       #The environmental temperature value
temp_final = 80     #The final desired temperature of the coffee

cons_k = 0.1        #A constant inlcuded in Newtons Law of Cooling Eqn.

time_step = 1       #The initial time step
delta_t = 10        #The change in timesteps?

tolerance = abs(temp_curr - temp_env)

#Section 3: Make a function for the cooling equation

def cooling_eqn():
    

if 




#within a if and else statment
temp_curr = #copy of initial

temp_new = #equation

temp_curr = #copy of new 



