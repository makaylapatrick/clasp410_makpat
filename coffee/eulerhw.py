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

temp_coffee = 145   #The initial temperature value for the coffee
temp_env = 22       #The environmental temperature value
temp_final = 80     #The final desired temperature of the coffee

cons_k = 0.1        #A constant inlcuded in Newtons Law of Cooling Eqn.

time_start = 0     #The initial time step
delta_t = 0.1     #The change in timesteps

#Create empty lists to save the output of the elif loop for plotting later
cup_history  = []
time_history = []

#Section 3: Run the loop

#Create a while loop that goes through the cooling equation
while temp_coffee > temp_final:
    #Run the cooling equation to calculate the new value
    temp_change  = - cons_k * delta_t * (temp_coffee - temp_env) 
    #Calculate the next values
    temp_coffee = temp_coffee + temp_change #Add the calculated temperature to the inital tempature
    time_start = delta_t + time_start #Add the time step to the initial time
    #Need to append the values before calculating new ones to not overwrite
    cup_history.append(temp_coffee)
    time_history.append(time_start)
    #Print the values to ensure the updates are successful
    print(cup_history)
    print(time_history)

print("Yay!! You can drink the coffee!")

plt.figure(figsize=(5,5))
plt.plot(time_history, cup_history, color='pink', linewidth='2')

plt.show()
