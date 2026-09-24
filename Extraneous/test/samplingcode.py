
'''


#Define constants for the different states of the forest
burnt_bare = 1 #depicted as a tan color in the figures to show bare ground
alive_forest = 2 #depicted as a green color in the figures to show trees
on_fire = 3 #depicted as a red color in the figures to show flames

p_spread = 1.0 # Chance to spread to adjacent cells.
p_bare = 0.0 # Chance of cell to start as a bare patch.
p_ignite = 0.0 # Chance of cell to start on fire.

#Implement time as a varaible for the model
time_step = 4

nx, ny = 3, 3

def make_forest():
    forest = np.zeros([ny,nx])
    return forest

ori_forest = make_forest()

print(ori_forest)

ori_forest[:, 1] =1
ori_forest[2, :] =1

print(ori_forest)


nx, ny = 10, 10 # Number of cells in X and Y direction

a, b = [ny-2], [nx-2] # Grid cell for initial burning

#creating the main forest array that the rest of the model will build from
ori_forest = np.zeros([ny, nx]) + 2

#setting the outside grid to be 1, so the fire will not spread into the ghost cells

ori_forest[[ny-1], :] = 1
ori_forest[:, [nx-1]] = 1
ori_forest[0, :] = 1
ori_forest[:, 0] = 1

print (ori_forest)

#Randomly generate forest with bare and forested cells
for i in range(1,nx-1): #Set range to avoid ghost nodes
    for j in range(1,ny-1):
        # Roll our "dice" to see if we get a bare spot:
        if np.random.rand() < p_bare:
            ori_forest[j, i] = burnt_bare # 1 is a bare spot randomly generated
        elif np.random.rand() < p_ignite:
            ori_forest[j,i] = on_fire

#Create a copy of the original forest to begin burning

ori_forest[a,b] = on_fire

burn_forest = np.copy(ori_forest)

print(burn_forest)

'''

'''
# Define time steps (e.g., 0 to 9)
time_steps = range(10)

# Create an empty list to save outputs
outputs = []

# Initial state or value
current_value = 1

for t in time_steps:
  # Simulate a step: update the value (e.g., double it each step)
  current_value = current_value * 2

  # Save the output of this step
  outputs.append(current_value)

print(outputs)
# Output: [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

'''
'''
import numpy as np

no_outbreak = np.zeros([ny, nx]) + 2

deceased = 0 #depicted as a black color in the figures to show those who perished
immune = 1 #depicted as a tan color in the figures to show those who begin immune
healthy = 2 #depicted as a green color in the figures to show healthy individuals
infected = 3 #depicted as a red color in the figures to show where the virus has spread

p_spread = 1.0 # Chance to spread to adjacent cells
p_immune = 0.2 # Chance of person to start as immune
p_zombie = 0.05 # Chance of infection to begin in a cell
p_fatal = 0.1 # Chance of infected to perish instead of recover

def chance_perish():
    if no_outbreak[j, i] == infected:
        if np.random.rand() < p_fatal:
            infected_pop[j, i] = deceased
            print("You Died :(")
'''     


import numpy as np
import matplotlib.pyplot as plt

# --- SIMULATED DATA SETUP ---
# Let's assume you have a 3D numpy array representing your simulation history.
# Shape: (num_time_steps, grid_rows, grid_cols)
num_steps = 50
rows, cols = 100, 100
total_cells = rows * cols

# Creating dummy history data where the virus spreads over time
# In your real code, this 'grid_history' comes directly from your simulation loop.
grid_history = np.zeros((num_steps, rows, cols))
for t in range(num_steps):
    # Simulating a growing infected zone over time for demonstration
    infected_radius = int((t / num_steps) * (rows // 2))
    grid_history[t, 50-infected_radius:50+infected_radius, 50-infected_radius:50+infected_radius] = 1

# --- STEP 1: CALCULATE PERCENTAGE OVER TIME ---
# Define what value represents "infected" in your grid (e.g., 1)
INFECTED_VALUE = 1

percentages = []
for t in range(num_steps):
    # Count how many cells equal the infected value at time step t
    infected_count = np.sum(grid_history[t] == INFECTED_VALUE)
    
    # Calculate percentage
    pct = (infected_count / total_cells) * 100
    percentages.append(pct)

# --- STEP 2: PLOT THE LOGISTIC/EXPONENTIAL CURVE ---
time_steps = np.arange(num_steps)

plt.figure(figsize=(9, 5))
plt.plot(time_steps, percentages, color='crimson', linewidth=2.5, label='Infected %')

# Formatting the Chart
plt.title('Zombie Virus Spread Dynamics Over Time', fontsize=14, fontweight='bold')
plt.xlabel('Time Step (t)', fontsize=12)
plt.ylabel('Percentage of Population Infected (%)', fontsize=12)
plt.ylim(-5, 105) # Keeps the scale bounded between 0% and 100%
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper left')

pct_susceptible = [np.sum(grid_history[t] == 0) / total_cells * 100 for t in range(num_steps)]
pct_infected    = [np.sum(grid_history[t] == 1) / total_cells * 100 for t in range(num_steps)]
pct_zombie      = [np.sum(grid_history[t] == 2) / total_cells * 100 for t in range(num_steps)]

plt.plot(time_steps, pct_susceptible, color='royalblue', label='Susceptible Humans')
plt.plot(time_steps, pct_infected, color='orange', label='Infected (Incubating)')
plt.plot(time_steps, pct_zombie, color='crimson', label='Active Zombies')

# Display the plot
plt.tight_layout()
plt.show()


