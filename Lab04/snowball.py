#!/usr/bin/env python3

'''
This file contains the tools and scripts needed to complete Lab 4
in CLIMATE410: Climate Modeling. To reproduce all plots from the
report follow the steps below!!

The script will be broken up into a series of sections that compiles
each segment of the analysis into one identifiable area.
'''

'''
Section 1: Import Libraries
'''

import numpy as np
import matplotlib.pyplot as plt

'''
Section 2: Define Constants
'''

#Debug statment
debug = True

#Define time in years
time_step = 1
max_time = 10000 #How long the simulation will run for


#Define all constant values
sol_flux = 1370 #Incident solar flux in W/m^2
heat_sea = 4.2e6 #Heat capacity of seawater in J/m^3K^1
dens_sea = 1020 #Density of seawater in kg/m^3
e_rad = 6378 #Earth radius in km
t_diff = 100 #Thermal diffusivity in m^2/s
oc_depth = 50 #Depth of the ocean in m
stef_bol = 5.67e-8 #Stefan-Boltzman constant in J/m^2sK^4
emis = 1 #Emissivity for range [0,1]
a_ground = 0.3 #Albedo of the water/ground
a_snow = 0.6 #Albedo of ice/snow

'''
Section 3: Functions that define different temperature scenarios
'''

#Function from class
def temp_warm(lats_in):
    
    '''
    Create a temperature profile for modern day "warm" earth.
    Parameters
    ----------
    lats_in : Numpy array
    Array of latitudes in degrees where temperature is required.
    0 corresponds to the south pole, 180 to the north.
    Returns
    -------
    temp : Numpy array
    Temperature in Celcius.
    '''

    # Set initial temperature curve
    T_warm = np.array([-47, -19, -11, 1, 9, 14, 19, 23, 25, 25,
    23, 19, 14, 9, 1, -11, -19, -47])
    # Get base grid:
    npoints = T_warm.size
    dlat = 180 / npoints # Latitude spacing.
    lats = np.linspace(dlat/2., 180-dlat/2., npoints) # Lat cell centers.
    # Fit a parabola to the above values
    coeffs = np.polyfit(lats, T_warm, 2)
    # Now, return fitting sampled at "lats".
    temp = coeffs[2] + coeffs[1]*lats_in + coeffs[0] * lats_in**2
    return temp

#Function for creating a very hot world
def temp_hot(lats):
    return np.full(lats.size, 60.0)

#Function for creating a very cold world
def temp_cold(lats):
    return np.full(lats.size, -60.0)

#Funciton from class
def insolation(S0, lats):

    '''
    Given a solar constant (`S0`), calculate average annual, longitude-averaged
    insolation values as a function of latitude.
    Insolation is returned at position `lats` in units of W/m^2.
    Parameters
    ----------
    S0 : float
    Solar constant (1370 for typical Earth conditions.)
    lats : Numpy array
    Latitudes to output insolation. Following the grid standards set in
    the diffusion program, polar angle is defined from the south pole.
    In other words, 0 is the south pole, 180 the north.
    Returns
    -------
    insolation : numpy array
    Insolation returned over the input latitudes.
    '''

    # Constants:
    max_tilt = 23.5 # tilt of earth in degrees
    # Create an array to hold insolation:
    insolation = np.zeros(lats.size)
    # Daily rotation of earth reduces solar constant by distributing the sun
    # energy all along a zonal band
    dlong = 0.01 # Use 1/100 of a degree in summing over latitudes
    angle = np.cos(np.pi/180. * np.arange(0, 360, dlong))
    angle[angle < 0] = 0
    total_solar = S0 * angle.sum()
    S0_avg = total_solar / (360/dlong)
    # Accumulate normalized insolation through a year.
    # Start with the spin axis tilt for every day in 1 year:
    tilt = [max_tilt * np.cos(2.0*np.pi*day/365) for day in range(365)]
    # Apply to each latitude zone:
    for i, lat in enumerate(lats):
        # Get solar zenith; do not let it go past 180. Convert to latitude.
        zen = lat - 90. + tilt
        zen[zen > 90] = 90
        # Use zenith angle to calculate insolation as function of latitude.
        insolation[i] = S0_avg * np.sum(np.cos(np.pi/180. * zen)) / 365.
    # Average over entire year; multiply by S0 amplitude:
    insolation = S0_avg * insolation / 365
    return insolation


'''
Section 4: Matrix functions
'''

def tridiag_matrix(n_lats, dy):
    '''
    Build the K matrix to solve the diffusivity equations

    Inputs
        n_lats: the amount of latitude bands defined across the 0 - 180

        dy: need to account for the distance over each band as seen in class

    Returns
        matrix_k_dy: the tridiagonal matrix K defined in the lecture notes
    '''
  
    #Set the values for the centerline and the diagonal above and below the center
    matrix_k = (
          np.diag(2.0 * np.ones(n_lats))
          + np.diag(-np.ones(n_lats - 1), k=1)
          + np.diag(-np.ones(n_lats - 1), k=-1)
            )

    #Set boundary conditions at north and south poles
    matrix_k[0, 0] = 1.0
    matrix_k[-1, -1] = 1.0

    matrix_k_dy = matrix_k / dy ** 2
            
    return matrix_k_dy    

def first_derivative_matrix(n_lats, dy):
    '''
    Build the B matrix for solving the Snowball Earth model

    Inputs
        n_lats: the amount of latitude bands defined across the 0 - 180
        
        dy: need to account for the distance over each band as seen in class

    Returns
        matrix_by: the matrix B used for solving the spherical correction
    '''

    #Set the values of the matrix
    matrix_b = (
        np.diag(np.ones(n_lats - 1), k=1)
        - np.diag(np.ones(n_lats - 1), k=-1)
        ) 

    #Set boundary conditions
    matrix_b[0, :] = 0.0
    matrix_b[-1, :] = 0.0

    matrix_by = matrix_b / (2.0 * dy)

    return matrix_by


'''
Section 5: Set up the base for the model
'''

def solver_setup(n_lats, temperature_profile, diffusivity):
    '''
    Prepare the grid, initial temperatures, and diffusion matrix.

    Inputs
        n_lats: the amount of latitude bands defined across the 0 - 180
    
        temperature_profile: calls the function for the input conditions as defined by above functions

        diffusivity: needed for solving the L matrix and defined as t_diff above
    
    Returns
        lats: the basic model shape with the values being measured at the center of the lat band

        dy: the grid spacing conversion so that the diffusivity will be accurate across the cell

        temperature: the plot of the temperature profile of the Snowball Earth

        matrix_L: the L matrix used for solving equation 2 seen in the lab assignment

        time_sec: the conversion of time from years into seconds

        n_steps: the total number of time steps taken during the model run
    '''

    #Generate the model of the earth with the layers of each band
    dlat = 180 / n_lats #Latitude spacing
    lats = np.linspace(dlat/2.0, 180-dlat/2.0, n_lats) #Lat cell centers

    #Convert earth radius to m for continuity of units
    e_rad_m = e_rad * 1000
    
    #Calculate grid spacing for accurate diffusivity across cells
    dy = e_rad_m * np.deg2rad(dlat)

    #Define the temperature profile
    temperature = np.array(temperature_profile(lats), dtype=float)

    #Calculate time in seconds for continuity of units
    time_sec = time_step * 3.154e7

    #Total number of temperature updates needed to simulate max_time
    n_steps = int(round(max_time / time_step))

    #Create the identity matrix
    identity = np.diag(np.ones(n_lats))

    #Solve for the matrix of L as shown with equation 2 in the lab assignment
    matrix_L = (identity + diffusivity * time_sec * tridiag_matrix(n_lats, dy))

    return lats, dy, temperature, matrix_L, time_sec, n_steps


'''
Section 6: Solve for the diffusion equation
'''

def diffusion(matrix_L, temperature):
    '''
    Solve the diffusion equation!

    Inputs
        matrix_L: defined in the solver_setup function and needed to solve equation 2

        temperature: the plot of the temperature profile of the Snowball Earth

    Returns
        the solution for matrix_l and temperature
    '''

    return np.linalg.solve(matrix_L, temperature)


'''
Section 7: Adding spherical correction
'''


def spherical(temperature, matrix_by, area_xz, t_diff, dt_sec):
    '''
    Add spherical correction to the model!

    Input
        temperature: the plot of the temperature profile of the Snowball Earth


        matrix_B: the first derivative matrix calculated earlier

        area_gradient: the spherical area to apply temperature over

    Returns:
        the result of adding the spherical term to the model
    '''

    #Solve for spherical term
    temperature_gradient = matrix_by @ temperature
    area_gradient = matrix_by @ area_xz

    return (t_diff * dt_sec * temperature_gradient * area_gradient / area_xz)


'''
Section 8: Adding radiative forcing
'''

def radiative(temperature, insol, emis, a_ground, dt_sec):
    '''
    Solve for the radiative forcing!

    Inputs
        temperature: the plot of the temperature profile of the Snowball Earth
    
        insol: the insolation value

        emis: the emissivity value as defined above

        a_ground: the albedo value of the surface
    
        dt_sec: the time_step duration in seconds

    Return
        the result of calculating the radiative forcing term
    '''

    #Need to use Kelvin to calculate the outgoing radiation
    temperature_kelvin = temperature + 273.15

    #The amount of solar radiation absorbed by the surface
    absorbed_solar = insol * (1.0 - a_ground)

    #The amount of radiation emitted
    outgoing_radiation = (emis * stef_bol * temperature_kelvin**4)

    #heat_sea is volumetric heat capacity, J/(m^3 K)
    #Multiplying by depth gives heat capacity per surface area
    surface_heat_capacity = dens_sea * heat_sea * oc_depth

    return (dt_sec / surface_heat_capacity * (absorbed_solar - outgoing_radiation))

'''
Section 9: Lab Question 1
'''

def question_1(n_lats=18, temperature_profile=temp_warm):
    '''
    Calculate all four of the model options as seen in the lab assignment model output plot

    Inputs
        n_lats=18: to have 18 latitude bands to calculate across

        temperature_profile: solving given the warm earth scenario provided from class notes

    Returns
        basic_temperature: the value for the diffusivity only plot

        spherical_temperature: the value for the spherical correction plot

        full_temperature: the value for the radiative forcing plot
    '''

    #Use the setup function already defined in Section 5
    lats, dy, initial_temperature, matrix_L, dt_sec, n_steps = (solver_setup(n_lats, temperature_profile, t_diff))

    #Prepare matrix_b
    matrix_B = first_derivative_matrix(n_lats, dy)

    #Convert earth radius to m for continuity of units
    e_rad_m = e_rad * 1000.0

    #Need to ensure solving is being done along the cross sectional area
    area_xz = (2.0 * np.pi * e_rad_m * oc_depth * np.sin(np.deg2rad(lats)))

    #Calculate solar forcing
    insol = insolation(sol_flux, lats)

    #Make arrays for all data output to be saved for plotting
    basic_temperature = initial_temperature.copy()
    spherical_temperature = initial_temperature.copy()
    full_temperature = initial_temperature.copy()

    for step in range(n_steps):

        #Red line: basic diffusion only
        basic_temperature = diffusion(matrix_L,basic_temperature)

        # Gold line: diffusion + spherical correction
        spherical_increment = spherical(spherical_temperature, matrix_B, area_xz, t_diff, dt_sec)

        spherical_temperature = diffusion(matrix_L, spherical_temperature + spherical_increment)

        # Green line: all three processes

        full_spherical_increment = spherical(full_temperature,matrix_B,area_xz,t_diff,dt_sec)

        radiation_increment = radiative(full_temperature,insol,emis,a_ground,dt_sec)

        full_temperature = diffusion(matrix_L,full_temperature + full_spherical_increment + radiation_increment)


    # Convert the internal 0–180 coordinate to -90–90 latitude.
    plot_lats = lats - 90.0

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(plot_lats,initial_temperature,"--",color="blue",label="Initial Condition")
    ax.plot(plot_lats,basic_temperature,color="tab:red",label="Basic Diffusion")
    ax.plot(plot_lats,spherical_temperature,"-.",color="goldenrod",label="Diffusion + Spherical Correction")
    ax.plot(plot_lats,full_temperature,":",color="tab:green",linewidth=2,label="Diffusion + Spherical Correction + Radiation")
    ax.set_xlabel("Latitude (degrees)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("Question 1: Snowball Earth Model Validation")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    plt.legend()
    plt.grid()
    plt.show()

    return (basic_temperature, spherical_temperature, full_temperature)

'''
Section 10: Full Model for Questions 2, 3, and 4
'''

def full_model(n_lats=18,temperature_profile=temp_warm,initial_temperature=None,diffusivity=t_diff,
               emissivity=emis,dynamic_albedo=False,fixed_albedo=a_ground,gamma=1.0):
    '''
    Run the full Snowball Earth model with all components

    Inputs
        n_lats=18: creating 18 latitude bands across the model earth

        temperature_profile: which temperature condiiton being used for solving

        initial_temperature: optional array from a previous simulation. If supplied, it overrides temperature_profile

        diffusivity: defined as t_diff in section 2

        emissivity: defined as emis in section 2

        dynamic_albedo: needed for when albedo values can be fluctuating

        fixed_albedo: defined as the ground albedo value from section 2

        gamma: needed for calculating

    Return
        lats: latitude-band center

        temperature: final temperatures in Celsius

        final_change: largest temperature change over the final approximately 100 years
    '''

    lats, dy, temperature, matrix_L, dt_sec, n_steps = (solver_setup(n_lats, temperature_profile, diffusivity))

    #Previous-run temperatures take precedence if supplied
    if initial_temperature is not None:
        temperature = np.array(initial_temperature, dtype=float)

    if temperature.shape != (n_lats,):
        raise ValueError("Provide one temperature per latitude band.")

    #Set up matrix_B for solving
    matrix_B = first_derivative_matrix(n_lats, dy)

    #Convert earth radius to m for continuity of units
    e_rad_m = e_rad * 1000

    #Need to ensure solving is being done along the cross sectional area
    area_xz = (2.0 * np.pi * e_rad_m * oc_depth * np.sin(np.deg2rad(lats)))

    #Calculate solar forcing
    insol = gamma * insolation(sol_flux, lats)

    #Precompute the implicit update operator to avoid repeated solving
    diffusion_operator = diffusion(matrix_L, np.eye(n_lats))
    check_steps = min(n_steps, max(1, int(round(100.0 / time_step))))
    reference_temperature = temperature.copy()


    for step in range(n_steps):

        #Making sure the right temperature is being used
        if step == n_steps - check_steps:
            reference_temperature = temperature.copy()

        #Allow for albedo to be updated from current conditions
        if dynamic_albedo:
            albedo = np.where(temperature <= -10.0, a_snow, a_ground)

        else:
            albedo = fixed_albedo

        #Calculte spherical correction
        spherical_increment = spherical(temperature, matrix_B, area_xz, diffusivity, dt_sec)

        #Calculate radiative forcing 
        radiation_increment = radiative(temperature, insol, emissivity, albedo, dt_sec)

        #Calculate the final model with all parts included
        temperature = diffusion_operator @ (temperature + spherical_increment + radiation_increment)

    #Calculate the final temperature change
    final_change = np.max(np.abs(temperature - reference_temperature))

    return lats, temperature, final_change


def latitude_weights(lats):
    '''
    Need to have relative surface areas of equal-angle latitude bands
    
    Inputs
        lats: the latitude band center calculated from full_model

    Returns
        weighted latitude solution
    '''
    #Total latitude spacing
    dlat = 180.0 / lats.size

    #Define the edge cases
    lower_edges = np.deg2rad(lats - dlat / 2.0)
    upper_edges = np.deg2rad(lats + dlat / 2.0)

    #Calculate the weights
    weights = np.cos(lower_edges) - np.cos(upper_edges)

    return weights / weights.sum()


def global_temperature(temperature, lats):
    '''
    Calculate area-weighted global mean temperature

    Inputs
        temperature: the temperature as calculated from the full_model function

        lats: the latitude band center as calculated from the full_model function

    Returns
        the area_weighted temperature
    '''

    return np.sum(latitude_weights(lats) * temperature)


def temperature_error(temperature, target, lats):
    '''
    This funciton is used for assessing how closely the model follows the warm Earth conditions

    Inputs
        temperature: the temperature as calculated from the full_model function

        target: the target temperature

        lats: the latitude band center as calculated from the full_model function 

    Returns
        a numerical score for comparing outputs in question 2
    
    '''

    return np.sqrt(np.sum(latitude_weights(lats)* (temperature - target)**2))


'''
Section 11: Lab Question 2
'''

def question_2(n_lats=18):
    '''
    Inputs
        n_lats=18: defined latitude bands for the snowball Earth model

    Returns
        best_diffusivity: the diffusivity value that most closely aligns with the warm Earth conditions

        best_emissivity: the emissivity value that most closely aligns with the warm Earth conditions
    '''

    #Define diffusivity and emissivity values to calculate over
    diffusivity_values = np.linspace(0.0, 150.0, 16)
    emissivity_values = np.linspace(0.0, 1.0, 21)

    #Want to have only a less than 0.01 change
    steady_tolerance = 0.01

    #Lists to hold the error outputs
    diffusivity_errors = []
    emissivity_errors = []
    nonsteady_runs = []

    # -------------------------------------------------------
    # A: Vary diffusivity; hold emissivity at 1.
    # -------------------------------------------------------

    for value in diffusivity_values:

        #Run the full model to calculate lats, temperature, and change based on constant emis
        lats, temperature, change = full_model(n_lats=n_lats,diffusivity=value,
                                               emissivity=1.0,fixed_albedo=0.3)

        #Append calculations to diffusivity error list
        diffusivity_errors.append(temperature_error(temperature, temp_warm(lats), lats))

        #Allow for the changes that are above tolerance to be saved
        if change >= steady_tolerance:
            nonsteady_runs.append(("diffusivity", float(value)))

    # -------------------------------------------------------
    # B: Vary emissivity; hold diffusivity at 100.
    # -------------------------------------------------------

    for value in emissivity_values:

        #Run the full model to calculate lats, temperature, and change based on constant t_diff
        lats, temperature, change = full_model(n_lats=n_lats,diffusivity=100.0,
                                               emissivity=value,fixed_albedo=0.3)

        #Append calculations to emissivity error list
        emissivity_errors.append(temperature_error(temperature, temp_warm(lats), lats))

        #Allow for the changes that are above tolerance to be saved
        if change >= steady_tolerance:
            nonsteady_runs.append(("emissivity", float(value)))


    #Plot to show the change in values over time
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    axes[0].plot(diffusivity_values,diffusivity_errors,"o-")
    axes[0].set_xlabel("Diffusivity (m²/s)")
    axes[0].set_ylabel("Warm-profile RMSE (°C)")
    axes[0].set_title("Emissivity fixed at 1")

    axes[1].plot(emissivity_values,emissivity_errors,"o-")
    axes[1].set_xlabel("Emissivity")
    axes[1].set_ylabel("Warm-profile RMSE (°C)")
    axes[1].set_title("Diffusivity fixed at 100 m²/s")

    for ax in axes:
        ax.grid(True, alpha=0.3)

    fig.tight_layout()

    #Plot to show actual profiles to help interpret the effects
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    #Select value increments to plot changing diffusivity over
    for value in [0.0, 50.0, 100.0, 150.0]:
        lats, temperature, _ = full_model(n_lats=n_lats,diffusivity=value,emissivity=1.0)
        axes[0].plot(lats - 90.0,temperature,label=f"Diffusivity = {value:g}")

    #Select value increments to plot changing emissivity over
    for value in [0.25, 0.5, 0.75, 1.0]:
        lats, temperature, _ = full_model(n_lats=n_lats,diffusivity=100.0,emissivity=value)
        axes[1].plot(lats - 90.0,temperature,label=f"Emissivity = {value:g}")

    #Plot the target condition values to make the comparison easier to see
    for ax in axes:
        ax.plot(lats - 90.0,temp_warm(lats),"k--",label="Warm-Earth target")
        ax.set_xlabel("Latitude (degrees)")
        ax.set_ylabel("Temperature (°C)")
        ax.grid(True, alpha=0.3)
        ax.legend()

    axes[0].set_title("Effect of diffusivity")
    axes[1].set_title("Effect of emissivity")
    fig.tight_layout()

    # -------------------------------------------------------
    # C: Isolate the curve that best fits the warm Earth
    # -------------------------------------------------------

    #Used to find the best fitting lines
    best_error = np.inf
    best_diffusivity = None
    best_emissivity = None
    best_temperature = None

    joint_emissivities = np.linspace(0.05, 1.0, 20)

    #Loop thorugh the diffusivity values and the emissivity values to find the best fitting
    for diffusivity in diffusivity_values:
        for emissivity in joint_emissivities:

            lats, temperature, change = full_model(n_lats=n_lats,diffusivity=diffusivity,
                                                   emissivity=emissivity,fixed_albedo=0.3)

            #Ensure the value follows the tolerance conditions
            if change >= steady_tolerance:
                continue

            #Calculate the level of error
            error = temperature_error(temperature,temp_warm(lats),lats)

            #Identify an error value that meets the criteria defined above
            if error < best_error:
                best_error = error
                best_diffusivity = float(diffusivity)
                best_emissivity = float(emissivity)
                best_temperature = temperature.copy()


    #Create a plot to compare the best error value to the initial warm Earth conditions
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(lats - 90.0,temp_warm(lats),"k--",label="Warm-Earth target")
    ax.plot(lats - 90.0,best_temperature,label="Best coarse-grid fit")

    ax.set_xlabel("Latitude (degrees)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("Question 2: Tuning the Warm-Earth Model")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    return best_diffusivity, best_emissivity

'''
Section 12: Lab Question 3
'''

def question_3(diffusivity, emissivity, n_lats=18):
    '''
    Create and compare hot, cold, and  flash freeze scenario

    Inputs
        diffusivity: as defined from t_diff above

        emissivity: as defined from emis above

        n_lats=18: the 18 latitude bands needed for the model configuration

    Returns
        hot_temperature: the valuea calculated for the hot earth model

        cold_temperature: the values calculated for the cold earth model

        frozen_temperature: the values calculated for the frozen earth model
    
    '''

    #Hot start with dynamic albedo
    lats, hot_temperature, hot_change = full_model(n_lats=n_lats,temperature_profile=temp_hot,
                                                   diffusivity=diffusivity,emissivity=emissivity,dynamic_albedo=True)

    #Cold start with dynamic albedo
    #Use _ as a throwaway variable
    _, cold_temperature, cold_change = full_model(n_lats=n_lats,temperature_profile=temp_cold,
                                                  diffusivity=diffusivity,emissivity=emissivity,dynamic_albedo=True)

    #Warm profile with ice albedo held fixed
    #Use _ as a throwaway variable
    _, frozen_temperature, frozen_change = full_model(n_lats=n_lats,temperature_profile=temp_warm,
                                                      diffusivity=diffusivity,emissivity=emissivity,
                                                      dynamic_albedo=False,fixed_albedo=0.6)


    #Plot all the outputs together to compare
    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(lats - 90.0,temp_warm(lats),"k--",label="Warm-Earth reference")
    ax.plot(lats - 90.0,hot_temperature,label="Hot start; dynamic albedo")
    ax.plot(lats - 90.0,cold_temperature,"-.",label="Cold start; dynamic albedo")
    ax.plot(lats - 90.0,frozen_temperature,":",linewidth=2,label="Warm start; fixed ice albedo")

    ax.set_xlabel("Latitude (degrees)")
    ax.set_ylabel("Temperature (°C)")
    ax.set_title("Question 3: Ice-Albedo Feedback")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    return hot_temperature, cold_temperature, frozen_temperature

'''
Section 13: Lab Question 4
'''

def question_4(diffusivity, emissivity, n_lats=18):
    '''
    Increase gamma from 0.4 to 1.4, then decrease it to 0.4.
    Need to reuse the previous final temperature profile at each step.

    Inputs
        diffusivity: as defined by t_diff above

        emissivity: as defined by emis above

        n_lats=18: the 18 latitude bands for the model grid

    Returns
        gamma_up: values calculated while gamma is increasing

        warming_means: the mean temperature values during increasing gamma

        gamma_down: values calculated while gamma is decreasing

        cooling_means: the mean temperature values during decreaing gamma
    '''

    #Need to use integer to avoid endpoint rounding error
    gamma_up = np.arange(8, 29) * 0.05
    gamma_down = gamma_up[::-1]

    #Empty lists to save calculated mean values
    warming_means = []
    cooling_means = []
    nonsteady_runs = []

    #First run starts from cold Earth
    previous_temperature = None

    #Solve for scenarios where gamma is going up
    for gamma in gamma_up:
        lats, temperature, change = full_model(n_lats=n_lats,temperature_profile=temp_cold,
                                               initial_temperature=previous_temperature,diffusivity=diffusivity,
                                               emissivity=emissivity,dynamic_albedo=True,gamma=float(gamma))

        #Append the warming mean values
        warming_means.append(global_temperature(temperature, lats))

        #Need to start from previous_temperature for each iteration
        previous_temperature = temperature.copy()

        #Only want changes that are greater than 0.01
        if change >= 0.01:
            nonsteady_runs.append(("increasing", float(gamma)))

    #Append the cooling mean values
    cooling_means.append(warming_means[-1])

    #Solve for scenarios where gamma is going down
    for gamma in gamma_down[1:]:
        lats, temperature, change = full_model(n_lats=n_lats,initial_temperature=previous_temperature,diffusivity=diffusivity,
                                               emissivity=emissivity,dynamic_albedo=True,gamma=float(gamma))

        #Append cooling values
        cooling_means.append(global_temperature(temperature, lats))

        #Need to start from previous_temperature for each iteration
        previous_temperature = temperature.copy()

        #Only want changes that are greater than 0.01
        if change >= 0.01:
            nonsteady_runs.append(("decreasing", float(gamma)))

    #Compare values against each other
    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(gamma_up,warming_means,"o-",label="Increasing solar forcing")
    ax.plot(gamma_down,cooling_means,"s--",label="Decreasing solar forcing")

    ax.set_xlabel("Solar multiplier γ")
    ax.set_ylabel("Area-weighted global mean temperature (°C)")
    ax.set_title("Question 4: Solar Forcing")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()


    return gamma_up, warming_means, gamma_down, cooling_means


'''
Section 14: Actually running the experiments
'''
#Question 1: validate the three model stages.
basic, spherical_result, full = question_1(18)

#Question 2: choose tuned parameters.
tuned_diffusivity, tuned_emissivity = question_2(18)

#Question 3: use the tuned parameters.
question_3(tuned_diffusivity,tuned_emissivity,n_lats=18)

#Question 4: use the same tuned parameters.
question_4(tuned_diffusivity,tuned_emissivity,n_lats=18)

plt.show()