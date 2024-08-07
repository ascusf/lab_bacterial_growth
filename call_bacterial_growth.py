#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 22:06:25 2024
This script provides two options for simulating bacterial growth: one using 
user-provided data and another using default parameters. It calls the main script
bacterial_growth.py and runs it.

Option 1: Enter Parameters Without Data
If no data is provided (test == None), users need to input the required parameters:

Option 2: Enter Parameters with Data
If data is provided (test != None), users need to input the required parameters:

Summary of User Inputs
Initial OD600 Reading (initial_read)
The initial optical density at 600 nm, which reflects the concentration of the 
bacterial culture.
Stationary OD600 Reading (stationary_od)
Steady State OD600 before declining in micrograms/mL
Sample Volume (sample_volume)
The volume of the bacterial sample taken for the experiment, should be 10% of the
Total volume (diluted volume and sample volume).
Dilute Volume (dilute_volume)
The volume of diluent added to the sample to achieve the desired dilution.
Total Volume (total_volume)
The total volume of the culture, including the sample and diluent.
Target Time (target_time)
The desired time to achieve a specific colony-forming unit (CFU) size.
Lagged Time (lagged_time)
The time it takes for the bacteria to stabilize in their environment before 
exponential growth begins.
Double Time (double_time or cal_double_time)
The time required for the bacterial population to double in size.
Final Time (final_time)
The end time for the simulation, which is the total duration of the experiment.
Polynomial Function Validity (poly_fun_is_valid)
Determines if the double time can be calculated automatically ('yes') or needs 
to be input manually ('no').
These inputs enable users to customize the bacterial growth simulation to match their experimental conditions or use collected data for more accurate modeling.


@author: andrewchadwell
"""
from bacterial_growth_copy import *

def main():
    '''
    
    Call bacterial_growth simulation script
    
    '''
    # Input data
    test = eppendorf_data # dictionary of the data transfer from ecoli_data.py
    #test = None # dictionary of the data transfer from ecoli_data.py
    
    # Input Parameters
    sample_volume= 200 # Volume of the sample in microliters
    dilute_volume= 1800 # Volume of diluent added in microliters
    target_time = 120 # Time for the desired CFU size in minutes
    final_time = 1000 # Time for the desired end time in minutes 
    cell_mass = 1e-6 # mass of bacteria in micrograms
    cell_per_volume = 1e9 # liquid concentration in cell/mL
    
    # ---------- Option 1 ----------------------------------------------------
    # Input for user inputs if no data is provided
    if not test:    
        # Data attained from inputs in the Eppendorf BioSpectrometer kinetic 
        initial_read =  0.5 # Initial average OD600 reading in micrograms/mL
        stationary_od = 4.3 # Steady State OD600 before declining in micrograms/mL
        total_volume = 10 # Total volume from liquid sample in micrograms
        lagged_time = 120 # Time for bacteria to stabilize to its environment at the beginning in minutes
        double_time = 40 # Time for bacterial population to double in minutes
    # ---------- Option 2 ----------------------------------------------------
    # Input for user inputs if data is provided
    else: 
        # enter parameters with data:
        # Data attained from inputs in the Eppendorf BioSpectrometer kinetic 
        initial_read = average_initial_od(data=test) # Initial average OD600 reading in micrograms/mL
        stationary_od = calculate_stationary_od(data=test) # Steady State OD600 before declining in micrograms/mL 
        # poly_fun_is_valid = 'yes': when the double time can be calculated
        # poly_fun_is_valid = 'no': when the double time will be enter as a 
        # parameter by the user.
        poly_fun_is_valid = 'no'

    #-----------------Calculations---------------------------------------------
    if not test:                                            
        # E.coli test
        volume_info = VolumeInformation(total_volume=total_volume, sample_volume=sample_volume,
                                dilute_volume=dilute_volume)
        time_info = TimeInformation(lagged_time=lagged_time,target_time=target_time,
                            double_time=double_time, final_time=final_time)
        cell_info = CellInformation(cell_mass=cell_mass, cell_per_volume=cell_per_volume)
        config = GrowthConfiguration(volume_information=volume_info,
                                    time_information=time_info,
                                    cell_information=cell_info)
        # This subplots the concentration micrograms/mL vs time in minututes and Population
        # size in colony-forming units vs time in minutes using the model equations.
        # A second is plotted using collected data and the model equations. It returns printed
        # inputs and outputs.
        ecoli_read, original_cfu = BacterialGrowth(initial_read=initial_read,
           data=test,stationary_od=stationary_od,
           configuration=config).bacterial_growth()
        # If data is provided, the plot_data function in Plot class returns a plot of
        # all the data runs with their average line.
        Plot(data=test).plot_data(bacterial_title='Bacterial Growth')
        
    else:
        total_volume = test['total_volume'][0]*1000 # total volume from the sample
        initial_read = average_initial_od(data=test) # average initial read 
        stationary_od = calculate_stationary_od(data=test)
        
        # If data is provided, the plot_data function in Plot class returns a plot of
        # all the data runs with their average line.
        Plot(data=test, final_time=final_time).plot_data(bacterial_title='Bacterial Growth')
        while True:
            # with collected data use calculated_double_time function in class Plot. It will
            # give a plotted line that correspondes with the data by changing the initial
            # initial read and degree of a polynomial.
            if poly_fun_is_valid == 'yes':
                # Initial degree of a polynomial
                lagged_time = int(input(f"Input the lagged time in minutes: ")) # Time in minutes for bacteria to stabilize to its environment at the beginning
                cal_double_time = Plot(data=test, initial_read=initial_read).calculate_double_time()
            else:
                lagged_time = int(input(f"Input the lagged time in minutes: ")) # Time in minutes for bacteria to stabilize to its environment at the beginning
                cal_double_time = int(input(f"Input the double time in minutes: ")) # Time in minutes when the bacteria is doubling in size.
               
            volume_info = VolumeInformation(total_volume=total_volume, sample_volume=sample_volume,
                                    dilute_volume=dilute_volume)
            time_info = TimeInformation(lagged_time=lagged_time,target_time=target_time,
                                double_time=cal_double_time, final_time=final_time)
            cell_info = CellInformation(cell_mass=1e-6, cell_per_volume=1e9)
            config = GrowthConfiguration(volume_information=volume_info,
                                        time_information=time_info,
                                        cell_information=cell_info)
            
            # This subplots the concentration micrograms/mL vs time in minututes and Population
            # size in colony-forming units vs time in minutes using the model equations.
            # A second is plotted using collected data and the model equations. It returns printed
            # inputs and outputs.
            ecoli_read, original_cfu = BacterialGrowth(initial_read=initial_read,
               stationary_od=stationary_od, data=test,configuration=config).bacterial_growth()
            # Ask user if the current double time is suitable
            user_input = input("Is the current double time suitable? (yes/no): ")
            if user_input.lower() == 'yes':
                break  # Exit loop if suitable

if __name__ == "__main__":
    main()