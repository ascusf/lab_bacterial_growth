#!/usr/bin/env python3
'''
Bacterial Growth Simulation Script

This script simulates bacterial growth using parameters for volume, time,
and cell information.
It includes functionality to plot the growth curves and validate inputs to
ensure they are within acceptable ranges.

Do list:
    1. death rate, death time
    2. How does different Broth impact growth rate?
    3. create a separate script for different bacteria
    4. create a separate script for different media

Classes:
    VolumeInformation: Holds volume-related parameters.
    TimeInformation: Holds time-related parameters.
    CellInformation: Holds cell-related parameters.
    Plot: Handles plotting of the simulation results.
    Print: Prints inputs and the results
    BacterialGrowth: Performs the bacterial growth simulation and analysis.

Inputs:
    initial_read: Initial read from spectrometer in micrograms/mL
    statuionary_od: Maximum reading from spectrometer in micrograms/mL 
    data: data collected using BioSpectrometer kinetic instrument
    total_volume: Total volume of sample in mL
    sample_volume: Volume of the sample in microliters
    dilute_volume: Volume of diluent added in microliters
    lagged_time: Time in minutes for bacteria to stabilize to its environment at the beginning
    target_time: Time in minutes for the desired CFU size
    double_time: Time in minutes for bacterial population to double
    final_time: Final simulation time in minutes
    cell_mass: Mass of a single bacterial cell in micrograms
    cell_per_volume: Number of cells per mL when saturated

Calculations:
    - The target concentration in micrograms/mL corresponding with target time.
    - The target colony-forming units corresponding with target time.
    - The average data is calculated and plotted.
    - The calculation accounts for the lagged time and exponential growth of
        bacteria restricted by the original total volume.
    - The estimated double time for bacteria to double its size in minutes.
    - dilute_conversion: Corrected concentration considering dilution factor (unitless)
    - time: Time intervals in minutes
    - total_volume: Total volume of sample and diluent in mL
    - rate: Growth rate in 1/minute
    - read: Spectrometer read in micrograms/mL
    - 1 OD600 = saturation
    - 1 OD600 = 10 microgram/mL

Outputs:
    - Prints all input parameters
    - Print target time with the predicted reading in micrograms/mL
        and colony-forming unit
    - Plots population size scaling up to the original total volume the sample was
        collected from in Colony Forming Units vs time in minutes
    - Plots population size in diluted sample size
        in micrograms/mL vs time in minutes
    - Using linear interpolation plots, it plots the desired target time and
        target amount on both plots.
    - Plots data if given which can plotted with the model or alone.
    - linear x-axis and linear y-axis
    - Returns read array and dilute_conversion array

Example:


Google Search:
- Information about Escherichia coli and Saccharomyces yeast

- Escherichia coli is a typical gram-negative rod bacterium. Its dimensions
are those of a cylinder 1.0-2.0 micrometers long, with radius about
0.5 micrometers. They double in population every 20 minutes. Their mass is
1e-6 micro_grams. They become saturated per volume for 1e9 cell/mL.
Cell dry weight	3 x 10-13 g

- Yeast cell has an average diameter between 3 and 4 micrometers (μm),The largest
yeast cells can be as big as 40 μm. Yeast cells can also be oval, spherical,
elongated, or rectangular. They double in population every 90 minutes.
Saccharomyces yeast cell is approximately 47.65 ± 1.05 pg. They become saturated
per volume for (1e7 to 1e8) cells/mL.

Andrew June 28, 2024
'''

# -*- coding: utf-8 -*-
import copy, math
from pylab import plot,xlabel,ylabel,yscale,grid,title,show,figure,subplot,legend,xticks
from numpy import mean,array,arange,zeros,argmax,exp,poly1d,polyfit,diff,log,std
from scipy.stats import linregress
from ecoli_data import *

def average_initial_od(data=None):
    '''
    
    Calculate the average initial optical density (OD) from the provided data.

    Parameters
    ----------
    data : dict
        Data collected using BioSpectrometer kinetic instrument. 
        Should contain 'minute' and other keys with lists of measurements.

    Returns
    -------
    float
        The average initial OD if data is available, otherwise None.
    
    '''
    
    if not data:
        print("No data available to plot.")
        return
   
    # Extract Data
    minutes = copy.copy(data['minute'])
    data_vectors = {key: values for key, values in data.items()
                    if key not in ('minute', 'total_volume')}

    # Calculate the initial OD average point.
    for i in range(len(minutes)):
        data_points = [data[i] for data in data_vectors.values() if data[i] is not None]
        if data_points:
            return mean(data_points)

def calculate_lagged():
    '''
    
    ??????

    Returns
    -------
    None.

    '''
    
    return           

def linear_growth(data=None,time=None):
    '''

    Estimates the linear growth of bacterial culture based on OD600 measurements.

    This function calculates the growth rate by fitting a line to the points with
    the steepest slope in the OD600 data. It returns the estimated OD600 value at 
    a given time using the fitted line equation.

    Parameters
    ----------
    data : dict, optional
        A dictionary containing the experimental data. It should include a key 'minute'
        that corresponds to the time points and other keys for OD600 values. The default 
        is None.
    time : float, optional
        The time in minutes for which the OD600 value is to be estimated based on 
        the linear growth model. The default is None.

    Returns
    -------
    float
        The estimated OD600 value at the specified time.

    Notes
    -----
    - The function first calculates the average OD600 values and filters out any None 
      values.
    - It computes the slopes between consecutive time points and identifies the maximum 
      slope to determine the linear growth phase.
    - A linear fit (y = mx + b) is performed using the points around the maximum slope, 
      and the fitted line equation is used to estimate the OD600 value at the specified 
      time.

    '''
    
    # Retrieving data
    copy_minute = copy.copy(data['minute'])
    average = Plot(data=data).calculate_data_average()

    # Filter out None values
    find_none = [i for i, item in enumerate(average) if item is None]
    average = [item for i, item in enumerate(average) if i not in find_none]
    copy_minute = [item for item in copy_minute if item not in find_none]

    # Calculate the slopes between each pair of points
    differences = diff(array(average)) / diff(array(copy_minute))

    # Find the index of the maximum slope
    max_index = argmax(differences)

    # Use the points around the maximum slope to fit a line
    coefficients = polyfit(array(copy_minute)[max_index:max_index + 2],
                           array(average)[max_index:max_index + 2], 1)

    # Calculate the time to double the OD600
    polynomial = poly1d(coefficients)
    
    return polynomial[1] * time + polynomial[0]  

def poly3_growth(data=None,time=None):
    '''
    
    Estimates the polynomial degree of 3 growth of bacterial culture based on OD600 measurements.
    
    '''
    # Retrieving data
    copy_minute = copy.copy(data['minute'])
    average = Plot(data=data).calculate_data_average()

    # Filter out None values
    find_none = [i for i, item in enumerate(average) if item is None]
    average = [item for i, item in enumerate(average) if i not in find_none]
    copy_minute = [item for item in copy_minute if item not in find_none]

    # Calculate the slopes between each pair of points
    differences = diff(array(average)) / diff(array(copy_minute))

    # Find the index of the maximum slope
    max_index = argmax(differences)

    # Use the points around the maximum slope to fit a line
    coefficients = polyfit(array(copy_minute)[max_index:max_index + 2],
                           array(average)[max_index:max_index + 2], 3)

    # Calculate the time to double the OD600
    polynomial = poly1d(coefficients)
    
    return polynomial[1] * time + polynomial[0]  

def calculate_stationary_od(data=None):
    '''
    
    Broth determines the stationary phase?  
    Calculate the average steady OD from the provided data.

    Parameters
    ----------
    data : dict
        Data collected using BioSpectrometer kinetic instrument. 
        Should contain 'minute' and other keys with lists of measurements.

    Returns
    -------
    float
        The average steady state OD if data is available, otherwise None.
    
    '''
    
    # Return if no data is provided
    if not data:
        print("No data available to plot.")
        return
    
    # Extract Data
    data_copy =  copy.copy(data)
    last_od = {key: values[-1] for key, values in data_copy.items()
                   if key not in ('minute', 'total_volume')}
    
    # Calculate the final OD average point.
    last_od_vector = array(list(last_od.values()))
    
    return mean(last_od_vector)

class VolumeInformation:
    '''
    
    Class to hold volume information for bacterial growth simulation.
    
    '''
    
    def __init__(self,total_volume=25, sample_volume=200, dilute_volume=1800):
        '''
        
        Initialize the VolumeInformation class with the provided volumes.

        Parameters
        ----------
        total_volume : int
            Total volume for the simulation (default is 25).
        sample_volume : int
            Sample volume for the simulation (default is 200).
        dilute_volume : int
            Dilute volume for the simulation (default is 1800).
        Returns
        -------
        None
        
        '''
        
        self.total_volume = total_volume
        self.sample_volume = sample_volume
        self.dilute_volume = dilute_volume

    def get_volumes(self):
        '''
       
        Get the volume information.

       Returns
       -------
       dict
           A dictionary containing the total volume, sample volume, and dilute volume.
       
        '''
        
        return {
    'total_volume': self.total_volume,
    'sample_volume': self.sample_volume,
    'dilute_volume': self.dilute_volume
                }

class TimeInformation:
    '''
    
    Class to hold time information for bacterial growth simulation.
    
    '''
    
    def __init__(self,target_time=20,lagged_time=2*60, double_time=20,
                 stationary_time=420, final_time=400):
        '''
        
        Initialize the TimeInformation class with parameters for bacterial growth simulation.

        Parameters
        ----------
        target_time : int or float, optional
            The time in minutes at which the target cell count is desired. Default is 20 minutes.
        lagged_time : int or float, optional
            The time in minutes during which the bacterial population remains
            constant before starting exponential growth. Default is 120 minutes (2 hours).
        double_time : int or float, optional
            The time in minutes required for the bacterial population to double
            in size. Default is 20 minutes.
        final_time : int or float, optional
            The final time in minutes up to which the simulation is run. Default is 400 minutes.

        Returns
        -------
        None
        
        '''

        self.target_time = target_time
        self.lagged_time = lagged_time
        self.double_time = double_time
        self.stationary_time = stationary_time
        self.final_time = final_time

    def get_time(self):
        '''
       
        Get the time information.

       Returns
       -------
       dict
           A dictionary containing the target_time, lagged_time, double_time, final_time.
       
        '''
        
        return {
    'target_time': self.target_time,
    'lagged_time': self.lagged_time,
    'double_time': self.double_time,
    'stationary_time': self.stationary_time,
    'final_time': self.final_time
                }

class CellInformation:
    '''
    
    Class to hold cell information for bacterial growth simulation.
    
    '''
    
    def __init__(self,cell_mass=1e-6, cell_per_volume=1E9):
        '''
        
        Initialize the CellInformation class with the provided cell mass and cell count per volume.

        Parameters
        ----------
        cell_mass : float
            Mass of a single cell in grams (default is 1e-6).
        cell_per_volume : float
            Number of cells per unit volume (default is 1E9).

        Returns
        -------
        None
        
        '''
        
        self.cell_mass = cell_mass
        self.cell_per_volume = cell_per_volume

    def get_cell(self):
        '''
       Get the cell information.

       Returns
       -------
       dict
           A dictionary containing the cell_mass and cell_per_volume.
       '''
        return {
    'cell_mass': self.cell_mass,
    'cell_per_volume': self.cell_per_volume
                }

class Plot:
    '''
    Class to handle plotting and linear interpolation of bacterial growth
    simulation results.
    '''
    def __init__(self,data=None,target_time=None, initial_read=None, final_time=None):
        '''
        Initialize the Plot class with the provided parameters.

        Parameters
        ----------
        data : dict
        Data for plotting bacterial growth. Default is None.
        target_time : int
        The time in minutes at which the target cell count is desired. Default None
        initial_read: int
        Initial read from spectrometer in micrograms/mL. Default None

        Returns
        -------
        None
        '''
        self.data = data
        self.target_time = target_time
        self.initial_read = initial_read
        self.final_time = final_time ### ADD

    def validate_inputs(self):
        '''
        Validate inputs to ensure they are within acceptable ranges or types.

        Returns
        -------
        None
        '''
        if (not isinstance(self.target_time, (int, float))
            or self.target_time <= 0):
            raise ValueError("Target time must be a positive number.")

        if "minute" not in self.data:
            raise ValueError("Data dictionary must contain the key 'minute'.")

        minutes = self.data["minute"]

        if not all(isinstance(m, (int, float)) and m >= 0 for m in minutes):
            raise ValueError("'minute' list must contain non-negative numbers.")

        data_length = len(minutes)

        for key, values in self.data.items():
            if key != "minute":
                if not isinstance(values, list):
                    raise ValueError(f"""All values in the data
                  dictionary must be lists. Error with key: {key}""")
                if len(values) != data_length:
                    raise ValueError(f"""All lists in the data
                dictionary must be of the same length as 'minute'
                list. Error with key: {key}""")
                if not all(isinstance(v, (int, float, type(None))) for v in values):
                    raise ValueError(f"""All elements in the data lists
                    must be numbers or None. Error with key: {key}""")

        if (not isinstance(self.initial_read, (int, float))
            or self.initial_read < 0):
            raise ValueError("Initial read must be a non-negative number.")

    def plot_result(self,read,dilute_conversion,time):
        '''
        
        Plot and calculation linear interpolations the results of bacterial
        growth simulation

        Parameters
        ----------
        read : Array of float64
            Reading from Eppendorf BioSpectrometer kinetic in micrograms/mL with respect to time.
        dilute_conversion : Array of float64
            Converted read array [micrograms/mL] to colony forming unit
        time : Array of float64
            Time intervals in minutes starting with 0 to the final time.

        Returns
        -------
        target_amount : float64
            The target concentration in micrograms/mL corresponding with target time.
        target_original_count : float64
            The target colony forming units corresponding with target time.
            
        '''
        
        target_time = self.target_time
        if target_time % 1 == 0:
            target_amount = read[int(target_time)]
            target_original_count = dilute_conversion[int(target_time)]
        else:
            x_1= int(math.ceil(target_time))
            x_0 = int(math.floor(target_time))
            y_1 = read[x_1]
            y_0 = read[x_0]
            slope = (y_1 - y_0)/(x_1- x_0)
            y_int = y_0 - slope*x_0
            target_amount = slope*target_time + y_int
            y_1 = dilute_conversion[x_1]
            y_0 = dilute_conversion[x_0]
            slope = (y_1 - y_0)/(x_1- x_0)
            y_int = y_0 - slope*x_0
            target_original_count = slope*target_time + y_int

        # Subplots
        figure(figsize=(10,5))
        subplot(1,2,1)
        plot(time, read, linestyle='-',color='k')
        plot(target_time, target_amount,'*k')
        title('Eppendorf BioSpectrometer kinetic read')
        xlabel('Time, Minutes')
        ylabel('Concentration, micrograms/mL')
        grid(visible=True)
        subplot(1,2,2)
        plot(time, dilute_conversion, linestyle='-',color='k')
        plot(target_time, target_original_count,'*k')
        title('From Original Sample vs Time')
        yscale('log')
        xlabel('Time, Minutes')
        ylabel('Population Size, in Colony Forming Units')
        grid(visible=True)
        show()
        return target_amount, target_original_count

    def plot_data(self, bacterial_title='Bacterial Growth'):
        '''
        
        Plot the bacterial growth data representing symbols with all runs
        and calculates and plots the averages representing a line. The x-axis
        represents time (in hours), and the y-axis represents bacterial growth.

        Parameters
        ----------
        self.data : dict
            Data dictionary containing time points and corresponding bacterial
            growth measurements.
        bacterial_title : str, optional
            Title for the bacterial growth plot. Default is 'Bacterial Growth'.

        Returns
        -------
        None
        
        '''
        
        data = self.data
        # extract data
        if not data:
            print("No data available to plot.")
            return
        minutes = data['minute']
        data_vectors = {key: values for key, values in data.items()
                        if key not in ('minute', 'total_volume')}

        # calculate the data average points.
        average = []
        for i in range(len(minutes)):
            data_points = [data[i] for data in data_vectors.values() if data[i] is not None]
            if data_points:
                average.append(mean(data_points))
            else:
                average.append(None)

        # Plot the data
        markers = ['.','o','*','+','x','>','<']  # Different markers for different datasets
        for i, (key, values) in enumerate(data_vectors.items()):
            plot(minutes, values, label=f'Data {key}',
                 marker=markers[i % len(markers)], linestyle='', color='k')

        plot(minutes, average, label='Average', linestyle='-', color='k')
        xlabel('Time (minutes)')
        ylabel('OD_600 (Micro_g/mL)')
        title(bacterial_title)
        legend()
        grid(True)
        xticks(arange(0, math.ceil(minutes[-1]+(minutes[-1]-minutes[-2])),20), rotation='vertical')
        show()
        return

    def calculate_data_average(self):
        '''
        Calculates the average of a set of data runs
       
        Returns
        -------
        average : list
            Average points of data runs in micrograms/mL.

        '''
        data = self.data
        # extract data
        if not data:
            return print("No data available to plot.")
        minutes = data['minute']
        data_vectors = {key: values for key, values in data.items()
                        if key not in ('minute','total_volume')}
        average = []
        # calculate averages
        for i in range(len(minutes)):
            data_points = [data[i] for data in data_vectors.values() if data[i] is not None]
            if data_points:
                average.append(mean(data_points))
            else:
                average.append(None)
        return average

    def plot_data_model(self,read,time):
        '''
        
        Plot the bacterial growth data average points and the Model.

        Parameters
        ----------
        data : dict
            Data dictionary containing time points and corresponding bacterial growth measurements.
        bacterial_title : str, optional
            Title for the bacterial growth plot. Default is 'Bacterial Growth'.

        Returns
        -------
        None
        
        '''
        
        # Retrieving data
        data = self.data
        
        # Subplots model
        figure(figsize=(10,5))
        subplot(1,2,1)
        plot(time, read,linestyle='-',color='k',label='Model')
        title('Eppendorf BioSpectrometer kinetic read')
        xlabel('Time, Minutes')
        ylabel('Concentration, micrograms/mL')
        grid(visible=True)

        # extract data
        if not data:
            return print("No data available to plot.")
        minutes = data['minute']

        # calculate the data average points.
        average = Plot.calculate_data_average(self)

        # Plot the data
        subplot(1,2,1)
        plot(minutes, average, label='Data Average',linestyle='', marker='.', color='k')
        xlabel('Time (minutes)')
        ylabel('OD_600 (Micro_g/mL)')
        legend()
        grid(True)
        show()
        
        return print('Data available to plot')

    def calculate_double_time(self):
        '''
        
        Estimates the bacterial doubling time based on OD600 measurements. For
        instance, it calculates the time it takes for the OD600 to increase from 1 to 2.
        The initial OD600 times by 2 and returns the time it takes
        the bacteria to grow twice itself using polyval. The data is collected.
        The code calculates the slopes between each point, finds the index of the 
        maximum slope and uses linear equation to estimate the doubling time.
        Parameters
        ----------
        None

        Returns
        -------
        float64
            The time in minutes it takes the bacteria to double in size.
        
        '''
        # Retrieving data
        copy_minute = copy.copy(self.data['minute'])
        average = Plot.calculate_data_average(self)
        
        # Filter out None values
        find_none = [i for i, item in enumerate(average) if item is None]
        average = [item for i, item in enumerate(average) if i not in find_none]
        copy_minute = [item for item in copy_minute if item not in find_none]
        
        x_axis = array(copy_minute)
        y_axis = array(average)
        
        # Calculate the slopes between each pair of points
        differences = diff(y_axis) / diff(x_axis)
        
        # Calculate mean and standard deviation of the growth rates
        mean_growth = mean(differences)
        std_growth = std(differences)
 
        # find growth phase of the data
        xx = []
        yy = []
        for index, ele in enumerate(differences):
            if ele > mean_growth + 0.05 * std_growth: # Define the threshold value (mean + n * std)
                xx.append(x_axis[index+1])
                yy.append(y_axis[index+1])
        
        # Ensure there are enough points to perform regression
        if len(xx) < 2:
            print("Insufficient data points in the exponential phase.")
            return None    
        
        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = linregress(xx, yy)
        
        return log(2)/slope

class GrowthConfiguration:
    '''
    
    Class to hold all configuration information for bacterial growth simulation.

        Parameters
        ----------
        volume_information : VolumeInformation
            An instance of the VolumeInformation class containing details about
            the volumes involved in the simulation nsuch as total_volume,
            sample_volume, and dilute_volume.
        time_information : TimeInformation
            An instance of the TimeInformation class containing details about
            the timing of the simulation, such as lag time, target time,
            doubling time, and final time.
        cell_information : CellInformation
            An instance of the CellInformation class containing details about
            the cells, including cell mass and cells per volume.

        Returns
        -------
        None
        
    '''
    
    def __init__(self, volume_information, time_information, cell_information):
        '''
        Initialize the Plot class with the provided parameters.

        Parameters
        ----------
        volume_information : VolumeInformation
            An instance of the VolumeInformation class containing details about
            the volumes involved in the simulation nsuch as total_volume,
            sample_volume, and dilute_volume.
        time_information : TimeInformation
            An instance of the TimeInformation class containing details about
            the timing of the simulation, such as lag time, target time,
            doubling time, and final time.
        cell_information : CellInformation
            An instance of the CellInformation class containing details about
            the cells, including cell mass and cells per volume.

        Returns
        -------
        None.

        '''
        
        self.volume_information = volume_information
        self.time_information = time_information
        self.cell_information = cell_information

class BacterialGrowth:
    '''
    
    Class to simulate bacterial growth using given parameters.
    
    '''
    
    def __init__(self,initial_read=1,stationary_od=5,data=None, configuration=None):
        '''
        
        Initialize the BacterialGrowth class with the provided parameters.

        Parameters
        ----------
        initial_read : int or float
        Initial read from spectrometer in micrograms/mL. Default 1
        data : dict
        data collected using BioSpectrometer kinetic instrument. Default is None.
        volume_information : VolumeInformation
        Instance of VolumeInformation class. Default is None.
        time_information : TimeInformation
        Instance of TimeInformation class. Default is None.
        cell_information : CellInformation
        Instance of CellInformation class. Default is None.

        Returns
        -------
        None
        
        '''
        
        self.initial_read = initial_read
        self.stationary_od = stationary_od
        self.data = data
        self.volume_information = configuration.volume_information
        self.time_information = configuration.time_information
        self.cell_information = configuration.cell_information
        self.validate_inputs()

    def calculate_growth_rate(self):
        '''
        
    Calculate the growth rate value based on the doubling time in minutes.

    The growth rate (k) is calculated using the formula:

        k = log(2) / doubling_time

    where:
    - log(2) is the natural logarithm of 2.
    - doubling_time is the time it takes for the bacterial population to double,
    in minutes.

    Parameters
    ----------
    None

    Returns
    -------
    float
        The growth rate in units of 1/minutes.
        
        '''
        
        return log(2)/self.time_information.double_time #1/minutes

    def calculation_measurements(self):
        '''

    Calculate the growth curve for bacterial culture based on given parameters.

    This method simulates the bacterial growth curve over time and returns the
    calculated values for dilute conversion, readings, and time points. The
    logistic equation was derived from 
    https://en.wikipedia.org/wiki/Logistic_function .

    The growth curve is modeled using the initial readings, dilution volumes,
    and other parameters provided during the initialization of the BacterialGrowth
    class instance.

    Parameters
    ----------
    None

    Returns
    -------
    dilute_conversion : Array of float64
        Array containing the calculated dilute conversion values over time.
        These values represent the colony-forming units adjusted for dilution.

    read : Array of float64
        Array containing the calculated optical density (OD) readings over time.
        These values represent the bacterial growth measured as OD_600 in
        micrograms/mL.

    time : Array of float64
        Array containing the time points (in minutes) at which the readings
        and dilute conversions were calculated.

        '''
        # Inputs
        lagged_time = self.time_information.lagged_time
        stationary_time =  self.time_information.stationary_time
        final_time = self.time_information.final_time
        initial_read = self.initial_read
        stationary_od = self.stationary_od
        sample_volume = self.volume_information.sample_volume
        dilute_volume = self.volume_information.dilute_volume
        total_volume = self.volume_information.total_volume
        cell_mass = self.cell_information.cell_mass

        # Prepare vectors for For Loop
        dilute_conversion = zeros(final_time+1) #colony forming unit
        read = zeros(final_time+1) #colony forming unit
        time = arange(0, final_time+1) #MINUTE
        
        # Lagged OD with respect with time.
        read[0:round(lagged_time)+1] = initial_read #micrograms/mL 
        dilute_conversion[0:round(lagged_time)+1] = (
            initial_read*((sample_volume + dilute_volume)/sample_volume)
            *total_volume/cell_mass 
            ) #colony forming unit
        
        # Growth OD with respect with time
        for index in range(1, final_time+1):
            if round(lagged_time) + index - 1 >= len(read):
                break
            
            # Logistic Equation
            read[round(lagged_time)+index-1] = (
                stationary_od / (1 + ((stationary_od - initial_read) / initial_read) 
                    * exp(-BacterialGrowth.calculate_growth_rate(self) * index)) 
                ) #micrograms/mL
            
            # Conversion from micrograms/mL to CFU
            dilute_conversion[round(lagged_time)+index-1] = (read[index]*((sample_volume + dilute_volume)
                            /sample_volume)*total_volume/cell_mass
                                    )  #unitless

        return dilute_conversion,read,time

    def print_information(self,target_amount,target_original_count):
        '''
        Print relative information
        
        Parameters
        ----------
        print_variables : dict
            A dictionary containing the following keys:
    - initial_read : float
    - stationary_od : float
    - total_volume : float
    - lagged_time : float
    - double_time : float
    - cell_mass : float
    - cell_per_volume : float
    - sample_volume : float
    - dilute_volume : float
    - target_time : float
    - target_amount : float
    - target_original_count : float

        Returns
        -------
        None.

        '''
        print_variables = {
            'initial_read': self.initial_read,
            'stationary_od': self.stationary_od,
            'total_volume': self.volume_information.total_volume,
            'lagged_time': self.time_information.lagged_time,
            'double_time': self.time_information.double_time,
            'cell_mass': self.cell_information.cell_mass,
            'cell_per_volume': self.cell_information.cell_per_volume,
            'sample_volume': self.volume_information.sample_volume,
            'dilute_volume': self.volume_information.dilute_volume,
            'target_time': self.time_information.target_time,
            'target_amount': target_amount,
            'target_original_count': target_original_count
        }
        print(f'''INPUTS:\n-----------
initial read: {print_variables['initial_read']:.2f} micrograms/mL
stationary OD: {print_variables['stationary_od']:.2f} micrograms/mL
original sample volume: {print_variables['total_volume']/1000} mL
lagged time: {print_variables['lagged_time']} minutes
double time: {print_variables['double_time']:.2f} minutes
cell mass: {print_variables['cell_mass']} micrograms
cell count per volume: {print_variables['cell_per_volume']:.2e} CFU/mL
sample volume: {print_variables['sample_volume']/1000} mL
dilute volume: {print_variables['dilute_volume']/1000} mL

OUTPUTS at target time {print_variables['target_time']} minutes
------------------------------------------
target cell concentration read: {print_variables['target_amount']:.4f} micrograms/mL
target cell size in colony-forming Units: {print_variables['target_original_count']:.2e} CFU
''')

    def validate_inputs(self):
        '''
        Validate inputs to ensure they are within acceptable ranges or types.

        Returns
        -------
        None
        '''
        if (not isinstance(self.initial_read, (int, float))
            or self.initial_read < 0):
            raise ValueError("Initial read must be a non-negative number.")
            
        if (not isinstance(self.stationary_od, (int, float))
            or self.stationary_od < 0):
            raise ValueError("Stationary OD must be a non-negative number.")    

        if (not isinstance(self.volume_information.total_volume, (int, float))
            or self.volume_information.total_volume <= 0):
            raise ValueError("Total volume must be a positive number.")

        if (not isinstance(self.volume_information.sample_volume, (int, float))
            or self.volume_information.sample_volume <= 0):
            raise ValueError("Sample volume must be a positive number.")

        if (not isinstance(self.volume_information.dilute_volume, (int, float))
            or self.volume_information.dilute_volume < 0):
            raise ValueError("Dilute volume must be a non-negative number.")

        if (not isinstance(self.time_information.target_time, (int, float))
            or self.time_information.target_time <= 0):
            raise ValueError("Target time must be a positive number.")

        if (not isinstance(self.time_information.lagged_time, (int, float))
            or self.time_information.lagged_time < 0):
            raise ValueError("Lagged time must be a non-negative number.")

        if (not isinstance(self.time_information.double_time, (int, float))
            or self.time_information.double_time <= 0):
            raise ValueError("Double time must be a positive number.")

        if (not isinstance(self.time_information.final_time, (int, float))
            or self.time_information.final_time <= 0):
            raise ValueError("Final time must be a positive number.")

        if (not isinstance(self.cell_information.cell_mass, (int, float))
            or self.cell_information.cell_mass <= 0):
            raise ValueError("Cell mass must be a positive number.")

        if (not isinstance(self.cell_information.cell_per_volume, (int, float))
            or self.cell_information.cell_per_volume <= 0):
            raise ValueError("Cell per volume must be a positive number.")

    def bacterial_growth(self):

        """
        Perform bacterial growth simulation and analysis using the Eppendorf
        BioSpectrometer kinetic. The calculation account for the lagged time
        and expontial growth restricted by the original total volume.

        Parameters
        ----------
        None

        Returns
        -------

        """
        
        # calculate the logistic equation
        dilute_conversion,read,time = BacterialGrowth.calculation_measurements(self)
        
        # Subplot
        target_amount, target_original_count = (
        Plot(target_time=self.time_information.target_time)
        .plot_result(read, dilute_conversion, time)
        )
        
        # Plot data with model
        Plot.plot_data_model(self,read,time)
        
        # print the information
        BacterialGrowth.print_information(self,target_amount,target_original_count)

        return read, dilute_conversion
    