#!/usr/bin/env python3
'''
Bacterial Growth Simulation Script

This script simulates bacterial growth using parameters for volume, time,
and cell information.
It includes functionality to plot the growth curves and validate inputs to
ensure they are within acceptable ranges.
Do list:
    1. collected growth data
    2. lagged time, death rate, death time

Classes:
    VolumeInformation: Holds volume-related parameters.
    TimeInformation: Holds time-related parameters.
    CellInformation: Holds cell-related parameters.
    Plot: Handles plotting of the simulation results.
    Print: Prints inputs and the results
    BacterialGrowth: Performs the bacterial growth simulation and analysis.

Example:
    volume_info = VolumeInformation(total_volume=25, sample_volume=200, dilute_volume=1800)
    time_info = TimeInformation(lagged_time=2*60, target_time=60*7,
                                double_time=60*2.5, final_time=60*10)
    cell_info = CellInformation(cell_mass=47.65e-6, cell_per_volume=1E7)

    bacterial_growth_simulation = BacterialGrowth(
        initial_read=0.48,
        volume_info=volume_info,
        time_info=time_info,
        cell_info=cell_info
    )
    array_1, array_2 = bacterial_growth_simulation.bacterial_growth()

Andrew Chadwell June 28, 2024
'''
# -*- coding: utf-8 -*-

from numpy import math, zeros, log, exp
from pylab import plot,xlabel,ylabel,yscale,grid,title,show,figure,subplot

class VolumeInformation:
    '''
    Class to hold volume information for bacterial growth simulation.
    '''
    def __init__(self,total_volume=25, sample_volume=200, dilute_volume=1800):
        '''
        Parameters
        ----------
        total_volume : int
            The default is 25.
        sample_volume : int
            The default is 200.
        dilute_volume : int
            The default is 1800.
        Returns
        -------
        None.

        '''
        self.total_volume = total_volume
        self.sample_volume = sample_volume
        self.dilute_volume = dilute_volume
    def get_volumes(self):
        '''
        Returns volumes variables
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
                 final_time=400):
        self.target_time = target_time
        self.lagged_time = lagged_time
        self.double_time = double_time
        self.final_time = final_time
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
    
        Attributes
        ----------
        target_time : int or float
            Stores the target time in minutes.
        lagged_time : int or float
            Stores the lagged time in minutes.
        double_time : int or float
            Stores the doubling time in minutes.
        final_time : int or float
            Stores the final time in minutes.
        '''
class CellInformation:
    '''
    Class to hold cell information for bacterial growth simulation.
    '''
    def __init__(self,cell_mass=1e-6, cell_per_volume=1E9):
        self.cell_mass = cell_mass
        self.cell_per_volume = cell_per_volume
class Plot:
    '''
    Class to handle plotting and linear interpolation of bacterial growth
    simulation results.
    '''
    def plot_result(self, target_time,read,dilute_conversion,time):
        '''
        Plot and calculation linear interpolations the results of bacterial
        growth simulation
        '''
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
        plot(time, read, '--')
        plot(target_time, target_amount,'*r')
        title('Eppendorf BioSpectrometer kinetic read')
        xlabel('Time, Minutes')
        ylabel('Population Size, micrograms/mL')
        grid(visible=True)
        subplot(1,2,2)
        plot(time, dilute_conversion, '--')
        plot(target_time, target_original_count,'*r')
        title('From Original Sample vs Time')
        yscale('log')
        xlabel('Time, Minutes')
        ylabel('Population Size, in Colony Forming Units')
        grid(visible=True)
        show()
        return target_amount, target_original_count

class Print:
    '''
    Print information on inputs and target amounts
    '''
    def print_information(self, print_variables):
        '''
        Parameters
        ----------
        print_variables : dict
            A dictionary containing the following keys:
    - initial_read : float
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
        print(f'''INPUTS:\n-----------
initial read: {print_variables['initial_read']} micrograms/mL
original sample volume: {print_variables['total_volume']} mL
lagged time: {print_variables['lagged_time']} minutes
double time: {print_variables['double_time']} minutes
cell mass: {print_variables['cell_mass']} micrograms
cell count per volume: {print_variables['cell_per_volume']:.2e} CFU/mL
sample volume: {print_variables['sample_volume']} mL
dilute volume: {print_variables['dilute_volume']} mL

OUTPUTS at target time {print_variables['target_time']} minutes
------------------------------------------
target cell count read: {print_variables['target_amount']:.4f} micrograms/mL
original count at {print_variables['total_volume']} mL: {print_variables['target_original_count']:.2e} CFU
''')
class BacterialGrowth:
    '''
    Class to simulate bacterial growth using given parameters.
    '''
    def __init__(self,initial_read=1,volume_information=None,
                   time_information=None, cell_information=None):
        self.initial_read = initial_read
        self.volume_information = volume_information or VolumeInformation()
        self.time_information = time_information or TimeInformation()
        self.cell_information = cell_information or CellInformation()
        self.validate_inputs()

    def validate_inputs(self):
        '''
        Validate inputs to ensure they are within acceptable ranges or types.
        '''
        if (not isinstance(self.initial_read, (int, float))
            or self.initial_read < 0):
            raise ValueError("Initial read must be a non-negative number.")

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
        and expontial growth restricted by the original sample volume.

        Inputs:
        initial_read: Initial read from spectrometer in micrograms/mL
        total_volume: Total volume of sample in mL
        lagged_time: Time in minutes for the bacterial population to stay constant
        target_time: Time in minutes for the desired CFU count
        double_time: Time in minutes for bacterial population to double
        final_time: Final simulation time in minutes
        cell_mass: Mass of a single bacterial cell in micrograms
        cell_per_volume: Number of cells per mL when saturated
        sample_volume: Volume of the sample in microliters
        dilute_volume: Volume of diluent added in microliters

        Calculations:
        dilute_conversion: Corrected concentration considering dilution factor
        (unitless)
        time: Time intervals in minutes
        VOLUME: Total volume of sample and diluent in mL
        rate: Growth rate in 1/minute
        max_count: Maximum number of cells per volume in micrograms/mL
        read: Spectrometer read in micrograms/mL

        Outputs:
        - Prints all input parameters
        - Print target time with the predicted reading and colony forming unit
        using linear interpolation
        - Plots population size in from original sample
            in Colony Forming Units vs time in minutes
        - Plots population size in diluted sample size
            in micrograms/mL vs time in minutes
        - linear x-axis and log y-axis
        - Returns read array and dilute_conversion array

        Google Search:
        - Information about Escherichia coli and Saccharomyces yeast

        - Escherichia coli is a typical gram-negative rod bacterium. Its dimensions
        are those of a cylinder 1.0-2.0 micrometers long, with radius about
        0.5 micrometers. They double in population every 20 minutes. Their mass is
        1e-6 micr_g. They become saturated per volume for 1e9 cell/mL.

        - Yeast cell has an average diameter between 3 and 4 micrometers (μm),The largest
        yeast cells can be as big as 40 μm. Yeast cells can also be oval, spherical,
        elongated, or rectangular. They double in population every 90 minutes.
        Saccharomyces yeast cell is approximately 47.65 ± 1.05 pg. They become saturated
        per volume for (1e7 to 1e8) cells/mL.
        """
        initial_read = self.initial_read
        total_volume = self.volume_information.total_volume  # Accessing through VolumeInformation
        target_time = self.time_information.target_time
        lagged_time = self.time_information.lagged_time
        double_time = self.time_information.double_time
        final_time = self.time_information.final_time
        cell_mass = self.cell_information.cell_mass
        cell_per_volume = self.cell_information.cell_per_volume
        sample_volume = self.volume_information.sample_volume
        dilute_volume = self.volume_information.dilute_volume

        dilute_conversion = zeros(final_time+1) #colony forming unit
        read = zeros(final_time+1) #colony forming unit
        time = zeros(final_time+1) #MINUTE
        rate = log(2)/double_time #1/minutes
        max_count = total_volume*cell_per_volume*cell_mass  #micro_g/mL the maximum cell per volume
        dilute_conversion[0:lagged_time] = (
            initial_read*((sample_volume + dilute_volume)/sample_volume)
            *total_volume/cell_mass #colony forming unit
            )
        read[0:lagged_time] = initial_read

        for index in range(lagged_time, final_time+1, 1):
            read[index] = (
                max_count / (1 + (max_count - initial_read)/initial_read
                             *exp(-rate * (index-lagged_time+1))) #micr_g/mL
                )
            dilute_conversion[index] = (read[index]*((sample_volume + dilute_volume)
                            /sample_volume)*total_volume/cell_mass #unitless
                                    )
            time[index] = index #minutes

        target_amount, target_original_count = (
        Plot.plot_result(self, target_time, read, dilute_conversion, time)
        )

        print_variables = {
            'initial_read': initial_read,
            'total_volume': total_volume,
            'lagged_time': lagged_time,
            'double_time': double_time,
            'cell_mass': cell_mass,
            'cell_per_volume': cell_per_volume,
            'sample_volume': sample_volume,
            'dilute_volume': dilute_volume,
            'target_time': target_time,
            'target_amount': target_amount,
            'target_original_count': target_original_count
        }
        Print.print_information(self, print_variables)

        return read, dilute_conversion
 
    
# ECOLI = bacterial_growth(initial_read=1, total_volume=1000, lagged_time=2*60,
   #                         target_time=40.0, double_time=20,
   #                       final_time=400, cell_mass=1e-6, cell_per_volume=1E9,
   #                       sample_volume=200, dilute_volume=1800)
# volume_info = VolumeInformation(total_volume=25, sample_volume=200,
#                                 dilute_volume=1800)
# time_info = TimeInformation(lagged_time=2*60,target_time=60*7,
#                             double_time=60*2.5, final_time=60*10)
# cell_info = CellInformation(cell_mass=47.65e-6, cell_per_volume=1E7)
# yeast_read, original_cfu = BacterialGrowth(initial_read=0.48,
#                     volume_information=volume_info,time_information=time_info,
#                     cell_information=cell_info).bacterial_growth()
