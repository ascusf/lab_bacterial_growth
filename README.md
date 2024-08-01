# Bacterial-Growth 
# laboratories - Bacterial Growth Simulation 

## Overview 
  
This project contains a script for simulating bacterial growth using various parameters for volume, time, and cell information. It includes functionality to plot growth curves and validate inputs to ensure they are within acceptable ranges.

The bacterial_growth.py is a comprehensive bacterial growth simulation tool designed to model and analyze the growth of bacterial cultures, specifically using E. coli. The script includes classes to manage volume, time, and cell information, and handles plotting and calculations related to bacterial growth.
 
The call_bacterial_growth.py simulates bacterial growth using either predefined data or user input. It calculates and plots bacterial concentration and population size over time, with options for inputting experimental parameters directly or using provided data sets. 

The ecoli_data.py collected data that can be call into scripts to analyze.

I’m currently in the early stages of programming, and I have a lot of work ahead. My goal is to become familiar with Python and GitHub.
  

## Features 

- **Volume Information:** Manage and validate volume-related parameters. 

- **Time Information:** Manage and validate time-related parameters. 

- **Cell Information:** Manage and validate cell-related parameters. 

- **Plotting:** Visualize bacterial growth simulation results with plotting functionality. 

- **Validation:** Ensure inputs are valid and handle errors gracefully. 

- **Simulation:** Perform bacterial growth simulation and analyze results. 

## Installation 

To use this script, you need to have Python 3, packages numpy and pylab installed and imports:  

  
- numpy (import mean,array,arange,zeros,argmax,exp,poly1d,polyfit,diff,log)

- pylab (import plot,xlabel,ylabel,yscale,grid,title,show,figure,subplot,legend)

- ecoli_data (have all the data)

- import copy 

You can install numpy and pylab packages using pip: 

## Define volume information 

volume_info = VolumeInformation(total_volume=25, sample_volume=200, dilute_volume=1800) 

# Functions
init(self, total_volume, sample_volume, dilute_volume): Initializes the volume parameters.
get_total_volume(self): Returns the total volume of the sample.
get_sample_volume(self): Returns the volume of the sample.
get_dilute_volume(self): Returns the volume of the diluent.
  

## Define time information 

time_info = TimeInformation(lagged_time=2*60, target_time=60*7, double_time=60*2.5, final_time=60*10) 

# Functions
init(self, lagged_time, target_time, double_time, final_time): Initializes the time parameters.
get_lagged_time(self): Returns the lagged time.
get_target_time(self): Returns the target time.
get_double_time(self): Returns the doubling time.
get_final_time(self): Returns the final simulation time.
  

## Define cell information 

cell_info = CellInformation(cell_mass=47.65e-6, cell_per_volume=1E7) 

# Functions
init(self, cell_mass, cell_per_volume): Initializes the cell parameters.
get_cell_mass(self): Returns the mass of a single bacterial cell.
get_cell_per_volume(self): Returns the number of cells per mL when saturated.
  

## Define Plot
average = Plot.calculate_data_average(self)
Plot.plot_data_model(self,read,time)
Plot(data=test).plot_data(bacterial_title='Bacterial Growth')

# Functions
init(self, data=None, initial_read=None): Initializes the plot parameters.
plot_growth_curve(self, time, concentration, title): Plots the growth curve of the bacteria.
calculate_double_time(self, deg=3): Calculates and plots the double time for bacterial growth.
plot_data(self, bacterial_title): Plots the data runs with their average line.

## Define Print
print_inputs(self, initial_read, data, volume_info, time_info, cell_info): Prints all input parameters.
print_results(self, target_time, predicted_reading, cfu): Prints the target time with the predicted reading and colony-forming units.

## Define BacterialGrowth
init(self, initial_read, data, configuration): Initializes the bacterial growth simulation parameters.
bacterial_growth(self): Performs the bacterial growth simulation and analysis.


# Run the bacterial growth simulation and plots
- with collected data use calculated_double_time function in class Plot.
- The code calculates the slopes between each point, finds the index of the maximum slope and uses linear equation to estimate the doubling time.
cal_double_time = Plot(data=test,initial_read=initial_read).calculate_double_time()

volume_info = VolumeInformation(total_volume=total_volume, sample_volume=200,
                        dilute_volume=1800)
time_info = TimeInformation(lagged_time=0*3,target_time=200,
                    double_time=cal_double_time, final_time=1500)
cell_info = CellInformation(cell_mass=1e-6, cell_per_volume=1e9)
config = GrowthConfiguration(volume_information=volume_info,
                            time_information=time_info,
                            cell_information=cell_info)
- subplots the concentration micrograms/mL vs time in minututes and Population size in colony-forming units vs time in minutes using the model equations.
- A second plot is plotted using collected data and the model equations.
- It returns printed inputs and outputs.
ecoli_read, original_cfu = BacterialGrowth(initial_read=initial_read,
   data=test,configuration=config).bacterial_growth()
- If data is provided, the plot_data function in Plot class returns a plot of
- all the data runs with their average line.
Plot(data=test).plot_data(bacterial_title='Bacterial Growth')

## Classes and Methods 

VolumeInformation: Holds volume-related parameters.

TimeInformation: Holds time-related parameters.

CellInformation: Holds cell-related parameters.

Plot: Handles plotting of the simulation results.

Print: Prints inputs and the results.

BacterialGrowth: Performs the bacterial growth simulation and analysis.
 
# VolumeInformation 

Manages volume-related parameters. 

Parameters 

total_volume (float): Total volume of the sample in mL. 

sample_volume (float): Volume of the sample in mL. 

dilute_volume (float): Volume of the diluent in mL. 

# TimeInformation 

Manages time-related parameters. 

Parameters 

target_time (int or float): Time in minutes for the desired CFU count. 

lagged_time (int or float): Time in minutes for the bacterial population to stay constant. 

double_time (int or float): Time in minutes for the bacterial population to double. 

final_time (int or float): Final simulation time in minutes. 

# CellInformation 

Manages cell-related parameters. 

Parameters 

cell_mass (float): Mass of a single bacterial cell in micrograms. 

cell_per_volume (float): Number of cells per mL when saturated. 

# BacterialGrowth 

Performs bacterial growth simulation and analysis. 

# Methods 

validate_inputs(): Validates the input parameters. 

bacterial_growth(): Performs the bacterial growth simulation and returns the read and dilute conversion arrays. 

Error Handling 

The script includes error handling to manage issues that may arise during execution, such as invalid input values. This ensures that the simulation runs smoothly and provides informative error messages when issues occur. 

# Inputs:  

initial_read: Initial read from spectrometer in micrograms/mL

data: Data collected using BioSpectrometer kinetic instrument

total_volume: Total volume of sample in mL

sample_volume: Volume of the sample in microliters

dilute_volume: Volume of diluent added in microliters

lagged_time: Time in minutes for the bacterial population to stay constant

target_time: Time in minutes for the desired Colony-Forming Units size

double_time: Time in minutes for bacterial population to double

final_time: Final simulation time in minutes

cell_mass: Mass of a single bacterial cell in micrograms

cell_per_volume: Number of cells per mL when saturated 

# Calculations:

The average data is calculated and plotted.

The calculation accounts for the lagged time and exponential growth of bacteria restricted by the original total volume.

The estimated double time for bacteria to double its size in minutes.

dilute_conversion: Corrected concentration considering dilution factor (unitless)

rate: Growth rate in 1/minute

max_od: Maximum concentration in micrograms/mL

1 OD600 = saturation
1 OD600 = 10 microgram/mL

# Outputs:

Prints all input parameters

Prints target time with the predicted reading in micrograms/mL and colony-forming units

Plots population size scaling up to the original total volume the sample was collected from in Colony Forming Units vs time in minutes

Plots population size in diluted sample size in micrograms/mL vs time in minutes

Using linear interpolation, plots the desired target time and target amount on both plots.

Plots data if given, which can be plotted with the model or alone.

Linear x-axis and linear y-axis

Returns read array and dilute conversion array

# Contributions 

Contributions are welcome! If you find any bugs or have suggestions for improvements, please open an issue or submit a pull request.  

  

## License 

  

  

## Acknowledgements 

Students and Professors at University of South Florida 

  

### Explanation TO-DO-LIST 

1. **Overview**: Brief introduction to the project and its purpose. 

2. **Features**: List of main features of the script. 

3. **Installation**: Instructions on how to install necessary dependencies. 

4. **Usage**: Example code demonstrating how to use the script. 

5. **Classes and Methods**: Detailed information about the main classes and methods. 

6. **Error Handling**: Information on how the script handles errors. 

7. **Contributions**: Guidelines for contributing to the project. 

8. **License**: License information. 

9. **Acknowledgements**: Credits and acknowledgments. 

Do list:
    1. lagged time, death rate, death time
    2. How does different Broth impact growth rate?
    3. Calculation of lagged, double, and stable time in def calculation_measurements
    4. create a separate script for different bacteria
    5. create a separate script for different media
    6. make most inputs calculated if data is provided. 
    7. * find a better method to calculate double time.
