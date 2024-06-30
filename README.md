# Bacterial-Growth 
# Stern Lab - Bacterial Growth Simulation 

## Overview 

  

This project contains a script for simulating bacterial growth using various parameters for volume, time, and cell information. It includes functionality to plot growth curves and validate inputs to ensure they are within acceptable ranges. 

  

## Features 

  

- **Volume Information:** Manage and validate volume-related parameters. 

- **Time Information:** Manage and validate time-related parameters. 

- **Cell Information:** Manage and validate cell-related parameters. 

- **Plotting:** Visualize bacterial growth simulation results with plotting functionality. 

- **Validation:** Ensure inputs are valid and handle errors gracefully. 

- **Simulation:** Perform bacterial growth simulation and analyze results. 

  

## Installation 

  

To use this script, you need to have Python 3 and the following packages installed: 

  

- numpy 

- pylab 

  

You can install these packages using pip: 

  

## Usage 

  

# Here is an example of how to use the script: 

from bacterial_growth_simulation import VolumeInformation, TimeInformation, CellInformation, BacterialGrowth 

  

# Define volume information 

volume_info = VolumeInformation(total_volume=25, sample_volume=200, dilute_volume=1800) 

  

# Define time information 

time_info = TimeInformation(lagged_time=2*60, target_time=60*7, double_time=60*2.5, final_time=60*10) 

  

# Define cell information 

cell_info = CellInformation(cell_mass=47.65e-6, cell_per_volume=1E7) 

  

# Run the bacterial growth simulation 

bacterial_growth_simulation = BacterialGrowth(initial_read=0.48, volume_information=volume_info, time_information=time_info, cell_information=cell_info) 

read, dilute_conversion = bacterial_growth_simulation.bacterial_growth() 

  

## Classes and Methods 

  

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

  

# Contributions 

Contributions are welcome! If you find any bugs or have suggestions for improvements, please open an issue or submit a pull request.  

  

## License 

  

  

## Acknowledgements 

Stern's Lab at University of South Florida 

  

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
