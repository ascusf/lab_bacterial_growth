#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 11:31:10 2024

- LB broth has carbincillin, pH 7.0 
- E.coli culture were examined from agar plates.
    - Wild Type CSF1-myc-aga, CSF1 NNK 'F' S.C., CSF1 NNV vector 'F' DH5alpha
- 37 degC at 280 RPM
- Eppendorf BioSpectrometer kinetic
- 1.5 mL cuvettes
- 12 mL tubes
- 200 micro mL were collected every measurement taken.

@author: andrewchadwell
"""
#---------7/2/24---------------------------------------------------------------
# ecoli_initial_37c: Ecoli data collected on 7/2/24 3 hour grow from 4degC LB 
# agar plate. Sample plate Wild type CSF1-myc-aga E. Coli culture 6/25/24. After
# three hours in a 37degC shake, the data was collected every 20 minutes. 12 mL
# tubes and 4 mL of sample are used.  
ecoli_initial_37c_7_2_24 = {
    "minute": [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 1440],
    "total_volume": [4],
    "1": [None, 0.050, 0.051, 0.040, 0.060, 0.120, 0.130, 0.140, 0.230, 0.230, 
          0.400, 0.390, 0.470, 4.210],
    "2": [None, 0.003, 0.200, 0.300, 0.200, 0.230, 0.230, 0.240, 0.270, 0.310,
          0.420, 0.480, 0.610, 4.260],
    "3": [None, 0.030, 0.056, 0.390, 0.220, 0.180, 0.180, 0.180, 0.210, 0.350, 
          0.410, 0.500, 0.570, 4.120],
}

# ecoli_initial_37c: Ecoli data collected on 7/2/24 3 hour grow from 4degC LB 
# agar plate. Sample plate Wild type CSF1-myc-aga E. Coli culture 6/25/24. After
# three hours in a 37degC shake. It was cooled to 4C then return to the 37C 
# shaker. The data was collected every 20 minutes. 12 mL tubes and 4 mL of 
# sample are used.   
ecoli_initial_4c_7_2_24 = {
    "minute": [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 1440],
    "total_volume": [4],
    "1": [None, .03, .015, .110, .13, .1, .18, .21, .23, .37, .53, .63, .72,
          4.44],
    "2": [None, .007, .28, .274, .3, .28, .3, .34, .37, .37, .9, 1.01,1.2, 4.78],
    "3": [None, .002, .02, .07, .12, .14, .17, .24, .28, .3, .45, .53, .66, 4.9],
}

#---------7/8/24---------------------------------------------------------------
# ecoli_initial_37c: Ecoli data collected on 7/8/24 3 hour grow from 4degC LB 
# agar plate. Sample plate Wild type CSF1-myc-aga E. Coli culture 6/25/24. After
# three hours in a 37degC shake, the data was collected. 12 mL tubes and 4 mL 
# of sample are used.  
ecoli_initial_37c_7_8_24 = {
    "minute": [0, 60, 80, 100, 120, 140, 160],
    "total_volume": [4],
    "1": [.06,.03,.313,.151,.25, None, None],
    "2": [.06,None,.05,.15,.202,.32,.4]
}

# ecoli_initial_4c: Ecoli data collected on 7/8/24 3 hour grow from 4degC LB 
# agar plate. Sample plate Wild type CSF1-myc-aga E. Coli culture 6/25/24. After
# three hours in a 37degC shake. It was cooled to 4C then return to the 37C 
# shaker. The data was collected every 20 minutes. 12 mL tubes and 4 mL of 
# sample are used.  
ecoli_initial_4c_7_8_24 = {
    "minute": [0, 60, 80, 100, 120, 140, 160],
    "total_volume": [4],
    "1": [.06,None,.01,.118,.272,.41,.40],
    "2": [.06,.18,.113,.24,.488,.54,.599]
}

#---------7/11/24---------------------------------------------------------------
# ecoli_initial_4c_7_11_24: Ecoli data collected on 7/11/24, grew from 4degC LB 
# agar plate. Put in shaker 37C and 280 RPM. Sample plate CSF1 NNK 'F' S.C. 
# The data was collected every 60 minutes. 12 mL tubes and 4.00 mL of sample are used.   
ecoli_initial_4c_7_11_24 = {
    "minute": [0, 60, 120, 180, 240, 300, 360, 420],
    "total_volume": [4],
    "1": [.70,.945,1.50,2.30,3.35,4.35,4.25,4.37],
    "2": [.62,.76,1.30,2.15,3.10,3.40,4.50,4.45]
}

#---------7/15/24---------------------------------------------------------------
# ecoli_initial_37c_7_15_24: Ecoli data collected on 7/15/24, grew from 4degC LB 
# agar plate. Sample plate CSF1 NNK 'F' S.C. After. Put in shaker 37C and 280 RPM.
# The data was collected every 60 minutes. 12 mL tubes and 11 mL of 
# sample are used 8 mL remain at end. The samples was diluted to 1 micrograms/mL. 
ecoli_7_15_24 = {
    "minute": [0, 60, 120, 180, 240, 300, 360, 420, 480, 540],
    "total_volume": [8],
    "1": [None,1.20,1.55,1.98,2.18,2.65,3.00,3.00,4.25,3.21],
    "2": [None,1.33,1.70,2.20,2.55,2.95,3.15,3.60,6.30,4.10]
}
  
# ecoli: Ecoli data collected on 7/8/24 8:22 grow from 4degC LB. Liquid in a 37C 
# shaker. Sample plate Wild type CSF1-myc-aga. The data was collected around every 
# hour. '1' 100 mL Bottle with initial 12 mL sample. '2' 250 mL flask with initial 
# 12 mL sample.
ecoli_7_17_24 = {
    "minute": [0, 60, 120, 180, 240, 300, 360, 433, 493],
    "total_volume": [10.20],
    "1": [.80,1.20,1.75,2.74,3.64,4.20,3.8,3.9,4.0],
    "2": [1.07,1.30,1.81,2.88,3.81,4.27,4,4.01,3.7]
}

