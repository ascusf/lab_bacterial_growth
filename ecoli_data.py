#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 11:31:10 2024

LB broth has carbincillin, pH 7.0 


@author: andrewchadwell
"""
#---------7/2/24---------------------------------------------------------------
# ecoli_initial_37c: Ecoli data collected on 7/2/24 3 hour grow from 4degC LB 
# agar plate. Sample plate Wild type CSF1-myc-aga E. Coli culture 6/25/24. After
# three hours in a 37degC shake, the data was collected every 20 minutes.
ecoli_initial_37c_7_2_24 = {
    "minute": [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 1440],
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
# shaker. The data was collected every 20 minutes. 
ecoli_initial_4c_7_2_24 = {
    "minute": [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 1440],
    "1": [None, .03, .015, .110, .13, .1, .18, .21, .23, .37, .53, .63, .72,
          4.44],
    "2": [None, .007, .28, .274, .3, .28, .3, .34, .37, .37, .9, 1.01,1.2, 4.78],
    "3": [None, .002, .02, .07, .12, .14, .17, .24, .28, .3, .45, .53, .66, 4.9],
}

#---------7/8/24---------------------------------------------------------------
# ecoli_initial_37c: Ecoli data collected on 7/8/24 3 hour grow from 4degC LB 
# agar plate. Sample plate Wild type CSF1-myc-aga E. Coli culture 6/25/24. After
# three hours in a 37degC shake, the data was collected.
ecoli_initial_37c_7_8_24 = {
    "minute": [0, 60, 80, 100, 120, 140, 160],
    "1": [.06,.03,.313,.151,.25, None, None],
    "2": [.06,None,.05,.15,.202,.32,.4]
}

# ecoli_initial_4c: Ecoli data collected on 7/8/24 3 hour grow from 4degC LB 
# agar plate. Sample plate Wild type CSF1-myc-aga E. Coli culture 6/25/24. After
# three hours in a 37degC shake. It was cooled to 4C then return to the 37C 
# shaker. The data was collected every 20 minutes. 
ecoli_initial_4c_7_8_24 = {
    "minute": [0, 60, 80, 100, 120, 140, 160],
    "1": [.06,None,.01,.118,.272,.41,.40],
    "2": [.06,.18,.113,.24,.488,.54,.599]
}

data_test = {
    "minute": [0, 60, 80, 100, 120, 140, 160],
    "1": [.06,None,.01,.118,.272,.41,.40]
}

# ecoli_initial_4c_7_11_24: Ecoli data collected on 7/11/24, grew from 4degC LB 
# agar plate. Sample plate Wild type CSF1-myc-aga E. Coli culture 6/25/24. After
# three hours in a 37degC shake. It was cooled to 4C then return to the 37C 
# shaker. The data was collected every 20 minutes. 
ecoli_initial_4c_7_11_24 = {
    "minute": [0, 60, 120, 180, 240, 300, 360, 420],
    "1": [.70,.945,1.50,2.30,3.35,4.35,4.25,4.37],
    "2": [.62,.76,1.30,2.15,3.10,3.40,4.50,4.45]
}