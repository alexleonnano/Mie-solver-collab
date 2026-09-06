#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday September 02 2026

Intened for a first pass of using minimisation of residual 
squares to find the best fitting spectra to experimentally measured data

@author: william.mcm.p
"""

from src import load_material as tools
from src import lambda_interpolate as interp
import numpy as np
import matplotlib.pyplot as plt
import miepython as mie
import pandas as pd

# - Helper functions
def normalised_data(data):
    """Noramlise the data between [0,1] - min-max method"""
    return (data - min(data))/(max(data)- min(data))

material = "materials/Silicon - crystalline - Green 2008.txt"
refrative_index_df = pd.read_csv(material,  sep=r'\s+', header = 0, names = ["wl", "n", "k"], dtype = float)

exp_data_path = 'experimental_data/SiNP_crystaline_scat.csv'
exp_data_df = pd.read_csv(exp_data_path, header=0,names=['wl', 'scat'])


# rnage of to be fitted particles
radius_range = np.linspace(5, 200, 200) * 1e-3 # in micron

# Pre allocated rms arrat
rms = np.zeros_like(radius_range)

# Loop over each particle size
#  -> 
