#--------------------------------------------------------------------------------------------------------------------------------------------------
# This code runs the TKE and TPE codes for each of the layers, and then calls on the energy_budget_transfers.py code to plot the energy budget
#
#--------------------------------------------------------------------------------------------------------------------------------------------------

import numpy as np
import TKE_func_new
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os.path

# Specify parameters
fluid_var = 'oc-coupled'
climate_var = 'p'
domain = 'partial domain'
yrs = [55,55]

if domain == 'fulldomain':
	y1 = 0
	y2 = 960
	x1 = 0
	x2 = 960
	datatitle = 'fulldomain'
else:
	y1 = 0
	y2 = 100
	x1 = 0
	x2 = 100

print_stuff = 1
spacetime = 'spacetime'
padding_fac = 1.0
kfac = 100 # number of wavenumbers desired

# Define constants to be used
dx = 5000. # meters
dy = 5000. # meters
H1 = 350.0 # meters
H2 = 750.0 # meters
H3 = 2900.0 # meters
Htot = H1 + H2 + H3

datapath = '/g/data/v45/pm2987/nco_and_output/ocean/'
dataname = 'p_Daily_1yr_dg2_output037_layer1_0_100_0_100_159_159_50daytest.nc'

#--------------------------------------------------------------------------------------------------------------------------------------------------
# Call TKE_func_new.py
TKE_func_new.main(datapath,dataname,print_stuff,spacetime,padding_fac,kfac)
