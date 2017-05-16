# buoyancy_func.py: Does calculation for buoyancy term in spectral transfer equation
#	- is called by run_Tkw.py
#	- uses xarray 
#	- calls transfer function for each individual row, then adds to current sum of rows 
#----------------------------------------------------------------------

import numpy as np
import math
import scipy.io

import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
from scipy.ndimage import gaussian_filter1d
import resource

import xarray as xr
from datetime import datetime
import window_func
import detrend_func
import calc_T_func_longTime
#-----------------------------------------------------------------------
def main(datapath,terms_dict,i,j):

	### Load pressure1
	#with xr.open_mfdataset(datapath) as ds:
	ds = xr.open_mfdataset(datapath)
	p1 = ds.p[:-1,0,i,j] # dimensions: [time, layer, x (lon), y (lat)]
	p2 = ds.p[:-1,1,i,j]
	e = ds.e[:-1,i,j] # dimensions: [time, x (lon), y (lat)]

	if terms_dict.get('print_stuff'):
		print 'p1.shape=',p1.shape
		print 'Mem usage after loading p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	# Take the difference of the two pressures
	p_diff = p2 - p1
	del p1,p2

	if terms_dict.get('print_stuff'):
		print 'Mem usage before transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	# Call the transfer function
	transfer_iso,ktiso = calc_T_func_longTime.main(p_diff,e,terms_dict)
	del p_diff,e

	# Multiply by the correct constants: 1 / ((f0) * Htot * gprime_layer)
	fac = -1.0/(terms_dict.get('H')[3])
	transfer_iso = fac * transfer_iso
	ds.close()
	if terms_dict.get('print_stuff'):
		print 'Mem usage after transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	# For correct units, I need to scale k and w
	if not terms_dict.get('spatial_flag'):
		kiso = 1000*kiso
	ktiso = 60*60*24*ktiso

	if terms_dict.get('save_data'):
		x1 = terms_dict.get('domain')[0]
		x2 = terms_dict.get('domain')[1]
		y1 = terms_dict.get('domain')[2]
		y2 = terms_dict.get('domain')[3]
		yrs = terms_dict.get('yrs')
		save_name = terms_dict.get('save_name')
		extra_name = terms_dict.get('extra_name')
		transfer_datapath = terms_dict.get('transfer_datapath')

		# Open, add, and resave to previous timeseries (for all but first i)
		if i > x1 or j > 0:

			with xr.open_dataset(transfer_datapath+'buoyancy_longTime_xarray'+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc') as ds_oldT:

				transfer_iso = transfer_iso +  ds_oldT.transfer[:]
		
			if terms_dict.get('print_stuff'):
				print 'transfer_iso.shape after addition = ',transfer_iso.shape

			# If on last i, divide by i
			if i == (x2-1) and j == (y2-1):
				transfer_iso = transfer_iso[:]/float(i)
				if terms_dict.get('print_stuff'):
					print 'i and j = ',i,j
			if terms_dict.get('print_stuff'):
				print 'i and j = ',i,j

		ds_new = xr.Dataset({'transfer':('w',transfer_iso)},coords={'w':ktiso})
		ds_new.to_netcdf(transfer_datapath+'buoyancy_longTime_xarray'+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
		ds_new.close()











	
