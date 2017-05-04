# buoyancy.py: Does calculation for buoyancy term in spectral transfer equation
#	- is called by run_Tkw.py
#----------------------------------------------------------------------

import numpy as np
import math
import scipy.io

import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
from scipy.ndimage import gaussian_filter1d
import resource

from netCDF4 import Dataset
from datetime import datetime
import window_func
import detrend_func
import calc_T_func_longTime
#-----------------------------------------------------------------------
def main(datapath,terms_dict,i,j):

	### Load pressure1
	p_data = Dataset(datapath+'ocpo.nc')
	p1 = p_data.variables['p'][:,0,j,i] # layer 1
	del p_data
	#p1 = np.transpose(np.squeeze(p1,(2,1,0))) # dimensions: (x,y,time)
	if terms_dict.get('print_stuff'):
		print 'p1.shape=',p1.shape
		print 'Mem usage after loading p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	### Load pressure2
	p_data = Dataset(datapath+'ocpo.nc')
	p2 = p_data.variables['p'][:,1,j,i] # layer 2
	del p_data
	#p2 = np.transpose(np.squeeze(p2,(2,1,0))) # dimensions: (x,y,time)

	# Take the difference of the two pressures
	p_diff = p2 - p1
	del p1,p2

	### Load entrainment
	e_data = Dataset(datapath+'ocpo.nc')
	e = e_data.variables['e'][:,j,i]
	del e_data
	#e = np.transpose(np.squeeze(e,(2,1,0))) # dimensions: (x,y,time)

	if terms_dict.get('print_stuff'):
		print 'Mem usage before transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	# Call the transfer function
	transfer_iso,ktiso = calc_T_func_longTime.main(p_diff,e,terms_dict)
	del p_diff,e

	# Multiply by the correct constants: 1 / ((f0) * Htot * gprime_layer)
	fac = -1.0/(terms_dict.get('H')[3])
	transfer_iso = fac * transfer_iso

	if terms_dict.get('print_stuff'):
		print 'Mem usage after transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	# For correct units, I need to scale k and w
	if not terms_dict.get('spatial_flag'):
		kiso_plot = 1000*kiso
		del kiso
	ktiso_plot = 60*60*24*ktiso
	del ktiso

	if terms_dict.get('save_data'):
		x1 = terms_dict.get('domain')[0]
		x2 = terms_dict.get('domain')[1]
		y1 = terms_dict.get('domain')[2]
		y2 = terms_dict.get('domain')[3]
		yrs = terms_dict.get('yrs')
		save_name = terms_dict.get('save_name')
		extra_name = terms_dict.get('extra_name')
		'''
		# Save first timeseries
		if i==0:
			Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/buoyancy_spatial'+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp.createDimension('w',transfer_iso.shape[0])
			T = Tgrp.createVariable('T','f4',('w'))
			T[:] = transfer_iso
			Tgrp.createDimension('ktiso_dim',len(ktiso_plot))
			ktiso = Tgrp.createVariable('ktiso','f4',('ktiso_dim'))
			ktiso[:] = ktiso_plot
			Tgrp.close()
		'''
		# Open, add, and resave to previous timeseries (for all but first i)
		if i > x1 or j > 0:

			data_open = Dataset(terms_dict.get('transfer_datapath')+'buoyancy_longTime'+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
			old_data = data_open.variables['T']
			del data_open
			if terms_dict.get('print_stuff'):
				print 'old_data.shape = ',old_data.shape
			transfer_iso = old_data + transfer_iso
			del old_data
			# If on last i, divide by i
			if i == (x2-1) and j == (y2-1):
				transfer_iso = transfer_iso/float(i)
				if terms_dict.get('print_stuff'):
					print 'i and j = ',i,j
			if terms_dict.get('print_stuff'):
				print 'i and j = ',i,j

		Tgrp = Dataset(terms_dict.get('transfer_datapath')+'buoyancy_longTime'+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc', 'w', format='NETCDF3_CLASSIC')
		Tgrp.createDimension('w',transfer_iso.shape[0])
		T = Tgrp.createVariable('T','f4',('w'))
		T[:] = transfer_iso
		Tgrp.createDimension('ktiso_dim',len(ktiso_plot))
		ktiso = Tgrp.createVariable('ktiso','f4',('ktiso_dim'))
		ktiso[:] = ktiso_plot
		Tgrp.close()

	'''
	if terms_dict.get('spatial_flag'):
		# Average over spatial dimensions
		transfer_iso = np.mean(np.mean(transfer_iso,axis=0),axis=0)
		print 'transfer_iso.shape after averaging = ',transfer_iso.shape
	
	# Plot (temporarily to test how it looks)
	# Define dk and dw
	if not terms_dict.get('spatial_flag'):
		dk = kiso_plot[-1] - kiso_plot[-2]
	dw = ktiso_plot[-1] - ktiso_plot[-2]

	# Specify smaller font size for screen viewing
	font = {'family' : 'normal',
	        'size'   : 10}
	matplotlib.rc('font', **font)

	if terms_dict.get('spatial_flag'):
		print 'Mem usage after transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	print 'Before TKE plots ',datetime.now().time()

	if terms_dict.get('spatial_flag'):

		# Create figure 2
		plt.figure(num=2, figsize=(15,12))
		# KEs
		plt.plot(ktiso_plot,transfer_iso/dw,linewidth=5.0,color='m',label='PE12')
		plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
		plt.xscale('log')
		plt.axis('tight')
		plt.title('Buoyancy, integrated over wavenumber')
		plt.xlabel('Frequency')
		plt.ylabel('(nW/kg)/(rad/day)')
		plt.legend()

	else:
		# Create figure 1
		plt.figure(num=1, figsize=(12,8))
		# KEs
		plt.plot(kiso_plot,np.sum(transfer_iso,axis=1)/dk,linewidth=5.0,color='m',label='PE12')
		plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
		plt.xscale('log')
		plt.axis('tight')
		plt.title('Buoyancy, integrated over frequency')
		plt.xlabel('Wavenumber')
		plt.ylabel('(nW/kg)/(rad/km)')
		plt.legend()

		# Create figure 2
		plt.figure(num=2, figsize=(15,12))
		# KEs
		plt.plot(ktiso_plot,np.sum(transfer_iso,axis=0)/dw,linewidth=5.0,color='m',label='PE12')
		plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
		plt.xscale('log')
		plt.axis('tight')
		plt.title('Buoyancy, integrated over wavenumber')
		plt.xlabel('Frequency')
		plt.ylabel('(nW/kg)/(rad/day)')
		plt.legend()

	if terms_dict.get('print_stuff'):
		print 'Mem usage before plotting =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
	'''
	#plt.show()










	
