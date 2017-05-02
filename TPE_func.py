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
import calc_T_func
import calc_Tspatial_func
#-----------------------------------------------------------------------
# Functions used in the code

# Take x-derivative and average y-axis
def ddx(var,dx):
	return (1./dx) * (var[:-1,:,:] - var[1:,:,:])
	
# Take y-derivative
def ddy(var,dy):
	return (1./dy) * (var[:,:-1,:] - var[:,1:,:]) 

# Average specified dimension(s)
def avg_dim(var,axis,number):
	for i in np.arange(number):
		if axis == 'x':
			var = 0.5 * (var[:-1,:,:] + var[1:,:,:])
		if axis == 'y':
			var = 0.5 * (var[:,:-1,:] + var[:,1:,:])
		if axis == 'xy':
			var = 0.5 * (var[:-1,:-1,:] + var[1:,1:,:])
	return var


#-----------------------------------------------------------------------

def main(datapath,dataname1,dataname2,layer,terms_dict):

	### Load QGCM data
	p_data = Dataset(datapath+dataname1)
	p1 = p_data.variables['p']
	del p_data
	p1 = np.transpose(np.squeeze(p1,(2,1,0))) # dimensions: (x,y,time)
	if terms_dict.get('print_stuff'):
		print 'p1.shape=',p1.shape
		print 'Mem usage after loading p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	### Load QGCM data
	p_data = Dataset(datapath+dataname2)
	p2 = p_data.variables['p']
	del p_data
	p2 = np.transpose(np.squeeze(p2,(2,1,0))) # dimensions: (x,y,time)

	# Define dx and dy
	dx = terms_dict.get('dx')[0]
	dy = terms_dict.get('dx')[1]

	# Take derivative of p
	p1_x = avg_dim(ddx(p1,dx),'y',1)
	p2_x = avg_dim(ddx(p2,dx),'y',1)
	p1_y = avg_dim(ddy(p1,dy),'x',1)
	p2_y = avg_dim(ddy(p2,dy),'x',1)
	
	# Take the difference of the two pressures
	p_diff = p1 - p2
	del p1,p2

	# Calculate Jacobian
	J = p1_x * p2_y - p1_y * p2_x
	del p1_x,p1_y,p2_x,p2_y

	if terms_dict.get('print_stuff'):
		#print 'p1.shape=',p1_x.shape
		print 'J.shape=',J.shape
		print 'Mem usage before transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	# Ensure that dimensions are the same
	p_diff = avg_dim(p_diff,'xy',1)

	if terms_dict.get('spatial_flag'):
		transfer_iso,kiso,ktiso = calc_T_func.main(p_diff,J,terms_dict)
		del kiso
	else:
		transfer_iso,kiso,ktiso = calc_T_func.main(p_diff,J,terms_dict)
	del p_diff,J

	if terms_dict.get('print_stuff'):
		print 'Mem usage after transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	# Multiply by the correct constants: 1 / ((f0) * Htot * gprime_layer)
	fac = 1.0/((terms_dict.get('f0')) * terms_dict.get('H')[3] * terms_dict.get('gprimes')[layer-1])
	transfer_iso = fac * transfer_iso

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
		if terms_dict.get('spatial_flag'):
			#Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/TPE_spatial_1yr_test_layer1_yr159_dg2_output037.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/TPE_spatial'+save_name+extra_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp.createDimension('x',transfer_iso.shape[0])
			Tgrp.createDimension('y',transfer_iso.shape[1])
			Tgrp.createDimension('w',transfer_iso.shape[2])
			T = Tgrp.createVariable('T','f4',('x','y','w'))
			T[:,:,:] = transfer_iso
			Tgrp.createDimension('ktiso_dim',len(ktiso_plot))
			ktiso = Tgrp.createVariable('ktiso','f4',('ktiso_dim'))
			ktiso[:] = ktiso_plot
			Tgrp.close()
		else:
			#Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/TPE_1yr_test_layer1_yr159_dg2_output037.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/TPE'+save_name+extra_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp.createDimension('k',transfer_iso.shape[0])
			Tgrp.createDimension('w',transfer_iso.shape[1])
			T = Tgrp.createVariable('T','f4',('k','w'))
			T[:,:] = transfer_iso
			Tgrp.createDimension('kiso_dim',len(kiso_plot))
			kiso = Tgrp.createVariable('kiso','f4',('kiso_dim'))
			kiso[:] = kiso_plot
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

	print 'Before TKE plots ',datetime.now().time()

	if terms_dict.get('spatial_flag'):

		# Create figure 2
		plt.figure(num=2, figsize=(15,12))
		# KEs
		plt.plot(ktiso_plot,transfer_iso/dw,linewidth=5.0,color='m',label='PE12')
		plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
		plt.xscale('log')
		plt.axis('tight')
		plt.title('PE Energy budget, integrated over wavenumber')
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
		plt.title('PE Energy budget, integrated over frequency')
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
		plt.title('PE Energy budget, integrated over wavenumber')
		plt.xlabel('Frequency')
		plt.ylabel('(nW/kg)/(rad/day)')
		plt.legend()

	if terms_dict.get('print_stuff'):
		print 'Mem usage before plotting =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
	'''
	#plt.show()










	
