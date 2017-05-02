# Bottom drag calculation
#-----------------------------------------------------------------------
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
# Constants used in this code

delta_ek = 2.0 #m

#-----------------------------------------------------------------------

def main(datapath,dataname,terms_dict):

	### Load QGCM data
	p_data = Dataset(datapath+dataname)
	p = p_data.variables['p']
	del p_data
	p = np.transpose(np.squeeze(p,(2,1,0))) # dimensions: (x,y,time)
	if terms_dict.get('print_stuff'):
		print 'p.shape=',p.shape

	# Take derivative of p
	p_x = avg_dim(ddx(p,terms_dict.get('dx')[0]),'y',1)
	p_y = avg_dim(ddy(p,terms_dict.get('dx')[1]),'x',1)

	
	# Take second derivative of p
	p_xx = avg_dim(ddx(p_x,terms_dict.get('dx')[0]),'y',1)
	p_yy = avg_dim(ddy(p_y,terms_dict.get('dx')[1]),'x',1)
	del p_x,p_y

	# Calculate del2(p)
	del2p = p_xx + p_yy
	del p_xx,p_yy

	# Ensure that dimensions are the same
	p = avg_dim(p,'xy',3)

	print 'Before calculating transfer ',datetime.now().time()

	if terms_dict.get('spatial_flag'):
		transfer_iso,kiso,ktiso = calc_T_func.main(p,del2p,terms_dict)
		del kiso
	else:
		transfer_iso,kiso,ktiso = calc_T_func.main(p,del2p,terms_dict)
	del p,del2p

	# Multiply by the correct constants: delek / (abs(f0) * Htot)
	fac = delta_ek/(abs(terms_dict.get('f0')*terms_dict.get('H')[3]))
	transfer_iso = fac * transfer_iso

	print 'After calculating transfer ',datetime.now().time()

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
		if terms_dict.get('spatial_flag'):
			#Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/bottomDrag_spatial_1yr_test_yr159_dg2_output037.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/bottomDrag_spatial'+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc', 'w', format='NETCDF3_CLASSIC')
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
			#Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/buoyancy_1yr_test_layer1_yr159_dg2_output037.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/bottomDrag_'+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc', 'w', format='NETCDF3_CLASSIC')
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
		plt.title('Bottom Drag, integrated over wavenumber')
		plt.xlabel('Frequency')
		plt.ylabel('(nW/kg)/(rad/day)')
		plt.legend()

	plt.show()
	'''








	
