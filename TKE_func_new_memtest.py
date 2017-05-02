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

import os

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]


def memory(since=0.0):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:') - since


def resident(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since


def stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since

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

def main(datapath,dataname,layer,terms_dict):

	if terms_dict.get('print_stuff'):
		print 'Mem usage before loading p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'memory=',memory()/1024/1024 #in MB
		print 'resdient=',resident()/1024/1024

	### Load QGCM data
	p_data = Dataset(datapath+dataname)
	p = p_data.variables['p']
	del p_data
	p = np.transpose(np.squeeze(p,(2,1,0))) # dimensions: (x,y,time)
	if terms_dict.get('print_stuff'):
		print 'p.shape=',p.shape
		print 'Mem usage after loading p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'memory=',memory()/1024/1024 #in MB
		print 'resdient=',resident()/1024/1024

	# Define dx and dy
	#dx = terms_dict.get('dx')[0]
	#dy = terms_dict.get('dx')[1]

	# Take derivative of p
	p_x = avg_dim(ddx(p,terms_dict.get('dx')[0]),'y',1)
	p_y = avg_dim(ddy(p,terms_dict.get('dx')[1]),'x',1)

	if terms_dict.get('print_stuff'):
		print 'Mem usage after first deriv of p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'memory=',memory()/1024/1024 #in MB
		print 'resdient=',resident()/1024/1024
	
	# Take second derivative of p
	p_xx = avg_dim(ddx(p_x,terms_dict.get('dx')[0]),'y',1)
	p_yy = avg_dim(ddy(p_y,terms_dict.get('dx')[1]),'x',1)

	if terms_dict.get('print_stuff'):
		print 'Mem usage after second deriv of p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'memory=',memory()/1024/1024 #in MB
		print 'resdient=',resident()/1024/1024

	# Calculate del2(p)
	del2p = p_xx + p_yy
	if terms_dict.get('print_stuff'):
		print 'Mem usage after del2p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'memory=',memory()/1024/1024 #in MB
		print 'resdient=',resident()/1024/1024

	del p_xx,p_yy

	if terms_dict.get('print_stuff'):
		print 'Mem usage after delete p_xx and p_yy =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'memory=',memory()/1024/1024 #in MB
		print 'resdient=',resident()/1024/1024

	# Calculate derivative of del2p
	del2p_x = avg_dim(ddx(del2p,terms_dict.get('dx')[0]),'y',1)
	del2p_y = avg_dim(ddy(del2p,terms_dict.get('dx')[1]),'x',1)
	if terms_dict.get('print_stuff'):
		print 'Mem usage after del2p_x =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'memory=',memory()/1024/1024 #in MB
		print 'resdient=',resident()/1024/1024

	del del2p
	if terms_dict.get('print_stuff'):
		print 'Mem usage after delete del2p =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'memory=',memory()/1024/1024 #in MB
		print 'resdient=',resident()/1024/1024

	# Calculate Jacobian
	J = del2p_x * avg_dim(p_y,'xy',2) - del2p_y * avg_dim(p_x,'xy',2)
	del p_x,p_y,del2p_x,del2p_y

	if terms_dict.get('print_stuff'):
		print 'p.shape=',p.shape
		print 'J.shape=',J.shape
		print 'Mem usage before transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	# Ensure that dimensions are the same
	p = avg_dim(p,'xy',3)

	print 'Before calculating transfer ',datetime.now().time()

	if terms_dict.get('spatial_flag'):
		transfer_iso,kiso,ktiso = calc_T_func.main(p,J,terms_dict)
		del kiso
	else:
		transfer_iso,kiso,ktiso = calc_T_func.main(p,J,terms_dict)
	del p,J

	print 'After calculating transfer ',datetime.now().time()

	# Multiply by the correct constants: - H_layer / ((f0**3) * Htot)
	fac = -(terms_dict.get('H')[layer-1])/((terms_dict.get('f0')**3) * (terms_dict.get('H')[3]))
	transfer_iso = fac * transfer_iso

	if terms_dict.get('print_stuff'):
		print 'Mem usage after transfer func =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000
		print 'ktiso.shape=',ktiso.shape

	# For correct units, I need to scale k and w
	if not terms_dict.get('spatial_flag'):
		kiso_plot = 1000*kiso
	ktiso_plot = 60*60*24*ktiso

	if terms_dict.get('save_data'):
		x1 = terms_dict.get('domain')[0]
		x2 = terms_dict.get('domain')[1]
		y1 = terms_dict.get('domain')[2]
		y2 = terms_dict.get('domain')[3]
		yrs = terms_dict.get('yrs')
		save_name = terms_dict.get('save_name')
		extra_name = terms_dict.get('extra_name')
		if terms_dict.get('spatial_flag'):
			#Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/TKE_spatial_1yr_test_layer1_yr159_dg2_output037.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/TKE_spatial'+save_name+extra_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc', 'w', format='NETCDF3_CLASSIC')
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
			#Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/TKE_1yr_test_layer1_yr159_dg2_output037.nc', 'w', format='NETCDF3_CLASSIC')
			Tgrp = Dataset('/g/data/v45/pm2987/netcdf_transfers/TKE_'+save_name+extra_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc', 'w', format='NETCDF3_CLASSIC')
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

	print 'After saved data ',datetime.now().time()
	'''
	# Plot (temporarily to test how it looks)
	# Define dk and dw
	if not terms_dict.get('spatial_flag'):
		dk = kiso_plot[-1] - kiso_plot[-2]
	dw = ktiso_plot[-1] - ktiso_plot[-2]

	if terms_dict.get('spatial_flag'):
		# Average over spatial dimensions
		transfer_iso = np.mean(np.mean(transfer_iso,axis=0),axis=0)
		print 'transfer_iso.shape after averaging = ',transfer_iso.shape

	# Specify smaller font size for screen viewing
	font = {'family' : 'normal',
	        'size'   : 10}
	matplotlib.rc('font', **font)

	if terms_dict.get('spatial_flag'):
		# Create figure 2
		plt.figure(num=1, figsize=(15,12))
		# KEs
		plt.plot(ktiso_plot,transfer_iso/dw,linewidth=5.0,color='m',label='KE1')
		plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
		plt.xscale('log')
		plt.axis('tight')
		plt.title('KE Energy budget, integrated over wavenumber')
		plt.xlabel('Frequency')
		plt.ylabel('(nW/kg)/(rad/day)')
		plt.legend()

	else:
		# Create figure 1
		plt.figure(num=1, figsize=(12,8))
		# KEs
		plt.plot(kiso_plot,np.sum(transfer_iso,axis=1)/dk,linewidth=5.0,color='m',label='KE1')
		plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
		plt.xscale('log')
		plt.axis('tight')
		plt.title('KE Energy budget, integrated over frequency')
		plt.xlabel('Wavenumber')
		plt.ylabel('(nW/kg)/(rad/km)')
		plt.legend()

		# Create figure 2
		plt.figure(num=2, figsize=(15,12))
		# KEs
		plt.plot(ktiso_plot,np.sum(transfer_iso,axis=0)/dw,linewidth=5.0,color='m',label='KE1')
		plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
		plt.xscale('log')
		plt.axis('tight')
		plt.title('KE Energy budget, integrated over wavenumber')
		plt.xlabel('Frequency')
		plt.ylabel('(nW/kg)/(rad/day)')
		plt.legend()

	if terms_dict.get('print_stuff'):
		print 'Mem usage before plots =',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000

	print 'Before TKE plots ',datetime.now().time()
	'''

	#plt.show()










	
