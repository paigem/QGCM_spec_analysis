import numpy as np
import math
import scipy.io

import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
from scipy.ndimage import gaussian_filter1d

from netCDF4 import Dataset
from datetime import datetime
import window_func
import detrend_func
import calc_T_func
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

def main(datapath,dataname,print_stuff,spacetime,padding_fac,kfac):

	### Load QGCM data
	p_data = Dataset(datapath+dataname)
	p = p_data['p']
	del p_data
	p = np.transpose(np.squeeze(p,(2,1,0))) # dimensions: (x,y,time)
	if print_stuff:
		print 'p.shape=',p.shape

	# Define dx and dy
	dx = 5000.
	dy = 5000.

	# Take derivative of p
	p_x = avg_dim(ddx(p,dx),'y',1)
	p_y = avg_dim(ddy(p,dy),'x',1)
	
	# Take second derivative of p
	p_xx = avg_dim(ddx(p_x,dx),'y',1)
	p_yy = avg_dim(ddy(p_y,dy),'x',1)

	# Calculate del2(p)
	del2p = p_xx + p_yy
	del p_xx,p_yy

	# Calculate derivative of del2p
	del2p_x = avg_dim(ddx(del2p,dx),'y',1)
	del2p_y = avg_dim(ddy(del2p,dy),'x',1)
	del del2p

	# Calculate Jacobian
	J = del2p_x * avg_dim(p_y,'xy',2) + del2p_y * avg_dim(p_x,'xy',2)
	del p_x,p_y,del2p_x,del2p_y

	if print_stuff:
		print 'p.shape=',p.shape
		print 'J.shape=',J.shape

	# Ensure that dimensions are the same
	p = avg_dim(p,'xy',3)

	calc_T_func.main(p,J,spacetime,padding_fac,kfac)

	print 'Got to the end!'
	











	
