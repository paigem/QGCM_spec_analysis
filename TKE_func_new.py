import numpy as np
import math
import scipy.io

import matplotlib
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
from scipy.ndimage import gaussian_filter1d

from netCDF4 import Dataset
from datetime import datetime
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
			return 0.5 * (var[:-1,:,:] + var[1:,:,:])
		if axis == 'y':
			return 0.5 * (var[:,:-1,:] + var[:,1:,:])
		if axis == 'xy'
			return 0.5 * (var[:-1,:-1,:] + var[1:,1:,:])


#-----------------------------------------------------------------------

def main(datapath,dataname,print_stuff):

	### Load QGCM data
	p_data = Dataset(datapath+dataname)
	p_grid = p_data['p']
	del p_data
	print 'p_grid.shape=',p_grid.shape

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
	J = del2p_x * p_y + del2p_y * p_x
	del p_x,p_y,del2p_x,del2p_y

	if print_stuff:
		print 'p.shape=',p.shape
		print 'J.shape=',J.shape

	# Window relevant terms (J and p)
	#p = window_func(avg_dim(p,'xy',3))
	#J = window_func(J)
	











	
