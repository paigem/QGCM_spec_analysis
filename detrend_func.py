import numpy as np
from scipy import signal

#--------------------------------------------------------------------------------------------
# Function for 2d spatial detrend
def spatialdetrend_2d_QGCM(var):

	# Loop through the time dimension
	for t in np.arange(var.shape[2]):
		
		# Define matrix A ('G' in Brian's code)
		
		x,y = np.meshgrid(np.arange(var.shape[0]),np.arange(var.shape[1]))
		
		# Flatten the arrays
		x = x.flatten()
		y = y.flatten()

		A = np.column_stack((x,y,np.ones(x.shape)))
		z = var[:,:,t]
		z = z.flatten()
	
		# Calculate the least squares coefficients
		c = np.linalg.lstsq(A,z)[0]
		
		zfit = c[0]*x + c[1]*y +c[2]
		
		# Reshape zfit and z
		zfit = np.reshape(zfit,[var.shape[0],var.shape[1]])
		z = np.reshape(z,[var.shape[0],var.shape[1]])
		
		# Subtract off the least squares plane and any nonzero mean
		var_temp = z-zfit-np.mean(z-zfit)
		if t==0:
			var_2d_detrend = var_temp
		else:
			var_2d_detrend = np.dstack((var_2d_detrend,var_temp))
	del var_temp,c,zfit,x,y,z,A
	return var_2d_detrend
#--------------------------------------------------------------------------------------------


def main(var,spacetime):
	if spacetime == 'time' or spacetime == 'spacetime':
		print 'Detrending in time'
		# Detrend in time (i.e. subtract off time mean)
		var = signal.detrend(var,axis=2,type='linear') # dimensions still (time,lat,lon)
		print 'After time detrend ',var.shape
	
		
	if spacetime == 'space' or spacetime == 'spacetime':
		print 'Detrending in space'
		# Do 2d spatial detrend
		var = spatialdetrend_2d_QGCM(var) # Note that outputted dimensions are now: (lat,lon,time)
		print 'After spatial detrend var.shape=',var.shape

	return var


