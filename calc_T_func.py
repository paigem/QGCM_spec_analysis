import numpy as np
import math

import detrend_func
import window_func




#----------------------------------------------------------------------------------------------------
# This code takes two quantities in grid space and computes the spectral transfer:
#	- detrends both terms
#	- windows both terms
#	- takes FT of both terms
#	- multiplies complex conjugate of one by the other
#----------------------------------------------------------------------------------------------------
# Functions used later in the code

# Make a given number even
def make_odd(number):
	number = math.ceil(number)
	if not np.mod(number,2):
		number = number+1		
	return number
	
# Make psi_array odd if originally even
def make_var_odd(var):
	nx = var.shape[0] # Get size of input psi
	ny = var.shape[1]
	nt = var.shape[2]
	
	if not np.mod(nx,2):
		var = np.concatenate((var,np.zeros((1,ny,nt))),axis=0)
	
	# Redefine size after previous if statement
	nx = var.shape[0]
	
	if not np.mod(ny,2):
		var = np.concatenate((var,np.zeros((nx,1,nt))),axis=1)
	del nx,ny,nt
	return var #,nx,ny

# Create spectral domain
def spec_domain(L):
	if np.mod(L,2):
		return np.array(np.hstack([0,np.arange(1,(L-1)/2+1),-1*np.arange((L-1)/2,0,-1)]))
		print 'L is odd!'
	else: 
		return np.array(np.hstack([0,np.arange(1,L/2),-1*np.arange(L/2,0,-1)]))

# Pad given array in spectral space
def pad(field,Lx,Ly,nx,ny):
	xpad = int((Lx-nx)/2)
	ypad = int((Ly-ny)/2)
	
	print 'xpad=',xpad,'ypad=',ypad
	
	# Will need to add in the case if the domain is even!!!
	
	return np.fft.ifftshift(np.lib.pad(np.fft.fftshift(field),((xpad,xpad),(ypad,ypad),(0,0)),'constant'))


# Unpad an array
def unpad(field,Lx,Ly,nx,ny):
	xpad = int((Lx-nx)/2)
	ypad = int((Ly-ny)/2)

	out = np.fft.fftshift(field)
	#if xpad:
	#	out = out[xpad:-(xpad),ypad:-(ypad),:]
	#	print 'in the unpad if loop'
	return np.fft.ifftshift(out)

def make_iso(var1,var2,wv,kiso):

	# Make the frequencies isotropic
	field = make_wiso(var1,var2,wv)	
	del var2
	print 'field.shape=',field.shape

	# Portion adapted from Brian's code
	# Loop through time/frequency
	for m in np.arange(field.shape[2]):
	
		dummy = np.squeeze(field[:,:,m]) # squeeze matrix in 2d from 3d

		# Initialize e0,en
		e0 = 0
		en = np.zeros((len(kiso),1))

		# Loop through isotropic wavenumbers
		for n in np.arange(len(kiso)):

			cut = wv < (kiso[n])**2 # Need to add if statement for n not =1?
			e1 = sum(sum(dummy*cut.T))
			en[n] = (e1 - e0)
			e0 = e1
		
		#Try to append columns each loop
		if m==0:
			out = en
		else:
			out = np.hstack((out,en))

	print 'out.shape',out.shape
	del wv,field,e0,en
	return out


# Make frequencies isotropic (i.e. symmetric)
def make_wiso(var1,var2,wv):
	endw = float(var1.shape[2])

	# Define first and last entries of isotropic transfer
	out1 = (10**9)*np.real(np.conjugate(var1[:,:,0])*var2[:,:,0])
	out_end = (10**9)*np.real(np.conjugate(var1[:,:,0.5*endw])*var2[:,:,0.5*endw])
	
	# Define all other entries of isotropic transfer
	out = (10**9)*np.real(np.conjugate(var1[:,:,1:0.5*endw])*var2[:,:,1:0.5*endw] + np.conjugate(var1[:,:,-1:0.5*endw+1:-1])*var2[:,:,-1:0.5*endw+1:-1])

	# Concatenate the arrays together
	out_final = np.dstack((out1,out,out_end))

	return out_final

#----------------------------------------------------------------------------------------------------

def main(var1,var2,spacetime,padding_fac,kfac):


	# Define FFT normalization
	nx = var1.shape[0]
	ny = var1.shape[1]
	nt = var1.shape[2]
	dt = 1 # in days
	print 'nx = ',nx

	# Normalization factor for wavenumbers (translating from lat/lon to wavenumber)
	kx_units = 2.*np.pi/(10000.0*ny) # in rad/m
	ky_units = 2.*np.pi/(10000.0*nx)
	omega_units = 2.*np.pi/(nt*dt*24*60*60) # in rad/day

	# Length of array dimensions after padding (calls on make_odd at top of file)
	Lx = make_odd(padding_fac*nx) #+ np.mod(range_nx,2)-1 # Make the padded domain size the same evenness/oddness as range_nx
	Ly = make_odd(padding_fac*ny) #+ np.mod(range_ny,2)-1
	T = make_odd(padding_fac*nt) # time-domain needs no padding
	print 'Lx = ',Lx

	# Add dimensions to variables to make them odd
	var1 = make_var_odd(var1)
	var2 = make_var_odd(var2)
	nx = var1.shape[0]
	ny = var1.shape[1]
	print 'var1.shape after making odd = ',var1.shape
	print 'new odd nx,ny = ',nx, ny

	# Construct the spectral domain
	kx = spec_domain(Lx)  # spec_domain is a function defined at top
	ky = spec_domain(Ly)
	#kt = spec_domain(T)
	print 'kx.shape = ',kx.shape

	# Make isotropic k: smallest of kx or ky
	kx_, ky_ = np.meshgrid(kx_units*kx,ky_units*ky)
	wv = kx_**2 + ky_**2
	maxk = math.sqrt(np.amax((wv)))
	mink = math.sqrt(min(wv[np.nonzero(wv>0)]))
	del kx_,ky_


	print 'maxk,mink=',maxk,mink

	# Isotropic wavenumber and frequency
	kiso = np.hstack((0,np.linspace(mink,maxk,kfac+1)))
	ktiso = omega_units*np.arange(0,np.floor(T/2)+1)	#range_kt[:len_kt/2-1] DOUBLE CHECK THIS LINE!!!!!
	print 'kiso.shape=',kiso.shape
	print 'ktiso.shape=',ktiso.shape
	del mink,maxk



	# Detrend relevant terms
	var1 = detrend_func.main(var1,spacetime)
	var2 = detrend_func.main(var2,spacetime)
	print 'var1.shape = ',var1.shape

	# Window relevant terms (J and p)
	var1 = window_func.main(var1,spacetime)
	var2 = window_func.main(var2,spacetime)

	# Take FFT
	if spacetime == 'spacetime':
		var1_fft = (1.0/(nx*ny*nt)) * np.fft.fftn(var1)
		var2_fft = (1.0/(nx*ny*nt)) * np.fft.fftn(var2)
		del var1,var2

	# Pad the variables
	if padding_fac > 1:
		var1_fft = pad(var1_fft,Lx,Ly,nx,ny)
		var2_fft = pad(var2_fft,Lx,Ly,nx,ny)

	# Calculate the isotropic transfer
	transfer_iso = make_iso(var1_fft,var2_fft,wv,kiso)
	print 'transfer_iso.shape = ',transfer_iso.shape









