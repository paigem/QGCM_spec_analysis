import numpy as np
import math

import detrend_func_longTime
import window_func_longTime


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
	nt = var.shape[0]
	
	if not np.mod(nt,2):
		var = np.append(var,[0])
	del nt

	return var


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
	field = make_wiso(var1,var2)	
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

	del wv,field,e0,en
	return out


# Make frequencies isotropic (i.e. symmetric)
def make_wiso(var1,var2):
	endw = float(var1.shape[0])

	# Define first and last entries of isotropic transfer
	out1 = (10**9)*np.real(np.conjugate(var1[0])*var2[0])
	out_end = (10**9)*np.real(np.conjugate(var1[0.5*endw])*var2[0.5*endw])
	
	# Define all other entries of isotropic transfer
	out = (10**9)*np.real(np.conjugate(var1[1:0.5*endw])*var2[1:0.5*endw] + np.conjugate(var1[-1:0.5*endw+1:-1])*var2[-1:0.5*endw+1:-1])

	# Concatenate the arrays together
	out_final = np.hstack((out1,out,out_end))

	return out_final

#----------------------------------------------------------------------------------------------------

def main(var1,var2,terms_dict):


	# Define FFT normalization
	nt = var1.shape[0]
	dt = terms_dict.get('dt') # in days
	padding_fac = terms_dict.get('padding_fac')

	# Normalization factor for wavenumbers (translating from lat/lon to wavenumber)
	omega_units = 2.*np.pi/(nt*dt*24*60*60) # in rad/day

	# Length of array dimensions after padding (calls on make_odd at top of file)
	T = make_odd(padding_fac*nt) # time-domain needs no padding

	# Add dimensions to variables to make them odd
	var1 = make_var_odd(var1)
	var2 = make_var_odd(var2)
	nt = var1.shape[0]
	if terms_dict.get('print_stuff'):
		print 'var1.shape after making odd = ',var1.shape

	# Isotropic frequency
	ktiso = omega_units*np.arange(0,np.floor(nt/2)+1)
	if terms_dict.get('print_stuff'):
		print 'ktiso.shape=',ktiso.shape

	# Detrend relevant terms
	var1 = detrend_func_longTime.main(var1,terms_dict.get('spacetime'),terms_dict)
	var2 = detrend_func_longTime.main(var2,terms_dict.get('spacetime'),terms_dict)
	if terms_dict.get('print_stuff'):
		print 'var1.shape = ',var1.shape

	# Window relevant terms (J and p)
	var1 = window_func_longTime.main(var1,terms_dict.get('spacetime'),terms_dict)
	var2 = window_func_longTime.main(var2,terms_dict.get('spacetime'),terms_dict)

	# Take FFT
	var1_fft = (1.0/(nt)) * np.fft.fft(var1)
	var2_fft = (1.0/(nt)) * np.fft.fft(var2)
	del var1,var2

	# Pad the variables
	if terms_dict.get('padding_fac') > 1:
		var1_fft = pad(var1_fft,Lx,Ly,nx,ny)
		var2_fft = pad(var2_fft,Lx,Ly,nx,ny)

	# Calculate the isotropic transfer
	transfer_iso = make_wiso(var1_fft,var2_fft)
	if terms_dict.get('print_stuff'):
		print 'transfer_iso.shape = ',transfer_iso.shape
	del var1_fft,var2_fft

	return transfer_iso,ktiso









