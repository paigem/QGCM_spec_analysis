# This code takes in spectral transfer data and creates a KE energy budget plot
#
# This code is run from 'run_TKE_TPE_Ebudget.py'
#-----------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def plot_Ebudget(TKE1,TKE2,TKE3,TPE12,TPE23,windstress,buoyancy,bottomDrag,kiso,ktiso,terms_dict): #figpath,climate_var,save_name,

	# Define dk and dw
	if not terms_dict.get('spatial_flag'):
		dk = kiso[-1] - kiso[-2]
	dw = ktiso[-5] - ktiso[-6]
	print 'dw = ',dw

	# For saving Tdiss
	diss_savename = ''

	# Define from terms_dict
	gauss_smooth = terms_dict.get('gauss_smooth')
	area_preserv = terms_dict.get('area_preserv')
	savefigs = terms_dict.get('savefigs')
	figpath = terms_dict.get('figpath')

	# Specify smaller font size for screen viewing
	font = {'family' : 'normal',
	        'size'   : 10}
	matplotlib.rc('font', **font)

	# Filler for now
	x = 1

	if terms_dict.get('include_all_terms'):

		# Frequency plots
		if gauss_smooth and area_preserv:
	
			print 'shapes=',ktiso.shape,dk.shape,TKE1.shape
	
			# Create figure 1
			plt.figure(num=1, figsize=(12,8))
			# KEs
			plt.plot(ktiso,gaussian_filter1d(TKE1*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='m',label='KE1')
			plt.plot(ktiso,gaussian_filter1d(TKE2*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='m',linestyle='dashed',label='KE2')
			plt.plot(ktiso,gaussian_filter1d(TKE3*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='m',linestyle='dotted',label='KE3')
			# PEs
			plt.plot(ktiso,gaussian_filter1d(TPE12*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='b',label='PE12')
			plt.plot(ktiso,gaussian_filter1d(TPE23*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='b',linestyle='dashed',label='PE23')
			plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
			# Windstress/Buoyancy/Bottom Drag
			plt.plot(ktiso,gaussian_filter1d(windstress*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='c',label='windstress')
			plt.plot(ktiso,gaussian_filter1d(buoyancy*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='g',label='buoyancy')
			plt.plot(ktiso,gaussian_filter1d(bottomDrag*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='y',label='bottom drag')
			#if include_dissipation:
			#	plt.plot(ktiso,gaussian_filter1d(np.sum(Tdiss1,axis=0)*dk*ktiso,sigma=gauss_smooth),linewidth=2.0,color='#778899',label='Tdiss1')
			#	plt.plot(ktiso,gaussian_filter1d(np.sum(Tdiss2,axis=0)*dk*ktiso,sigma=gauss_smooth),linewidth=2.0,color='#778899',linestyle='dashed',label='Tdiss2')
			#	plt.plot(ktiso,gaussian_filter1d(np.sum(Tdiss3,axis=0)*dk*ktiso,sigma=gauss_smooth),linewidth=2.0,color='#778899',linestyle='dotted',label='Tdiss3')
			#	diss_savename = 'withTdiss_'

			plt.axvline(2*np.pi/500,color='k',linestyle='dashed',linewidth=3.0) # it takes a rossby wave 908 days to cross the basin (I think)
			plt.xscale('log')
			plt.axis('tight')
			plt.title('Energy budget TKE(w), TPE(w), TBe(w)')
			plt.xlabel('Frequency (rad/day)')
			plt.ylabel('(nW/kg)/(rad/km)')
			plt.legend()
			if savefigs:
				plt.savefig(figpath+'/'+terms_dict.get('fig_savename')+'_gauss_smooth'+str(gauss_smooth)+'area_preserv.png')
		elif area_preserv:
			# Create figure 1
			plt.figure(num=6, figsize=(15,12))
			# KEs
			plt.plot(ktiso,TKE1*ktiso/dw,linewidth=5.0,color='m',label='KE1')
			plt.plot(ktiso,TKE2*ktiso/dw,linewidth=5.0,color='m',linestyle='dashed',label='KE2')
			plt.plot(ktiso,TKE3*ktiso/dw,linewidth=5.0,color='m',linestyle='dotted',label='KE3')
			# PEs
			plt.plot(ktiso,TPE12*ktiso/dw,linewidth=5.0,color='b',label='PE12')
			plt.plot(ktiso,TPE23*ktiso/dw,linewidth=5.0,color='b',linestyle='dashed',label='PE23')
			# Windstress/Buoyancy/Bottom Drag
			plt.plot(ktiso,windstress*ktiso/dw,linewidth=5.0,color='c',label='windstress')
			plt.plot(ktiso,buoyancy*ktiso/dw,linewidth=5.0,color='g',label='buoyancy')
			plt.plot(ktiso,bottomDrag*ktiso/dw,linewidth=5.0,color='y',label='bottom drag')


			plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
			plt.xscale('log')
			plt.axis('tight')
			plt.title('Spectral transfer energy budget')
			plt.xlabel('Frequency')
			plt.ylabel('(nW/kg)/(rad/km)')
			plt.legend()
			if savefigs:	
				plt.savefig(figpath+'/'+terms_dict.get('fig_savename')+'_area_preserv.png')
	
		elif gauss_smooth:
			# Create figure 2
			plt.figure(num=6, figsize=(15,12))
			# KEs
			# KEs
			plt.plot(ktiso,gaussian_filter1d(TKE1/dw,sigma=gauss_smooth),linewidth=5.0,color='m',label='KE1')
			plt.plot(ktiso,gaussian_filter1d(TKE2/dw,sigma=gauss_smooth),linewidth=5.0,color='m',linestyle='dashed',label='KE2')
			plt.plot(ktiso,gaussian_filter1d(TKE3/dw,sigma=gauss_smooth),linewidth=5.0,color='m',linestyle='dotted',label='KE3')
			# PEs
			plt.plot(ktiso,gaussian_filter1d(TPE12/dw,sigma=gauss_smooth),linewidth=5.0,color='b',label='PE12')
			plt.plot(ktiso,gaussian_filter1d(TPE23/dw,sigma=gauss_smooth),linewidth=5.0,color='b',linestyle='dashed',label='PE23')
			plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
			# Windstress/Buoyancy/Bottom Drag
			plt.plot(ktiso,gaussian_filter1d(windstress/dw,sigma=gauss_smooth),linewidth=5.0,color='c',label='windstress')
			plt.plot(ktiso,gaussian_filter1d(buoyancy/dw,sigma=gauss_smooth),linewidth=5.0,color='g',label='buoyancy')
			plt.plot(ktiso,gaussian_filter1d(bottomDrag/dw,sigma=gauss_smooth),linewidth=5.0,color='y',label='bottom drag')

			plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
			plt.xscale('log')
			plt.axis('tight')
			plt.title('Energy budget TKE(w), TPE(w), TBe(w)')
			plt.xlabel('Frequency')
			plt.ylabel('(nW/kg)/(rad/day)')
			plt.legend()
			if savefigs:
				plt.savefig(figpath+'/'+terms_dict.get('fig_savename')+'_gauss_smooth'+str(gauss_smooth)+'.png')

	'''
	# Make plots with just TKE and TPE
	else:
		# Frequency plots of PE and KE only
		if gauss_smooth and area_preserv:
	
			# Create figure 1
			plt.figure(num=6, figsize=(15,12))
			# KEs
			plt.plot(ktiso,gaussian_filter1d(np.sum(TKE1,axis=0)*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='m',label='KE1')
			plt.plot(ktiso,gaussian_filter1d(np.sum(TKE2,axis=0)*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='m',linestyle='dashed',label='KE2')
			plt.plot(ktiso,gaussian_filter1d(np.sum(TKE3,axis=0)*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='m',linestyle='dotted',label='KE3')
			# PEs
			plt.plot(ktiso,gaussian_filter1d(np.sum(TPE12,axis=0)*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='b',label='PE12')
			plt.plot(ktiso,gaussian_filter1d(np.sum(TPE23,axis=0)*ktiso/dw,sigma=gauss_smooth),linewidth=5.0,color='b',linestyle='dashed',label='PE23')
			plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
			plt.axvline(2*np.pi/500,color='k',linestyle='dashed',linewidth=3.0) # it takes a rossby wave 908 days to cross the basin (I think)
			plt.xscale('log')
			plt.axis('tight')
			plt.title('Energy budget TKE(w) and TPE(w)')
			plt.xlabel('Frequency (rad/day)')
			plt.ylabel('(nW/kg)/(rad/km)')
			plt.legend()
			if savefigs:
				plt.savefig(figpath+'/'+terms_dict.get('fig_savename')+'_gauss_smooth'+str(gauss_smooth)+'area_preserv.png')

		elif area_preserv:
			# Create figure 1
			plt.figure(num=6, figsize=(15,12))
			# KEs
			plt.plot(ktiso,TKE1*ktiso/dw,linewidth=5.0,color='m',label='KE1')
			plt.plot(ktiso,TKE2*ktiso/dw,linewidth=5.0,color='m',linestyle='dashed',label='KE2')
			plt.plot(ktiso,TKE3*ktiso/dw,linewidth=5.0,color='m',linestyle='dotted',label='KE3')
			# PEs
			plt.plot(ktiso,TPE12*ktiso/dw,linewidth=5.0,color='b',label='PE12')
			plt.plot(ktiso,TPE23*ktiso/dw,linewidth=5.0,color='b',linestyle='dashed',label='PE23')
			plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
			plt.xscale('log')
			plt.axis('tight')
			plt.title('Energy budget TKE(w) and TPE(w)')
			plt.xlabel('Frequency')
			plt.ylabel('(nW/kg)/(rad/km)')
			plt.legend()
			if savefigs:
				plt.savefig(figpath+'/'+terms_dict.get('fig_savename')+'_area_preserv.png')
	
		elif gauss_smooth:
			# Create figure 2
			plt.figure(num=6, figsize=(15,12))
			# KEs
			plt.plot(ktiso,gaussian_filter1d(np.sum(TKE1,axis=0)*dk,sigma=gauss_smooth),linewidth=5.0,color='m',label='KE1')
			plt.plot(ktiso,gaussian_filter1d(np.sum(TKE2,axis=0)*dk,sigma=gauss_smooth),linewidth=5.0,color='m',linestyle='dashed',label='KE2')
			plt.plot(ktiso,gaussian_filter1d(np.sum(TKE3,axis=0)*dk,sigma=gauss_smooth),linewidth=5.0,color='m',linestyle='dotted',label='KE3')
			# PEs
			plt.plot(ktiso,gaussian_filter1d(np.sum(TPE12,axis=0)*dk,sigma=gauss_smooth),linewidth=5.0,color='b',label='PE12')
			plt.plot(ktiso,gaussian_filter1d(np.sum(TPE23,axis=0)*dk,sigma=gauss_smooth),linewidth=5.0,color='b',linestyle='dashed',label='PE23')
			plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
			plt.xscale('log')
			plt.axis('tight')
			plt.title('Energy budget TKE(w) and TPE(w)')
			plt.xlabel('Frequency')
			plt.ylabel('(nW/kg)/(rad/day)')
			plt.legend()
			if savefigs:
				plt.savefig(figpath+'/'+terms_dict.get('fig_savename')+'_gauss_smooth'+str(gauss_smooth)+'.png')
	'''


	print 'plot_Ebudget is done!'
	plt.show()


