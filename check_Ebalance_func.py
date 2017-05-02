# This code takes in spectral transfer data and creates a KE energy budget plot
#
# This code is run from 'run_TKE_TPE_Ebudget.py'
#-----------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def check_Ebalance(TKE1,TKE2,TKE3,TPE12,TPE23,windstress,buoyancy,bottomDrag,kiso,ktiso,terms_dict):

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

	# Add all terms together
	TotalE = TKE1 + TKE2 + TKE3 + TPE12 + TPE23 + windstress + buoyancy + bottomDrag
	
	if terms_dict.get('print_stuff'):
		print 'TotalE.shape = ',TotalE.shape
		print 'TotalE.shape[0] = ',TotalE.shape[0]

	cum_avg = np.zeros(TotalE.shape[0])
	for i in np.arange(1,TotalE.shape[0]):
		cum_avg[i] = np.sum(TotalE[:i])/i 

	if terms_dict.get('print_stuff'):
		print 'cum_avg.shape = ',cum_avg.shape

	# Create figure
	plt.figure(figsize=(12,8))
	plt.plot(cum_avg)

	if savefigs:
		plt.savefig(figpath+'/'+terms_dict.get('fig_savename')+'_gauss_smooth'+str(gauss_smooth)+'area_preserv.png')

	print 'check_Ebudget is done!'
	plt.show()


