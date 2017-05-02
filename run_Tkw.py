#--------------------------------------------------------------------------------------------------------------------------------------------------
# This code runs the TKE and TPE codes for each of the layers, and then calls on the energy_budget_transfers.py code to plot the energy budget
#
#--------------------------------------------------------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os.path
from datetime import datetime

import TKE_func_new_memtest
import TPE_func
import windstress_func
import buoyancy_func
import bottom_drag_func
import plot_func
import check_Ebalance_func

# Specify parameters
fluid_var = 'oc-coupled'
climate_var = 'p'
domain = 'fulldomain'
yrs = [1,1]
output_num = '013'
run_name = 'dg2'

if domain == 'fulldomain':
	y1 = 0
	y2 = 960
	x1 = 0
	x2 = 960
	datatitle = 'fulldomain'
else:
	y1 = 0
	y2 = 100
	x1 = 0
	x2 = 100

print_stuff = 1
save_data = 1
spatial_flag = 1
padding_fac = 1.0
kfac = 100 # number of wavenumbers desired
max_layer = 3

# Plotting flags
plot_flag = 1
include_all_terms = 1
area_preserv = 1
gauss_smooth = 0
savefigs = 0
check_Ebalance = 1

if spatial_flag:
	spacetime = 'time'
	spatial_name = '_spatial'
else:
	spacetime = 'spacetime'
	spatial_name = ''

# Define constants to be used
dx = 5000. # meters
dy = 5000. # meters
dt = 1 # in days
H1 = 350.0 # meters
H2 = 750.0 # meters
H3 = 2900.0 # meters
Htot = H1 + H2 + H3
f0 = 9.37456*(10**(-5)) #1/s (Coriolis parameter)
g1 = .015
g2 = .0075

datapath = '/g/data/v45/pm2987/nco_and_output/oc_spunup/'
#datapath = '/g/data/v45/pm2987/nco_and_output/ocean/'
transfer_datapath = '/g/data/v45/pm2987/netcdf_transfers/'
figpath = '~/Documents/Python/Figures/'
#save_name = '_output'+output_num
save_name = ''
extra_name = '_recalc'

if include_all_terms:
	allTerms = 'allTerms'
else:
	allTerms = 'KEandPE'
fig_savename = 'Tbudget'+spatial_name+save_name+extra_name+'_'+allTerms+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])

# Try constructing a dict of all these single quantities
terms_dict = {'spacetime':spacetime,'padding_fac':padding_fac,'kfac':kfac,'print_stuff':print_stuff,'save_data':save_data,'savefigs':savefigs,'figpath':figpath,'fig_savename':fig_savename,'area_preserv':area_preserv,'gauss_smooth':gauss_smooth,'include_all_terms':include_all_terms,'save_name':save_name,'extra_name':extra_name,'spatial_flag':spatial_flag,'dx':[dx,dy],'dt':dt,'H':[H1,H2,H3,Htot],'f0':f0,'gprimes':[g1,g2],'yrs':yrs,'domain':[x1,x2,y1,y2]}


#--------------------------------------------------------------------------------------------------------------------------------------------------
print 'Start time ',datetime.now().time()

for layer in np.arange(1,max_layer+1):
	# Call TKE_func_new.py
	if not os.path.exists(transfer_datapath+'TKE'+spatial_name+save_name+extra_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'):
		dataname = climate_var+save_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
		#print 'dataname = ',dataname
		#dataname = 'p_Daily_1yr_dg2_output037_layer1_0_100_0_100_159_159_50daytest.nc'
		TKE_func_new_memtest.main(datapath,dataname,layer,terms_dict)
	else:
		print 'TKE code of layer '+str(layer)+' already exists!'

print 'End TKE time ',datetime.now().time()
#--------------------------------------------------------------------------
# Call TPE_func.py
for layer in np.arange(1,max_layer):
	if not os.path.exists(transfer_datapath+'TPE'+spatial_name+save_name+extra_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'):
		dataname1 = climate_var+save_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
		dataname2 = climate_var+save_name+'_layer'+str(layer+1)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
		#dataname1 = 'p_Daily_1yr_dg2_output037_layer1_0_100_0_100_159_159_50daytest.nc'
		#dataname2 = 'p_Daily_1yr_dg2_output037_layer2_0_100_0_100_159_159_50daytest.nc'
		TPE_func.main(datapath,dataname1,dataname2,layer,terms_dict)

	else:
		print 'TPE code of layer '+str(layer)+' already exists!'

print 'End TPE time ',datetime.now().time()

#--------------------------------------------------------------------------
# Call windstress_func.py
if not os.path.exists(transfer_datapath+'windstress'+spatial_name+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'):
	dataname1 = 'p'+save_name+'_layer1_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
	dataname2 = 'wekt'+save_name+'_'+str(x1)+'_'+str(x2-1)+'_'+str(y1)+'_'+str(y2-1)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
	#dataname1 = 'p_Daily_1yr_dg2_output037_layer1_0_100_0_100_159_159_50daytest.nc'
	#dataname2 = 'wekt_Daily_1yr_dg2_output037_0_100_0_100_159_159_50daytest.nc'
	windstress_func.main(datapath,dataname1,dataname2,terms_dict)

else:
	print 'windstress code already exists!'

print 'End windstress time ',datetime.now().time()

#--------------------------------------------------------------------------
# Call buoyancy_func.py
if not os.path.exists(transfer_datapath+'buoyancy'+spatial_name+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'):
	dataname1 = 'p'+save_name+'_layer1_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
	dataname2 = 'p'+save_name+'_layer2_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
	e_dataname = 'e'+save_name+'_layer1_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
	#dataname1 = 'p_Daily_1yr_dg2_output037_layer1_0_100_0_100_159_159_50daytest.nc'
	#dataname2 = 'p_Daily_1yr_dg2_output037_layer2_0_100_0_100_159_159_50daytest.nc'
	#e_dataname = 'e_Daily_1yr_dg2_output037_layer2_0_100_0_100_159_159_50daytest.nc'
	buoyancy_func.main(datapath,dataname1,dataname2,e_dataname,terms_dict)

else:
	print 'buoyancy code already exists!'

print 'End buoyancy time ',datetime.now().time()

#--------------------------------------------------------------------------
# Call bottom_drag_func.py
if not os.path.exists(transfer_datapath+'bottomDrag'+spatial_name+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'):
	dataname = 'p'+save_name+'_layer3_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
	#dataname = 'p_Daily_1yr_dg2_output037_layer2_0_100_0_100_159_159_50daytest.nc'
	bottom_drag_func.main(datapath,dataname,terms_dict)

else:
	print 'bottom drag code already exists!'

print 'End bottom drag time ',datetime.now().time()
#--------------------------------------------------------------------------

# Plot

if plot_flag:
	# Cycle through layers to pull up all three TKEs
	if terms_dict.get('spatial_flag'):
		transfer_name = 'TKE_spatial'
	else:
		transfer_name = 'TKE'

	for layer in np.arange(1,max_layer+1):
		KE_data = Dataset(transfer_datapath+transfer_name+save_name+extra_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
		T = KE_data.variables['T'] # shape is (x,y,time)

		# Average over spatial dimensions
		T = np.mean(np.mean(T,axis=0),axis=0)

		# Trick to name TKE variable names according to the layer
		temp = 'TKE'+str(layer)
		vars()[temp] = T

		# Define kiso and ktiso
		if not terms_dict.get('spatial_flag'):
			kiso = KE_data.variables['kiso']
		else:
			kiso = 0
		ktiso = KE_data.variables['ktiso']
		print 'in TKE',ktiso.shape, TKE1.shape
		print 'ktiso[6]=',ktiso[6]
		del KE_data

	# Cycle through layers to pull up all two TPEs
	if terms_dict.get('spatial_flag'):
		transfer_name = 'TPE_spatial'
	else:
		transfer_name = 'TPE'

	for layer in np.arange(1,3):
		PE_data = Dataset(transfer_datapath+transfer_name+save_name+extra_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
		T = PE_data.variables['T']
	
		# Average over spatial dimensions
		T = np.mean(np.mean(T,axis=0),axis=0)

		# Trick to name TKE variable names according to the layer
		temp = 'TPE'+str(layer)+str(layer+1)
		vars()[temp] = T
		del PE_data

	# Load windstress
	if terms_dict.get('spatial_flag'):
		transfer_name = 'windstress_spatial'
	else:
		transfer_name = 'windstress'
	print transfer_datapath+transfer_name+save_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc'
	wind_data = Dataset(transfer_datapath+transfer_name+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
	windstress = wind_data.variables['T']
	windstress = np.mean(np.mean(windstress,axis=0),axis=0)
	del wind_data

	# Load buoyancy
	if terms_dict.get('spatial_flag'):
		transfer_name = 'buoyancy_spatial'
	else:
		transfer_name = 'buoyancy'
	buoy_data = Dataset(transfer_datapath+transfer_name+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
	buoyancy = buoy_data.variables['T']
	buoyancy = np.mean(np.mean(buoyancy,axis=0),axis=0)
	del buoy_data

	# Load bottom drag
	if terms_dict.get('spatial_flag'):
		transfer_name = 'bottomDrag_spatial'
	else:
		transfer_name = 'bottomDrag'
	bd_data = Dataset(transfer_datapath+transfer_name+save_name+extra_name+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
	bottomDrag = bd_data.variables['T']
	bottomDrag = np.mean(np.mean(bottomDrag,axis=0),axis=0)
	del bd_data

	plot_func.plot_Ebudget(TKE1,TKE2,TKE3,TPE12,TPE23,windstress,buoyancy,bottomDrag,kiso,ktiso,terms_dict)

	
	if check_Ebalance:
		check_Ebalance_func.check_Ebalance(TKE1,TKE2,TKE3,TPE12,TPE23,windstress,buoyancy,bottomDrag,kiso,ktiso,terms_dict)



'''
# Plot

plot_datapath = '/g/data/v45/pm2987/netcdf_transfers/'
layer = 1

KE_data = Dataset(plot_datapath+'TKE_spatial_'+save_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
TKE_spatial = KE_data.variables['TKE']
del KE_data

# Average over spatial dimensions
TKE_spatial = np.mean(np.mean(TKE_spatial,axis=0),axis=0)

PE_data = Dataset(plot_datapath+'TPE_spatial_'+save_name+'_layer'+str(layer)+'_'+str(x1)+'_'+str(x2)+'_'+str(y1)+'_'+str(y2)+'_'+str(yrs[0])+'_'+str(yrs[1])+'.nc')
TPE_spatial = PE_data.variables['TKE']

# Define kiso and ktiso
ktiso = PE_data.variables['ktiso']
dw = ktiso[-1] - ktiso[-2]
del PE_data

# Average over spatial dimensions
TPE_spatial = np.mean(np.mean(TPE_spatial,axis=0),axis=0)
		
# Trick to name TKE variable names according to the layer
#temp = 'TPE'+str(layer)+str(layer+1)
#vars()[temp] = a.variables['TPE']

dw = ktiso[-1] - ktiso[-2]

if print_stuff:
	print 'ktiso.shape=',ktiso.shape
	print 'TKE.shape = ', TPE12.shape
	print 'TKE_spatial.shape = ',TKE1.shape
	print 'windstress',windstress.shape
	print 'buoyancy',buoyancy.shape

dw = ktiso[-1] - ktiso[-2]

# Create figure 1
plt.figure(num=1, figsize=(12,8))
plt.plot(ktiso,windstress/dw,linewidth=5.0,color='g',label='wind')
plt.plot(ktiso,buoyancy/dw,linewidth=5.0,color='y',label='buoyancy')
plt.axhline(0,color='k',linestyle='dotted',linewidth=3.0)
plt.xscale('log')
plt.axis('tight')
plt.title('Energy terms (calculated spatially)')
plt.xlabel('Frequency')
plt.ylabel('(nW/kg)/(rad/day)')
plt.legend()

plt.show()

'''
print 'End of code!'

