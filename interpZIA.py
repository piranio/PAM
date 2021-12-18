
#interpolate zia models to aeronet wavelengths

import csv
import json
import requests
import os
import pandas as pd
import numpy as np
from datetime import datetime as dt 
import sys
from more_itertools import unique_everseen
import math
import glob
import re
from scipy import interpolate
import matplotlib.pyplot as plt


##*************************************
##Author: Martin A. Montes
##SSAI, April 29, 2020
##Current version derived from interpZIA3.py
##*************************************

def interp_ref(fullname):
	df = pd.read_csv(fullname,header=None,error_bad_lines=False) 
	df=df.values
	
	am=df[:,0]
	rm=df[:,2]
	sigma=df[:,3]
	nr=df[:,4]
	ni=df[:,5]
	ndf=df[:,6]
	
	rm=rm.reshape(len(rm),1)
	sigma=sigma.reshape(len(sigma),1)
	ndf=ndf.reshape(len(ndf),1)
	
	nr=nr.reshape(len(nr),1)
	ni=ni.reshape(len(ni),1)
	am=am.reshape(len(am),1)
	
	kmax=80
	ix1=0
	ix2=34
	allmodels=np.empty([34,7])
	for k in range(1,kmax+1):
		print('processing model...',k)
		
		am_one=am[ix1:ix2] #34 rows=17 lambdas x 2 psd modes
		rm_one=rm[ix1:ix2]
		sigma_one=sigma[ix1:ix2]
		nr_one=nr[ix1:ix2]
		ni_one=ni[ix1:ix2]
		ndf_one=ndf[ix1:ix2]
	
	#separate in fine and coarse modes---------17 values for nr and ni each
		nr_f=nr_one[0::2]
		nr_c=nr_one[1::2]
	
		ni_f=ni_one[0::2].T
		ni_c=ni_one[1::2].T
			
	
	###-FINE MODE-############## WE INTERPOLATE TWICE FINE/COARSE MODES #################### FOR EACH ZIA MODEL	
		wl=[0.412,0.443,0.469,0.488,0.531,0.547,0.551,0.555,0.645,0.667,0.678,0.748,0.859,0.869,1.240,1.640,2.130]
		wi=[0.440,0.469,0.488,0.531,0.547,0.551,0.555,0.645,0.667,0.675,0.678,0.748,0.859,0.870,0.98,0.99,1.020]
		
	#####CURBIC INTERPOLATION AT AERONET OR 6S WAVELENGTHS----------------------------------------------	
	#beyond 1020 nm the interpolation as expected is bad!!
		s = interpolate.InterpolatedUnivariateSpline(wl, nr_f) #numpy array
		nri=s(wi)
		s = interpolate.InterpolatedUnivariateSpline(wl, ni_f) #numpy array
		nii=s(wi)
	
		nri=nri.reshape(len(nri),1) #17x1
		nii=nii.reshape(len(nii),1)
	
	#######-COARSE MODE-##################
		wl=[0.412,0.443,0.469,0.488,0.531,0.547,0.551,0.555,0.645,0.667,0.678,0.748,0.859,0.869,1.240,1.640,2.130]
		wi=[0.440,0.469,0.488,0.531,0.547,0.551,0.555,0.645,0.667,0.675,0.678,0.748,0.859,0.870,0.98,0.99,1.020]
		
	#####CURBIC INTERPOLATION coarse mode (in theory the same for the current inversion)----------------------------------------------

		s = interpolate.InterpolatedUnivariateSpline(wl, nr_c) #numpy array
		nri2=s(wi)
		s = interpolate.InterpolatedUnivariateSpline(wl, ni_c) #numpy array
		nii2=s(wi)
	

		nri2=nri2.reshape(len(nri2),1) #17x1
		nii2=nii2.reshape(len(nii2),1)
	
	###overriding values#########
		nr_one[0::2]=nri
		nr_one[1::2]=nri2
		ni_one[0::2]=nii
		ni_one[1::2]=nii2
	#print(ni_one)
	
	#merge as a 34 x 1 vector------------
	######################################
		wi=np.asarray(wi)
		wi=wi.reshape(len(wi),1)
		wlam_one=np.zeros([34,1])
		wlam_one[0::2]=wi
		wlam_one[1::2]=wi
	### merge output file for Mie simulations #########
	#concatenate interpolated variables plus those without like rm_n and sigma
		
		onemodel=np.empty([34,1])
		onemodel=np.concatenate((onemodel,am_one),axis=1)
		onemodel=np.concatenate((onemodel,wlam_one),axis=1)
		onemodel=np.concatenate((onemodel,rm_one),axis=1)
		onemodel=np.concatenate((onemodel,sigma_one),axis=1)
		onemodel=np.concatenate((onemodel,nr_one),axis=1)
		onemodel=np.concatenate((onemodel,ni_one),axis=1)
		onemodel=np.concatenate((onemodel,ndf_one),axis=1)
	
		onemodel=onemodel[:,1:]
	
		ix1=ix2
		ix2=ix2+34
		
		allmodels=np.concatenate((allmodels,onemodel),axis=0)
		
	
	
	
	allmodels=allmodels[34:,:]

	###re-arrangement of rows for Mie input
	m=1
	ixi=m-1
	ixe=ixi+10
	allmodels2=np.zeros([2720,7])
	m=0 #model index
	c=0
	
	
	for lam in range(0,34,2): #lambda index for each wavelength	
		#2720 x7
		for m in range(0,80): #80 models!!
			ix=lam+m*34
			allmodels2[c:c+2,:]=allmodels[ix:ix+2,:]
			
			c+=2
			
	
	allmodels2=pd.DataFrame(data=allmodels2)
	allmodels2.to_csv(path_out+'zia80interpMieInp.csv', header=False,index=False) 
	
	
	
	
#**********MAIN******************************	

if __name__ == '__main__':

	
	fpath='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/INTERP/in/' #other aerosol properties, *.all ext		
	path_out='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/INTERP/out/'	
	filename='zia80inp.csv'	#17 wavelengths
	
	fullname=fpath+filename

	interp_ref(fullname)