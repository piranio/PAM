import csv
import json
import requests
import os
import pandas as pd
import numpy as np
from datetime import datetime as dt 
import sys
import math
import glob
import PAMbinder as rbind

##*************************************
##Author: Martin A. Montes
##SSAI, April 29, 2020
#version derived from merge_csv5.py, call_Rpy.py and PAM4.R
##*************************************

	######################################################################
	#MAIN------------
	######################################################################	

if __name__ == '__main__':
	
	path_in='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/in/'	
	path_out='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/out/'	
	path_r_in='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/R/in/'
		
	####----Reading AERONET data including SSA values-------
	fullname=path_in+'globe_aeronet_SSArecordsID.csv'
	data2= pd.read_csv(fullname, header=0,error_bad_lines=False)
	data2=data2.values	
	row,col=data2.shape
	
	###--Split the data in even and uneven record ids-----
	odd=data2[0::2,:]
	even=data2[1::2,:]	
	even1=odd
	ix=np.round(len(even)/2)
	ix=int(ix)
	even2=even[0:ix,:]
	even3=even[ix:,:]
		
	even1=pd.DataFrame(even1) #around50% of dataset--
	even2=pd.DataFrame(even2)
	even3=pd.DataFrame(even3)
	
	####-create two subdatasets with 75% and 25% of data for clustering and validation, respectively
	data75=pd.concat([even1, even2],axis=0) #adding rows
	data25=even3
		
	data75.to_csv(path_out+'globe_aeronet_SSAtrain2'+'.csv', index=False)
	data25.to_csv(path_out+'globe_aeronet_SSAtest2'+'.csv', index=False)
	
	
	
	
	###-We start the clustering process with a the training dataset-------

	fullname=path_out+'globe_aeronet_SSAtrain2.csv'
	traindata= pd.read_csv(fullname, header=0,error_bad_lines=False)
	traindata=traindata.values #header off
	
	####-------Data Stratification-------Set aod(440) levels	
	aod440c=traindata[:,11]
	aod440t=traindata[:,15]
	
	####-Selecting clustering parameters-----	
	ssa440=traindata[:,19] #single scattering albedo at 440 nm
	#ssa870=traindata[:,22] #single scattering albedo at 870 nm
	AEe=traindata[:,37]   #angstrom coefficient for extinction for spectral range 440-870 nm
	AEa=traindata[:,39]    #angstrom coefficient for absorption for spectral range 440-870 nm
	
	aod440m=[]
	##We use mean values of aod440coincident and non-coincident for each record
	for i in range(0,len(aod440c)):
		temp=[aod440c[i],aod440t[i]]
		aod440m.append(np.mean(temp))
	
	####--Filtering aerosol parameter indices per aod level--------
	aod440m=np.asarray(aod440m)
	
	ix=np.where((aod440m>0.) & (aod440m<=0.4))
	
	
	
	data1=traindata[ix,:]
	aod1=aod440m[ix]
	#ssa870=ssa870[ix]
	ssa440=ssa440[ix]
	AEe=AEe[ix]
	AEa=AEa[ix]
	
	#concatenate and save-----
	ix=np.asarray(ix,dtype=int)
	ix=ix.T
	
	
	#ssa870=ssa870.reshape((len(ix),1))
	ssa440=ssa440.reshape((len(ix),1))
	AEe=AEe.reshape((len(ix),1))
	AEa=AEa.reshape((len(ix),1))
	
	dataclus=np.concatenate((ssa440,AEe),axis=1)
	dataclus=np.concatenate((dataclus,AEa),axis=1)
	#dataclus=np.concatenate((dataclus,ssa870),axis=1)
	
	
	#count how many clustering params
	nrec,nvar=dataclus.shape
	print('number of clustering variables...',nvar)
	
	dataclus=pd.DataFrame(dataclus)
	dataclus.to_csv(path_r_in+'trainglobe75_3params2'+'.csv', index=False,header=False)
	#dataclus.to_csv(path_r_in+'trainglobe75_4params2'+'.csv', index=False,header=False)
	
	
	data1=np.squeeze(data1)
	dataclus2=pd.DataFrame(data1)
	dataclus2.to_csv(path_r_in+'trainglobe75_3params_allprops2'+'.csv', index=False,header=False)
	
	
	
	
	
	####---Calling python-R binder------
	print('calling PAM clustering wrapper...')
	
	###---kmaxi=number of tentative clusters--initial guess
	kmaxi=30 #first estimate of maximum number of clusters
	flag_run=1 #1 = initial clustering, 2 = final clustering
	print('kmaxi-----initial number of max clusters..',kmaxi)
	
	##--PAM clustering
	rbind.pam_clus2(kmaxi,flag_run)

	nsamples=10 #number of records per cluster with a Silhouette threshold
	sil_thre=0.51 #0.51 # 0.71 #0.51 #0.26 #0.71
	
	
	####--Opening clustering results---medoids, silouette-and cluster ids----
	####--Determining number of clusters for each aod and clustering settings (nvars)
	kmax=rbind.check_kmax(kmaxi,sil_thre,nsamples)
	
	###---Maximum number of clusters-final estimate
	print('kmax..maximum number of clusters final estimate..',kmax)
	
	
	kmaxi=kmax
	flag_run=2 
	##---Optimization of Mahalanobis distances with the right number of clusters--------
	rbind.pam_clus2(kmaxi,flag_run)
	
	
	
	

	
	
	
	
	
	