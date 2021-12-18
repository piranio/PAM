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
import matplotlib.pyplot as plt

##*************************************
##Author: Martin A. Montes
##SSAI, April 29, 2020
##Current version derived from clus_stat5.py 
##*************************************

if __name__ == '__main__':
	#for each site and each cluster count number of records and nsites x 29 clusters array for each month
	#compute most frequent cluster per month, mean and std of aerosol properties

	path_in='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/out/'
	path_out='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/out'
	
	filename='train75withclus.csv'
	fullname=path_in+filename	
	#fields
	#sta	obs_id	cluster id	year	month	day	hour	min	sec	lat	lon	elev	aod1c	aod2c	aod3c	aod4c	aod1t	aod2t	aod3t	aod4t	ssa1	ssa2	ssa3	ssa4	aod1a	aod2a	aod3a	aod4a	g1	g2	g3	g4	volcf	vmrf	stdf	volcc	vmrc	stdc	AEe1	AEe2	AEa	nr1	nr2	nr3	nr4	ni1	ni2	ni3	ni4
	data= pd.read_csv(fullname, header=0,error_bad_lines=False)
	data=data.values
	row,col=data.shape
	month=data[:,4]
		
	#load sitenames file
	filename2='aeronetSitesGlobalTill92019.csv'
	#filename2='fake.csv'
	f2=path_in+filename2
	snames=pd.read_csv(f2, header=0,error_bad_lines=False)
	snames=snames.values
	sitenames=snames[:,0]
		
	si=np.shape(sitenames)
	#big array nsites x 29 clust for each month
	bigarray=np.zeros([si[0],12,29])
	#print(np.shape(jan))
	#sys.exit()
	
	for c,i in enumerate(sitenames):
		#find i matches
		matches=np.where(i==data[:,0])		
		print(c,i)
	
		onesta=data[matches,:]	
		nd=onesta.ndim
		mm=np.shape(onesta)
		if mm[0]==1 and nd==3 and mm[1]==1:
			onesta=onesta.reshape(1,49)
		else:
			onesta=np.squeeze(onesta)
		
		#for each month screen all clusters---
		for k in range(1,13):						
			match_month=np.where(k==onesta[:,4])
			
			if any(map(len, match_month)):
				onemonth=onesta[match_month,:]
				
				nd=onemonth.ndim
				mm=np.shape(onemonth)
				if mm[0]==1 and nd==3 and mm[1]==1:					
					onemonth=onemonth.reshape(1,49)
				else:
					onemonth=np.squeeze(onemonth)
									
				#for each cluster id ---compute frequencies mean and std and per month
				for j in range(1,30):
					match_clus=np.where(j==onemonth[:,2]) #check clus id
					if any(map(len, match_clus)):						
						mm=onemonth[match_clus,2]	
						row,col=np.shape(mm)						
						col=int(col)										
						bigarray[c][k-1][j-1]=bigarray[c][k-1][j-1]+col #site, month, clus
							
	
	nc_per_month=np.zeros([si[0],12]) #number of clusters per site
	cmaxn_per_month=np.zeros([si[0],12])#number of recrods for dominant cluster for that month
	cmax_per_month=np.zeros([si[0],12]) #dominant cluster for that month
	cmax_frac_per_month=np.zeros([si[0],12]) #dominant cluster fraction per month
	
	
	for i in range(0,si[0]):
		for k in range(0,12):
			mm=np.sum(bigarray[i,k,:]) #this is ok when we only we no repetitive clusters per month
			mm2=np.max(bigarray[i,k,:])								
			mm3=mm2/mm
						
			if mm2>0:#index method is for list!!
				ix=np.where(bigarray[i,k,:]==mm2)[0]
				#print(bigarray[i,k,:])
				print(i)
				print(k)
				print(ix)
				print(mm2)
				#if there is a problem when two clusters have the same N!! they are tied
				if len(ix)>1: #co-dominant clus
					ix=49
					
				cmax_per_month[i,k]=ix+1
				
			nc_per_month[i,k]=mm
			cmaxn_per_month[i,k]=mm2
			
			cmax_frac_per_month[i,k]=mm3
	
	#**********OUTPUT*********************
	cmax=pd.DataFrame(data=cmax_per_month) 		
	cmax.to_csv(path_out+'cmax_globalOK.csv', index=False) 
	cmaxn=pd.DataFrame(data=cmaxn_per_month) 		
	cmaxn.to_csv(path_out+'cmaxn_globalOK.csv', index=False) 
			
			
		
	
	