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
from scipy import interpolate
import math

##*************************************
##Author: Martin A. Montes
##SSAI, April 29, 2020
##Current version derived from make_fortran_inp18.py 
##*************************************



#**********************************************************************************
def nd_frac2(vo_r,mun,sigma_v): #vof and voc are volume concentration fractions----
	#Method to compute number number density fractions 
	#vo_r volumen concentration of aerosol particles (integrated across the log-Normal distribution)
	#mun is the ln of rn
	#sigma_v is the width of PSD = sigma_n, with v and n representing volume and number density space
	
	fudge_n=1000000 #this is an arbitrary number of particles--fudge factor-
	
	vo_f=vo_r[0::2] #fine mode Vo, these values are derived from 10 Cvf:Cvc cases suggested by A2010 for Mie simulations
	vo_c=vo_r[1::2] #coarse mode Vo
	
	mun_f=mun[0::2]
	mun_c=mun[1::2]
	
	sigma_f=sigma_v[0::2] #log log-normal distribution
	sigma_c=sigma_v[1::2]
			
	Cv_f=vo_f*fudge_n
	Cv_c=vo_c*fudge_n
	Cv_t=Cv_f+Cv_c
	
	
	##Calculations below based on 
	#Vo=(4/3)*pi No exp(3 mu +4.5 sigma^2)
	
	#---eq 77 Grainger, R. G.: Some Useful Formulae for Aerosol Size Distributions 
	#and Optical Properties, available at: http://eodg.atm.ox.ac.uk/user/grainger/research/aerosols.pdf (last access: 9 April 2020), 2017.
	#this is based on infinitum  integration plus/minus!!	
	
	expterm_f=np.exp(3*mun_f+4.5*np.power(sigma_f,2))
	No_f=(3./4.)*Cv_f*(1./math.pi)*(1./expterm_f)	
	expterm_c=np.exp(3*mun_c+4.5*np.power(sigma_c,2))
	No_c=(3./4.)*Cv_c*(1./math.pi)*(1./expterm_c)
	No_t=No_f+No_c
			
	fracn_f=No_f/No_t
	fracn_c=No_c/No_t
	return fracn_f,fracn_c

#**********************************************************************************
def rv_to_rm2(rv,sigma_v):
	#Method to convert rv from AERONET (volume density space) to rn of A2010 (number density space)
	
	muv=[math.log(i) for i in rv]	
	muv=np.asarray(muv)
				
	mun=muv-3*np.power(sigma_v,2) #from Grainger 2017 eq. 78
	
	rm_n=np.exp(mun)

	return rm_n,mun,muv




#**********************************************************************************
def inp_miefile():
	#Method for gluing routines and concatenating files for MODTRAN and FORTRAN mie input files

	#path for cluster ids derived from PAM-----
	path_in='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/R/out/'
	
	#path for training dataset
	path_in2='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/out/'	
		
	#path to keep intermediate Mie input files
	path_out2='/Users/mmontes/Documents/aeronet/globe/scripts/fortran/out/'
	
	
	filename1='clus_idinit.csv' #sil values
	fullname1=path_in+filename1
	
	filename2='trainglobe75_3params_allprops.csv' #all aero attributes
	fullname2=path_in2+filename2
	
	df1= pd.read_csv(fullname1, header=0,error_bad_lines=False)
	df1=df1.values
	df2= pd.read_csv(fullname2,header=None,error_bad_lines=False)
	df2=df2.values
	
	#merging labels with full database-------	
	df3=np.concatenate((df1,df2),axis=1)

	clus_lab=df3[:,1]#cluster id
	rf=df3[:,34] #median radius fine
	stdf=df3[:,35] #sigma
	rc=df3[:,37] #median radius coarse
	stdc=df3[:,38] #sigma coarse	
	nr=df3[:,42:46] #real refractive index
	ni=df3[:,46:50] #imaginary refractive index
		
	#interpolate AERONET wavelengths to A2010 wavelengths---I picked 17 rather than 170 lambdas from ZIA--
	wl=[0.44,0.675,0.87,1.020]
	wi=[0.440,0.469,0.488,0.531,0.547,0.551,0.555,0.645,0.667,0.675,0.678,0.748,0.859,0.870,0.98,0.99,1.020] #17 lambdas
	
	nri=[]
	nii=[]
	for i in nr:
		s = interpolate.InterpolatedUnivariateSpline(wl, i) #numpy array
		nri.append(s(wi))	
	for i in ni:
		s = interpolate.InterpolatedUnivariateSpline(wl, i) #numpy array
		nii.append(s(wi))	
	
	
	nri=np.asarray(nri)
	nii=np.asarray(nii)
	
	clus_lab=clus_lab.reshape((len(clus_lab),1))
	
	rf=rf.reshape((len(rf),1))
	stdf=stdf.reshape((len(stdf),1))			
	rc=rc.reshape((len(rc),1))
	stdc=stdc.reshape((len(stdc),1))

	#merging properties		
	df4=np.concatenate((clus_lab,rf),axis=1)
	df4=np.concatenate((df4,stdf),axis=1)
	df4=np.concatenate((df4,rc),axis=1)
	df4=np.concatenate((df4,stdc),axis=1)
	df4=np.concatenate((df4,nri),axis=1)
	df4=np.concatenate((df4,nii),axis=1)

	#############COMPUTING MEAN VALUES FOR AERO PROPERTIES##########################
	var=[]
	row,col=df4.shape
	allm=np.zeros([col,1])
	
	kmax=np.max(clus_lab) #number of clusters
	
	for i in range(1,kmax+1):	
		ix=np.where(df4[:,0]==i)
	
		for j in range(0,39): #number of aerosol parameters
			mm=np.mean(df4[ix,j])
			var.append(mm)
	
		var=np.asarray(var)
		var=var.reshape((len(var),1)) #going from 0-d to 1-d
		allm=np.concatenate((allm,var),axis=1) #clus1, clus2........clus29
		var=[]
	
	
	#Create input txt file for fortran runs----
	#################################################################################
	#**this is the structure************************
	#row1=lambda1 in microns clus1 fine mode	
	#row2=lambda1 clus1 fine coarse
	#row3=lambda1 clus2 fine
	#row4=lambda1 clus2 coarse
	#and so on.......
	#row x lambda2 clus1 fine
	#row x lambda 2 clus1 coarse
	#and so on....
	
	#create an 1-d array of ZIA lambdas 
	########################################################
	wi=np.asarray(wi)
	col1=wi.T
	col1=col1.reshape((len(col1),1)) #going from 0-d to 1-d
	
	#repeat each wavelength  = nclus x 2 modes = 29 x 2 = 58
	arr_lam=np.tile(wi,(kmax*2,1)) #repeat each lambda kmax*2 times rows
	
	#size parameters--------------------------------------------------
	sizeall=np.zeros([2,2])
	for k in range(1,kmax+1):
		fine=allm[1:3,k].T
		coarse=allm[3:5,k].T	
		size=[fine,coarse]
		sizeall=np.concatenate((sizeall,size),axis=0)
		
	#volume fraction values-----for fine and coarse modes-----
	vf=[0.00,0.01,0.02,0.05,0.1,0.20,0.30,0.50,0.80,0.95]	
	vfall=np.zeros([2,1])
	
	for k in range(1,kmax+1):
		vfone=[0.50,0.50] #default Cvf:Cfvc value will be updated later
		vfone=np.asarray(vfone)
		vfone=vfone.reshape(1,2)
		vfone=vfone.T #concatenating nr and ni for each cluster (1 wavelength)
		vfall=np.concatenate((vfall,vfone),axis=0)
	
	#integration intervals for computing density number fractions for fine and coarse modes----
	dxall=np.zeros([2,1])
	for k in range(1,kmax+1):
		dxone=[0.05,0.15]
		dxone=np.asarray(dxone)
		dxone=dxone.reshape(1,2)
		dxone=dxone.T #concatenating nr and ni for each cluster (1 wavelength)
		dxall=np.concatenate((dxall,dxone),axis=0)
	
	
	#refractive indices----------------------------------------------
	#same nr and ni for each psd mode-----
	clus1all=np.zeros([2,8]) #Mie JC input file-Geany, 2 rows for fine and coarse
	clus2all=np.zeros([2,8])
	clus3all=np.zeros([2,8])
	clus4all=np.zeros([2,8])
	clus5all=np.zeros([2,8])
	clus6all=np.zeros([2,8])
	clus7all=np.zeros([2,8])
	clus8all=np.zeros([2,8])
	clus9all=np.zeros([2,8])
	clus10all=np.zeros([2,8])
	clus11all=np.zeros([2,8])
	clus12all=np.zeros([2,8])
	clus13all=np.zeros([2,8])
	clus14all=np.zeros([2,8])
	clus15all=np.zeros([2,8])
	clus16all=np.zeros([2,8])
	clus17all=np.zeros([2,8])
	clus18all=np.zeros([2,8])
	clus19all=np.zeros([2,8])
	clus20all=np.zeros([2,8])
	clus21all=np.zeros([2,8])
	clus22all=np.zeros([2,8])
	clus23all=np.zeros([2,8])
	clus24all=np.zeros([2,8])
	clus25all=np.zeros([2,8])
	clus26all=np.zeros([2,8])
	clus27all=np.zeros([2,8])
	clus28all=np.zeros([2,8])
	clus29all=np.zeros([2,8])
	
	c6=np.zeros([58,1]) #29 PAM models x 2, aerosol model-cluster id or label starting with one
	for i in range(1,len(wi)+1): #for each lambda!
		refall=np.zeros([2,2])
		for k in range(1,kmax+1): #cluster id or name label is already ordered row-wise
			fine=[allm[4+i,k],allm[21+i,k]] #lambda1----
			fine=np.asarray(fine)
			fine=fine.T #concatenating nr and ni for each cluster (1 wavelength)
			coarse=fine		
			ref=[fine,coarse]
			refall=np.concatenate((refall,ref),axis=0)	
				
	#Creating concatenated arrays for each lambda	
		c1=arr_lam[:,i-1]	#lambdax 1-16
		c1=c1.reshape(len(c1),1) #adding 1 dimension!
		c2=refall[2:,:]	#refractive index real and imaginary
	
		wave=str(c1[0])
		c3=sizeall[2:,:] #median volumen radius	and sig=log(std)		
		rv=c3[:,0]
		sigma_v=c3[:,1]
		rv=rv.reshape(len(rv),1)
		sigma_v=sigma_v.reshape((len(sigma_v),1))
		c3=np.concatenate((rv,sigma_v),axis=1)
		
		no_f=0.99 #default values for NDF will be updated later
		no_c=0.01
		c4=np.zeros([len(rv),1])
		c4[0::2]=no_f #just a fudge number
		c4[1::2]=no_c
		c5=dxall[2:,:] #integration step for density number calculations
		
		aeromod=[]
		for i in range(1,30):
			aeromod.append(i)
		
		aeromod=np.asarray(aeromod)
		aeromod=aeromod.reshape(29,1)
		c6[0::2]=aeromod
		c6[1::2]=aeromod
	
		####merging variables for each wavelength ALL CLUSTERS##############################
		###################################################
		miefile=np.concatenate((c1,c3),axis=1) 
		miefile=np.concatenate((miefile,c2),axis=1)
		miefile=np.concatenate((miefile,c4),axis=1)
		miefile=np.concatenate((miefile,c5),axis=1)
		miefile=np.concatenate((miefile,c6),axis=1)
		

		#####################SPLIT FOR EACH CLUSTER#########################
		clus1all=np.concatenate((clus1all,miefile[0:2,:]),axis=0)
		clus2all=np.concatenate((clus2all,miefile[2:4,:]),axis=0)
		clus3all=np.concatenate((clus3all,miefile[4:6,:]),axis=0)
		clus4all=np.concatenate((clus4all,miefile[6:8,:]),axis=0)
		clus5all=np.concatenate((clus5all,miefile[8:10,:]),axis=0)
		clus6all=np.concatenate((clus6all,miefile[10:12,:]),axis=0)
		clus7all=np.concatenate((clus7all,miefile[12:14,:]),axis=0)
		clus8all=np.concatenate((clus8all,miefile[14:16,:]),axis=0)
		clus9all=np.concatenate((clus9all,miefile[16:18,:]),axis=0)
		clus10all=np.concatenate((clus10all,miefile[18:20,:]),axis=0)
		clus11all=np.concatenate((clus11all,miefile[20:22,:]),axis=0)
		clus12all=np.concatenate((clus12all,miefile[22:24,:]),axis=0)
		clus13all=np.concatenate((clus13all,miefile[24:26,:]),axis=0)
		clus14all=np.concatenate((clus14all,miefile[26:28,:]),axis=0)
		clus15all=np.concatenate((clus15all,miefile[28:30,:]),axis=0)
		clus16all=np.concatenate((clus16all,miefile[30:32,:]),axis=0)
		clus17all=np.concatenate((clus17all,miefile[32:34,:]),axis=0)
		clus18all=np.concatenate((clus18all,miefile[34:36,:]),axis=0)
		clus19all=np.concatenate((clus19all,miefile[36:38,:]),axis=0)
		clus20all=np.concatenate((clus20all,miefile[38:40,:]),axis=0)
		clus21all=np.concatenate((clus21all,miefile[40:42,:]),axis=0)
		clus22all=np.concatenate((clus22all,miefile[42:44,:]),axis=0)
		clus23all=np.concatenate((clus23all,miefile[44:46,:]),axis=0)
		clus24all=np.concatenate((clus24all,miefile[46:48,:]),axis=0)
		clus25all=np.concatenate((clus25all,miefile[48:50,:]),axis=0)
		clus26all=np.concatenate((clus26all,miefile[50:52,:]),axis=0)
		clus27all=np.concatenate((clus27all,miefile[52:54,:]),axis=0)
		clus28all=np.concatenate((clus28all,miefile[54:56,:]),axis=0)
		clus29all=np.concatenate((clus29all,miefile[56:58,:]),axis=0)

	##Backup each cluster input params in csv files------
	
	clus1all=pd.DataFrame(data=clus1all) 		
	clus1all.to_csv(path_out2+'clus_global_1'+'.csv', index=False,header=False) 
	clus2all=pd.DataFrame(data=clus2all) 		
	clus2all.to_csv(path_out2+'clus_global_2'+'.csv', index=False,header=False) 
	clus3all=pd.DataFrame(data=clus3all) 		
	clus3all.to_csv(path_out2+'clus_global_3'+'.csv', index=False,header=False) 
	clus4all=pd.DataFrame(data=clus4all) 		
	clus4all.to_csv(path_out2+'clus_global_4'+'.csv', index=False,header=False) 
	clus5all=pd.DataFrame(data=clus5all) 		
	clus5all.to_csv(path_out2+'clus_global_5'+'.csv', index=False,header=False) 
	clus6all=pd.DataFrame(data=clus6all) 		
	clus6all.to_csv(path_out2+'clus_global_6'+'.csv', index=False,header=False) 
	clus7all=pd.DataFrame(data=clus7all) 		
	clus7all.to_csv(path_out2+'clus_global_7'+'.csv', index=False,header=False) 
	clus8all=pd.DataFrame(data=clus8all) 		
	clus8all.to_csv(path_out2+'clus_global_8'+'.csv', index=False,header=False) 
	clus9all=pd.DataFrame(data=clus9all) 		
	clus9all.to_csv(path_out2+'clus_global_9'+'.csv', index=False,header=False) 
	clus10all=pd.DataFrame(data=clus10all) 		
	clus10all.to_csv(path_out2+'clus_global_10'+'.csv', index=False,header=False) 
	clus11all=pd.DataFrame(data=clus11all) 		
	clus11all.to_csv(path_out2+'clus_global_11'+'.csv', index=False,header=False) 
	clus12all=pd.DataFrame(data=clus12all) 		
	clus12all.to_csv(path_out2+'clus_global_12'+'.csv', index=False,header=False) 
	clus13all=pd.DataFrame(data=clus13all) 		
	clus13all.to_csv(path_out2+'clus_global_13'+'.csv', index=False,header=False) 
	clus14all=pd.DataFrame(data=clus14all) 		
	clus14all.to_csv(path_out2+'clus_global_14'+'.csv', index=False,header=False) 
	clus15all=pd.DataFrame(data=clus15all) 		
	clus15all.to_csv(path_out2+'clus_global_15'+'.csv', index=False,header=False) 
	clus16all=pd.DataFrame(data=clus16all) 		
	clus16all.to_csv(path_out2+'clus_global_16'+'.csv', index=False,header=False) 
	clus17all=pd.DataFrame(data=clus17all) 		
	clus17all.to_csv(path_out2+'clus_global_17'+'.csv', index=False,header=False) 
	clus18all=pd.DataFrame(data=clus18all) 		
	clus18all.to_csv(path_out2+'clus_global_18'+'.csv', index=False,header=False) 
	clus19all=pd.DataFrame(data=clus19all) 		
	clus19all.to_csv(path_out2+'clus_global_19'+'.csv', index=False,header=False) 
	clus20all=pd.DataFrame(data=clus20all) 		
	clus20all.to_csv(path_out2+'clus_global_20'+'.csv', index=False,header=False) 
	clus21all=pd.DataFrame(data=clus21all) 		
	clus21all.to_csv(path_out2+'clus_global_21'+'.csv', index=False,header=False) 
	clus22all=pd.DataFrame(data=clus22all) 		
	clus22all.to_csv(path_out2+'clus_global_22'+'.csv', index=False,header=False) 
	clus23all=pd.DataFrame(data=clus23all) 		
	clus23all.to_csv(path_out2+'clus_global_23'+'.csv', index=False,header=False) 
	clus24all=pd.DataFrame(data=clus24all) 		
	clus24all.to_csv(path_out2+'clus_global_24'+'.csv', index=False,header=False) 
	clus25all=pd.DataFrame(data=clus25all) 		
	clus25all.to_csv(path_out2+'clus_global_25'+'.csv', index=False,header=False) 
	clus26all=pd.DataFrame(data=clus26all) 		
	clus26all.to_csv(path_out2+'clus_global_26'+'.csv', index=False,header=False) 
	clus27all=pd.DataFrame(data=clus27all) 		
	clus27all.to_csv(path_out2+'clus_global_27'+'.csv', index=False,header=False) 
	clus28all=pd.DataFrame(data=clus28all) 		
	clus28all.to_csv(path_out2+'clus_global_28'+'.csv', index=False,header=False) 
	clus29all=pd.DataFrame(data=clus29all) 		
	clus29all.to_csv(path_out2+'clus_global_29'+'.csv', index=False,header=False) 	

		
	print('creating outfile...........')
	#defining columns
	row1=''
	mm=str(17)+' '+' '+' '+str(80)	
	chunk=[mm,',',',','nlamb',',','nsd']

	row2=""
	for ele in chunk:  
		row2 += ele 
	mm=str(0)+str(1)+' '+' '+' '+str(17)
	chunk=[mm,',',',','ilim1',',','ilim2']
	
	row3=""
	for ele in chunk:  
		row3 += ele 
	mm=str(0)+str(1)+' '+' '+' '+str(80)
	chunk=[mm,',',',','isd1',',','isd2']
	
	row4=""
	for ele in chunk:  
		row4 += ele 
	
	
	############## WE NEED TO LOOP HERE Cvf:Cvc 10 levels#######################################################
	####0.95:0.05, 0.80:0.20, 0.50:0.50, 0.30:0.70, 0.20:0.80, 0.10:0.90, 0.05:0.95, 0.02:0.98, 0.01:0.90, 0:1
	#### WE NEED TO OPEN MEAN AEROSOL PROPERTIES AND COMPUTE ND FRACTIONS--------------
	####RH IS 60 IN AVERAGE ALWAYS FOR AERONET CLUSTERS---IDEALLY SHOULD BE MEANR RH OF DATA FOR THAT SPECIFIC CLUSTER-----------------
	for f2 in sorted(glob.glob(path_out2+'*.csv')): 
		print(f2)
		df = pd.read_csv(f2,header=None,error_bad_lines=False) 
		df=df.values[2:,:]
		clusname=f2[71:]
		print(clusname)
			
	##################################
	#variables that we do not modify for fine/coarse modes---
		wlam=df[:,0] #wavelength
		nr=df[:,3]
		ni=df[:,4]
		delx=df[:,6]
		am=df[:,7] #aerosol model
		
		wlam=wlam.reshape(len(wlam),1)
		am=am.reshape(len(am),1)	
		nr=nr.reshape(len(nr),1)	
		ni=ni.reshape(len(ni),1)	
		delx=delx.reshape(len(delx),1)	
	###################################
	
	#converting rv to rn-------
		rv=df[:,1] #34 (17 wavelengths x 2 psd modes) x 1
		sigma_v=df[:,2]
		rm_n,mun,muv=rv_to_rm2(rv,sigma_v) #sigma_n = sigma_v, rv=aeronet volume space, rn=zia N space
		rm_n=rm_n.reshape(len(rm_n),1)	
		mun=mun.reshape(len(mun),1)
		sigma_v=sigma_v.reshape(len(sigma_v),1)
			
		voall_f=[0.95,0.80,0.50,0.30,0.20,0.10,0.05,0.02,0.01,0.0]
		voall_c=[0.05,0.20,0.50,0.70,0.80,0.90,0.95,0.98,0.99,1.0]
	
	########## each iteration will be 1 file with a specific cvf:cfc 
		for k in range(0,len(voall_f)): #10 cvf:cfc volume concentration fractions levels	
			vo_f=voall_f[k]
			vo_c=voall_c[k]
			vo_r=[vo_f,vo_c]
			vo_r=np.asarray(vo_r)
			vo_r=vo_r.reshape(len(vo_r),1)
			vo_r=np.tile(vo_r,(17,1))
		
			RH=60 #this is used by default like a mean RH between 30 and 95%
		
		#####################################
		#computing number density fractions-----		
			no_f,no_c=nd_frac2(vo_r,mun,sigma_v)
			no_f=no_f.reshape(len(no_f),1)
			no_c=no_c.reshape(len(no_c),1)
			no_r=np.zeros([len(no_f)*2,1])
			no_r[0::2]=no_f
			no_r[1::2]=no_c
			ndfcase=np.zeros([len(no_f)*2,1])
			ndfcase[0::2]=k+1
			ndfcase[1::2]=k+1
		
		####--MERGINg MODTRAN INPUT VARIABLES-----------------------------
			df2=np.concatenate((wlam,nr),axis=1)
			df2=np.concatenate((df2,ni),axis=1)
			df2=np.concatenate((df2,rm_n),axis=1)
			df2=np.concatenate((df2,sigma_v),axis=1)
			df2=np.concatenate((df2,no_r),axis=1)
			df2=np.concatenate((df2,delx),axis=1)
		
		
		####--MERGIN FORTRAN INPUT VARIABLES-----------------------------
			df3=np.concatenate((wlam,rm_n),axis=1)
			df3=np.concatenate((df3,sigma_v),axis=1)
			df3=np.concatenate((df3,nr),axis=1)
			df3=np.concatenate((df3,ni),axis=1)
			df3=np.concatenate((df3,no_r),axis=1)
			df3=np.concatenate((df3,am),axis=1)
			df3=np.concatenate((df3,ndfcase),axis=1)
					
			mm="{:.5f}".format(df[0,5])
			ndf_f=str(mm)
		
			volc_frac='(CvF:CvC='+str(vo_f)+':'+str(vo_c)+')'
			ndfrac_f='(WT='+ndf_f
	
			RHst='5 RH_'+str(RH)+' '+'percent)'
	
	#size distribution: (	WT=0.99994	5 RH_30 pe	rcent) (Cv	F:CvC=0.95	0:0.050)	
	#chunk=['size distribution: ',',','(WT=0.00000',',','5 RH_30 percent)',',','(CvF:CvC=0.95	0:0.050)']
			chunk=['size distribution: ',ndfrac_f,',',RHst,volc_frac]
			row5=""
			for ele in chunk:  #from list to string
				row5 += ele 	
	
	
			mm=str(4)+' '+' '+' '+' '+str(0)
			chunk=[mm,',',',','ifunc i',',','ifunc_sub']
			row6=""		
			for ele in chunk:  #from list to string
				row6 += ele 
	
			mm=str(1)+' '+' '+' '+str(RH)
			chunk=[mm,',',str(3),',','irgm   R',',','H    cset']
			row7=""		
			for ele in chunk:  #from list to string
				row7 += ele
	#header--------------------------------------
			chunk=['lamb',',','n1',',','n2',',','rm',',','sig',',','n',',','delx']
			row8=""		
			for ele in chunk:  #from list to string
				row8 += ele
	
	
			print('processing cluster...',clusname)
			print('volume concentration fraction...',volc_frac)
			outfile=path_out2+'MODTRAN_input_globe'+'_'+clusname+volc_frac+'.txt'
			with open(outfile, "w+") as f: 
					f.write(row1+"\n")
					f.write(row2+"\n")#number of wavelengths and number of legendre polynomials =80 default
					f.write(row3+"\n")#limits for wavelength, 17 default
					f.write(row4+"\n") #limits for of legendre polynomial 01 80 default
					f.write(row5+"\n")#size distr RH WT and CvfCvc
					f.write(row6+"\n")# ifunc-----integration function step?, delx 0.05 fine 0.15 coarse---default
					f.write(row7+"\n")#RH limits ??
					f.write(row8+"\n")#header
			
			#open one aerosol model and write properties for each wavelength---
			#for each line--------------------
			#separate with ','
			
					for i in df2:
						chunk=i.tolist()
						row_data=""		
						for ele in chunk:  #from list to string
							row_data += str(ele)+','			
						f.write(row_data+"\n")#header				
			
    			#f.write(row9+"\n")#monochr
    			#f.write(row10+"\n")#lambda
    			#f.write(row11+"\n")#homo surf
    			#f.write(row12+"\n")#direct effect
    			#f.write(row13+"\n")#ocean
    			#f.write(row14+"\n")#wind water
    			#f.write(row15+"\n")#no atmo corre flag   				  		
			f.close()
			
			#### Write FORTRAIN INPUT---#############################
			print('processing FORTRAN input file..',clusname)
			print('volume concentration fraction...',volc_frac)
			outfile=path_out2+'FORTRAN_input_globe'+'_'+clusname+volc_frac+'.txt'
			with open(outfile, "w+") as f3: 
			
					for i in df3:
						chunk=i.tolist()			
						row_data=""		
						for ele in chunk:  #from list to string
							row_data += str(ele)+','
				
						f3.write(row_data+"\n")#header
			
    			#f.write(row9+"\n")#monochr
    			#f.write(row10+"\n")#lambda
    			#f.write(row11+"\n")#homo surf
    			#f.write(row12+"\n")#direct effect
    			#f.write(row13+"\n")#ocean
    			#f.write(row14+"\n")#wind water
    			#f.write(row15+"\n")#no atmo corre flag   				  		
			f3.close()
			
	

if __name__ == '__main__':

	inp_miefile()