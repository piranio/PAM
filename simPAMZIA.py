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
from sklearn.metrics import mean_squared_error as mse
from sklearn.metrics import mean_absolute_error as mae
from numpy import linalg as LA

##*************************************
##Author: Martin A. Montes
##SSAI, April 29, 2020
#version derived from cosimZIAPAM7.py
##*************************************


def cosim(a,b):

# manually compute cosine similarity
	dot = np.dot(a, b)
	norma = np.linalg.norm(a)
	normb = np.linalg.norm(b)
	cosind = dot / (norma * normb) #manual way to compute the cosine similarity
 
# use library that operates on sets of vectors
	#aa = a.reshape(1,3)
	#ba = b.reshape(1,3)
	#cos_lib = cosine_similarity(aa, ba)
	return cosind

#*************MAIN**************--------------------------------------------------------	

if __name__ == '__main__':
	
	#these are radiative aerosol properties for PAM models and all NDFs variations per cluster 
	#computed from Mie code, col.1 wavelength, 2.SSA, 3.g, 4. Kext. 5. Ka, 6. NDFfine and NDFcoarse	
	path_pam='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/SIM/in/PAM/'
	#same as above but for A2010 models
	path_zia='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/SIM/in/ZIA/'
	
	#Sorted output Mie files described in path_pam and based on descending NDFfine
	path_out='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/SIM/out/'
	
	#Metrics results
	path_out2='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/SIM/results/'
	
	for f in sorted(glob.glob(path_pam+'*.csv')):
		fname=f[-10:-2:]
		print('processing model..#..',fname)
		file_pam= pd.read_csv(f,header=None,error_bad_lines=False)
		file_pam=file_pam.values
		file_pam=file_pam[9:,:]
			
		sorted_pam=sorted(file_pam, reverse=True,key=lambda x: x[5]) #will create a sorted list		
		sorted_pam=pd.DataFrame(sorted_pam)
		sorted_pam.to_csv(path_out+'pam_'+'sortedMod'+fname+'.csv', header=False,index=False)
	
	
	####### OPEN ZIA MODELS (all of them in 1 file)###############
	for f2 in sorted(glob.glob(path_zia+'*.csv')):
		
		file_zia= pd.read_csv(f2,header=None,error_bad_lines=False)
		file_zia=file_zia.values
		file_zia=file_zia[80:,:]
		
		#print(file_zia.shape) #1360 x 8=17lam x 8 mod x 10 ndf,col1:lambda col2:ssa col3:g col4:kext col5:ksca col6:ndff col7:ndfc col8:modnum 1-80
		
		#ONE NDF AT THE TIME!!!!
		nd2=0
		########## WE NEED TO PRESELECT WHICH AEROSOL PROPERTY WE WANT TO COMPARE ##############
		op=3   #op1: ssa, op2:g, op3: kext, op4: ksca
		cosim1=[] #to backup 29 pam per 1 zia
		rmse1=[]
		nrmse1=[]
		mae1=[]
		nmae1=[]
		#################################################################################
		#Initalizations of arrays holding the resulting metric results---		
		for nd in range(0,9): #ndf 1 to 10, ndf=10 is ignored cause ndff=0 so nan in Mie simulations
			cosim2=np.zeros([29,8]) #29 pam x 8 zia models per ndf
			rmse2=np.zeros([29,8])	#root mean square error	
			mae2=np.zeros([29,8]) #mean absolute error
			nmae2=np.zeros([29,8]) #normalized mae
			nrmse2=np.zeros([29,8]) #normalized rmse
			
			ndf=nd+1 #1 to 10 number density fractions
			print('ndf case study number ....#',ndf)
			print('aerosol property #..op1:ssa,op2:g,op3:kext.................',op)
			c=0 #zia model 1-8
			zia_mods=[]
			for ziamod in range(ndf,81,10):
				print('ziamod number...................#',ziamod)
				zia_mod=str(ziamod)
				zia_mods.append(ziamod)
				izia=np.where(file_zia[:,7]==ziamod) #find all ziamod with ndff=?? (1 to 10 ndff descending order)
				one_zia=file_zia[izia,:]
				one_zia=np.squeeze(one_zia) ##17 x 8 params
				one_zia_p=one_zia[:,op] #17 X 1 parameter
				
		
		#CALLING PAM models for the above ndf case study----
		#####################################################
		####### OPEN 1 PAM model at the time (10 ndf levels)#################
				pam_mods=[]
				for f in sorted(glob.glob(path_out+'*.csv')): #FOR EACH PAM MODEL!! not read in order though for name differences..
					fname=f[-15:-10:]
					pam_mods.append(fname)
					print('reading sorted PAM model...',fname)
					file_pamo= pd.read_csv(f,header=None,error_bad_lines=False)
					file_pamo=file_pamo.values			
					one_pam=file_pamo[nd2:nd2+17:,:]
					
		########## COMPUTE COS SIMILARITY INDEX: FOR EACH NDF I DO  29 PAM X 8 ZIA ###############
		########## WE NEED TO PRESELECT WHICH AEROSOL PROPERTY WE WANT TO COMPARE ##############
					one_pam_p=one_pam[:,op] #op1: ssa, op2:g, op3: kext, op4: ksca					
					mm=cosim(one_zia_p,one_pam_p) #row-wise and dimension zero-cosine similarity
					one_zia_pt=one_zia_p.T
					one_pam_pt=one_pam_p.T					
					mm2=mse(one_zia_pt,one_pam_pt)  #col-wise and dimens zero--root mean square error
					mm3=mae(one_zia_pt,one_pam_pt)
					nn=LA.norm(one_zia_pt-one_pam_pt)
					mm4=mm3/nn
					mm5=mm2/nn
				
	##*******************METRICS CALC***********************************************		
		#Similarity and accuracy performance metrics calculation			
					cosim1.append(mm) #1 zia x 29 pam--cosine similarity
					rmse1.append(mm2) #root mean square root
					mae1.append(mm3)#mae or mean absolute error
					nmae1.append(mm4)#normalized mae
					nrmse1.append(mm5)#normalized mae
				
				
				mm=np.asarray(cosim1)
				mm2=np.asarray(rmse1)
				mm3=np.asarray(mae1)
				mm4=np.asarray(nmae1)
				mm5=np.asarray(nrmse1)
								
				cosim2[:,c]=mm		#8 zia x 29 pam models
				rmse2[:,c]=mm2
				mae2[:,c]=mm3
				nmae2[:,c]=mm4
				nrmse2[:,c]=mm5
								
				cosim1=[]
				rmse1=[]
				mae1=[]
				nmae1=[]
				nrmse1=[]
				
				c+=1 #zia model 8 per ndf
			nd2=nd2+17 #increases when nd increases
			ndff_d=nd+1
			ndff_d=str(ndff_d)
			
	##**********BACKUP OF RESULTS*************************				
			zia_mods=pd.DataFrame(zia_mods)
			zia_mods=zia_mods.T #8 columns
			pam_mods=pd.DataFrame(pam_mods) #29 rows
		
		#cosine similarity
			cosim2=pd.DataFrame(cosim2)
			merged1 = pd.concat([zia_mods,cosim2], axis=0)
			merged1 = pd.concat([merged1,pam_mods], axis=1)
			#merged1.to_csv(path_out2+'COSIM_pamzia_'+ndff_d+'_'+'_ssa'+'.csv', header=False,index=False)
			#merged1.to_csv(path_out2+'COSIM_pamzia_'+ndff_d+'_'+'_g'+'.csv', header=False,index=False)
			merged1.to_csv(path_out2+'COSIM_pamzia_'+ndff_d+'_'+'_kext'+'.csv', header=False,index=False)
		#rmse
			rmse2=pd.DataFrame(rmse2)
			merged2 = pd.concat([zia_mods,rmse2], axis=0)
			merged2 = pd.concat([merged2,pam_mods], axis=1)			
			#merged2.to_csv(path_out2+'RMSE_pamzia_'+ndff_d+'_'+'_ssa'+'.csv', header=False,index=False)
			#merged2.to_csv(path_out2+'RMSE_pamzia_'+ndff_d+'_'+'_g'+'.csv', header=False,index=False)
			merged2.to_csv(path_out2+'RMSE_pamzia_'+ndff_d+'_'+'_kext'+'.csv', header=False,index=False)
		#mae							
			mae2=pd.DataFrame(mae2)
			merged3 = pd.concat([zia_mods,mae2], axis=0)
			merged3 = pd.concat([merged3,pam_mods], axis=1)	
			#merged3.to_csv(path_out2+'MAE_pamzia_'+ndff_d+'_'+'_ssa'+'.csv', header=False,index=False)
			#merged3.to_csv(path_out2+'MAE_pamzia_'+ndff_d+'_'+'_g'+'.csv', header=False,index=False)		
			merged3.to_csv(path_out2+'MAE_pamzia_'+ndff_d+'_'+'_kext'+'.csv', header=False,index=False)
		#name
			nmae2=pd.DataFrame(nmae2)
			merged4 = pd.concat([zia_mods,nmae2], axis=0)
			merged4 = pd.concat([merged4,pam_mods], axis=1)
			#merged4.to_csv(path_out2+'NMAE_pamzia_'+ndff_d+'_'+'_ssa'+'.csv', header=False,index=False)
			#merged4.to_csv(path_out2+'NMAE_pamzia_'+ndff_d+'_'+'_g'+'.csv', header=False,index=False)			
			merged4.to_csv(path_out2+'NMAE_pamzia_'+ndff_d+'_'+'_kext'+'.csv', header=False,index=False)
		#nrsme
			nrmse2=pd.DataFrame(nrmse2)
			merged5 = pd.concat([zia_mods,nrmse2], axis=0)
			merged5 = pd.concat([merged5,pam_mods], axis=1)
			#merged5.to_csv(path_out2+'NRMSE_pamzia_'+ndff_d+'_'+'_ssa'+'.csv', header=False,index=False)
			#merged5.to_csv(path_out2+'NRMSE_pamzia_'+ndff_d+'_'+'_g'+'.csv', header=False,index=False)			
			merged5.to_csv(path_out2+'NRMSE_pamzia_'+ndff_d+'_'+'_kext'+'.csv', header=False,index=False)
			
			
			print('finishing with ndff descending ..#',ndff_d)
			
		
			
		
		
		
		
		
		
		
		
		
		
		
		
		