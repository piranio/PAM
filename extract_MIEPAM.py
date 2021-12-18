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
##Current version derived from extract_Mie_propPAM.py
##*************************************

#****************MAIN************************--------------------------------------------------------	

if __name__ == '__main__':
	
	#initialization of arrays
	ph=[]
	qh=[]
	uh=[]
	g=[]
	ssa=[]
	kext=[]
	ksca=[]
	sca=[]
	wl=[]
	nr=[]
	ni=[]
	nr2=[]
	ni2=[]
	rm=[]	
	std=[]
	rm2=[]
	std2=[]
	volf=[] 
	volf2=[]

	path_in='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/FORTRAN/in/PAM/'
	path_out='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/FORTRAN/out/'
	
	#### we loop each model in PAM ###################
	for f in sorted(glob.glob(path_in+'*ok')):
		fname=f[-9:-2:]
		print('processing model..#..',fname)	
		an=0
		c=1	
		while c<=170: #while not eof 16 = 4 cluster x 4 wavelengths		
			line=file1.readline()
			mm=line.split()
			if mm[2]==':Log-Normal':
				for i in range(4): #skip these lines
					next(file1)				
		#refractive index----FINE MODE (in AERONET fine and coarse nr and ni are the same)
			line=file1.readline()
			mm=line.split()
			nr.append(mm[7])
			ni.append(mm[8])
				
			#radius-std----line7----fine mode
			line=file1.readline()
			mm=line.split()
			rm.append(mm[5])
			std.append(mm[10])
			
		#vol fract-----------------
			line=file1.readline()
			mm=line.split()
			volf.append(mm[4])
							
			for i in range(2): #skip these lines
				next(file1)	
			line=file1.readline()
			mm=line.split()
		#refractive index----COARSE MODE
			nr2.append(mm[7])
			ni2.append(mm[8])
			
		#radius-std----line12----coarse mode
			line=file1.readline()
			mm=line.split()
			rm2.append(mm[5])
			std2.append(mm[10])
			
		#vol fract-----------------
			line=file1.readline()
			mm=line.split()
			volf2.append(mm[4])
		
			line=file1.readline()
			line=file1.readline()
			mm=line.split()
			
			if mm[0]=='180.000': #scattering angle			
				while an<99:
					sca.append(mm[0])
					ph.append(mm[1])
					qh.append(mm[2])
					uh.append(mm[3])
					line=file1.readline()
					mm=line.split()					
					an+=1
	
				an=0
				wl.append(mm[2])
				line=file1.readline()
				mm=line.split()
				#print(mm)
				#sys.exit()				
				ssa.append(mm[2])
				line=file1.readline()
				mm=line.split()
				kext.append(mm[2])
				line=file1.readline()
				mm=line.split()
				ksca.append(mm[2])
				line=file1.readline()
				mm=line.split()
				g.append(mm[4])
				
				if c<=170:			
					for i in range(18):
						line=file1.readline()
			
				mm=line.split()
								
				c+=1
					
	#### CREATE RADIATIVE PROPS TABLE ###############
	
	#merge all ndf cases in one table------
	#we need to split in chunks of 10 rows, each chunk is 10 ndfs and 1 lambda
	
	#lambda 1----------------
		i=0 
		wl1=wl[0+i*10:10+i*10] #first lambda for 10 ndfs
		ssa1=ssa[0+i*10:10+i*10] 
		g1=g[0+i*10:10+i*10] 
		ksca1=ksca[0+i*10:10+i*10] 
		kext1=kext[0+i*10:10+i*10] 
		rf=rm[0+i*10:10+i*10] 
		sigmaf=std[0+i*10:10+i*10] 
		rc=rm2[0+i*10:10+i*10] 
		sigmac=std2[0+i*10:10+i*10]	
		nrf1=nr[0+i*10:10+i*10] 
		nif1=ni[0+i*10:10+i*10] 
		nrc1=nr2[0+i*10:10+i*10] 
		nic1=ni2[0+i*10:10+i*10]
		ndff=volf[0+i*10:10+i*10] 
		ndfc=volf2[0+i*10:10+i*10]
	
	#lambda 2----------------
		i=1 
		wl2=wl[0+i*10:10+i*10] 
		ssa2=ssa[0+i*10:10+i*10] 
		g2=g[0+i*10:10+i*10] 
		ksca2=ksca[0+i*10:10+i*10] 
		kext2=kext[0+i*10:10+i*10] 
		nrf2=nr[0+i*10:10+i*10] 
		nif2=ni[0+i*10:10+i*10] 
		nrc2=nr2[0+i*10:10+i*10] 
		nic2=ni2[0+i*10:10+i*10]
	
	#lambda 3----------------
		i=2 
		wl3=wl[0+i*10:10+i*10] 
		ssa3=ssa[0+i*10:10+i*10] 
		g3=g[0+i*10:10+i*10] 
		ksca3=ksca[0+i*10:10+i*10] 
		kext3=kext[0+i*10:10+i*10] 
		nrf3=nr[0+i*10:10+i*10] 
		nif3=ni[0+i*10:10+i*10] 
		nrc3=nr2[0+i*10:10+i*10] 
		nic3=ni2[0+i*10:10+i*10]
		
	#lambda 4----------------
		i=3 
		wl4=wl[0+i*10:10+i*10] 
		ssa4=ssa[0+i*10:10+i*10] 
		g4=g[0+i*10:10+i*10] 
		ksca4=ksca[0+i*10:10+i*10] 
		kext4=kext[0+i*10:10+i*10] 
		nrf4=nr[0+i*10:10+i*10] 
		nif4=ni[0+i*10:10+i*10] 
		nrc4=nr2[0+i*10:10+i*10] 
		nic4=ni2[0+i*10:10+i*10]
	
	#lambda 5----------------
		i=4 
		wl5=wl[0+i*10:10+i*10] 
		ssa5=ssa[0+i*10:10+i*10] 
		g5=g[0+i*10:10+i*10] 
		ksca5=ksca[0+i*10:10+i*10] 
		kext5=kext[0+i*10:10+i*10] 
		nrf5=nr[0+i*10:10+i*10] 
		nif5=ni[0+i*10:10+i*10] 
		nrc5=nr2[0+i*10:10+i*10] 
		nic5=ni2[0+i*10:10+i*10]
	
	
	#lambda 6----------------
		i=5 
		wl6=wl[0+i*10:10+i*10]
		ssa6=ssa[0+i*10:10+i*10] 
		g6=g[0+i*10:10+i*10] 
		ksca6=ksca[0+i*10:10+i*10] 
		kext6=kext[0+i*10:10+i*10] 
		nrf6=nr[0+i*10:10+i*10] 
		nif6=ni[0+i*10:10+i*10] 
		nrc6=nr2[0+i*10:10+i*10] 
		nic6=ni2[0+i*10:10+i*10]
	
	
	#lambda 7----------------
		i=6 
		wl7=wl[0+i*10:10+i*10] 
		ssa7=ssa[0+i*10:10+i*10] 
		g7=g[0+i*10:10+i*10] 
		ksca7=ksca[0+i*10:10+i*10] 
		kext7=kext[0+i*10:10+i*10] 
		nrf7=nr[0+i*10:10+i*10] 
		nif7=ni[0+i*10:10+i*10] 
		nrc7=nr2[0+i*10:10+i*10] 
		nic7=ni2[0+i*10:10+i*10]
	
	
	#lambda 8----------------
		i=7 
		wl8=wl[0+i*10:10+i*10] 
		ssa8=ssa[0+i*10:10+i*10] 
		g8=g[0+i*10:10+i*10] 
		ksca8=ksca[0+i*10:10+i*10] 
		kext8=kext[0+i*10:10+i*10] 
		nrf8=nr[0+i*10:10+i*10] 
		nif8=ni[0+i*10:10+i*10] 
		nrc8=nr2[0+i*10:10+i*10] 
		nic8=ni2[0+i*10:10+i*10]
	
	#lambda 9----------------
		i=8 
		wl9=wl[0+i*10:10+i*10] 
		ssa9=ssa[0+i*10:10+i*10] 
		g9=g[0+i*10:10+i*10] 
		ksca9=ksca[0+i*10:10+i*10] 
		kext9=kext[0+i*10:10+i*10] 
		nrf9=nr[0+i*10:10+i*10] 
		nif9=ni[0+i*10:10+i*10] 
		nrc9=nr2[0+i*10:10+i*10] 
		nic9=ni2[0+i*10:10+i*10]
	
	#lambda 10----------------
		i=9 
		wl10=wl[0+i*10:10+i*10]
		ssa10=ssa[0+i*10:10+i*10] 
		g10=g[0+i*10:10+i*10] 
		ksca10=ksca[0+i*10:10+i*10] 
		kext10=kext[0+i*10:10+i*10] 
		nrf10=nr[0+i*10:10+i*10] 
		nif10=ni[0+i*10:10+i*10] 
		nrc10=nr2[0+i*10:10+i*10] 
		nic10=ni2[0+i*10:10+i*10]
	
	#lambda 11----------------
		i=10 
		wl11=wl[0+i*10:10+i*10] 
		ssa11=ssa[0+i*10:10+i*10] 
		g11=g[0+i*10:10+i*10] 
		ksca11=ksca[0+i*10:10+i*10] 
		kext11=kext[0+i*10:10+i*10] 
		nrf11=nr[0+i*10:10+i*10] 
		nif11=ni[0+i*10:10+i*10] 
		nrc11=nr2[0+i*10:10+i*10] 
		nic11=ni2[0+i*10:10+i*10]
	
	#lambda 12----------------
		i=11 
		wl12=wl[0+i*10:10+i*10] 
		ssa12=ssa[0+i*10:10+i*10] 
		g12=g[0+i*10:10+i*10] 
		ksca12=ksca[0+i*10:10+i*10] 
		kext12=kext[0+i*10:10+i*10] 
		nrf12=nr[0+i*10:10+i*10] 
		nif12=ni[0+i*10:10+i*10] 
		nrc12=nr2[0+i*10:10+i*10] 
		nic12=ni2[0+i*10:10+i*10]
	
	#lambda 13----------------
		i=12 
		wl13=wl[0+i*10:10+i*10] 
		ssa13=ssa[0+i*10:10+i*10] 
		g13=g[0+i*10:10+i*10] 
		ksca13=ksca[0+i*10:10+i*10] 
		kext13=kext[0+i*10:10+i*10] 
		nrf13=nr[0+i*10:10+i*10] 
		nif13=ni[0+i*10:10+i*10] 
		nrc13=nr2[0+i*10:10+i*10] 
		nic13=ni2[0+i*10:10+i*10]
	
	#lambda 14----------------
		i=13 
		wl14=wl[0+i*10:10+i*10] 
		ssa14=ssa[0+i*10:10+i*10] 
		g14=g[0+i*10:10+i*10] 
		ksca14=ksca[0+i*10:10+i*10] 
		kext14=kext[0+i*10:10+i*10] 
		nrf14=nr[0+i*10:10+i*10] 
		nif14=ni[0+i*10:10+i*10] 
		nrc14=nr2[0+i*10:10+i*10] 
		nic14=ni2[0+i*10:10+i*10]
	
	#lambda 15----------------
		i=14 
		wl15=wl[0+i*10:10+i*10] 
		ssa15=ssa[0+i*10:10+i*10] 
		g15=g[0+i*10:10+i*10] 
		ksca15=ksca[0+i*10:10+i*10] 
		kext15=kext[0+i*10:10+i*10] 
		nrf15=nr[0+i*10:10+i*10] 
		nif15=ni[0+i*10:10+i*10] 
		nrc15=nr2[0+i*10:10+i*10] 
		nic15=ni2[0+i*10:10+i*10]
	
	#lambda 16----------------
		i=15 
		wl16=wl[0+i*10:10+i*10] 
		ssa16=ssa[0+i*10:10+i*10] 
		g16=g[0+i*10:10+i*10] 
		ksca16=ksca[0+i*10:10+i*10] 
		kext16=kext[0+i*10:10+i*10] 
		nrf16=nr[0+i*10:10+i*10] 
		nif16=ni[0+i*10:10+i*10] 
		nrc16=nr2[0+i*10:10+i*10] 
		nic16=ni2[0+i*10:10+i*10]
	
	#lambda 17----------------
		i=16 
		wl17=wl[0+i*10:10+i*10] 
		ssa17=ssa[0+i*10:10+i*10] 
		g17=g[0+i*10:10+i*10] 
		ksca17=ksca[0+i*10:10+i*10] 
		kext17=kext[0+i*10:10+i*10] 
		nrf17=nr[0+i*10:10+i*10] 
		nif17=ni[0+i*10:10+i*10] 
		nrc17=nr2[0+i*10:10+i*10] 
		nic17=ni2[0+i*10:10+i*10]
		
		table1=np.zeros([10,7])
		for k in range(0,10): #10 different number density fractions!!
			#print(k)
			#print('processing ndf number##....',k+1)	
			wl_ndf=[wl1[k],wl2[k],wl3[k],wl4[k],wl5[k],wl6[k],wl7[k],wl8[k],wl9[k],wl10[k],wl11[k],wl12[k],wl13[k],wl14[k],wl15[k],wl16[k],wl17[k]]
			ssa_ndf=[ssa1[k],ssa2[k],ssa3[k],ssa4[k],ssa5[k],ssa6[k],ssa7[k],ssa8[k],ssa9[k],ssa10[k],ssa11[k],ssa12[k],ssa13[k],ssa14[k],ssa15[k],ssa16[k],ssa17[k]]
			g_ndf=[g1[k],g2[k],g3[k],g4[k],g5[k],g6[k],g7[k],g8[k],g9[k],g10[k],g11[k],g12[k],g13[k],g14[k],g15[k],g16[k],g17[k]]
			kext_ndf=[kext1[k],kext2[k],kext3[k],kext4[k],kext5[k],kext6[k],kext7[k],kext8[k],kext9[k],kext10[k],kext11[k],kext12[k],kext13[k],kext14[k],kext15[k],kext16[k],kext17[k]]
			ksca_ndf=[ksca1[k],ksca2[k],ksca3[k],ksca4[k],ksca5[k],ksca6[k],ksca7[k],ksca8[k],ksca9[k],ksca10[k],ksca11[k],ksca12[k],ksca13[k],ksca14[k],ksca15[k],ksca16[k],ksca17[k]]
			ndff_ndf=[ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k],ndff[k]]
			ndfc_ndf=[ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k],ndfc[k]]
		
		#convert to numpy---
			wl_ndf=np.asarray(wl_ndf)
			ssa_ndf=np.asarray(ssa_ndf)
			g_ndf=np.asarray(g_ndf)
			kext_ndf=np.asarray(kext_ndf)
			ksca_ndf=np.asarray(ksca_ndf)
			ndff_ndf=np.asarray(ndff_ndf)
			ndfc_ndf=np.asarray(ndfc_ndf)
		
		#reshape from zero 1 -d
			wl_ndf=wl_ndf.reshape(len(wl_ndf),1)
			ssa_ndf=ssa_ndf.reshape(len(ssa_ndf),1)
			g_ndf=g_ndf.reshape(len(g_ndf),1)
			kext_ndf=kext_ndf.reshape(len(kext_ndf),1)
			ksca_ndf=ksca_ndf.reshape(len(ksca_ndf),1)
			ndff_ndf=ndff_ndf.reshape(len(ndff_ndf),1)
			ndfc_ndf=ndfc_ndf.reshape(len(ndfc_ndf),1)
			
		#---concatenate
			one_ndf=np.concatenate((wl_ndf,ssa_ndf),axis=1)
			one_ndf=np.concatenate((one_ndf,g_ndf),axis=1)
			one_ndf=np.concatenate((one_ndf,kext_ndf),axis=1)
			one_ndf=np.concatenate((one_ndf,ksca_ndf),axis=1)
			one_ndf=np.concatenate((one_ndf,ndff_ndf),axis=1)
			one_ndf=np.concatenate((one_ndf,ndfc_ndf),axis=1)
		
			table1=np.concatenate((table1,one_ndf),axis=0)
		
	
		table1=pd.DataFrame(table1)
		table1.to_csv(path_out+'pam_'+fname+'_rad'+'.csv', header=False,index=False)
		print('finish extracting model #..',fname)
		
		file1.close()
		
		#initialize the arrays again for the next cluster associated file-----
		ph=[]
		qh=[]
		uh=[]
		g=[]
		ssa=[]
		kext=[]
		ksca=[]
		sca=[]
		wl=[]
		nr=[]
		ni=[]
		nr2=[]
		ni2=[]
		rm=[]	
		std=[]
		rm2=[]
		std2=[]
		volf=[] 
		volf2=[]	
	
	
	
	
		
		
