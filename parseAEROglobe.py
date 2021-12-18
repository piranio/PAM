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

##*************************************
##Author: Martin A. Montes
##SSAI, April 29, 2020
##Current version derived from parseAEROglobe_nopf3.py
##*************************************


def parse_line(oneline,hh):		
	tok=oneline.split(",")
	header=hh.split(",")
	flag=1

	if float(tok[23])>0: #only values with ssa and absorbing properties
	#if flag==1: this is used in case we want values including those records without SSA
		sta.append(tok[0])
		
		#parsing date and time		
		mm=str(tok[1])
		ddmmyyyy = mm.split(':')
		mm2=str(tok[2])
		hhmmss= mm2.split(':')		
		day.append(ddmmyyyy[0])
		month.append(ddmmyyyy[1])
		year.append(ddmmyyyy[2])		
		hour.append(hhmmss[0])
		min.append(hhmmss[1])
		sec.append(hhmmss[2])
		
		#parsing geospatial data
		lat.append(tok[236])
		lon.append(tok[237])
		elev.append(tok[238])
		
		#parsing aerosol properties	
		aod1c.append(tok[5]) #coincident aod at 4 wavelengths aeronet
		aod2c.append(tok[6])
		aod3c.append(tok[7])
		aod4c.append(tok[8])
		aod1t.append(tok[10]) #aod extinction 440
		aod2t.append(tok[11])#675
		aod3t.append(tok[12])#870
		aod4t.append(tok[13])#1020
	
		ssa1.append(tok[23])#ssa443
		ssa2.append(tok[24])#ssa667
		ssa3.append(tok[25])#ssa865
		ssa4.append(tok[26])#ssa1020
		
		aod1a.append(tok[27]) #aod absorption440
		aod2a.append(tok[28])#675
		aod3a.append(tok[29])#870
		aod4a.append(tok[30])#1020
		g1.append(tok[40])#g443
		g2.append(tok[41])#g667
		g3.append(tok[42])#g865
		g4.append(tok[43])#g440
		
		nr1.append(tok[32])#nr443
		nr2.append(tok[33])#nr443
		nr3.append(tok[34])#nr443
		nr4.append(tok[35])#nr443
		ni1.append(tok[36])#nr443
		ni2.append(tok[37])#nr443
		ni3.append(tok[38])#nr443
		ni4.append(tok[39])#nr443
		
		volct.append(tok[76]) 
		vmrt.append(tok[78])
		stdt.append(tok[79])		
		volcf.append(tok[80]) 
		vmrf.append(tok[82])
		stdf.append(tok[83])	
		volcc.append(tok[84]) 
		vmrc.append(tok[86])
		stdc.append(tok[87])
			
		AEe1.append(tok[22])		
		AEe2.append(np.log(float(tok[12])/float(tok[13]))/np.log(1.02/0.87)) #870-1020 #based on two wavelengths		
		AEa.append(tok[31])	
	
	#outuput parsed results----
		path_out='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/L2V3/out/'
		aero1 = pd.DataFrame({'sta': sta,'year': year,'month': month,'day': day,'hour': hour,'min': min,'sec': sec,'lat':lat,'lon':lon,'elev':elev,'aod1c': aod1c,'aod2c': aod2c,'aod3c': aod3c,'aod4c': aod4c,'aod1t': aod1t,'aod2t': aod2t,'aod3t': aod3t,'aod4t': aod4t},columns=['sta','year','month','day','hour','min','sec','lat','lon','elev','aod1c','aod2c','aod3c','aod4c','aod1t','aod2t','aod3t','aod4t'])		
		aero2 = pd.DataFrame({'ssa1': ssa1,'ssa2': ssa2,'ssa3': ssa3,'ssa4': ssa4,'aod1a': aod1a,'aod2a': aod2a,'aod3a': aod3a,'aod4a': aod4a,'g1': g1,'g2': g2,'g3': g3,'g4': g4,'volcf':volcf,'vmrf':vmrf,'stdf':stdf,'volcc':volcc,'vmrc':vmrc,'stdc':stdc},columns=['ssa1','ssa2','ssa3','ssa4','aod1a','aod2a','aod3a','aod4a','g1','g2','g3','g4','volcf','vmrf','stdf','volcc','vmrc','stdc'])
		aero3 = pd.DataFrame({'nr1':nr1,'nr2':nr2,'nr3':nr3,'nr4':nr4,'ni1':ni1,'ni2':ni2,'ni3':ni3,'ni4':ni4,'AEe1': AEe1,'AEe2': AEe2,'AEa': AEa},columns=['AEe1','AEe2','AEa','nr1','nr2','nr3','nr4','ni1','ni2','ni3','ni4'])

		aero_all=pd.concat([aero1, aero2],axis=1)
		aero_all=pd.concat([aero_all, aero3],axis=1)	
		aero_all.to_csv(path_out+'globe_aeronet_extract_nopfALLfullAOD-part4.csv', index=False) 
	
	
	

def parse_aero(fullname,flagv):
	c=0
	f=open(fullname, 'r')
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	line = f.readline()
	#line = f.readline()
	hh=line #header
	
	while line:
		oneline=f.readline()
		if not oneline:
			print('done extracting file...')
			break		
		parse_line(oneline,hh) #flag_data =1 no pf, 2 = pf
		c+=1		
	f.close()
        	


#*******************MAIN****************************
	
if __name__ == '__main__':
		
	tok=[]
	sta=[]
	date=[]
	month=[]
	day=[]
	year=[]
	hour=[]
	min=[]
	sec=[]
	utc=[]
	lat=[]
	lon=[]
	elev=[]
	aod1c=[] #coincident aod at 440 nm
	aod2c=[] #coincident aod at 675 nm
	aod3c=[] #coincident aod at 870 nm
	aod4c=[] #coincident aod at 1020 nm
	aod1t=[]#total aod extinction at 440 nm
	aod2t=[] #total aod extinction at 675 nm
	aod3t=[] #total aod extinction at 870 nm
	aod4t=[] #total aod extinction at 1020 nm
	
	ssa1=[]#single scattering albedo at 440 nm
	ssa2=[] #single scattering albedo at 675 nm
	ssa3=[] #single scattering albedo at 870 nm
	ssa4=[] #single scattering albedo at 1020 nm
	
	aod1a=[] #total aod absorption at 440 nm
	aod2a=[]#total aod absorption at 675 nm
	aod3a=[]#total aod absorption at 870 nm
	aod4a=[]#total aod absorption at 1020 nm
	g1=[]#g440
	g2=[]#g675
	g3=[]#g870
	g4=[]#g1020
	
	volct=[] #total volumen concentration 
	vmrt=[]#total median log radius 
	stdt=[]#sigma
	volcf=[] #volumen concentration fine mode
	vmrf=[]#total median log radius 
	stdf=[]#sigma fine mode
	volcc=[] #volumen concentration coarse mode
	vmrc=[]#total median log radius coarse mode
	stdc=[]#sigma coarse mode
	
	AEe1=[] #Angstrom exponent for extinction 440-870 nm
	AEe2=[]#Angstrom exponent for extinction 870-1020 nm
	AEa=[] #Angstrom exponent for absorption 440-870 nm
	nr1=[] #refractive index 440 nm
	nr2=[] #refractive index 675 nm
	nr3=[] #refractive index 870 nm
	nr4=[] #refractive index 1020 nm
	ni1=[]
	ni2=[]
	ni3=[]
	ni4=[]
	
#*********READING AERONET SITE-SPECIFIC FILES*********************
	fpath='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/L2V3/' #other aerosol properties, *.all ext
	flagv=1
	
	for f in sorted(glob.glob(fpath+'*.all')):
		print(f)		 
		data=parse_aero(f,flagv)
	
			
			


