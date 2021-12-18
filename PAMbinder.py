import subprocess
import sys
import pandas as pd
import numpy as np

##*************************************
##Author: Martin A. Montes
##SSAI, April 29, 2020
##*************************************




def check_kmax(nclus,sil_thre,nsamples):
##Method to compute maximum number of clusters based on Silhouette index and number of records per cluster

#open file with the silhouette results--------
#for each initial cluster check if we have a cluster based on the sil_thre and number of samples 
	#Kaufman and Rousseux 1990 criterion--
	#-1 to 0.25 (no substantial structure)
	#0.26-0.5 (weak structure)
	#0.51-0.70 (reasonable structure)
	#>0.71 (strong structure))
	
	path_in='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/R/out/'
	filename='clus_silinit.csv'
	fullname=path_in+filename
	sildata = pd.read_csv(fullname)
	
	##--retrieving the Silhouette values
	sildata=sildata.values
	arr=sildata[:,3:5]
	arrmax=[]
	
	for i in arr:
		arrmax.append(np.max(i))
		
	
	arrmax=np.asarray(arrmax,dtype=float)
	###--Compute Silhouette ------------------------
	sil=(sildata[:,4]-sildata[:,3])/arrmax
	####----Get values with the correct silhouette
	ix=np.where(sil>=sil_thre) #all clusters with correct silouette
	
	dataok=sildata[ix,:]
	clusid=sildata[ix,1]	
	
	oix=[]
	oix2=[]
	####----for each clusterid check if number of records is >10
	for j in range(1,nclus+1):

		ix=np.where(clusid==j)
		ix=np.asarray(ix,dtype=int)

		row,col=ix.shape
		oix.append(col)
	
	#do we have a cluster?----------
	nc=0
	for i in oix:
		if i>=nsamples:
			#count that partition as a cluster
			nc+=1
	kmax_opt=nc
	return kmax_opt

	

def pam_clus2(kmaxi,flag_run):
####--Binder calling cluster algorithm in R----------
	command = 'Rscript'
	path2script = '/Users/mmontes/Documents/aeronet/globe/scripts/PAM/PAMalgo.R'
	
	if flag_run==1:
		flag_run='init'
	elif flag_run==2:
		flag_run='optm'
	
	kmaxi=str(kmaxi)
	#print(kmaxi)
	#args=list(kmaxi.split())
	args1=list(kmaxi.split())		
	args2=list(flag_run.split())
	args=[args1,args2]
	
# Build subprocess command--first argument is kmaxi or kmax and second argument if the processing flag
	cmd = [command, path2script] + args1 + args2
	subprocess.Popen(cmd)
	
#Output results
	outvar = subprocess.check_output(cmd)
	print('done running PAM clustering:', outvar)