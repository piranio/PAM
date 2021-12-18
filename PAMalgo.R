library(ClusterR)
##*************************************
##Author: Martin A. Montes
##SSAI, April 29, 2020
##*************************************

#replace input traiing file path by your local folder-----
raw = read.csv("/Users/mmontes/Documents/aeronet/globe/scripts/PAM/R/in/trainglobe75_3params2.csv", header = FALSE)

path_r_out='/Users/mmontes/Documents/aeronet/globe/scripts/PAM/R/out/'


myArgs <- commandArgs(trailingOnly = TRUE)
#Convert to numerics
kmaxi = as.numeric(myArgs[1])
flag_run=myArgs[2]
print(kmaxi)
print(flag_run)

dat = center_scale(raw)
cm = Cluster_Medoids(dat, clusters = kmaxi, distance_metric = 'mahalanobis', swap_phase = TRUE)


nums='Clustering processing ----OK!!'
cat(nums)

###---output indices (ci),medoid centers (cc) and silhouette (si) values as csv files---

ci=cm[5]
str1=paste(path_r_out,"clus_id",flag_run,".csv",sep="")
write.csv(ci,str1)
print("saving ci done")

cc=cm[1]
str2=paste(path_r_out,"clus_medoids",flag_run,".csv",sep="")
write.csv(cc,str2)
print("saving cc done")


si=cm[6]
str3=paste(path_r_out,"clus_sil",flag_run,".csv",sep="")
write.csv(si,str3)
print("saving si done")