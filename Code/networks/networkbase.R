#constructing bipartite networks and calculating network metrics for all sites
filenets = list.files(path="data/Localities", pattern="*.csv", full.names=T) #dir for all sites csv data
filenets2 = list.files(path="data/Localities", pattern="*.csv", full.names=F) #dir for all sites csv data

library(dplyr)
library(bipartite)

networklev=NULL
higherlev=NULL
plantlev=NULL

for (i in 1:length(filenets)){
  nam=paste(filenets[i])
  netcurr=read.csv(nam, header=T)
  colnames(netcurr)[1]="ID"
  netcurr2=netcurr %>% group_by(ID) %>%  summarise_if(is.numeric, sum, na.rm = TRUE)
  netmat = data.matrix(netcurr2[,-c(1)], rownames.force = NA)
  rownames(netmat)= netcurr2$ID
  
  
  
  nam2=strsplit(filenets2[i],".csv")[[1]]
  plotweb(net) #(to plot the networks)
  networklev=rbind(networklev,c(nam2,networklevel(netmat)))
  X=specieslevel(netmat)
  X1=X$`higher level`
  X2=X$`lower level`
  
  X1$web=nam2
  X2$web=nam2
  
  higherlev=rbind(higherlev,X1)
  plantlev=rbind(plantlev,X2)
  #print(i) #marker
}

