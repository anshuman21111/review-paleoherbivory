##Calculation of network properties with bootstrapping at 300 leaves without replacement

#names of files -- put the location of the csv files here
filenets = list.files(path="data/Localities_3/select", pattern="*.csv", full.names=T)
filenets2 = list.files(path="data/Localities_3/select", pattern="*.csv", full.names=F)

#seeding
set.seed(1)

#initializing lists to measure properties and their variances
networklevlist=NULL
networklevmean=NULL
networklevelvar=NULL

library(Rfast)
library(bipartite)

for (i in 1:length(filenets)){
  nam=paste(filenets[i])
  netcurr=read.csv(nam, header=T)
  colnames(netcurr)[1]="ID"
  netcurr2=netcurr %>% group_by(ID) %>%  summarise_if(is.numeric, sum, na.rm = TRUE)
  netmat = data.matrix(netcurr2[,-c(1)], rownames.force = NA)
  rownames(netmat)= netcurr2$ID
  
  minrow=300 #minimum number of leaves
  
  netpro=NULL
  allnames=cbind(rownames(netmat),0)
  colnames(allnames)=c("ID","Nul")
  allnames=as.data.frame(allnames)
  
  for (j in 1:500){
    q2=sample(nrow(netcurr), minrow, replace=F)
    datarand=netcurr[q2,]
    netrand=datarand %>% group_by(ID) %>%  summarise_if(is.numeric, sum, na.rm = TRUE)
    netrandfinal=left_join(allnames,netrand, by="ID")
    netrandfinal[is.na(netrandfinal)] <- 0
    netrandfinal=netrandfinal[,-c(1:2)]
    netrandfinal=data.matrix(netrandfinal)
    netpro=rbind(netpro,networklevel(netrandfinal))
    
    if (j %% 100 == 0){print(paste("Sequence", " ", j))}
  }
  
  
  
  nam2=strsplit(filenets2[i],".csv")[[1]]
  
  networklevmean=rbind(networklevmean,c(nam2,colmeans(netpro)))
  networklevelvar=rbind(networklevelvar,c(nam2,colVars(netpro, std = F)))
  
  networklevlist[[i]]=netpro
  
  print(i)
}

colnames(networklevmean)=c("Web",colnames(netpro))
colnames(networklevelvar)=c("Web",colnames(netpro))
