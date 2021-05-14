#calculating co-occurrence matrix for all sites
filenets = list.files(path="data/Localities", pattern="*.csv", full.names=T)
filenets2 = list.files(path="data/Localities", pattern="*.csv", full.names=F)

codata=NULL

for (i in 1:length(filenets)){
  nam=paste(filenets[i])
  netcurr=read.csv(nam, header=T)
  colnames(netcurr)[1]="ID"
  netcurr2=netcurr %>% group_by(ID) %>%  summarise_if(is.numeric, sum, na.rm = TRUE)
  netmat = data.matrix(netcurr2[,-c(1)], rownames.force = NA)
  rownames(netmat)= netcurr2$ID
  
  netmat=as.data.frame(netmat)
  
  #print(ncol(netmat))
  
  # net=net[,-1]
  nam2=strsplit(filenets2[i],".csv")[[1]]
  
  netmat$web=nam2
  if (ncol(netmat)==224){
    codata=rbind(codata,netmat)
  }
  print(i)
}

alldata2=codata[,1:223]
#rownames(alldata2)=alldata[,1]

alldata3=(t(alldata2))
alldata3[alldata3>0]=1

rownames(alldata3)=colnames(alldata2)



library(cooccur)

paleo3DT.co=cooccur(mat=alldata3,type = "spp_site", thresh = TRUE, true_rand_classifier = 0.1, spp_names = TRUE)

library(ggplot2)
library(ggthemes)
expobs=function (mod) 
{
  ptab <- mod$results
  ptab$signs <- ifelse(ptab$p_gt >= 0.05, 0, 1) + ifelse(ptab$p_lt >= 
                                                           0.05, 0, -1)
  exp_cooccur <- ptab$exp_cooccur
  obs_cooccur <- ptab$obs_cooccur
  signs <- ptab$signs
  p <- ggplot(ptab, aes(x = exp_cooccur, y = obs_cooccur)) + 
    geom_point(aes(fill = factor(signs, levels = c(-1, 0, 
                                                   1))),  pch = 21, size = 4)+theme_few()
  p <- p + scale_fill_manual(values = c("lightcoral", "honeydew", 
                                        "greenyellow"), name = "", labels = c("negative", "random", 
                                                                              "positive"), drop = FALSE)
  p <- p + theme(plot.title = element_text(vjust = 12, size = 15, 
                                           face = "bold"), legend.text = element_text(size = 12), 
                 axis.title = element_text(size = 15), axis.text = element_text(size = 15), 
                 axis.text.x = element_text(hjust = 0.5, vjust = 0)) + xlab("Expected Co-occurrences") + 
    ylab("Observed Co-occurrences")
  p <- p + ggtitle("Observed-Expected Plot") + geom_abline(color = "dark gray")
  p
}

expobs(paleo3DT.co)


files=paleo3DT.co[[2]]
n=223
vals=ifelse(files$p_gt >= 0.05, 0, 1)+ifelse(files$p_lt >= 0.05, 0, -1)
p3DTmatrix=matrix(0,n,n)
for (z in 1:nrow(files)){
  i=files$sp1[z]
  j=files$sp2[z]
  p3DTmatrix[i,j]=vals[z]
  p3DTmatrix[j,i]=p3DTmatrix[i,j]
}

#contdata=as.data.frame(contmatrix)

colnames(p3DTmatrix)=colnames(codata[,1:223])
rownames(p3DTmatrix)=colnames(p3DTmatrix)

datum=as.data.frame(cbind(files$sp1,files$sp2, vals))

mat2=abs(p3DTmatrix)

ss=which(rowSums(mat2)==0)

p3DTmatF=p3DTmatrix[-ss,-ss]

heatmap(p3DTmatF, symm = T, col = c("honeydew", "yellowgreen"))



filenets = list.files(path="data/Localities", pattern="*.csv", full.names=T)
filenets2 = list.files(path="data/Localities", pattern="*.csv", full.names=F)

coDTs=list()
allDTmat=matrix(0,223,223)
rownames(allDTmat)=colnames(alldata2)
colnames(allDTmat)=rownames(allDTmat)

for (k in 57:58){
  nam=paste(filenets[k])
  netcurr=read.csv(nam, header=T)
  colnames(netcurr)[1]="ID"
  
  alldata2=netcurr[,2:224]
  #rownames(alldata2)=alldata[,1]
  
  alldata3=(t(alldata2))
  alldata3[alldata3>0]=1
  
  midDT.co=cooccur(mat=alldata3,type = "spp_site", thresh = TRUE, true_rand_classifier = 0.1, spp_names = TRUE)
  
  coDTs[[k]]=midDT.co
  
  files=midDT.co[[2]]
  n=223
  vals=ifelse(files$p_gt >= 0.05, 0, 1)+ifelse(files$p_lt >= 0.05, 0, -1)
  newmatrix=matrix(0,n,n)
  for (z in 1:nrow(files)){
    i=files$sp1[z]
    j=files$sp2[z]
    newmatrix[i,j]=vals[z]
    newmatrix[j,i]=newmatrix[i,j]
  }
  
  allDTmat=allDTmat+newmatrix
  
  print(k)
}


