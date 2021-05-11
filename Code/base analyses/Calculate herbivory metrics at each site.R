#Calculation of herbivory metrics from raw data
#Data input matrix: presence-absence of each DT (columns) on each leaf at a site (rows)
#The analytic rarefaction solution for herbivory data was provided by Torsten Wappler and was first described and published in Gunkel & Wappler 2015.

library(vegan)
mask <- read.csv("mask.csv", header=T) #Read mask- this file defines each DT- herbivory, specialized, gall, mine
ffg_def <- read.csv("ffgdef.csv", header=T)$x #file defining FFGs: Boring, Fungal, Seed Predation, Incertae, Oviposition =0, Galling=3, Hole Feeding=4, Mining =5, Margin Feeding =6,  Piercing and Sucking=7, Skeletonization=1, Surface feeding =2, 

################################################################################
#READ IN DATA AND FORMAT FOR ANALYSES
inputfilename <- "filename.csv" #change "filename" to the name of your matrix
dataname <- "dataname" #free to choose 
cutoff <- 100 #cutoff for resampling data for a single species

dataall   <- read.csv(inputfilename, header=T) #read in data
datatable <- dataall[,2:ncol(dataall)]
species   <- dataall[,1]

ffg_data  <- matrix (ncol=max(ffg_def), nrow=nrow(datatable), data=0)
for (i in 1:max(ffg_def)) ffg_data[,i]<-as.numeric(as.logical(rowSums(datatable[,ffg_def==i])))

################################################################################
#CONSTRUCT A MATRIX OF % OF LEAVES WITH EACH FFG
ffg_comp <- matrix(ncol = 7, nrow=1)
colnames(ffg_comp) <- c("Skeletonization", "Surface feeding", "Gall", "Hole Feeding", "Mining", "Margin Feeding", "Piercing")
ffg_comp[1,] <- colSums(ffg_data)/nrow(dataall)*100

################################################################################
#GUNKEL AND WAPPLER 2015 CODE
#calculates total, specialized, mine, and gall richness
#The summary .csv file output includes total and specialized damage frequency
counta <- function (dat, classa) {
 v <- dat[,classa]
 count <- sum(v)
 count}

countb <- function (dat, classa, classb) {
 va <- dat[,classa]
 vb <- dat[,classb]
 va <- 1-va
 vb <- 1-vb
 v  <- va*vb
 v <- 1-v
 count <- sum(v)
count}

rarefy <- function(data, dname) {

if (!(is.vector(data[,colSums(data)>0]))) data <-data[,colSums(data)>0]

N <- nrow(data) # find the number of individuals
m <- ncol(data) # find the number of classes

k <- array(0, c(m))
l <- array(0, c(m,m))

rarefaction <- matrix(data=NA, nrow=N, ncol=3) #This is your final data matrix 
barefaction <- matrix(data=NA, nrow=N, ncol=3) #Intermediate results

colnames(rarefaction) <- c("E(X)", "Var(X)", "SD")

#Main Calculations

lister <- array (0, c(N,m)) #Start: 1st Moment
for (i in 1:m) k[i] <- counta (data, i)
for (j in 1:m) lister[1,j] <- ((N-k[j])/N)
for (s in 2:N) for (j in 1:m) if (N-s-k[j]>0) lister[s,j] <- lister[s-1,j]*(N-s-k[j])/(N-s) else lister [s,j]=0
for (s in 1:N) barefaction [s,1] <- sum(lister[s,])
for (s in 1:N) rarefaction[s,1] <- m-barefaction [s,1] #End: 1st Moment

listerb <- array (0, c(N,m,m)) #Start: 2nd Moment
for (i in 2:m) for (j in 1:(i-1)) l[i,j] <- countb (data, i, j)
for (i in 2:m) for (j in 1:(i-1)) listerb[1,i,j] <- ((N-l[i,j])/N)
for (s in 2:N) for (i in 2:m) for (j in 1:(i-1)) if (N-s-l[i,j]>0) listerb[s,i,j] <- listerb[s-1,i,j]*(N-s-l[i,j])/(N-s) else listerb[s,i,j] <- 0
for (s in 1:N) barefaction [s,2] <- sum(listerb[s,,])
for (s in 1:N) rarefaction[s,2] <- (barefaction [s,1] +2*barefaction[s,2]-barefaction[s,1]**2) #End: 2nd Moment
for (s in 1:N) rarefaction[s,3] <- rarefaction[s,2]**(1/2)

write.table(rarefaction, paste(dname,".csv"), sep=",")
}

rarefy (datatable[,mask$All],paste(dataname,"_All"))
rarefy (datatable[,mask$Spec],paste(dataname,"_Spec"))
rarefy (datatable[,mask$Gall],paste(dataname,"_Gall"))
rarefy (datatable[,mask$Mine],paste(dataname,"_Mine"))
rarefy (ffg_data, paste(dataname,"_FFG"))

specnames<-levels(species)
species<-as.numeric(species)

for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(datatable[(species==i), mask$All])>0) rarefy (datatable[(species==i),mask$All],paste(dataname, "_", specnames[i], "_All"))
for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(datatable[(species==i), mask$Spec])>0) rarefy (datatable[(species==i),mask$Spec],paste(dataname, "_", specnames[i], "_Spec"))
for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(datatable[(species==i), mask$Gall])>0) rarefy (datatable[(species==i),mask$Gall],paste(dataname, "_", specnames[i], "_Gall"))
for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(datatable[(species==i), mask$Mine])>0) rarefy (datatable[(species==i),mask$Mine],paste(dataname, "_", specnames[i], "_Mine"))
for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(ffg_data[species==i,])>0) rarefy (ffg_data[(species==i),], paste(dataname, "_", specnames[i],"_FFG"))

specdata<-matrix(nrow=length(species), ncol=max(species), data=0)
for (i in 1:length(species)) specdata[i,(species[i])]<-1
rarefy (specdata, paste(dataname, "_Species"))

overview<-matrix (nrow=length(specnames)+1, ncol=11)
colnames(overview)<-c("Species", "#leaves", "%DMG", "%Spec", "%Gall", "%Mine", "DTs", "SpecDTs", "GDTs", "MDTs", "#FFGs")
overview[1:length(specnames),1]<-specnames
overview[1:length(specnames),2]<-colSums(specdata)

for (i in 1:length(specnames)) overview[i,3]  <-mean(rowSums(datatable[species==i, mask$All])>0)
for (i in 1:length(specnames)) overview[i,4]  <-mean(rowSums(datatable[species==i, mask$Spec])>0)
for (i in 1:length(specnames)) overview[i,5]  <-mean(rowSums(datatable[species==i, mask$Gall])>0)
for (i in 1:length(specnames)) overview[i,6]  <-mean(rowSums(datatable[species==i, mask$Mine])>0)

for (i in 1:length(specnames)) overview[i,7]  <-sum(colSums(datatable[species==i, mask$All])>0)
for (i in 1:length(specnames)) overview[i,8]  <-sum(colSums(datatable[species==i, mask$Spec])>0)
for (i in 1:length(specnames)) overview[i,9]  <-sum(colSums(datatable[species==i, mask$Gall])>0)
for (i in 1:length(specnames)) overview[i,10] <-sum(colSums(datatable[species==i, mask$Mine])>0)

for (i in 1:length(specnames)) overview[i,11] <-sum(colSums(as.matrix(ffg_data[species==i, ]))>0)

overview[length(specnames)+1,1]<-"Total"
overview[length(specnames)+1,2]<-sum(as.numeric(overview[1:length(specnames),2]))
overview[length(specnames)+1,3]<-mean(rowSums(datatable[,mask$All])>0)
overview[length(specnames)+1,4]<-mean(rowSums(datatable[,mask$Spec])>0)
overview[length(specnames)+1,5]<-mean(rowSums(datatable[,mask$Gall])>0)
overview[length(specnames)+1,6]<-mean(rowSums(datatable[,mask$Mine])>0)
overview[length(specnames)+1,7]<-sum(colSums(datatable[,mask$All])>0)
overview[length(specnames)+1,8]<-sum(colSums(datatable[,mask$Spec])>0)
overview[length(specnames)+1,9]<-sum(colSums(datatable[,mask$Gall])>0)
overview[length(specnames)+1,10]<-sum(colSums(datatable[,mask$Mine])>0)
overview[length(specnames)+1,11]<-sum(colSums(ffg_data)>0)


write.table(overview, paste(dataname,"_summary.csv"), sep=",")

################################################################################
#CALCULATE PLANT DIVERSITY METRICS
PlantCounts <- matrix(data = NA, nrow = nrow(overview)-1, ncol=1)
PlantCounts[,1] <- as.numeric(overview[1:(nrow(overview)-1),2])
shannon <- diversity(PlantCounts, index = "shannon")
J <- (shannon)/log(nrow(PlantCounts))

################################################################################
#Compute the % DTOs that are specialized for all datasets
DTO.data <- colSums(datatable) # finds the number of occurrences of each DT
spec.DTOs <- subset(DTO.data[,c(mask$Spec==T)]) #subset the data so only Spec. DTs are included
perc.spec.DTOs <- sum(spec.DTOs)/sum(DTO.data) #calculate the percent of DT occurrences that are specialized