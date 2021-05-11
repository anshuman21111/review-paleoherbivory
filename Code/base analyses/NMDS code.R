#NMDS plots for the Paleoherbivory Project
#Author: EDC
#Date: February 2021

#load packages
library(ggplot2)
library(tidyverse)
library(MASS)
library(vegan)

#load data file
full_data <- read.csv("Full database.csv", header=T, sep=",", 
                         na.strings=  "NA", dec=".", strip.white=TRUE)

#separate out FFG data & only include sites with data
ffg.data <- full_data[,c(1:6,18,20,22:26)]
ffg.data <- subset(ffg.data, Perc.HF>0) #eliminates rows with NA

#NMDS of sites by FFGs
#construct a matrix with only the proportional species data for each site
prop.ffg.data <- ffg.data[,7:13]/100
row.names(prop.ffg.data)<- ffg.data$Flora

#arcsin sqrt transform
prop.ffg.asin <- asin(sqrt(prop.ffg.data))

#NMS
prop.ffg.nms <- metaMDS(prop.ffg.asin, distance="bray")
s1 <- scores(prop.ffg.nms, display=c("sites"), choices=1)
s2 <- scores(prop.ffg.nms, display=c("sites"), choices=2)
v1 <- scores(prop.ffg.nms, display=c("species"), choices=1)
v2 <- scores(prop.ffg.nms, display=c("species"), choices=2)

#Write NMS data to a spreadsheet.
nms.coords <- tibble(ffg.data$Flora, ffg.data$Abbreviation, ffg.data$Region, ffg.data$Latitude, ffg.data$Epoch.Period, s1, s2)
write.table(nms.coords,"NMS Coordinates.txt", sep="\t")

#Create Plot of NMDS ordination with points coded by region and sized based on diversity
nms.coords <- read.csv("NMS Coordinates.csv", header=T)
nms.region <- split(nms.coords, nms.coords$Region)
s1.split <- nms.region$s1

plot(nms.coords$s1, nms.coords$s2, type="n", xlab="NMDS 1", ylab="NMDS 2", xlim = c(-1,1), ylim = c(-0.5,0.5))
points(nms.region$Patagonia$s1,nms.region$Patagonia$s2, pch=16, col="#F5DD90", cex=nms.region$Patagonia$NMDS.scale)
points(nms.region$`North America`$s1,nms.region$`North America`$s2, pch=16, col="#177e89", cex=nms.region$`North America`$NMDS.scale)
points(nms.region$Eurasia$s1,nms.region$Eurasia$s2, pch=8, col="#177e89", cex=nms.region$Eurasia$NMDS.scale)
points(nms.region$Africa$s1,nms.region$Africa$s2, pch=8, col="#F76C5E", cex=nms.region$Africa$NMDS.scale)
points(nms.region$`Tropical South America`$s1,nms.region$`Tropical South America`$s2, pch=16, col="#F76C5E", cex=nms.region$`Tropical South America`$NMDS.scale)
points(nms.region$`New Zealand`$s1,nms.region$`New Zealand`$s2, pch=8, col="#F5DD90", cex=nms.region$`New Zealand`$NMDS.scale)
points(nms.region$Iceland$s1,nms.region$Iceland$s2, pch=8, col="#a6808c", cex=nms.region$Iceland$NMDS.scale)
points(nms.region$Spitsbergen$s1,nms.region$Spitsbergen$s2, pch=4, col="#a6808c", cex=nms.region$Spitsbergen$NMDS.scale)
points(nms.region$Alaska$s1,nms.region$Alaska$s2, pch=16, col="#a6808c", cex=nms.region$Alaska$NMDS.scale)
points(nms.region$Antarctica$s1, nms.region$Antarctica$s2, pch=17, col="#586BA4", cex=nms.region$Antarctica$NMDS.scale)
points(v1,v2, col="black", pch=1, cex=1.7)
leg.text <- c("Iceland", "Spitsbergen", "Alaska", "North America", "Eurasia", "Africa", "Tropical South America", "Patagonia", "New Zealand", "Antarctica")
legend("bottomleft", leg.text, col = c("#a6808c","#a6808c","#a6808c", "#177e89", "#177e89", "#F76C5E","#F76C5E", "#F5DD90","#F5DD90","#586BA4"), pch = c(8,4,16,16,8,8,16,16,8,17), cex=1)
title("NMDS of Sites by FFG, coded by region")


anosim.cont <- anosim(prop.ffg.asin, continent, distance="bray")
summary(anosim.cont)