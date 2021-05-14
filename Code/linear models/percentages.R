---
  title: "Gussie's LME"
author: "Gussie and Anshuman"
date: "3/2/2021"
output: html_document


#MixedEffectsModel

library(ggplot2)
library(bootstrap)
library(fossil)
library(iNEXT)
library(plyr)
library(lme4)
library(lmerTest)
library(rsq)

setwd("/Users/gussiemaccracken/Desktop/Review Paper- Ellen Lauren Anshuman")

x = read.csv("Full database.csv", header=T)
head(x)
summary(x)

x$DTs.at.300<-as.numeric(x$DTs.at.300)
x$MAT.C<-as.numeric(x$MAT.C)
x$MAP.mm.yr<-as.numeric(x$MAP.mm.yr)
x$PubDate<-as.numeric(x$PubDate)
x$Shannon<-as.numeric(x$Shannon)
x$Pielou.s.J<-as.numeric(x$Pielou.s.J)
x$Plant.div.300.leaves<-as.numeric(x$Plant.div.300.leaves)
x$Spec.DTs.at.300<-as.numeric(x$Spec.DTs.at.300)
x$Mine.DTs.at.300<-as.numeric(x$Mine.DTs.at.300)
x$Gall.DTs.at.300<-as.numeric(x$Gall.DTs.at.300)
x$Perc.Damage<-as.numeric(x$Perc.Damage)
x$Perc.Galls<-as.numeric(x$Perc.Galls)
x$Perc.Mines<-as.numeric(x$Perc.Mines)
x$Perc.Spec.Dam<-as.numeric(x$Perc.Spec.Dam)
x$Perc.Spec.DTOs<-as.numeric(x$Perc.Spec.DTOs)


#SINGLE LINNEAR MODELS
#1 How does MAT affect DT richness?
#Plot the data
hist(x$Perc.Damage)
plot(x$MAT.C,x$Perc.Damage, pch=21, bg="lightblue", 
     xlab="Mean Annual Temperature (C)", ylab="Perc.Damage", bty="l")

#Fit a simple linear model to explore the relationship between MAT and DTs
lmod1<-glm(Perc.Damage~MAT.C, family=gaussian(link="identity"),dat=x)
coef(lmod1)
abline(lmod1)

#Review the results
summary(lmod1) 

# Significant

#2 How does precipitation affect DT richness?
plot(x$MAP.mm.yr,x$Perc.Damage, pch=21, bg="lightblue", 
     xlab="Mean Annual Precip (C)", ylab="Perc.Damage", bty="l")
lmod2<-glm(Perc.Damage~MAP.mm.yr, family=gaussian(link="identity"),dat=x)
coef(lmod2)
abline(lmod2)
summary(lmod2) 
#Not significant

#3 How does date of publication affect DT richness?
plot(x$PubDate,x$Perc.Damage, pch=21, bg="#a6808c", 
     xlab="Publication Date", ylab="Perc.Damage", bty="l")
lmod3<-glm(Perc.Damage~PubDate, family=gaussian(link="identity"),dat=x)
coef(lmod3)
summary(lmod3) 
#Not significant

#4 How does geologic age affect DT richness?
plot(x$Age,x$Perc.Damage, pch=21, bg="#177e89", 
     xlab="Age (MY)", ylab="Perc.Damage", bty="l")
lmod4<-glm(Perc.Damage~Age, family=gaussian(link="identity"),dat=x)
coef(lmod4)
summary(lmod4) 
#Not Significant

#5 How does plant diversity (Shannon's) affect DT richness?
plot(x$Shannon,x$Perc.Damage, pch=21, bg="lightblue", 
     xlab="Plant Div Shannon", ylab="Perc.Damage", bty="l")
lmod5<-glm(Perc.Damage~Shannon, family=gaussian(link="identity"),dat=x)
coef(lmod5)
abline(lmod5)
summary(lmod5) 
#NOT Significant

#6 How does plant diversity (Pielou's) affect DT richness?
plot(x$Pielou.s.J,x$Perc.Damage, pch=21, bg="lightblue", 
     xlab="Plant Div Pielou", ylab="Perc.Damage", bty="l")
lmod6<-glm(Perc.Damage~Pielou.s.J, family=gaussian(link="identity"),dat=x)
coef(lmod6)
abline(lmod6)
summary(lmod6) 
#Not significant


#7 How does plant diversity (at 300) affect DT richness?
plot(x$Plant.div.300.leaves,x$Perc.Damage, pch=21, bg="lightblue", 
     xlab="Plant div at 300", ylab="Perc.Damage", bty="l")
lmod7<-glm(Perc.Damage~Plant.div.300.leaves, family=gaussian(link="identity"),dat=x)
coef(lmod7)
abline(lmod7)
summary(lmod7) 
#Not significant


#8 How does MAT affect specialist DT richness?
hist(x$Perc.Spec.Dam)
plot(x$MAT.C,x$Perc.Spec.Dam, pch=21, bg="lightblue", 
     xlab="Mean Annual Temperature (C)", ylab="Perc.Spec.Dam", bty="l")
lmod8<-glm(Spec.DTs.at.300~MAT.C, family=gaussian(link="identity"),dat=x)
coef(lmod8)
abline(lmod8)
summary(lmod8) 
#Significant


#9 How does MAP affect Spec DT richness?
hist(x$Perc.Spec.Dam)
plot(x$MAP.mm.yr,x$Perc.Spec.Dam, pch=21, bg="lightblue", 
     xlab="Mean Annual precip (mm/yr)", ylab="Perc.Spec.Dam", bty="l")
lmod10<-glm(Perc.Spec.Dam~MAP.mm.yr, family=gaussian(link="identity"),dat=x)
coef(lmod10)
abline(lmod10)
summary(lmod10) 
#Not significant


#10 How does date of publication affect SPEC DT richness?
plot(x$PubDate,x$Perc.Spec.Dam, pch=21, bg="#a6808c", 
     xlab="Publication Date", ylab="Perc.Spec.Dam", bty="l")
lmod3<-glm(Perc.Spec.Dam~PubDate, family=gaussian(link="identity"),dat=x)
coef(lmod3)
summary(lmod3) 
#Not significant

#11 How does geologic age affect SPEC DT richness?
plot(x$Age,x$Perc.Spec.Dam, pch=21, bg="#177e89", 
     xlab="Age (MY)", ylab="Perc.Spec.Dam", bty="l")
lmod4<-glm(Perc.Spec.Dam~Age, family=gaussian(link="identity"),dat=x)
coef(lmod4)
summary(lmod4) 
#Not Significant

#12 How does plant diversity (Shannon's) affect SPEC DT richness?
plot(x$Shannon,x$Perc.Spec.Dam, pch=21, bg="lightblue", 
     xlab="Plant Div Shannon", ylab="Perc.Spec.Dam", bty="l")
lmod5<-glm(Perc.Spec.Dam~Shannon, family=gaussian(link="identity"),dat=x)
coef(lmod5)
abline(lmod5)
summary(lmod5) 
#Not Significant

#13 How does plant diversity (Pielou's) affect SPEC DT richness?
plot(x$Pielou.s.J,x$Perc.Spec.Dam, pch=21, bg="lightblue", 
     xlab="Plant Div Pielou", ylab="Perc.Spec.Dam", bty="l")
lmod6<-glm(Perc.Spec.Dam~Pielou.s.J, family=gaussian(link="identity"),dat=x)
coef(lmod6)
abline(lmod6)
summary(lmod6) 
#Not significant

#14 How does plant diversity (at 300) affect SPEC DT richness?
plot(x$Plant.div.300.leaves,x$Perc.Spec.Dam, pch=21, bg="lightblue", 
     xlab="Plant div at 300", ylab="Perc.Spec.Dam", bty="l")
lmod7<-glm(Perc.Spec.Dam~Plant.div.300.leaves, family=gaussian(link="identity"),dat=x)
coef(lmod7)
abline(lmod7)
summary(lmod7) 
#Not significant


#15 How does MAT affect Galling DT richness?
hist(x$Gall.DTs.at.300)
plot(x$MAT.C,x$Perc.Galls, pch=21, bg="lightblue", 
     xlab="Mean Annual T", ylab=" Perc.Galls", bty="l")
lmod12<-glm(Perc.Galls~MAT.C, family=gaussian(link="identity"),dat=x)
coef(lmod12)
abline(lmod12)
summary(lmod12) 
# NOT Significant


#16 How does MAP affect GALL DT richness?
hist(x$Perc.Galls)
plot(x$MAP.mm.yr,x$Perc.Galls, pch=21, bg="lightblue", 
     xlab="Mean Annual Temperature (C)", ylab=" Perc.Galls", bty="l")
lmod9<-glm(Perc.Galls~MAP.mm.yr, family=gaussian(link="identity"),dat=x)
coef(lmod9)
abline(lmod9)
summary(lmod9) 
#Not Significant

#17 How does date of publication affect GALL DT richness?
plot(x$PubDate,x$Perc.Galls, pch=21, bg="#a6808c", 
     xlab="Publication Date", ylab="Perc.Galls", bty="l")
lmod3<-glm(Perc.Galls~PubDate, family=gaussian(link="identity"),dat=x)
coef(lmod3)
summary(lmod3) 
# Significant

#18 How does geologic age affect GALL DT richness?
plot(x$Age,x$Perc.Galls, pch=21, bg="#177e89", 
     xlab="Age (MY)", ylab="Perc.Galls", bty="l")
lmod4<-glm(Perc.Galls~Age, family=gaussian(link="identity"),dat=x)
coef(lmod4)
summary(lmod4) 
#Not Significant

#19 How does plant diversity (Shannon's) affect GALL DT richness?
plot(x$Shannon,x$Perc.Galls, pch=21, bg="lightblue", 
     xlab="Plant Div Shannon", ylab="Perc.Galls", bty="l")
lmod5<-glm(Perc.Galls~Shannon, family=gaussian(link="identity"),dat=x)
coef(lmod5)
abline(lmod5)
summary(lmod5) 
# NOT significant

#20 How does plant diversity (Pielou's) affect GALL DT richness?
plot(x$Pielou.s.J,x$Perc.Galls, pch=21, bg="lightblue", 
     xlab="Plant Div Pielou", ylab="Perc.Galls", bty="l")
lmod6<-glm(Perc.Galls~Pielou.s.J, family=gaussian(link="identity"),dat=x)
coef(lmod6)
abline(lmod6)
summary(lmod6) 
#Not significant

#21 How does plant diversity (at 300) affect GALL DT richness?
plot(x$Plant.div.300.leaves,x$Perc.Galls, pch=21, bg="lightblue", 
     xlab="Plant div at 300", ylab="Perc.Galls", bty="l")
lmod7<-glm(Perc.Galls~Plant.div.300.leaves, family=gaussian(link="identity"),dat=x)
coef(lmod7)
abline(lmod7)
summary(lmod7) 
#Not significant

#22 How does MAT affect Mining DT richness?
hist(x$Perc.Mines)
plot(x$MAT.C,x$Perc.Mines, pch=21, bg="lightblue", 
     xlab="Mean Annual T", ylab="Perc.Mines", bty="l")
lmod13<-glm(Perc.Mines~MAT.C, family=gaussian(link="identity"),dat=x)
coef(lmod13)
abline(lmod13)
summary(lmod13) 
#Significant


#23 How does MAP affect mining MINE DT richness?
hist(x$Perc.Mines)
plot(x$MAP.mm.yr,x$Perc.Mines, pch=21, bg="lightblue", 
     xlab="Mean Annual precip (mm/yr)", ylab="Perc.Miness", bty="l")
lmod10<-glm(Perc.Mines~MAP.mm.yr, family=gaussian(link="identity"),dat=x)
coef(lmod10)
abline(lmod10)
summary(lmod10) 
#Not significant

#24 How does date of publication affect MINE DT richness?
plot(x$PubDate,x$Perc.Mines, pch=21, bg="#a6808c", 
     xlab="Publication Date", ylab="Perc.Mines", bty="l")
lmod3<-glm(Perc.Mines~PubDate, family=gaussian(link="identity"),dat=x)
coef(lmod3)
summary(lmod3) 
#Not significant

#25 How does geologic age affect MINE DT richness?
plot(x$Age,x$Perc.Mines, pch=21, bg="#177e89", 
     xlab="Age (MY)", ylab="Perc.Mines", bty="l")
lmod4<-glm(Perc.Mines~Age, family=gaussian(link="identity"),dat=x)
coef(lmod4)
summary(lmod4) 
#Not Significant

#26 How does plant diversity (Shannon's) affect MINE DT richness?
plot(x$Shannon,x$Perc.Mines, pch=21, bg="lightblue", 
     xlab="Plant Div Shannon", ylab="Perc.Mines", bty="l")
lmod5<-glm(Perc.Mines~Shannon, family=gaussian(link="identity"),dat=x)
coef(lmod5)
abline(lmod5)
summary(lmod5) 
#significant

#27 How does plant diversity (Pielou's) affect MINE DT richness?
plot(x$Pielou.s.J,x$Perc.Mines, pch=21, bg="lightblue", 
     xlab="Plant Div Pielou", ylab="Perc.Mines", bty="l")
lmod6<-glm(Perc.Mines~Pielou.s.J, family=gaussian(link="identity"),dat=x)
coef(lmod6)
abline(lmod6)
summary(lmod6) 
#Not significant

#28 How does plant diversity (at 300) affect MINE DT richness?
plot(x$Plant.div.300.leaves,x$Perc.Mines, pch=21, bg="lightblue", 
     xlab="Plant div at 300", ylab="Perc.Mines", bty="l")
lmod7<-glm(Perc.Mines~Plant.div.300.leaves, family=gaussian(link="identity"),dat=x)
coef(lmod7)
abline(lmod7)
summary(lmod7) 
#Not significant