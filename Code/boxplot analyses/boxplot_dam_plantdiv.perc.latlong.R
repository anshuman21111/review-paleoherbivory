#Paleoherbivory reveiew paper
#Purpose: Damage diveristy and percentage by lat/long
          #plant diversity by lat/long
#Date started: 02.03.2021
#R version:3.6.1
#Script author: LAS

#load packages
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(agricolae)
library(EnvStats)

#load data file----
dam_data <- read.csv("LAS analyses/Paleoherb_paper/allflora.csv", header=T, sep=",", 
                     na.strings=  "NA", dec=".", strip.white=TRUE)
dam_data$Flora <- NULL
#dam_data <-dam_data[-c(64:982), ] 
#for some reason lots of blank rows are being added to the bottom of 
#the data frame. removes them
#orgnize data----
total.perc <- dam_data[,c(2,6)]#percent total damage
total.div <- drop_na(dam_data[,c(2,17)]) #total diversity

#following FFG same organization as above
spec.perc <- drop_na(dam_data[,c(2,7)]) 
spec.div <- drop_na(dam_data[,c(2,18)]) 
mine.perc <- drop_na(dam_data[,c(2,8)])
mine.div <- drop_na(dam_data[,c(2,19)])
gall.perc <- drop_na(dam_data[,c(2,9)])
gall.div <- drop_na(dam_data[,c(2,20)])
HF.perc <- drop_na(dam_data[,c(2,10)])
sk.perc <- drop_na(dam_data[,c(2,12)])
surf.perc <- drop_na(dam_data[,c(2,13)])
marg.perc <- drop_na(dam_data[,c(2,11)])
P.perc <- drop_na(dam_data[,c(2,14)])
plant.div <- drop_na(dam_data[,c(2,22)])
st.dt.ocur <- drop_na(dam_data[,c(2,23)])
#############################################################################
#Cleaning up data and making color palette----
#Ordering legend properly

#Tukey's test for significance across time
#Asking does Xdamage change across time? 
#total percent damage tukey test----
total.lm <- lm(perc.damage ~ Latitude, data = dam_data)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'Latitude')
total_test

ros.totalperc <- rosnerTest(dam_data$perc.damage, k=10)
ros.totalperc
#no outliers 

#Total percent damage
avgtotal.perc <- total.perc %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(Latitude)%>% #grouping by epoch
  summarise(mean.all= mean(perc.damage), #average
            sd.all = sd(perc.damage), #standard deviation
            count=n())

#Bar graph Total %Dam total----
total.perc$Latitude <- factor(total.perc$Latitude,
                                 levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
total.perc <- arrange(total.perc, Latitude)

total.perc.lat <- ggplot(total.perc, aes(fill=Latitude, y=perc.damage, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(avgtotal.perc), label = c("b","a","ab","b","b"), 
            aes(y = c(-0.1,-0.1,-0.1,-0.1,-0.1), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values=col5.bar)+ #color palette created above 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("% Total Damage")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
total.perc.lat

#total damage diversity tukey test----
totaldiv.lm <- lm(stand.dt ~ Latitude, data = dam_data)
totaldiv.av <- aov(totaldiv.lm)
summary(totaldiv.av)
totaldiv_test <- HSD.test(totaldiv.av, 'Latitude')
totaldiv_test

ros.totaldiv <- rosnerTest(total.div$stand.dt, k=3)
ros.totaldiv
#remove outlier
total.div <- total.div[-8,]

#rerun tukey
totaldiv.lm <- lm(stand.dt ~ Latitude, data = dam_data)
totaldiv.av <- aov(totaldiv.lm)
summary(totaldiv.av)
totaldiv_test <- HSD.test(totaldiv.av, 'Latitude')
totaldiv_test

#Total damage diversity 
avgtotal.div <- total.div %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(Latitude)%>% #grouping by epoch
  summarise(mean.all= mean(stand.dt), #average
            sd.all = sd(stand.dt), #standard deviation
            count=n())

#Bar graph Total Damage Diversity----
total.div$Latitude <- factor(total.div$Latitude,
                             levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
total.div <- arrange(total.div, Latitude)

avgtotal.div$Latitude <- factor(avgtotal.div$Latitude,
                             levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
avgtotal.div <- arrange(avgtotal.div, Latitude)

total.div.lat <- ggplot(total.div, aes(fill=Latitude, y=stand.dt, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(avgtotal.div), label = c("a","ab","b","b"), 
            aes(y = c(-0.1,-0.1,-0.1,-0.1), size = 2))+
  scale_fill_manual(values=c("#F5DD90","#F76C5E","#177e89","#a6808c"))+ #color palette created above 
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("Total Damage \n(standarized to 300 leaves)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
total.div.lat

#percent gall damage tukey test----
gall.lm <- lm(perc.gall ~ Latitude, data = dam_data)
gall.av <- aov(gall.lm)
summary(gall.av)
gall_test <- HSD.test(gall.av, trt = 'Latitude')
gall_test

#test for outliers
ros.gallperc <- rosnerTest(gall.perc$perc.gall,
                           k = 10)
ros.gallperc
#remove outliers
gall.perc <- gall.perc[c(-50,-14,-8,-1),]
#rerun tukey without outliers
gall.lm <- lm(perc.gall ~ Latitude, data = dam_data)
gall.av <- aov(gall.lm)
summary(gall.av)
gall_test <- HSD.test(gall.av, trt = 'Latitude')
gall_test

#gall diversity tukey test----
galldiv.lm <- lm(stand.gall ~ Latitude, data = dam_data)
galldiv.av <- aov(galldiv.lm)
summary(galldiv.av)
galldiv_test <- HSD.test(galldiv.av, trt = 'Latitude')
galldiv_test

#test for outliers
ros.galldiv <- rosnerTest(gall.div$stand.gall,
                           k = 10)
ros.galldiv
#remove outliers
gall.div <- gall.div[c(-8,-23),]
#rerun tukey without outliers
galldiv.lm <- lm(stand.gall ~ Latitude, data = dam_data)
galldiv.av <- aov(galldiv.lm)
summary(galldiv.av)
galldiv_test <- HSD.test(galldiv.av, trt = 'Latitude')
galldiv_test

#percent hole damage tukey test----
hole.lm <- lm(perc.hole ~ Latitude, data = HF.perc)
hole.av <- aov(hole.lm)
summary(hole.av)
hole_test <- HSD.test(hole.av, trt = 'Latitude')
hole_test

#test for outliers
ros.hole <- rosnerTest(HF.perc$perc.hole,
                          k = 10)
ros.hole
#no outliers 

#Total hole feeding 
avghole <- HF.perc %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(Latitude)%>% #grouping by epoch
  summarise(mean.all= mean(perc.hole), #average
            sd.all = sd(perc.hole), #standard deviation
            count=n())

#Bar graph Total Damage Diversity----
HF.perc$Latitude <- factor(HF.perc$Latitude,
                             levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
HF.perc <- arrange(HF.perc, Latitude)

avghole$Latitude <- factor(avghole$Latitude,
                                levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
avghole <- arrange(avghole, Latitude)

hole.lat <- ggplot(HF.perc, aes(fill=Latitude, y=perc.hole, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(avghole), 
            label = c("b","a","ab","ab","b"), 
            aes(y = c(-0.1,-0.1,-0.1,-0.1,-0.1), size = 2))+
  scale_fill_manual(values=col5.bar)+ #color palette created above 
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("% Hole Feeding Damage")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
hole.lat

#percent specilized damage tukey test----
spec.lm <- lm(perc.spec ~ Latitude, data = spec.perc)
spec.av <- aov(spec.lm)
summary(spec.av)
spec_test <- HSD.test(spec.av, trt = 'Latitude')
spec_test

ros.spec <- rosnerTest(spec.perc$perc.spec,
                       k = 10)
ros.spec
#remove outlier
spec.perc <- spec.perc[-8,]

spec.lm <- lm(perc.spec ~ Latitude, data = spec.perc)
spec.av <- aov(spec.lm)
summary(spec.av)
spec_test <- HSD.test(spec.av, trt = 'Latitude')
spec_test

#percent specialized damage
avgspec.perc <- spec.perc %>%
  drop_na(perc.spec)%>%
  group_by(Latitude)%>%
  summarise(mean.spec = mean(perc.spec),
            sd.spec = sd(perc.spec),
            count=n()) #number of samples/sites per epoch)

#% specialized dam ----
spec.perc$Latitude <- factor(spec.perc$Latitude,
                              levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
spec.perc <- arrange(spec.perc, Latitude)

avgspec.perc$Latitude <- factor(avgspec.perc$Latitude,
                             levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
avgspec.perc <- arrange(avgspec.perc, Latitude)

spec.perc.lat <- ggplot(spec.perc, aes(fill=Latitude, y=perc.spec, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(avgspec.perc), label = c("a","ab","ab","b"), 
            aes(y = c(-0.1,-0.1,-0.1,-0.1), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values=c("#F5DD90","#F76C5E","#177e89","#a6808c"))+ #color palette created above 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("% Specialized Feeding")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") 
spec.perc.lat

#specialized damage diversity tukey test----
specdiv.lm <- lm(stand.spec ~ Latitude, data = spec.div)
specdiv.av <- aov(specdiv.lm)
summary(specdiv.av)
specdiv_test <- HSD.test(specdiv.av, trt = 'Latitude')
specdiv_test

#test for outliers
ros.specdiv <- rosnerTest(spec.div$stand.spec,
                       k = 10)
ros.specdiv
#remove outlier
spec.div <- spec.div[-8,]
specdiv.lm <- lm(stand.spec ~ Latitude, data = spec.div)
specdiv.av <- aov(specdiv.lm)
summary(specdiv.av)
specdiv_test <- HSD.test(specdiv.av, trt = 'Latitude')
specdiv_test

#specialized damage diversity 
avgspec.div <- spec.div %>%
  drop_na(stand.spec)%>%
  group_by(Latitude)%>%
  summarise(mean.spec = mean(stand.spec),
            sd.spec = sd(stand.spec),
            count=n()) #number of samples/sites per epoch)

#Specialized damage diveristy
spec.div$Latitude <- factor(spec.div$Latitude,
                                levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
spec.div <- arrange(spec.div, Latitude)

avgspec.div$Latitude <- factor(avgspec.div$Latitude,
                            levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
avgspec.div <- arrange(avgspec.div, Latitude)

spec.div.lat <- ggplot(spec.div, aes(fill=Latitude, y=stand.spec, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(avgspec.div), label = c("a","ab","b","b"), 
            aes(y = c(-0.1,-0.1,-0.1,-0.1), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values=c("#F5DD90","#F76C5E","#177e89","#a6808c"))+ #color palette created above 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("Specialized Damage \n(standarized to 300 leaves)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") 
spec.div.lat

#percent mine damage tukey test----
mine.lm <- lm(perc.mine ~ Latitude, data = mine.perc)
mine.av <- aov(mine.lm)
summary(mine.av)
mine_test <- HSD.test(mine.av, trt = 'Latitude')
mine_test

#test for outliers
ros.mineperc <- rosnerTest(mine.perc$perc.mine,
                          k = 10)
ros.mineperc
#remove outliers
mine.perc <- mine.perc[c(-8,-61),]
mine.lm <- lm(perc.mine ~ Latitude, data = mine.perc)
mine.av <- aov(mine.lm)
summary(mine.av)
mine_test <- HSD.test(mine.av, trt = 'Latitude')
mine_test
#removing outliers makes mine percentage insignificant across lats. 

#mine diversity tukey tes----
minediv.lm <- lm(stand.mine ~ Latitude, data = mine.div)
minediv.av <- aov(minediv.lm)
summary(minediv.av)
minediv_test <- HSD.test(minediv.av, trt = 'Latitude')
minediv_test

#check for outliers
ros.minediv <- rosnerTest(mine.div$stand.mine,
                           k = 10)
ros.minediv
mine.div <- mine.div[-8,]
minediv.lm <- lm(stand.mine ~ Latitude, data = mine.div)
minediv.av <- aov(minediv.lm)
summary(minediv.av)
minediv_test <- HSD.test(minediv.av, trt = 'Latitude')
minediv_test

#skeletonization tukey test----
sk.lm <- lm(perc.skel ~ Latitude, data = sk.perc)
sk.av <- aov(sk.lm)
summary(sk.av)
sk_test <- HSD.test(sk.av, trt = 'Latitude')
sk_test

#test for outliers
ros.skel <- rosnerTest(sk.perc$perc.skel,
                          k = 10)
ros.skel
sk.perc <- sk.perc[c(-35,-33,-49,-34,-3),]
sk.lm <- lm(perc.skel ~ Latitude, data = sk.perc)
sk.av <- aov(sk.lm)
summary(sk.av)
sk_test <- HSD.test(sk.av, trt = 'Latitude')
sk_test

#piercing/sucking tukey test----
p.lm <- lm(perc.pierce ~ Latitude, data = P.perc)
p.av <- aov(p.lm)
summary(p.av)
p_test <- HSD.test(p.av, trt = 'Latitude')
p_test

#outliers
ros.p <- rosnerTest(P.perc$perc.pierce,
                       k = 10)
ros.p
P.perc <- P.perc[c(-51,-56,-41,-63,-34,-58,-8,-50,-37,-59),]
ros.p <- rosnerTest(P.perc$perc.pierce,
                    k = 10)
ros.p

p.lm <- lm(perc.pierce ~ Latitude, data = P.perc)
p.av <- aov(p.lm)
summary(p.av)
p_test <- HSD.test(p.av, trt = 'Latitude')
p_test

#percent piercing/sucking damage
avgP.perc <- P.perc %>%
  drop_na(perc.pierce)%>%
  group_by(Latitude)%>%
  summarise(mean.P = mean(perc.pierce),
            sd.P = sd(perc.pierce),
            count=n()) #number of samples/sites per epoch)

#% piercing damage----
P.perc$Latitude <- factor(P.perc$Latitude,
                             levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
P.perc <- arrange(P.perc, Latitude)

avgP.perc$Latitude <- factor(avgP.perc$Latitude,
                          levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
avgP.perc <- arrange(avgP.perc, Latitude)

P.perc.lat <- ggplot(P.perc, aes(fill=Latitude, y=perc.pierce, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(avgP.perc), label = c("b","ab","a","ab","b"), 
            aes(y = c(-0.1,-0.1,-0.1,-0.1,-0.1), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values=col5.bar)+ #color palette created above 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("% Piercing and \nSucking Feeding")+ #axis labels
  xlab("")+
  theme_light()+ #light ackground makes bar plot easier to read
  theme(legend.position = "none") 
P.perc.lat

#perc. margin tukey test----
marg.lm <- lm(perc.marg ~ Latitude, data = marg.perc)
marg.av <- aov(marg.lm)
summary(marg.av)
m_test <- HSD.test(marg.av, trt = 'Latitude')
m_test

#outliers
ros.marg <- rosnerTest(marg.perc$perc.marg,
                    k = 10)
ros.marg

#percent margin damage
avgmarg.perc <- marg.perc %>%
  drop_na(perc.marg)%>%
  group_by(Latitude)%>%
  summarise(mean.P = mean(perc.marg),
            sd.P = sd(perc.marg),
            count=n()) #number of samples/sites per epoch)

#% margin damage----
marg.perc$Latitude <- factor(marg.perc$Latitude,
                          levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
marg.perc <- arrange(marg.perc, Latitude)

avgmarg.perc$Latitude <- factor(avgmarg.perc$Latitude,
                             levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
avgmarg.perc <- arrange(avgmarg.perc, Latitude)

marg.perc.lat <- ggplot(marg.perc, aes(fill=Latitude, y=perc.marg, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(avgmarg.perc), label = c("b","a","ab","b","b"), 
            aes(y = c(-0.1,-0.1,-0.1,-0.1,-0.1), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values=col5.bar)+ #color palette created above 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("% Margine Feeding")+ #axis labels
  xlab("")+
  theme_light()+ #light ackground makes bar plot easier to read
  theme(legend.position = "none") 
marg.perc.lat

#perc. surface tukey test----
surf.lm <- lm(perc.surface ~ Latitude, data = surf.perc)
surf.av <- aov(surf.lm)
summary(surf.av)
srf_test <- HSD.test(surf.av, trt = 'Latitude')
srf_test

ros.surf <- rosnerTest(surf.perc$perc.surface,
                       k = 10)
ros.surf
surf.perc <- surf.perc[c(-8,-11,-61,-54,-53,-39,-57),]
surf.lm <- lm(perc.surface ~ Latitude, data = surf.perc)
surf.av <- aov(surf.lm)
summary(surf.av)
srf_test <- HSD.test(surf.av, trt = 'Latitude')
srf_test

#percent surface damage
avgsurf.perc <- surf.perc %>%
  drop_na(perc.surface)%>%
  group_by(Latitude)%>%
  summarise(mean.P = mean(perc.surface),
            sd.P = sd(perc.surface),
            count=n()) #number of samples/sites per epoch)

#% surface damage----
surf.perc$Latitude <- factor(surf.perc$Latitude,
                          levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
surf.perc <- arrange(surf.perc, Latitude)

avgsurf.perc$Latitude <- factor(avgsurf.perc$Latitude,
                             levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
avgsurf.perc <- arrange(avgsurf.perc, Latitude)

surf.perc.lat <- ggplot(surf.perc, aes(fill=Latitude, y=perc.surface, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(avgsurf.perc), label = c("b","b","a","b","b"), 
            aes(y = c(-0.1,-0.1,-0.1,-0.1,-0.1), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values=col5.bar)+ #color palette created above 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("% Surface Feeding")+ #axis labels
  xlab("")+
  theme_light()+ #light ackground makes bar plot easier to read
  theme(legend.position = "none") 
surf.perc.lat

#plant diversity tukey test----
plant.lm <- lm(stand.plant.div ~ Latitude, data= plant.div)
plant.av <- aov(plant.lm)
summary(plant.av)
plant_test <- HSD.test(plant.av, trt="Latitude")
plant_test

ros.plant <- rosnerTest(plant.div$stand.plant.div,
                       k = 10)
ros.plant
plant.div <- plant.div[c(-53,-3),]
plant.lm <- lm(stand.plant.div ~ Latitude, data= plant.div)
plant.av <- aov(plant.lm)
summary(plant.av)
plant_test <- HSD.test(plant.av, trt="Latitude")
plant_test

#percent of specialized dt occurances----
occur.lm <- lm(Perc.Spec.DTOs ~ Latitude, data= st.dt.ocur)
occur.av <- aov(occur.lm)
summary(occur.av)
occur_test <- HSD.test(occur.av, trt="Latitude")
occur_test

ros.occur <- rosnerTest(st.dt.ocur$Perc.Spec.DTOs,
                        k = 10)
ros.occur
st.dt.ocur <- st.dt.ocur[-46,]
occur.lm <- lm(Perc.Spec.DTOs ~ Latitude, data= st.dt.ocur)
occur.av <- aov(occur.lm)
summary(occur.av)
occur_test <- HSD.test(occur.av, trt="Latitude")
occur_test

#putting plots together into superplot----
super.spec <- plot_grid(spec.div.lat,
                        spec.perc.lat, 
                        nrow = 1, 
                        ncol = 2)
super.spec

super.total <- plot_grid(total.div.lat,
                         total.perc.lat,
                         nrow = 1,
                         ncol = 2)
super.total

ffgsuper <- plot_grid(surf.perc.lat, 
                      marg.perc.lat,
                      P.perc.lat,
                      hole.lat,
                      nrow = 1, 
                      ncol = 4)
ffgsuper
#saving figures----
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/surf.perc.lat.pdf", surf.perc.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/marg.perc.lat.pdf", marg.perc.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/P.perc.lat.pdf", P.perc.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/spec.div.lat.pdf", spec.div.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/spec.perc.lat.pdf", spec.perc.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/total.div.lat.pdf", total.div.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/total.perc.lat.pdf", total.perc.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/hole.lat.pdf", hole.lat)

#super figures
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/super.spec.pdf", super.spec)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/super.total.pdf", super.total)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~lat/ffgsuper.pdf", ffgsuper)
