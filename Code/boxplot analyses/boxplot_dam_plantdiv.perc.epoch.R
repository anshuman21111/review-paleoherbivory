#Paleoherbivory reveiew paper
#Purpose: Damage diveristy and percentage by epoch
#plant diversity by epoch
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
dam_data <- read.csv("LAS analyses/Paleoherb_paper/allflora_v2.csv", header=T, sep=",", 
                     na.strings=  "NA", dec=".", strip.white=TRUE)
dam_data$Flora <- NULL
#dam_data <- dam_data[c(-20,-66),]
#dam_data <-dam_data[-c(64:982), ] 
#for some reason lots of blank rows are being added to the bottom of 
#the data frame. removes them
#orgnize data----
total.perc <- dam_data[,c(1,6)]#percent total damage
total.div <- drop_na(dam_data[,c(1,17)]) #total diversity
#following FFG same organization as above
spec.perc <- drop_na(dam_data[,c(1,7)]) 
spec.div <- drop_na(dam_data[,c(1,18)]) 
mine.perc <- drop_na(dam_data[,c(1,8)])
mine.div <- drop_na(dam_data[,c(1,19)])
gall.perc <- drop_na(dam_data[,c(1,9)])
gall.div <- drop_na(dam_data[,c(1,20)])
HF.perc <- drop_na(dam_data[,c(1,10)])
sk.perc <- drop_na(dam_data[,c(1,12)])
surf.perc <- drop_na(dam_data[,c(1,13)])
marg.perc <- drop_na(dam_data[,c(1,11)])
P.perc <- drop_na(dam_data[,c(1,14)])
plant.div <- drop_na(dam_data[,c(1,22)])
st.dt.ocur <- drop_na(dam_data[,c(1,23)])
#----
#boxplots
col7 <- c("#F2CC8F","#81B29A","#5F797B","#3D405B","#8F5D5D","#E07A5F","#c89f9c")

#Total percent damage
#tot.perc.data <-compare_means(perc.damage~epoch, data = total.perc)
#write.csv(tot.perc.data, "export/tot.perc.data.csv")

#Tukey's test for significance across time
#Asking does Xdamage change across time?
#total percent damage tukey test----
total.lm <- lm(perc.damage ~ epoch, data = total.perc)
total.av <- aov(total.lm)
summary(total.av)
total_test <- HSD.test(total.av, 'epoch')
total_test

#outliers
ros.totalperc <- rosnerTest(total.perc$perc.damage,
                        k = 10)
ros.totalperc

#Total percent damage
avgtotal.perc <- total.perc %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(epoch)%>% #grouping by epoch
  summarise(mean.all= mean(perc.damage), #average
            sd.all = sd(perc.damage), #standard deviation
            count=n())

total.perc$epoch <- factor(total.perc$epoch,
                              levels = c("Cretaceous", "Paleocene", "Eocene",
                                         "Oligocene","Miocene","Pliocene",
                                         "Pleistocene"))
total.perc <- arrange(total.perc, epoch)

box.total.perc <- ggplot(total.perc, aes(fill=epoch, y=perc.damage, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgtotal.perc),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("% Total Damage")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.total.perc

#Total diversity box plot
#total damage diversity tukey test----
totaldiv.lm <- lm(stand.dt ~ epoch, data = total.div)
totaldiv.av <- aov(totaldiv.lm)
summary(totaldiv.av)
totaldiv_test <- HSD.test(totaldiv.av, 'epoch')
totaldiv_test

ros.totaldiv <- rosnerTest(total.div$stand.dt,
                        k = 10)
ros.totaldiv
total.div <- total.div[-8,]
totaldiv.lm <- lm(stand.dt ~ epoch, data = total.div)
totaldiv.av <- aov(totaldiv.lm)
summary(totaldiv.av)
totaldiv_test <- HSD.test(totaldiv.av, 'epoch')
totaldiv_test

#Total damage diversity
avgtotal.div <- total.div %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(epoch)%>% #grouping by epoch
  summarise(mean.all= mean(stand.dt), #average
            sd.all = sd(stand.dt), #standard deviation
            count=n())
total.div$epoch <- factor(total.div$epoch,
                           levels = c("Cretaceous", "Paleocene", "Eocene",
                                      "Oligocene","Miocene","Pliocene",
                                      "Pleistocene"))
total.div <- arrange(total.div, epoch)

box.total.div <- ggplot(total.div, aes(fill=epoch, y=stand.dt, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  #stat_compare_means(comparisons = my_comparisons)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgtotal.div),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("Number of DTs (standarized to 300 leaves)")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.total.div

#Gall diversity box plot
#gall diversity tukey test----
galldiv.lm <- lm(stand.gall ~ epoch, data = dam_data)
galldiv.av <- aov(galldiv.lm)
summary(galldiv.av)
galldiv_test <- HSD.test(galldiv.av, trt = 'epoch')
galldiv_test

ros.galldiv <- rosnerTest(gall.div$stand.gall,
                           k = 10)
ros.galldiv
gall.div <- gall.div[c(-22,-8),]
galldiv.lm <- lm(stand.gall ~ epoch, data = dam_data)
galldiv.av <- aov(galldiv.lm)
summary(galldiv.av)
galldiv_test <- HSD.test(galldiv.av, trt = 'epoch')
galldiv_test

#gall damage diversity
avggall.div <-gall.div %>%
  drop_na(stand.gall)%>%
  group_by(epoch)%>%
  summarise(mean.gall = mean(stand.gall),
            sd.gall = sd(stand.gall),
            count=n()) #number of samples/sites per epoch)

# gall.div.data <- compare_means(stand.gall~epoch, data = gall.div)
# write.csv(gall.div.data, "export/gall.div.data.csv")
gall.div$epoch <- factor(gall.div$epoch,
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
gall.div <- arrange(gall.div, epoch)

box.gall.div <- ggplot(gall.div, aes(fill=epoch, y=stand.gall, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avggall.div),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("Number of Galling DTs (standarized to 300 leaves)")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.gall.div

#percent gall damage tukey test----
gall.lm <- lm(perc.gall ~ epoch, data = gall.perc)
gall.av <- aov(gall.lm)
summary(gall.av)
gall_test <- HSD.test(gall.av, trt = 'epoch')
gall_test

ros.gallperc <- rosnerTest(gall.perc$perc.gall,
                          k = 10)
ros.gallperc

gall.perc <- gall.perc[c(-49,-14,-8,-1),]
gall.lm <- lm(perc.gall ~ epoch, data = gall.perc)
gall.av <- aov(gall.lm)
summary(gall.av)
gall_test <- HSD.test(gall.av, trt = 'epoch')
gall_test

#percent gall damage
avggall.perc <-gall.perc %>%
  drop_na(perc.gall)%>%
  group_by(epoch)%>%
  summarise(mean.gall = mean(perc.gall),
            sd.gall = sd(perc.gall),
            count=n()) #number of samples/sites per epoch)


#%Gall damage box plot
#compare_means(perc.gall~epoch, data = gall.perc)
gall.perc$epoch <- factor(gall.perc$epoch,
                         levels = c("Cretaceous", "Paleocene", "Eocene",
                                    "Oligocene","Miocene","Pliocene",
                                    "Pleistocene"))
gall.perc <- arrange(gall.perc, epoch)

box.gall.perc <- ggplot(gall.perc, aes(fill=epoch, y=perc.gall, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avggall.perc),
            label = c("a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("% Gall Damage")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.gall.perc

#percent hole damage tukey test----
hole.lm <- lm(perc.hole ~ epoch, data = HF.perc)
hole.av <- aov(hole.lm)
summary(hole.av)
hole_test <- HSD.test(hole.av, trt = 'epoch')
hole_test

ros.hole <- rosnerTest(HF.perc$perc.hole,
                           k = 10)
ros.hole

#percent hole feeding damage
avgHF.perc <- HF.perc %>%
  drop_na(perc.hole)%>%
  group_by(epoch)%>%
  summarise(mean.HF = mean(perc.hole),
            sd.HF = sd(perc.hole),
            count=n()) #number of samples/sites per epoch)

#%Hole Feeding damage box plot
# perc.HF.data <- compare_means(perc.hole~epoch, data = HF.perc)
# write.csv(perc.HF.data, "export/perc.HF.data.csv")
HF.perc$epoch <- factor(HF.perc$epoch,
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
HF.perc <- arrange(HF.perc, epoch)

box.HF.perc <- ggplot(HF.perc, aes(fill=epoch, y=perc.hole, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgHF.perc),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("% Hole Feeding Damage")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.HF.perc

#perc. margin tukey test----
marg.lm <- lm(perc.marg ~ epoch, data = marg.perc)
marg.av <- aov(marg.lm)
summary(marg.av)
m_test <- HSD.test(marg.av, trt = 'epoch')
m_test

ros.marg <- rosnerTest(marg.perc$perc.marg,
                       k = 10)
ros.marg

#percent margin damage
avgmarg.perc <- marg.perc %>%
  drop_na(perc.marg)%>%
  group_by(epoch)%>%
  summarise(mean.P = mean(perc.marg),
            sd.P = sd(perc.marg),
            count=n()) #number of samples/sites per epoch)

#%Margin Feeding damage box plot
#compare_means(perc.marg~epoch, data = marg.perc)
marg.perc$epoch <- factor(marg.perc$epoch,
                        levels = c("Cretaceous", "Paleocene", "Eocene",
                                   "Oligocene","Miocene","Pliocene",
                                   "Pleistocene"))
marg.perc <- arrange(marg.perc, epoch)

box.marg.perc <- ggplot(marg.perc, aes(fill=epoch, y=perc.marg, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgmarg.perc),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("% Margin Feeding Damage")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.marg.perc

#percent mine damage tukey test----
mine.lm <- lm(perc.mine ~ epoch, data = mine.perc)
mine.av <- aov(mine.lm)
summary(mine.av)
mine_test <- HSD.test(mine.av, trt = 'epoch')
mine_test

ros.mineperc <- rosnerTest(mine.perc$perc.mine,
                       k = 10)
ros.mineperc

mine.perc <- mine.perc[c(-60,-8),]
mine.lm <- lm(perc.mine ~ epoch, data = mine.perc)
mine.av <- aov(mine.lm)
summary(mine.av)
mine_test <- HSD.test(mine.av, trt = 'epoch')
mine_test

#percent mine damage
avgmine.perc <- mine.perc %>%
  drop_na(perc.mine)%>%
  group_by(epoch)%>%
  summarise(mean.mine = mean(perc.mine),
            sd.mine = sd(perc.mine),
            count=n()) #number of samples/sites per epoch)

#%mine damage box plot
# perc.mine.data <-compare_means(perc.mine~epoch, data = mine.perc)
# write.csv(perc.mine.data, "export/perc.mine.data.csv")
mine.perc$epoch <- factor(mine.perc$epoch,
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
mine.perc <- arrange(mine.perc, epoch)

box.mine.perc <- ggplot(mine.perc, aes(fill=epoch, y=perc.mine, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgmine.perc),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("% Mine Damage")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.mine.perc

#mine diversity tukey tes----
minediv.lm <- lm(stand.mine ~ epoch, data = mine.div)
minediv.av <- aov(minediv.lm)
summary(minediv.av)
minediv_test <- HSD.test(minediv.av, trt = 'epoch')
minediv_test

ros.minediv <- rosnerTest(mine.div$stand.mine,
                           k = 10)
ros.minediv
mine.div <- mine.div[-8,]
minediv.lm <- lm(stand.mine ~ epoch, data = mine.div)
minediv.av <- aov(minediv.lm)
summary(minediv.av)
minediv_test <- HSD.test(minediv.av, trt = 'epoch')
minediv_test

#mine damage diversity
avgmine.div <- mine.div %>%
  drop_na(stand.mine)%>%
  group_by(epoch)%>%
  summarise(mean.mine = mean(stand.mine),
            sd.mine = sd(stand.mine),
            count=n()) #number of samples/sites per epoch)

#Mine diversity box plot
#compare_means(stand.mine~epoch, data = mine.div)
mine.div$epoch <- factor(mine.div$epoch,
                         levels = c("Cretaceous", "Paleocene", "Eocene",
                                    "Oligocene","Miocene","Pliocene",
                                    "Pleistocene"))
mine.div <- arrange(mine.div, epoch)

box.mine.div <- ggplot(mine.div, aes(fill=epoch, y=stand.mine, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgmine.div),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("Number of Mining DTs (standarized to 300 leaves)")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.mine.div

#percent specilized damage tukey test----
spec.lm <- lm(perc.spec ~ epoch, data = spec.perc)
spec.av <- aov(spec.lm)
summary(spec.av)
spec_test <- HSD.test(spec.av, trt = 'epoch')
spec_test

ros.specperc <- rosnerTest(spec.perc$perc.spec,
                          k = 10)
ros.specperc
spec.perc <- spec.perc[-8,]

spec.lm <- lm(perc.spec ~ epoch, data = spec.perc)
spec.av <- aov(spec.lm)
summary(spec.av)
spec_test <- HSD.test(spec.av, trt = 'epoch')
spec_test

#percent specialized damage
avgspec.perc <- spec.perc %>%
  drop_na(perc.spec)%>%
  group_by(epoch)%>%
  summarise(mean.spec = mean(perc.spec),
            sd.spec = sd(perc.spec),
            count=n()) #number of samples/sites per epoch)

#%spec damage box plot
spec.perc$epoch <- factor(spec.perc$epoch, 
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
spec.perc <- arrange(spec.perc, epoch)
avgspec.perc$epoch <- factor(avgspec.perc$epoch, 
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
avgspec.perc <- arrange(avgspec.perc, epoch)

box.spec.perc <- ggplot(spec.perc, aes(fill=epoch, y=perc.spec, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  geom_text(data= data.frame(avgspec.perc), 
            label = c("a","ab","ab","ab","ab","ab","a"), 
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("% Specialized Damage")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.spec.perc 

#specialized damage diversity tukey test----
specdiv.lm <- lm(stand.spec ~ epoch, data = spec.div)
specdiv.av <- aov(specdiv.lm)
summary(specdiv.av)
specdiv_test <- HSD.test(specdiv.av, trt = "epoch")
specdiv_test

ros.specdiv <- rosnerTest(spec.div$stand.spec,
                           k = 10)
ros.specdiv
spec.div <- spec.div[-8,]
specdiv.lm <- lm(stand.spec ~ epoch, data = spec.div)
specdiv.av <- aov(specdiv.lm)
summary(specdiv.av)
specdiv_test <- HSD.test(specdiv.av, trt = "epoch")
specdiv_test

#specialized damage diversity
avgspec.div <- spec.div %>%
  drop_na(stand.spec)%>%
  group_by(epoch)%>%
  summarise(mean.spec = mean(stand.spec),
            sd.spec = sd(stand.spec),
            count=n()) #number of samples/sites per epoch)

#Spec diversity box plot
# spec.div.data <- compare_means(stand.spec~epoch, data = spec.div)
# write.csv(spec.div.data, "export/spec.div.data.csv")
spec.div$epoch <- factor(spec.div$epoch,
                         levels = c("Cretaceous", "Paleocene", "Eocene",
                                    "Oligocene","Miocene","Pliocene",
                                    "Pleistocene"))
spec.div <- arrange(spec.div, epoch)

box.spec.div <- ggplot(spec.div, aes(fill=epoch, y=stand.spec, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgspec.div),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("Number of Specialized DTs (standarized to 300 leaves)")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.spec.div

#plant diversity tukey test----
plant.lm <- lm(stand.plant.div ~ epoch, data= plant.div)
plant.av <- aov(plant.lm)
summary(plant.av)
plant_test <- HSD.test(plant.av, trt="epoch")
plant_test

ros.plant <- rosnerTest(plant.div$stand.plant.div,
                           k = 10)
ros.plant
plant.div <- plant.div[c(-52,-3),]
plant.lm <- lm(stand.plant.div ~ epoch, data= plant.div)
plant.av <- aov(plant.lm)
summary(plant.av)
plant_test <- HSD.test(plant.av, trt="epoch")
plant_test

#plant diveristy
avgplant.div <- plant.div %>%
  drop_na()%>%
  group_by(epoch)%>%
  summarise(mean.plant = mean(stand.plant.div),
            sd.plant = sd(stand.plant.div),
            count=n())

#Plant diversity box plot
# plant.div.data <- compare_means(stand.plant.div~epoch, data = plant.div)
# write.csv(plant.div.data, "export/plant.div.data.csv")
plant.div$epoch <- factor(plant.div$epoch,
                         levels = c("Cretaceous", "Paleocene", "Eocene",
                                    "Oligocene","Miocene","Pliocene",
                                    "Pleistocene"))
plant.div <- arrange(plant.div, epoch)

box.plant.div <- ggplot(plant.div, aes(fill=epoch, y=stand.plant.div, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgplant.div),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("Plant Diversity (standarized to 300 leaves)")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.plant.div

#piercing/sucking tukey test----
p.lm <- lm(perc.pierce ~ epoch, data = P.perc)
p.av <- aov(p.lm)
summary(p.av)
p_test <- HSD.test(p.av, trt = 'epoch')
p_test

ros.p <- rosnerTest(P.perc$perc.pierce,
                        k = 10)
ros.p
P.perc <- P.perc[c(-50,-55,-40,-62,-33,-57,-8,-49,-36,-58),]
p.lm <- lm(perc.pierce ~ epoch, data = P.perc)
p.av <- aov(p.lm)
summary(p.av)
p_test <- HSD.test(p.av, trt = 'epoch')
p_test

#percent piercing/sucking damage
avgP.perc <- P.perc %>%
  drop_na(perc.pierce)%>%
  group_by(epoch)%>%
  summarise(mean.P = mean(perc.pierce),
            sd.P = sd(perc.pierce),
            count=n()) #number of samples/sites per epoch)


#%Piercing damage box plot
#compare_means(perc.pierce~epoch, data = P.perc)
P.perc$epoch <- factor(P.perc$epoch,
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
P.perc <- arrange(P.perc, epoch)

box.pierce.perc <- ggplot(P.perc, aes(fill=epoch, y=perc.pierce, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgP.perc),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(0,0,0,0,0,0,0), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("% Piercing Damage")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.pierce.perc

#skeletonization tukey test----
sk.lm <- lm(perc.skel ~ epoch, data = sk.perc)
sk.av <- aov(sk.lm)
summary(sk.av)
sk_test <- HSD.test(sk.av, trt = 'epoch')
sk_test
ros.sk <- rosnerTest(sk.perc$perc.skel,
                    k = 10)
ros.sk
sk.perc <- sk.perc[c(-34,-32,-48,-33,-3),]
sk.lm <- lm(perc.skel ~ epoch, data = sk.perc)
sk.av <- aov(sk.lm)
summary(sk.av)
sk_test <- HSD.test(sk.av, trt = 'epoch')
sk_test

#percent skeletonization damage
avgsk.perc <-sk.perc %>%
  drop_na(perc.skel)%>%
  group_by(epoch)%>%
  summarise(mean.sk = mean(perc.skel),
            sd.sk = sd(perc.skel),
            count=n()) #number of samples/sites per epoch)


#%skeletonization damage box plot
# perc.skel.data <- compare_means(perc.skel~epoch, data = sk.perc)
# write.csv(perc.skel.data, "export/perc.skel.data.csv")
sk.perc$epoch <- factor(sk.perc$epoch,
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
sk.perc <- arrange(sk.perc, epoch)

box.skel.perc <- ggplot(sk.perc, aes(fill=epoch, y=perc.skel, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgsk.perc),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("% Skeletonization Damage")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.skel.perc

#perc. surface tukey test----
surf.lm <- lm(perc.surface ~ epoch, data = surf.perc)
surf.av <- aov(surf.lm)
summary(surf.av)
srf_test <- HSD.test(surf.av, trt = 'epoch')
srf_test
ros.srf <- rosnerTest(surf.perc$perc.surface,
                     k = 10)
ros.srf
surf.perc <- surf.perc[c(-8,-11,-60,-53,-52,-38,-56),]
surf.lm <- lm(perc.surface ~ epoch, data = surf.perc)
surf.av <- aov(surf.lm)
summary(surf.av)
srf_test <- HSD.test(surf.av, trt = 'epoch')
srf_test

#percent surface damage
avgsurf.perc <- surf.perc %>%
  drop_na(perc.surface)%>%
  group_by(epoch)%>%
  summarise(mean.P = mean(perc.surface),
            sd.P = sd(perc.surface),
            count=n()) #number of samples/sites per epoch)

#%surface damage box plot
#compare_means(perc.surface~epoch, data = surf.perc)
surf.perc$epoch <- factor(surf.perc$epoch,
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
surf.perc <- arrange(surf.perc, epoch)

box.surf.perc <- ggplot(surf.perc, aes(fill=epoch, y=perc.surface, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgsurf.perc),
            label = c("a","a","a","a","a","a","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  ggtitle("% Surface Feeding Damage")+ #title
  ylab("")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.surf.perc

#percent of specialized dt occurances----
occur.lm <- lm(Perc.Spec.DTOs ~ epoch, data= st.dt.ocur)
occur.av <- aov(occur.lm)
summary(occur.av)
occur_test <- HSD.test(occur.av, trt="epoch")
occur_test

ros.occur <- rosnerTest(st.dt.ocur$Perc.Spec.DTOs,
                        k = 10)
ros.occur
st.dt.ocur <- st.dt.ocur[-45,]
occur.lm <- lm(Perc.Spec.DTOs ~ epoch, data= st.dt.ocur)
occur.av <- aov(occur.lm)
summary(occur.av)
occur_test <- HSD.test(occur.av, trt="epoch")
occur_test

#occurence data
avgoccur <- st.dt.ocur %>%
  drop_na()%>%
  group_by(epoch)%>%
  summarise(mean.P = mean(Perc.Spec.DTOs),
            sd.P = sd(Perc.Spec.DTOs),
            count=n()) #number of samples/sites per epoch)

st.dt.ocur$epoch <- factor(st.dt.ocur$epoch,
                          levels = c("Cretaceous", "Paleocene", "Eocene",
                                     "Oligocene","Miocene","Pliocene",
                                     "Pleistocene"))
st.dt.ocur<- arrange(st.dt.ocur, epoch)

avgoccur$epoch <- factor(avgoccur$epoch,
                           levels = c("Cretaceous", "Paleocene", "Eocene",
                                      "Oligocene","Miocene","Pliocene",
                                      "Pleistocene"))
avgoccur<- arrange(avgoccur, epoch)

occur.epoch <- ggplot(st.dt.ocur, aes(fill=epoch, y=Perc.Spec.DTOs, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgoccur),
            label = c("ab","b","b","b","b","ab","a"),
            aes(y = c(-1,-1,-1,-1,-1,-1,-1), size = 2))+
  scale_fill_manual(values=col7)+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("Percent Occurence of Specialized DTs")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
occur.epoch

#composite figure
superfig <- plot_grid(box.spec.perc,
                      occur.epoch,
                      nrow = 1,
                      ncol = 2)
superfig

#Creating non-sig supplementary figure----
nonsig.super <- plot_grid(box.total.div, box.total.perc,
                          box.gall.div, box.gall.perc,
                          box.mine.div, box.mine.perc,
                          box.spec.div,NULL,
                          box.HF.perc, box.marg.perc,
                          box.pierce.perc,box.skel.perc,
                          box.surf.perc, NULL,
                          box.plant.div,
                          nrow = 8,
                          ncol = 2)

nonsig.super
ggsave("Figures/nonsig.super.pdf", nonsig.super, height = 24, units = "in")
#Creating figures with only Cretaceous-Eocene boxes----
specperc.KE <- filter(spec.perc, epoch != "Oligocene", epoch!="Miocene",
                      epoch !="Pliocene", epoch!="Pleistocene")
occur.KE <- filter(st.dt.ocur, epoch != "Oligocene", epoch!="Miocene",
                                  epoch !="Pliocene", epoch!="Pleistocene")

#Percentage of Specialized damage K-Eocene
occur.lm <- lm(Perc.Spec.DTOs ~ epoch, data= occur.KE)
occur.av <- aov(occur.lm)
summary(occur.av)
occur_test <- HSD.test(occur.av, trt="epoch")
occur_test

ros.occur <- rosnerTest(occur.KE$Perc.Spec.DTOs,
                        k = 10)
ros.occur

avgoccur <- occur.KE %>%
  drop_na()%>%
  group_by(epoch)%>%
  summarise(mean.P = mean(Perc.Spec.DTOs),
            sd.P = sd(Perc.Spec.DTOs),
            count=n()) #number of samples/sites per epoch)

occur.KE$epoch <- factor(occur.KE$epoch,
                           levels = c("Cretaceous", "Paleocene", "Eocene"))
occur.KE<- arrange(occur.KE, epoch)

avgoccur$epoch <- factor(avgoccur$epoch,
                         levels = c("Cretaceous", "Paleocene", "Eocene"))
avgoccur<- arrange(avgoccur, epoch)
occur.epochKE <- ggplot(occur.KE, aes(fill=epoch, y=Perc.Spec.DTOs, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl
  geom_text(data= data.frame(avgoccur),
            label = c("a","b","b"),
            aes(y = c(-1,-1,-1), size = 2))+
  scale_fill_manual(values=c("#F2CC8F","#81B29A","#5F797B"))+ #color palette created above
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("Percent Occurence of Specialized DTs")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
occur.epochKE

#percent specialized damage K-Eocene
spec.lm <- lm(perc.spec ~ epoch, data = specperc.KE)
spec.av <- aov(spec.lm)
summary(spec.av)
spec_test <- HSD.test(spec.av, trt = 'epoch')
spec_test

ros.specperc <- rosnerTest(specperc.KE$perc.spec,
                           k = 10)
ros.specperc

avgspec.perc <- specperc.KE %>%
  drop_na(perc.spec)%>%
  group_by(epoch)%>%
  summarise(mean.spec = mean(perc.spec),
            sd.spec = sd(perc.spec),
            count=n()) #number of samples/sites per epoch)

#%spec damage box plot
specperc.KE$epoch <- factor(specperc.KE$epoch, 
                          levels = c("Cretaceous", "Paleocene", "Eocene"))
specperc.KE <- arrange(specperc.KE, epoch)
avgspec.perc$epoch <- factor(avgspec.perc$epoch, 
                             levels = c("Cretaceous", "Paleocene", "Eocene"))
avgspec.perc <- arrange(avgspec.perc, epoch)

box.specpercKE <- ggplot(specperc.KE, aes(fill=epoch, y=perc.spec, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  geom_text(data= data.frame(avgspec.perc), 
            label = c("a","a","a"), 
            aes(y = c(-1,-1,-1), size = 2))+
  scale_fill_manual(values=c("#F2CC8F","#81B29A","#5F797B"))+ #color palette created above 
  #ggtitle("Average Percent Damage by Epoch")+ #title
  ylab("% Specialized Damage")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.specpercKE

#Saving figures----
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~epoch/spec.perc.pdf", box.spec.perc)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~epoch/occur.pdf", occur.epoch)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~epoch/super.pdf", superfig)

#only K-Eocene figures
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~epoch/box.specpercKE.pdf", box.specpercKE)
ggsave("LAS analyses/Paleoherb_paper/figures/boxplots~epoch/occur.epochKE.pdf", occur.epochKE)

