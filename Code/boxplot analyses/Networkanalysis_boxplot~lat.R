#load packages
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(agricolae)
library(EnvStats)

#load data file----
data <- read.csv("LAS analyses/Paleoherb_paper/networklevelAll_v2.csv", header=T, sep=",", 
                 na.strings=  "NA", dec=".", strip.white=TRUE)
#data <- data[c(-24:-27),]

conn <- drop_na(data[,c(4,7)])
ndo <- drop_na(data[,c(4,14)])
cco <- drop_na(data[,c(4,12)])
nest <- drop_na(data[,c(4,13)])
int.even <- drop_na(data[,c(4,23)])
al.int.even <- drop_na(data[,c(4,24)])
h2 <- drop_na(data[,c(4,25)])
c.HL <- drop_na(data[,c(4,38)])
c.LL <- drop_na(data[,c(4,39)])
niche.HL <- drop_na(data[,c(4,34)])
niche.LL <- drop_na(data[,c(4,35)])

#Tukey's test for connectance across lat----
conn.lm <- lm(connectance ~ Latitude, data = conn)
conn.av <- aov(conn.lm)
summary(conn.av)
conn_test <- HSD.test(conn.av, 'Latitude')
conn_test

ros.conn <- rosnerTest(conn$connectance,
                           k = 10)
ros.conn
#remove outliers
conn <- conn[c(-61,-28,-15,-21),]

#rerun tukey test
conn.lm <- lm(connectance ~ Latitude, data = conn)
conn.av <- aov(conn.lm)
summary(conn.av)
conn_test <- HSD.test(conn.av, 'Latitude')
conn_test

#Tukey's test for NDO across lat----
ndo.lm <- lm(NODF ~ Latitude, data = ndo)
ndo.av <- aov(ndo.lm)
summary(ndo.av)
ndo_test <- HSD.test(ndo.av, 'Latitude')
ndo_test

ros.ndo <- rosnerTest(ndo$NODF,
                       k = 10)
ros.ndo

#Tukey's test for cluster coeff. across lat----
cco.lm <- lm(cluster.coefficient ~ Latitude, data = cco)
cco.av <- aov(cco.lm)
summary(cco.av)
cco_test <- HSD.test(cco.av, 'Latitude')
cco_test

ros.cco <- rosnerTest(cco$cluster.coefficient,
                      k = 10)
ros.cco

cco <- cco[c(-15,-28,-21,-61),]

cco.lm <- lm(cluster.coefficient ~ Latitude, data = cco)
cco.av <- aov(cco.lm)
summary(cco.av)
cco_test <- HSD.test(cco.av, 'Latitude')
cco_test

#Tukey's test for interaction evenness across lat----
ie.lm <- lm(interaction.evenness ~ Latitude, data = int.even)
ie.av <- aov(ie.lm)
summary(ie.av)
ie_test <- HSD.test(ie.av, 'Latitude')
ie_test

ros.ie <- rosnerTest(int.even$interaction.evenness,
                      k = 10)
ros.ie
int.even <- int.even[-40,]

ie.lm <- lm(interaction.evenness ~ Latitude, data = int.even)
ie.av <- aov(ie.lm)
summary(ie.av)
ie_test <- HSD.test(ie.av, 'Latitude')
ie_test

#Avg interaction eveneness 
ieavg <- int.even %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(Latitude)%>% #grouping by epoch
  summarise(mean.all= mean(interaction.evenness), #average
            sd.all = sd(interaction.evenness), #standard deviation
            count=n())

#Boxplot interaction evenness----
int.even$Latitude <- factor(int.even$Latitude,
                       levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
int.even <- arrange(int.even, Latitude)
ieavg$Latitude <- factor(ieavg$Latitude,
                          levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
ieavg <- arrange(ieavg, Latitude)

ie.lat <- ggplot(int.even, aes(fill=Latitude, y=interaction.evenness, x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(ieavg), 
            label = c("a","ab","b","ab"), 
            aes(y = c(0.45,0.45,0.45,0.45), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = c("#F5DD90","#F76C5E","#177e89","#a6808c"))+ #color palette created above 
  #ggtitle("")+ #title
  ylab("Interaction Evenness")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
ie.lat

#Tukey's test for nestedness across lat----
nest.lm <- lm(nestedness ~ Latitude, data = nest)
nest.av <- aov(nest.lm)
summary(nest.av)
nest_test <- HSD.test(nest.av, 'Latitude')
nest_test

ros.nest <- rosnerTest(nest$nestedness,
                     k = 10)
ros.nest
nest <- nest[-62,-24]

nest.lm <- lm(nestedness ~ Latitude, data = nest)
nest.av <- aov(nest.lm)
summary(nest.av)
nest_test <- HSD.test(nest.av, 'Latitude')
nest_test

#Avg nestedness  
nestavg <- nest %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(Latitude)%>% #grouping by epoch
  summarise(mean.all= mean(nestedness), #average
            sd.all = sd(nestedness), #standard deviation
            count=n())

#Boxplot nestedness----
nest$Latitude <- factor(nest$Latitude,
                               levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
nest <- arrange(nest, Latitude)
nestavg$Latitude <- factor(nestavg$Latitude,
                         levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
nestavg <- arrange(nestavg, Latitude)

nest.lat <- ggplot(nest, aes(fill=Latitude, y=nestedness, 
                                x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(nestavg), 
            label = c("b","b","b","a"), 
            aes(y = c(0.25,0.25,0.25,0.25), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = c("#F5DD90","#F76C5E","#177e89","#a6808c"))+ #color palette created above 
  #ggtitle("")+ #title
  ylab("Nestedness")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
nest.lat

#Tukey's test for nestedness across lat----
al.lm <- lm(Alatalo.interaction.evenness ~ Latitude, data = al.int.even)
al.av <- aov(al.lm)
summary(al.av)
al_test <- HSD.test(al.av, 'Latitude')
al_test

ros.al <- rosnerTest(al.int.even$Alatalo.interaction.evenness,
                       k = 10)
ros.al

#Tukey's test for H2 lat----
h2.lm <- lm(H2 ~ Latitude, data = h2)
h2.av <- aov(h2.lm)
summary(h2.av)
h2_test <- HSD.test(h2.av, 'Latitude')
h2_test

ros.h2 <- rosnerTest(h2$H2,
                     k = 10)
ros.h2
h2 <- h2[-62,]

h2.lm <- lm(H2 ~ Latitude, data = h2)
h2.av <- aov(h2.lm)
summary(h2.av)
h2_test <- HSD.test(h2.av, 'Latitude')
h2_test

#Avg H2 
h2avg <- h2 %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(Latitude)%>% #grouping by epoch
  summarise(mean.all= mean(H2), #average
            sd.all = sd(H2), #standard deviation
            count=n())

#Boxplot alatalo int.even.----
h2$Latitude <- factor(h2$Latitude,
                               levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
h2 <- arrange(h2, Latitude)
h2avg$Latitude <- factor(h2avg$Latitude,
                         levels = c("High-S","Mid-S","Low","Mid-N","High-N"))
h2avg <- arrange(h2avg, Latitude)

h2.lat <- ggplot(h2, aes(fill=Latitude, y=H2, 
                                  x=Latitude)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(h2avg), 
            label = c("ab","ab","a","a"), 
            aes(y = c(0,0,0,0), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = c("#F5DD90","#F76C5E","#177e89","#a6808c"))+ #color palette created above 
  #ggtitle("")+ #title
  ylab("H2")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
h2.lat

#Tukey's test for C.score.HL lat----
chl.lm <- lm(C.score.HL ~ Latitude, data = c.HL)
chl.av <- aov(chl.lm)
summary(chl.av)
chl_test <- HSD.test(chl.av, 'Latitude')
chl_test

ros.chl <- rosnerTest(c.HL$C.score.HL,
                     k = 10)
ros.chl

#Tukey's test for C.score.LL lat----
cll.lm <- lm(C.score.LL ~ Latitude, data = c.LL)
cll.av <- aov(cll.lm)
summary(cll.av)
cll_test <- HSD.test(cll.av, 'Latitude')
cll_test

ros.cll <- rosnerTest(c.LL$C.score.LL,
                      k = 10)
ros.cll

#Tukey's test for niche overlap HL lat----
nhl.lm <- lm(niche.overlap.HL ~ Latitude, data = niche.HL)
nhl.av <- aov(nhl.lm)
summary(nhl.av)
nhl_test <- HSD.test(nhl.av, 'Latitude')
nhl_test

ros.nhl <- rosnerTest(niche.HL$niche.overlap.HL,
                      k = 10)
ros.nhl

#Tukey's test for niche overlap LL lat----
nll.lm <- lm(niche.overlap.LL ~ Latitude, data = niche.LL)
nll.av <- aov(nll.lm)
summary(nll.av)
nll_test <- HSD.test(nll.av, 'Latitude')
nll_test

ros.nll <- rosnerTest(niche.LL$niche.overlap.LL,
                      k = 10)
ros.nll

#composite figure: 
netbox.lat.final <- plot_grid(ie.lat,
                              h2.lat,
                              nest.lat,
                              nrow = 1,
                              ncol = 3)
netbox.lat.final

#saving figures
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~lat/ie.lat.pdf", ie.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~lat/h2.lat.pdf", h2.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~lat/nest.lat.pdf", nest.lat)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~lat/netbox.lat.pdf", netbox.lat.final)
