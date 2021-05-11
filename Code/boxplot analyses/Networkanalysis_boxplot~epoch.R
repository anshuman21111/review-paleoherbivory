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
data <- data[c(-1,-55),]

conn <- drop_na(data[,c(3,7)])
ndo <- drop_na(data[,c(3,14)])
cco <- drop_na(data[,c(3,12)])
nest <- drop_na(data[,c(3,13)])
int.even <- drop_na(data[,c(3,23)])
al.int.even <- drop_na(data[,c(3,24)])
h2 <- drop_na(data[,c(3,25)])
c.HL <- drop_na(data[,c(3,38)])
c.LL <- drop_na(data[,c(3,39)])
niche.HL <- drop_na(data[,c(3,34)])
niche.LL <- drop_na(data[,c(3,35)])

#Tukey's test for connectance across lat----
conn.lm <- lm(connectance ~ epoch, data = conn)
conn.av <- aov(conn.lm)
summary(conn.av)
conn_test <- HSD.test(conn.av, 'epoch')
conn_test

ros.conn <- rosnerTest(data$connectance,
                           k = 10)
ros.conn
#remove outliers
conn <- conn[c(-59,-27,-14),]


#rerun tukey test
conn.lm <- lm(connectance ~ epoch, data = conn)
conn.av <- aov(conn.lm)
summary(conn.av)
conn_test <- HSD.test(conn.av, 'epoch')
conn_test

#Avg connectance ~ epoch 
connavg <- conn %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(epoch)%>% #grouping by epoch
  summarise(mean.all= mean(connectance), #average
            sd.all = sd(connectance), #standard deviation
            count=n())

#Boxplot alatalo int.even.----
conn$epoch <- factor(conn$epoch, 
                         levels = c("Cretaceous", "Paleocene", "Eocene",
                                    "Oligocene","Miocene","Pliocene",
                                    "Pleistocene"))
conn <- arrange(conn, epoch)
connavg$epoch <- factor(connavg$epoch, 
                       levels = c("Cretaceous", "Paleocene", "Eocene",
                                  "Oligocene","Miocene","Pliocene",
                                  "Pleistocene"))
connavg <- arrange(connavg, epoch)

conn.epoch <- ggplot(conn, aes(fill=epoch, y=connectance, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(connavg), 
            label = c("b","a","b","ab","ab","b","b"), 
            aes(y = c(0,0,0,0,0,0,0), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = col7)+ #color palette created above 
  #ggtitle("")+ #title
  ylab("Connectance")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
conn.epoch

#Tukey's test for NDO across lat----
ndo.lm <- lm(NODF ~ epoch, data = ndo)
ndo.av <- aov(ndo.lm)
summary(ndo.av)
ndo_test <- HSD.test(ndo.av, 'epoch')
ndo_test

ros.ndo <- rosnerTest(ndo$NODF,
                       k = 10)
ros.ndo

#Tukey's test for cluster coeff. across lat----
cco.lm <- lm(cluster.coefficient ~ epoch, data = cco)
cco.av <- aov(cco.lm)
summary(cco.av)
cco_test <- HSD.test(cco.av, 'epoch')
cco_test

ros.cco <- rosnerTest(cco$cluster.coefficient,
                      k = 10)
ros.cco

cco <- cco[c(-14,-27,-20,-59),]

cco.lm <- lm(cluster.coefficient ~ epoch, data = cco)
cco.av <- aov(cco.lm)
summary(cco.av)
cco_test <- HSD.test(cco.av, 'epoch')
cco_test

#Tukey's test for interaction evenness across lat----
ie.lm <- lm(interaction.evenness ~ epoch, data = int.even)
ie.av <- aov(ie.lm)
summary(ie.av)
ie_test <- HSD.test(ie.av, 'epoch')
ie_test

ros.ie <- rosnerTest(int.even$interaction.evenness,
                      k = 10)
ros.ie
int.even <- int.even[-39,]

ie.lm <- lm(interaction.evenness ~ epoch, data = int.even)
ie.av <- aov(ie.lm)
summary(ie.av)
ie_test <- HSD.test(ie.av, 'epoch')
ie_test

#Tukey's test for nestedness across lat----
nest.lm <- lm(nestedness ~ epoch, data = nest)
nest.av <- aov(nest.lm)
summary(nest.av)
nest_test <- HSD.test(nest.av, 'epoch')
nest_test

ros.nest <- rosnerTest(nest$nestedness,
                     k = 10)
ros.nest
nest <- nest[c(-60,-23),]


nest.lm <- lm(nestedness ~ epoch, data = nest)
nest.av <- aov(nest.lm)
summary(nest.av)
nest_test <- HSD.test(nest.av, 'epoch')
nest_test

#Tukey's test for nestedness across lat----
al.lm <- lm(Alatalo.interaction.evenness ~ epoch, data = al.int.even)
al.av <- aov(al.lm)
summary(al.av)
al_test <- HSD.test(al.av, 'epoch')
al_test

ros.al <- rosnerTest(al.int.even$Alatalo.interaction.evenness,
                       k = 10)

ros.al

#Tukey's test for H2 lat----
h2.lm <- lm(H2 ~ epoch, data = h2)
h2.av <- aov(h2.lm)
summary(h2.av)
h2_test <- HSD.test(h2.av, 'epoch')
h2_test

ros.h2 <- rosnerTest(h2$H2,
                     k = 10)
ros.h2
h2 <- h2[-60,]


h2.lm <- lm(H2 ~ epoch, data = h2)
h2.av <- aov(h2.lm)
summary(h2.av)
h2_test <- HSD.test(h2.av, 'epoch')
h2_test

#Tukey's test for C.score.HL lat----
chl.lm <- lm(C.score.HL ~ epoch, data = c.HL)
chl.av <- aov(chl.lm)
summary(chl.av)
chl_test <- HSD.test(chl.av, 'epoch')
chl_test

ros.chl <- rosnerTest(c.HL$C.score.HL,
                     k = 10)
ros.chl

#Tukey's test for C.score.LL epoch----
cll.lm <- lm(C.score.LL ~ epoch, data = c.LL)
cll.av <- aov(cll.lm)
summary(cll.av)
cll_test <- HSD.test(cll.av, 'epoch')
cll_test

ros.cll <- rosnerTest(c.LL$C.score.LL,
                      k = 10)
ros.cll

#Avg C score LL ~ epoch 
cllavg <- c.LL %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(epoch)%>% #grouping by epoch
  summarise(mean.all= mean(C.score.LL), #average
            sd.all = sd(C.score.LL), #standard deviation
            count=n())

#Boxplot C.Score LL.----
c.LL$epoch <- factor(c.LL$epoch, 
                     levels = c("Cretaceous", "Paleocene", "Eocene",
                                "Oligocene","Miocene","Pliocene",
                                "Pleistocene"))
c.LL <- arrange(c.LL, epoch)
cllavg$epoch <- factor(cllavg$epoch, 
                        levels = c("Cretaceous", "Paleocene", "Eocene",
                                   "Oligocene","Miocene","Pliocene",
                                   "Pleistocene"))
cllavg <- arrange(cllavg, epoch)

cll.epoch <- ggplot(c.LL, aes(fill=epoch, y=C.score.LL, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(cllavg), 
            label = c("ab","b","ab","a","ab","ab","ab"), 
            aes(y = c(0,0,0,0,0,0,0), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = col7)+ #color palette created above 
  #ggtitle("")+ #title
  ylab("C Score Plants (LL)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
cll.epoch

#Tukey's test for niche overlap HL epoch----
nhl.lm <- lm(niche.overlap.HL ~ epoch, data = niche.HL)
nhl.av <- aov(nhl.lm)
summary(nhl.av)
nhl_test <- HSD.test(nhl.av, 'epoch')
nhl_test

ros.nhl <- rosnerTest(niche.HL$niche.overlap.HL,
                      k = 10)
ros.nhl

#Tukey's test for niche overlap LL epoch----
nll.lm <- lm(niche.overlap.LL ~ epoch, data = niche.LL)
nll.av <- aov(nll.lm)
summary(nll.av)
nll_test <- HSD.test(nll.av, 'epoch')
nll_test

ros.nll <- rosnerTest(niche.LL$niche.overlap.LL,
                      k = 10)
ros.nll

#Avg niche overlap LL ~ epoch 
nllavg <- niche.LL %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(epoch)%>% #grouping by epoch
  summarise(mean.all= mean(niche.overlap.LL), #average
            sd.all = sd(niche.overlap.LL), #standard deviation
            count=n())

#Boxplot niche overlap LL.----
niche.LL$epoch <- factor(niche.LL$epoch, 
                      levels = c("Cretaceous", "Paleocene", "Eocene",
                                 "Oligocene","Miocene","Pliocene",
                                 "Pleistocene"))
niche.LL <- arrange(niche.LL, epoch)
nllavg$epoch <- factor(nllavg$epoch, 
                         levels = c("Cretaceous", "Paleocene", "Eocene",
                                    "Oligocene","Miocene","Pliocene",
                                    "Pleistocene"))
nllavg <- arrange(nllavg, epoch)

nll.epoch <- ggplot(niche.LL, aes(fill=epoch, y=niche.overlap.LL, 
                         x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(nllavg), 
            label = c("ab","a","b","b","ab","b","b"), 
            aes(y = c(0,0,0,0,0,0,0), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = col7)+ #color palette created above 
  #ggtitle("")+ #title
  ylab("Niche Overlap Plants (LL)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
nll.epoch

#making composite figure:
net.box.super <- plot_grid(conn.epoch, 
                           cll.epoch, nll.epoch,
                           nrow=1,
                           ncol = 3)
net.box.super

##############################################################################
#only Cretaceous-Eocene----
conn.KE <- filter(conn, epoch != "Oligocene", epoch!="Miocene",
                  epoch !="Pliocene", epoch!="Pleistocene")
nichell.KE <- filter(niche.LL, epoch != "Oligocene", epoch!="Miocene",
                     epoch !="Pliocene", epoch!="Pleistocene")

cll.KE <- filter(c.LL, epoch != "Oligocene", epoch!="Miocene",
                 epoch !="Pliocene", epoch!="Pleistocene")

conn.lm <- lm(connectance ~ epoch, data = conn.KE)
conn.av <- aov(conn.lm)
summary(conn.av)
conn_test <- HSD.test(conn.av, 'epoch')
conn_test

ros.conn <- rosnerTest(conn.KE$connectance,
                       k = 10)
ros.conn
#remove outliers
conn.KE <- conn.KE[-13,]

#rerun tukey test
conn.lm <- lm(connectance ~ epoch, data = conn.KE)
conn.av <- aov(conn.lm)
summary(conn.av)
conn_test <- HSD.test(conn.av, 'epoch')
conn_test

#Avg connectance ~ epoch 
connavg <- conn.KE %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(epoch)%>% #grouping by epoch
  summarise(mean.all= mean(connectance), #average
            sd.all = sd(connectance), #standard deviation
            count=n())

#Boxplot connectence K-Eocene.----
conn.KE$epoch <- factor(conn.KE$epoch, 
                     levels = c("Cretaceous", "Paleocene", "Eocene"))
conn.KE <- arrange(conn.KE, epoch)
connavg$epoch <- factor(connavg$epoch, 
                        levels = c("Cretaceous", "Paleocene", "Eocene"))
connavg <- arrange(connavg, epoch)

box.conn.KE <- ggplot(conn.KE, aes(fill=epoch, y=connectance, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(connavg), 
            label = c("b","a","b"), 
            aes(y = c(0,0,0), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = c("#F2CC8F","#81B29A","#5F797B"))+ #color palette created above 
  #ggtitle("")+ #title
  ylab("Connectance")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.conn.KE

#Tukey's test for C.score.LL epoch----
cll.lm <- lm(C.score.LL ~ epoch, data = cll.KE)
cll.av <- aov(cll.lm)
summary(cll.av)
cll_test <- HSD.test(cll.av, 'epoch')
cll_test

ros.cll <- rosnerTest(cll.KE$C.score.LL,
                      k = 10)
ros.cll

#Avg C score LL ~ epoch 
cllavg <- cll.KE %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(epoch)%>% #grouping by epoch
  summarise(mean.all= mean(C.score.LL), #average
            sd.all = sd(C.score.LL), #standard deviation
            count=n())

#Boxplot C.Score LL.----
cll.KE$epoch <- factor(cll.KE$epoch, 
                     levels = c("Cretaceous", "Paleocene", "Eocene"))
cll.KE <- arrange(cll.KE, epoch)
cllavg$epoch <- factor(cllavg$epoch, 
                       levels = c("Cretaceous", "Paleocene", "Eocene"))
cllavg <- arrange(cllavg, epoch)

box.cll.KE <- ggplot(cll.KE, aes(fill=epoch, y=C.score.LL, x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(cllavg), 
            label = c("a","a","a"), 
            aes(y = c(0,0,0), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = c("#F2CC8F","#81B29A","#5F797B"))+ #color palette created above 
  #ggtitle("")+ #title
  ylab("C Score Plants (LL)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.cll.KE

#Tukey's test for niche overlap LL epoch----
nll.lm <- lm(niche.overlap.LL ~ epoch, data = nichell.KE)
nll.av <- aov(nll.lm)
summary(nll.av)
nll_test <- HSD.test(nll.av, 'epoch')
nll_test

ros.nll <- rosnerTest(nichell.KE$niche.overlap.LL,
                      k = 10)
ros.nll

#Avg niche overlap LL ~ epoch 
nllavg <- nichell.KE %>%
  drop_na()%>% #dropping NAs at bottom of .csv
  group_by(epoch)%>% #grouping by epoch
  summarise(mean.all= mean(niche.overlap.LL), #average
            sd.all = sd(niche.overlap.LL), #standard deviation
            count=n())

#Boxplot niche overlap LL.----
nichell.KE$epoch <- factor(nichell.KE$epoch, 
                         levels = c("Cretaceous", "Paleocene", "Eocene"))
nichell.KE <- arrange(nichell.KE, epoch)
nllavg$epoch <- factor(nllavg$epoch, 
                       levels = c("Cretaceous", "Paleocene", "Eocene"))
nllavg <- arrange(nllavg, epoch)

box.nll.KE <- ggplot(nichell.KE, aes(fill=epoch, y=niche.overlap.LL, 
                                  x=epoch)) +
  geom_boxplot()+
  geom_jitter(color="black", size=0.4, alpha=0.9)+
  geom_text(data= data.frame(nllavg), 
            label = c("ab","a","b"), 
            aes(y = c(0,0,0), size = 2))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2),drop=FALSE)+ #offsets the xlabesl 
  scale_fill_manual(values = c("#F2CC8F","#81B29A","#5F797B"))+ #color palette created above 
  #ggtitle("")+ #title
  ylab("Niche Overlap Plants (LL)")+ #axis labels
  xlab("")+
  theme_light()+ #light background makes bar plot easier to read
  theme(legend.position = "none") #removes legend
box.nll.KE

#composite figure
KE.super <- plot_grid(box.conn.KE, 
                      box.nll.KE,
                      nrow = 1,
                      ncol = 2)
KE.super

#saving figures----
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~epoch/nll.epoch.pdf", nll.epoch)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~epoch/connectance.pdf", conn.epoch)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~epoch/cll.epoch.pdf", cll.epoch)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~epoch/net.box.final.pdf", net.box.super)

#saving K-Eocene figures
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~epoch/box.nll.KE.pdf", box.nll.KE)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~epoch/box.cll.KE.pdf", box.cll.KE)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~epoch/box.conn.KE.pdf", box.conn.KE)
ggsave("LAS analyses/Paleoherb_paper/figures/network_boxplot~epoch/super.KE.pdf", KE.super)
