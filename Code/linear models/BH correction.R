#title: "Gussie's P value correction "
#author: "Gussie and Anshuman"
#date: "4/2/2021"



#BH corrrection


setwd("/Users/gussiemaccracken/Desktop/Review Paper- Ellen Lauren Anshuman")

pvalues <- read.csv("Table 1.csv", header=F)

x <- as.numeric(unlist(pvalues)) 

t1 <- p.adjust(x, method = "BH")
View(t1)

pvalues <- read.csv("Table 2.csv", header=F)

y <- as.numeric(unlist(pvalues)) 

t2 <- p.adjust(y, method = "BH")
View(t2)

