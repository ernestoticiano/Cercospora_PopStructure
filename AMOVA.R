# set wd
setwd("~/Desktop/GWS")

#Packages
library("vcfR")
library("adegenet")
library("tidyr")
library("poppr")
library("ape")

#Opening and examining the vcf file
cercospora.VCF <- read.vcfR("cflag.filt2.homo.renamed.07.05.vcf")
cercospora.VCF

#Read information about population
pop.data <- read.csv("Label.csv")
all(colnames(cercospora.VCF@gt)[-1] == pop.data$AccessID)

#subset
subset.filt1 <- sample(size = 50000, x= c(1:nrow(cercospora.VCF)))
cercospora.VCF.sub <- cercospora.VCF[subset.filt1,]

#genlight object
cercospora.sub.gi <- vcfR2genind(cercospora.VCF.sub)
ploidy(cercospora.sub.gi) <- 1

#Create a genclone object from a genind object
cercospora.sub.gc <- poppr::as.genclone(cercospora.sub.gi)
ploidy(cercospora.sub.gc) <- 1

#add strata
strata(cercospora.sub.gc) <- data.frame(pop.data)
strata(cercospora.sub.gi) <- data.frame(pop.data)

table(strata(cercospora.sub.gc, ~State))
table(strata(cercospora.sub.gc, ~Tree))
table(strata(cercospora.sub.gc, ~State/Tree, combine = FALSE))

#AMOVA
#State~Tree
#no clone correction
cercospora.sub.state.amova <- poppr.amova(cercospora.sub.gc, ~State/Tree, cutoff = 0.4)
cercospora.sub.state.amova
write.table(cercospora.sub.state.amova$componentsofcovariance, sep = ",", file = "~/Desktop/GWS/cercospora.sub.state.amova.csv")

#clone correction
cercospora.sub.state.amovacc <- poppr.amova(cercospora.sub.gc, ~State/Tree, cutoff = 0.4, clonecorrect = TRUE)
cercospora.sub.state.amovacc
write.table(cercospora.sub.state.amovacc$componentsofcovariance, sep = ",", file = "~/Desktop/GWS/cercospora.sub.state.amovacc.csv")

cercospora.sub.state.amova.signif   <- randtest(cercospora.sub.state.amova, nrepet = 999)
cercospora.sub.state.amova.signif 


set.seed(1999)
randtest(cercospora.sub.state.amova, nrepet = 999)


#Tree~State
cercospora.sub.tree.amova <- poppr.amova(cercospora.sub.gc, ~Tree/State, cutoff = 0.4, clonecorrect = TRUE)
cercospora.sub.tree.amova

cercospora.sub.tree.amova.signif   <- randtest(cercospora.sub.tree.amova, nrepet = 999)
cercospora.sub.tree.amova.signif 

set.seed(1999)
randtest(cercospora.sub.tree.amova, nrepet = 999)

##Full data
#Opening and examining the vcf file
cercospora.VCF <- read.vcfR("cflag.filt2.homo.renamed.07.05.vcf")
cercospora.VCF

#Read information about population
pop.data <- read.csv("Label.csv")
all(colnames(cercospora.VCF@gt)[-1] == pop.data$AccessID)

#genlight object
cercospora.gi <- vcfR2genind(cercospora.VCF)
ploidy(cercospora.sub.gi) <- 1

#Create a genclone object from a genind object
cercospora.gc <- poppr::as.genclone(cercospora.gi)
ploidy(cercospora.gc) <- 1

#add strata
strata(cercospora.gc) <- data.frame(pop.data)
strata(cercospora.gi) <- data.frame(pop.data)

table(strata(cercospora.gc, ~State))
table(strata(cercospora.gc, ~Tree))
table(strata(cercospora.gc, ~State/Tree, combine = FALSE))

#AMOVA
#State~Tree
cercospora.state.amova <- poppr.amova(cercospora.gc, ~State/Tree, cutoff = 0.4, clonecorrect = TRUE)
cercospora.state.amova
#write.table(cercospora.state.amova$componentsofcovariance, sep = ",", file = "~/Desktop/GWS.csv")
cercospora.state.amova.signif   <- randtest(cercospora.sub.state.amova, nrepet = 999)
cercospora.state.amova.signif 

set.seed(1999)
randtest(cercospora.state.amova, nrepet = 999)


#Tree~State
cercospora.tree.amova <- poppr.amova(cercospora.sub.gc, ~Tree/State, cutoff = 0.4, clonecorrect = TRUE)
cercospora.tree.amova

cercospora.tree.amova.signif   <- randtest(cercospora.sub.tree.amova, nrepet = 999)
cercospora.tree.amova.signif 

set.seed(1999)
randtest(cercospora.tree.amova, nrepet = 999)







 


