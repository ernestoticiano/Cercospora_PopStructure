library("vcfR")
library("poppr")
library("ape")
library("RColorBrewer")
library("ggplot2")

# set wd
setwd("~/Desktop/GWS")

#Opening and examining the vcf file
cercospora.VCF <- read.vcfR("cflag.filt2.homo.renamed.07.05.vcf")
cercospora.VCF

pop.data <- read.csv("Label.csv", header = TRUE)
all(colnames(cercospora.VCF@gt)[-1] == pop.data$AccessID)

#Converting the dataset to a genlight object
gl.cercospora <- vcfR2genlight(cercospora.VCF)
ploidy(gl.cercospora) <- 2
pop(gl.cercospora) <- pop.data$State

#Population genetic analyses for SNP data
#Distance matrices
x <- vcfR2genlight(cercospora.VCF)
x.dist <- dist(x)
x.dist <- poppr::bitwise.dist(x)

#Distance tree
tree <- aboot(gl.cercospora, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

cols <- brewer.pal(n = nPop(gl.cercospora), name = "Dark2")
plot.phylo(tree, cex = 0.15, font =10, adj =0.5, tip.color =  cols[pop(gl.cercospora)])
nodelabels(tree$node.label, adj = c(1.0, -0.5), frame = "n", cex = 0.3,font = 2, xpd = TRUE)

legend(35,10,c("Mississippii", "Louisiana", "Arkansas", "Texas", "Missouri", "Tennessee", "Alabama"),cols, border = FALSE, bty = "n")
legend('topleft', legend = c("Mississippii", "Louisiana", "Arkansas", "Texas", "Missouri", "Tennessee", "Alabama"), fill = cols, border = FALSE, bty = "n", cex = 1)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")

#Minimum spanning networks
library(igraph)

cercospora.dist <- bitwise.dist(gl.cercospora)
cercospora.msn <- poppr.msn(gl.cercospora, cercospora.dist, showplot = FALSE, include.ties = T)

node.size <- rep(2, times = nInd(gl.cercospora))
names(node.size) <- indNames(gl.cercospora)
vertex.attributes(cercospora.msn$graph)$size <- node.size

set.seed(9)
plot_poppr_msn(gl.cercospora, cercospora.msn , palette = brewer.pal(n = nPop(gl.cercospora), name = "Dark2"), gadj = 70)


#Principal components analysis
#PCA Eingenvalues
cercospora.pca <- glPca(gl.cercospora, nf = 3)
barplot(100*cercospora.pca$eig/sum(cercospora.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

##PCA
cercospora.pca.scores <- as.data.frame(cercospora.pca$scores)
cercospora.pca.scores$pop <- pop(gl.cercospora)

library(ggplot2)
set.seed(9)
p <- ggplot(cercospora.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p


#DAPCA
pnw.dapc <- dapc(gl.cercospora, n.pca = 3, n.da = 2)

scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)

#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc,col = cols, posi = 'top')

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.rubi)
dapc.results$indNames <- rownames(dapc.results)

# library(reshape2)
# dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

head(dapc.results, n = 6)










cols <- brewer.pal(n = nPop(gl.subset.cercospora), name = "Dark2")
##Subset to 50,000 random variants
setwd("~/Desktop/GWS")

#Opening and examining the vcf file
cercospora.VCF <- read.vcfR("cflag.filt2.homo.renamed.07.05.vcf")
cercospora.VCF


set.seed(123)
subset_indices <- sample.int(nrow(cercospora.VCF), 50000)
cercospora.VCF_subset <- cercospora.VCF[subset_indices, ]

#Converting the dataset to a genlight object
gl.subset.cercospora <- vcfR2genlight(cercospora.VCF_subset)
ploidy(gl.subset.cercospora) <- 2
pop(gl.subset.cercospora) <- pop.data$State

#Population genetic analyses for SNP data
#Distance matrices
x <- vcfR2genlight(cercospora.VCF_subset)
x.dist <- dist(x)
x.dist <- poppr::bitwise.dist(x)

#PCA_States
gl.subset.cercospora.pca <- glPca(gl.subset.cercospora, nf = 3)
gl.subset.cercospora.pca.scores <- as.data.frame(gl.subset.cercospora.pca$scores)
gl.subset.cercospora.pca.scores$pop <- pop(gl.subset.cercospora)

library(ggplot2)
set.seed(9)
p <- ggplot(gl.subset.cercospora.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, linewidth = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p

#DAPCA
pnw.dapc <- dapc(gl.subset.cercospora, n.pca = 3, n.da = 2)

scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)

#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc,col = cols, posi = 'top')

dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.subset.cercospora)
dapc.results$indNames <- rownames(dapc.results)

# library(reshape2)
# dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))

head(dapc.results, n = 6)


p <- ggplot(dapc.results, aes(x=indNames, y=value, fill=name))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
