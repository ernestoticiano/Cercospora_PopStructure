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
***** Object of Class vcfR *****
96 samples
19 CHROMs
1,020,094 variants
Object size: 944.8 Mb
0 percent missing data
*****        *****         *****

#Read information about population
pop.data <- read.csv("Label.csv")
all(colnames(cercospora.VCF@gt)[-1] == pop.data$AccessID)
[1] TRUE

#subset
subset.filt1 <- sample(size = 50000, x= c(1:nrow(cercospora.VCF)))
cercospora.VCF.sub <- cercospora.VCF[subset.filt1,]
***** Object of Class vcfR *****
96 samples
17 CHROMs
50,000 variants
Object size: 46.7 Mb
0 percent missing data
*****        *****         *****

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
State
     Alabama     Arkansas    Louisiana Mississippii     Missouri   Tennessee         Texas 
          10           18           27           10           10           11           10 

table(strata(cercospora.sub.gc, ~Tree))
Tree
 2  1  3 
40 52  4 

table(strata(cercospora.sub.gc, ~State/Tree, combine = FALSE))

#AMOVA
#State~Tree
#no clone correction
cercospora.sub.state.amova <- poppr.amova(cercospora.sub.gc, ~State, cutoff = 0.4)
cercospora.sub.state.amova
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                Df    Sum Sq  Mean Sq
Between samples  6  298972.5 49828.75
Within samples  89 3394612.4 38141.71
Total           95 3693584.9 38879.84

$componentsofcovariance
                                 Sigma          %
Variations  Between samples   880.8865   2.257375
Variations  Within samples  38141.7119  97.742625
Total variations            39022.5983 100.000000

$statphi
                         Phi
Phi-samples-total 0.02257375

#pvalue
cercospora.sub.state.amova.signif   <- randtest(cercospora.sub.state.amova, nrepet = 999)
cercospora.sub.state.amova.signif
Monte-Carlo test
Call: as.randtest(sim = res, obs = sigma[1])

Observation: 880.8865 

Based on 999 replicates
Simulated p-value: 0.054 
Alternative hypothesis: greater 

      Std.Obs   Expectation      Variance 
     1.823301    -22.132708 245288.580231 

#clone correction
cercospora.sub.state.amovacc <- poppr.amova(cercospora.sub.gc, ~State, cutoff = 0.4, clonecorrect = TRUE)
cercospora.sub.state.amovacc
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                Df    Sum Sq  Mean Sq
Between samples  6  298972.5 49828.75
Within samples  89 3394612.4 38141.71
Total           95 3693584.9 38879.84

$componentsofcovariance
                                 Sigma          %
Variations  Between samples   880.8865   2.257375
Variations  Within samples  38141.7119  97.742625
Total variations            39022.5983 100.000000

$statphi
                         Phi
Phi-samples-total 0.02257375

#pvalue
cercospora.sub.state.amova.signifcc   <- randtest(cercospora.sub.state.amovacc, nrepet = 999)
cercospora.sub.state.amova.signifcc
Monte-Carlo test
Call: as.randtest(sim = res, obs = sigma[1])

Observation: 880.8865 

Based on 999 replicates
Simulated p-value: 0.057 
Alternative hypothesis: greater 

     Std.Obs  Expectation     Variance 
1.743736e+00 5.147341e+00 2.522250e+05 

#Tree~State
#no clone correction
cercospora.sub.tree.amova <- poppr.amova(cercospora.sub.gc, ~Tree, cutoff = 0.4)
cercospora.sub.tree.amova
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                Df     Sum Sq  Mean Sq
Between samples  2   71542.42 35771.21
Within samples  93 3622042.44 38946.69
Total           95 3693584.86 38879.84

$componentsofcovariance
                                 Sigma          %
Variations  Between samples  -124.5287  -0.320767
Variations  Within samples  38946.6929 100.320767
Total variations            38822.1642 100.000000

$statphi
                          Phi
Phi-samples-total -0.00320767

#pvalue
cercospora.sub.tree.amova.signif   <- randtest(cercospora.sub.tree.amova, nrepet = 999)
cercospora.sub.tree.amova.signif 
Monte-Carlo test
Call: as.randtest(sim = res, obs = sigma[1])

Observation: -124.5287 

Based on 999 replicates
Simulated p-value: 0.537 
Alternative hypothesis: greater 

      Std.Obs   Expectation      Variance 
-2.938267e-01  1.355420e+01  2.208499e+05 

#clone correction
cercospora.sub.tree.amovacc <- poppr.amova(cercospora.sub.gc, ~Tree, cutoff = 0.4, clonecorrect = TRUE)
cercospora.sub.tree.amovacc
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                Df     Sum Sq  Mean Sq
Between samples  2   71542.42 35771.21
Within samples  93 3622042.44 38946.69
Total           95 3693584.86 38879.84

$componentsofcovariance
                                 Sigma          %
Variations  Between samples  -124.5287  -0.320767
Variations  Within samples  38946.6929 100.320767
Total variations            38822.1642 100.000000

$statphi
                          Phi
Phi-samples-total -0.00320767

#pvalue
cercospora.sub.tree.amova.signifcc   <- randtest(cercospora.sub.tree.amova, nrepet = 999)
cercospora.sub.tree.amova.signifcc 
Monte-Carlo test
Call: as.randtest(sim = res, obs = sigma[1])

Observation: -124.5287 

Based on 999 replicates
Simulated p-value: 0.553 
Alternative hypothesis: greater 

      Std.Obs   Expectation      Variance 
-3.015535e-01  1.005121e+01  1.991735e+05 
