# Load the vcfR package
library(vcfR)
library(adegenet)

# set wd
setwd("~/Desktop/GWS")

# Read in the VCF file
myvcf <- read.vcfR("cflag.filt1.homo.07.05.vcf")

# Get the total number of variants in the VCF file
total_variants <- length(myvcf)

# Subset to 50,000 random variants
set.seed(123)
subset_indices <- sample.int(nrow(myvcf), 50000)
myvcf_subset <- myvcf[subset_indices, ]

# check if vcfR object has a 'genotypes' slot
slotNames(myvcf_subset)
class(slot(myvcf_subset, "gt"))

# Convert to genlight object
mygl <- vcfR2genlight(myvcf_subset)

#ID names
indNames(mygl)
