##### This script was designed to perform quality control of ddRADseq samples  ######

## loading packages

library(adegenet)
library(ade4)
library(car)
library(canadamaps)
library(data.table)
library(dartRverse)
library(ecodist)
library(hierfstat)
library(LEA)
library(lfmm)
library(maps)
library(mapplots)
library(mapproj)
library(rnaturalearth)
library(pegas)
library(poppr)
library(qvalue)
library(sf)
library(SeqArray)
library(SeqVarTools)
library(SNPRelate)
library(stringr)
library(vcfR)
library(xtable)


## reading .vcf file and producing genind object

my_vcf <- read.vcfR("/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/ipyrad-reference/aRanSyl_outfiles/aRanSyl_copy_reheadered_snps.vcf.gz")

aRanSyl_genind <- vcfR2genind(my_vcf)
aRanSyl_genind

## Access coordinates from the vcfR object's fix slot

aRanSyl_genind@other$chromosome <- my_vcf@fix[ , 1]
aRanSyl_genind@other$position <- my_vcf@fix[ , 2]
aRanSyl_genind

## removing samples that did not pass the adapter content check qc from the multiqc_report.html file

ind <- c("s_272_S121_R1_001", "s_346_S186_R1_001")
aRanSyl_genind <- aRanSyl_genind[!indNames(aRanSyl_genind) %in% ind, ]

## calculate the percentage of complete genotypes per loci in the aRanSyl SNP data set

## remove loci with > 20% missing data

aRanSyl_genind <- missingno(aRanSyl_genind, type = "loc", cutoff = 0.20)

## calculate the percentage of complete genotypes per individual in the aRanSyl SNP data set

## remove individuals with > 20% missing genotypes

aRanSyl_genind <- missingno(aRanSyl_genind, type = "geno", cutoff = 0.20)

## filtering out loci with > 2 allelles

n <- names(which(nAll(aRanSyl_genind) == 2))
aRanSyl_genind <- aRanSyl_genind[loc = n]
aRanSyl_genind

## print the number of multilocus genotypes

mlg(aRanSyl_genind)
isPoly(aRanSyl_genind) %>% summary
poly_loci_aRanSyl <- names(which(isPoly(aRanSyl_genind) == TRUE))
aRanSyl_genind <- aRanSyl_genind[loc = poly_loci_aRanSyl]
isPoly(aRanSyl_genind) %>% summary

## convert genind object to a data frame

save(aRanSyl_genind, file = "/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/data/aRanSyl_genind_chr_pos.RData")

## identify duplicated genotypes

dups_aRanSyl <- mlg.id(aRanSyl_genind)
for (i in dups_aRanSyl){ # for each element in the list object
  if (length(dups_aRanSyl[i]) > 1){ # if the length is greater than 1
    print(i) # print individuals that are duplicates
  }
}


## preparing .geno file for pop structure analysis

aRanSyl_gl <- gi2gl(aRanSyl_genind)
save(aRanSyl_gl, file = "/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/data/aRanSyl_gl_chr_pos.RData")