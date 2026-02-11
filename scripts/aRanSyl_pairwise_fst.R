##### This script was designed to compute pairwise Fst values across populations of aRanSyl  ######

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

load("/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/data/aRanSyl_genind_pops.RData")
load("/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/data/aRanSyl_gl_pops.RData")


## Computing pairwise Fst test


#aRanSyl_fst <- genet.dist(aRanSyl_genind_pops, method = "WC84") %>% round(digits = 3)
#save(aRanSyl_fst, file = "/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/data/aRanSyl_fst.RData") ## saving here just in case the commands below fail
#aRanSyl_fst_mat <- as.matrix(aRanSyl_fst)
#aRanSyl_fst_mat[aRanSyl_fst_mat < 0] <- 0
#aRanSyl_fst_mat[1:4, 1:4] ## I ran this script using 1:10 in this command and it gave me an error because there are only 9 columns, not 10
#save(aRanSyl_fst_mat, file = "/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/data/aRanSyl_fst_mat.RData")


aRanSyl_pval_fst <- gl.fst.pop(aRanSyl_gl_pops, nboots = 1000, percent = 95, nclusters = 1, verbose = NULL)
save(aRanSyl_pval_fst, file = "/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/data/aRanSyl_pval_fst.RData") ## saving here just in case the commands below fail
aRanSyl_fst <- as.matrix(as.dist(aRanSyl_pval_fst$Fsts))
aRanSyl_fst[aRanSyl_fst < 0] <- 0 ## I ran this script using aRanSyl_pval_fst here instead, and it gave me an error because the matrix that should have been used is aRanSyl_fst
aRanSyl_fst
save(aRanSyl_fst, file = "/gpfs/gibbs/project/skelly/dp996/Wood-frog_ddRADseq/data/aRanSyl_pval_mat.RData")


