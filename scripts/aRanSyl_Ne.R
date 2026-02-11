##### This script was designed to estimate Ne across populations of aRanSyl  ######

## libraries

library(adegenet)
library(ade4)
library(car)
library(canadamaps)
library(data.table)
library(dartRverse)
library(dartR.popgen)
library(ecodist)
library(geohippos)
library(ggplot2)
library(hierfstat)
library(LEA)
library(lfmm)
library(maps)
library(mapplots)
library(mapproj)
library(rnaturalearth)
library(pegas)
library(poppr)
library(prettymapr)
library(qvalue)
library(sf)
library(scales)
library(SeqArray)
library(SeqVarTools)
library(shape)
library(SNPRelate)
library(stringr)
library(vcfR)
library(xtable)

## set paths for different versions

path.binaries <- "/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/binaries/"
load("/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/data/aRanSyl_gl_list.RData")
aRanSyl_gl_list

aRanSyl_Ne_list <- vector("list", 9)

for(i in 1:9){

    ## current effective population size

    aRanSyl_gl <- gl.filter.allna(aRanSyl_gl_list[[i]])
    aRanSyl_gl

    aRanSyl_Ne <- gl.LDNe(aRanSyl_gl,  outfile = "popsLD.txt",
                      neest.path = path.binaries,
                      critical = c(0.01, 0.02, 0.05),
                      singleton.rm = TRUE, mating = "random", plot.out = FALSE)
    aRanSyl_Ne_list[[i]] <- aRanSyl_Ne
    names(aRanSyl_Ne_list)[[i]] <- paste("aRanSyl_Ne_pop","_",i, sep = "")
    save(aRanSyl_Ne_list, file = "/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/outputs/aRanSyl_Ne/aRanSyl_Ne_list.RData")
}




