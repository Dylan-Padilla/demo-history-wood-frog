##### This script was designed to estimate population structure of the aRanSyl  ######

## loading packages

library(adegenet)
library(ade4)
library(car)
library(canadamaps)
library(data.table)
library(dartRverse)
library(ecodist)
library(fields)
library(hierfstat)
library(LEA)
library(lfmm)
library(maps)
library(mapplots)
library(mapproj)
library(MASS)
library(rnaturalearth)
library(pegas)
library(poppr)
library(prettymapr)
library(qvalue)
library(RColorBrewer)
library(raster)
library(sf)
library(scales)
library(SeqArray)
library(SeqVarTools)
library(shape)
library(SNPRelate)
library(stringr)
library(viridis)
library(vcfR)
library(xtable)


## loading project

load("/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/outputs/aRanSyl_Ne/aRanSyl_Ne_list.RData")

Nes <- data.frame("Ne_0.05" = rep(NA, length(aRanSyl_Ne_list)),
                           "Ne_0.02" = rep(NA, length(aRanSyl_Ne_list)),
                           "Ne_0.01" = rep(NA, length(aRanSyl_Ne_list)),
                           "Pops" = paste(rep("Pop", 9),seq(1, 9, 1), sep = ""))

Nes

## Getting the Ne values

for(i in 1:length(aRanSyl_Ne_list)){

    obj <- data.frame(aRanSyl_Ne_list[[i]])
    Nes[i, c(1, 2, 3)] <-  as.numeric(obj[6, c(2, 3, 4)])
}
str(Nes)
Nes


## Loading metadata

meta <- read.csv("/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/data/aRanSyl_metadata.csv")
colnames(meta)[17] <- "Pops"
str(meta)

meta$Pops <- as.factor(meta$Pops)
levels(meta$Pops) <- paste(rep("Pop", 9), seq(1, 9, 1), sep = "")

tapply(meta$id, meta$Pops, length)

submeta <- data.frame("longitude" = meta$longitude, "latitude" = meta$latitude,
                      "Pops" = meta$Pops)
data <- merge(submeta, Nes, by = "Pops")
str(data)

files <- list.files("/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/data/wc2.1_30s_bio", full.names = TRUE)
files

sd.temp <- stack(files[14])


png("/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/imgs/Ne_map.png",
width = 7, height = 7, units = 'in', res = 360)


## Plot the world map of annual mean temperature

plot(sd.temp, axes = FALSE, xlim = c(-180, -50), ylim = c(25, 70),
     las = 1, ylab = "", xlab = "", col = "white", legend = FALSE)

grid()
par(new = TRUE)

plot(sd.temp, axes = TRUE, xlim = c(-180, -50), ylim = c(25, 70),
     legend = FALSE, las = 1, xlab = "Longitude", ylab = "Latitude",
     col = brewer.pal(6, "Blues"))

plot(sd.temp, 
     legend.only = TRUE, 
     horizontal = TRUE,
     legend.args = list(text = "Temperature (standard deviation × 100)", side = 3,
                        line = 0.5, cex = 0.6),
     axis.args = list(cex.axis = 0.6, mgp = c(3, 0.2, 0)), ## Adjust axis label size
     smallplot = c(0.13, 0.37, 0.2, 0.21), col = brewer.pal(6, "Blues"))

addnortharrow(pos = "topleft", scale = 0.4)

## Plot the points with colors corresponding to temperature

with(data, points(longitude, latitude, pch = 21,
                  bg = RColorBrewer::brewer.pal(9,"Paired")[Pops],
                  col = alpha("black", 0.3),
                  cex = log10(Ne_0.05)))

legend("bottomright", legend = levels(data$Pops), pch = 16,
       col = RColorBrewer::brewer.pal(9,"Paired"), bty = "n",
       cex = 0.6)

dev.off()
