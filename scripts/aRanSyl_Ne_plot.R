##### This script was designed to estimate population structure of the aRanSyl  ######

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
library(prettymapr)
library(qvalue)
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

barplot.data <- data.frame("Ne_0.05" = rep(NA, length(aRanSyl_Ne_list)),
                           "Ne_0.02" = rep(NA, length(aRanSyl_Ne_list)),
                           "Ne_0.01" = rep(NA, length(aRanSyl_Ne_list)),
                           "Pops" = paste(rep("Pop", 9),seq(1, 9, 1), sep = ""))

barplot.data

## Getting the Ne values

for(i in 1:length(aRanSyl_Ne_list)){

    obj <- data.frame(aRanSyl_Ne_list[[i]])
    barplot.data[i, c(1, 2, 3)] <-  as.numeric(obj[6, c(2, 3, 4)])
}
str(barplot.data)
rownames(barplot.data) <- barplot.data[ , 4]
barplot.data <- as.matrix(barplot.data[ , -4])
barplot.data


## getting the CIs

ci.data <- data.frame("CI_lower_0.05" = rep(NA, length(aRanSyl_Ne_list)),
                      "CI_upper_0.05" = rep(NA, length(aRanSyl_Ne_list)),
                      "CI_lower_0.02" = rep(NA, length(aRanSyl_Ne_list)),
                      "CI_upper_0.02" = rep(NA, length(aRanSyl_Ne_list)),
                      "CI_lower_0.01" = rep(NA, length(aRanSyl_Ne_list)),
                      "CI_upper_0.01" = rep(NA, length(aRanSyl_Ne_list)),
                      "Pops" = paste(rep("Pop", 9),seq(1, 9, 1), sep = ""))


ci.data

for(i in 1:length(aRanSyl_Ne_list)){

    obj <- data.frame(aRanSyl_Ne_list[[i]])
    ci.data[i, 1:6] <-  c(as.numeric(obj[7:8, 2]),
                          as.numeric(obj[7:8, 3]),
                          as.numeric(obj[7:8, 4]))
}
str(ci.data)
rownames(ci.data) <- ci.data[ , 7]
#ci.data <- as.matrix(ci.data[ , -7])
ci.data


png("/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/imgs/Pops_Ne.png",
    width = 7, height = 7, units = "in", res = 360)


bp <- barplot(t(barplot.data), beside = TRUE, col = "white", ylim = c(0, 130),
        xlab = "", ylab = "", las = 1, xaxt = "n", yaxt = "n", border = FALSE)

grid()
par(new = TRUE)

barplot(t(barplot.data), beside = TRUE,
        col = viridis(3),
        ylim = c(0, 130) ,
        ylab = expression("Effective Population size"~italic((N[e]))),
        border = FALSE, las = 1)
box()

for(i in seq(2, 27, 3)){
    axis(1, at = bp[i], labels = "")
}

ci <- aggregate(cbind(CI_lower_0.05, CI_upper_0.05, CI_lower_0.02, CI_upper_0.02, CI_lower_0.01, CI_upper_0.01)~Pops , data=ci.data , print)
rownames(ci) <- ci[,1]
ci <- ci[, -1]
ci

CIs <- rbind(abs(barplot.data[ , 1]-ci[ , 1]),
             abs(barplot.data[ , 1]-ci[ , 2]),
             abs(barplot.data[ , 2]-ci[ , 3]),
             abs(barplot.data[ , 2]-ci[ , 4]),
             abs(barplot.data[ , 3]-ci[ , 5]),
             abs(barplot.data[ , 3]-ci[ , 6]))

rownames(CIs) <- colnames(ci)
CIs

Nes <- rbind(barplot.data[ , 1],
             barplot.data[ , 2],
             barplot.data[ , 3])
Nes

bp

idx <- 1
idx2 <- 2
for(i in 1:nrow(Nes)){

    arrows(x0 = bp[i, ], y0 = abs(Nes[i, ]-CIs[idx, ]),
           x1 = bp[i, ], y1 = abs(Nes[i, ]+CIs[idx2, ]),
           code = 3, angle = 90, length = 0.04, col = "orange")
    idx = idx + 2
    idx2 = idx2 + 2
}

legend("topright", legend = c("0.05", "0.02", "0.01"), fill = viridis(3),
       title = "Critical value", bty = "n", cex = 0.8, border = FALSE)


dev.off()
