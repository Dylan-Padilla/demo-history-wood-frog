##### This script was designed to plot Fst values across populations of the aRanSyl  ######

## loading packages

library(adegenet)
library(ade4)
library(car)
library(canadamaps)
library(data.table)
library(dartRverse)
library(ecodist)
library(gplots)
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


## loading project

load("/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/data/aRanSyl_fst_mat.RData")

fst <- as.matrix(aRanSyl_fst_mat)
fst[fst < 0] <- 0
colnames(fst) <- paste(rep("Pop", 9), seq(1, 9, 1), sep = "")
rownames(fst) <- paste(rep("Pop", 9), seq(1, 9, 1), sep = "")

get_lower_tri<-function(cor_matrix){
  cor_matrix[upper.tri(cor_matrix)] <- NA
  return(cor_matrix)
}

lower_tri_matrix <- get_lower_tri(round(fst, 2))

png("/Users/dpadil10/ASU Dropbox/Dylan Padilla/Yale/aRanSyl-PopGen/imgs/Pops_Fst.png",
    width = 7, height = 7, units = "in", res = 360)



my_palette <- colorRampPalette(c("lightblue", "red"))(n = 100)

heatmap.2(lower_tri_matrix,
          cellnote = lower_tri_matrix,
          trace = "none",
          na.color = "white",
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          col = my_palette,
          key = TRUE,
          offsetRow = -30,
#          symbreaks = TRUE,
          key.title = "Fst",
          density.info = "none"
          )

dev.off()



# 1. Setup the Data
values <- c(
  0,
  0.45, 0,
  0.53, 0.49, 0,
  0.30, 0.42, 0.70, 0,
  0.32, 0.66, 0.76, 0.48, 0,
  0.55, 0.47, 0.24, 0.67, 0.75, 0,
  0.16, 0.53, 0.64, 0.30, 0.20, 0.66, 0,
  0.06, 0.45, 0.55, 0.27, 0.30, 0.58, 0.09, 0,
  0.66, 0.64, 0.48, 0.81, 0.83, 0.25, 0.78, 0.70, 0
)

# Create matrix
fst_matrix <- matrix(NA, nrow=9, ncol=9)
fst_matrix[lower.tri(fst_matrix, diag=TRUE)] <- values
labs <- paste0("Pop", 1:9)

# 2. Prepare the Plotting Area
# This is the key fix: par(mar=...) sets the margins manually.
# Order: Bottom, Left, Top, Right. 
# We give extra space (5 lines) to the Left for your labels.
par(mar = c(5, 5, 2, 2))

# 3. Create the Image
# We transform the matrix so it looks like the table (Row 1 at top)
# t() transposes it, and we reverse the columns to flip the y-axis
grid_data <- t(fst_matrix[nrow(fst_matrix):1, ])

# Generate the heatmap colors
my_colors <- colorRampPalette(c("#8ccddb", "#bf6c6c", "red"))(100)

image(
  1:9, 1:9,           # X and Y coordinates
  grid_data,          # The data
  axes = FALSE,       # Turn off default axes (we will draw our own)
  col = my_colors,    # Colors
  xlab = "", ylab = ""
)

# 4. Add Custom Axes
# Axis 2 is the Left side. 
# We define positions 1:9, but label them in reverse (Pop9 at bottom, Pop1 at top)
axis(2, at = 1:9, labels = rev(labs), las = 2, tick = FALSE, cex.axis = 1)

# Axis 1 is the Bottom side.
axis(1, at = 1:9, labels = labs, las = 2, tick = FALSE, cex.axis = 1)

# 5. Add Text Values
# We loop through the original matrix to place the numbers
n <- 9
for (row in 1:n) {
  for (col in 1:n) {
    val <- fst_matrix[row, col]
    
    if (!is.na(val)) {
      # Logic: 
      # x coordinate = column index
      # y coordinate = n + 1 - row index (because y=1 is the bottom)
      text(x = col, y = n + 1 - row, labels = val, col = "cyan", cex = 0.9)
    }
  }
}

# --- 3. ADD LEGEND (TOP RIGHT) ---
# We draw the legend in the empty white space (Coordinate x=6 to x=8 approx)

# Legend settings
leg_x <- 6.5       # X starting position
leg_y <- 6.0       # Y starting position (bottom of legend)
leg_w <- 0.5       # Width of the bar
leg_h <- 2.5       # Height of the bar
num_cols <- length(my_colors)

# Draw the Gradient Bar
# We stack tiny rectangles on top of each other
for (i in 1:num_cols) {
  y_start <- leg_y + (i - 1) * (leg_h / num_cols)
  y_end   <- leg_y + (i) * (leg_h / num_cols)
  rect(
    xleft = leg_x, 
    ybottom = y_start, 
    xright = leg_x + leg_w, 
    ytop = y_end, 
    col = my_colors[i], 
    border = NA
  )
}

# Add a black border around the legend bar
rect(leg_x, leg_y, leg_x + leg_w, leg_y + leg_h, border = "black")

# Add Legend Text Labels (0, 0.5, 1.0)
# We map the data range (0 to 0.83) to the legend height
min_val <- 0
max_val <- max(values, na.rm=TRUE) # approx 0.83

# Label position 1 (Bottom)
text(x = leg_x + leg_w + 0.2, y = leg_y, labels = min_val, adj = 0, cex = 0.8)

# Label position 2 (Middle)
text(x = leg_x + leg_w + 0.2, y = leg_y + (leg_h / 2),
     labels = round((max_val/2), 2), adj = 0, cex = 0.8)

# Label position 3 (Top)
text(x = leg_x + leg_w + 0.2, y = leg_y + leg_h, labels = max_val, adj = 0, cex = 0.8)

# Add Title above legend
text(x = leg_x + (leg_w/2), y = leg_y + leg_h + 0.3, labels = "Fst", font = 2, cex = 1)

# --- 4. ADD FIGURE LABEL "B" (Top Left) ---
# side=3 (Top), line=1 (distance from plot), adj implies alignment (-0.1 moves it into left margin)
mtext("B", side = 3, at = -1, line = 1, adj = -0.15, cex = 1, font = 1)
