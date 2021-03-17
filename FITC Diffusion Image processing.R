### FITC Diffusion Image Processing -- Andrew R Gross -- 2020/03/01
### The script reads in timeseries data in TIFF format, averages intensity across the length of a channel,
### And plots the intensity based on it's position and time.

### 1.0 - Header
#####################################################################
require(raster)             # Load raster package
require(rgdal)
library(plotly)
library(ggplot2)
library(zoo)

### 2.0 - Import Raster layers and stacks.
#####################################################################
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/")
assay.list <- c()

setwd("E245 - FITC/E245_FITC_02-24/stacktoprocess/"); time.step = 20 ; title = "FITC 2-24"
setwd("E245 - FITC/E245_FITC_02-21/stacktoprocess/"); time.step = 20 ; title = "FITC 2-21"
setwd("E245 - FITC/E245_FITC_02-17and19/stacktoprocess/"); time.step = 20
setwd("E245 - FITC/E245_FITC_02-04/stacktoprocess/"); time.step = 20
setwd("E245 - FITC/E245_6-15 iEC seeding assessment/stacktoprocess/"); time.step = 20

setwd("E300 - iEC FITC Assay/E300_06-27_07-15/sequencetoprocess-6-27_7-15_A"); time.step = 10 ; title = "FITC 7-15 A"
setwd("E300 - iEC FITC Assay/E300_06-27_07-15/sequencetoprocess-6-27_07-15_B"); time.step = 10 ; title = "FITC 7-15 B"
setwd("E300 - iEC FITC Assay/E300_07-23_07-27/sequence_07-23_07-27"); time.step = 10 ; title = "FITC 7-27"
setwd("E300 - iEC FITC Assay/E300_HU_07-24/sequence-HU-empty"); time.step = 10
setwd("E300 - iEC FITC Assay/E300_HU_07-24/sequence-HU"); time.step = 10
setwd("E300 - iEC FITC Assay/E300_HU_07-24/sequence-HU"); time.step = 10

setwd("E300 - iEC FITC Assay/07-25_07-29/sequence-07-25_07-29 cells/"); time.step = 10 ; title = "FITC 7-25 cells"
setwd("E300 - iEC FITC Assay/07-25_07-29/sequence-07-25_07-29 no cells//"); time.step = 10 ; title = "FITC 7-25 no cells"     #

setwd("E352 - FITC Analysis/s12-21/sequence_s12-21_C_no_cells/"); time.step = 10 ; title = "FITC s12-21 no cells"     #
setwd("E352 - FITC Analysis/s12-21/sequence_s12-21_D_cells/"); time.step = 10 ; title = "FITC s12-21 w cells"     #
#setwd("E352 - FITC Analysis/s12-21/seq_s12-21_C_no_cells_raw/"); time.step = 10 ; title = "FITC s12-21 no cells"
#setwd("E352 - FITC Analysis/s12-21/seq_s12-21_D_cells_raw/"); time.step = 10 ; title = "FITC s12-21 w cells"

#setwd("E352 - FITC Analysis/c02-18/seq 0.5 001/"); time.step = 20 ; title = "FITC-dex_c02-17_nc"
#setwd("E352 - FITC Analysis/c02-18/seq 0.5 002/"); time.step = 20 ; title = "FITC-dex_c02-17_nc 5-10min"     
setwd("E352 - FITC Analysis/c02-18/seq B 0.5 00001/"); time.step = 20 ; title = "FITC-dex_c02-17_B_nc"     #
setwd("E352 - FITC Analysis/c02-18/seq 0.5 001_002/"); time.step = 20 ; title = "FITC-dex_c02-17_nc_joined"     #

setwd("E352 - FITC Analysis/c02-19/E352-c02-19A_timeseries_reg/"); time.step = 20 ; title = "FITC-dex_c09MAR_A_nc"    #
setwd("E352 - FITC Analysis/c02-19/E352-c02-19B_timeseries_reg/"); time.step = 20 ; title = "FITC-dex_c09MAR_B_nc"    #

stack.to.process <- stack(list.files())

### Report data info
print(names(stack.to.process)[1])
(assayID = title)
print(paste('Dimensions: height =', dim(stack.to.process)[1], ', width =', dim(stack.to.process)[2]))
print(paste('Time stop =', time.step, 'seconds'))
print(paste(dim(stack.to.process)[3], 'time points. Total run = ', dim(stack.to.process)[3]*time.step/60, 'minutes'))


### 3.0 - Automated generation of intensity profiles for each layer
#####################################################################
### Generates a matrix of the averaged intensity of each frame
row.averages.matrix <- matrix(0, nrow = dim(stack.to.process)[3], ncol = dim(stack.to.process)[1] - 8)
nlayers = dim(stack.to.process)[3]

for(stacklayer in 1:nlayers){                           #$$$$$$$$$$$ ~0.5 sec/slice $$$$$$$$$$$
  layer.current = stack.to.process[[stacklayer]]        ### Select the current layer data
  layer.matrix <- as.matrix(layer.current)+1            ### Convert to Matrix
  layer.matrix <- log(layer.matrix)+1
  row.averages <- apply(layer.matrix,1,mean)            ### Calcuate row averages
  row.averages <- rollmean(row.averages, k = 9)
  row.averages.matrix[stacklayer,] <- row.averages      ### Assign to matrix
}

dim(row.averages.matrix)


### 4.0 - Plot multiple intensity profiles
#####################################################################
### Define the timesteps of the timeseries and the plot
timestepInSeconds = time.step
umPerPixel = 2.5
(lengthinminutes = timestepInSeconds * dim(row.averages.matrix)[1] / 60)
(plottingTimestepInMinutes = floor(lengthinminutes/5))
(frameJump = plottingTimestepInMinutes * 60 / timestepInSeconds)

### Subsample matrix and populated dataframe
position = umPerPixel * 1:dim(row.averages.matrix)[2]
diffusion.df <- data.frame(position, t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0)         # Generate a dataframe with the position as rows

for (timepoint in 1:6) {
  row = 1 + (timepoint - 1) * frameJump
  newColumn = row.averages.matrix[row,]
  #diffusion.df <- cbind(diffusion.df, newColumn)
  diffusion.df[timepoint+1] = newColumn
}

### Plot multiple intensity profiles
diff.front <- ggplot(data=diffusion.df, aes(x = position), label = c("t0","t1","t2", "t3", "t4", "t5")) + 
  geom_line(aes(y=t0, col = 't0'), size = 2) + 
  geom_line(aes(y=t1, col = 't1'), size = 1) +
  geom_line(aes(y=t2, col = 't2'), size = 1) + 
  geom_line(aes(y=t3, col = 't3'), size = 1) +
  geom_line(aes(y=t4, col = 't4'), size = 1) +
  geom_line(aes(y=t5, col = 't5'), size = 1) +
  xlim(000, 6200) +
  labs(title="Diffusion range over Time", y="Relative Fluorescent Intensity", x="Position (px)",
       subtitle = assayID) +
  theme(panel.background = element_rect(fill = "black", linetype = "blank", colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin(0, 0, 0, 0)) ,
        axis.text=element_text(size=14),
        legend.text=element_text(size=14),
        axis.title=element_text(size=16, margin(1, 1, 1, 1))) +
  guides(guide_legend()) +
  scale_color_manual(name="Time", 
                     labels = c("t0"="time 0", 
                                "t1" = paste(plottingTimestepInMinutes * 1, 'min'), 
                                "t2" = paste(plottingTimestepInMinutes * 2, 'min'), 
                                "t3" = paste(plottingTimestepInMinutes * 3, 'min'), 
                                "t4" = paste(plottingTimestepInMinutes * 4, 'min'),
                                "t5" = paste(plottingTimestepInMinutes * 5, 'min')), 
  #                   values = c("t0"="#08519c","t1"="#3182bd","t2"="#6baed6","t3"="#9ecae1", "t4"="#c6dbef","t5"="#eff3ff"),  # Blues for W/O cells
                     values = c("t0"="#006d2c","t1"="#31a354","t2"="#74c476","t3"="#a1d99b","t4"="#c7e9c0","t5"="#edf8e9"),  # Greens for WITH cells
                     limits = c("t0", "t1", "t2", "t3", "t4", "t5")) +
  labs(title="Diffusion range over Time", y="Log Relative Fluorescent Intensity", x="Position (um)")

### Display Plot
diff.front

### 4.1 - Identify channel edges and fluorescence cutoffs
### 
max.intensity = max(diffusion.df$t0)
max.background = max(diffusion.df$t0[1:100])
above.threshold = which(diffusion.df$t0 > max.background*2)
left.channel.edge <- diffusion.df$position[above.threshold[1]]
right.channel.edge <- diffusion.df$position[above.threshold[length(above.threshold)]]

(boundary.lines <- diff.front + geom_hline(yintercept = max.background * 2, linetype="dashed", color = "red", size=1) +
  geom_vline(xintercept = left.channel.edge, color = "blue") +
  geom_vline(xintercept = right.channel.edge, color = "blue"))

### Export plot
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")

png(filename= paste0('Diff ranges over time/', assayID,'_Diff range over time.png'), width = 800, height = 800, units = "px", pointsize = 12)
plot(boundary.lines)
dev.off()

#####################################################################
#####################################################################
### 6.0 - Automated plotting of X-positions above cutoff
#####################################################################
#####################################################################

### Generate a dataframe of diffusion-front positions over time
#time.step = 10
umPerPixel = 2.5
channelWidth = 600
#cutoff = 1000                                             # For log data
cutoff = max.background * 2                                               # For linear data


### Method 2
edge.pos.df <- data.frame(Layer = c(), Time = c(),Left.edge = c(), Right.edge = c(), Width = c(), dist.left = c(), dist.right = c(), Distance = c(), Max_Intensity = c())
for (layer.num in 1:nlayers(stack.to.process)){           # Begin loop through the number of slices in the stack
  row.averages <- row.averages.matrix[layer.num,]
  ### Calculate positions where intensity exceeds cutoff
  dye.location <- diffusion.df$position[which(row.averages > cutoff)]            # Calculate the locations over the fluorescences
  
  ### Fill data frame
  time = time.step*(layer.num-1)                                 # Assign time value
  left.edge = dye.location[1]                             # Assign left pos
  right.edge = dye.location[length(dye.location)]         # Assign right pos
  if (left.edge == 0) {right.edge = 0}                  # Error handling
  width = (right.edge - left.edge) * umPerPixel           # Calculate width of flourescent front
  dist.left = left.channel.edge - left.edge
  dist.right = right.edge - right.channel.edge
  distance = mean(c(dist.left, dist.right))                      # Calcualte distance of flourescent front from starting point (channel edge)
  maxIntensity = max(row.averages)/100                    # Calculate max intensity
  new.row <- data.frame(layer.num, time, left.edge, right.edge, width, dist.left, dist.right, distance, maxIntensity)
  edge.pos.df <- rbind(edge.pos.df, new.row)
}

edge.pos.df$distance  <- rollmean(edge.pos.df$distance, k = 3, fill = c(edge.pos.df$distance[1], 0, edge.pos.df$distance[nrow(edge.pos.df)]))

### Examining diffusion front data frame
str(edge.pos.df)
head(edge.pos.df)

### Plot diffusion front over time
dist.over.time <- ggplot(data = edge.pos.df, aes(x = time/60, y = distance)) +
  geom_line(aes(y = distance), color = 'blue', size = 2) +
  geom_point(aes(y = dist.left), color = 'grey50', pch = '-', size = 7) +
  geom_point(aes(y = dist.right), color = 'grey50', pch = '-', size = 7) +
  scale_x_continuous(name="Time (Min)", breaks = seq(0,10000, 5)) +
  scale_y_continuous(name="Distance from channel (um)", breaks = seq(0,10000, 100)) +
  labs(title="Diffusion front Distance from Channel Edge Over Time", 
       y = "Distance from Channel (um)", x="Time (min)", subtitle = assayID) +
  theme(panel.background = element_rect(fill = "grey95", colour = "black"), 
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin(2, 2, 2, 2)) ,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14, margin(2, 2, 2, 2)))

### Show plot
dist.over.time

### Save plot
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")
png(filename= paste0('Diff fronts over time/', assayID,'_front dist over time.png'), width = 800, height = 800, units = "px", pointsize = 12)
plot(dist.over.time)
dev.off()

### Save Diffusion Front Data
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")
write.csv(edge.pos.df, paste0('Edge_position_tables/',title,'_',toupper(format(Sys.Date(),"%d%b%y")),'.csv'),row.names = FALSE)








### 6.0 - Comparison of multiple x-positions over cutoff
#####################################################################

#names(edge.pos.df)[6] <- assayID
assay.list <- c(assay.list, assayID)

multi.edge.df <- 1
  
  
### Generate first columns
multi.edge.df <- edge.pos.df[c(2,6)]
head(multi.edge.df)

### Generate additional columns
multi.edge.df <- cbind(multi.edge.df, edge.pos.df[6])
head(multi.edge.df)
names(multi.edge.df)[2] <- 'distance.2'

ggplot(data = multi.edge.df, aes(x = time/60)) +
  geom_point(aes(y = distance),color = "navyblue", size = 4) + 
  scale_x_continuous(name="Time (Min)", breaks = seq(0,10000, 5)) +
  scale_y_continuous(name="Distance from channel (um)", breaks = seq(0,10000, 100))


### Plot diffusion front over time
multi.dist.over.time <- ggplot(data = multi.edge.df, aes(x = time/60)) +
  geom_point(aes(y = distance), color = "navyblue", size = 4) +
  geom_point(aes(y = distance.2), color = 'red', size = 4) +
  scale_x_continuous(name="Time (Min)", breaks = seq(0,10000, 5)) +
  scale_y_continuous(name="Distance from channel (um)", breaks = seq(0,10000, 100)) +
  labs(title="Diffusion front Distance from Channel Edge Over Time", 
       y = "Distance from Channel (um)", 
       x="Time (min)",
       subtitle = assayID) +
  theme(panel.background = element_rect(fill = "grey95", linetype = "blank", colour = "black"), 
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin(2, 2, 2, 2)) ,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14, margin(2, 2, 2, 2)))

### Show plot
multi.dist.over.time
assay.list

### Save plot
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")
png(filename= paste0('multi.dist.over.time.png'), width = 800, height = 800, units = "px", pointsize = 12)
plot(multi.dist.over.time)
dev.off()

### 7.0 - Automated plotting of intensity at a position over time
#####################################################################
### Generate a dataframe of diffusion-front positions over time
#
#
#

### 3.0 - Manual plotting of intensity profile
#####################################################################
### Select stack layer
stacklayer = 1
layer.current <- stack.to.process[[stacklayer]]           # Pull the specified layer
#layer.current <- subset(stack.to.process,stacklayer)     # An alternative means of calling a layer

### Convert to Matrix
layer.matrix <- as.matrix(layer.current)+1                # Convert to matrix

### Subset image
range = dim(layer.matrix)[[1]]
offset = 0
#range = 1800
#offset = 250
(xpos = nrow(layer.matrix)/2 + offset - (range/2))
layer.matrix.ss <- layer.matrix[xpos:(xpos+range),1:500]
dim(layer.matrix.ss)
layer.matrix.ss <- log(layer.matrix.ss)+1

### Intensity distribution check
#hist(layer.matrix, breaks = 50)
#hist(layer.matrix.ss, breaks = 50)
#quantile(as.vector(layer.matrix), seq(0,1,0.1))

### Plot as heatmap
#colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
my_palette <- colorRampPalette(c("black", "green"))(n = 299)
heatmap(layer.matrix.ss, Rowv = NA, Colv = NA, col = my_palette, scale = "none")
#heatmap(layer.matrix, Rowv = NA, Colv = NA, col = my_palette)

### Calcuate row averages
row.averages <- apply(layer.matrix.ss,1,mean)
#row.averages <- apply(layer.matrix,1,mean)

plot.20 <- plot(row.averages)



### TS-2 - Manual plotting of positions above cuttoff
#####################################################################
umPerPixel = 2.5
channelWidth = 600

stacklayer = 20

### Calculate cutoff positions
#(cutoff = max(row.averages)/100)                                     # Define the cutoff as 1/20th the max value
cutoff = 1000
row.averages <- row.averages.matrix[stacklayer,]
cutoff.vector <- row.averages > cutoff                              # Report which rows are above the cutoff value
dye.location <- which(row.averages > cutoff)                       # 

### Create a table of edge positions
edge.pos.df <- data.frame(Layer = c(), Time = c(),Left.edge = c(), Right.edge = c(), Width = c(), Distance = c(), Max_Intensity = c())

(time = 20*stacklayer)
(left.edge = dye.location[1])
(right.edge = dye.location[length(dye.location)])
(width = (right.edge - left.edge) * umPerPixel)
(distance = width/2 - channelWidth)
(maxIntensity = max(row.averages)/100)                    # Calculate max intensity

new.row <- data.frame(stacklayer, time, left.edge, right.edge, width, distance, maxIntensity)
(edge.pos.df <- rbind(edge.pos.df, new.row))


plot(row.averages.matrix[10,]);abline(h = cutoff, col = "red")   
abline(v = 1370); abline(v = 1371+(600/2.5))

#plot(log(row.averages.matrix[2,]));abline(h = cutoff, col = "blue")  

ggplot(data = row.averages.matrix, aes()) +
  geom_point(aes(x = row.averages.matrix[10,])) + 
  scale_x_continuous(name="Layer", breaks = seq(0,10000, 10)) +
  scale_y_continuous(name="Distance from channel", breaks = seq(0,10000, 100))

p + scale_y_continuous(breaks=seq(0,40,5))

### Plotting the curve front
diffusion.profile.list <- list()
layer.num = 10
layer.of.interest <- stack.to.process[[layer.num]]      # Call the current layer in the stack
layer.matrix <- as.matrix(layer.of.interest)+1
#layer.matrix <- log(layer.matrix)
### Calculate row averages
dim(layer.matrix)
row.averages <- apply(layer.matrix,1,mean)
plot(row.averages)

### Add to list
diffusion.profile.list[[2]] <- row.averages
diffusion.df <- data.frame(x = 1:2458, "time0" = diffusion.profile.list[[1]], "time4min" = diffusion.profile.list[[2]], 
                           "time8min" = diffusion.profile.list[[3]], "time12min" = diffusion.profile.list[[4]],
                           "time16min" = diffusion.profile.list[[5]], "time20min" = diffusion.profile.list[[6]] )

### Plotting a 3D contour surface 

### Plot fluorescent intensity at three points over time
# Poins of measure: 340, 300, 250, 200
intensity.at.pos.df <- data.frame(time=c(),Int.340 = c(), Int.300 = c(), Int.250 = c(), Int.200 = c())

### Full loop
for (layer.num in 1:nlayers(stack.to.process)){           # Begin loop through the number of slices in the stack
  layer.of.interest <- stack.to.process[[layer.num]]      # Call the current layer in the stack
  layer.matrix <- as.matrix(layer.of.interest)+1
  #layer.matrix <- log(layer.matrix)
  ### Calculate row averages
  row.averages <- apply(layer.matrix,1,mean)

  new.row = data.frame(time = 20*(layer.num-1), Int.340 = row.averages[340], Int.300 = row.averages[300], 
                       Int.250 = row.averages[250], Int.200 = row.averages[200])
  intensity.at.pos.df <- rbind(intensity.at.pos.df, new.row)

}

ggplot(data = intensity.at.pos.df, aes(x = time)) +
  geom_point(aes(y = Int.340, color = 'pos340')) +
  geom_point(aes(y = Int.300, color = 'pos300')) +
  geom_point(aes(y = Int.250, color = 'pos250')) +
  geom_point(aes(y = Int.200, color = 'pos200')) +
  theme(panel.background = element_rect(fill = "black", linetype = "blank", colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(guide_legend()) +
  scale_color_manual(name="Position", 
                     labels = c("pos340" = "Position 340", 
                                "pos300" = "Position 300", 
                                "pos250" = "Position 250",
                                "pos200" = "Position 200"), 
                     values = c("pos340" = "red", 
                                "pos300" = "orange", 
                                "pos250" = "yellow",
                                "pos200" = "green"), 
                     limits = c("pos340", "pos300", "pos250", "pos200")) +
  labs(title="Fluorescent Intensity over Time at a position", y="Fluorescent Intensity", x="Time (Sec)")


### Method 1 for creating a dataframe of diffusion front positiosn
edge.pos.df <- data.frame(Layer = c(), Time = c(),Left.edge = c(), Right.edge = c(), Width = c(), Distance = c(), Max_Intensity = c())
for (layer.num in 1:nlayers(stack.to.process)){           # Begin loop through the number of slices in the stack
  ### Prepare the row averages vector
  layer.of.interest <- stack.to.process[[layer.num]]      # Call the current layer in the stack
  layer.matrix <- as.matrix(layer.of.interest)+1          # Convert to matrix
  layer.matrix <- log(layer.matrix)                       # Log transform
  row.averages <- apply(layer.matrix,1,mean)              # Calculate row averages
  
  ### Calculate positions where intensity exceeds cutoff
  dye.location <- which(row.averages > cutoff)            # Calculate the locations over the fluorescences
  
  ### Fill data frame
  time = time.step*(layer.num-1)                                 # Assign time value
  left.edge = dye.location[1]                             # Assign left pos
  right.edge = dye.location[length(dye.location)]         # Assign right pos
  if (left.edge == 0) {right.edge = 0}                  # Error handling
  width = (right.edge - left.edge) * umPerPixel           # Calculate width of flourescent front
  distance = width/2 - channelWidth                       # Calcualte distance of flourescent front from starting point (channel edge)
  maxIntensity = max(row.averages)/100                    # Calculate max intensity
  new.row <- data.frame(layer.num, time, left.edge, right.edge, width, distance, maxIntensity)
  edge.pos.df <- rbind(edge.pos.df, new.row)
}