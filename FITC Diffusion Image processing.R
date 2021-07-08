### FITC Diffusion Image Processing -- Andrew R Gross -- 2020/03/01
### The script reads in timeseries data in TIFF format, averages intensity across the length of a channel,
### And plots the intensity based on its position and time.
##########################################################################################################################################
### 1.0 - Header
##########################################################################################################################################
require(raster)
require(rgdal)
library(plotly)
library(ggplot2)
library(zoo)

##########################################################################################################################################
### 2.0 - Import Raster layers and stacks.
##########################################################################################################################################
path.prefix = "C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/"

### 2020
setwd("E245 - FITC/E245_FITC_02-24/stacktoprocess/"); time.step = 20 ; title = "FITC-BSA 2-24_B"
setwd("E245 - FITC/E245_FITC_02-21/stacktoprocess/"); time.step = 20 ; title = "FITC-Dex 2-21_A"
setwd("E245 - FITC/E245_FITC_02-17and19/stacktoprocess/"); time.step = 20; title = "FITC-BSA 2-17_A"
setwd("E245 - FITC/E245_FITC_02-04/stacktoprocess/"); time.step = 20; title = "FITC 02-04_A"
setwd("E245 - FITC/E245_6-15 iEC seeding assessment/stacktoprocess/"); time.step = 20; title = "FITC 06-15"
## June/July
setwd("E300 - iEC FITC Assay/06-27_07-15/sequencetoprocess-6-27_7-15_A/"); time.step = 10 ; title = "FITC 7-15 A"     #
setwd("E300 - iEC FITC Assay/06-27_07-15/sequencetoprocess-6-27_07-15_B"); time.step = 10 ; title = "FITC 7-15 B"     #
setwd("E300 - iEC FITC Assay/07-23_07-27/sequence_07-23_07-27"); time.step = 10 ; title = "FITC 7-27"
setwd("E300 - iEC FITC Assay/E300_HU_07-24/sequence-HU-empty"); time.step = 10
setwd("E300 - iEC FITC Assay/E300_HU_07-24/sequence-HU"); time.step = 10
setwd("E300 - iEC FITC Assay/E300_HU_07-24/sequence-HU"); time.step = 10
setwd("E300 - iEC FITC Assay/07-25_07-29/sequence-07-25_07-29 cells/"); time.step = 10 ; title = "FITC 7-25 cells"
setwd("E300 - iEC FITC Assay/07-25_07-29/sequence-07-25_07-29 no cells//"); time.step = 10 ; title = "FITC 7-25 no cells"     #
## Dec.
setwd("E352 - FITC Analysis/s12-21/sequence_s12-21_C_no_cells/"); time.step = 10 ; title = "FITC s12-21 no cells"     #
setwd("E352 - FITC Analysis/s12-21/sequence_s12-21_D_cells/"); time.step = 10 ; title = "FITC s12-21 w cells"     #
#setwd("E352 - FITC Analysis/s12-21/seq_s12-21_C_no_cells_raw/"); time.step = 10 ; title = "FITC s12-21 no cells"
#setwd("E352 - FITC Analysis/s12-21/seq_s12-21_D_cells_raw/"); time.step = 10 ; title = "FITC s12-21 w cells"

### 2021, Jan
setwd("E352 - FITC Analysis/c01-18_s01-18/seq_chA_cells/"); time.step = 20; title = "FITC-A-BSA cells"            #
setwd("E352 - FITC Analysis/c01-18_s01-18/seq_chA_cells_non-norm/"); time.step = 20; title = "FITC-BSA, c01-18_A (cells)"
setwd("E352 - FITC Analysis/c01-18_s01-18/seq_chB_nocells/"); time.step = 20; title = "FITC-BSA, c01-18_B (nocells)"        
## Feb
#setwd("E352 - FITC Analysis/c02-18/seq 0.5 001/"); time.step = 20 ; title = "FITC-dex_c02-17_nc"
#setwd("E352 - FITC Analysis/c02-18/seq 0.5 002/"); time.step = 20 ; title = "FITC-dex_c02-17_nc 5-10min"     
setwd("E352 - FITC Analysis/c02-18/seq B 0.5 00001/"); time.step = 20 ; title = "FITC-dex_c02-17_B_nc"     #
setwd("E352 - FITC Analysis/c02-18/seq 0.5 001_002/"); time.step = 20 ; title = "FITC-dex_c02-17_nc_joined"     #
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/c02-19/E352-c02-19A_timeseries_reg/"); time.step = 20 ; title = "FITC-dex_c19FEB_A_nc"    
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/c02-19/E352-c02-19B_timeseries_reg/"); time.step = 20 ; title = "FITC-dex_c19FEB_B_nc"    
## Mar
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/03-March/c04MAR/seq - B/"); time.step = 20 ; title = 'FITC-dex_c04MAR_B - dead cells'
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/03-March/c17MAR/sequence-c17MAR_B/"); time.step = 20 ; title = 'FITC-dex_c17MAR_B - cells'
## Apr
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/04-April/c02APR/seq_A/"); time.step = 20 ; title = 'FITC-dex_c02APR_A - dead cells'
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/04-April/c02APR/seq_B/"); time.step = 20 ; title = 'FITC-dex_c02APR_B - dead cells'
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/04-April/c05APR/Seq_A"); time.step = 20 ; title = 'FITC-dex_c05APR_A - cells'
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/04-April/c09APR/seq_A_01/"); time.step = 20 ; title = 'FITC-dex_c09APR_A - cells'
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/04-April/c09APR/seq_B"); time.step = 20 ; title = 'FITC-dex_c09APR_B - nc'
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/04-April/c15APR/seq_A_nc/"); time.step = 20 ; title = 'FITC-dex_c15APR_A - nc'
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/04-April/c15APR/seq_B_cells/"); time.step = 20 ; title = 'FITC-dex_c15APR_B - cells'

#####################################################################
## 2.1 - Import selected image set
stack.to.process <- stack(list.files())

#####################################################################
## 2.2 - Report dataset infomation
print(names(stack.to.process)[1])                                    # Print the name of the first image
(assayID = title)                                                    # Print the name assigned to the dataset
print(paste('Dimensions: height =', dim(stack.to.process)[1], ', width =', dim(stack.to.process)[2]))  # Print the dimensions
print(paste('Time stop =', time.step, 'seconds'))                                                      # Print the time step
print(paste(dim(stack.to.process)[3], 'time points. Total run = ', dim(stack.to.process)[3]*time.step/60, 'minutes'))  # Print the length of the time series

##########################################################################################################################################
### 3.0 - Automated generation of intensity profiles for each layer
##########################################################################################################################################
### The following generates a matrix of the averaged intensity of each frame
row.averages.matrix <- matrix(0, nrow = dim(stack.to.process)[3], ncol = dim(stack.to.process)[1] - 8) # Generate an empty matrix with a row for each time point as long as the image height
nlayers = dim(stack.to.process)[3]                                  # Define the number of images to loop through

for(stacklayer in 1:nlayers){                           #$$$$$$$$$$$ ~0.5 sec/slice $$$$$$$$$$$
  layer.current = stack.to.process[[stacklayer]]        ### Select the current layer data
  layer.matrix <- as.matrix(layer.current)+1            ### Convert to Matrix
  layer.matrix <- log(layer.matrix)+1                   ### Log transform the intensity
  row.averages <- apply(layer.matrix,1,mean)            ### Calcuate row averages
  row.averages <- rollmean(row.averages, k = 9)         ### Smooth the data
  row.averages.matrix[stacklayer,] <- row.averages      ### Assign to the matrix
}

##########################################################################################################################################
### 4.0 - Plot five intensity profiles summarizing the time series
##########################################################################################################################################
### Define the timesteps of the timeseries and the plot
timestepInSeconds = time.step                                        # Define the time step in seconds
umPerPixel = 2.5                                                     # Define the distance each pixel represents in um
lengthinminutes = timestepInSeconds * dim(row.averages.matrix)[1]/60 # Redefine the length in minutes
(plottingTimestepInMinutes = floor(lengthinminutes/5))               # Select a timestep to represent the diffusion using five evenly spaced time points 
(frameJump = plottingTimestepInMinutes * 60 / timestepInSeconds)

#####################################################################
## 4.1 - Subsample matrix based on the selected time points and populate a dataframe
position = umPerPixel * 1:dim(row.averages.matrix)[2]
diffusion.df <- data.frame(position, t0=0, t1=0, t2=0, t3=0, t4=0, t5=0)  # Generate a dataframe with the position as rows

for (timepoint in 1:6) {                                             
  row = 1 + (timepoint - 1) * frameJump
  newColumn = row.averages.matrix[row,]
  diffusion.df[timepoint+1] = newColumn
}

#####################################################################
## 4.2 - Plot multiple intensity profiles using the data frame of subsampled time points
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

#####################################################################
## 4.3 - Display Plot of multiple intensity profiles
diff.front +coord_flip()

#####################################################################
## 4.4 - Identify channel edges and fluorescence cutoffs
max.background = max(diffusion.df$t0[1:100])                        # Find the maximum intensity in the first 100 pixels of the first image
cutoff = max.background * 2                                         # Define the cutoff as twice the maximum background intensity
#max.intensity = max(diffusion.df$t0)
above.threshold = which(diffusion.df$t0 > cutoff)                   # Identify which positions are above the cutoff
left.channel.edge <- diffusion.df$position[above.threshold[1]]      # Find the left-most and right-most positions above the cutoff
right.channel.edge <- diffusion.df$position[above.threshold[length(above.threshold)]]

## Plot the intensity profiles with the cutoff and boundary lines superimposed
(boundary.lines <- diff.front + geom_hline(yintercept = cutoff, linetype="dashed", color = "white", size=1) +
  geom_vline(xintercept = left.channel.edge, color = "blue") +
  geom_vline(xintercept = right.channel.edge, color = "blue")) 

boundary.lines + coord_flip()


##########################################################################################################################################
### 5.0 - Export plot
##########################################################################################################################################
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")

png(filename= paste0('Diff ranges over time/', assayID,'_Diff range over time.png'), width = 800, height = 800, units = "px", pointsize = 12)
plot(boundary.lines)
dev.off()

##########################################################################################################################################
### 6.0 - Automated plotting of X-positions above cutoff
##########################################################################################################################################
### Generate a dataframe listing the distance from the channel edge to the diffusion front for both sides of the channel at each timepoint

edge.pos.df <- data.frame(Layer = c(), Time = c(),Left.edge = c(), Right.edge = c(), Width = c(), dist.left = c(), dist.right = c(), Distance = c(), Max_Intensity = c())
                                                                    # Generate an empty data frame
for (layer.num in 1:nlayers(stack.to.process)){                     # Begin loop through the number of slices in the stack
  row.averages <- row.averages.matrix[layer.num,]                   # Define the current row average vector
  
  dye.location <- diffusion.df$position[which(row.averages > cutoff)] # Calculate the locations where intensity exceeds the cutoff value
  
  ### Fill data frame
  time = time.step*(layer.num-1)                                    # Assign time value
  left.edge = dye.location[1]                                       # Assign left position
  right.edge = dye.location[length(dye.location)]                   # Assign right position
  if (left.edge == 0) {right.edge = 0}                              # Error handling
  width = (right.edge - left.edge) * umPerPixel                     # Calculate width from one diffusion front to the other
  dist.left = left.channel.edge - left.edge                         # Calculate the distance from the channel to the left front
  dist.right = right.edge - right.channel.edge                      # Calculate the distance to the right front
  distance = mean(c(dist.left, dist.right))                         # Average the distances of the two diffusion fronts
  maxIntensity = max(row.averages)/100                              # Calculate max intensity
  new.row <- data.frame(layer.num, time, left.edge, right.edge, width, dist.left, dist.right, distance, maxIntensity)
  edge.pos.df <- rbind(edge.pos.df, new.row)                        # Assign new values to the data frame
}

edge.pos.df$distance  <- rollmean(edge.pos.df$distance, k = 3, fill = c(edge.pos.df$distance[1], 0, edge.pos.df$distance[nrow(edge.pos.df)]))
                                                                    # Apply a smoothing function the average of diffusion front distances

#####################################################################
## 6.1 - Examining the data frame of diffusion front positions over time
str(edge.pos.df)
head(edge.pos.df)

#####################################################################
## 6.2 -  Plot the diffusion front distance over the time series
dist.over.time <- ggplot(data = edge.pos.df, aes(x = time/60, y = distance)) +
  geom_line(aes(y = distance), color = 'blue', size = 2) +
  geom_point(aes(y = dist.left), color = 'grey50', pch = '-', size = 7) +
  geom_point(aes(y = dist.right), color = 'grey50', pch = '-', size = 7) +
  scale_x_continuous(name="Time (Min)", limits = c(0,30), breaks = seq(0,10000, 5)) +
  scale_y_continuous(name="Distance from channel (um)", limits = c(0,1800), breaks = seq(0,10000, 100)) +
  labs(title="Diffusion front Distance from Channel Edge Over Time", 
       y = "Distance from Channel (um)", x="Time (min)", subtitle = assayID) +
  theme(panel.background = element_rect(fill = "grey95", colour = "black"), 
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin(2, 2, 2, 2)) ,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14, margin(2, 2, 2, 2)))

#####################################################################
## 6.3 - Display the plot
dist.over.time

#####################################################################
## 6.4 Export the plot
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")
png(filename= paste0('Diff DISTANCE over time/', assayID,'_front dist over time.png'), width = 800, height = 800, units = "px", pointsize = 12)
plot(dist.over.time)
dev.off()

#####################################################################
## 6.5 - Export the diffusion front data frame
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")
write.csv(edge.pos.df, paste0('Edge_position_tables/',title,'_',toupper(format(Sys.Date(),"%d%b%y")),'.csv'),row.names = FALSE)



#
##
###
#####
########
############
#############
############
########
#####
###
##
#




### 7.0 - Scratchwork
#####################################################################
# Comparison of multiple x-positions over cutoff
  
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


### 3.0 - Manual plotting of intensity profile
#####################################################################
### Select stack layer
stacklayer = 10
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

