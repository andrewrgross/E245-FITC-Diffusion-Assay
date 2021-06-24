### FITC Diffusion Data Comparison -- Andrew R Gross -- 2021/03/10
### The script reads in tables generated in the FITC Diffusion Image Processing script
### And plots the positions releative to one antoher for comparison purposes

### 1.0 - Header
#####################################################################
#require(raster)             # Load raster package
#require(rgdal)
#library(plotly)
library(ggplot2)
library(zoo)
#library('plyr')


### 2.0 - Import Diffusion Data Tables
#####################################################################

setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/Edge_position_tables/")
diffusion.data.list <- list()
(file.list = list.files())

for(file.name in file.list){
  print(paste('Opening', file.name))
  input.data <- read.csv(file.name)
  diffusion.data.list[[length(diffusion.data.list)+1]] <- input.data
  names(diffusion.data.list)[length(diffusion.data.list)] = file.name
}

str(diffusion.data.list)
names(diffusion.data.list)

#(metadata.df <- read.csv('../metadata.csv'))
metadata.df <- data.frame(file.list = file.list, cells = c('dead', 'dead', 'dead', 'cells', 'cells', 'nc', 'nc', 'cells', 'cells', 'nc', 'nc', NA))

## 2.1 - Report data info
##########################

### 3.0 - Plot distance over time of individual diffusion fronts
#####################################################################
### Plot the left and right side distance over time at the cutoff

print(data.frame(file.list))
print(metadata.df)

table.number <- 6
(assayID <- names(diffusion.data.list)[table.number])
edge.pos.df <- diffusion.data.list[[table.number]]

### Plot diffusion front over time
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

### Show plot
dist.over.time

### 4.0 - Generate a table of fronts to plot
#####################################################################
### Combine columns from multiple tables in a single table

## 4.05 - List available datasets
metadata.df
#metadata.df[metadata.df$Dye == 'Dex',]
#metadata.df[metadata.df$Seeded == 'nc',]

## 4.1 - List the assay files to include

#assays.to.include <- c(21,15,11,13,12,8);colors.to.plot = c('grey', 'grey', 'grey', 'grey', 'green', 'green')
#assays.to.include <- c(8, 9, 11, 12, 13, 6); colors.to.plot = c('green', 'green', 'green', 'grey', 'grey', 'grey')
#assays.to.include <- c(4, 5, 8, 9, 6, 7) ;colors.to.plot = c('green', 'green', 'green', 'green', 'grey', 'grey')
assays.to.include <- c(4, 5, 8, 9, 6, 7, 1, 2, 3) ;colors.to.plot = c('green2', 'green2', 'green2', 'green2', 'grey', 'grey', 'grey', 'grey', 'grey')
assays.to.include <- c(10, 11, 4, 5, 8, 9, 6, 7, 1, 2, 3) ;colors.to.plot = c('grey', 'grey', 'green', 'green', 'green', 'green', 'grey60', 'grey60', 'grey60', 'grey60', 'grey60')
assays.to.include <- c(5, 6, 5, 6, 5, 6, 5, 6, 5, 5, 5) ;colors.to.plot = c('green2', 'grey', 'green', 'grey', 'green', 'green', 'grey60', 'grey60', 'grey60', 'grey60', 'grey60')


print(data.frame(file.list[assays.to.include]))

## 4.2 - Create a list of all times to include
times.to.include = c() 
for(assay in assays.to.include){
  times.to.include <- c(times.to.include,diffusion.data.list[[assay]]$time)
}
times.to.include <- sort(unique(times.to.include))

edge.pos.multi <- data.frame('time' = times.to.include)

## 4.3 - For each assay file, merge the column of interest to the full edge.pos.multi df

for(assay in assays.to.include){
  current.df = current.df <- diffusion.data.list[[assay]][c(2,8)]
  current.df$distance <- rollmean(current.df$distance, k = 5, fill = c(current.df$distance[1], 0, current.df$distance[nrow(current.df)]))
  edge.pos.multi = merge(edge.pos.multi, current.df, by = 'time', all = TRUE, suffixes = c('', assay))
}

head(edge.pos.multi)
#names(edge.pos.multi) <- c('time', 'distance', 'distance2', 'distance3', 'distance4', 'distance5' , 'distance6')
#names(edge.pos.multi) <- c('time', 'distance', 'distance2', 'distance3', 'distance4', 'distance5' , 'distance6', 'distance7')
names(edge.pos.multi) <- c('time', 'distance', 'distance2', 'distance3', 'distance4', 'distance5' , 'distance6', 'distance7', 'distance8', 'distance9', 'distance10', 'distance11')
metadata.df[assays.to.include,]

## 4.4 - Smooth lines

#edge.pos.multi$distance  <- rollmean(edge.pos.df$distance, k = 3, fill = c(edge.pos.df$distance[1], 0, edge.pos.df$distance[nrow(edge.pos.df)]))


## 4.5 - Create label data to plot
labels.to.plot <- metadata.df[assays.to.include,]
labels.to.plot$text <- substr(labels.to.plot$file.list, 10,17)
labels.to.plot$posx <- rep(32, 11)
labels.to.plot$posy <- unlist(edge.pos.multi[91,][,2:12])
labels.to.plot$color = colors.to.plot
labels.to.plot

### 4.0 - Plot distance over time of multiple fronts
#####################################################################
### Overlay diffusion front distance over time plots to compare them

multi.track.plot <- ggplot(edge.pos.multi, aes(x = time/60)) +
  geom_point(aes(y = distance), color = 'green4', size = 2) +
  geom_point(aes(y = distance2), color = colors.to.plot[2], size = 2) +
#  geom_point(aes(y = distance3), color = colors.to.plot[3], size = 2) +
#  geom_point(aes(y = distance4), color = colors.to.plot[4], size = 2) +
#  geom_point(aes(y = distance5), color = colors.to.plot[5], size = 2) +
#  geom_point(aes(y = distance6), color = colors.to.plot[6], size = 2) +
#  geom_point(aes(y = distance7), color = colors.to.plot[7], size = 2) +
#  geom_point(aes(y = distance8), color = colors.to.plot[8], size = 2) +
#  geom_point(aes(y = distance9), color = colors.to.plot[9], size = 2) +
#  geom_point(aes(y = distance10), color = colors.to.plot[10], size = 2) +
#  geom_point(aes(y = distance11), color = colors.to.plot[11], size = 2) +
  scale_x_continuous(name="Time (Min)", breaks = seq(0,100, 5), limits = c(0,30)) +
  scale_y_continuous(name="Distance from channel (um)", breaks = seq(0,10000, 100), limits = c(0,1400)) +
  labs(title="Diffusion Front Distance from Channel Edge Over Time", 
       y = "Distance from Channel (um)", x="Time (min)") +
  theme(panel.background = element_rect(fill = "grey98", colour = "black"), 
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin(2, 2, 2, 2)) ,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14, margin(2, 2, 2, 2))) 

multi.track.plot

labels.to.plot <- labels.to.plot[1:7,]

multi.track.plot + geom_text(data = labels.to.plot, aes(x = posx, y = posy, label = text, color = color),hjust = 0, size = 4) + 
  scale_color_manual(values = c('green2', 'grey20', 'grey20')) +
  theme(legend.position = "none") 



ggplot(edge.pos.multi, aes(x = time/60)) +
  geom_point(aes(y = distance), color = colors.to.plot[1], size = 2) +#  geom_point(aes(y = distance2), color = colors.to.plot[2], size = 2) +
  geom_point(aes(y = distance2), color = colors.to.plot[2], size = 2) +
  geom_point(aes(y = distance3), color = colors.to.plot[3], size = 2) +
  geom_point(aes(y = distance4), color = colors.to.plot[4], size = 2) +
  geom_point(aes(y = distance5), color = colors.to.plot[5], size = 2) +
  geom_point(aes(y = distance6), color = colors.to.plot[6], size = 2) +
  geom_point(aes(y = distance7), color = colors.to.plot[7], size = 2) +
  geom_point(aes(y = distance8), color = colors.to.plot[8], size = 2) +
  geom_point(aes(y = distance9), color = colors.to.plot[9], size = 2) +
  geom_point(aes(y = distance10), color = colors.to.plot[10], size = 2) +
  geom_point(aes(y = distance11), color = colors.to.plot[11], size = 2) +
  scale_x_continuous(name="Time (Min)", breaks = seq(0,100, 5), limits = c(0,30)) +
  scale_y_continuous(name="Distance from channel (um)", breaks = seq(0,10000, 100), limits = c(0,1600)) +
  labs(title="Diffusion Front Distance from Channel Edge Over Time", 
       y = "Distance from Channel (um)", x="Time (min)") +
  theme(panel.background = element_rect(fill = "grey98", colour = "black"), 
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin(2, 2, 2, 2)) ,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14, margin(2, 2, 2, 2))) 

############################################################################## 
### Scratchwork
############################################################################## 

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

ggplot(edge.pos.multi, aes(x = time/60)) +
  geom_point(aes(y = distance), color = 'green', size = 2) +
  geom_point(aes(y = distance2), color = 'Orange', size = 2) +
  geom_point(aes(y = distance3), color = 'Red', size = 2) +
  geom_point(aes(y = distance4), color = 'Blue', size = 2) +
  geom_point(aes(y = distance5), color = 'Violet', size = 2) +
  geom_point(aes(y = distance7), color = 'Red', size = 2) +
  labs(title="Diffusion front Distance from Channel Edge Over Time", 
       y = "Distance from Channel (um)", x="Time (min)", subtitle = assayID) +
  theme(panel.background = element_rect(fill = "grey95", colour = "black"), 
        panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin(2, 2, 2, 2)) ,
        axis.text=element_text(size=12),
        axis.title=element_text(size=14, margin(2, 2, 2, 2)))

#####################################################################
### 5.0 - Calculate fit to quantify diffusion rate
#####################################################################





### Save plot
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")
png(filename= paste0(assayID,'_front dist over time.png'), width = 800, height = 800, units = "px", pointsize = 12)
plot(dist.over.time)
dev.off()

### Save Diffusion Front Data
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E352 - FITC Analysis/")
#write.csv(edge.pos.df, paste0('Edge_position_tables/',title,'_',toupper(format(Sys.Date(),"%d%b%y")),'.csv'),row.names = FALSE)




