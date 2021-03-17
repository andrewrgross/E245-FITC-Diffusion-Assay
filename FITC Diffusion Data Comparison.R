### FITC Diffusion Data Comparison -- Andrew R Gross -- 2021/03/10
### The script reads in tables generated in the FITC Diffusion Image Processing script
### And plots the positions releative to one antoher for comparison purposes

### 1.0 - Header
#####################################################################
#require(raster)             # Load raster package
#require(rgdal)
#library(plotly)
library(ggplot2)
#library(zoo)
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

## 2.1 - Report data info
##########################

### 3.0 - Plot distance over time of individual diffusion fronts
#####################################################################
### Plot the left and right side distance over time at the cutoff

table.number <- 2
(assayID <- names(diffusion.data.list)[table.number])
edge.pos.df <- diffusion.data.list[[table.number]]

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

### 4.0 - Generate a table of fronts to plot
#####################################################################
### Combine columns from multiple tables in a single table

## 4.1 - List the assay files to include
assays.to.include <- c(1,2)
print(data.frame(file.list[assays.to.include]))

## 4.2 - Create a list of all times to include
times.to.include = c() 
for(assay in assays.to.include){
  times.to.include <- c(times.to.include,diffusion.data.list[[assay]]$time)
}
times.to.include <- sort(unique(times.to.include))

edge.pos.multi <- data.frame('time' = times.to.include)

## 4.3 - For each assay file, merge the column of interest to the full edge.pos.multi df

assay = 2
current.df = current.df <- diffusion.data.list[[assay]][c(2,8)]
edge.pos.multi = merge(edge.pos.multi, current.df, by = 'time', all = TRUE, suffixes = c('', assay))

head(edge.pos.multi)

for(assay in assays.to.include){
  print(names(diffusion.data.list)[assay])
  current.df <- diffusion.data.list[[assay]]
  times.to.add <- match(times.to.include, current.df$time)
  distances.to.add <- current.df$distance[]
  distances.to.add.df <- data.frame(distances = distances.to.add)
  #names(distances.to.add.df) <- names(diffusion.data.list)[assay]
  edge.pos.multi[ncol(edge.pos.multi)+1] <- distances.to.add.df
}

head(edge.pos.multi)


times.to.add <- match(current.df$time,times.to.include)

current.df <- diffusion.data.list[[assay]]
distances.to.add <- match(times.to.include, current.df$time)
distances.to.add.df <- data.frame(distances = distances.to.add)
names(distances.to.add.df) <- names(diffusion.data.list)[assay]
edge.pos.multi[ncol(edge.pos.multi)+1] <- distances.to.add.df


## 4.4 - Join all prepared columns with the universal time column





### 4.0 - Plot distance over time of multiple fronts
#####################################################################
### Overlay diffusion front distance over time plots to compare them

test <- dist.over.time

ggplot(edge.pos.multi, aes(x = time/60)) +
  geom_line(aes(y = distance), color = 'Blue', size = 2) +
  geom_line(aes(y = distance2), color = 'Red', size = 2)



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




