### raster and rgdal demo -- Andrew R Gross -- 2020/05/28
### The script documents basic use of the raster and rgdal packages.

### Header
require(raster)             # Load raster package
require(rgdal)

### Reading Data

setwd("C:/Users/grossar/Pictures/")                                       # Set working directory
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E245 - FITC")   # Set alternative working directory 
list.files()                                                              # List files in working directory 
#list.dirs()                                                              # Alternative commands to list contents
#dir()                                                                    # Alternative commands to list contents

### Testing read commands (readGDAL, raster, stack)
test1 = readGDAL(system.file("pictures/Rlogo.jpg", package = "rgdal")[1]) # Using readGDAL to read raster (if in multi-band stack format) 
test1 = readGDAL("cat-test.tif")
test2 <- stack( raster(system.file("pictures/Rlogo.jpg", package = "rgdal")[1], band=1), # Using raster to create stack from individual bands and coerce to SpatialGridDataFrame
                raster(system.file("pictures/Rlogo.jpg", package = "rgdal")[1], band=2),
                raster(system.file("pictures/Rlogo.jpg", package = "rgdal")[1], band=3))

### Assess objects
str(test1)
str(test1@data)
head(test1)
head(test1@data)
class(test1)
class(test1@data)

test1 <- as(test1, "SpatialGridDataFrame")

# do something to the data and write raster to disk  
test1@data <- test1@data * 0.01  
writeGDAL(y, "corrected.tif", drivername="GTiff", type="Float32") 



### Import Data
#####################################################################
setwd("C:/Users/grossar/Pictures/")                                       # Set working directory
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E245 - FITC")   # Set alternative working directory 
list.files()                                                              # List files in working directory 
test2 <- raster("cat-test.jpg")
test2 <- as.matrix(test2)
heatmap(test2, Rowv = NA, Colv = NA)                                      # A heatmap of the raster data
test2 <- apply(test2, 2, rev)                                             # Flips matrix along y axis


my_palette <- colorRampPalette(c("black", "yellow", "white"))(n = 299)
my_palette <- colorRampPalette(c("black", "blue", "yellow", "white"))(n = 299)

heatmap(test2, Rowv = NA, Colv = NA, col = my_palette)

summary(as.vector(test2))
hist(as.vector(test2), breaks = 50)

as.vector(test2[1,])
plot(dist(test2[1,]))



### Import TIFF file
## TESTING
setwd("C://Users/grossar/Pictures/")
test2 <- raster("cat-test.jpg")
test2 <- as.matrix(test2)
heatmap(test2, Rowv = NA, Colv = NA)                                      # A heatmap of the raster data
test2 <- apply(test2, 2, rev)                                             # Flips matrix along y axis


my_palette <- colorRampPalette(c("black", "yellow", "white"))(n = 299)
heatmap(test2, Rowv = NA, Colv = NA, col = my_palette)

summary(as.vector(test2))
hist(as.vector(test2), breaks = 50)
as.vector(test2[1,])
plot(dist(test2[1,]))


### Import method 2
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E245 - FITC/E245_FITC_02-04/")
channel <- stack("E245_FITC_02-04_A_004.tiff")
channel1 <- raster("E245_FITC_02-04_A_004 kept stack.tiff")
test.stack <- stack("test-stack2.tif")
test.stack2 <- stack("E245_02-4_004/004.ome.tif")

stack.to.process <- channel1






setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E245 - FITC/E245_FITC_02-04/")
channel <- stack("E245_FITC_02-04_A_004.tiff")
channel1 <- raster("E245_FITC_02-04_A_004 kept stack.tiff")
test.stack <- stack("test-stack2.tif")
test.stack2 <- stack("E245_02-4_004/004.ome.tif")

stack.to.process <- channel1

### New Stack loading method
setwd("C://Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E245 - FITC/E245_FITC_02-21/stacktoprocess/")
stack.to.process <- stack(list.files())


### Subset stack
layer.of.interest <- subset(stack.to.process,170)
layer.of.interest <- stack.to.process[[2]]

### This will be the beginning of the loop
### Convert to Matrix
#test3 <- as.matrix(test2)
#channel1m <- as.matrix(channel1)
layer.matrix <- as.matrix(layer.of.interest)+1

### Subset
#test3 <- test3[20:120,20:120]
(dim(layer.matrix))
range = 2400
offset = 0
xpos = nrow(layer.matrix)/2 + offset - (range/2)
layer.matrix.ss <- layer.matrix[xpos:(xpos+range),50:90]
layer.matrix.ss <- log(layer.matrix.ss)+1

dim(layer.matrix.ss)

### Plot as heatmap
#heatmap(test3, Rowv = NA, Colv = NA)

colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))

my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
my_palette <- colorRampPalette(c("black", "green", "black"))(n = 299)

heatmap(layer.matrix, Rowv = NA, Colv = NA, col = my_palette)
heatmap(layer.matrix.ss, Rowv = NA, Colv = NA, col = my_palette, breaks = colors)

### Calcuate row averages
dim(layer.matrix)

row.averages <- apply(layer.matrix.ss,1,mean)

plot(row.averages)

### Calculate cutoff positions
max(row.averages)
(cutoff = max(row.averages)/100)
abline(h = cutoff, col = "blue")

cutoff.vector <- row.averages > cutoff
dye.location <- which(cutoff.vector > cutoff)

### Create a table of edge positions
#edge.pos.df <- data.frame(Time = c(),Left.edge = c(), Right.edge = c(), Left.travel = c(), Right.travel = c(), Avg.travel = c())
edge.pos.df <- data.frame(Time = c(),Left.edge = c(), Right.edge = c())

layer.number = 1
time = 20*layer.number
left.edge = dye.location[1]
right.edge = dye.location[length(dye.location)]

new.row <- data.frame(time, left.edge, right.edge)
edge.pos.df <- rbind(edge.pos.df,data.frame(time, left.edge, right.edge))

### Full loop
layer.num = 17
for (layer.num in 1:nlayers(stack.to.process)){           # Begin loop through the number of slices in the stack
  layer.of.interest <- stack.to.process[[layer.num]]      # Call the current layer in the stack
  layer.matrix <- as.matrix(layer.of.interest)+1
  layer.matrix <- log(layer.matrix)
  ### Calculate row averages
  #dim(layer.matrix)
  row.averages <- apply(layer.matrix,1,mean)
  #plot(row.averages)
  
  ### Calculate cutoff positions
  max(row.averages)
  (cutoff = max(row.averages)/100)
  abline(h = cutoff, col = "blue")
  
  cutoff.vector <- row.averages > cutoff
  dye.location <- which(cutoff.vector > cutoff)
  
  ### Fill data frame
  time = 20*(layer.num-1)
  left.edge = dye.location[1]
  right.edge = dye.location[length(dye.location)]
  if (left.edge == 0) {right.edge = 0}
  new.row <- data.frame(time, left.edge, right.edge)
  edge.pos.df <- rbind(edge.pos.df,data.frame(time, left.edge, right.edge))
  
}

### Plotting the curve front
diffusion.profile.list <- list()
layer.num = 61
layer.of.interest <- stack.to.process[[layer.num]]      # Call the current layer in the stack
layer.matrix <- as.matrix(layer.of.interest)+1
layer.matrix <- log(layer.matrix)
### Calculate row averages
dim(layer.matrix)
row.averages <- apply(layer.matrix,1,mean)
plot(row.averages)

### Add to list
diffusion.profile.list[[6]] <- row.averages
diffusion.df <- data.frame(x = 1:1212, "time0" = diffusion.profile.list[[1]], "time4min" = diffusion.profile.list[[2]], 
                           "time8min" = diffusion.profile.list[[3]], "time12min" = diffusion.profile.list[[4]],
                           "time16min" = diffusion.profile.list[[5]], "time20min" = diffusion.profile.list[[6]] )

### Plotting a 3D contour surface 

library(plotly)
# volcano is a numeric matrix that ships with R
p <- plot_ly(z = ~volcano) %>% add_surface()

p

ggplot(data=diffusion.df, aes(x=x), label = c("t0","t4","t8", "t12", "t16", "t20")) + 
  geom_line(aes(y=time0, col = 't0'), size = 1) + 
  geom_line(aes(y=time4min, col = 't4'), size = 1) +
  geom_line(aes(y=time8min, col = 't8'), size = 1) + 
  geom_line(aes(y=time12min, col = 't12'), size = 1) +
  geom_line(aes(y=time16min, col = 't16'), size = 1) +
  geom_line(aes(y=time20min, col = 't20'), size = 1) +
  theme(panel.background = element_rect(fill = "black", linetype = "blank", colour = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(guide_legend()) +
  scale_color_manual(name="Time", 
                     labels = c("t0"="time0", 
                                "t4" = "4 min", 
                                "t8" = "8 min", 
                                "t12" = "12 min", 
                                "t16" = "16 min",
                                "t20" = "20 min"), 
                     values = c("t0"="red", 
                                "t4"="orange", 
                                "t8"="yellow", 
                                "t12"="green", 
                                "t16"="blue",
                                "t20"="purple"),
                     limits = c("t0", "t4", "t8", "t12", "t16", "t20")) +
  labs(title="Diffusion range over Time", y="Relative Fluorescent Intensity", x="Position (px)")

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

library(plotly)
# volcano is a numeric matrix that ships with R
p <- plot_ly(z = ~volcano) %>% add_surface()

p
