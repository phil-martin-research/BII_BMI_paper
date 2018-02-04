#replicable example for stackoverflow

#load packages
library(raster)
library(colorplaner)
library(ggplot2)


#here's my dummy data
r1<- raster(ncol=10, nrow=10)
set.seed(0)
values(r1) <- runif(ncell(r1))

r2<- raster(ncol=10, nrow=10)
values(r2) <- runif(ncell(r2))

plot(r1)
plot(r2)


grid<-raster(ncol=10, nrow=10)
grid[] <- 1:ncell(grid)
grid.pdf<-as(grid, "SpatialPixelsDataFrame")

grid.pdf$r1<-(extract(r1,grid.pdf))
grid.pdf$r2<-(extract(r2,grid.pdf))

grid.df<-as.data.frame(grid.pdf)

#plot a map  of this

ggplot(data=grid.df,aes(x,y,fill=r1,fill2=r2))+geom_raster()+scale_fill_colourplane("")
ggsave("figures/bivariate_eg2.png",width = 6,height = 4,dpi = 200,units = "in")

#This colourscale doesn't really suit my needs - give example of colour scheme that I'd like

#however I'm finding it tricky to modify the colourscheme in the function scale_fill_colourplane

#the closest I can get to the colourscale I want is this
ggplot(data=grid.df,aes(x,y,fill=r1,fill2=r2))+
  geom_raster()+
  scale_fill_colourplane(name = "",na.color = "white",
                         color_projection = "interpolate",
                         vertical_color = "#FAE30C",
                         horizontal_color = "#0E91BE",
                         zero_color = "#E8E6F2",
                         limits_y = c(0,1),
                         limits=c(0,1))
