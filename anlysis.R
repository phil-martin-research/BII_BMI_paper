#this is a script to look at the differences in biomass intactness and the biodiversity intactness index

rm(list = ls())

#use this function to produce a hsv bivariate palette for the map
#set hue as 0.15 and 0.62 for yellow to blue pallete
#set hue as 0.07 and 0.62 for orange to blue pallete
#set hue as 1 and 0.62 for red to blue palette
#set hue as 0.10 and 0.29 for orange to green palette
#set hue as 0.82 and 0.29 for pink to green palette

col_func <- function(x, y){
  
  x[x == 0] <- 0.000001
  y[y == 0] <- 0.000001
  x[x == 1] <- 0.999999
  y[y == 1] <- 0.999999
  
  # upper or lower triangle?
  u <- y > x
  
  # Change me for different hues.
  hue <- ifelse(u, 0.82, 0.29)
  
  
  # distance from (0,0) to (x,y)
  hyp <- sqrt(x^2 + y^2) 
  
  # Angle between x axis and line to our point
  theta <- asin(y / hyp)
  
  # Angle between 45 degree line and (x,y)
  phi <- ifelse(u, theta - pi/4, pi/4 - theta)
  phi <- ifelse(phi < 0, 0, phi)
  
  # Distance from 45 degree line and (x,y)
  s <- hyp * sin(phi) / sqrt(2)
  
  # Draw line from (x, y) to 45 degree line that is at right angles.
  # How far along 45 degree line, does that line join.
  v <- 0.9 - hyp * cos(phi) / sqrt(2)
  
  # Get hsv values.
  sapply(seq_along(x), function(i) hsv(hue[i], s[i], v[i],alpha = 1))
  
}


#load packages to be used
library(raster)
library(rgdal)
library(ggplot2)
library(raster)
library(maps)
library(sp)
library(colorplaner)
library(plyr)

#load in biomes  data
shape <- readOGR(dsn = "data", layer = "Ecoregions2017")
#load in coastlien data
coast<-readOGR("data",layer="ne_10m_coastline")

#subset biomes data so that it only includes tropical forests
trop_shape<-subset(shape,BIOME_NAME=="Tropical & Subtropical Moist Broadleaf Forests")

#subset biome data to only include forest biomes
forest_shape<-subset(shape,BIOME_NAME=="Tropical & Subtropical Moist Broadleaf Forests"|
                       BIOME_NAME=="Mediterranean Forests, Woodlands & Scrub"|
                       BIOME_NAME=="Boreal Forests/Taiga"|
                       BIOME_NAME=="Temperate Conifer Forests"|
                       BIOME_NAME=="Temperate Conifer Forests"|
                       BIOME_NAME=="Temperate Broadleaf & Mixed Forests"|
                       BIOME_NAME=="Mangroves"|
                       BIOME_NAME=="Tropical & Subtropical Dry Broadleaf Forests"|
                       BIOME_NAME=="Tropical & Subtropical Coniferous Forests")

#load in BII map data
BII_map<-raster("data/lbii.asc")

#aggregate this map to a coarser resolution and project
BII_agg<-aggregate(BII_map, fact=10, fun=mean, expand=TRUE)
proj4string(BII_agg) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


#mask data to only include BII data from tropical forest biomes
bii_trop<-mask(BII_agg,trop_shape)

#mask data to only include BII data from forest biomes
bii_forest<-mask(BII_agg,forest_shape)

#load in biomass map
biomass<-raster("data/NATURE_DATADOWNLOAD/C_reduction_perc.tif")

#clean up biomass data so that anything with values <0 and >1 get a a value of NA
values(biomass)[values(biomass) <0] = NA
values(biomass)[values(biomass) >=1] = NA

#standardise values for BII so max values=1 and min=0
values(BII_agg)[values(BII_agg)>1]=1

proj4string(biomass) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
biomass<-crop(biomass,BII_agg)

biomass_inv<-raster.invert(biomass)
biom_trop<-mask(biomass_inv,trop_shape)
biom_forest<-mask(biomass_inv,forest_shape)

plot(biomass)
plot(biomass_inv)
plot(biom_forest)


#create a grid to extract data to
r <- raster(extent(matrix( c(-180, -90, 180,  90), nrow=2)), nrow=1800, ncol=1800, 
            crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")            
r[] <- 1:ncell(r)
summary(r)
print(r)
plot(r)
sp.r <- as(r, "SpatialPixelsDataFrame")
class(sp.r)   
 

#extract values from biomass raster
sp.r$biomass<-(extract(biomass_inv,sp.r))
hist(sp.r$biomass)

#extract values from BII data
sp.r$bii<-(extract(BII_agg,sp.r))
sp.r$map<-1
sp.r$map2<-ifelse(sp.r$bii>0,1,0)

spr_df<-as.data.frame(sp.r)

#fortify coastline shapefile for plotting
coast@data$id<-rownames(coast@data)
coast.points<-fortify(coast,region="id")
coast.df<-join(coast.points, coast@data, by="id")


head(world)
biv_map1<-ggplot()
biv_map2<-biv_map1+geom_raster(data=spr_df,aes(x=x,y=y,fill=biomass,fill2=bii))+
         scale_fill_colourplane(name = "",na.color = "NA",
                         color_projection = col_func,
                         limits_y = c(0.5,1),
                         limits=c(0,1),
                         axis_title = "biomass \nintactness",
                         axis_title_y = "biodiversity\n intactness \nindex",
                         breaks = c(0,0.5,1),
                         breaks_y = c(0.5,0.75,1))+
  theme_bw()+
  coord_equal()+
  theme(axis.line=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())+
  guides(fill = guide_colorplane(
    title.theme = element_text(size = 13),
    label.theme = theme_gray(),
    label.y.theme = theme_gray()))+
  theme(legend.position = "bottom",
  panel.background = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank())


biv_map2+geom_path(data=coast.points,aes(x = long, y = lat, group = group),colour="black",size=0.1)

ggsave("figures/bivariate_map5_1.png",height = 10,width=14,dpi = 800,units = "in")



#plot map with a different colour scheme
ggplot(spr_df,aes(x=x,y=y,fill=biomass,fill2=bii))+
  geom_raster()+
  scale_fill_colourplane(name = "",na.color = "white",
                         color_projection = "interpolate",
                         vertical_color = "#E6A3D0",
                         horizontal_color = "#8BE2AF",
                         zero_color = " light grey",
                         axis_title = "biomass \nintactness",
                         axis_title_y = "biodiversity\n intactness \nindex",
                         limits_y = c(0.5,1),
                         limits=c(0,1),
                         breaks = c(0,0.5,1),
                         breaks_y = c(0.5,0.75,1))+
  theme_bw()+
  coord_equal()+
  theme(axis.line=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("figures/bivariate_map_colour_scheme_2.png",height = 10,width=14,dpi = 400,units = "in")




#define new colour scheme
# Define a custom projection function that uses HSV color space by
# mapping the data to hue and saturation with a fixed level for value
hsv_projection <- function(x, y, v) {
  # Convert y value from a position to a hue angle
  h <- atan(scales::rescale(1, from = c(0,1), to = c(-1, 1)) / 
              scales::rescale(2, from =  c(0, 1), to = c(-1, 1)))
  
  
  ?rescale
  
  
  h <- atan(scales::rescale(y, from = c(0,1), to = c(-1, 1)) / 
              scales::rescale(x, from =  c(0, 1), to = c(-1, 1)))
  # There are no missing values in the input, but atan can create some
  h <- ifelse(is.na(h), 0, h)
  h <- scales::rescale(h, from = c(-pi / 2, pi / 2), to = c(0, 1))
  # hsv takes inputs on a scale of [0, 1] and returns an RGB color string
  grDevices::hsv(h, x, v)
}


ggplot(spr_df,aes(x=x,y=y,fill=biomass,fill2=bii))+
  geom_raster()+
  scale_fill_colorplane(color_projection = hsv_projection,v=1)+
  theme_bw()+
  coord_equal()+
  theme(axis.line=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())




#minus one raster from another
BII_diff<-biomass-BII_agg


#function to produce a colour matrix for maps
colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}


col.matrix<-colmat(nquantiles=6,xlab = "biomass intactness",ylab = "biodiversity intactness index")

#function to produce bivariate map
rastery<-biomass_inv
bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-seq(0,1,by=0.1)
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-seq(0.5,1,by=0.1)
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}

#plot of BII map
my.colors = colorRampPalette(c("white","lightblue", "yellow","orangered", "red"))
plot(BII_agg,frame.plot=F,axes=F,box=F,add=F,legend.width=1,legend.shrink=1,col=my.colors(255)) 
map(add=T)

#plot of biomass map
my.colors = colorRampPalette(c("white","lightblue", "yellow","orangered", "red"))
plot(biomass,frame.plot=F,axes=F,box=F,add=F,legend.width=1,legend.shrink=1,col=my.colors(255)) 
map(add=T)


bivmap<-bivariate.map(biomass_inv,BII_agg, colormatrix=col.matrix, nquantiles=10)

plot(bivmap,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
plot(BII_agg)
plot(biomass_inv)

col.matrix
