#this is a script to look at the differences in biomass intactness and the biodiversity intactness index

rm(list = ls())

#use this function to produce a hsv bivariate palette for the map
#set hue as 0.16 and 0.62 for yellow to blue pallete
#set hue as 0.07 and 0.62 for orange to blue pallete
#set hue as 1 and 0.62 for red to blue palette
#set hue as 0.10 and 0.29 for orange to green palette
#set hue as 0.82 and 0.29 for pink to green palette

col_func <- function(x, y, hue1 = 0.3, hue2 = 0.8, light_grey = 0, dark_grey = 1){
  
  x[x == 0] <- 0.000001
  y[y == 0] <- 0.000001
  x[x == 1] <- 0.999999
  y[y == 1] <- 0.999999
  
  # upper or lower triangle?
  u <- y > x
  hue <- ifelse(u, hue1, hue2)
  
  
  # distace from (0,0) to (x,y)
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
  v <- (hyp * cos(phi) / sqrt(2))
  scaled_v <- v * (dark_grey - light_grey) + light_grey
  # Get hsv values.
  sapply(seq_along(x), function(i) hsv(hue[i], s[i], scaled_v[i]))
  
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
library(spatialEco)
library(cowplot)

###################################################################
#1 - Processsing of spatial data###################################
###################################################################

#load in ecoregion  dataset
eco_regions <- readOGR(dsn = "data", layer = "Ecoregions2017")
#load in coastline dataset
coast<-readOGR("data",layer="ne_10m_coastline")
#fortify coastline shapefile for plotting in ggplot
coast@data$id<-rownames(coast@data)
coast.points<-fortify(coast,region="id")
#load in BII map data
BII_map<-raster("data/lbii.asc")
#load in biomass map
biomass<-raster("data/NATURE_DATADOWNLOAD/C_reduction_perc.tif")
#load in human footprint index
hfp<-raster("data/HumanFootprintv2/Dryadv3/Maps/HFP2009.tif")

#aggregate the BII map to a coarser resolution and project
BII_agg<-aggregate(BII_map, fact=10, fun=mean, expand=TRUE)
proj4string(BII_agg) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
BII_agg_corr<-BII_agg
values(BII_agg_corr)[values(BII_agg_corr) >=1] = 1
#save this altered version of the BII map
writeRaster(BII_agg_corr, "data/BII_corrected", format = "GTiff")

#clean up biomass data so that anything with values <0 and >1 get a a value of NA
values(biomass)[values(biomass) <0] = NA
values(biomass)[values(biomass) >=1] = NA
#set biomass extent and projection to be the sames as that for the BII_agg map
proj4string(biomass) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
biomass<-crop(biomass,BII_agg_corr)
#invert values for raster so that biomass raster represents the intactness of biomass rather than it's loss
biomass_inv<-raster.invert(biomass)
#save this altered version of the biomass map
writeRaster(biomass_inv, "data/biomass_corrected", format = "GTiff",overwrite=T)

#clean up hfp data and reproject
hfp_coord<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
projected_hfp<-projectRaster(hfp, crs = hfp_coord)
projected_hfp_inv<-raster.invert(projected_hfp)

#create a grid to extract data to with a resolution of 0.1 degrees - may need to change this to 0.083333 degrees
r <- raster(extent(matrix( c(-180, -90, 180,  90), nrow=2)), nrow=1800, ncol=3600, 
            crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")            
r[] <- 1:ncell(r)
sp.r <- as(r, "SpatialPixelsDataFrame")

#extract values from biomass raster
sp.r$biomass<-(extract(biomass_inv,sp.r))
#extract values from BII data
sp.r$bii<-(extract(BII_agg_corr,sp.r))
#extract hfp values to grid
sp.r$hfp<-(extract(projected_hfp,sp.r))
sp.r$hfp_inv<-(extract(projected_hfp_inv,sp.r))
sp.r$hfp_norm<-sp.r$hfp/50
sp.r$hfp_inv_norm<-sp.r$hfp_inv/50
#convert grid to a dataframe for plotting in  ggplot
spr_df<-as.data.frame(sp.r)
#save this data as a csv file
write.csv(spr_df,"data/BII_HF_BMI.csv")

#####################################################################
#2 - Plotting of maps################################################
#####################################################################

#plot bivariate map of biomass intactness vs biodivesity intactness
#alter function 'col_func' to change the tones used in the map
biv_map1<-ggplot()
biv_map2<-biv_map1+geom_raster(data=spr_df,aes(x=x,y=y,fill=biomass,fill2=bii))+
         scale_fill_colourplane(name = "",na.color = "NA",color_projection = col_func,
         limits_y = c(0,1),limits=c(0,1),axis_title = "biomass \nintactness",
         axis_title_y = "biodiversity\n intactness \nindex",breaks = c(0,0.5,1),breaks_y = c(0,0.5,1),
         hue1 = 1, hue2 = 0.62,
         dark_grey = 0.5,
         light_grey = 0.9)+
         theme_bw()+coord_equal(ylim=c(-55,75))+
         theme(axis.line=element_blank(),panel.background=element_blank(),panel.border=element_blank(),
         panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
         guides(fill = guide_colorplane(title.theme = element_text(size = 13),label.theme = theme_gray(),label.y.theme = theme_gray()))+
         theme(legend.position=c(.1,.3),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.key.size=unit(0.8,"cm"))
biv_map3<-biv_map2+geom_path(data=coast.points,aes(x = long, y = lat, group = group),colour="black",size=0.1)

ggsave("figures/bivariate_map_new_scale2.png",height = 10,width=14,dpi = 800,units = "in")

#plot bivariate map of human footprint vs biodiversity intactness
biv_map1<-ggplot()
biv_map2<-biv_map1+geom_raster(data=spr_df,aes(x=x,y=y,fill=hfp_norm,fill2=bii))+
  scale_fill_colourplane(name = "",na.color = "NA",color_projection = col_func,
                         limits_y = c(0,1),limits=c(0,1),axis_title = "normalized\nhuman\nfootprint",
                         axis_title_y = "biodiversity\n intactness \nindex",breaks = c(0,0.5,1),breaks_y = c(0,0.5,1),
                         hue1 = 1, hue2 = 0.62,
                         dark_grey = 0.5,
                         light_grey = 0.9)+
  theme_bw()+coord_equal(ylim=c(-55,75))+
  theme(axis.line=element_blank(),panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
  guides(fill = guide_colorplane(title.theme = element_text(size = 13),label.theme = theme_gray(),label.y.theme = theme_gray()))+
  theme(legend.position=c(.1,.3),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.key.size=unit(0.8,"cm"))
biv_map3<-biv_map2+geom_path(data=coast.points,aes(x = long, y = lat, group = group),colour="black",size=0.1)

ggsave("figures/bivariate_map_hfp.png",height = 10,width=14,dpi = 800,units = "in")


#plot bivariate map of human footprint vs biomass intactness
biv_map1<-ggplot()
biv_map2<-biv_map1+geom_raster(data=spr_df,aes(x=x,y=y,fill=hfp_norm,fill2=biomass))+
  scale_fill_colourplane(name = "",na.color = "NA",color_projection = col_func,
                         limits_y = c(0,1),limits=c(0,1),axis_title = "normalized\nhuman\nfootprint",
                         axis_title_y = "biomass\n intactness",breaks = c(0,0.5,1),breaks_y = c(0,0.5,1),
                         hue1 = 1, hue2 = 0.62,
                         dark_grey = 0.5,
                         light_grey = 0.9)+
  theme_bw()+coord_equal(ylim=c(-55,75))+
  theme(axis.line=element_blank(),panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
  guides(fill = guide_colorplane(title.theme = element_text(size = 13),label.theme = theme_gray(),label.y.theme = theme_gray()))+
  theme(legend.position=c(.1,.3),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.key.size=unit(0.8,"cm"))
biv_map3<-biv_map2+geom_path(data=coast.points,aes(x = long, y = lat, group = group),colour="black",size=0.1)

ggsave("figures/bivariate_map_hfp_biomass.png",height = 10,width=14,dpi = 800,units = "in")

head(spr_df)

ggplot(data=spr_df,aes(x=hfp,y=bii))+geom_point(shape=1,alpha=0.01)+geom_smooth(se=F,size=4)
ggplot(data=spr_df,aes(x=hfp,y=biomass))+geom_point(shape=1,alpha=0.01)+geom_smooth(se=F,size=4)
ggplot(data=spr_df,aes(x=biomass,y=bii))+geom_point(shape=1,alpha=0.01)+geom_smooth(se=F,size=4)


cor.test(spr_df$hfp,spr_df$bii)
cor.test(spr_df$biomass,spr_df$bii)



#####################################################################################
#extract information from BII and biomass intactness index for individual ecoregions#
#highlighting the ones that are contained within biodiversity hotspots###############
#####################################################################################

#load in hotspots  dataset
hotspots <- readOGR(dsn = "data/commondata/data0/", layer = "hotspots_revisited_2004_polygons.shp")

#calculate mean  values of bii in ecoregions
BII_vals<-extract(BII_agg,eco_regions)
BII_mean <- lapply(BII_vals, FUN=mean)
eco_regions@data<- data.frame(eco_regions@data, bii=BII_mean)
str(eco_regions@data)

head(eco_regions@data)
#calculate mean values of biomass in ecoregions
biomass_vals<-over(eco_regions,biomass_inv,fn = mean)
biomass_mean <- lapply(r.vals, FUN=mean)
eco_regions@data <- data.frame(eco_regions@data, bii=biomass_mean)
