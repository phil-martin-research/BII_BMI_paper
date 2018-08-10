#this script is to calculate statistics from rasters using a grid

#load packages to be used
library(raster)
library(rgdal)
library(ggplot2)
library(maps)
library(sp)
library(colorplaner)
library(plyr)
library(spatialEco)
library(cowplot)
library(viridis)


#create a grid
r <- raster(extent(matrix( c(-180, -90, 180,  90), nrow=2)), nrow=1800, ncol=3600, 
            crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")            
r[] <- 1:ncell(r)
sp.r <- as(r, "SpatialPixelsDataFrame")
summary(sp.r)

sp.r_sub<-subset(sp.r,layer<100)

plot(sp.r_sub)

#load in raster datasets

#load in BII map data
BII_map<-raster("data/bii/lbii.asc")
#load in biomass map
biomass<-raster("data/biomass_intactness/biomass_corrected.tif")
#load in human footprint index
hfp<-raster("data/HumanFootprintv2/Dryadv3/Maps/HFP2009.tif")
plot(hfp)
#load in nightlights, pasture, croplands, and human population density data
pop<-raster("data/population_density/gpw_v4_population_density_rev10_2010_2pt5_min.tif")
lights<-raster("data/lights/F182013.v4c_web.stable_lights.avg_vis.tif")
pasture<-raster("data/agriculture/pasture.tif")
croplands<-raster("data/agriculture/cropland.tif")
agriculture<-pasture+croplands
#stack rasters
var_stack<-stack(BII_map,biomass,hfp,pop,lights,pasture,croplands,agriculture)

res(BII_map)
res(biomass)
res(hfp)
res(pop)
res(lights)
res(pasture)
res(croplands)
res(agriculture)

#extract raster values

sp.r$biomass<-(extract(biomass,sp.r))
#extract values from BII data
sp.r$bii<-(extract(BII_map,sp.r))
#extract hfp values to grid
sp.r$hfp<-(extract(projected_hfp,sp.r))
#extract croplands to grid
sp.r$crops<-(extract(croplands,sp.r))
#extract pastures to grid
sp.r$pasture<-(extract(pasture,sp.r))
#extract agriculture to grid
sp.r$agriculture<-(extract(agriculture,sp.r))
# extract lights to grid
sp.r$lights<-(extract(lights,sp.r))
#convert grid to a dataframe for plotting in  ggplot
spr_df<-as.data.frame(sp.r)
#save this data as a csv file
write.csv(spr_df,"data/Grid_analysis.csv")

ggplot(data=spr_df,aes(x=agriculture,y=biomass))+geom_point()+geom_smooth()
