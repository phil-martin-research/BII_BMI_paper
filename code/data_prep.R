#this script is for data formatting and extraction from raster layers

rm(list = ls())

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
library(snow)

beginCluster(n=4,type="SOCK") 

###################################################################
#1 - Processsing of spatial data###################################
###################################################################

#load in BII map data
BII_map<-raster("data/bii/lbii.asc")
#load in biomass map
biomass<-raster("data/biomass_intactness/C_reduction_perc.tif")
#load in nightlights, pasture, croplands, and human population density data
pop<-raster("data/population_density/gpw_v4_population_density_rev10_2010_2pt5_min.tif")
lights<-raster("data/lights/F182013.v4c_web.stable_lights.avg_vis.tif")
pasture<-raster("data/agriculture/pasture.tif")
croplands<-raster("data/agriculture/cropland.tif")

#project all data to mollweide equal area projection

#vegetation biomass
biomass_ex <- projectExtent(biomass, CRS("+proj=moll +lon_0=0 +ellps=WGS84"))
biomass_proj <- projectRaster(biomass, biomass_ex)
#clean up biomass data so that anything with values <0 and >1 get a a value of NA
values(biomass_proj)[values(biomass_proj) <0] = NA
values(biomass_proj)[values(biomass_proj) >=1] = NA
#invert values for raster so that biomass raster represents the intactness of biomass rather than it's loss
biomass_inv<-raster.invert(biomass_proj)
#save this altered version of the biomass map
writeRaster(biomass_inv, "data/biomass_intactness/biomass_corrected_moll", format = "GTiff",overwrite=T)

#human population density
pop_ex <- projectExtent(pop, CRS("+proj=moll +lon_0=0 +ellps=WGS84"))
pop_proj <- projectRaster(pop, pop_ex)
writeRaster(pop_proj, "data/population_density/pop__moll", format = "GTiff",overwrite=T)

#night light data
light_resample<-resample(lights,biomass)
lights_ex<-projectExtent(light_resample, CRS("+proj=moll +lon_0=0 +ellps=WGS84"))
lights_proj<-projectRaster(lights, lights_ex)
writeRaster(lights_proj, "data/lights/lights_moll", format = "GTiff",overwrite=T)

#pasture data
pasture_ex <- projectExtent(pasture, CRS("+proj=moll +lon_0=0 +ellps=WGS84"))
pasture_proj <- projectRaster(pasture, pasture_ex)
writeRaster(pasture_proj, "data/agriculture/pasture_moll", format = "GTiff",overwrite=T)

#croplands data
croplands_ex <- projectExtent(croplands, CRS("+proj=moll +lon_0=0 +ellps=WGS84"))
croplands_proj <- projectRaster(croplands, croplands_ex)
writeRaster(croplands_proj, "data/agriculture/croplands_moll", format = "GTiff",overwrite=T)

#aggregate the BII map to a coarser resolution and project
BII_agg<-resample(BII_map,biomass)
BII_ex <- projectExtent(BII_agg, CRS("+proj=moll +lon_0=0 +ellps=WGS84"))
BII_proj <- projectRaster(BII_agg, BII_ex)
writeRaster(BII_proj, "data/BII_moll", format = "GTiff")
values(BII_proj)[values(BII_proj) >=1] = 1
#save this altered version of the BII map
writeRaster(BII_proj, "data/BII_corrected", format = "GTiff",overwrite=T)

######################################################################
#2. - EXTRACT DATA FROM RASTERS TO A GRID FOR ANALYSIS################
######################################################################

rm(list = ls())
gc()

beginCluster(n=4,type="SOCK") 

biomass<-raster("data/biomass_intactness/biomass_corrected_moll.tif")
BII<-raster("data/bii/BII_moll.tif")
pop<-raster("data/population_density/pop__moll.tif")
lights<-raster("data/lights/lights_moll.tif")
pasture<-raster("data/agriculture/pasture_moll.tif")
croplands<-raster("data/agriculture/croplands_moll.tif")


#create a grid to extract data to with a resolution of 0.083333 degrees
r <- raster(extent(biomass), nrow=2169, ncol=4337,
            crs = "+proj=moll +lon_0=0 +ellps=WGS84 +units=m +no_defs")
r[] <- 1:ncell(r)
sp.r <- as(r, "SpatialPixelsDataFrame")

#extract values from biomass raster
sp.r$biomass<-(extract(biomass,sp.r))
#extract values from BII data
sp.r$bii<-(extract(BII,sp.r))
#extract croplands to grid
sp.r$crops<-(extract(croplands,sp.r))
#extract pastures to grid
sp.r$pasture<-(extract(pasture,sp.r))
# extract lights to grid
sp.r$lights<-(extract(lights,sp.r))
#convert grid to a dataframe for plotting in  ggplot
spr_df<-as.data.frame(sp.r)
#save this data as a csv file
write.csv(spr_df,"data/BII_BMI_data.csv")


