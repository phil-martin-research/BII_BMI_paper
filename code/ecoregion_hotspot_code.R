#script to do the gis analysis to look at mean bii and mean biomass in 
#different ecoregions as well as overlap beween ecoregions and hotspots

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
library(rgdal)

#load in ecoregion  dataset
eco_regions <- readOGR(dsn = "data", layer = "Ecoregions2017")
#load in hotspots dataset
hotspots <- shapefile("data/hotspots_unzipped/commondata/data0/hotspots_revisited_2004_polygons.shp")
#load in coastline dataset
coast<-readOGR("data",layer="ne_10m_coastline")
#load in bii dataset
bii_corr<-raster("data/BII_corrected.tif")
#load in biomass dataset
biomass_corr<-raster("data/biomass_corrected.tif")

#subset hotspot data to only include extent
hotspots <- hotspots[hotspots$TYPE=="hotspot_area",]

#work out area of hotspot that overlaps with ecoregions
head(eco_regions$ECO_NAME)
#do this in a loop
ecoregions_unique<-eco_regions@data$ECO_NAME
ecoregion_overlap<-data.frame(ECO_NAME=character(),hotspot_area=numeric())
for (i in 1:length(ecoregions_unique)){
  tryCatch({
  ecoregion_sub <-eco_regions[eco_regions$ECO_NAME==ecoregions_unique[i],]
  area_intersect<-(area(intersect(ecoregion_sub,hotspots))/area(ecoregion_sub))*100
  ecoregion_subset<-data.frame(ECO_NAME=ecoregions_unique[i],hotspot_area=ifelse(exists("area_intersect"),area_intersect,0))
  ecoregion_overlap<-rbind(ecoregion_overlap,ecoregion_subset)
  remove(area_intersect)
  print(paste("percent done:",round((i/length(ecoregions_unique)*100),digits = 1)))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.csv(ecoregion_overlap,file="data/ecoregion_overlap.csv")

#attach this data to the ecoregion shapefile data
eco_regions2<-merge(eco_regions,ecoregion_overlap,by="ECO_NAME",all=T)

#extract information on BII values in each ecoregion
eco_regions_subset <- eco_regions[eco_regions2@data$OBJECTID<20,]
eco_regions_subset@data$ECO_NAME

bii_ecoregions<-extract(bii_corr,eco_regions_subset,fun=mean,na.rm=TRUE)
#extract information on biomass values in each ecoregion
biomass_ecoregions<-extract(biomass_corr,eco_regions_subset,fun=mean,na.rm=TRUE)

#remove nas
eco_regions2@data$hotspot_area2<-ifelse(is.na(eco_regions2@data$hotspot_area),0,eco_regions2@data$hotspot_area)
#put all this in a dataframe
hotspot_overlap_df<-data.frame(ECO_NAME=eco_regions2@data$ECO_NAME,hotspot_area=eco_regions2@data$hotspot_area2)


write.csv(hotspot_overlap_df,"data/hotspot_overlap.csv")
#plot some exploratory plots
ggplot(hotspot_overlap_df,aes(x=hotspot_area))+geom_histogram()


