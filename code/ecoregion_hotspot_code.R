########################################################################
#this script calculates the mean bii, biomass, and human footprint of###
#different ecoregions as well as the overlap between ecoregions and#####
#hotspots###############################################################
########################################################################


rm(list = ls())

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
library(cowplot)
library(reshape2)
library(dplyr)

#load in ecoregion  dataset
eco_regions <- shapefile("data/ecoregions/Ecoregions2017.shp")
#load in hotspots dataset
hotspots <- shapefile("data/hotspots/hotspots_revisited_2004_polygons.shp")
#load in coastline dataset
coast<-readOGR("data",layer="ne_10m_coastline")
#load in BII map data
BII_map<-raster("data/bii/lbii.asc")
#aggregate the BII map to a coarser resolution and project
BII_agg<-aggregate(BII_map, fact=10, fun=mean, expand=TRUE)
proj4string(BII_agg) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#extract information on BII values in each ecoregion
bii_ecoregions<-extract(BII_agg,eco_regions,fun=mean,na.rm=TRUE)

#load in gridded population data
pop_dens<-raster("data/population_density/gpw_v4_population_density_rev10_2010_2pt5_min.tif")
plot(pop_dens)

#extact population values in each ecoreion
pop_ecoregions<-extract(pop_dens,eco_regions,fun=mean,na.rm=TRUE)

#load in agriculture data
pasture<-raster("data/agriculture/pasture.tif")
cropland<-raster("data/agriculture/cropland.tif")
pasture_ecoregions<-extract(pasture,eco_regions,fun=mean,na.rm=TRUE)
cropland_ecoregions<-extract(cropland,eco_regions,fun=mean,na.rm=TRUE)




hotspot_overlap_df<-read.csv("data/hotspot_overlap.csv")
head(hotspot_overlap_df)
hotspot_overlap_df$bii<-bii_ecoregions
hotspot_overlap_df$pop_den<-pop_ecoregions
hotspot_overlap_df$pasture<-pasture_ecoregions
hotspot_overlap_df$cropland<-cropland_ecoregions

write.csv(hotspot_overlap_df,"data/hotspot_overlap.csv")
require(scales)
ggplot(hotspot_overlap_df,aes(x=pop_den,y=bii))+geom_point()+scale_x_continuous(trans='log10',labels = comma)+geom_smooth(method="glm")
ggplot(hotspot_overlap_df,aes(x=pop_den,y=biomass))+geom_point()+scale_x_continuous(trans='log10',labels = comma)+geom_smooth(method="glm")

ggplot(hotspot_overlap_df,aes(x=pasture,y=biomass))+geom_point()+geom_smooth(method="glm")
ggplot(hotspot_overlap_df,aes(x=cropland,y=biomass))+geom_point()+geom_smooth(method="glm")
ggplot(hotspot_overlap_df,aes(x=cropland+pasture,y=biomass,size=ecoregion_area))+geom_point()+geom_smooth()
ggplot(hotspot_overlap_df,aes(x=cropland+pasture,y=bii))+geom_point()+geom_smooth()


#load in biomass dataset
biomass_corr<-raster("data/biomass_corrected.tif")
#load in human footprint data
hfp<-raster("data/HumanFootprintv2/Dryadv3/Maps/HFP2009.tif")
#load elements of HF layer
lights<-raster("data/HumanFootprintv2/Dryadv3/Maps/Lights2009.tif")
pasture<-raster("data/HumanFootprintv2/Dryadv3/Maps/Pasture2009.tif")
croplands<-raster("data/HumanFootprintv2/Dryadv3/Maps/croplands2005.tif.tif")
pop<-raster("data/HumanFootprintv2/Dryadv3/Maps/Popdensity2010.tif")


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

########################################################################
#2 - Extract data from BII, biomass, and human footprint rasters to ####
#using shapefile of ecoregions##########################################
########################################################################

#aggregate the BII map to a coarser resolution and project
BII_agg<-aggregate(BII_map, fact=10, fun=mean, expand=TRUE)
proj4string(BII_agg) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#extract information on BII values in each ecoregion
bii_ecoregions<-extract(BII_agg,eco_regions,fun=mean,na.rm=TRUE)
#extract information on biomass values in each ecoregion
biomass_ecoregions<-extract(biomass_corr,eco_regions,fun=mean,na.rm=TRUE)
#extract information on human footprint in each ecoregion
hf_ecoregions<-extract(hfp,eco_regions,fun=mean,na.rm=TRUE)


#remove nas
eco_regions2@data$hotspot_area2<-ifelse(is.na(eco_regions2@data$hotspot_area),0,eco_regions2@data$hotspot_area)
#put all this in a dataframe
hotspot_overlap_df<-data.frame(ECO_NAME=eco_regions2@data$ECO_NAME,hotspot_area=eco_regions2@data$hotspot_area2,bii=bii_ecoregions,biomass=biomass_ecoregions,hf=hf_ecoregions)
hotspot_overlap_df$hotspot<-ifelse(hotspot_overlap_df$hotspot_area>=50,"hotspot","not hotspot")
#save dataframe to csv
write.csv(hotspot_overlap_df,"data/hotspot_overlap.csv")

#get statistics for hotspots
#get mean bii and biomass for each hotspot
hotspot_biomass<-extract(biomass_corr,hotspots,fun=mean,na.rm=TRUE)
hotspot_bii<-extract(bii_corr,hotspots,fun=mean,na.rm=TRUE)
#check the different that using the uncorrected version of bii makes
hotspot_bii_unc<-extract(BII_agg,hotspots,fun=mean,na.rm=TRUE)


hotspot_stats<-data.frame(hotspot_name=hotspots@data$NAME,biomass=hotspot_biomass,bii=hotspot_bii,bii_unc=hotspot_bii_unc)
ggplot(hotspot_stats,aes(x=hotspot_bii,y=hotspot_bii_unc))+geom_point()+geom_abline(lty=2)
write.csv(hotspot_stats,"data/hotspot_stats.csv")

#plot some exploratory plots


median(hotspot_overlap_df$bii,na.rm=T)
median(hotspot_overlap_df$biomass,na.rm=T)
quantile(hotspot_overlap_df$biomass,na.rm=T)
quantile(hotspot_overlap_df$bii,probs=0.25,na.rm=T)


#plot alternative version with line of unity and mean values for hotspots and non-hotspots

#first calculate the median and interquartile ranges
hotspot_median<-ddply(hotspot_overlap_df,.(hotspot),summarize,m_bii=median(bii,na.rm=T),m_biomass=median(biomass,na.rm=T),
                      biom_25=quantile(biomass,probs=0.25,na.rm=T),biom_75=quantile(biomass,probs=0.75,na.rm=T),
                      bii_25=quantile(bii,probs=0.25,na.rm=T),bii_75=quantile(bii,probs=0.75,na.rm=T))


hot_plot1_v2<-ggplot(hotspot_overlap_df,aes(x=biomass,y=bii,colour=hotspot))+geom_point(shape=16,alpha=0.4)+geom_abline(lty=2)+coord_cartesian(xlim=c(0,1.1),ylim=c(0,1.1),expand = F)
hot_plot2_v2<-hot_plot1_v2+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biodiversity intactness index")+xlab("biomass intactness")
hot_plot3_v2<-hot_plot2_v2+geom_point(data=hotspot_median,aes(x=m_biomass,y=m_bii,colour=hotspot),shape=15,size=2,stroke=1)+ theme(legend.position="none")
hot_plot4_v2<-hot_plot3_v2+geom_errorbarh(data=hotspot_median,aes(x=m_biomass,y=m_bii,xmin=biom_25,xmax=biom_75,colour=hotspot),size=0.5,height=0.005)+
  geom_errorbar(data=hotspot_median,aes(x=m_biomass,y=m_bii,ymin=bii_25,ymax=bii_75,colour=hotspot),size=0.5,width=0.005)
save_plot(hot_plot4_v2,filename = "figures/ecoregion_scatter_line_2.png",base_aspect_ratio = 1.2,dpi=600)


#plot hf vs bii for ecoregions
hf_bii_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=bii,colour=hotspot))+geom_point(shape=1,alpha=0.8)
hf_bii_plot2<-hf_bii_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biodiversity intactness index")+xlab("human footprint")
hf_bii_plot3<-hf_bii_plot2+theme(legend.position="none")
save_plot(hf_bii_plot3,filename = "figures/hf_bii_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#plot hf vs biomass for ecoregions
hf_bmi_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=biomass,colour=hotspot))+geom_point(shape=1,alpha=0.8)
hf_bmi_plot2<-hf_bmi_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biomass intactness")+xlab("human footprint")
hf_bmi_plot3<-hf_bmi_plot2+theme(legend.position="none")
save_plot(hf_bmi_plot3,filename = "figures/hf_bmi_scatter.png",base_aspect_ratio = 1.2,dpi=600)
