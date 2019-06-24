#####################################################################################
#this script is to process raster and vector data to  look at differences in  #######
#the Biodiversity Intactness Index and the biomass intactness index and compare them#
#to the human footprint index########################################################
#####################################################################################


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
  v <- (hyp * cos(phi) / sqrt(2))
  scaled_v <- v * (dark_grey - light_grey) + light_grey
  # Get hsv values.
  sapply(seq_along(x), function(i) hsv(hue[i], s[i], scaled_v[i]))
  
}

#load packages to be used

libPaths("C:/Program Files/R/R-3.6.0/library")

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


#####################################################################
#2 - Plotting of maps################################################
#####################################################################

spr_df<-read.csv("data/analysis_output/BII_HF_BMI.csv")

summary(spr_df)
spr_df_cc<-spr_df[complete.cases(spr_df), ]

#load in coastline dataset
coast<-readOGR(dsn="data/coastline/ne_10m_coastline.shp",layer="ne_10m_coastline")
#fortify coastline shapefile for plotting in ggplot
coast@data$id<-rownames(coast@data)
coast.points<-fortify(coast,region="id")

#plot bivariate map of biomass intactness vs biodiversity intactness
#alter function 'col_func' to change the tones used in the map

biv_map1<-ggplot()
biv_map2<-biv_map1+geom_raster(data=spr_df,aes(x=x,y=y,fill=biomass,fill2=bii))+
         scale_fill_colourplane(name = "",na.color = "NA",color_projection = col_func,
         limits_y = c(0,1),limits=c(0,1),axis_title = "biomass \nintactness",
         axis_title_y = "Biodiversity\n Intactness \nIndex",breaks = c(0,0.5,1),breaks_y = c(0,0.5,1),
         hue1 = 1, hue2 = 0.62,
         dark_grey = 0.5,
         light_grey = 0.9)+
         theme_bw()+coord_equal(xlim = c(-180,180),ylim = c(-50,80))+
         theme(axis.line=element_blank(),panel.background=element_blank(),panel.border=element_blank(),
         panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
         guides(fill = guide_colorplane(title.theme = element_text(size = 13),label.theme = theme_gray(),label.y.theme = theme_gray()))+
         theme(legend.position=c(.1,.3),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.key.size=unit(0.8,"cm"))
biv_map2<-biv_map1+
  theme_bw()+coord_equal(xlim = c(-180,180),ylim = c(-50,80))+
  theme(axis.line=element_blank(),panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())+
  guides(fill = guide_colorplane(title.theme = element_text(size = 13),label.theme = theme_gray(),label.y.theme = theme_gray()))+
  theme(legend.position=c(.1,.3),panel.background = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),legend.key.size=unit(0.8,"cm"))
biv_map3<-biv_map2+geom_path(data=coast.points,aes(x = long, y = lat, group = group),colour="black",size=0.1)

ggsave("figures/bivariate_map_new_scale_final.png",height = 10,width=14,dpi = 800,units = "in")
ggsave("figures/bivariate_map_blank_presentation.png",height = 10,width=14,dpi = 800,units = "in")




########################################################################
#3- calculate the mean bii, biomass, and human footprint of############
#different ecoregions as well as the overlap between ecoregions and#####
#hotspots###############################################################
########################################################################

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


############################################
#produce final figure#######################
############################################

hotspot_overlap_df<-read.csv("data/hotspot_overlap.csv")
head(ecoregion_analysis)

#first calculate the median and interquartile ranges
hotspot_median<-ddply(hotspot_overlap_df,.(hotspot),summarize,m_bii=median(bii,na.rm=T),m_biomass=median(biomass,na.rm=T),
                      biom_25=quantile(biomass,probs=0.25,na.rm=T),biom_75=quantile(biomass,probs=0.75,na.rm=T),
                      bii_25=quantile(bii,probs=0.25,na.rm=T),bii_75=quantile(bii,probs=0.75,na.rm=T))


hot_plot1_v2<-ggplot(hotspot_overlap_df,aes(x=biomass,y=bii,colour=hotspot))+geom_point(shape=16,alpha=0.4)+geom_abline(lty=2)+coord_cartesian(xlim=c(0,1.1),ylim=c(0,1.1),expand = F)
hot_plot2_v2<-hot_plot1_v2+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("Biodiversity Intactness Index")+xlab("biomass intactness")
hot_plot3_v2<-hot_plot2_v2+geom_point(data=hotspot_median,aes(x=m_biomass,y=m_bii,colour=hotspot),shape=15,size=2,stroke=1)+ theme(legend.position="none")
hot_plot4_v2<-hot_plot3_v2+geom_errorbarh(data=hotspot_median,aes(x=m_biomass,y=m_bii,xmin=biom_25,xmax=biom_75,colour=hotspot),size=0.5,height=0.005)+
  geom_errorbar(data=hotspot_median,aes(x=m_biomass,y=m_bii,ymin=bii_25,ymax=bii_75,colour=hotspot),size=0.5,width=0.005)
save_plot(hot_plot4_v2,filename = "figures/ecoregion_scatter_line_2.png",base_aspect_ratio = 1.2,dpi=600)


#plot hf vs bii for ecoregions
hf_bii_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=bii,colour=hotspot))+geom_point(shape=1,alpha=0.8)
hf_bii_plot2<-hf_bii_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("Biodiversity Intactness Index")+xlab("Human Footprint")
hf_bii_plot3<-hf_bii_plot2+theme(legend.position="none")
save_plot(hf_bii_plot3,filename = "figures/hf_bii_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#plot hf vs biomass for ecoregions
hf_bmi_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=biomass,colour=hotspot))+geom_point(shape=1,alpha=0.8)
hf_bmi_plot2<-hf_bmi_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biomass intactness")+xlab("Human Footprint")
hf_bmi_plot3<-hf_bmi_plot2+theme(legend.position="none")
save_plot(hf_bmi_plot3,filename = "figures/hf_bmi_scatter.png",base_aspect_ratio = 1.2,dpi=600)


bottom_row<-plot_grid(hf_bmi_plot3,hf_bii_plot3,hot_plot4_v2,labels=c("b","c","d"),align='h',ncol=3)
combined_plots<-plot_grid(biv_map3,bottom_row,labels=c('a',''),ncol=1,rel_heights=c(1.6,1))
save_plot(combined_plots,filename = "figures/combined_plot_new.png",base_aspect_ratio = 1.2,dpi=600,base_height = 10,base_width = 12)
