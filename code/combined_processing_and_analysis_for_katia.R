#####################################################################################
#this script is to process raster and vector data to  look at differences in  #######
#the Biodiversity Intactness Index and the biomass intactness index and compare them#
#to the human footprint index########################################################
#####################################################################################


rm(list = ls())
.libPaths("C:/R/Library")
.Library<-"C:/R/Library"


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


beginCluster(n=6,type="SOCK") 

###################################################################
#1 - Processsing of spatial data###################################
###################################################################

BII_map<-raster("data/lbii.asc")#load in BII map data
biomass<-raster("data/C_reduction_perc.tif")#load in biomass map
hfp<-raster("data/HFP2009.tif")


#project all data to mollweide equal area projection#
biomass_ex <- projectExtent(biomass, CRS("+proj=moll +lon_0=0 +ellps=WGS84"))#vegetation biomass
biomass_proj <- projectRaster(biomass, biomass_ex)
values(biomass_proj)[values(biomass_proj) <0] = NA #set biomass values <0 and >1 to NA
values(biomass_proj)[values(biomass_proj) >=1] = NA
biomass_inv<-raster.invert(biomass_proj)#invert values for raster to represent intactness

writeRaster(biomass_inv, #save this altered version of the biomass map
            "data/biomass_corrected_moll", 
            format = "GTiff",overwrite=T)

BII_agg<-resample(BII_map,biomass)#aggregate the BII map to a coarser resolution
BII_ex <- projectExtent(BII_agg, CRS("+proj=moll +lon_0=0 +ellps=WGS84")) #project BII data
BII_proj <- projectRaster(BII_agg, BII_ex)
writeRaster(BII_proj, "data/BII_moll", format = "GTiff")#save projected BII
BII_proj_capped<-BII_proj
values(BII_proj_capped)[values(BII_proj_capped) >=1] = 1#set maximum BII value as 1
writeRaster(BII_proj_capped, "data/BII_corrected", format = "GTiff",overwrite=T) #save BII map with max value of 1

#create grid to extract data to
r <- raster(extent(biomass), nrow=2169, ncol=4337, #create a grid with resolution 0.083 degrees to extract data to
            crs = "+proj=moll +lon_0=0 +ellps=WGS84 +units=m +no_defs")
r[] <- 1:ncell(r)
sp.r <- as(r, "SpatialPixelsDataFrame")

sp.r$biomass<-(extract(biomass_inv,sp.r))#extract values from biomass raster
sp.r$bii<-(extract(BII_proj,sp.r))#extract values from BII data
spr_df<-as.data.frame(sp.r)#convert grid to a dataframe for plotting in  ggplot
write.csv(spr_df,"data/output_data//BII_BMI_data.csv")#save this data as a csv file
head(spr_df)

########################################################################
#2- calculate the mean bii, biomass, and human footprint of############
#different ecoregions as well as the overlap between ecoregions and#####
#hotspots###############################################################
########################################################################

#load ecoregion, hotspot, and coast shapefiles
eco_regions <-shapefile("data/Ecoregions2017.shp")#load in ecoregion  dataset
hotspots<-shapefile("data/hotspots_revisited_2004_polygons.shp")#load in hotspots dataset

hotspots<-hotspots[hotspots$TYPE=="hotspot_area",]#subset hotspot data to only include polygons of hotspot extent

#work out area of hotspot that overlaps with ecoregions
ecoregions_unique<-eco_regions@data$ECO_NAME#produce list of unique hotspot names
ecoregion_overlap<-data.frame(ECO_NAME=character(),hotspot_area=numeric())#create df with ecoregion name and area of hotspot
for (i in 1:length(ecoregions_unique)){
  tryCatch({
    ecoregion_sub <-eco_regions[eco_regions$ECO_NAME==ecoregions_unique[i],]#subet to ecoregion 'i'
    area_intersect<-(area(intersect(ecoregion_sub,hotspots))
                     /area(ecoregion_sub))*100#calculate percentage of ecoregion occupied by hotspot
    ecoregion_subset<-data.frame(ECO_NAME=ecoregions_unique[i],#if ecoregion overlaps with hotspot area>0
                                 hotspot_area=ifelse(exists("area_intersect"),
                                                     area_intersect,0))
    ecoregion_overlap<-rbind(ecoregion_overlap,ecoregion_subset)#bind this interation to df
    remove(area_intersect)#remove the area calculation
    print(paste("percent done:",round((i/length(ecoregions_unique)*100),digits = 1)))#print progress message
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



#attach this data to the ecoregion shapefile data
eco_regions2<-merge(eco_regions,ecoregion_overlap,by="ECO_NAME",all=T)

#remove nas
eco_regions2@data$hotspot_area2<-ifelse(is.na(eco_regions2@data$hotspot_area),0,eco_regions2@data$hotspot_area)


#extract mean raster values for each ecoregion
bii_ecoregions<-extract(BII_proj,eco_regions,fun=mean,na.rm=TRUE)#extract mean BII for each ecoregion
biomass_ecoregions<-extract(biomass_inv,eco_regions,fun=mean,na.rm=TRUE)#extract mean biomass for each ecoregion
hf_ecoregions<-extract(hfp,eco_regions,fun=mean,na.rm=TRUE) #extract mean human footprint for each ecoregion

#put all this in a dataframe
ecoregion_hotspot_stats<-data.frame(ECO_NAME=eco_regions2@data$ECO_NAME,hotspot_area=eco_regions2@data$hotspot_area2,bii=bii_ecoregions,biomass=biomass_ecoregions,hf=hf_ecoregions)
ecoregion_hotspot_stats$hotspot<-ifelse(ecoregion_hotspot_stats$hotspot_area>=50,"hotspot","not hotspot")
write.csv(ecoregion_hotspot_stats,"data/output_data/ecoregion_hotspot_stats.csv")

#####################################################################
#3 - Plotting of final figure########################################
#####################################################################

spr_df<-read.csv("data/output_data/BII_BMI_data.csv")#read in data extracted to grid
summary(spr_df)#summary of data extracted to grid
mean(spr_df$bii,na.rm=T)

spr_df_cc<-spr_df[complete.cases(spr_df), ]#remove rows with missing data


coast<-readOGR(dsn="data/coastline/ne_10m_coastline.shp",layer="ne_10m_coastline")#load in coastline dataset
coast@data$id<-rownames(coast@data)#fortify coastline shapefile for plotting in ggplot
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
biv_map3<-biv_map2+geom_path(data=coast.points,aes(x = long, y = lat, group = group),colour="black",size=0.1)


hotspot_overlap_df<-read.csv("data/hotspot_overlap.csv")#load ecoregion statistics


#calculate median and quartiles for BII and BMI for ecoregions with >50% overlap with hotspots
#and ecoregions with <50% overlap with hotspots
hotspot_median<-ddply(hotspot_overlap_df,.(hotspot),  
                      summarize,m_bii=median(bii,na.rm=T),
                      m_biomass=median(biomass,na.rm=T),
                      biom_25=quantile(biomass,probs=0.25,na.rm=T),
                      biom_75=quantile(biomass,probs=0.75,na.rm=T),
                      bii_25=quantile(bii,probs=0.25,na.rm=T),
                      bii_75=quantile(bii,probs=0.75,na.rm=T))



#plot figure of relationship between BMI and BII
hs_BMI_BII_1<-ggplot(hotspot_overlap_df,aes(x=biomass,y=bii,colour=hotspot))+
  geom_point(shape=16,alpha=0.4)+
  geom_abline(lty=2)+coord_cartesian(xlim=c(0,1.1),ylim=c(0,1.1),expand = F)
hs_BMI_BII_2<-hs_BMI_BII_1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+
  ylab("Biodiversity Intactness Index")+
  xlab("biomass intactness")
hs_BMI_BII_3<-hs_BMI_BII_2+geom_point(data=hotspot_median,aes(x=m_biomass,y=m_bii,
             colour=hotspot),shape=15,size=2,stroke=1)+theme(legend.position="none")
hs_BMI_BII_4<-hs_BMI_BII_3+geom_errorbarh(data=hotspot_median,aes(x=m_biomass,y=m_bii,
              xmin=biom_25,xmax=biom_75,colour=hotspot),size=0.5,height=0.005)+
  geom_errorbar(data=hotspot_median,aes(x=m_biomass,y=m_bii,
  ymin=bii_25,ymax=bii_75,colour=hotspot),size=0.5,width=0.005)
save_plot(hs_BMI_BII_4,filename = "figures/ecoregion_scatter_line_2.png",base_aspect_ratio = 1.2,dpi=600)


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
