#####################################################################################
#this script is to process raster and vector data to  look at differences in  #######
#the Biodiversity Intactness Index and the biomass intactness index and compare them#
#to the human footprint index########################################################
#####################################################################################

.libPaths("C:/R/Library")
.Library<-"C:/R/Library"


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
library(rgeos)
library(devtools)
library(archive)
library(rdryad)
library(dplyr)
library(stringr)


beginCluster(n=6,type="SOCK") 

###################################################################
#1 - Processsing of spatial data###################################
###################################################################

BII_map<-raster("data/lbii.asc")#load in BII map data downloaded from https://doi.org/10.5519/0009936

#download biomass data from Erb et al 2017 https://doi.org/10.1038/nature25138
download.file("https://boku.ac.at/fileadmin/data/H03000/H73000/H73700/Data_Download/Data/NatureDataDownload.zip", 
              destfile = 'NatureDataDownload.zip')
unzip(zipfile = "data/NatureDataDownload.zip",#unzip biomass dataset 
      exdir = 'data')
unzip(zipfile = "data/NATURE_DATADOWNLOAD.zip",#unzip biomass dataset 
      exdir = 'data')
file.rename(from="data/NEW/C_reduction_perc.tif",to="data/C_reduction_perc.tif")#move file from subdirectory
unlink("data/NEW/",recursive = TRUE)#remove biomass folder with unwanted files
file.remove("data/NATURE_DATADOWNLOAD.zip")#remove biomass zips with unwanted files
file.remove("data/NatureDataDownload.zip")#remove biomass zips with unwanted files
biomass<-raster("data/C_reduction_perc.tif")#load in biomass map
getwd()
#download human footprint data from Venter et al 2015 https://doi.org/10.1038/sdata.2016.67 (data location http://dx.doi.org/10.5061/dryad.052q5)
tf <- tempfile()
td <- tempdir()
file.path <-"https://datadryad.org/bitstream/handle/10255/dryad.131447/HumanFootprintv2.7z?sequence=2"
download.file( file.path , tf , mode = "wb" )
hfp_files<-a%>% 
  select(path, size,date) %>% 
  filter(str_detect(string=tolower(path), pattern = "hfp2009.tif"))
hfp_files
hfp<-raster(hfp_files[1,]) # need to fix here


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
writeRaster(BII_proj_capped, "data/BII_capped", format = "GTiff",overwrite=T) #save BII map with max value of 1

#create a grid to extract data to with a resolution of 0.083333 degrees
r <- raster(extent(BII_proj_capped), nrow=1804, ncol=3607,
            crs = "+proj=moll +lon_0=0 +ellps=WGS84")
r[] <- 1:ncell(r)
sp.r <- as(r, "SpatialPixelsDataFrame")


sp.r$biomass<-(extract(biomass_inv,sp.r))#extract values from biomass raster
sp.r$bii<-(extract(BII_proj_capped,sp.r))#extract values from BII data
spr_df<-as.data.frame(sp.r)#convert grid to a dataframe for plotting in  ggplot

write.csv(spr_df,"data/output_data//BII_BMI_data.csv")#save this data as a csv file
summary(spr_df)

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
spr_df_cc<-spr_df[complete.cases(spr_df),]#remove rows with missing data

#download coastline dataset
download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", 
              destfile = 'coastlines.zip')
unzip(zipfile = "coastlines.zip",#unzip coastline dataset 
      exdir = 'data')
coast<-readOGR("data/ne_10m_coastline.shp")
coast_moll<- spTransform(coast, crs(biomass_inv))
coast_simp<-gSimplify(coast_moll,tol = 0.1,topologyPreserve = TRUE)#simplify coastline to speed up plotting
coast_df<-SpatialLinesDataFrame(coast_simp,coast@data) 


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
         theme_bw()+coord_equal(ylim = c(-6500000,8500000))+
         guides(fill = guide_colorplane(title.theme = element_text(size = 13),label.theme = theme_gray(),
         label.y.theme = theme_gray(),axis.ticks = element_blank()))+
         theme(axis.line=element_blank(),panel.background=element_blank(),panel.border=element_blank(),
         panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank(),
         axis.ticks = element_blank(),
         legend.position=c(.1,.3),legend.key.size=unit(0.8,"cm"))
biv_map3<-biv_map2+geom_path(data=coast_df,aes(x = long, y = lat, group = group),colour="black",size=0.1)
biv_map3


hotspot_overlap_df<-read.csv("data/output_data/ecoregion_hotspot_stats.csv")#load ecoregion statistics
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

#plot hf vs bii for ecoregions
hf_bii_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=bii,colour=hotspot))+geom_point(shape=1,alpha=0.8)
hf_bii_plot2<-hf_bii_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("Biodiversity Intactness Index")+xlab("Human Footprint")
hf_bii_plot3<-hf_bii_plot2+theme(legend.position="none")

#plot hf vs biomass for ecoregions
hf_bmi_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=biomass,colour=hotspot))+geom_point(shape=1,alpha=0.8)
hf_bmi_plot2<-hf_bmi_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biomass intactness")+xlab("Human Footprint")
hf_bmi_plot3<-hf_bmi_plot2+theme(legend.position="none")

bottom_row<-plot_grid(hf_bmi_plot3,hf_bii_plot3,hs_BMI_BII_4,labels=c("b","c","d"),align='h',ncol=3)
combined_plots<-plot_grid(biv_map3,bottom_row,labels=c('a',''),ncol=1,rel_heights=c(1.6,1))
save_plot(combined_plots,filename = "figures/combined_plot_new.png",base_aspect_ratio = 1.2,dpi=300,base_height = 10,base_width = 12)
#this plot was then modified in GIMP to change contrast for the map and remove axis 
#labels that I couldn't remove programatically
