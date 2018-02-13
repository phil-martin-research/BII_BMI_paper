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
library(cowplot)
library(reshape2)

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
bii_ecoregions<-extract(bii_corr,eco_regions,fun=mean,na.rm=TRUE)
#extract information on biomass values in each ecoregion
biomass_ecoregions<-extract(biomass_corr,eco_regions,fun=mean,na.rm=TRUE)

#remove nas
eco_regions2@data$hotspot_area2<-ifelse(is.na(eco_regions2@data$hotspot_area),0,eco_regions2@data$hotspot_area)
#put all this in a dataframe
hotspot_overlap_df<-data.frame(ECO_NAME=eco_regions2@data$ECO_NAME,hotspot_area=eco_regions2@data$hotspot_area2,bii=bii_ecoregions,biomass=biomass_ecoregions)
hotspot_overlap_df$hotspot<-ifelse(hotspot_overlap_df$hotspot_area>=50,"hotspot","not hotspot")


write.csv(hotspot_overlap_df,"data/hotspot_overlap.csv")
#plot some exploratory plots


median(hotspot_overlap_df$bii,na.rm=T)
median(hotspot_overlap_df$biomass,na.rm=T)


#plot data on bii vs biomass
hot_plot1<-ggplot(hotspot_overlap_df,aes(x=biomass,y=bii,colour=hotspot))+geom_hline(yintercept = 0.885,lty=2)+geom_vline(xintercept = 0.508,lty=2)+geom_point(shape=1)
hot_plot2<-hot_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="grey"))+ylab("biodiversity intactness index")+xlab("biomass intactness")

hot_plot1_v2<-ggplot(hotspot_overlap_df,aes(x=biomass,y=bii,colour=hotspot))+geom_point(shape=1)+geom_abline()+coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = T)
hot_plot2_v2<-hot_plot1_v2+scale_color_manual("",values=c("hotspot"="red","not hotspot"="grey"))+ylab("biodiversity intactness index")+xlab("biomass intactness")


save_plot(hot_plot2,filename = "figures/ecoregion_scatter.png",base_aspect_ratio = 1.5)


#plot the same kind of data as a violin plot
head(hotspot_overlap_df)
hotspot_melt<-melt(hotspot_overlap_df,id.vars = .("ECO_NAME","hotspot_area","hotspot"))
summary(hotspot_melt)
hot_violin1<-ggplot(hotspot_melt,aes(x=variable,y=value,fill=hotspot),colour="black")+geom_violin(draw_quantiles = c(0.5),trim = T)
hot_violin2<-hot_violin1+scale_fill_manual("",values=c("hotspot"="red","not hotspot"="grey"))+ylab("index value")
save_plot(hot_violin2,filename = "figures/ecoregion_violin.png",base_aspect_ratio = 1.5)

#boxplot
hot_box1<-ggplot(hotspot_melt,aes(x=variable,y=value,fill=hotspot),colour="black")+geom_boxplot(notch = TRUE,position = position_dodge(1))
hot_box2<-hot_box1+scale_fill_manual("",values=c("hotspot"="red","not hotspot"="grey"))+ylab("index value")
save_plot(hot_box2,filename = "figures/ecoregion_boxplot.png",base_aspect_ratio = 1.5)

#the same but as a dotplot
summary(hotspot_melt)
hot_dot1<-ggplot(hotspot_melt,aes(x=variable,y=value,colour=hotspot,fill=hotspot,group=interaction(variable, hotspot)))+geom_dotplot(shape=1,alpha=0.5,dotsize = 0.5,stackdir = "center",binaxis = "y",position="dodge")
hot_dot2<-hot_dot1+scale_fill_manual("",values=c("hotspot"="red","not hotspot"="grey"))+ylab("index value")+scale_colour_manual("",values=c("hotspot"="red","not hotspot"="grey"))

save_plot(hot_dot2,filename = "figures/ecoregion_dotplot.png",base_aspect_ratio = 4)
