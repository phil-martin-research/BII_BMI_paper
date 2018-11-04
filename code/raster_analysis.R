#this script is to calculate statistics from rasters using a grid

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
library(reghelper)
library(MuMIn)
library(mgcv)

##########################
#load in raster datasets##
##########################

#load in coastline dataset
coast<-readOGR(dsn="data/coastline/ne_10m_coastline.shp",layer="ne_10m_coastline")
#fortify coastline shapefile for plotting in ggplot
coast@data$id<-rownames(coast@data)
coast.points<-fortify(coast,region="id")


#load in BII map data
BII_map<-raster("data/bii/lbii.asc")

#load in biomass map
biomass<-raster("data/biomass_intactness/biomass_corrected.tif")

#load in nightlights, pasture, croplands, and human population density data
pop<-raster("data/population_density/gpw_v4_population_density_rev10_2010_2pt5_min.tif")
lights<-raster("data/lights/F182013.v4c_web.stable_lights.avg_vis.tif")
pasture<-raster("data/agriculture/pasture.tif")
croplands<-raster("data/agriculture/cropland.tif")

#reproject BII and HFP data to the  projection used by biomass data
proj4string(BII_map) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#resample all rasters based on biomass raster - the one with the lowest resolution
BII_resample<-resample(BII_map,biomass)
hfp_resample<-resample(hfp,biomass)
pop_resample<-resample(pop,biomass)
lights_resample<-resample(lights,biomass)
pasture_resample<-resample(pasture,biomass)
croplands_resample<-resample(croplands,biomass)
agriculture_resample<-resample(agriculture,biomass)


#get lat long coordinates for all rasters

#stack rasters
var_stack<-stack(BII_resample,biomass,pop_resample,lights_resample,pasture_resample,croplands_resample,agriculture_resample)

#convert raster to df
stack_df<- data.frame(na.omit(values(var_stack)))
stack_df$bii_trans<-stack_df$bii/max(stack_df$bii)
names(stack_df)<-c("bii","bmi","pop","light","pasture","cropland","agriculture","bii_trans")
write.csv(stack_df,file = "data/bii_bmi_stack.csv")

stack_df<-read.csv(file = "data/bii_bmi_stack.csv")

m1<- glm(bii_trans~log(pop+1)+log(light+1)+pasture+cropland, data=stack_df,family = 'binomial',na.action = "na.fail")
m2<- glm(bmi~log(pop+1)+log(light+1)+pasture+cropland, data=stack_df,family = 'binomial',na.action = "na.fail")

BII_dredge<-dredge(global.model = m1,evaluate = T,rank = "AICc",trace = T)
BMI_dredge<-dredge(global.model = m2,evaluate = T,rank = "AICc",trace = T)


bii_sc<-data.frame(std.coef(m1,partial.sd=T))
bii_sc$metric<-"BII"
bii_sc$variable<-c("intercept","log human\npopulation\ndensity","log nightlight","pasture","cropland")
bmi_sc<-data.frame(std.coef(m2,partial.sd=T))
bmi_sc$metric<-"BMI"
bmi_sc$variable<-c("intercept","log human\npopulation\ndensity","log nightlight","pasture","cropland")

sc_comb<-rbind(bii_sc,bmi_sc)
sc_comb_sub<-subset(sc_comb,variable!="intercept")

plot1<-ggplot(sc_comb_sub,aes(x=Estimate.,xmin=Estimate.-(1.96*Std..Error.),xmax=Estimate.+(1.96*Std..Error.),y=variable,colour=metric))+geom_point()+geom_errorbarh(height=0.1)+geom_vline(xintercept = 0,lty=2)
plot1+xlab("standardised coefficient")
ggsave(filename = "figures/standardised_coeff.png")
