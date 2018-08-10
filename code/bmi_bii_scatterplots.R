#script to produce scatterplot figures of BII and BMI vs various aspects of HF

rm(list = ls())

library(ggplot2)
library(lme4)
library(MuMIn)
library(rsq)
library(cowplot)
library(dplyr)
library(scales)


#load data
hotspot_overlap_df<-read.csv("data/hotspot_overlap.csv",stringsAsFactors = F)
hotspot_overlap_df$agriculture<-hotspot_overlap_df$pasture+hotspot_overlap_df$cropland
str(hotspot_overlap_df)

eco_unique<-unique(hotspot_overlap_df$biome)
eco_biome<-hotspot_overlap_df$biome

hotspot_overlap_df$forest<-case_when(
  eco_biome==eco_unique[2]~"Forest",
  eco_biome==eco_unique[3]~"Forest",
  eco_biome==eco_unique[5]~"Forest",
  eco_biome==eco_unique[6]~"Forest",
  eco_biome==eco_unique[7]~"Forest",
  eco_biome==eco_unique[8]~"Forest",
  eco_biome==eco_unique[10]~"Forest",
  eco_biome==eco_unique[13]~"Forest",
  eco_biome==eco_unique[14]~"Forest")
hotspot_overlap_df$forest<-case_when(is.na(hotspot_overlap_df$forest)~"Non-forest",
                                     hotspot_overlap_df$forest=="Forest"~"Forest")

ggplot(hotspot_overlap_df,aes(x=biomass,y=bii,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)
ggplot(hotspot_overlap_df,aes(x=hf,y=bii,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")
ggplot(hotspot_overlap_df,aes(x=hf,y=biomass,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")
ggplot(hotspot_overlap_df,aes(x=cropland,y=biomass,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")+geom_smooth(method="glm")
ggplot(hotspot_overlap_df,aes(x=cropland,y=bii,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")+geom_smooth(method="glm")
ggplot(hotspot_overlap_df,aes(x=pasture,y=biomass,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")+geom_smooth(method="glm")
ggplot(hotspot_overlap_df,aes(x=pasture,y=bii,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")
ggplot(hotspot_overlap_df,aes(x=pasture+croplands,y=biomass,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")+geom_smooth(method="glm")
ggplot(hotspot_overlap_df,aes(x=pasture+croplands,y=bii,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")+geom_smooth(method="glm")
ggplot(hotspot_overlap_df,aes(x=pop_den,y=biomass,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")+geom_smooth(method="glm")+scale_x_continuous(trans = "log10",labels=comma,breaks = c(0.01,1,100,1000))
ggplot(hotspot_overlap_df,aes(x=pop_den,y=bii,size=ecoregion_area,colour=forest))+geom_point(shape=1,alpha=0.8)+theme(legend.position="none")+geom_smooth(method="glm")+scale_x_continuous(trans = "log10",labels=comma,breaks = c(0.01,1,100,1000))



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


#model relationship between hf and bii
hotspot_overlap_df$bii_trans<-(hotspot_overlap_df$bii)/max(hotspot_overlap_df$bii,na.rm = T)
M0<-glm(bii_trans~1,data=hotspot_overlap_df,family = 'binomial')
M1<-glm(bii_trans~hf,data=hotspot_overlap_df,family = 'binomial')
AICc(M0,M1)
#null model is more likely



#plot hf vs bii for ecoregions
hf_bii_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=bii,size=ecoregion_area,colour=hotspot))+geom_point(shape=1,alpha=0.8)
hf_bii_plot2<-hf_bii_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biodiversity intactness index")+xlab("human footprint")
hf_bii_plot3<-hf_bii_plot2+theme(legend.position="none")
save_plot(hf_bii_plot3,filename = "figures/hf_bii_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#model relationship between hf and bmi
M0<-glm(biomass~1,data=hotspot_overlap_df,family = 'binomial')
M1<-glm(biomass~hf,data=hotspot_overlap_df,family = 'binomial')
AICc(M0,M1)
summary(M1)
rsq(M1)
bmi_hf_preds<-data.frame(hf=seq(0,40,0.001))
bmi_hf_preds$preds<-predict(M1,newdata=bmi_hf_preds,type="response")

#plot hf vs biomass for ecoregions
hf_bmi_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=biomass,colour=hotspot,size=ecoregion_area))+geom_point(shape=1,alpha=0.8)
hf_bmi_plot2<-hf_bmi_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biomass intactness")+xlab("human footprint")
hf_bmi_plot3<-hf_bmi_plot2+theme(legend.position="none")+geom_line(data=bmi_hf_preds,aes(x=hf,y=preds),size=1,colour="black")
save_plot(hf_bmi_plot3,filename = "figures/hf_bmi_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#exploratory plots for elements of hf against bmi and bii


#produce models for agriculture
AG_BII_M0<-glm
AG_BII_M1<-glm(bii_trans~agriculture,data=hotspot_overlap_df,family="binomial")
rsq(AG_BII_M1)
bmi_ag_preds<-data.frame(agriculture=seq(0,1,0.001))
bmi_ag_preds$preds<-predict(AG_BII_M1,newdata=bmi_ag_preds,type="response")


AICc(AG_BII_M0,AG_BII_M1)
#plot these results
ag_bii_plot1<-ggplot(hotspot_overlap_df,aes(x=agriculture*100,y=bii,colour=hotspot,size=ecoregion_area))+geom_point(shape=1,alpha=0.8)
ag_bii_plot2<-ag_bii_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biodiversity intactness index")+xlab("Percentage agricultural land")
ag_bii_plot3<-ag_bii_plot2+theme(legend.position="none")+geom_line(data=bmi_ag_preds,aes(x=agriculture*100,y=preds),size=1,colour="black")
save_plot(ag_bii_plot3,filename = "figures/agriculture_bii_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#models for biomass and agriculture
AG_BMI_M0<-glm(biomass~1,data=hotspot_overlap_df,family="binomial")
AG_BMI_M1<-glm(biomass~agriculture,data=hotspot_overlap_df,family="binomial")
AG_BMI_M2<-glm(biomass~agriculture*forest,data=hotspot_overlap_df,family="binomial")
rsq(AG_BMI_M2)
AICc(AG_BMI_M0,AG_BMI_M1,AG_BMI_M2)
bmi_ag_preds<-expand.grid(agriculture=seq(0,1,0.001),forest=c("Forest","Non-forest"))

bmi_ag_preds$preds<-predict(AG_BMI_M2,newdata=bmi_ag_preds,type="response")
#plot these results
ag_bmi_plot1<-ggplot(hotspot_overlap_df,aes(x=agriculture*100,y=biomass,colour=forest,size=ecoregion_area))+geom_point(shape=1,alpha=0.8)
ag_bmi_plot2<-ag_bmi_plot1+ylab("biomass intactness")+xlab("Percentage agricultural land")
ag_bmi_plot3<-ag_bmi_plot2+theme(legend.position="none")+geom_line(data=bmi_ag_preds,aes(x=agriculture*100,y=preds,colour=forest),size=1)+scale_color_manual("",values=c("Forest"="green4","Non-forest"="dark grey"))
save_plot(ag_bmi_plot3,filename = "figures/agriculture_bmi_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#now plot relationship between population density and abii and bmi


#produce models population vs bii
POP_BII_M0<-glm(bii_trans~1,data=hotspot_overlap_df,family="binomial")
POP_BII_M1<-glm(bii_trans~log(pop_den),data=hotspot_overlap_df,family="binomial")
rsq(POP_BII_M1)
AICc(POP_BII_M0,POP_BII_M1)
#null model is better

pop_bii_plot1<-ggplot(hotspot_overlap_df,aes(x=pop_den,y=bii,colour=forest,size=ecoregion_area))+geom_point(shape=1,alpha=0.8)
pop_bii_plot2<-pop_bii_plot1+ylab("biodiversity intactness index")+xlab("Human population density")
pop_bii_plot3<-pop_bii_plot2+theme(legend.position="none")+scale_x_continuous(trans="log",labels=comma,breaks=c(0.1,1,10,100,1000))+scale_color_manual("",values=c("Forest"="green4","Non-forest"="dark grey"))
save_plot(pop_bii_plot3,filename = "figures/pop_bii_scatter.png",base_aspect_ratio = 1.2,dpi=600)


#produce models population vs bmi
POP_BMI_M0<-glm(biomass~1,data=hotspot_overlap_df,family="binomial")
POP_BMI_M1<-glm(biomass~log(pop_den),data=hotspot_overlap_df,family="binomial")
POP_BMI_M2<-glm(biomass~log(pop_den)*forest,data=hotspot_overlap_df,family="binomial")
rsq(POP_BMI_M1)
AICc(POP_BMI_M0,POP_BMI_M1,POP_BMI_M2)
bmi_pop_preds<-expand.grid(pop_den=seq(0.001,2600,0.01),forest=c("Forest","Non-forest"))
bmi_pop_preds$preds<-predict(POP_BMI_M2,newdata=bmi_pop_preds,type="response")
#null model is better

pop_bmi_plot1<-ggplot(hotspot_overlap_df,aes(x=pop_den,y=biomass,colour=forest,size=ecoregion_area))+geom_point(shape=1,alpha=0.8)
pop_bmi_plot2<-pop_bmi_plot1+ylab("biodiversity intactness index")+xlab("Human population density")
pop_bmi_plot3<-pop_bmi_plot2+theme(legend.position="none")+geom_line(data=bmi_pop_preds,aes(x=pop_den,y=preds,colour=forest),size=1)+scale_x_continuous(trans="log",labels=comma,breaks=c(0.1,1,10,100,1000))+scale_color_manual("",values=c("Forest"="green4","Non-forest"="dark grey"))
save_plot(pop_bmi_plot3,filename = "figures/pop_bmi_scatter.png",base_aspect_ratio = 1.2,dpi=600)

