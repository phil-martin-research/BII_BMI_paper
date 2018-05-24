#script to produce scatterplot figures of BII and BMI vs various aspects of HF

library(ggplot2)
library(lme4)
library(MuMIn)
library(rsq)

#load data
hotspot_overlap_df<-read.csv("data/hotspot_overlap.csv")


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
M0<-glm(bii~1,data=hotspot_overlap_df,family = 'binomial',weights=ecoregion_area)
M1<-glm(bii~hf,data=hotspot_overlap_df,family = 'binomial',weights=ecoregion_area)
AICc(M0,M1)
summary(M1)
rsq(M1)
bii_hf_preds<-data.frame(hf=seq(0,40,0.011))
bii_hf_preds$preds<-plogis(predict(M1,newdata=bii_hf_preds))


#plot hf vs bii for ecoregions
hf_bii_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=bii,colour=hotspot,size=ecoregion_area))+geom_point(shape=1,alpha=0.8)
hf_bii_plot2<-hf_bii_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biodiversity intactness index")+xlab("human footprint")
hf_bii_plot3<-hf_bii_plot2+theme(legend.position="none")+geom_line(data=bii_hf_preds,aes(x=hf,y=preds),size=1,colour="black")
save_plot(hf_bii_plot3,filename = "figures/hf_bii_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#model relationship between hf and bmi
M0<-glm(biomass~1,data=hotspot_overlap_df,family = 'binomial',weights=ecoregion_area)
M1<-glm(biomass~hf,data=hotspot_overlap_df,family = 'binomial',weights=ecoregion_area)
AICc(M0,M1)
summary(M1)
rsq(M1)
bmi_hf_preds<-data.frame(hf=seq(0,40,0.001))
bmi_hf_preds$preds<-plogis(predict(M1,newdata=bmi_hf_preds))

#plot hf vs biomass for ecoregions
hf_bmi_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=biomass,colour=hotspot,size=ecoregion_area))+geom_point(shape=1,alpha=0.8)
hf_bmi_plot2<-hf_bmi_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biomass intactness")+xlab("human footprint")
hf_bmi_plot3<-hf_bmi_plot2+theme(legend.position="none")+geom_line(data=bmi_hf_preds,aes(x=hf,y=preds),size=1,colour="black")
save_plot(hf_bmi_plot3,filename = "figures/hf_bmi_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#exploratory plots for elements of hf against bmi and bii
head(hotspot_overlap_df)
ggplot(data=hotspot_overlap_df,aes(x=pop_den,y=bii,size=ecoregion_area))+geom_point()
ggplot(data=hotspot_overlap_df,aes(x=pop_den,y=biomass,size=ecoregion_area))+geom_point()
ggplot(data=hotspot_overlap_df,aes(x=light,y=bii,size=ecoregion_area))+geom_point()
ggplot(data=hotspot_overlap_df,aes(x=light,y=biomass,size=ecoregion_area))+geom_point()
ggplot(data=hotspot_overlap_df,aes(x=pasture,y=bii,size=ecoregion_area))+geom_point()
ggplot(data=hotspot_overlap_df,aes(x=pasture,y=biomass,size=ecoregion_area))+geom_point()
ggplot(data=hotspot_overlap_df,aes(x=croplands,y=bii,size=ecoregion_area))+geom_point()
ggplot(data=hotspot_overlap_df,aes(x=croplands,y=biomass,size=ecoregion_area))+geom_point()

