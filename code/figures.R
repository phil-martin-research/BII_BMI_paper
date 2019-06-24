hotspot_overlap_df<-read.csv("data/hotspot_overlap.csv")
head(ecoregion_analysis)

#first calculate the median and interquartile ranges
hotspot_median<-ddply(hotspot_overlap_df,.(hotspot),summarize,m_bii=median(bii,na.rm=T),m_biomass=median(biomass,na.rm=T),
                      biom_25=quantile(biomass,probs=0.25,na.rm=T),biom_75=quantile(biomass,probs=0.75,na.rm=T),
                      bii_25=quantile(bii,probs=0.25,na.rm=T),bii_75=quantile(bii,probs=0.75,na.rm=T))


hot_plot1_v2<-ggplot(hotspot_overlap_df,aes(x=biomass,y=bii,colour=hotspot))+geom_point(shape=16,alpha=0.8)+geom_abline(lty=2)+coord_cartesian(xlim=c(0,1.1),ylim=c(0,1.1),expand = F)
hot_plot1_v2_pres<-hot_plot1_v2+theme(text=element_text(size=20))+xlab("biomass intactness")+ylab("biodiversity intactness")
save_plot(hot_plot1_v2_pres,filename = "figures/ecoregion_scatter_line_presentation.png",base_aspect_ratio = 1.6,dpi=600)

hot_plot2_v2<-hot_plot1_v2+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("Biodiversity Intactness Index")+xlab("biomass intactness")
hot_plot3_v2<-hot_plot2_v2+geom_point(data=hotspot_median,aes(x=m_biomass,y=m_bii,colour=hotspot),shape=15,size=2,stroke=1)+ theme(legend.position="none")
hot_plot4_v2<-hot_plot2_v2+geom_errorbarh(data=hotspot_median,aes(x=m_biomass,y=m_bii,xmin=biom_25,xmax=biom_75,colour=hotspot),size=0.5,height=0.005)+
  geom_errorbar(data=hotspot_median,aes(x=m_biomass,y=m_bii,ymin=bii_25,ymax=bii_75,colour=hotspot),size=0.5,width=0.005)+ theme(legend.position="none")
save_plot(hot_plot4_v2,filename = "figures/ecoregion_scatter_line_2.png",base_aspect_ratio = 1.2,dpi=600)


#plot hf vs bii for ecoregions
hf_bii_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=bii))+geom_point(shape=16,alpha=0.8)
hf_bii_plot1_pres<-hf_bii_plot1+theme(text=element_text(size=20))+xlab("human footprint index")+ylab("biodiversity intactness")
save_plot(hf_bii_plot1_pres,filename = "figures/hf_bii_presentation.png",base_aspect_ratio = 1.2,dpi=600)

hf_bii_plot2<-hf_bii_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("Biodiversity Intactness Index")+xlab("Human Footprint")
hf_bii_plot3<-hf_bii_plot2+theme(legend.position="none")
save_plot(hf_bii_plot3,filename = "figures/hf_bii_scatter.png",base_aspect_ratio = 1.2,dpi=600)

#plot hf vs biomass for ecoregions
hf_bmi_plot1<-ggplot(hotspot_overlap_df,aes(x=hf,y=biomass))+geom_point(shape=16,alpha=0.8)+geom_smooth(se=F,method = "glm",size=2,method.args = list(family = "binomial"))
hf_bmi_plot1_pres<-hf_bmi_plot1+theme(text=element_text(size=20))+xlab("human footprint index")+ylab("biomass intactness")+xlim(0,32)
save_plot(hf_bmi_plot1_pres,filename = "figures/hf_bmi_presentation.png",base_aspect_ratio = 1.2,dpi=600)

hf_bmi_plot2<-hf_bmi_plot1+scale_color_manual("",values=c("hotspot"="red","not hotspot"="dark grey"))+ylab("biomass intactness")+xlab("Human Footprint")
hf_bmi_plot3<-hf_bmi_plot2+theme(legend.position="none")
save_plot(hf_bmi_plot3,filename = "figures/hf_bmi_scatter.png",base_aspect_ratio = 1.2,dpi=600)


bottom_row<-plot_grid(hf_bmi_plot3,hf_bii_plot3,hot_plot4_v2,labels=c("b","c","d"),align='h',ncol=3)
combined_plots<-plot_grid(biv_map3,bottom_row,labels=c('a',''),ncol=1,rel_heights=c(1.6,1))
save_plot(combined_plots,filename = "figures/combined_plot_new.png",base_aspect_ratio = 1.2,dpi=600,base_height = 10,base_width = 12)
