
########### Input files and libraries
library(ggplot2)
library(foreach)

attach(Final.daily.data)

NEEmax<-max(-Final.daily.data$NEE,na.rm = T)
LAI.NEE<-4*((-Final.daily.data$NEE)/NEEmax)
LAI.Xin<-LAI.seasonal.Xin.method(LAIs.guess = LAI.NEE,GPP.reported = "Yes") # make a new function for Xin for tower derived using moving average
No.kinetic.sets<-1  # enter No.sets of kinetics you intent to test in one go (or run of script) e.g. 1 if testing just one set
Kc<-18.6    # enter Kc field from dataframe or enter manually here
Kcat<-2  # enter Kcat field from dataframe or enter manually here
Sc.o<-101  # enter Sc.o field from dataframe or enter manually here
VmaxHA<-51.2
KcHA<-57
GamHA<-19.8
Ca<-360 # if reported specify field here, or give ambient atmospheric CO2 i.e. 360
LMA<-0.087 # chosen plant function type leaf mass area for C3 species
Sensor.height<-2.7 # example: change to match site flux tower height (Reported in site detail files)


##### Validation model

Canopy<-Two.big.leaf.concept(File.length = nrow(Final.daily.data),No.kinetic.sets=No.kinetic.sets,LAI=LAI.Xin,crop.type="C3",leaf.angle="Spherical",pathway="C3",Tair = Tair,wind = wind,humidity = RH/100,Precip = Precip
                             ,SWC_1 = SWC_1,Rn=Rn,H=H,LE=LE,LMA=LMA,pressure = pressure,SZA = SZA.zenith,SW_down = SW_down
                             ,Ca=Ca,Sensor.height = Sensor.height,Kcat =Kcat,Kc =Kc,Sc.o = Sc.o,KcHA=KcHA,VmaxHA=VmaxHA,GamHA=GamHA)

Canopy.data.frame<-as.data.frame(unlist(Canopy))

ggplot()+geom_line(data=Canopy.data.frame,aes(x=0:152,y=Canopy.data.frame$`unlist(Canopy)`,color="Modelled"),size=1.25)+geom_line(data=Final.daily.data,aes(x=0:152,y=-Final.daily.data$NEE,color="Measured"),size=1.25)+ylim(0,20)+ylab(expression("Net CO"[2]~"assimilation"~mu*"mol"~"m"^-2~"s"^-1))+xlab("")+scale_color_manual(values = c(
  'Modelled' = 'black',
  'Measured' = 'darkgray'))+labs(color="")+theme(axis.title.x=element_blank(),axis.title=element_text(size=22,face="bold"),axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size=22),panel.background = element_rect(fill = 'white', colour = 'black'))

Stats.RMSE.MAE.R2(Predict=Canopy,True=-Final.daily.data$NEE,No.obs = nrow(Final.daily.data))
