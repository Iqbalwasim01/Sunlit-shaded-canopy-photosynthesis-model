
########### Input files and libraries
library(ggplot2)
library(foreach)

attach(Final.daily.data)

#### 2005
data2005<-Final.daily.data[0:152,]
NEEmax<-max(-data2005$NEE,na.rm = T)
LAI.NEE<-4.21*((-data2005$NEE)/NEEmax)
LAI.Xin<-LAI.seasonal.Xin.method(LAIs.guess = LAI.NEE,GPP.reported = "Yes") # make a new function for Xin for tower derived using moving average
No.kinetic.sets<-1  # enter No.sets of kinetics you intent to test in one go (or run of script) e.g. 1 if testing just one set
Kc<-42    # enter Kc field from dataframe or enter manually here
Kcat<-4.1  # enter Kcat field from dataframe or enter manually here
Sc.o<-87.3  # enter Sc.o field from dataframe or enter manually here
VmaxHA<-31
KcHA<-44.5
GamHA<-24.3
Ca<-360 # if reported specify field here, or give ambient atmospheric CO2 i.e. 360
LMA<-0.097 # chosen plant function type leaf mass area for C4 species
Sensor.height<-3.05 # example: change to match site flux tower height (Reported in site detail files)

##### Validation model LAI

### 2005 Growing season
GS.2005<-Two.big.leaf.concept(File.length = nrow(data2005),No.kinetic.sets=No.kinetic.sets,LAI=LAI.Xin,crop.type="Maize",leaf.angle="Corn",pathway="C4",Tair = data2005$Tair,wind = data2005$wind,humidity = data2005$RH/100,Precip = data2005$Precip
                             ,SWC_1 = data2005$SWC_2,Rn=data2005$Rn,H=data2005$H,LE=data2005$LE,LMA=LMA,pressure = data2005$pressure,SZA = data2005$SZA,SW_down = data2005$SW_down
                             ,Ca=Ca,Sensor.height = Sensor.height,Kcat =Kcat,Kc =Kc,Sc.o = Sc.o,KcHA=KcHA,VmaxHA=VmaxHA,GamHA=GamHA)

Canopy.data.frame.MODIS<-as.data.frame(unlist(GS.2005))

ggplot()+geom_line(data=Canopy.data.frame.MODIS,aes(x=0:151,y=`unlist(GS.2005)`,color="Modelled"),size=1.25)+ylim(0,40)+geom_line(data=data2005,aes(x=0:151,y=-NEE,color="Measured"),size=1.25)+ylab(expression("Net CO"[2]~"assimilation"~mu*"mol"~"m"^-2~"s"^-1))+xlab("")+scale_color_manual(values = c(
  'Modelled'='black',
  'Measured'='darkgray'))+labs(color="")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size=22),panel.background = element_rect(fill = 'white', colour = 'black'))

Stats.RMSE.MAE.R2(Predict=GS.2005,True=-data2005$NEE,No.obs = nrow(data2005))

### 2007 Growing season
data2007<-Final.daily.data[153:285,]
NEEmax<-max(-data2007$NEE,na.rm = T)
LAI.NEE<-4.21*((-data2007$NEE)/NEEmax)
LAI.Xin<-LAI.seasonal.Xin.method(LAIs.guess = LAI.NEE,GPP.reported = "Yes")

GS.2007<-Two.big.leaf.concept(File.length = nrow(data2007),No.kinetic.sets=No.kinetic.sets,LAI=LAI.Xin,crop.type="Maize",leaf.angle="Corn",pathway="C4",Tair = data2007$Tair,wind = data2007$wind,humidity = data2007$RH/100,Precip = data2007$Precip
                              ,SWC_1 = data2007$SWC_2,Rn=data2007$Rn,H=data2007$H,LE=data2007$LE,LMA=LMA,pressure = data2007$pressure,SZA = data2007$SZA,SW_down = data2007$SW_down
                              ,Ca=Ca,Sensor.height = Sensor.height,Kcat =Kcat,Kc =Kc,Sc.o = Sc.o,KcHA=KcHA,VmaxHA=VmaxHA,GamHA=GamHA)



Canopy.data.frame.MODIS<-as.data.frame(unlist(GS.2007))

ggplot()+geom_line(data=Canopy.data.frame.MODIS,aes(x=0:132,y=`unlist(GS.2007)`,color="Modelled"),size=1.25)+ylim(0,40)+geom_line(data=data2007,aes(x=0:132,y=-NEE,color="Measured"),size=1.25)+ylab(expression("Net CO"[2]~"assimilation"~mu*"mol"~"m"^-2~"s"^-1))+xlab("")+scale_color_manual(values = c(
  'Modelled'='black',
  'Measured'='darkgray'))+labs(color="")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),text = element_text(size=22),panel.background = element_rect(fill = 'white', colour = 'black'))

Stats.RMSE.MAE.R2(Predict=GS.2007,True=-data2007$NEE,No.obs = nrow(data2007))
