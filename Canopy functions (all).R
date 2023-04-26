
Sunlit.shaded.parameter<-function(Param=17,Kb=0.7,LAI=3,Kn=0.3) { 
  # Kb calculated from two stream approximation of canopy radiation
  Kb=Kb
  Param=Param
   # fraction of green leaf area index
  
  # Calculate parameter for sunlit region and wuthout any sunlit scaling  
  # Scale leaf-level Vcmax to canopy level using LAI   
  Top.canopy<-Param
  
  
  Paramsunlit=Top.canopy*((1-exp(-(Kn+Kb)*LAI))*(1/(Kn+Kb)))
 
  # Calculate parameter for shaded region

  Paramshaded=Top.canopy*(((1-exp(-Kn*LAI))*(1/Kn))-((1-exp(-(Kn+Kb)*LAI))*(1/(Kn+Kb))))
  
  Paramsunlit<-as.data.frame(Paramsunlit)
  Paramshaded<-as.data.frame(Paramshaded)
  Data<-cbind(Paramsunlit,Paramshaded)
  
  return(Data)
  }
  
  
Ci.quadratic<-function(Jmax25=93,Vcmax25=Vcmax25.1,gs=GS,Ca=Ca) {
  
  K.c=404.9
  K.o = 278.4
  gamma=42.75
  O.i=0.21
  Rd<-0.1
  Jb=2*gamma
  Vb=K.c*(1 + O.i/K.o)
  
  # ax^2+bx+c=0
  V1=gs
  V2=(Vcmax25-Rd+(Vb-Ca)*gs)
  V3=-(Vcmax25*gamma+(Ca*gs+Rd)*Vb)
  
  J1=gs
  J2=((Jmax25/4)-Rd+(Jb-Ca)*gs)
  J3=-((Jmax25/4)*gamma+(Ca*gs+Rd)*Jb)
 
  Cic<-(-V2+sqrt(V2^2-4*V1*V3))/(2*V1)
  Cij<-(-J2+sqrt(J2^2-4*J1*J3))/(2*J1)
  
  Ci<-min(Cic,Cij)
  
  
  return(Ci) }

Esat<-function(Temp=285,pressure=100){
  
  if(is.na(pressure)){
    pressure=100
  }else{
    pressure=pressure
  }
  
  T.final<-Temp-273.15
  pressure1<-pressure*1000

  Esat.1<-ifelse(T.final<0,100*(6.11123516+T.final*(5.03109514e-01+T.final*(1.88369801e-02+T.final*(4.20547422e-04+T.final*(6.14396778e-06+T.final*(6.02780717e-08+T.final*(3.87940929e-10+T.final*(1.49436277e-12+T.final*2.62655803e-15)))))))),
                    100*(6.11213476+T.final*(4.44007856e-01+T.final*(1.43064234e-02+T.final*(2.64461437e-04+T.final*(3.05903558e-06+T.final*(1.96237241e-08+T.final*(8.92344772e-11+T.final*(-3.73208410e-13+T.final*2.09339997e-16)))))))))  
  
  Esat.derivative<-ifelse(T.final<0,100*(5.03277922e-01+T.final*(3.77289173e-02+T.final*(1.26801703e-03+T.final*(2.49468427e-05+T.final*(3.13703411e-07+T.final*(2.57180651e-09+T.final*(1.33268878e-11+T.final*(3.94116744e-14+T.final*4.98070196e-17)))))))),
                          100*(4.44017302e-01+T.final*(2.86064092e-02+T.final*(7.94683137e-04+T.final*(1.21211669e-05+T.final*(1.03354611e-07+T.final*(4.04125005e-10+T.final*(-7.88037859e-13+T.final*(-1.14596802e-14+T.final*3.81294516e-17)))))))))  
  
  qsat.derivative.air<-Esat.derivative/pressure1
  qsat.air<-Esat.1/pressure1
  qsat.derivative.air<-as.data.frame(qsat.derivative.air)
  qsat.air<-as.data.frame(qsat.air)
  Esat.pascals<-as.data.frame(Esat.1)
  colnames(Esat.pascals)<-"Esat(Pascals)"
  colnames(qsat.derivative.air)<-"qsat.derivative"
  colnames(qsat.air)<-"qsat"
  
  Data<-cbind(qsat.derivative.air,qsat.air,Esat.pascals)
  
  return(Data)
}



Stomatal.quadratic.equation<-function(Stress=0.9,Gb=2,An=0.1,g0=0.01,g1=1.79,Cs=360,ea=450,esat=500){
  
  Stress=Stress
  Gb=Gb
  An=An
  g0=g0
  g1=g1
  Cs=Cs
  ea=ea
  esat=esat
  
  a<-1
  b<-Gb-g0-g1*An/Cs
  c<--Gb*(g0+((g1*An/Cs)*(ea/esat)))  
  
  Gs1<-(-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  Gs2<-(-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  Gs<-pmax(Gs1,Gs2)
  return(Gs)
}


Ci.iteration<-function(max.iterations=20,Ci.old=241,fxCi=118.72,Ci.intermediate=3.933,Ca=360,ea=0.01264362790,esat=0.01456363309,Gb=2.416,g1=1.79,Gs.stress=9,Stress=0.7
                       ,Vcmax=78,PPAR=500,Kc25=269,pathway="c3",Sc.o = 167,Temp=273,VmaxHA=20,GamHA=40,KcHA=50){
  count<-0
  max.iterations<-max.iterations
  Stress<-Stress
  Ci.old<-Ci.old
  Ci.bisec.guess<-Ci.old
  FxCi.old<-fxCi
  Ci.intermediate.old<-Ci.intermediate
  Gb<-Gb
  Gs.stress<-Gs.stress
  Ca<-Ca
  ea<-ea
  esat<-esat
  Vcmax=Vcmax
  PPAR=PPAR
  Kc25=Kc25
  Sc.o = Sc.o
  Temp=Temp
  VmaxHA=VmaxHA
  GamHA=GamHA
  KcHA=KcHA
  pathway=pathway

  repeat{
    
    Ci.new<-Ci.old-(FxCi.old/((Ci.intermediate.old-FxCi.old)/FxCi.old))
    
  
    An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = Ci.new,
                           Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
    Cs<-Ca-(1.4/Gb)*An
    Gs<-Stomatal.quadratic.equation(Stress,Gb,An,0.01,Gs.stress,Cs,ea,esat)
    FxCi.old<-(Cs-(1.65/Gs)*An)-Ci.new # f(x)
    
   
    An.i<-Farquhar.version.3(Kc25 = Kc25,Sc.o = Sc.o,c.i = FxCi.old+Ci.new,
                             Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate f'(x)
    Cs.i<-Ca-(1.4/Gb)*An.i
    Gs.i<-Stomatal.quadratic.equation(Stress,Gb,An.i,0.01,Gs.stress,Cs.i,ea,esat)
    Ci.i<-(Cs.i-(1.65/Gs.i)*An.i)-(FxCi.old+Ci.new)
    Ci.intermediate.old<-Ci.i
   
    Conv<-abs(Ci.new-Ci.old)
    ifelse(Conv<0.0001,break,Ci.new)
    Ci.old<-Ci.new
    
    count<-count+1
    
    ifelse(count>max.iterations,break,NA) 
  }
      ### Bisection backup method
      Conv<-Conv # break next loop if first loop is successful
      count<-0 # Reset the count
      a<-5 # min interval guess
      An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = a,
                             Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
      Cs<-Ca-(1.4/Gb)*An
      Gs<-Stomatal.quadratic.equation(Stress,Gb,An,0.01,Gs.stress,Cs,ea,esat)
      Fx.a<-(Cs-(1.65/Gs)*An)-a
      
      b<-400 # max interval guess
      
      repeat{
        ifelse(Conv<0.0001,break,NA)
        c=(a+b)/2
        
        An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = c,
                               Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
        Cs<-Ca-(1.4/Gb)*An
        Gs<-Stomatal.quadratic.equation(Stress,Gb,An,0.01,Gs.stress,Cs,ea,esat)
        Fx.c<-(Cs-(1.65/Gs)*An)-c
        
        ifelse(abs(Fx.c)<0.0001,break,c)

          a=ifelse((Fx.a*Fx.c)<0,a,c)
          b=ifelse((Fx.a*Fx.c)<0,c,400)
          Fx.a=ifelse((Fx.a*Fx.c)>0,Fx.c,Fx.a)
        
        Ci.new<-a
        count<-count+1
        ifelse(count==max.iterations,break,NA)
        
      }
      return(Ci.new)
      options(digits = 20)
    }
 


Air.density<-function(Tair=285,pressure=100){
  g=8.314
  pressure=pressure*1000
  
  Pm<-pressure/(Tair*g)
  
  return(Pm)
}



Farquhar.version.3<-function(pathway="C3",Kc25=10.8,Sc.o=99.9,c.i=300,Vcmax25=20,KcHA=63.4
                             ,VmaxHA=41.5,GamHA=23.6,PPAR=700,Temp=298.15){
  Jmax25=Vcmax25*1.65
  Rd25=0.015*Vcmax25
  Kc25=Kc25/0.0334
  O.i=210
  Gam25=(0.5*(O.i/Sc.o))/0.0334
  Kc_Ha=KcHA
  Gam_Ha=GamHA
  Vmax_Ha=VmaxHA
  S_Vmax=0.66839
  S_Jmax=0.65970
  S_Rd=0.490
  Vmax_Hd=200
  Rd_Ha=46.390
  Rd_Hd=150.650
  Jmax_Ha=50
  Jmax_Hd=200

  ### C3 temperature adjustments
  K.c<-Kc25 * exp(Kc_Ha/298.15*8.314*(1-(298.15/Temp))) 
  
  gamma<-Gam25 * exp(Gam_Ha/298.15*8.314*(1-(298.15/Temp)))
  
  Rd<-ifelse(Temp>298.15,Rd25 * exp(Rd_Ha/298.15*8.314*(1-(298.15/Temp))) * (1+exp((298.15*S_Rd-Rd_Ha)/(298.15*0.001*8.314)))/(1+exp((S_Rd*Temp-Rd_Hd)/(0.001*8.314*Temp))),Rd25 * exp(Rd_Ha/298.15*8.314*(1-(298.15/Temp))))
  
  Vcmax<-ifelse(Temp>298.15,Vcmax25 * exp(Vmax_Ha/298.15*8.314*(1-(298.15/Temp))) * (1+exp((298.15*S_Vmax-Vmax_Hd)/(298.15*0.001*8.314)))/(1+exp((S_Vmax*Temp-Vmax_Hd)/(0.001*8.314*Temp))),Vcmax25 * exp(Vmax_Ha/298.15*8.314*(1-(298.15/Temp))))
    
  Tp<-0.167*Vcmax # CLM 5.0 'Tp' uses same temperature adjustment parameters as Vcmax
  
  Jmax<-ifelse(Temp>298.15,Jmax25 * exp(Jmax_Ha/298.15*8.314*(1-(298.15/Temp))) * (1+exp((298.15*S_Jmax-Jmax_Hd)/(298.15*0.001*8.314)))/(1+exp((S_Jmax*Temp-Jmax_Hd)/(0.001*8.314*Temp))),Jmax25 * exp(Jmax_Ha/298.15*8.314*(1-(298.15/Temp)))) 

  
  # Jmax adjusted for light conditions
  IPSII<- 0.5*0.85*PPAR
  a <- 0.7
  b <- -(IPSII+Jmax)
  c <-IPSII*Jmax
  J.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  J.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  Jmax <-min(J.1, J.2)
  
  if(pathway=="C3"){

  Ac <- Vcmax * (c.i - gamma) / (c.i + K.c)
  
  Aj <- Jmax * (c.i - gamma) / (4 * (c.i + 2 * gamma))
  
  Ap <- 3*Tp
  
  a=0.98
  b=-(Ac+Aj)
  c=Ac*Aj
  A.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  Ai <-min(A.1, A.2)
  
  a=0.95
  b=-(Ai+Ap)
  c=Ai*Ap
  A.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A <-min(A.1, A.2)
  

  A.n<-A-Rd
  
  
  }else{
    ### C4 temperature adjustments
 

  Vcmax=Vcmax25*(exp(Vmax_Ha/298.15*8.314*(1-(298.15/Temp)))/((1+exp(0.3*(Temp-313.15)))*(1+exp(0.2*(288.15-Temp)))))
   
  Rd=Vcmax*0.025
   
  Kp=20000*Vcmax # Houborg et al 2013

  Ap<-Kp*c.i
  Aj<-0.05*PPAR
  Ac<-Vcmax
  
  a=0.80
  b=-(Ac+Aj)
  c=Ac*Aj
  A.i.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A.i.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  Ai <-min(A.i.1, A.i.2) # solution given by the smaller root
  
  a=0.95
  b=-(Ai+Ap)
  c=Ai*Ap
  A.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A <-min(A.1, A.2)
  
  A.n<-A-Rd
  
}  
  
  return(A.n)
  
}

  

Soil.stress.parameter<-function(Soil.moist=0.5,soil.type="Silty clay loam"){
  
  if(soil.type=="Silty clay loam"){
    Soil.moist.sat<-0.477 
  }else if(soil.type=="Sand"){
    Soil.moist.sat<-0.395
  }else if(soil.type=="Loamy sand"){
    Soil.moist.sat<-0.410
  }else if(soil.type=="Sandy loam"){
    Soil.moist.sat<-0.435
  }else if(soil.type=="Silt loam"){
    Soil.moist.sat<-0.485
  }else if(soil.type=="Loam"){
    Soil.moist.sat<-0.451
  }else if(soil.type=="Sandy clay loam"){
    Soil.moist.sat<-0.420
  }else if(soil.type=="Clay loam"){
    Soil.moist.sat<-0.476
  }else if(soil.type=="Sandy clay"){
    Soil.moist.sat<-0.426
  }else if(soil.type=="Silty clay"){
    Soil.moist.sat<-0.492
  }else{
    Soil.moist.sat<-0.482
  } # Clay

  MinBw<-0.01
    
  if(is.na(Soil.moist)){
    Bw<-1
  }else{
    Bw<-ifelse(Soil.moist>Soil.moist.sat,1,ifelse(Soil.moist>0.01 & Soil.moist<0.477,((Soil.moist-MinBw)/(Soil.moist.sat-MinBw)),0))
    
  }

  return(Bw)
}


  
  
  Two.big.leaf.concept<-function(File.length=1,No.kinetic.sets=1,LAI=4,crop.type="C3",leaf.angle="Spherical",pathway="C3",Tair=25,wind=2,humidity=90/100,Precip=0.5
                                 ,SWC_1=30,Rn=NA,H=NA,LE=NA,LMA=2.2,pressure=100,SZA=30,SW_down=500
                                 ,Ca=360,Sensor.height=2.7,Kcat=2,Kc=10,Sc.o=101,KcHA=20,VmaxHA=30,GamHA=40){
              Count<-0

          An<-foreach(Kcat=Kcat,Kc=Kc,Sc.o=Sc.o,KcHA=KcHA,VmaxHA=VmaxHA,GamHA=GamHA)%:%
              foreach(Tair=Tair,LAI=LAI,wind=wind,humidity=humidity,Precip=Precip
              ,SWC_1=SWC_1,Rn=Rn,H=H,LE=LE,pressure=pressure,SZA=SZA,SW_down=SW_down
              ,.errorhandling = "remove",.combine = "rbind")%dopar%{
                
                ############################ 
                # MODEL INPUTS AND CONSTANTS
                
                  
  
                # Constants for Leaf nitrogen taken from CLM4.5 & 5
                if(crop.type=="Wheat"|crop.type=="Sugar beet"){
                SLA<-0.07 # conversion factor from LAI to SLA
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.79
                Ci.ratio<-0.87
                }else if(crop.type=="Soybean"){
                SLA<-0.07 
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.79
                Ci.ratio<-0.87
                }else if(crop.type=="Rice"){
                Gs.max<-5.79
                SLA<-0.07
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.87
                }else if(crop.type=="Maize"|crop.type=="Test 1"|crop.type=="Test 5"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.1758
                CN<-25 
                Ci.ratio<-0.67
                }else if(crop.type=="C3"){ # any other crop type
                SLA<-0.07
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.25
                Ci.ratio<-0.87
                }else if(crop.type=="Test 2"){
                  Gs.max<-1.79
                  SLA<-0.08
                  FLnr<-0.1758
                  CN<-25
                  Ci.ratio<-0.67
                }else if(crop.type=="Test 3"){
                  Gs.max<-1.79
                  SLA<-0.09
                  FLnr<-0.1758
                  CN<-25
                  Ci.ratio<-0.67
                }else if(crop.type=="Test 4"){
                  Gs.max<-1.79
                  SLA<-0.1
                  FLnr<-0.1758
                  CN<-25
                  Ci.ratio<-0.67
                }else if(crop.type=="Test 6"){
                  Gs.max<-1.79
                  SLA<-0.07
                  FLnr<-0.16
                  CN<-25
                  Ci.ratio<-0.67
                }else if(crop.type=="Test 7"){
                  Gs.max<-1.79
                  SLA<-0.07
                  FLnr<-0.15
                  CN<-25
                  Ci.ratio<-0.67
                }else if(crop.type=="Test 8"){
                  Gs.max<-1.79
                  SLA<-0.07
                  FLnr<-0.14
                  CN<-25
                  Ci.ratio<-0.67
                }else{
                SLA<-0.07
                FLnr<-0.1758
                CN<-25
                Gs.max<-1.62
                Ci.ratio<-0.67
                }
                
                # Compute Vcmax same as CLM4.5 & CLM 5 but using Kcat
             
                Na<-(1/(SLA*CN))/14.0057  # CLM uses SLA derived from top canopy to approximate Na.aver for entire canopy
                Co=7/5500 # Houborg et al (2013) and Friend (1995) factor mole Rubisco site per mole N

                Vcmax<-Co*Kcat*FLnr*Na*(10^6) # Houborg et al (2013) and Friend (1995)
              
                
                # Vcmax<-Vcmax+(0.2*Vcmax) # modified Vcmax
                
                Kn<-0.3 # assumed 0.3 in CLM 5 and 4.5 for all canopy types (except forest and tropical canopies)
                
                Tair1<-273.15+Tair
                Tair.atm<-Tair1+(0.0098*Sensor.height)
                Esat.air<-Esat(Temp=Tair.atm,pressure = pressure)
                E<-(Esat.air$`Esat(Pascals)`*humidity)
                Pm<-29.2
                
                # Two stream approximation for canopy radiation
                Canopy.radiation.visible<-Two.stream.approximation.version.3(Pl=0.10,Tl=0.05,LAI=LAI,SZA=SZA,Crop.type=leaf.angle,CI=1,SW.down = SW_down,Pg=0.1)
                
                
                # Soil moisture 
                Stress<-Soil.stress.parameter(Soil.moist = (SWC_1/100),soil.type = "Silty clay loam")
                Gs.stress<-Stress*Gs.max

                # Heat capacity of leaf/canopy
                Cwat<-4188 # specific heat of water
                Cdry<-Cwat*0.8 
                CL<-(LMA/(1-0.8))*((1-0.8)*Cdry+0.8*Cwat)
                
                if(is.na(Rn)|is.na(H)|is.na(LE)){
                  Leaf.temp.guess=Tair.atm
                }else{
                  # leaf temperature
                  Leaf.flux<-Rn-H-LE  
                  Leaf.temp.guess<-Tair.atm+(Leaf.flux/CL) # Initial leaf temp guess
                  # Leaf/canopy temp energy balance iteration
                }

                                 # Partition  kinetics across canopy
                Vcmax_modelled<-Sunlit.shaded.parameter(Param = Vcmax,Kb=Canopy.radiation.visible$Kb,LAI=LAI,Kn=Kn)
                
                
                
                # Leaf boundary layer conductance
                Dh<-18.9e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
                Dm<-13.3e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
                Pr<-Dm/Dh
                a<-(1.5*0.66*Pm*Dh*(Pr^0.33))/(Dm^0.5)
                Gb<-a*((wind/0.05)^0.5) # leaf width must be expressed in mm
                
                # Saturated vapor pressure of leaf/canopy
                Esat.leaf<-Esat(Temp=Leaf.temp.guess,pressure = pressure)
                
                # Initial intercellular CO2 guess
                Ci.guess<-Ci.ratio*Ca  # assumed ambient C02 concentrations at US-Bondville site
                
                ################## Shaded fraction
                An.shaded<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                              ,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA,Temp = Leaf.temp.guess,c.i =Ci.guess)
                Cs.shaded<-Ca-(1.4/Gb)*An.shaded
                Gs.shaded<-Stomatal.quadratic.equation(Stress,Gb,An.shaded,0.01,Gs.stress,Cs.shaded,E,Esat.leaf$`Esat(Pascals)`)
                Ciguess.shaded<-(Cs.shaded-(1.65/Gs.shaded)*An.shaded)-Ci.guess
                An.i.shaded<-Farquhar.version.3(pathway = pathway,Kc25=Kc,Sc.o = Sc.o,c.i=Ciguess.shaded+Ci.guess,Vcmax25=Vcmax_modelled$Paramshaded,PPAR=Canopy.radiation.visible$Shaded
                                                ,Temp=Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
                Cs.i.shaded<-Ca-(1.4/Gb)*An.i.shaded
                Gs.i.shaded<-Stomatal.quadratic.equation(Stress,Gb,An.i.shaded,0.01,Gs.stress,Cs.i.shaded,E,Esat.leaf$`Esat(Pascals)`)
                Ci.intermediate.shaded<-(Cs.i.shaded-(1.65/Gs.i.shaded)*An.i.shaded)-(Ciguess.shaded+Ci.guess)
                Ci.final.shaded<-Ci.iteration(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ciguess.shaded,Ci.intermediate=Ci.intermediate.shaded,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress = Stress,
                                              Vcmax =Vcmax_modelled$Paramshaded,Temp=Leaf.temp.guess,PPAR = Canopy.radiation.visible$Shaded,Kc25=Kc,Sc.o = Sc.o,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
                
                
                An.shaded.final<-Farquhar.version.3(pathway=pathway,c.i=Ci.final.shaded,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                                    ,Temp = Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
                
                ################### Sunlit fraction
                An.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                              ,Temp = Leaf.temp.guess,c.i =Ci.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
                Cs.sunlit<-Ca-(1.4/Gb)*An.sunlit
                Gs.sunlit<-Stomatal.quadratic.equation(Stress,Gb,An.sunlit,0.01,Gs.stress,Cs.sunlit,E,Esat.leaf$`Esat(Pascals)`)
                Ciguess.sunlit<-(Cs.sunlit-(1.65/Gs.sunlit)*An.sunlit)-Ci.guess
                An.i.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,c.i=Ciguess.sunlit+Ci.guess,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR=Canopy.radiation.visible$Sunlit
                                                ,Temp=Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
                Cs.i.sunlit<-Ca-(1.4/Gb)*An.i.sunlit
                Gs.i.sunlit<-Stomatal.quadratic.equation(Stress,Gb,An.i.sunlit,0.01,Gs.stress,Cs.i.sunlit,E,Esat.leaf$`Esat(Pascals)`)
                Ci.intermediate.sunlit<-(Cs.i.sunlit-(1.65/Gs.i.sunlit)*An.i.sunlit)-(Ciguess.sunlit+Ci.guess)
                Ci.final.sunlit<-Ci.iteration(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ciguess.sunlit,Ci.intermediate=Ci.intermediate.sunlit,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress=Stress,
                                              Vcmax =Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit,Kc25=Kc,Sc.o = Sc.o,Temp = Leaf.temp.guess,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
                
                
                An.sunlit.final<-Farquhar.version.3(pathway=pathway,c.i = Ci.final.sunlit,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                                    ,Temp = Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
 
                Count<-Count+1

                Progress<-(Count/(File.length*No.kinetic.sets))*100
                
                print(Progress)
                

                as.data.frame(An.sunlit.final+An.shaded.final)
              
                
              }
          
          return(An)
}
  
  
  
  Stats.sum<-function(Final.An=Final.An){
    Data<-as.data.frame(Final.An)
    Data[Data>100]<-NA # remove extreme values 
    Data[Data<0]<-NA   # remove extreme values
    Sum<-Data %>% summarise_each(funs(sum(.,na.rm=T)))
    Sumcol<-data.frame(Sum=unlist(Sum, use.names = FALSE))
    return(Sumcol)
  }
  
  Stats.RMSE.MAE.R2<-function(Predict=Predict.An,True=An,No.obs=328){
    Predict<-as.data.frame(unlist(Predict))
    Predict[Predict>100]<-NA
    Predict[Predict<0]<-NA
    True<-as.data.frame(True)
    True[True<0]<-NA
    
    Sum.diff<-colSums((abs(Predict-True))^2,na.rm = T)
    RMSE<-sqrt(((1/No.obs)*Sum.diff))
    Sum.dif<-colSums((abs(Predict-True)),na.rm = T)
    MAE<-(1/No.obs)*Sum.dif
    R2<-cor(Predict,True,use = "pairwise.complete.obs")^2
    data<-data.frame(MAE=MAE,RMSE=RMSE,R2=R2)

    return(data)
  }
  
  
  
  
  LAI.seasonal.Xin.method<-function(LAIs.guess=0.01,GPP.reported="Yes"){
    
    if(GPP.reported=="Yes")
    {
      LAI<-(lag(LAIs.guess)+LAIs.guess)/2
    }  
  return(LAI)
    }
  

  
  Two.stream.approximation.version.3<-function(Pl=0.10,Tl=0.05,LAI=6,SZA_=30
                                               ,Crop.type="Spherical",CI=1,SW.down=500,Pg=0.1){
    Pg=Pg # This is P_gd
    CI=CI # This is clumping index
    Pl=Pl # This is ??l
    Tl=Tl # This is ??l
    SZA=SZA_*pi/180 # convert to radians
    L=LAI # This is LAI
    
    # Use Ross index quantifies leaf angle departure from spherical distribution 
    if(Crop.type=="Crop"){   # CLM 4.5 
      XL=-0.30
    }else if(Crop.type=="Cereal"){
      XL=-0.30 
    }else if(Crop.type=="Soybean"){
      XL=-0.5 
    }else if(Crop.type=="Corn"){
      XL=-0.5
    }else{# spherical
      XL=0
    } 
    
    # Direct beam coefficient calculated assuming absence of azimuth preference
    O1=0.5-0.633*XL-0.33*(XL^2) # This is ???1
    O2=0.877*(1-2*O1) # This is ???2
    G=O1+O2*cos(SZA) # This is G(Z)
    Kb<-G/cos(SZA)
    
    # Diffuse beam coefficient for 30 degree incriment sky hemisphere
    G30=O1+O2*cos(0.5235988)
    Kb30<-G30/cos(0.5235988)
    Td30<-2*exp(-Kb30*L)*sin(0.5235988)*cos(0.5235988)*0.5235988 # 0.5235988 refers to 
    
    G60=O1+O2*cos(1.047198)
    Kb60<-G60/cos(1.047198)
    Td60<-2*exp(-Kb60*L)*sin(1.047198)*cos(1.047198)*0.5235988
    
    G90=O1+O2*cos(1.570796)
    Kb90<-G90/cos(1.570796)
    Td90<-2*exp(-Kb90*L)*sin(1.570796)*cos(1.570796)*0.5235988
    
    if(is.infinite(Td90)|is.na(Td90)){
      Td=Td30+Td60
    }else if(is.infinite(Td60)|is.na(Td60) & is.infinite(Td90)|is.na(Td90)){
      Td=Td30
    }else if(is.infinite(Td60)|is.na(Td60)&is.infinite(Td90)|is.na(Td90)&is.infinite(Td30)|is.na(Td30)){
      Td=0  
    }else{
      Td=Td30+Td60+Td90 
    }
    Kd<--(log(Td)/L)
    
    # Two stream approximation by Liou et al (2002)
    W=Pl+Tl
    Cos.angle=((1+XL)/2)^2
    B=((1/2)*(Pl+Tl+(Pl-Tl)*Cos.angle))/W # This is ????l
    a_s=(W/2)*(G/(G+cos(SZA)*O2) * (1 - (cos(SZA)*O1/(G+cos(SZA)*O2)) * log((G+cos(SZA)*O1+cos(SZA)*O2)/(cos(SZA)*O1))))# This is a_s (??)
    Bo=(((Kb+Kd)/Kb)*a_s)/W  # This is ??_0 ??l
    
    # Partition Direct and diffuse solar radiation
    SW.diffuse=SW.down*0.2 # This gives I_(sky,d)^???
    SW.direct=SW.down-SW.diffuse # This gives I_(sky,b)^???
    
    # Coefficients
    b=(1-(1-B)*W)*Kd
    c=B*W*Kd
    h.a=(b^2)-(c^2)
    h=sqrt((b^2)-(c^2))
    u=h-b-c/(2*h)
    v=h+b+c/(2*h)
    
    
    # Direct beam flux 
    y1=((Bo*Kb-b*Bo-c*(1-Bo))*W*Kb*SW.direct)/((h^2)-(Kb^2))
    y2=(((1-Bo)*Kb+c*Bo+b*(1-Bo))*W*Kb*SW.direct)/((h^2)-(Kb^2))
    
    n2=(v*(y1+y2*Pg+Pg*SW.direct)*exp(-Kb*CI*L)-y2*(u+v*Pg)*exp(-h*CI*L))/(v*(v+u*Pg)*exp(h*CI*L)-u*(u+v*Pg)*exp(-h*CI*L))
    n1=y2-n2*u/v
    
    SW.up.direct.above=-(y1*exp(-Kb*CI*0))+(n1*u*exp(-h*CI*0))+(n2*v*exp(h*CI*0)) # upward flux of direct beam solar radiation above the leaves with L=0
    SW.down.direct.below=(y2*exp(-Kb*CI*L))-(n1*v*exp(-h*CI*L))-(n2*u*exp(h*CI*L)) # downward flux of direct beam solar radiation below the leaves at depth equal to the cumulative leaf area index (L) with L=LAI 
    SW.up.direct.L=-(y1*exp(-Kb*CI*L))+(n1*u*exp(-h*CI*L))+(n2*v*exp(h*CI*L)) #  scattered flux of direct beam solar radiation upward in the canopy at depth equal to the cumulative leaf area index (L) with L=LAI
    SW.canopy.b=(1-exp(-Kb*CI*L))*SW.direct-SW.up.direct.above+SW.up.direct.L-SW.down.direct.below # The direct beam flux absorbed by the canopy is obtained from the radiative balance by considering fluxes above and below the canopy. 
    
    # Diffuse beam flux      
    n2.diffuse=SW.diffuse*(u+v*Pg)*exp(-h*CI*L)/(v*(v+u*Pg)*exp(h*CI*L)-u*(u+v*Pg)*exp(-h*CI*L))
    n1.diffuse=-SW.diffuse-n2.diffuse*u/v
    
    SW.up.diffuse.above=(n1.diffuse*u*exp(-h*CI*0))+(n2.diffuse*v*exp(h*CI*0)) # upward flux of diffuse beam solar radiation above the leaves with L=0
    SW.down.diffuse.below=(-n1.diffuse*v*exp(-h*CI*L))-(n2.diffuse*u*exp(h*CI*L)) # downward flux of diffuse beam solar radiation below the leaves at depth equal to the cumulative leaf area index (L) with L=LAI 
    SW.up.diffuse.L=(n1.diffuse*u*exp(-h*CI*L))+(n2.diffuse*v*exp(h*CI*L)) #  scattered flux of diffuse beam solar radiation upward in the canopy at depth equal to the cumulative leaf area index (L) with L=LAI
    SW.canopy.d=SW.diffuse-SW.up.diffuse.above+SW.up.diffuse.L-SW.down.diffuse.below # # The diffuse beam flux absorbed by the canopy is obtained from the radiative balance by considering fluxes above and below the canopy. 
    
    
    # Calculate sunlit/shaded fractions of LAI
    Lsun=((1-exp(-Kb*L)))/Kb 
    Lshad=L-Lsun
    
    
    
    # Shaded solar radiation per leaf area
    Shaded_=SW.canopy.d/Lshad # Overall having difficulties correcting sunlit/shaded bug. 
    # Used basic assumption shaded leaves recieve diffusion only and sunlit recieves both direct and diffusion.
    
    # Sunlit solar radiation per leaf area
    Sunlit_=Shaded_+(SW.canopy.b/Lsun) # This is 
    
    #Scale to canopy
    Shaded<-as.data.frame(Shaded_*4.6*Lshad) #Scale to canopy and conver to umol photons 
    colnames(Shaded)<-"Shaded"
    Sunlit<-as.data.frame(Sunlit_*4.6*Lsun) #Scale to canopy and conver to umol photons 
    colnames(Sunlit)<-"Sunlit"
    Total.canopy<-as.data.frame(Sunlit+Shaded)
    colnames(Total.canopy)<-"Total"
    
    Kb<-as.data.frame(Kb)
    colnames(Kb)<-"Kb"
    
    Data<-cbind(Kb,Total.canopy,Sunlit,Shaded)
    
    return(Data)
  }
  
  
  



VPD_GSI_Jolly=function(VPD=c(1,5)){
  
  VPDmax=2.5 # adjusted to Ball & Berry response
  VPDmin=0.5 # adjusted to Ball & Berry response 
  
  GSI_VPD=ifelse(VPD>VPDmax,0.01,ifelse(VPDmax>VPD&VPD>VPDmin,1-((VPD-VPDmin)/(VPDmax-VPDmin)),1))
  
  return(GSI_VPD)
  
}
  




Two.big.leaf.concept.with.yield<-function(File.length=1,No.kinetic.sets=1,LAIMax=2.2,crop.type="Wheat",leaf.angle="Wheat",pathway="C3",Temperature=25,wind=2,humidity=80/100,Precip=0.2
                                          ,VPD=0.1,SWC_1=0.3,pressure=100,SZA=67,SW_down=500
                                          ,Ca=360,Sensor.height=2.7,Kcat=2,Kc=10,Sc.o=101,KcHA=20,VmaxHA=30,GamHA=40){
  Count<-0
  
  
  # Initial parameters for GDD phenology
  fvn=0 # cumulative daily vernalization rate
  GDD_O<-10 # wheat and sugarbeet initial GDDmin based on 20-year running mean growing degree days for planting tracked from April to September
  GDD_8<-5 # maize initial GDD based on 20-year running mean growing degree days for planting tracked from April to September
  i_Carbon=3 # Initial carbon in seed prior to leef emergence assumed to be 3 g Carbon m^-2
  LAI=0.001
  Temperature<-273.15+Temperature
  LAIMax=LAIMax
  VPD_GSI_stress=VPD_GSI_Jolly(VPD=VPD)
  
  
  An<-foreach(Kcat=Kcat,Kc=Kc,Sc.o=Sc.o,KcHA=KcHA,VmaxHA=VmaxHA,GamHA=GamHA)%:%
    foreach(Tair=Temperature,wind=wind,humidity=humidity,Precip=Precip
            ,VPD_stress=VPD_GSI_stress,SWC_1=SWC_1,pressure=pressure,SZA=SZA,SW_down=SW_down
            ,.errorhandling = "remove",.combine = "rbind")%dopar%{
              
              
              ############################ 
              # MODEL INPUTS AND CONSTANTS
              
              # Constants for Leaf nitrogen taken from CLM4.5 & 5
              if(crop.type=="Wheat"|crop.type=="Sugar beet"){
                SLA<-0.035 # conversion factor from LAI to SLA
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.25
                Ci.ratio<-0.87
              }else if(crop.type=="Maize"|crop.type=="Test 1"|crop.type=="Test 5"){
                SLA<-0.09
                FLnr<-0.1758
                Gs.max<-1.79
                Ci.ratio<-0.67
              }else if(crop.type=="C3"){ # any other crop type
                SLA<-0.035
                FLnr<-0.4120
                Gs.max<-5.25
                Ci.ratio<-0.87
              }else{
                SLA<-0.07
                FLnr<-0.1758
                Gs.max<-5.25
                Ci.ratio<-0.67
              }
              
              #### GROWING DEGREE DAYS AND LAI CALCULATION
              
              if(crop.type=="Maize"){
                GDD_8=GDD_8+Tair-273.15-8 # Unlike CLM5.0 we have chosen to use daily mean Tair instead of total Tair in day.
              }else{
                GDD_O=GDD_O+Tair-273.15 # Unlike CLM5.0 we have chosen to use daily mean Tair instead of total Tair in day. 
              }
              
              
              
              if(crop.type=="Maize"){ # MAX GDDmat values or can modify to stop once reached a specific leaf area index
                 
                GDD_mat<-max(950,min(0.85*5,1700)) # Peng et al (2017)
                GDD_mat<-max(950,min(GDD_mat+150),1850) # Peng et al (2017)
                
              }
              
              if(crop.type=="Sugar beet"){
                GDD_mat<-1700
                
              }
              
              if(crop.type=="Wheat"){
                GDD_mat<-2000
                
              }
              
              if(crop.type=="Wheat"&Tair<273.15){
                Tcrown=2+(Tair-273.15)*(0.4+0.0018*(min(0*100,15)-15)^2)
              }else{
                Tcrown=Tair-273.15
              }
              
  
              if(crop.type=="Potato"){
                GDD_mat=1400
              }
        
              
              
              # Maize Leaf emergence stage
              if(LAI<LAIMax&crop.type=="Maize"&(GDD_mat*0.65)>GDD_8&GDD_8>(0.05*GDD_mat)){
                
                CN=25# Vary CN ratio for development stages
                
                # The 'h' heat unit threshold is the same for leaf emergence and grain filling stage
                CRM=max(min((GDD_mat+53.683)/13.882,135),73)  # Peng et al (2017)
                frac=min(max(-0.002*(CRM-73)+0.65,0.65-0.1),0.65)  # Peng et al (2017)
                h=frac*GDD_mat # Peng et al 2017
                
                arep=0 #### Carbon allocation
                aroot=0.1-(0.1-0.05)*(GDD_8/GDD_mat) 
                aleaf3=(1-aroot)*0.6*(exp(-0.1)-exp(-0.1*(GDD_8/h)))/(exp(-0.1)-1)
                alivestem3=1-arep-aroot-aleaf3
                
                Cleaf=i_Carbon*aleaf3*VPD_stress
                LAI=LAI+(Cleaf*SLA)
              }else{
                
                CN=25  # Vary CN ratio for development stages
                
                
                CRM=max(min((GDD_mat+53.683)/13.882,135),73)  # Peng et al (2017)
                frac=min(max(-0.002*(CRM-73)+0.65,0.65-0.1),0.65)  # Peng et al (2017)
                h=frac*GDD_mat # Peng et al 2017
                
                arep=0 #### Carbon allocation
                aroot=0.1-(0.1-0.05)*(GDD_8/GDD_mat) 

                Croot=i_Carbon*aroot # if max LAI met then carbon is directed exclusively to roots during leaf emergence
              }
              
              
              
              # C3 Leaf emergence stage
              if(LAI<LAIMax&crop.type=="Sugar beet"&(GDD_mat*0.6)>GDD_O&GDD_O>100){
                
                CN<-25 # Vary CN ratio for development stages
                
                h=0.6*GDD_mat
                
                arep=0 #### Carbon allocation
                aroot=0.05-(0.05-0)*(GDD_O/GDD_mat)
                aleaf3=(1-aroot)*0.9*(exp(-0.1)-exp(-0.1*(GDD_O/(h))))/(exp(-0.1)-1)
                alivestem3=1-arep-aroot-aleaf3 
                
                Cleaf=i_Carbon*aleaf3*VPD_stress
                LAI=LAI+(Cleaf*SLA) # leaf area index calculation using SLA 0.035 m^2 leaf g-1 Carbon
                
              }else{
                
                CN<-25 # Vary CN ratio for development stages
                
                h=0.6*GDD_mat
                
                arep=0 #### Carbon allocation
                aroot=0.05-(0.05-0)*(GDD_O/GDD_mat)
                
                Croot=i_Carbon*aroot # if max LAI met then carbon is directed exclusively to roots during leaf emergence
                
              }
              
              if(crop.type=="Wheat"){
                
                Tmin=-1.3 # -1.3 degrees from Streck et al (2003)
                Tmax=15.7 # 15.7 degrees from Streck et al (2003)
                Topt=4.9 # 4.9 degrees from Streck et al (2003)
                
                if(Tcrown>Tmax | Tcrown<Tmin){
                  fvn=fvn+0
                }else if(Tcrown>Tmin & Tcrown<Tmax){
                  
                  a=log(2)/log((Tmax-Tmin)/(Topt-Tmin))
                  b=(2*(Tcrown-Tmin)^a*(Topt-Tmin)^a-(Tcrown-Tmin)^(2*a))/((Topt-Tmin)^(2*a))
                  fvn=fvn+b
                  
                }else{
                  fvn=fvn+1
                }
                VF=(fvn^5)/(22.5^5+fvn^5) # Vernalization factor revised for mean daily air temperature
                
                
              }
              
              
              
              if(LAI<LAIMax&crop.type=="Wheat"&(GDD_mat*0.6)>GDD_O&GDD_O>60){
                
                CN<-50 # Vary CN ratio for development stages
                
                GDD_O=GDD_O*VF # adjusted for Winter wheat vernalization
                
                h=0.6*GDD_mat
                
                arep=0 #### Carbon allocation
                aroot=0.05-(0.05-0)*(GDD_O/GDD_mat)
                aleaf3=(1-aroot)*0.9*(exp(-0.1)-exp(-0.1*(GDD_O/(h))))/(exp(-0.1)-1)
                alivestem3=1-arep-aroot-aleaf3 
                
                Cleaf=i_Carbon*aleaf3*VPD_stress
                LAI=LAI+(Cleaf*SLA) # leaf area index calculation using SLA 0.035 m^2 leaf g-1 Carbon
                
              }else{
                
                CN<-50 # Vary CN ratio for development stages
                
                h=0.6*GDD_mat
                
                arep=0 #### Carbon allocation
                aroot=0.05-(0.05-0)*(GDD_O/GDD_mat)
                
                Croot=i_Carbon*aroot # if max LAI met then carbon is directed exclusively to roots during leaf emergence
                
              }
              
              
              # Wheat grain fill stage
              if(crop.type=="Wheat"&GDD_O>(GDD_mat*0.6)&GDD_O<GDD_mat){
                
                CN=50 # Vary CN ratio for development stages
                
                h=0.6*GDD_mat
                
                aleaf=aleaf3*(1-((GDD_O-(h))/((GDD_mat*1.05)-(h))))^3 # something wrong in this step!!!! with 'h' 
                alivestem=alivestem3*(1-((GDD_O-(h))/((GDD_mat*1.05)-(h))))^1
                agrain=(1-0-alivestem-aleaf)*VF # adjusted by vernalization factor if winter wheat is not fully vernalized. 
                
                
                Cleaf=i_Carbon*aleaf*VPD_stress
                LAI=LAI # Also included leaflitter decline function. Once MaxLAI reached in leaf emergence, LAI begins to decline
              }
              

              # Maize grain fill stage
              if(crop.type=="Maize"&GDD_8>(GDD_mat*0.65)&GDD_8<GDD_mat){
                
                CN<-60 # Vary CN ratio for development stages
                
                
                CRM=max(min((GDD_mat+53.683)/13.882,135),73)  # Peng et al (2017)
                frac=min(max(-0.002*(CRM-73)+0.65,0.65-0.1),0.65)  # Peng et al (2017)
                h=frac*GDD_mat # Peng et al (2017)
                
                aleaf=aleaf3*(1-((GDD_8-(h))/((GDD_mat*1.05)-(h))))^5
                alivestem=alivestem3*(1-((GDD_8-(h))/((GDD_mat*1.05)-(h))))^2
                agrain=1-0.05-alivestem-aleaf
                
                Cleaf=i_Carbon*aleaf*VPD_stress 
                LAI=LAI-(LAI*(1/(0.1*365))) # Also included leaflitter decline function
                
              }
              
              
              
              # C3 grain fill stage
              if(crop.type=="Sugar beet"&GDD_O>(GDD_mat*0.6)){
                
                CN=25 # Vary CN ratio for development stages
                
                
                h=0.6*GDD_mat
                
                aleaf=aleaf3*(1-((GDD_O-(h))/((GDD_mat*1.05)-(h))))^3 # something wrong in this step!!!! with 'h' 
                alivestem=alivestem3*(1-((GDD_O-(h))/((GDD_mat*1.05)-(h))))^1
                agrain=1-0-alivestem-aleaf
                
                
                Cleaf=i_Carbon*aleaf*VPD_stress
                LAI=LAI-(LAI*(1/(0.2*365))) # Also included leaflitter decline function. Once MaxLAI reached in leaf emergence, LAI begins to decline
              }
              
              
              # Harvest for C3
              # if(crop.type=="Sugar beet"&GDD_O>GDD_mat){ # If max days after planting or GDD_mat reached harvest occurs
                # LAI=0.001

              # }
              
              if(crop.type=="Wheat"&GDD_O>GDD_mat){
                LAI=0.001
              }
              
              # Harvest for C4
              if(crop.type=="Maize"&GDD_8>GDD_mat){
                LAI=0.001

              }
              
              
              # Compute Vcmax same as CLM4.5 & CLM 5 but using Kcat
              
              Na<-(1/(SLA*CN))/14.0057  # CLM uses SLA derived from top canopy to approximate Na.aver for entire canopy
              Co=7/5500 # Houborg et al (2013) and Friend (1995) factor mole Rubisco site per mole N
              
              Vcmax<-Co*Kcat*FLnr*Na*(10^6) # Houborg et al (2013) and Friend (1995)
              
              Kn<-0.3 # assumed 0.3 in CLM 5 and 4.5 for all canopy types (except forest and tropical canopies)
              
              
              Esat.air<-Esat(Temp=Tair,pressure = pressure)
              E<-(Esat.air$`Esat(Pascals)`*humidity)
              Pm<-29.2
              
              
              
              # Two stream approximation for canopy radiation
              Canopy.radiation.visible<-Two.stream.approximation.version.3(Pl=0.10,Tl=0.05,LAI=LAI,SZA=SZA,Crop.type=leaf.angle,CI=1,SW.down = SW_down,Pg=0.1)
              
              # Soil moisture 
              
              Stress<-Soil.stress.parameter(Soil.moist = (SWC_1/100),soil.type = "Silty clay loam")
              Gs.stress<-Stress*Gs.max
              
              
              
              # Partition  kinetics across canopy
              Vcmax_modelled<-Sunlit.shaded.parameter(Param = Vcmax,Kb=Canopy.radiation.visible$Kb,LAI=LAI,Kn=Kn)
              
              if(is.na(wind)){
                wind=1
                
              }
              
              # Leaf boundary layer conductance
              Dh<-18.9e-06*(101325/(pressure*1000))*((Tair/273.15)^1.81)
              Dm<-13.3e-06*(101325/(pressure*1000))*((Tair/273.15)^1.81)
              Pr<-Dm/Dh
              a<-(1.5*0.66*Pm*Dh*(Pr^0.33))/(Dm^0.5)
              Gb<-a*((wind/0.05)^0.5) # leaf width must be expressed in mm
              
              # Saturated vapor pressure of leaf/canopy
              Esat.leaf<-Esat(Temp=Tair,pressure = pressure)
              
              # Initial intercellular CO2 guess
              Ci.guess<-Ci.ratio*Ca  # assumed ambient C02 concentrations at US-Bondville site
              
              ################## Shaded fraction
              An.shaded<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                            ,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA,Temp = Tair,c.i =Ci.guess)
              Cs.shaded<-Ca-(1.4/Gb)*An.shaded
              Gs.shaded<-Stomatal.quadratic.equation(Stress,Gb,An.shaded,0.01,Gs.stress,Cs.shaded,E,Esat.leaf$`Esat(Pascals)`)
              Ciguess.shaded<-(Cs.shaded-(1.65/Gs.shaded)*An.shaded)-Ci.guess
              An.i.shaded<-Farquhar.version.3(pathway = pathway,Kc25=Kc,Sc.o = Sc.o,c.i=Ciguess.shaded+Ci.guess,Vcmax25=Vcmax_modelled$Paramshaded,PPAR=Canopy.radiation.visible$Shaded
                                              ,Temp=Tair,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
              Cs.i.shaded<-Ca-(1.4/Gb)*An.i.shaded
              Gs.i.shaded<-Stomatal.quadratic.equation(Stress,Gb,An.i.shaded,0.01,Gs.stress,Cs.i.shaded,E,Esat.leaf$`Esat(Pascals)`)
              Ci.intermediate.shaded<-(Cs.i.shaded-(1.65/Gs.i.shaded)*An.i.shaded)-(Ciguess.shaded+Ci.guess)
              Ci.final.shaded<-Ci.iteration(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ciguess.shaded,Ci.intermediate=Ci.intermediate.shaded,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress = Stress,
                                            Vcmax =Vcmax_modelled$Paramshaded,Temp=Tair,PPAR = Canopy.radiation.visible$Shaded,Kc25=Kc,Sc.o = Sc.o,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              
              An.shaded.final<-Farquhar.version.3(pathway=pathway,c.i=Ci.final.shaded,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                                  ,Temp = Tair,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              ################### Sunlit fraction
              An.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                            ,Temp = Tair,c.i =Ci.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              Cs.sunlit<-Ca-(1.4/Gb)*An.sunlit
              Gs.sunlit<-Stomatal.quadratic.equation(Stress,Gb,An.sunlit,0.01,Gs.stress,Cs.sunlit,E,Esat.leaf$`Esat(Pascals)`)
              Ciguess.sunlit<-(Cs.sunlit-(1.65/Gs.sunlit)*An.sunlit)-Ci.guess
              An.i.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,c.i=Ciguess.sunlit+Ci.guess,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR=Canopy.radiation.visible$Sunlit
                                              ,Temp=Tair,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
              Cs.i.sunlit<-Ca-(1.4/Gb)*An.i.sunlit
              Gs.i.sunlit<-Stomatal.quadratic.equation(Stress,Gb,An.i.sunlit,0.01,Gs.stress,Cs.i.sunlit,E,Esat.leaf$`Esat(Pascals)`)
              Ci.intermediate.sunlit<-(Cs.i.sunlit-(1.65/Gs.i.sunlit)*An.i.sunlit)-(Ciguess.sunlit+Ci.guess)
              Ci.final.sunlit<-Ci.iteration(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ciguess.sunlit,Ci.intermediate=Ci.intermediate.sunlit,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress=Stress,
                                            Vcmax =Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit,Kc25=Kc,Sc.o = Sc.o,Temp = Tair,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              An.sunlit.final<-Farquhar.version.3(pathway=pathway,c.i = Ci.final.sunlit,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                                  ,Temp = Tair,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              Count<-Count+1
              
              
              Progress<-(Count/(File.length*No.kinetic.sets))*100
              
              i_Carbon=An.sunlit.final+An.shaded.final
              
              if(abs(i_Carbon) > 30){
                i_Carbon=0.1
              }else{
                i_Carbon=i_Carbon # Total carbon uptake for current time step
              }
           
              
              print(Progress)
              
              as.data.frame(An.sunlit.final+An.shaded.final)
              
            }
  
  return(An)
}





Mesophyll_conductance_equation<-function(Gm_25=0.55,Temp=298.15){
  S_Jmax=0.65970 # Gm uses same temperature function as Jmax
  Jmax_Ha=50 # Gm uses same temperature function as Jmax
  Jmax_Hd=200 # Gm uses same temperature function as Jmax
  
  
  Gm<-ifelse(Temp>298.15,Gm_25 * exp(Jmax_Ha/298.15*8.314*(1-(298.15/Temp))) * (1+exp((298.15*S_Jmax-Jmax_Hd)/(298.15*0.001*8.314)))/(1+exp((S_Jmax*Temp-Jmax_Hd)/(0.001*8.314*Temp))),Gm_25 * exp(Jmax_Ha/298.15*8.314*(1-(298.15/Temp)))) 
  
  return(Gm)
}



Ci.iteration_mesophyll<-function(max.iterations=20,Ci.old=241,fxCi=118.72,Ci.intermediate=3.933,Ca=360,ea=0.01264362790,esat=0.01456363309,Gb=2.416,g1=1.79,Gs.stress=9,Stress=0.7
                                 ,Vcmax=78,PPAR=500,Kc25=269,pathway="c3",Sc.o = 167,Temp=273,VmaxHA=20,GamHA=40,KcHA=50){
  count<-0
  max.iterations<-max.iterations
  Stress<-Stress
  Ci.old<-Ci.old
  Ci.bisec.guess<-Ci.old
  FxCi.old<-fxCi
  Ci.intermediate.old<-Ci.intermediate
  Gb<-Gb
  Gs.stress<-Gs.stress
  Ca<-Ca
  ea<-ea
  esat<-esat
  Vcmax=Vcmax
  PPAR=PPAR
  Kc25=Kc25
  Sc.o = Sc.o
  Temp=Temp
  VmaxHA=VmaxHA
  GamHA=GamHA
  KcHA=KcHA
  pathway=pathway
  # Mesophyll conductance adjusted for temperature. Assume 0.55 mol CO2 m-2 s-1 for C3 wheat at 25oC
  Gm<-Mesophyll_conductance_equation(Gm_25=0.55,Temp=Temp)
  
  repeat{
    
    Ci.new<-Ci.old-(FxCi.old/((Ci.intermediate.old-FxCi.old)/FxCi.old))
    
    
    An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = Ci.new,
                           Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
    Cs<-Ca-(1.4/Gb)*An
    Gs<-Stomatal.quadratic.equation(Stress,Gb,An,0.01,Gs.stress,Cs,ea,esat)
    Ci<-(Cs-(1.65/Gs)*An)-Ci.new # f(x)
    FxCi.old<-Ci-(An/Gm)
    
    An.i<-Farquhar.version.3(Kc25 = Kc25,Sc.o = Sc.o,c.i = FxCi.old+Ci.new,
                             Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate f'(x)
    Cs.i<-Ca-(1.4/Gb)*An.i
    Gs.i<-Stomatal.quadratic.equation(Stress,Gb,An.i,0.01,Gs.stress,Cs.i,ea,esat)
    Ci.i<-(Cs.i-(1.65/Gs.i)*An.i)-(FxCi.old+Ci.new)
    Cc.i<-Ci.i-(An.i/Gm)
    Ci.intermediate.old<-Cc.i
    
    Conv<-abs(Ci.new-Ci.old)
    ifelse(Conv<0.0001,break,Ci.new)
    Ci.old<-Ci.new
    
    count<-count+1
    
    ifelse(count>max.iterations,break,NA) 
  }
  ### Bisection backup method
  Conv<-Conv # break next loop if first loop is successful
  count<-0 # Reset the count
  a<-5 # min interval guess
  An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = a,
                         Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
  Cs<-Ca-(1.4/Gb)*An
  Gs<-Stomatal.quadratic.equation(Stress,Gb,An,0.01,Gs.stress,Cs,ea,esat)
  Fx<-(Cs-(1.65/Gs)*An)-a
  Fx.a<-Fx-(An/Gm)
  
  b<-400 # max interval guess
  
  repeat{
    ifelse(Conv<0.0001,break,NA)
    c=(a+b)/2
    
    An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = c,
                           Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
    Cs<-Ca-(1.4/Gb)*An
    Gs<-Stomatal.quadratic.equation(Stress,Gb,An,0.01,Gs.stress,Cs,ea,esat)
    Fx<-(Cs-(1.65/Gs)*An)-c
    Fx.c<-Fx-(An/Gm)
    
    ifelse(abs(Fx.c)<0.0001,break,c)
    
    a=ifelse((Fx.a*Fx.c)<0,a,c)
    b=ifelse((Fx.a*Fx.c)<0,c,400)
    Fx.a=ifelse((Fx.a*Fx.c)>0,Fx.c,Fx.a)
    
    Ci.new<-a
    count<-count+1
    ifelse(count==max.iterations,break,NA)
    
  }
  return(Ci.new)
  options(digits = 20)
}



Two.big.leaf.concept_mesophyll_stomatal<-function(File.length=1,No.kinetic.sets=1,LAI=4,crop.type="Soybean",leaf.angle="Soybean",pathway="C3",Tair=25,wind=2,humidity=80/100,Precip=0.2
                                                  ,SWC_1=0.3,Rn=NA,H=NA,LE=NA,LMA=2.2,pressure=100,SZA=67,SW_down=500
                                                  ,Ca=360,Sensor.height=2.7,Kcat=2,Kc=10,Sc.o=101,KcHA=20,VmaxHA=30,GamHA=40){
  Count<-0
  
  An<-foreach(Kcat=Kcat,Kc=Kc,Sc.o=Sc.o,KcHA=KcHA,VmaxHA=VmaxHA,GamHA=GamHA)%:%
    foreach(Tair=Tair,LAI=LAI,wind=wind,humidity=humidity,Precip=Precip
            ,SWC_1=SWC_1,Rn=Rn,H=H,LE=LE,pressure=pressure,SZA=SZA,SW_down=SW_down
            ,.errorhandling = "remove",.combine = "rbind")%dopar%{
              
              ############################ 
              # MODEL INPUTS AND CONSTANTS
              
              
              
              # Constants for Leaf nitrogen taken from CLM4.5 & 5
              if(crop.type=="Wheat"|crop.type=="Sugar beet"){
                SLA<-0.07 # conversion factor from LAI to SLA
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.79
                Ci.ratio<-0.87
              }else if(crop.type=="Soybean"){
                SLA<-0.07 
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.79
                Ci.ratio<-0.87
              }else if(crop.type=="Rice"){
                Gs.max<-5.79
                SLA<-0.07
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.87
              }else if(crop.type=="Maize"|crop.type=="Test 1"|crop.type=="Test 5"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.1758
                CN<-25 
                Ci.ratio<-0.67
              }else if(crop.type=="C3"){ # any other crop type
                SLA<-0.07
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.25
                Ci.ratio<-0.87
              }else if(crop.type=="Test 2"){
                Gs.max<-1.79
                SLA<-0.08
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 3"){
                Gs.max<-1.79
                SLA<-0.09
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 4"){
                Gs.max<-1.79
                SLA<-0.1
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 6"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.16
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 7"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.15
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 8"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.14
                CN<-25
                Ci.ratio<-0.67
              }else{
                SLA<-0.07
                FLnr<-0.1758
                CN<-25
                Gs.max<-1.62
                Ci.ratio<-0.67
              }
              
              # Compute Vcmax same as CLM4.5 & CLM 5 but using Kcat
              
              Na<-(1/(SLA*CN))/14.0057  # CLM uses SLA derived from top canopy to approximate Na.aver for entire canopy
              Co=7/5500 # Houborg et al (2013) and Friend (1995) factor mole Rubisco site per mole N
              
              Vcmax<-Co*Kcat*FLnr*Na*(10^6) # Houborg et al (2013) and Friend (1995)
              
              
              # Vcmax<-Vcmax+(0.2*Vcmax) # modified Vcmax
              
              Kn<-0.3 # assumed 0.3 in CLM 5 and 4.5 for all canopy types (except forest and tropical canopies)
              
              Tair1<-273.15+Tair
              Tair.atm<-Tair1+(0.0098*Sensor.height)
              Esat.air<-Esat(Temp=Tair.atm,pressure = pressure)
              E<-(Esat.air$`Esat(Pascals)`*humidity)
              Pm<-29.2
              
              # Two stream approximation for canopy radiation
              Canopy.radiation.visible<-Two.stream.approximation.version.3(Pl=0.10,Tl=0.05,LAI=LAI,SZA=SZA,Crop.type=leaf.angle,CI=1,SW.down = SW_down,Pg=0.1)
              
              
              # Soil moisture 
              Stress<-Soil.stress.parameter(Soil.moist = (SWC_1/100),soil.type = "Silty clay loam")
              Gs.stress<-Stress*Gs.max
              
              # Heat capacity of leaf/canopy
              Cwat<-4188 # specific heat of water
              Cdry<-Cwat*0.8 
              CL<-(LMA/(1-0.8))*((1-0.8)*Cdry+0.8*Cwat)
              
              if(is.na(Rn)|is.na(H)|is.na(LE)){
                Leaf.temp.guess=Tair.atm
              }else{
                # leaf temperature
                Leaf.flux<-Rn-H-LE  
                Leaf.temp.guess<-Tair.atm+(Leaf.flux/CL) # Initial leaf temp guess
                # Leaf/canopy temp energy balance iteration
              }
              
              # Partition  kinetics across canopy
              Vcmax_modelled<-Sunlit.shaded.parameter(Param = Vcmax,Kb=Canopy.radiation.visible$Kb,LAI=LAI,Kn=Kn)
              
              
              
              # Leaf boundary layer conductance
              Dh<-18.9e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
              Dm<-13.3e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
              Pr<-Dm/Dh
              a<-(1.5*0.66*Pm*Dh*(Pr^0.33))/(Dm^0.5)
              Gb<-a*((wind/0.05)^0.5) # leaf width must be expressed in mm
              
              # Saturated vapor pressure of leaf/canopy
              Esat.leaf<-Esat(Temp=Leaf.temp.guess,pressure = pressure)
              
              # Initial intercellular CO2 guess
              Ci.guess<-Ci.ratio*Ca  # assumed ambient C02 concentrations at US-Bondville site
              
              # Mesophyll conductance adjusted for temperature. Assume 0.55 mol CO2 m-2 s-1 for C3 wheat at 25oC
              Gm<-Mesophyll_conductance_equation(Gm_25=0.55,Temp=Leaf.temp.guess)
              
              ################## Shaded fraction
              An.shaded<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                            ,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA,Temp = Leaf.temp.guess,c.i =Ci.guess)
              Cs.shaded<-Ca-(1.4/Gb)*An.shaded
              Gs.shaded<-Stomatal.quadratic.equation(Stress,Gb,An.shaded,0.01,Gs.stress,Cs.shaded,E,Esat.leaf$`Esat(Pascals)`)
              Ciguess.shaded<-(Cs.shaded-(1.65/Gs.shaded)*An.shaded)-Ci.guess
              Ccguess.shaded<-Ciguess.shaded-(An.shaded/Gm) # account for mesophyll
              An.i.shaded<-Farquhar.version.3(pathway = pathway,Kc25=Kc,Sc.o = Sc.o,c.i=Ccguess.shaded+Ci.guess,Vcmax25=Vcmax_modelled$Paramshaded,PPAR=Canopy.radiation.visible$Shaded
                                              ,Temp=Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
              Cs.i.shaded<-Ca-(1.4/Gb)*An.i.shaded
              Gs.i.shaded<-Stomatal.quadratic.equation(Stress,Gb,An.i.shaded,0.01,Gs.stress,Cs.i.shaded,E,Esat.leaf$`Esat(Pascals)`)
              Ci.intermediate.shaded<-(Cs.i.shaded-(1.65/Gs.i.shaded)*An.i.shaded)-(Ciguess.shaded+Ci.guess)
              Cc.intermediate.shaded<-Ci.intermediate.shaded-(An.i.shaded/Gm)
              Cc.final.shaded<-Ci.iteration_mesophyll(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ccguess.shaded,Ci.intermediate=Cc.intermediate.shaded,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress = Stress,
                                                      Vcmax =Vcmax_modelled$Paramshaded,Temp=Leaf.temp.guess,PPAR = Canopy.radiation.visible$Shaded,Kc25=Kc,Sc.o = Sc.o,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              
              An.shaded.final<-Farquhar.version.3(pathway=pathway,c.i=Cc.final.shaded,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                                  ,Temp = Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              ################### Sunlit fraction
              An.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                            ,Temp = Leaf.temp.guess,c.i =Ci.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              Cs.sunlit<-Ca-(1.4/Gb)*An.sunlit
              Gs.sunlit<-Stomatal.quadratic.equation(Stress,Gb,An.sunlit,0.01,Gs.stress,Cs.sunlit,E,Esat.leaf$`Esat(Pascals)`)
              Ciguess.sunlit<-(Cs.sunlit-(1.65/Gs.sunlit)*An.sunlit)-Ci.guess
              Ccguess.sunlit<-Ciguess.sunlit-(An.sunlit/Gm) # account for mesophyll
              An.i.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,c.i=Ccguess.sunlit+Ci.guess,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR=Canopy.radiation.visible$Sunlit
                                              ,Temp=Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
              Cs.i.sunlit<-Ca-(1.4/Gb)*An.i.sunlit
              Gs.i.sunlit<-Stomatal.quadratic.equation(Stress,Gb,An.i.sunlit,0.01,Gs.stress,Cs.i.sunlit,E,Esat.leaf$`Esat(Pascals)`)
              Ci.intermediate.sunlit<-(Cs.i.sunlit-(1.65/Gs.i.sunlit)*An.i.sunlit)-(Ciguess.sunlit+Ci.guess)
              Cc.intermediate.sunlit<-Ci.intermediate.sunlit-(An.i.sunlit/Gm)
              
              Cc.final.sunlit<-Ci.iteration_mesophyll(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ccguess.sunlit,Ci.intermediate=Cc.intermediate.sunlit,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress=Stress,
                                                      Vcmax =Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit,Kc25=Kc,Sc.o = Sc.o,Temp = Leaf.temp.guess,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              
              An.sunlit.final<-Farquhar.version.3(pathway=pathway,c.i = Cc.final.sunlit,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                                  ,Temp = Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              Count<-Count+1
              
              Progress<-(Count/(File.length*No.kinetic.sets))*100
              
              print(Progress)
              
              
              as.data.frame(An.sunlit.final+An.shaded.final)
              
              
            }
  
  return(An)
}





Stomatal.Medlyn.equation<-function(Stress=0.9,An=10,g0=0,g1=5.75,Cs=350,ea=2535,esat=3169.217){
  
  Ds=(esat-ea)/1000 # convert from pa to kpa
  
  Gs<-g0+1.6*(1+(g1/sqrt(Ds)))*An/Cs
  
  return(Gs)
  
}




Ci.iteration_medlyn_stomatal<-function(max.iterations=20,Ci.old=241,fxCi=118.72,Ci.intermediate=3.933,Ca=360,ea=0.01264362790,esat=0.01456363309,Gb=2.416,g1=1.79,Gs.stress=9,Stress=0.7
                                       ,Vcmax=78,PPAR=500,Kc25=269,pathway="c3",Sc.o = 167,Temp=273,VmaxHA=20,GamHA=40,KcHA=50){
  count<-0
  max.iterations<-max.iterations
  Stress<-Stress
  Ci.old<-Ci.old
  Ci.bisec.guess<-Ci.old
  FxCi.old<-fxCi
  Ci.intermediate.old<-Ci.intermediate
  Gb<-Gb
  Gs.stress<-Gs.stress
  Ca<-Ca
  ea<-ea
  esat<-esat
  Vcmax=Vcmax
  PPAR=PPAR
  Kc25=Kc25
  Sc.o = Sc.o
  Temp=Temp
  VmaxHA=VmaxHA
  GamHA=GamHA
  KcHA=KcHA
  pathway=pathway
  
  repeat{
    
    Ci.new<-Ci.old-(FxCi.old/((Ci.intermediate.old-FxCi.old)/FxCi.old))
    
    
    An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = Ci.new,
                           Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
    Cs<-Ca-(1.4/Gb)*An
    Gs<-Stomatal.Medlyn.equation(Stress,An,0.01,Gs.stress,Cs,ea,esat)
    FxCi.old<-(Cs-(1.65/Gs)*An)-Ci.new # f(x)
    
    
    An.i<-Farquhar.version.3(Kc25 = Kc25,Sc.o = Sc.o,c.i = FxCi.old+Ci.new,
                             Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate f'(x)
    Cs.i<-Ca-(1.4/Gb)*An.i
    Gs.i<-Stomatal.Medlyn.equation(Stress,An.i,0.01,Gs.stress,Cs.i,ea,esat)
    Ci.i<-(Cs.i-(1.65/Gs.i)*An.i)-(FxCi.old+Ci.new)
    Ci.intermediate.old<-Ci.i
    
    Conv<-abs(Ci.new-Ci.old)
    ifelse(Conv<0.0001,break,Ci.new)
    Ci.old<-Ci.new
    
    count<-count+1
    
    ifelse(count>max.iterations,break,NA) 
  }
  ### Bisection backup method
  Conv<-Conv # break next loop if first loop is successful
  count<-0 # Reset the count
  a<-5 # min interval guess
  An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = a,
                         Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
  Cs<-Ca-(1.4/Gb)*An
  Gs<-Stomatal.Medlyn.equation(Stress,An,0.01,Gs.stress,Cs,ea,esat)
  Fx.a<-(Cs-(1.65/Gs)*An)-a
  
  b<-400 # max interval guess
  
  repeat{
    ifelse(Conv<0.0001,break,NA)
    c=(a+b)/2
    
    An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = c,
                           Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
    Cs<-Ca-(1.4/Gb)*An
    Gs<-Stomatal.Medlyn.equation(Stress,An,0.01,Gs.stress,Cs,ea,esat)
    Fx.c<-(Cs-(1.65/Gs)*An)-c
    
    ifelse(abs(Fx.c)<0.0001,break,c)
    
    a=ifelse((Fx.a*Fx.c)<0,a,c)
    b=ifelse((Fx.a*Fx.c)<0,c,400)
    Fx.a=ifelse((Fx.a*Fx.c)>0,Fx.c,Fx.a)
    
    Ci.new<-a
    count<-count+1
    ifelse(count==max.iterations,break,NA)
    
  }
  return(Ci.new)
  options(digits = 20)
}









Two.big.leaf.concept_medlyn_stomatal<-function(File.length=1,No.kinetic.sets=1,LAI=4,crop.type="Soybean",leaf.angle="Soybean",pathway="C3",Tair=25,wind=2,humidity=80/100,Precip=0.2
                                               ,SWC_1=0.3,Rn=NA,H=NA,LE=NA,LMA=2.2,pressure=100,SZA=67,SW_down=500
                                               ,Ca=360,Sensor.height=2.7,Kcat=2,Kc=10,Sc.o=101,KcHA=20,VmaxHA=30,GamHA=40){
  Count<-0
  
  An<-foreach(Kcat=Kcat,Kc=Kc,Sc.o=Sc.o,KcHA=KcHA,VmaxHA=VmaxHA,GamHA=GamHA)%:%
    foreach(Tair=Tair,LAI=LAI,wind=wind,humidity=humidity,Precip=Precip
            ,SWC_1=SWC_1,Rn=Rn,H=H,LE=LE,pressure=pressure,SZA=SZA,SW_down=SW_down
            ,.errorhandling = "remove",.combine = "rbind")%do%{
              
              ############################ 
              # MODEL INPUTS AND CONSTANTS
              
              
              
              # Constants for Leaf nitrogen taken from CLM4.5 & 5
              if(crop.type=="Wheat"|crop.type=="Sugar beet"){
                SLA<-0.07 # conversion factor from LAI to SLA
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.79
                Ci.ratio<-0.87
              }else if(crop.type=="Soybean"){
                SLA<-0.07 
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.79
                Ci.ratio<-0.87
              }else if(crop.type=="Rice"){
                Gs.max<-5.79
                SLA<-0.07
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.87
              }else if(crop.type=="Maize"|crop.type=="Test 1"|crop.type=="Test 5"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.1758
                CN<-25 
                Ci.ratio<-0.67
              }else if(crop.type=="C3"){ # any other crop type
                SLA<-0.07
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.25
                Ci.ratio<-0.87
              }else if(crop.type=="Test 2"){
                Gs.max<-1.79
                SLA<-0.08
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 3"){
                Gs.max<-1.79
                SLA<-0.09
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 4"){
                Gs.max<-1.79
                SLA<-0.1
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 6"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.16
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 7"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.15
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 8"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.14
                CN<-25
                Ci.ratio<-0.67
              }else{
                SLA<-0.07
                FLnr<-0.1758
                CN<-25
                Gs.max<-1.62
                Ci.ratio<-0.67
              }
              
              # Compute Vcmax same as CLM4.5 & CLM 5 but using Kcat
              
              Na<-(1/(SLA*CN))/14.0057  # CLM uses SLA derived from top canopy to approximate Na.aver for entire canopy
              Co=7/5500 # Houborg et al (2013) and Friend (1995) factor mole Rubisco site per mole N
              
              Vcmax<-Co*Kcat*FLnr*Na*(10^6) # Houborg et al (2013) and Friend (1995)
              
              
              # Vcmax<-Vcmax+(0.2*Vcmax) # modified Vcmax
              
              Kn<-0.3 # assumed 0.3 in CLM 5 and 4.5 for all canopy types (except forest and tropical canopies)
              
              Tair1<-273.15+Tair
              Tair.atm<-Tair1+(0.0098*Sensor.height)
              Esat.air<-Esat(Temp=Tair.atm,pressure = pressure)
              E<-(Esat.air$`Esat(Pascals)`*humidity)
              Pm<-29.2
              
              # Two stream approximation for canopy radiation
              Canopy.radiation.visible<-Two.stream.approximation.version.3(Pl=0.10,Tl=0.05,LAI=LAI,SZA=SZA,Crop.type=leaf.angle,CI=1,SW.down = SW_down,Pg=0.1)
              
              
              # Soil moisture 
              Stress<-Soil.stress.parameter(Soil.moist = (SWC_1/100),soil.type = "Silty clay loam")
              Gs.stress<-Stress*Gs.max
              
              # Heat capacity of leaf/canopy
              Cwat<-4188 # specific heat of water
              Cdry<-Cwat*0.8 
              CL<-(LMA/(1-0.8))*((1-0.8)*Cdry+0.8*Cwat)
              
              if(is.na(Rn)|is.na(H)|is.na(LE)){
                Leaf.temp.guess=Tair.atm
              }else{
                # leaf temperature
                Leaf.flux<-Rn-H-LE  
                Leaf.temp.guess<-Tair.atm+(Leaf.flux/CL) # Initial leaf temp guess
                # Leaf/canopy temp energy balance iteration
              }
              
              # Partition  kinetics across canopy
              Vcmax_modelled<-Sunlit.shaded.parameter(Param = Vcmax,Kb=Canopy.radiation.visible$Kb,LAI=LAI,Kn=Kn)
              
              
              
              # Leaf boundary layer conductance
              Dh<-18.9e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
              Dm<-13.3e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
              Pr<-Dm/Dh
              a<-(1.5*0.66*Pm*Dh*(Pr^0.33))/(Dm^0.5)
              Gb<-a*((wind/0.05)^0.5) # leaf width must be expressed in mm
              
              # Saturated vapor pressure of leaf/canopy
              Esat.leaf<-Esat(Temp=Leaf.temp.guess,pressure = pressure)
              
              # Initial intercellular CO2 guess
              Ci.guess<-Ci.ratio*Ca  # assumed ambient C02 concentrations at US-Bondville site
              
              ################## Shaded fraction
              An.shaded<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                            ,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA,Temp = Leaf.temp.guess,c.i =Ci.guess)
              Cs.shaded<-Ca-(1.4/Gb)*An.shaded
              Gs.shaded<-Stomatal.Medlyn.equation(Stress,An.shaded,0.01,Gs.stress,Cs.shaded,E,Esat.leaf$`Esat(Pascals)`)
              Ciguess.shaded<-(Cs.shaded-(1.65/Gs.shaded)*An.shaded)-Ci.guess
              An.i.shaded<-Farquhar.version.3(pathway = pathway,Kc25=Kc,Sc.o = Sc.o,c.i=Ciguess.shaded+Ci.guess,Vcmax25=Vcmax_modelled$Paramshaded,PPAR=Canopy.radiation.visible$Shaded
                                              ,Temp=Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
              Cs.i.shaded<-Ca-(1.4/Gb)*An.i.shaded
              Gs.i.shaded<-Stomatal.Medlyn.equation(Stress,An.i.shaded,0.01,Gs.stress,Cs.i.shaded,E,Esat.leaf$`Esat(Pascals)`)
              Ci.intermediate.shaded<-(Cs.i.shaded-(1.65/Gs.i.shaded)*An.i.shaded)-(Ciguess.shaded+Ci.guess)
              Ci.final.shaded<-Ci.iteration_medlyn_stomatal(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ciguess.shaded,Ci.intermediate=Ci.intermediate.shaded,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress = Stress,
                                                             Vcmax =Vcmax_modelled$Paramshaded,Temp=Leaf.temp.guess,PPAR = Canopy.radiation.visible$Shaded,Kc25=Kc,Sc.o = Sc.o,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              
              An.shaded.final<-Farquhar.version.3(pathway=pathway,c.i=Ci.final.shaded,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                                  ,Temp = Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              ################### Sunlit fraction
              An.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                            ,Temp = Leaf.temp.guess,c.i =Ci.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              Cs.sunlit<-Ca-(1.4/Gb)*An.sunlit
              Gs.sunlit<-Stomatal.Medlyn.equation(Stress,An.sunlit,0.01,Gs.stress,Cs.sunlit,E,Esat.leaf$`Esat(Pascals)`)
              Ciguess.sunlit<-(Cs.sunlit-(1.65/Gs.sunlit)*An.sunlit)-Ci.guess
              An.i.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,c.i=Ciguess.sunlit+Ci.guess,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR=Canopy.radiation.visible$Sunlit
                                              ,Temp=Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
              Cs.i.sunlit<-Ca-(1.4/Gb)*An.i.sunlit
              Gs.i.sunlit<-Stomatal.Medlyn.equation(Stress,An.i.sunlit,0.01,Gs.stress,Cs.i.sunlit,E,Esat.leaf$`Esat(Pascals)`)
              Ci.intermediate.sunlit<-(Cs.i.sunlit-(1.65/Gs.i.sunlit)*An.i.sunlit)-(Ciguess.sunlit+Ci.guess)
              Ci.final.sunlit<-Ci.iteration_medlyn_stomatal(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ciguess.sunlit,Ci.intermediate=Ci.intermediate.sunlit,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress=Stress,
                                                             Vcmax =Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit,Kc25=Kc,Sc.o = Sc.o,Temp = Leaf.temp.guess,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              
              An.sunlit.final<-Farquhar.version.3(pathway=pathway,c.i = Ci.final.sunlit,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                                  ,Temp = Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              Count<-Count+1
              
              Progress<-(Count/(File.length*No.kinetic.sets))*100
              
              print(Progress)
              
              
              as.data.frame(An.sunlit.final+An.shaded.final)
              
              
            }
  
  return(An)
}




Stomatal.leuning.equation<-function(Stress=0.9,An=10,g0=0.01,g1=9,Sc.o=100,Cs=350,ea=2535,esat=3169.217){
  Gamma=(0.5*(210/Sc.o))/0.0334
  Ds=(esat-ea)/1000
  Do=1.5
  Gs<-g0+g1*An/(Cs-Gamma)*(1+(Ds/Do))^(-1)
  
  return(Gs)
}





Ci.iteration_leuning_stomatal<-function(max.iterations=20,Ci.old=241,fxCi=118.72,Ci.intermediate=3.933,Ca=360,ea=0.01264362790,esat=0.01456363309,Gb=2.416,g1=1.79,Gs.stress=9,Stress=0.7
                                        ,Vcmax=78,PPAR=500,Kc25=269,pathway="c3",Sc.o = 167,Temp=273,VmaxHA=20,GamHA=40,KcHA=50){
  count<-0
  max.iterations<-max.iterations
  Stress<-Stress
  Ci.old<-Ci.old
  Ci.bisec.guess<-Ci.old
  FxCi.old<-fxCi
  Ci.intermediate.old<-Ci.intermediate
  Gb<-Gb
  Gs.stress<-Gs.stress
  Ca<-Ca
  ea<-ea
  esat<-esat
  Vcmax=Vcmax
  PPAR=PPAR
  Kc25=Kc25
  Sc.o = Sc.o
  Temp=Temp
  VmaxHA=VmaxHA
  GamHA=GamHA
  KcHA=KcHA
  pathway=pathway
  
  repeat{
    
    Ci.new<-Ci.old-(FxCi.old/((Ci.intermediate.old-FxCi.old)/FxCi.old))
    
    
    An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = Ci.new,
                           Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
    Cs<-Ca-(1.4/Gb)*An
    Gs<-Stomatal.leuning.equation(Stress,An,0.01,Gs.stress,Sc.o,Cs,ea,esat)
    FxCi.old<-(Cs-(1.65/Gs)*An)-Ci.new # f(x)
    
    
    An.i<-Farquhar.version.3(Kc25 = Kc25,Sc.o = Sc.o,c.i = FxCi.old+Ci.new,
                             Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate f'(x)
    Cs.i<-Ca-(1.4/Gb)*An.i
    Gs.i<-Stomatal.leuning.equation(Stress,An.i,0.01,Gs.stress,Sc.o,Cs.i,ea,esat)
    Ci.i<-(Cs.i-(1.65/Gs.i)*An.i)-(FxCi.old+Ci.new)
    Ci.intermediate.old<-Ci.i
    
    Conv<-abs(Ci.new-Ci.old)
    ifelse(Conv<0.0001,break,Ci.new)
    Ci.old<-Ci.new
    
    count<-count+1
    
    ifelse(count>max.iterations,break,NA) 
  }
  ### Bisection backup method
  Conv<-Conv # break next loop if first loop is successful
  count<-0 # Reset the count
  a<-5 # min interval guess
  An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = a,
                         Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
  Cs<-Ca-(1.4/Gb)*An
  Gs<-Stomatal.leuning.equation(Stress,An,0.01,Gs.stress,Sc.o,Cs,ea,esat)
  Fx.a<-(Cs-(1.65/Gs)*An)-a
  
  b<-400 # max interval guess
  
  repeat{
    ifelse(Conv<0.0001,break,NA)
    c=(a+b)/2
    
    An<-Farquhar.version.3(pathway=pathway,Kc25 = Kc25,Sc.o = Sc.o,c.i = c,
                           Vcmax25 = Vcmax,PPAR =PPAR,Temp = Temp,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
    Cs<-Ca-(1.4/Gb)*An
    Gs<-Stomatal.leuning.equation(Stress,An,0.01,Gs.stress,Sc.o,Cs,ea,esat)
    Fx.c<-(Cs-(1.65/Gs)*An)-c
    
    ifelse(abs(Fx.c)<0.0001,break,c)
    
    a=ifelse((Fx.a*Fx.c)<0,a,c)
    b=ifelse((Fx.a*Fx.c)<0,c,400)
    Fx.a=ifelse((Fx.a*Fx.c)>0,Fx.c,Fx.a)
    
    Ci.new<-a
    count<-count+1
    ifelse(count==max.iterations,break,NA)
    
  }
  return(Ci.new)
  options(digits = 20)
}



Two.big.leaf.concept_leuning_stomatal<-function(File.length=1,No.kinetic.sets=1,LAI=4,crop.type="Soybean",leaf.angle="Soybean",pathway="C3",Tair=25,wind=2,humidity=80/100,Precip=0.2
                                                ,SWC_1=0.3,Rn=NA,H=NA,LE=NA,LMA=2.2,pressure=100,SZA=67,SW_down=500
                                                ,Ca=360,Sensor.height=2.7,Kcat=2,Kc=10,Sc.o=101,KcHA=20,VmaxHA=30,GamHA=40){
  Count<-0
  
  An<-foreach(Kcat=Kcat,Kc=Kc,Sc.o=Sc.o,KcHA=KcHA,VmaxHA=VmaxHA,GamHA=GamHA)%:%
    foreach(Tair=Tair,LAI=LAI,wind=wind,humidity=humidity,Precip=Precip
            ,SWC_1=SWC_1,Rn=Rn,H=H,LE=LE,pressure=pressure,SZA=SZA,SW_down=SW_down
            ,.errorhandling = "remove",.combine = "rbind")%do%{
              
              ############################ 
              # MODEL INPUTS AND CONSTANTS
              
              
              
              # Constants for Leaf nitrogen taken from CLM4.5 & 5
              if(crop.type=="Wheat"|crop.type=="Sugar beet"){
                SLA<-0.07 # conversion factor from LAI to SLA
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.79
                Ci.ratio<-0.87
              }else if(crop.type=="Soybean"){
                SLA<-0.07 
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.79
                Ci.ratio<-0.87
              }else if(crop.type=="Rice"){
                Gs.max<-5.79
                SLA<-0.07
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.87
              }else if(crop.type=="Maize"|crop.type=="Test 1"|crop.type=="Test 5"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.1758
                CN<-25 
                Ci.ratio<-0.67
              }else if(crop.type=="C3"){ # any other crop type
                SLA<-0.07
                FLnr<-0.4120
                CN<-25
                Gs.max<-5.25
                Ci.ratio<-0.87
              }else if(crop.type=="Test 2"){
                Gs.max<-1.79
                SLA<-0.08
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 3"){
                Gs.max<-1.79
                SLA<-0.09
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 4"){
                Gs.max<-1.79
                SLA<-0.1
                FLnr<-0.1758
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 6"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.16
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 7"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.15
                CN<-25
                Ci.ratio<-0.67
              }else if(crop.type=="Test 8"){
                Gs.max<-1.79
                SLA<-0.07
                FLnr<-0.14
                CN<-25
                Ci.ratio<-0.67
              }else{
                SLA<-0.07
                FLnr<-0.1758
                CN<-25
                Gs.max<-1.62
                Ci.ratio<-0.67
              }
              
              # Compute Vcmax same as CLM4.5 & CLM 5 but using Kcat
              
              Na<-(1/(SLA*CN))/14.0057  # CLM uses SLA derived from top canopy to approximate Na.aver for entire canopy
              Co=7/5500 # Houborg et al (2013) and Friend (1995) factor mole Rubisco site per mole N
              
              Vcmax<-Co*Kcat*FLnr*Na*(10^6) # Houborg et al (2013) and Friend (1995)
              
              
              # Vcmax<-Vcmax+(0.2*Vcmax) # modified Vcmax
              
              Kn<-0.3 # assumed 0.3 in CLM 5 and 4.5 for all canopy types (except forest and tropical canopies)
              
              Tair1<-273.15+Tair
              Tair.atm<-Tair1+(0.0098*Sensor.height)
              Esat.air<-Esat(Temp=Tair.atm,pressure = pressure)
              E<-(Esat.air$`Esat(Pascals)`*humidity)
              Pm<-29.2
              
              # Two stream approximation for canopy radiation
              Canopy.radiation.visible<-Two.stream.approximation.version.3(Pl=0.10,Tl=0.05,LAI=LAI,SZA=SZA,Crop.type=leaf.angle,CI=1,SW.down = SW_down,Pg=0.1)
              
              
              # Soil moisture 
              Stress<-Soil.stress.parameter(Soil.moist = (SWC_1/100),soil.type = "Silty clay loam")
              Gs.stress<-Stress*Gs.max
              
              # Heat capacity of leaf/canopy
              Cwat<-4188 # specific heat of water
              Cdry<-Cwat*0.8 
              CL<-(LMA/(1-0.8))*((1-0.8)*Cdry+0.8*Cwat)
              
              if(is.na(Rn)|is.na(H)|is.na(LE)){
                Leaf.temp.guess=Tair.atm
              }else{
                # leaf temperature
                Leaf.flux<-Rn-H-LE  
                Leaf.temp.guess<-Tair.atm+(Leaf.flux/CL) # Initial leaf temp guess
                # Leaf/canopy temp energy balance iteration
              }
              
              # Partition  kinetics across canopy
              Vcmax_modelled<-Sunlit.shaded.parameter(Param = Vcmax,Kb=Canopy.radiation.visible$Kb,LAI=LAI,Kn=Kn)
              
              
              
              # Leaf boundary layer conductance
              Dh<-18.9e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
              Dm<-13.3e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
              Pr<-Dm/Dh
              a<-(1.5*0.66*Pm*Dh*(Pr^0.33))/(Dm^0.5)
              Gb<-a*((wind/0.05)^0.5) # leaf width must be expressed in mm
              
              # Saturated vapor pressure of leaf/canopy
              Esat.leaf<-Esat(Temp=Leaf.temp.guess,pressure = pressure)
              
              # Initial intercellular CO2 guess
              Ci.guess<-Ci.ratio*Ca  # assumed ambient C02 concentrations at US-Bondville site
              
              ################## Shaded fraction
              An.shaded<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                            ,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA,Temp = Leaf.temp.guess,c.i =Ci.guess)
              Cs.shaded<-Ca-(1.4/Gb)*An.shaded
              Gs.shaded<-Stomatal.leuning.equation(Stress,An.shaded,0.01,Gs.stress,Sc.o,Cs.shaded,E,Esat.leaf$`Esat(Pascals)`)
              Ciguess.shaded<-(Cs.shaded-(1.65/Gs.shaded)*An.shaded)-Ci.guess
              An.i.shaded<-Farquhar.version.3(pathway = pathway,Kc25=Kc,Sc.o = Sc.o,c.i=Ciguess.shaded+Ci.guess,Vcmax25=Vcmax_modelled$Paramshaded,PPAR=Canopy.radiation.visible$Shaded
                                              ,Temp=Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
              Cs.i.shaded<-Ca-(1.4/Gb)*An.i.shaded
              Gs.i.shaded<-Stomatal.leuning.equation(Stress,An.i.shaded,0.01,Gs.stress,Sc.o,Cs.i.shaded,E,Esat.leaf$`Esat(Pascals)`)
              Ci.intermediate.shaded<-(Cs.i.shaded-(1.65/Gs.i.shaded)*An.i.shaded)-(Ciguess.shaded+Ci.guess)
              Ci.final.shaded<-Ci.iteration_leuning_stomatal(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ciguess.shaded,Ci.intermediate=Ci.intermediate.shaded,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress = Stress,
                                                             Vcmax =Vcmax_modelled$Paramshaded,Temp=Leaf.temp.guess,PPAR = Canopy.radiation.visible$Shaded,Kc25=Kc,Sc.o = Sc.o,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              
              An.shaded.final<-Farquhar.version.3(pathway=pathway,c.i=Ci.final.shaded,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded
                                                  ,Temp = Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              ################### Sunlit fraction
              An.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                            ,Temp = Leaf.temp.guess,c.i =Ci.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              Cs.sunlit<-Ca-(1.4/Gb)*An.sunlit
              Gs.sunlit<-Stomatal.leuning.equation(Stress,An.sunlit,0.01,Gs.stress,Sc.o,Cs.sunlit,E,Esat.leaf$`Esat(Pascals)`)
              Ciguess.sunlit<-(Cs.sunlit-(1.65/Gs.sunlit)*An.sunlit)-Ci.guess
              An.i.sunlit<-Farquhar.version.3(pathway=pathway,Sc.o = Sc.o,Kc25=Kc,c.i=Ciguess.sunlit+Ci.guess,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR=Canopy.radiation.visible$Sunlit
                                              ,Temp=Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA) #intermediate ".i" values for derivative approximate
              Cs.i.sunlit<-Ca-(1.4/Gb)*An.i.sunlit
              Gs.i.sunlit<-Stomatal.leuning.equation(Stress,An.i.sunlit,0.01,Gs.stress,Sc.o,Cs.i.sunlit,E,Esat.leaf$`Esat(Pascals)`)
              Ci.intermediate.sunlit<-(Cs.i.sunlit-(1.65/Gs.i.sunlit)*An.i.sunlit)-(Ciguess.sunlit+Ci.guess)
              Ci.final.sunlit<-Ci.iteration_leuning_stomatal(max.iterations = 20,Ci.old =Ci.guess,fxCi =Ciguess.sunlit,Ci.intermediate=Ci.intermediate.sunlit,Ca=Ca,ea=E,esat=Esat.leaf$`Esat(Pascals)`,Gb=Gb,Gs.stress = Gs.stress,Stress=Stress,
                                                             Vcmax =Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit,Kc25=Kc,Sc.o = Sc.o,Temp = Leaf.temp.guess,pathway=pathway,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              
              An.sunlit.final<-Farquhar.version.3(pathway=pathway,c.i = Ci.final.sunlit,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit
                                                  ,Temp = Leaf.temp.guess,VmaxHA=VmaxHA,GamHA=GamHA,KcHA=KcHA)
              
              Count<-Count+1
              
              Progress<-(Count/(File.length*No.kinetic.sets))*100
              
              print(Progress)
              
              
              as.data.frame(An.sunlit.final+An.shaded.final)
              
              
            }
  
  return(An)
}




Yield_wheat_APSIM_lite<-function(Total_CO2=100,plant_density=100){
  
  max_yield=12166
  adjusted_totalc=Total_CO2*(200/plant_density)
  
  yield=pmin(-874.47+57.59*adjusted_totalc,12166)
  
  return(yield)

}










