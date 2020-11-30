## Wasim's Farquhar function for calculating leaf-level photosynthetic rate
Farquharforassay<-function(Jmax=100,C.i=200,Sc.o=300,kc=400,Kcat=3.2) {

# use constant O.i provided by Sharwood et al 2016
O.i<-253/1000
# use constant Rd provided by Sharwood et al 2016
Rd<-1.04
# use constant TPU provided by Sharwood et al 2016
TPU<-11.40


# Rubisco-limited rate of photosynthesis
Ac=Kcat*30*(C.i-0.5*O.i/Sc.o)/(C.i+kc)

# Light-limited rate of photosynthesis
Aj=Jmax*(C.i-0.5/Sc.o)/4*(C.i+O.i/Sc.o)

# Sink-limited rate of photosynthesis  
At=3*TPU

# Net assimilation 
A.n<-min(Ac, Aj)

return(A.n)

}

Farquhar.version.2<-function(Jmax25=100,c.i=300,Vcmax25=50,PPAR=100,Temp=273.15){
Jmax25=Jmax25
c.i=c.i
Vcmax25=Vcmax25
PPAR=PPAR

# Jmax adjusted for light conditions
  a <- 0.7
  b <- -(Jmax25 + 0.385 * PPAR)
  c <- 0.385 * Jmax25 * PPAR
  J.1 <- (-b + sqrt(b^2 - 4 * a * c) ) / (2 * a)
  J.2 <- (-b - sqrt(b^2 - 4 * a * c) ) / (2 * a)
  Jmax25 <-pmin(J.1, J.2)
  
  Kc25=404.9
  Ko25= 278.4
  Gam25=42.75
  O.i=0.21
  Kc_Ha=59.36
  Ko_Ha=35.94
  Gam_Ha=23.4
  Vmax_Ha=58.52
  Jmax_Ha=37

  
  K.c<-Kc25 * exp(Kc_Ha/298.15*8.314*(1-(298.15/Temp)))
  
  K.o<-Ko25 * exp(Ko_Ha/298.15*8.314*(1-(298.15/Temp)))
  
  gamma<-Gam25 * exp(Gam_Ha/298.15*8.314*(1-(298.15/Temp)))
  
  Vcmax<-Vcmax25 * exp(Vmax_Ha/298.15*8.314*(1-(298.15/Temp)))
  
  Jmax<-Jmax25 * exp(Jmax_Ha/298.15*8.314*(1-(298.15/Temp)))
  
  Ac <- Vcmax * (c.i - gamma) / (c.i + K.c * (1 + O.i/K.o))
  
  Aj <- Jmax * (c.i - gamma) / (4 * (c.i + 2 * gamma))
  
  A.n<-pmin(Ac,Aj)
  
  return(A.n)
}




## Ci function
Ci.calc<-function(Ca=Ca,An=An,GsT=GsT) { 
  
Ci.=Ca-(An/GsT)
  
  return(Ci.)
}

## Ball&Berry function
BallBerry<-function(An=An,Cs=Cs,rh=rH_surf,g1=g1) {
  
g0<-0.01

G<-g0+g1*((An/Cs)*rh)

Gs<- G/1.6

return(Gs)
}

## Leuning function
Leuning<-function(An=An, Cs=Cs,Ds=VPD, Gamma=C02.compensation.point, g1=3.962) {

g0<-0.01

Do<-1.5

G<-g0+g1*(An/(Cs-Gamma))*(1+Ds/Do)

Gs<- G/1.6

return(Gs)

}


Vcmaxbigleaf<-function ( Temp=Tair, GPP = "GPP", Ci, 
                      Oi = 0.21, 
                        Kc25 = 404.9, Ko25 = 278.4, Gam25 = 42.75, Kc_Ha = 79.43, 
                        Ko_Ha = 36.38, Gam_Ha = 37.83, Vcmax_Ha = 65.33, Vcmax_Hd = 200, 
                        Vcmax_dS = 0.635,
                        constants = bigleaf.constants()) 
{
  Temp <- Temp + constants$Kelvin
  Tref <- 25 + constants$Kelvin
  Kc_Ha <- Kc_Ha * constants$kJ2J
  Ko_Ha <- Ko_Ha * constants$kJ2J
  Gam_Ha <- Gam_Ha * constants$kJ2J
  Kc <- Kc25 * exp(Kc_Ha * (Temp - Tref)/(Tref * constants$Rgas * 
                                            Temp))
  Ko <- Ko25 * exp(Ko_Ha * (Temp - Tref)/(Tref * constants$Rgas * 
                                            Temp))
  Gam <- Gam25 * exp(Gam_Ha * (Temp - Tref)/(Tref * constants$Rgas * 
                                               Temp))
  Ko <- Ko * constants$J2kJ
  Vcmax <- (GPP) * (Ci + Kc * (1 + Oi/Ko))/(Ci - Gam)
  
  Vcmax25 <- Arrhenius.temp.response(Vcmax, Temp - constants$Kelvin, 
                                     Ha = Vcmax_Ha, Hd = Vcmax_Hd, dS = Vcmax_dS, constants = constants)
 
  return(Vcmax25)
}

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


ZetaRiB<-function(RiB=0.167,sh=4,Zmo=0.34,d=1){ 
  
  ZetaRiB<-ifelse(is.na(RiB),NA,ifelse(0<RiB,RiB*log((sh-d)/Zmo)/(1-5*min(RiB,0.19)),RiB*log((sh-d)/Zmo)))
  ZetaRiB1<-ifelse(0<ZetaRiB,min(2,max(ZetaRiB,0.01)),max(-100,min(ZetaRiB,-0.01))) 
  return(ZetaRiB1)
}

Psi.mRiB<-function(ZetaRiB,x.m1){
  
  Psi<-ifelse(is.na(ZetaRiB),NA,ifelse(0<ZetaRiB,-5*ZetaRiB,2*log((1+x.m1)/2)+log((1+x.m1^2)/2)-2*tan(x.m1)+pi/2))
  
  return(Psi)
}

Psi.cRiB<-function(ZetaRiB,x.c1){
  Psi<-ifelse(is.na(ZetaRiB),NA,ifelse(0<ZetaRiB,-5*ZetaRiB,2*log((1+x.c1^2)/2)))
  
  return(Psi)
}


Psi.c<-function(Zeta,x.c){
  
  Psi<-ifelse(is.na(Zeta),NA,ifelse(0<Zeta,-5*Zeta,2*log((1+x.c^2)/2)))
  return(Psi)
}

Psi.m<-function(Zeta,x.m){
  Psi<-ifelse(is.na(Zeta),NA,ifelse(0<Zeta,-5*Zeta,2*log((1+x.m)/2)+log((1+x.m^2)/2)-2*tan(x.m)+pi/2))
  return(Psi)
}


Phi.m<-function(ZetaRiB){
  Phi<-ifelse(is.na(ZetaRiB),NA,ifelse(0<ZetaRiB,1+5*ZetaRiB,(1-16*ZetaRiB)^-(1/4)))
  return(Phi)
}

Phi.c<-function(ZetaRiB){
  Phi<-ifelse(is.na(ZetaRiB),NA,ifelse(0<ZetaRiB,1+5*ZetaRiB,(1-16*ZetaRiB)^-(1/2)))
  return(Phi)
}



Stability.parameter.iteration.version.2<-function(max.iterations=10,RiBguess=5,RiBderivativeguess=162.7344,ZetaRiBguess=2,sh=3.5,d=1.4,Zco=0.004,Zmo=0.04) {
  RiB.old<-RiBguess
  ZetaRiB.old<-ZetaRiBguess
  RiBderivative.old<-RiBderivativeguess
  max.iterations<-max.iterations
  counter<-0
  
  
  repeat{
    
    RiB.new<-RiB.old-(RiB.old/RiBderivative.old)
    
    
    ZetaRiB.old<-ZetaRiB(RiB.new,sh,Zmo,d) # Approximate for zeta calculated as original equation z-d/Lmo
    ZetaRiB2c<-((ZetaRiB.old*Zco)/(sh-d)) # Approximate for Zco/Lmo
    ZetaRiB2m<-((ZetaRiB.old*Zmo)/(sh-d)) # Approximate for Zmo/Lmo
    
    x.m1=(1-16*ZetaRiB.old)^-1/4
    x.c1=(1-16*ZetaRiB.old)^-1/2
    
    x.m2=(1-16*ZetaRiB2c)^-1/4
    x.c2=(1-16*ZetaRiB2m)^-1/2
    
    Psi.mRiB1.old<-Psi.mRiB(ZetaRiB.old,x.m1)
    Psi.cRiB1.old<-Psi.cRiB(ZetaRiB.old,x.c1)
    
    Psi.cRiB2.old<-Psi.cRiB(ZetaRiB2c,x.c2)
    Psi.mRiB2.old<-Psi.mRiB(ZetaRiB2m,x.m2)
    
    RiB.old<-ZetaRiB.old*(log((sh-d)/Zco)-Psi.cRiB1.old+Psi.cRiB2.old)*(log((sh-d)/Zmo)-Psi.mRiB1.old+Psi.mRiB2.old)^-2 # RiB=ZetaFc(Zeta)Fm(Zeta)
    
    Phi.m1<-Phi.m(ZetaRiB.old)
    Phi.c1<-Phi.c(ZetaRiB.old)
    
    Phi.m2<-Phi.m(ZetaRiB2m)
    Phi.c2<-Phi.c(ZetaRiB2c)
    
    Fc<-(log((sh-d)/Zco)-Psi.cRiB1.old+Psi.cRiB2.old)
    Fm<-(log((sh-d)/Zmo)-Psi.mRiB1.old+Psi.mRiB2.old)^-2
    
    RiBderivative.old<-1+(Phi.c1-Phi.c2)*Fc-2*(Phi.m1-Phi.m2)*Fm
    
    
    # Exit test
    counter<-counter+1
    
    ifelse(counter==max.iterations,break,NA)
    
  }
  return(RiB.new) 
}


eulers.function.soil.temp<-function(max.iterations=10,Tn.central=-0.3882,Tn.plus=-0.3804,Tn.minus=-0.396,time=0.1,space=-0.020,k1=2.92,Cv=29999){
  Tn.plus.old<-Tn.plus
  Tn.central.old<-Tn.central
  Tn.minus.old<-Tn.minus
  t<-time
  s<-space
  s.central<-space/2
  count<-0
  k1<-k1
  Cv<-Cv
  repeat{
    
    
    Tn.new<-Tn.central.old+((k1/Cv)*(t/s^2))*(Tn.minus.old-2*Tn.central.old+Tn.plus.old)
    
    Tn.central.old<--(k1/Cv)*(Tn.new/(s.central^2)) # current time point
    Tn.plus.old<-Tn.central.old+Tn.central.old*s     # Spatial forward approx
    Tn.minus.old<-Tn.central.old-Tn.central.old*s        # Spatial backward approx
    
    
    count<-count+1
    
    Conver<-abs(((k1/Cv)*(t/s^2)))
    
    ifelse(Conver>0.5,break,NA)
    
    ifelse(count==max.iterations,break,NA)
    
    
  }
  return(Tn.new)
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
 

Seasonal.LAI.GSI.estimate<-function(Lmax=5,NEE=300,NEEmax=500,NEEmin=50){
  # GSI stands for growing season index
  
  GSI<-((NEE-NEEmin)/(NEEmax-NEEmin))
  LAI<-GSI*Lmax  
  
  return(LAI)
}





Seasonal.NDVI.estimate<-function(PPAR.down=100,PPAR.up=300,SW.up=200,SW.down=300){
  Vis=0.45
  NIR=0.55
  TotalSW=SW.down-SW.up
  
  Rvis=PPAR.up/PPAR.down # calculate reflectance in visible band
  
  VISin=Vis*TotalSW # calculate reflectance in NIR band
  NIRin=NIR*TotalSW
  VISout=Rvis*VISin
  NIRout=SW.up-VISout
  Rnir=NIRout/NIRin
  
  NDVI=(Rnir-Rvis)/(Rnir+Rvis)
  
  return(NDVI)
}



Seasonal.LAI.estimate<-function(NDVI=0.6,NDVImax=0.9,NDVImin=0.2,k=2.40){
  
  NDVI.normalised=((NDVImax-NDVI)/(NDVImax-NDVImin)) # Winston
  
  LAI<--k*log(NDVI.normalised)
  
  return(LAI)
}



Two.stream.approximation.version.2<-function(Pl=0.53,Tl=0.34,XL=-0.30,LAI=4,SZA=27,Leaf.type="Plagiophile",CI=1,SW.down=500,Pg=0.1){
  Pg=Pg
  CI=CI # clumping index
  Pl=Pl
  Tl=Tl
  XL=XL
  SZA=SZA*pi/180 # convert to radians
  
  
  
  
  # leaf angles evaluated at each nine inclination angles/adjusted for probability of inclination
  if(Leaf.type=="Planophile"){   # Goel and Strebel (1984)
    Inclination.angle=0.4670501
  }else if(Leaf.type=="Erectophile"){
    Inclination.angle=1.103746 
  }else if(Leaf.type=="Plagiophile"){
    Inclination.angle=0.7853982 
  }else if(Leaf.type=="Uniform"){
    Inclination.angle=0.7853982
  }else{# spherical
    Inclination.angle=1.000074
  } 
  
  
  # Frequency of leaves in each leaf angle class
  if(Leaf.type=="Planophile"){
    Fi=0.472
  }else if(Leaf.type=="Erectophile"){
    Fi=0.707
  }else if(Leaf.type=="Plagiophile") {
    Fi=0.706
  }else if(Leaf.type=="Uniform"){
    Fi=1
  }else{
    Fi=0.1
  }
  
  
  # Direct beam coefficient calculated assuming absence of azimuth preference
  Zenithdiff=1.570796-SZA
  c=((sin(Inclination.angle)^2)-(cos(SZA)^2))^(1/2)
  a=cos(SZA)*cos(Inclination.angle)
  b=sin(SZA)*sin(Inclination.angle)
  a.b=a/b
  G.z<-ifelse(Inclination.angle<Zenithdiff,cos(SZA)*cos(Inclination.angle),(2/pi)*(c+a*(asin(a.b))))
  G<-G.z*Fi
  Kb<-G/cos(SZA)
  
  # Diffuse beam coefficient for 30 degree incriment sky hemisphere
  Zenithdiff=1.570796-0.5235988
  c=((sin(Inclination.angle)^2)-(cos(0.5235988)^2))^(1/2)
  a=cos(0.5235988)*cos(Inclination.angle)
  b=sin(0.5235988)*sin(Inclination.angle)
  a.b=a/b
  G.z30<-ifelse(Inclination.angle<Zenithdiff,cos(0.5235988)*cos(Inclination.angle),(2/pi)*(c+a*(asin(a.b))))
  G30<-G.z30*Fi
  Kb30<-G30/cos(0.5235988)
  Td30<-2*(exp(-(Kb30*LAI))*sin(0.5235988)*cos(0.5235988)*0.1745)
  Kd30<--((log(Td30))/LAI)
  
  Zenithdiff=1.570796-1.047198
  c=((sin(Inclination.angle)^2)-(cos(1.047198)^2))^(1/2)
  a=cos(1.047198)*cos(Inclination.angle)
  b=sin(1.047198)*sin(Inclination.angle)
  a.b=a/b
  G.z60<-ifelse(Inclination.angle<Zenithdiff,cos(1.047198)*cos(Inclination.angle),(2/pi)*(c+a*(asin(a.b))))
  G60<-G.z60*Fi
  Kb60<-G60/cos(1.047198)
  Td60<-2*exp(-Kb60*LAI)*sin(1.047198)*cos(1.047198)*0.1745
  Kd60<--(log(Td60)/LAI)
  
  
  Zenithdiff=1.570796-1.570796
  c=((sin(Inclination.angle)^2)-(cos(1.570796)^2))^(1/2)
  a=cos(1.570796)*cos(Inclination.angle)
  b=sin(1.570796)*sin(Inclination.angle)
  a.b=a/b
  G.z90<-ifelse(Inclination.angle<Zenithdiff,cos(1.570796)*cos(Inclination.angle),(2/pi)*(c+a*(asin(a.b))))
  G90<-G.z90*Fi
  Kb90<-G90/cos(1.570796)
  Td90<-2*(exp(-(Kb90*LAI))*sin(1.570796)*cos(1.570796)*0.1745)
  Kd90<--((log(Td90))/LAI)
  
  Kd<-Kd30+Kd60+Kd90
  Kd<-ifelse(Kd>1,1,Kd)
  
  Td=Td90+Td60+Td30 # transmission for solar down diffuse estimaste
  
  # Two stream approximation by Liou et al (2002)
  W=Pl+Tl
  O1=0.5-0.633*XL-0.33*(XL^2)
  O2=0.877*(1-2*O1)
  as=(W/2)*((1/Kd)/((1/Kd)+(cos(SZA)*O2)))
  Bo=(Kb+Kd/Kb)*as
  B=(1/2)*(Pl+Tl+(Pl-Tl)*(cos(Inclination.angle)^2)) # Portion of light scattered in the backward (upward) for diffuse beam
  
  
  b=(1-(1-B)*W)*Kd
  c=B*W*Kd
  h.a=(b^2)-(c^2)
  h=sqrt(h.a)
  u=(h-b-c)/(2*h)
  v=(h+b+c)/(2*h)
  y1=(Bo*Kb-b*Bo-c*(1-Bo))*((W*Kb*SW.down)/((h^2)-(Kb^2)))
  y2=((1-Bo)*Kb+c*Bo+b*(1-Bo))*((W*Kb*SW.down)/((h^2)-(Kb^2)))
  
  # Direct beam flux
  
  n2=((v*(y1+y2*Pg+Pg*SW.down)*exp(-Kb*CI*LAI))-(y2*(u+v*Pg)*exp(-h*CI*LAI)))/((v*(v+u*Pg)*exp(h*CI*LAI))-(u*(u+v*Pg)*exp(-h*CI*LAI)))
  n1=((y2-n2*u)/v)
  
  
  SW.up.direct.above=(-y1*exp(-Kb*CI))+(n1*u*exp(-h*CI))+(n2*v*exp(h*CI))
  SW.down.direct.below=(y2*exp(-Kb*CI*LAI))-(n1*v*exp(-h*CI*LAI))-(n2*u*exp(h*CI*LAI))
  SW.up.direct.LAI=(-y1*exp(-Kb*CI*LAI))+(n1*u*exp(-h*CI*LAI))+(n2*v*exp(h*CI*LAI))
  SW.canopy.b=((1-exp(-Kb*CI*LAI))*SW.down)-(SW.up.direct.above)+(SW.up.direct.LAI-SW.down.direct.below)
  
  # Diffuse beam flux
  SW.diffuse=SW.down*Td
  n2.diffuse=(SW.diffuse*(u+v*Pg)*exp(-h*CI*LAI))/((v*(v+u*Pg)*exp(-h*CI*LAI))-(u*(u+v*Pg)*exp(-h*CI*LAI)))
  n1.diffuse=((-SW.diffuse-n2*u)/v)
  
  SW.up.diffuse.above=(n1.diffuse*u*exp(-h*CI))+(n2.diffuse*v*exp(h*CI))
  SW.down.diffuse.below=(-n1.diffuse*v*(exp(-h*CI*LAI)))-(n2.diffuse*u*exp(h*CI*LAI))
  SW.up.diffuse.LAI=(n1.diffuse*u*exp(-h*CI*LAI))+(n2.diffuse*v*exp(h*CI*LAI))
  SW.canopy.d=SW.diffuse-SW.up.diffuse.above+SW.up.diffuse.LAI-SW.down.diffuse.below
  
  
  Total.canopy<-SW.canopy.d+SW.canopy.b
  Total.canopy<-as.data.frame(Total.canopy)
  PPAR.canopy<-Total.canopy*4.6
  
  # Incorporate sunlit and shaded fractions
  Solar_up<-SW.up.diffuse.above+SW.up.direct.above
  Solar_up<-as.data.frame(Solar_up)
  Sunlit=PPAR.canopy*exp(-Kb*LAI)
  Sunlit<-as.data.frame(Sunlit)
  colnames(Sunlit)<-"Sunlit"
  Shaded=PPAR.canopy*(1-exp(-Kb*LAI))
  Shaded<-as.data.frame(Shaded)
  colnames(Shaded)<-"Shaded"
  Kb<-as.data.frame(Kb)
  colnames(Kb)<-"Kb"
  
  Data<-cbind(Kb,PPAR.canopy,Sunlit,Shaded,Solar_up)
  
  return(Data)
}




Air.density<-function(Tair=285,pressure=100){
  g=8.314
  pressure=pressure*1000
  
  Pm<-pressure/(Tair*g)
  
  return(Pm)
}


General.euler.forward<-function(max.iterations=25,y.old=54.666,t=0.01){
  Count=0
  y.old=y.old
  t=t
  max.iterations=max.iterations
  
  repeat{
    y.new<-y.old+y.old*t
    
    Conv<-((y.new-y.old)/y.old)*100
    
    y.old<-y.new
    
    Count<-Count+1
    
    ifelse(Count==max.iterations,break,NA)
    
    options(digits = 5)
    
  }
  return(y.new)
}

General.euler.backward<-function(max.iterations=25,y.old=54.666,t=0.01){
  Count=0
  y.old=y.old
  t=t
  
  max.iterations=max.iterations
  
  repeat{
    y.new<-y.old-y.old*t
    
    y.old<-y.new
  
    Count<-Count+1
    
    ifelse(Count==max.iterations,break,NA)
    
    options(digits = 5)
    
  }
  return(y.new)
}


Farquhar.version.3<-function(pathway="C3",Kc25=10.8,Sc.o=99.9,c.i=300,Vcmax25=20,KcHA=63.4
                             ,VmaxHA=41.5,GamHA=23.6,PPAR=700,Temp=298.15){
  
  Kc25=Kc25/0.0334
  Ko25=278.4
  O.i=210
  Gam25=(0.5*(O.i/Sc.o))/0.0334
  Kc_Ha=KcHA
  Ko_Ha=35.94
  Gam_Ha=GamHA
  Vmax_Ha=VmaxHA
  S_Vmax=0.66839
  Vmax_Hd=200

  ### C3 temperature adjustments
  K.c<-Kc25 * exp(Kc_Ha/298.15*8.314*(1-(298.15/Temp))) 
  
  K.o<-Ko25 * exp(Ko_Ha/298.15*8.314*(1-(298.15/Temp)))
  
  gamma<-Gam25 * exp(Gam_Ha/298.15*8.314*(1-(298.15/Temp)))
  
  
  if(Temp>298.15){
    Vcmax<-Vcmax25 * exp(Vmax_Ha/298.15*8.314*(1-(298.15/Temp))) * (1+exp((298.15*S_Vmax-Vmax_Hd)/(298.15*0.001*8.314)))/(1+exp((S_Vmax*Temp-Vmax_Hd)/(0.001*8.314*Temp)))
  }else{
    Vcmax<-Vcmax25 * exp(Vmax_Ha/298.15*8.314*(1-(298.15/Temp)))
  }
  
  Jmax<-1.67*Vcmax
  
  # Jmax adjusted for light conditions
  IPSII<- 0.5*0.85*PPAR
  a <- 0.7
  b <- -(IPSII+Jmax)
  c <-IPSII*Jmax
  J.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  J.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  Jmax <-pmin(J.1, J.2)
  
  if(pathway=="C3"){

  Ac <- Vcmax * (c.i - gamma) / (c.i + K.c * (1 + O.i/K.o))
  
  Aj <- Jmax * (c.i - gamma) / (4 * (c.i + 2 * gamma))
  
  a=0.98
  b=-(Ac+Aj)
  c=Ac*Aj
  A.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A <-pmin(A.1, A.2)
  
  Rd=0.015*Vcmax
  
  A.n<-A-Rd
 
  
  }else{
    ### C4 temperature adjustments
 

  Vcmax=Vcmax25*(exp(Vmax_Ha/298.15*8.314*(1-(298.15/Temp)))/((1+exp(0.3*(Temp-313.15)))*(1+exp(0.2*(288.15-Temp)))))
   
  Rd=Vcmax*0.025
   
  Kp=20000*Vcmax
  
  Ap<-Kp*c.i
  Aj<- 0.05*PPAR
  Ac<-Vcmax
  
  a=0.80
  b=-(Ac+Aj)
  c=Ac*Aj
  A.i.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A.i.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  Ai <-pmin(A.i.1, A.i.2) # solution given by the smaller root
  
  a=0.95
  b=-(Ai+Ap)
  c=Ai*Ap
  A.1 <- (-b + sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A.2 <- (-b - sqrt((b^2) - 4 * a * c) ) / (2 * a)
  A <-pmin(A.1, A.2)
  
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
    
  
  Bw<-ifelse(Soil.moist>Soil.moist.sat,1,ifelse(Soil.moist>0.01 & Soil.moist<0.477,((Soil.moist-MinBw)/(Soil.moist.sat-MinBw)),0))
  
  return(Bw)
}




Moist.predictor.correct.method<-function(G.water=0.2){
  
  soil.type="Silty clay loam"
  max.increments.time=4
  t=500
  Count.1=0
  

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
  
  
  
  if(soil.type=="Silty clay loam"){
    Matric.sat<--35.6*2 
  }else if(soil.type=="Sand"){
    Matric.sat<--12.1*2
  }else if(soil.type=="Loamy sand"){
    Matric.sat<--9*2
  }else if(soil.type=="Sandy loam"){
    Matric.sat<--21.8*2
  }else if(soil.type=="Silt loam"){
    Matric.sat<--78.6*2
  }else if(soil.type=="Loam"){
    Matric.sat<--47.8*2
  }else if(soil.type=="Sandy clay loam"){
    Matric.sat<--29.9*2
  }else if(soil.type=="Clay loam"){
    Matric.sat<--63*2
  }else if(soil.type=="Sandy clay"){
    Matric.sat<--15.3*2
  }else if(soil.type=="Silty clay"){
    Matric.sat<--49*2
  }else{
    Matric.sat<--40.5*2
  } # Clay
  
  
  
  if(soil.type=="Silty clay loam"){
    K.sat<-0.61*2 
  }else if(soil.type=="Sand"){
    K.sat<-63.36*2
  }else if(soil.type=="Loamy sand"){
    K.sat<-56.28*2
  }else if(soil.type=="Sandy loam"){
    K.sat<-12.48*2
  }else if(soil.type=="Silt loam"){
    K.sat<-2.59*2
  }else if(soil.type=="Loam"){
    K.sat<-2.50*2
  }else if(soil.type=="Sandy clay loam"){
    K.sat<-2.27*2
  }else if(soil.type=="Clay loam"){
    K.sat<-0.88*2
  }else if(soil.type=="Sandy clay"){
    K.sat<-0.78*2
  }else if(soil.type=="Silty clay"){
    K.sat<-0.37*2
  }else {
    K.sat<-0.46*2
    
  } # Clay
  
  
  
  if(soil.type=="Silty clay loam"){
    b<-7.75 
  }else if(soil.type=="Sand"){
    b<-4.05
  }else if(soil.type=="Loamy sand"){
    b<-4.38
  }else if(soil.type=="Sandy loam"){
    b<-4.90
  }else if(soil.type=="Silt loam"){
    b<-5.30
  }else if(soil.type=="Loam"){
    b<-5.39
  }else if(soil.type=="Sandy clay loam"){
    b<-7.12
  }else if(soil.type=="Clay loam"){
    b<-8.52
  }else if(soil.type=="Sandy clay"){
    b<-10.4
  }else if(soil.type=="Silty clay"){
    b<-10.4
  }else{
    b<-11.4
  } # Clay
  
  Boundary<-ifelse(G.water==0|is.na(G.water),0.01,Matric.sat*(G.water/Soil.moist.sat)^(-b))
  Wilting.matric.1<-ifelse(G.water==0|is.na(G.water),-29990.59,Matric.sat*(G.water/Soil.moist.sat)^(-b))
  Wilting.matric.2<-ifelse(G.water==0|is.na(G.water),-29990.59,Matric.sat*(G.water/Soil.moist.sat)^(-b))
  Wilting.matric.3<-ifelse(G.water==0|is.na(G.water),-29990.59,Matric.sat*(G.water/Soil.moist.sat)^(-b))
  Wilting.matric.4<-ifelse(G.water==0|is.na(G.water),-29990.59,Matric.sat*(G.water/Soil.moist.sat)^(-b))
  Wilting.matric.5<-ifelse(G.water==0|is.na(G.water),-29990.59,Matric.sat*(G.water/Soil.moist.sat)^(-b))

  
  

  repeat{
    
    Moist.1<-Soil.moist.sat*(Wilting.matric.1/Matric.sat)^(-1/b)
    Moist.2<-Soil.moist.sat*(Wilting.matric.2/Matric.sat)^(-1/b)
    Moist.3<-Soil.moist.sat*(Wilting.matric.3/Matric.sat)^(-1/b)
    Moist.4<-Soil.moist.sat*(Wilting.matric.4/Matric.sat)^(-1/b)
    Moist.5<-Soil.moist.sat*(Wilting.matric.5/Matric.sat)^(-1/b)
    
    K1<-K.sat*(Moist.1/Soil.moist.sat)^(2*b+3)
    K2<-K.sat*(Moist.2/Soil.moist.sat)^(2*b+3)
    K3<-K.sat*(Moist.3/Soil.moist.sat)^(2*b+3)
    K4<-K.sat*(Moist.4/Soil.moist.sat)^(2*b+3)
    K5<-K.sat*(Moist.5/Soil.moist.sat)^(2*b+3)
    
    
    K1.half<-0.5*(K1) # layer 1 effective K is half K.cond of layer 1
    K1.for.half<-0.5*(K1+K2) # check that your halfing is correct!!
    K2.for.half<-0.5*(K2+K3) 
    K3.for.half<-0.5*(K3+K4)# k1.for/2 gives centre of next layer and k1 for gives the adjacent centre layer
    K4.for.half<-0.5*(K4+K5)
    K5.for.half<-0.5*(K5) # special case for last layer
    
    K2.back.half<-0.5*(K1+K2)
    K3.back.half<-0.5*(K2+K3)
    K4.back.half<-0.5*(K4+K3)
    K5.back.half<-0.5*(K5+K4)
    
    Ci.n.1<-(-Soil.moist.sat/(b*Matric.sat))*((Wilting.matric.1/Matric.sat)^(-1/(b-1)))
    Ci.n.2<-(-Soil.moist.sat/(b*Matric.sat))*((Wilting.matric.2/Matric.sat)^(-1/(b-1)))# Specific moist capacity assumed uniform for each layer
    Ci.n.3<-(-Soil.moist.sat/(b*Matric.sat))*((Wilting.matric.3/Matric.sat)^(-1/(b-1)))
    Ci.n.4<-(-Soil.moist.sat/(b*Matric.sat))*((Wilting.matric.4/Matric.sat)^(-1/(b-1)))
    Ci.n.5<-(-Soil.moist.sat/(b*Matric.sat))*((Wilting.matric.5/Matric.sat)^(-1/(b-1)))
    
    
    A1.ex<-0                # tridiagonal coefficients for explicit half time step
    A2.ex<--K2.back.half/2
    A3.ex<--K3.back.half/2
    A4.ex<--K4.back.half/2
    A5.ex<--K5.back.half/2
    
    A<-as.vector(c(A2.ex,A3.ex,A4.ex,A5.ex))
    
    
    B1.ex<-(Ci.n.1*2/(t/2))+(K1.half/2)+(K1.for.half/2)
    B2.ex<-(Ci.n.2*2/(t/2))+(K2.back.half/2)+(K2.for.half/2)
    B3.ex<-(Ci.n.3*2/(t/2))+(K3.back.half/2)+(K3.for.half/2)
    B4.ex<-(Ci.n.4*2/(t/2))+(K4.back.half/2)+(K4.for.half/2)
    B5.ex<-(Ci.n.5*2/(t/2))+(K5.back.half/2)
    
    B<-as.vector(c(B1.ex,B2.ex,B3.ex,B4.ex,B5.ex))
    
    C1.ex<--K1.for.half/2
    C2.ex<--K2.for.half/2
    C3.ex<--K3.for.half/2
    C4.ex<--K4.for.half/2
    C5.ex<-0   
    
    C<-as.vector(c(C1.ex,C2.ex,C3.ex,C4.ex))
    
    D1.ex<-Ci.n.1*2/(t/2)*Wilting.matric.1+(K1.half/2)*Boundary+K1.half-K1.for.half
    D2.ex<-Ci.n.2*2/(t/2)*Wilting.matric.2+K2.back.half-K2.for.half
    D3.ex<-Ci.n.3*2/(t/2)*Wilting.matric.3+K3.back.half-K3.for.half
    D4.ex<-Ci.n.4*2/(t/2)*Wilting.matric.4+K4.back.half-K4.for.half
    D5.ex<-Ci.n.5*2/(t/2)*Wilting.matric.5+K5.back.half-K5.for.half
    
    D<-as.vector(c(D1.ex,D2.ex,D3.ex,D4.ex,D5.ex))
    
    Matric.t.half<-trisolve(B,C,A,rhs = D)
    
    Moist.t.half<-Soil.moist.sat*(Matric.t.half/Matric.sat)^(-1/b)
    

    Ci.t.half<-(-Soil.moist.sat/(b*Matric.sat))*(Matric.t.half/Matric.sat)^(-1/(b-1))
    
    K.t.half<-K.sat*(Moist.t.half/Soil.moist.sat)^(2*b+3)
    
    ######## Now solve over a full time step using Crank-Nicolson method
    
    K1.half<-0.5*(K.t.half[1]) # layer 1 effective K is half K.cond of layer 1
    K1.for.half<-0.5*(K.t.half[1]+K.t.half[2])
    K2.for.half<-0.5*(K.t.half[2]+K.t.half[3])
    K3.for.half<-0.5*(K.t.half[3]+K.t.half[4]) 
    K4.for.half<-0.5*(K.t.half[4]+K.t.half[5])
    K5.for.half<-0.5*(K.t.half[5])
    
    K2.back.half<-0.5*(K.t.half[1]+K.t.half[2])
    K3.back.half<-0.5*(K.t.half[2]+K.t.half[3]) 
    K4.back.half<-0.5*(K.t.half[4]+K.t.half[3])
    K5.back.half<-0.5*(K.t.half[5]+K.t.half[4])
    
    
    A1<-0               # tridiagonal coefficients for Crank-Nicolson 
    A2<--K2.back.half/(2*2)
    A3<--K3.back.half/(2*2)
    A4<--K4.back.half/(2*2)
    A5<--K5.back.half/(2*2)
    
    A<-c(A2,A3,A4,A5)
    
    
    B1<-(Ci.t.half[1]*2/t)+(K1.half/(2*2))+(K1.for.half/(2*2))
    B2<-(Ci.t.half[2]*2/t)+(K2.back.half/(2*2))+(K2.for.half/(2*2))
    B3<-(Ci.t.half[3]*2/t)+(K3.back.half/(2*2))+(K3.for.half/(2*2))
    B4<-(Ci.t.half[4]*2/t)+(K4.back.half/(2*2))+(K4.for.half/(2*2))
    B5<-(Ci.t.half[5]*2/t)+(K5.back.half/(2*2))
    
    B<-c(B1,B2,B3,B4,B5)
    
    
    C1<--K1.for.half/(2*2)
    C2<--K2.for.half/(2*2)
    C3<--K3.for.half/(2*2)
    C4<--K4.for.half/(2*2)
    C5<-0  
    
    C<-c(C1,C2,C3,C4)
    
    D1<-(Ci.t.half[1]*2/t)*Wilting.matric.1+(K1.half/(2*2))*(Boundary-Wilting.matric.1)-(K1.for.half/(2*2))*(Wilting.matric.1-Wilting.matric.2)+K1.half-K1.for.half
    D2<-(Ci.t.half[2]*2/t)*Wilting.matric.2+(K2.back.half/(2*2))*(Wilting.matric.1-Wilting.matric.2)-(K2.for.half/(2*2))*(Wilting.matric.2-Wilting.matric.3)+K2.back.half-K2.for.half # UP TO HERE WITH CHECKING
    D3<-(Ci.t.half[3]*2/t)*Wilting.matric.3+(K3.back.half/(2*2))*(Wilting.matric.2-Wilting.matric.3)-(K3.for.half/(2*2))*(Wilting.matric.3-Wilting.matric.4)+K3.back.half-K3.for.half
    D4<-(Ci.t.half[4]*2/t)*Wilting.matric.4+(K4.back.half/(2*2))*(Wilting.matric.3-Wilting.matric.4)-(K4.for.half/(2*2))*(Wilting.matric.4-Wilting.matric.5)+K4.back.half-K4.for.half
    D5<-(Ci.t.half[5]*2/t)*Wilting.matric.5+(K5.back.half/(2*2))*(Wilting.matric.4-Wilting.matric.5)-(K5.for.half/(2*2))*(Wilting.matric.5-Wilting.matric.5)+K5.back.half-K5.for.half
    
    D<-c(D1,D2,D3,D4,D5)
    
    # Matric/2 is the initial condition assuming already half of the boundary condition
    
    
    Matric.full.t<-trisolve(B,C,A,rhs = D)
    
    Final.moist<-Soil.moist.sat*((Matric.full.t/Matric.sat)^(-(1/b)))
    
    Wilting.matric.1<-Matric.full.t[1]
    Wilting.matric.2<-Matric.full.t[2]
    Wilting.matric.3<-Matric.full.t[3]
    Wilting.matric.4<-Matric.full.t[4]
    Wilting.matric.5<-Matric.full.t[5]
    
    Count.1=Count.1+1
    
    ifelse(Final.moist[1]>=0.3,break,NA)
    
    ifelse(max.increments.time==Count.1,break,NA)
    
  }
  
  
  return(Final.moist[2])
  
}


  
  
  Soil.moist.iteration<-function(max.iterations=5,Moist.guess=0.4,Matric.guess=400,soil.type="Silty clay loam"){
    Matric.old<-Matric.guess
    Moist.old<-Moist.guess
    max.iterations=max.iterations
    soil.type<-soil.type
    Count<-0
    
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
    
    
    
    if(soil.type=="Silty clay loam"){
      Matric.sat<-35.6 
    }else if(soil.type=="Sand"){
      Matric.sat<-12.1
    }else if(soil.type=="Loamy sand"){
      Matric.sat<-9
    }else if(soil.type=="Sandy loam"){
      Matric.sat<-21.8
    }else if(soil.type=="Silt loam"){
      Matric.sat<-78.6
    }else if(soil.type=="Loam"){
      Matric.sat<-47.8
    }else if(soil.type=="Sandy clay loam"){
      Matric.sat<-29.9
    }else if(soil.type=="Clay loam"){
      Matric.sat<-63
    }else if(soil.type=="Sandy clay"){
      Matric.sat<-15.3
    }else if(soil.type=="Silty clay"){
      Matric.sat<-49
    }else{
      Matric.sat<-40.5
      
    } # Clay
    
    
    
    if(soil.type=="Silty clay loam"){
      K.sat<-0.61 
    }else if(soil.type=="Sand"){
      K.sat<-63.36
    }else if(soil.type=="Loamy sand"){
      K.sat<-56.28
    }else if(soil.type=="Sandy loam"){
      K.sat<-12.48
    }else if(soil.type=="Silt loam"){
      K.sat<-2.59
    }else if(soil.type=="Loam"){
      K.sat<-2.50
    }else if(soil.type=="Sandy clay loam"){
      K.sat<-2.27
    }else if(soil.type=="Clay loam"){
      K.sat<-0.88
    }else if(soil.type=="Sandy clay"){
      K.sat<-0.78
    }else if(soil.type=="Silty clay"){
      K.sat<-0.37
    }else{
      K.sat<-0.46
      
    } # Clay
    
    
    
    if(soil.type=="Silty clay loam"){
      b<-7.75 
    }else if(soil.type=="Sand"){
      b<-4.05
    }else if(soil.type=="Loamy sand"){
      b<-4.38
    }else if(soil.type=="Sandy loam"){
      b<-4.90
    }else if(soil.type=="Silt loam"){
      b<-5.30
    }else if(soil.type=="Loam"){
      b<-5.39
    }else if(soil.type=="Sandy clay loam"){
      b<-7.12
    }else if(soil.type=="Clay loam"){
      b<-8.52
    }else if(soil.type=="Sandy clay"){
      b<-10.4
    }else if(soil.type=="Silty clay"){
      b<-10.4
    }else{
      b<-11.4
      } # Clay
    
    
    repeat{
      Moist.new=Moist.old-((Matric.old/Matric.sat)^(-1/b))/(-Soil.moist.sat/b*Matric.sat*(Matric.old/Matric.sat)^(-1/b-1))
      
      Count<-Count+1
      
      ifelse(Count==max.iterations,break,NA)
      
      Matric.old<-(Moist.new/Soil.moist.sat)^(1/(-1/b))*Matric.sat
      
      Moist.old<-Moist.new
        
      
      
      print(Moist.new)
    }
    return(Moist.new)
  }
 
 
  

  K.conductivity.iteration<-function(max.iterations=2,Moist.guess=0.4,K.guess=0.001,soil.type="Silty clay loam"){
    Moist.old<-Moist.guess
    K.old<-K.guess
    Count<-0
    max.iterations=max.iterations
    
    if(soil.type=="Silty clay loam"){
      Soil.moist.sat<-0.477*10e5
    }
    else if(soil.type=="Sand"){
      Soil.moist.sat<-0.395*10e5
    }
    else if(soil.type=="Loamy sand"){
      Soil.moist.sat<-0.410*10e5
    }
    else if(soil.type=="Sandy loam"){
      Soil.moist.sat<-0.435*10e5
    }
    else if(soil.type=="Silt loam"){
      Soil.moist.sat<-0.485*10e5
    }
    else if(soil.type=="Loam"){
      Soil.moist.sat<-0.451*10e5
    }
    else if(soil.type=="Sandy clay loam"){
      Soil.moist.sat<-0.420*10e5
    }
    else if(soil.type=="Clay loam"){
      Soil.moist.sat<-0.476*10e5
    }
    else if(soil.type=="Sandy clay"){
      Soil.moist.sat<-0.426*10e5
    }
    else if(soil.type=="Silty clay"){
      Soil.moist.sat<-0.492*10e5
    }
    else{
      Soil.moist.sat<-0.482*10e5
      
    } # Clay
    
    
    
    if(soil.type=="Silty clay loam"){
      Matric.sat<-35.6 
    }
    else if(soil.type=="Sand"){
      Matric.sat<-12.1
    }
    else if(soil.type=="Loamy sand"){
      Matric.sat<-9
    }
    else if(soil.type=="Sandy loam"){
      Matric.sat<-21.8
    }
    else if(soil.type=="Silt loam"){
      Matric.sat<-78.6
    }
    else if(soil.type=="Loam"){
      Matric.sat<-47.8
    }
    else if(soil.type=="Sandy clay loam"){
      Matric.sat<-29.9
    }
    else if(soil.type=="Clay loam"){
      Matric.sat<-63
    }
    else if(soil.type=="Sandy clay"){
      Matric.sat<-15.3
    }
    else if(soil.type=="Silty clay"){
      Matric.sat<-49
    }
    else{
      Matric.sat<-40.5
      
    } # Clay
    
    
    
    if(soil.type=="Silty clay loam"){
      K.sat<-0.61 
    }
    else if(soil.type=="Sand"){
      K.sat<-63.36
    }
    else if(soil.type=="Loamy sand"){
      K.sat<-56.28
    }
    else if(soil.type=="Sandy loam"){
      K.sat<-12.48
    }
    else if(soil.type=="Silt loam"){
      K.sat<-2.59
    }
    else if(soil.type=="Loam"){
      K.sat<-2.50
    }
    else if(soil.type=="Sandy clay loam"){
      K.sat<-2.27
    }
    else if(soil.type=="Clay loam"){
      K.sat<-0.88
    }
    else if(soil.type=="Sandy clay"){
      K.sat<-0.78
    }
    else if(soil.type=="Silty clay"){
      K.sat<-0.37
    }
    else{
      K.sat<-0.46
      
    } # Clay
    
    
    
    if(soil.type=="Silty clay loam"){
      b<-7.75 
    }
    else if(soil.type=="Sand"){
      b<-4.05
    }
    else if(soil.type=="Loamy sand"){
      b<-4.38
    }
    else if(soil.type=="Sandy loam"){
      b<-4.90
    }
    else if(soil.type=="Silt loam"){
      b<-5.30
    }
    else if(soil.type=="Loam"){
      b<-5.39
    }
    else if(soil.type=="Sandy clay loam"){
      b<-7.12
    }
    else if(soil.type=="Clay loam"){
      b<-8.52
    }
    else if(soil.type=="Sandy clay"){
      b<-10.4
    }
    else if(soil.type=="Silty clay"){
      b<-10.4
    }
    else{
      b<-11.4
      
    } # Clay
    
    
    
    repeat{
      K.new=K.old-(K.sat*(Moist.old/Soil.moist.sat)^(2*b+3))/(K.sat*(2*b+3)/Soil.moist.sat*(Moist.old/Soil.moist.sat)^(2*b+2))
      Count<-Count+1
      Moist.old<-(K.new/K.sat)^(1/(2*b+3))*Soil.moist.sat
      
      ifelse(Count==max.iterations,break,NA)
      
      print(K.new)
    }
    return(K.new)
  }

  
  
  Two.big.leaf.concept<-function(File.length=1,No.kinetic.sets=1,LAI=4,crop.type="Soybean",leaf.angle="Soybean",pathway="C3",Tair=25,wind=2,humidity=80/100,Precip=0.2
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
                if(crop.type=="Wheat"){
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
  
  Stats.mean<-function(Final.An=Final.An){
    Data<-as.data.frame(Final.An)
    Data[Data<0]<-NA
    Mean<-Data %>% summarise_each(funs(mean(.,na.rm=T)))
    Meancol<-data.frame(Mean=unlist(Mean, use.names = FALSE))
    SD<-Data %>% summarise_each(funs(sd(.,na.rm=T)))
    SDcol<-data.frame(SD=unlist(SD, use.names = FALSE))
    Data<-cbind(Meancol,SDcol)
    
    return(Data)
  }
  
  Stats.sum<-function(Final.An=Final.An){
    Data<-as.data.frame(Final.An)
    Data[Data<0]<-NA
    Sum<-Data %>% summarise_each(funs(sum(.,na.rm=T)))
    Sumcol<-data.frame(Sum=unlist(Sum, use.names = FALSE))
    return(Sumcol)
  }
  
  Stats.RMSE.MAE.R2<-function(Predict=Predict.An,True=An,No.obs=328){
    Predict<-as.data.frame(unlist(Predict))
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
  
  
  
  Two.big.leaf.concept.Xin.LAI.method<-function(File.length=1,No.kinetic.sets=1,LAI=LAI,SLN=SLN,m.ratio=0.137,Tair=25,wind=7,humidity=80/100,VPD=0.6,Precip=0.2
                                 ,SWC_1=0.3,Rn=NA,H=NA,LE=NA,LMA=2.2,pressure=100,SZA=67,SW_down=300,PPFD_down=600
                                 ,NEE.or.GGP=27,NEEmax=30,NEEmin=1,Ca=360,leaf.type="Spherical",LAI.max=6
                                 ,Sensor.height=10,Kcat=1.5*30,Kc=9/0.0334,Sc.o=67){
    Count<-0
    
    foreach(Vcmax=Vcmax,Kc=Kc,Sc.o=Sc.o)%:%
      foreach(Tair=Tair,LAI=LAI,SLN=SLN,wind=wind,humidity=humidity,Precip=Precip
              ,SWC_1=SWC_1,Rn=Rn,H=H,LE=LE,pressure=pressure,SZA=SZA,SW_down=SW_down
              ,NEE.or.GGP=NEE.or.GGP,.combine="rbind",.errorhandling = "remove")%dopar%{
                
                ############################
                # MODEL INPUTS AND CONSTANTS
                
                
                # Allocation of nitrogen throughout the canopy: USE PFTs
                
              
                Jmax<-Vcmax*1.7 # assumed 1.7*Vcmax[v]
                Tair1<-273.15+Tair
                Tair.atm<-Tair1+(0.0098*Sensor.height)
                Esat.air<-Esat(Temp=Tair.atm,pressure = pressure)
                E<-(Esat.air$`Esat(Pascals)`*humidity)
                q.atm<-E/(pressure*1000)
                Pm<-29.2
                
                # Two stream approximation for canopy radiation
                Canopy.radiation.visible<-Two.stream.approximation.version.2(Pl=0.11,Tl=0.05,XL=-0.30,LAI=LAI,SZA=SZA,Leaf.type=leaf.type,CI=1,SW.down = SW_down,Pg=0.1)
                Canopy.radiation.NIR<-Two.stream.approximation.version.2(Pl=0.35,Tl=0.34,XL=-0.30,LAI=LAI,SZA=SZA,Leaf.type=leaf.type,CI=1,SW.down = SW_down,Pg=0.2)
                
                # Soil moisture 
                Fi<-1*tanh(LAI) # Fraction of leaf that collects water, default value of "1"
                WInter<-Fi*Precip
                WThru<-WInter*(1-Fi)
                G.water<-(WThru)+(SWC_1/100) # add Soil water content from Data.fileset as "initial conditions"
                Soil.moist<-Moist.predictor.correct.method(G.water = G.water)
                Stress<-Soil.stress.parameter(Soil.moist = Soil.moist,soil.type = "Silty clay loam")
                Gs.stress<-Stress*9
                
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
                Vcmax_modelled<-Sunlit.shaded.parameter(Param = Vcmax,Kb=Canopy.radiation.visible$Kb,LAI=LAI)
                Jmax_modelled<-Sunlit.shaded.parameter(Param = Jmax,Kb=Canopy.radiation.visible$Kb,LAI=LAI)
                
                
                # Leaf boundary layer conductance
                Dh<-18.9e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
                Dm<-13.3e-06*(101325/(pressure*1000))*((Leaf.temp.guess/273.15)^1.81)
                Pr<-Dm/Dh
                a<-(1.5*0.66*Pm*Dh*(Pr^0.33))/(Dm^0.5)
                Gb<-a*((wind/5)^0.5)
                
                # Saturated vapor pressure of leaf/canopy
                Esat.leaf<-Esat(Temp=Leaf.temp.guess,pressure = pressure)
                
                # Initial intercellular CO2 guess
                Ci.guess<-0.67*Ca  # assumed ambient C02 concentrations at US-Bondville site
                
                ################## Shaded fraction
                An.shaded<-Farquhar.version.3(Jmax25 =Vcmax_modelled$Paramshaded*1.7,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded,Temp = Leaf.temp.guess,c.i =Ci.guess)
                Cs.shaded<-Ca-(1.4/Gb)*An.shaded
                Gs.shaded<-Stomatal.quadratic.equation(Stress,Gb,An.shaded,0.01,Gs.stress,Cs.shaded,q.atm,Esat.leaf$qsat)
                Ciguess.shaded<-(Cs.shaded-(1.65/Gs.shaded)*An.shaded)-Ci.guess
                An.i.shaded<-Farquhar.version.3(Jmax25=Vcmax_modelled$Paramshaded*1.7,Kc25=Kc,Sc.o = Sc.o,c.i=Ciguess.shaded+Ci.guess,Vcmax25=Vcmax_modelled$Paramshaded,PPAR=Canopy.radiation.visible$Shaded,Temp=Leaf.temp.guess) #intermediate ".i" values for derivative approximate
                Cs.i.shaded<-Ca-(1.4/Gb)*An.i.shaded
                Gs.i.shaded<-Stomatal.quadratic.equation(Stress,Gb,An.i.shaded,0.01,Gs.stress,Cs.i.shaded,q.atm,Esat.leaf$qsat)
                Ci.intermediate.shaded<-(Cs.i.shaded-(1.65/Gs.i.shaded)*An.i.shaded)-(Ciguess.shaded+Ci.guess)
                Ci.final.shaded<-Ci.iteration(max.iterations = 50,Ci.old =Ci.guess,fxCi =Ciguess.shaded,Ci.intermediate=Ci.intermediate.shaded,Ca=Ca,ea=q.atm,esat=Esat.leaf$qsat,Gb=Gb,Gs.stress = Gs.stress,Stress = Stress,
                                              Vcmax =Vcmax_modelled$Paramshaded,Jmax = Vcmax_modelled$Paramshaded*1.7,PPAR = Canopy.radiation.visible$Shaded,Kc25=Kc,Sc.o = Sc.o)
                
                
                An.shaded.final<-Farquhar.version.3(c.i=Ci.final.shaded,Kc25=Kc,Sc.o = Sc.o,Jmax25 = Vcmax_modelled$Paramshaded*1.7,Vcmax25 = Vcmax_modelled$Paramshaded,PPAR = Canopy.radiation.visible$Shaded,Temp = Leaf.temp.guess)
                
                
                ################### Sunlit fraction
                An.sunlit<-Farquhar.version.3(Jmax25 =Vcmax_modelled$Paramsunlit*1.7,Sc.o = Sc.o,Kc25=Kc,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit,Temp = Leaf.temp.guess,c.i =Ci.guess)
                Cs.sunlit<-Ca-(1.4/Gb)*An.sunlit
                Gs.sunlit<-Stomatal.quadratic.equation(Stress,Gb,An.sunlit,0.01,Gs.stress,Cs.sunlit,q.atm,Esat.leaf$qsat)
                Ciguess.sunlit<-(Cs.sunlit-(1.65/Gs.sunlit)*An.sunlit)-Ci.guess
                An.i.sunlit<-Farquhar.version.3(Jmax25=Vcmax_modelled$Paramsunlit*1.7,Sc.o = Sc.o,Kc25=Kc,c.i=Ciguess.sunlit+Ci.guess,Vcmax25=Vcmax_modelled$Paramsunlit,PPAR=Canopy.radiation.visible$Sunlit,Temp=Leaf.temp.guess) #intermediate ".i" values for derivative approximate
                Cs.i.sunlit<-Ca-(1.4/Gb)*An.i.sunlit
                Gs.i.sunlit<-Stomatal.quadratic.equation(Stress,Gb,An.i.sunlit,0.01,Gs.stress,Cs.i.sunlit,q.atm,Esat.leaf$qsat)
                Ci.intermediate.sunlit<-(Cs.i.sunlit-(1.65/Gs.i.sunlit)*An.i.sunlit)-(Ciguess.sunlit+Ci.guess)
                Ci.final.sunlit<-Ci.iteration(max.iterations = 50,Ci.old =Ci.guess,fxCi =Ciguess.sunlit,Ci.intermediate=Ci.intermediate.sunlit,Ca=Ca,ea=q.atm,esat=Esat.leaf$qsat,Gb=Gb,Gs.stress = Gs.stress,Stress=Stress,
                                              Vcmax =Vcmax_modelled$Paramsunlit,Jmax = Vcmax_modelled$Paramsunlit*1.7,PPAR = Canopy.radiation.visible$Sunlit,Kc25=Kc,Sc.o = Sc.o,Temp = Leaf.temp.guess)
                
                An.sunlit.final<-Farquhar.version.3(c.i = Ci.final.sunlit,Jmax25 = Vcmax_modelled$Paramsunlit*1.7,Kc25=Kc,Sc.o = Sc.o,Vcmax25 = Vcmax_modelled$Paramsunlit,PPAR = Canopy.radiation.visible$Sunlit,Temp = Leaf.temp.guess)
                
                Final.An<-as.data.frame(An.sunlit.final+An.shaded.final)
                
                Count<-Count+1
                
                Progress<-(Count/(File.length*No.kinetic.sets))*100
                print(Progress)
                
                return(Final.An)
              }
  }
  
  Stats<-function(Final.An=Final.An){
    Data<-as.data.frame(Final.An)
    Data[Data<0]<-NA
    Mean<-Data %>% summarise_each(funs(mean(.,na.rm=T)))
    Meancol<-data.frame(Mean=unlist(Mean, use.names = FALSE))
    SD<-Data %>% summarise_each(funs(sd(.,na.rm=T)))
    SDcol<-data.frame(SD=unlist(SD, use.names = FALSE))
    Data<-cbind(Meancol,SDcol)
    
    return(Data)
  }
  
  
  LAI.seasonal.Xin.method<-function(LAIs.guess=0.01,GPP.reported="Yes"){
    
    if(GPP.reported=="Yes")
    {
      LAI<-(lag(LAIs.guess)+LAIs.guess)/2
    }  
  return(LAI)
    }
  

ft.MODIS<-function(Current=1.5){
  Tmin<--8.00
  Tmax<-12.02
  
  f<-pmax(pmin(((Current-Tmin)/(Tmax-Tmin)),1,na.rm=T),0)
  return(f)
}

fVPD.MODIS<-function(Current=1.5){
  VPDmin<-6.5
  VPDmax<-43
  
  f<-pmax(pmin((1-((Current-VPDmin)/(VPDmax-VPDmin))),1,na.rm=T),0)
  
  return(f)
}

Fixed.point.LAI.method.Xin.<-function(Kb=0.5,LAIs.guess=0.01
                             ,ft=1,fVPD=1,PPFD_down = 224
                             ,LUE=0.0330,m.ratio =0.21503315){
  
LAIs.new<-LAIs.guess

  repeat{
  

  FPARs.new<-1-exp((-Kb*LAIs.new))
  GPPs.new<-PPFD_down*FPARs.new*LUE*ft*fVPD
  LAIs.old<-LAIs.new
  LAIs.new<-m.ratio*GPPs.new
  print(LAIs.new)
  
  Conv<-abs(LAIs.new-LAIs.old)
  ifelse(Conv<0.000001,break,NA)
  
  }
  return(LAIs.new)
}

Gammastar<-function(Sc.o=100){
  Gammastar<-0.5*(210/Sc.o)
return(Gammastar)
}

HA.estimate<-function(T1=35,T2=25,V1=5.6,V2=3.2){
  T1=T1+273.15
  T2=T2+273.15
  HA<-log((V1/V2))*8.314/((1/T2)-(1/T1))
  return(HA)
}

Two.stream.approximation.version.3<-function(Pl=0.10,Tl=0.05,LAI=6,SZA=30
                                             ,Crop.type="Spherical",CI=1,SW.down=500,Pg=0.1){
  Pg=Pg
  CI=CI # clumping index
  Pl=Pl
  Tl=Tl
  SZA=SZA*pi/180 # convert to radians
  L=LAI
  
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
  O1=0.5-0.633*XL-0.33*(XL^2)
  O2=0.877*(1-2*O1)
  G=O1+O2*cos(SZA)
  Kb<-G/cos(SZA)
  
  # Diffuse beam coefficient for 30 degree incriment sky hemisphere
  G30=O1+O2*cos(0.5235988)
  Kb30<-G30/cos(0.5235988)
  Td30<-2*(exp(-(Kb30*L))*sin(0.5235988)*cos(0.5235988)*0.1745)

  G60=O1+O2*cos(1.047198)
  Kb60<-G60/cos(1.047198)
  Td60<-2*exp(-Kb60*L)*sin(1.047198)*cos(1.047198)*0.1745
  
  G90=O1+O2*cos(1.570796)
  Kb90<-G90/cos(1.570796)
  Td90<-2*(exp(-(Kb90*L))*sin(1.570796)*cos(1.570796)*0.1745)
  
  if(is.infinite(Td90)|is.na(Td90)){
    Td=Td30+Td60
  }else if(is.infinite(Td60)|is.na(Td60)&is.infinite(Td90)|is.na(Td90)){
    Td=Td30
  }else if(is.infinite(Td60)|is.na(Td60)&is.infinite(Td90)|is.na(Td90)&is.infinite(Td30)|is.na(Td30)){
  Td=0  
  }else{
    Td=Td30+Td60+Td90 
  }
  Kd<--(log(Td)/L)
  
  # Two stream approximation by Liou et al (2002)
  W=Pl+Tl
  Cos.angle=(1+XL/2)^2
  B=((1/2)*(Pl+Tl+(Pl-Tl)*Cos.angle))/W # Fraction of light scattered in the backward (upward) for diffuse beam
  Bo=((Kb+Kd/Kb)*((W/2)*G/G+cos(SZA)*O2*(1-cos(SZA)*O1/G+cos(SZA)*O2*log(G+cos(SZA)*O1+cos(SZA)*O2/cos(SZA)*O1))))/W ##### FINISH!!!!!!
  
  # Partition Direct and diffuse solar radiation
  SW.diffuse=SW.down*Td
  SW.down=SW.down-SW.diffuse
  
  # Coefficients
  b=(1-(1-B)*W)*Kd
  c=B*W*Kd
  h.a=(b^2)-(c^2)
  h=sqrt(h.a)
  u=(h-b-c)/(2*h)
  v=(h+b+c)/(2*h)
  y1=(Bo*Kb-b*Bo-c*(1-Bo))*((W*Kb*SW.down)/((h^2)-(Kb^2)))
  y2=((1-Bo)*Kb+c*Bo+b*(1-Bo))*((W*Kb*SW.down)/((h^2)-(Kb^2)))

  # Direct beam flux
  n2=((v*(y1+y2*Pg+Pg*SW.down)*exp(-Kb*CI*L))-(y2*(u+v*Pg)*exp(-h*CI*L)))/((v*(v+u*Pg)*exp(h*CI*L))-(u*(u+v*Pg)*exp(-h*CI*L)))
  n1=((y2-n2*u)/v)
  
  SW.up.direct.above=(-y1*exp(-Kb*CI*L))+(n1*u*exp(-h*CI*L))+(n2*v*exp(h*CI*L))
  SW.down.direct.below=(y2*exp(-Kb*CI*L))-(n1*v*exp(-h*CI*L))-(n2*u*exp(h*CI*L))
  SW.up.direct.L=(-y1*exp(-Kb*CI*L))+(n1*u*exp(-h*CI*L))+(n2*v*exp(h*CI*L))
  SW.canopy.b=((1-exp(-Kb*CI*L))*SW.down)-(SW.up.direct.above)+(SW.up.direct.L-SW.down.direct.below) # give as per leaf area
  
  # Diffuse beam flux
  n2.diffuse=(SW.diffuse*(u+v*Pg)*exp(-h*CI*L))/((v*(v+u*Pg)*exp(-h*CI*L))-(u*(u+v*Pg)*exp(-h*CI*L)))
  n1.diffuse=((-SW.diffuse-n2.diffuse*u)/v)
  
  SW.up.diffuse.above=(n1.diffuse*u*exp(-h*CI*L))+(n2.diffuse*v*exp(h*CI*L))
  SW.down.diffuse.below=(-n1.diffuse*v*(exp(-h*CI*L)))-(n2.diffuse*u*exp(h*CI*L))
  SW.up.diffuse.L=(n1.diffuse*u*exp(-h*CI*L))+(n2.diffuse*v*exp(h*CI*L))
  SW.canopy.d=SW.diffuse-SW.up.diffuse.above+SW.up.diffuse.L-SW.down.diffuse.below # give as per leaf area
  
  # Calculate per leaf area
  SW.leaf.d=SW.up.diffuse.above-SW.down.diffuse.below
  SW.leaf.b=SW.up.direct.above-SW.down.direct.below
  
  
  # Calculate sunlit/shaded fractios of LAI
  
  Lsun=(1-exp(-Kb*L))/Kb
  Lshad=L-Lsun
  
  # Partition PAR over Lsun and Lshad
  
  #Shaded
  Scatter.b.rad=((W*Kb*SW.down*exp(-Kb*L)))+SW.leaf.b
  Shaded=(SW.leaf.d+Scatter.b.rad)

  #Sunlit
  Sunlit=Shaded+((1-W)*Kb*SW.down)
  
  #Scale to canopy
  Shaded<-as.data.frame(Shaded*4.6*Lshad) #Scale to canopy and conver to umol photons 
  colnames(Shaded)<-"Shaded"
  Sunlit<-as.data.frame(Sunlit*4.6*Lsun) #Scale to canopy and conver to umol photons 
  colnames(Sunlit)<-"Sunlit"
  Total.canopy<-as.data.frame(Sunlit+Shaded)
  colnames(Total.canopy)<-"Total"
  
  Kb<-as.data.frame(Kb)
  colnames(Kb)<-"Kb"
  
  Data<-cbind(Kb,Total.canopy,Sunlit,Shaded)
  
  return(Data)
}

light.response<-function (data, NEE = "NEE", Reco = "Reco", PPFD = "PPFD", PPFD_ref = 2000, 
          ...) 
{
  
  mod <- nls(NEE ~ alpha * PPFD/(1 - (PPFD/PPFD_ref) + (alpha *PPFD/GPP_ref)) - Reco, start = list(alpha = 0.01, GPP_ref = 5), 
  control = nls.control(minFactor = 1/10096))
  return(mod)
}


