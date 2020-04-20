library(lubridate)
library(raster)
library(rworldmap)
library(ncdf4)
library(ggplot2)
library(reshape2)
library(tidyverse)

orders=read.csv("data/stayathomerestrictions.csv")
orders$Start.Date=as.Date(orders$Start.Date,format="%m/%d/%Y")
orders$doy=yday(orders$Start.Date)

#read in harvest and planting dates
datecrops=list.files("data/crop calendars/filled")

data(countriesCoarse)
r0=raster(nrow=2160,ncol=4320)
extent(r0)=c(-180,180,-90,90)
#merge in start date of orders
countriesCoarse=merge(countriesCoarse,orders[,3:4],by.x="ISO3",by.y="iso3",all.x=T)
countryraster=rasterize(countriesCoarse,r0,field=countriesCoarse$doy)

window=60 #number of days from start of stay at home order to be considered at risk

windowfun=function(r1,r2,date,w=window){
  ifelse(r1>date&r1<(date+w),1,ifelse(r2>date&r2<(date+w),1,0))
}

# windowfun_justplanting=function(r1,r2,date,w=window){
#   ifelse(r1>date&r1<(date+w),1,0)
# }

for(i in 1:length(datecrops)){
  plantdate=raster(paste0("data/crop calendars/filled/",datecrops[i]),varname="plant")
  harvestdate=raster(paste0("data/crop calendars/filled/",datecrops[i]),varname="harvest")
  atrisktemp=overlay(plantdate,harvestdate,countryraster,fun=function(x,y,z) windowfun(x,y,z))
  if(i==1) atrisk=atrisktemp
  if(i>1) atrisk=stack(atrisk,atrisktemp)
  print(i)
}

# names(atrisk)=datecrops;save(atrisk,file="data/productionatrisk_justplanting.Rdat")

#identify double cropping areas for barley, maize, oats, rice, sorghum, and wheat to allocate production across two growing seasons
doubles=list(c(1,2),c(5,6),c(8,9),c(13,14),c(16,17),c(22,23))
doublecropareas=function(r1,r2){
  ifelse(is.na(r1)|is.na(r2),NA,1)
}
datecropsdouble=list.files("data/crop calendars/unfilled")

for(i in 1:length(doubles)){
  season_small=raster(paste0("data/crop calendars/unfilled/",datecropsdouble[doubles[[i]][1]]),varname="plant")
  season_main=raster(paste0("data/crop calendars/unfilled/",datecropsdouble[doubles[[i]][2]]),varname="plant")
  doubletemp=overlay(season_small,season_main,fun=function(x,y) doublecropareas(x,y))
  if(i==1) doublecrop=doubletemp; if(i>1) doublecrop=stack(doublecrop,doubletemp)
}

doublesplit=0.3 #fraction production assumed to come in smaller season for double cropping areas

#area layers to match up to calendar layers
pulses=c("bambara","bean","broadbean","chickpea","cowpea","lentil","pigeonpea","pulsenes","vetch")

crops=c("barley","cassava","groundnut","maize","millet","oats","potato","pulses","rapeseed","rice","rye","sorghum","soybean","sugarbeet","sunflower","sweetpotato","wheat","yam")
names(doublecrop)=c("barley","maize","oats","rice","sorghum","wheat")
names(atrisk)=c("barley2","barley","cassava","groundnut","maize2","maize","millet","oats2","oats","potato","pulses","rapeseed","rice2","rice","rye","sorghum2","sorghum","soybean","sugarbeet","sunflower","sweetpotato","wheat2","wheat","yam")

doubles=c(1,4,6,10,12,17)

globalprod=data.frame(crop=crops,total_prod=rep(0,length(crops)),atrisk_prod=rep(0,length(crops)))

"%notin%"=Negate("%in%")

doublecropadjustment=function(r1,r2,r3,adj=doublesplit){
  #r1 is main season at risk, r2 is small season at risk, r3 is map of double cropped areas
  if(is.na(r3)){
    #no double cropping
    return(ifelse(r2==0,r1,ifelse(r1==0,r2,r1)))
  }
  if(r3==1){
    #with double cropping
    return(r1*(1-doublesplit)+r2*doublesplit)
  }
}

for(i in 1:length(crops)){
  crop=crops[i]
  if(crop=="pulses"){
    #for pulses need to loop through production of 9 different pulses and sum total
    cropprod=raster(paste0("data/crop areas/",pulses[1],"_HarvAreaYield_Geotiff/",pulses[1],"_Production.tif"))
    for(j in 2:length(pulses)){
      cropprod=cropprod+raster(paste0("data/crop areas/",pulses[j],"_HarvAreaYield_Geotiff/",pulses[j],"_Production.tif"))
    }
  }
  if(crop!="pulses") cropprod=raster(paste0("data/crop areas/",crop,"_HarvAreaYield_Geotiff/",crop,"_Production.tif"))
  
  globalprod$total_prod[i]=cellStats(cropprod,stat="sum")
  if(i==1) cropprodmap=cropprod; if(i>1) cropprodmap=stack(cropprodmap,cropprod)

  if(i%notin%doubles){
    atriskprodtemp=cropprod*atrisk[[crop]]
    globalprod$atrisk_prod[i]=cellStats(atriskprodtemp,stat="sum")
    if(i==1) atriskprod=atriskprodtemp; if(i>1) atriskprod=stack(atriskprod,atriskprodtemp)
  }
  
  if(i%in%doubles){
    atrisk_smallseason=cropprod*atrisk[[paste0(crop,2)]]
    atrisk_mainseason=cropprod*atrisk[[crop]]
    #adjust for double cropping
    atriskprodtemp=overlay(atrisk_mainseason,atrisk_smallseason,doublecrop[[crop]],fun=function(x,y,z) doublecropadjustment(x,y,z,adj=doublesplit))
    globalprod$atrisk_prod[i]=cellStats(atriskprodtemp,stat="sum")
    if(i==1) atriskprod=atriskprodtemp; if(i>1) atriskprod=stack(atriskprod,atriskprodtemp)
  }
}
globalprod$frac_atrisk=globalprod$atrisk_prod/globalprod$total_prod*100
write.csv(globalprod,"data/globalproductionatrisk.csv")

#check for potential nationally-important production short-falls

#aggreagate at-risk production to the national level to estimate national production short-falls
nationalprod=raster::extract(cropprodmap,countriesCoarse[which(countriesCoarse$ISO3%in%orders$iso3),],fun="sum",na.rm=TRUE)
atrisknationalprod=raster::extract(atriskprod,countriesCoarse[which(countriesCoarse$ISO3%in%orders$iso3),],fun="sum",na.rm=TRUE)
colnames(nationalprod)=crops; colnames(atrisknationalprod)=crops
nationalprod=as.data.frame(nationalprod);atrisknationalprod=as.data.frame(atrisknationalprod)
nationalprod$iso3=countriesCoarse$ISO3[which(countriesCoarse$ISO3%in%orders$iso3)]; atrisknationalprod$iso3=countriesCoarse$ISO3[which(countriesCoarse$ISO3%in%orders$iso3)]
nationalprod=melt(nationalprod,id.vars="iso3")
atrisknationalprod=melt(atrisknationalprod,id.vars="iso3")
colnames(nationalprod)=c("iso3","crop","nationalprod");colnames(atrisknationalprod)=c("iso3","crop","atriskprod")
nationalprod=merge(nationalprod,atrisknationalprod)
nationalprod$fractionatrisk=nationalprod$atriskprod/nationalprod$nationalprod*100

#aggregate to country level, weigting by caloric importance
fbs=read.csv("data/fao_foodbalancesheets.csv")
cropcrosswalk=read.csv("data/faocropcrosswalk.csv")
fbs=merge(fbs,cropcrosswalk[,2:3])
countrycrosswalk=read.csv("data/faoiso3crosswalk.csv")
fbs=merge(fbs,countrycrosswalk[,c(1,5)],by.x="Area.Code",by.y="FAOSTAT")

#crop weights by country
calweights=fbs%>%
  select(ISO3,crop,Value,Element)%>%
  group_by(ISO3,crop)%>%
  summarize(croptotalcal=sum(Value[Element=="Food supply (kcal/capita/day)"]),croptotalprod=sum(Value[Element=="Production"]),croptotalnetimports=sum((Value[Element=="Import Quantity"]-Value[Element=="Export Quantity"])))%>%
  group_by(ISO3)%>%
  mutate(calweights=croptotalcal/sum(croptotalcal),prodweights=croptotalprod/sum(croptotalprod+croptotalnetimports))%>%
  mutate(weights=calweights*prodweights/(sum(calweights*prodweights)))
  
nationalprod=merge(nationalprod,calweights[,c(1,2,8)],by.x=c("iso3","crop"),by.y=c("ISO3","crop"))
foodatrisk=nationalprod%>%
  group_by(iso3)%>%
  summarize(atriskfood=sum(fractionatrisk*weights,na.rm=T)/sum(weights))

#calculate value of production per agricultural worker
#2016 value ag production
agvalue=read.csv("data/fao_valueagproduction.csv")
agvalue=merge(agvalue,countrycrosswalk[,c(1,5)],by.x="Area.Code",by.y="FAOSTAT")
agvalue=agvalue%>%select("ISO3","Value")
colnames(agvalue)=c("iso3","valueagprod")

atriskiso=c("ARM","BGD","BOL","IND","IRQ","MAR","NPL","PAK","RWA")

temp=nationalprod%>%
  filter(iso3%in%atriskiso)%>%
  mutate(bigcrops=fractionatrisk*weights)%>%
  arrange(iso3,desc(bigcrops))

#2016 population
pop=read.csv("data/population_2016.csv")
pop$population=as.numeric(as.character(pop$population))
agvalue=merge(agvalue,pop)

#share labor force in agriculture
aglaborshare=read.csv("data/laborshare_ag.csv")
aglaborshare=aglaborshare%>%
  filter(Year==2016)%>%
  select("Code","X...of.total.employment.")
colnames(aglaborshare)=c("iso3","aglaborshare")
agvalue=merge(agvalue,aglaborshare)

agvalue=agvalue%>%
  mutate(agvalueperlabor=valueagprod/(population*aglaborshare/100))

foodatrisk=merge(foodatrisk,agvalue[,c(1,5)],all.x=TRUE,all.y=FALSE)
write.csv(foodatrisk,file="data/nationalfoodatrisk.csv")

#global stocks as fraction of consumption graph
psd=read.csv("data/globalcerealstocks.csv")
a=ggplot(psd,aes(x=Year,y=Frac_Cons*100,group=Crop,col=Crop))+theme_bw()+geom_line(lwd=1.3)
a=a+labs(y="Stocks as Fraction of Consumption (%)")+scale_color_manual(values=c("#1a2d32","#b0ac4b","#d49228"),labels=c("Maize","Rice","Wheat"))

