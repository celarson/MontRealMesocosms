#Analysis of spiny water flea eDNA in mesocosm experiments
#Datasets generated from physical enumeration and qPCR analysis

#########
#libraries
library(ggplot2)
library(plyr)
library(lme4)
library(sjPlot)
library(MASS)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)

######################
#Color vectors
yearcolvec<-c("deepskyblue4","#d95f02")
tankcolvec<-c("#E41A1C","#A24057","#606692","#3A85A8","#42977E", "#4AAA54", "#629363", "#7E6E85", "#9C509B",
              "#C4625D", "#EB751F", "#FF9709", "#FFC81D", "black", "#E1C62F", "#BF862B", "grey60", "#CC6A6F",
              "#EB7AA9", "#E086B5", "green", "grey30")

#################################
#Upload and standardize datasets

#Environmental datasets
bythoenv<-read.csv("EnvVarsBytho.csv")
#Create wide dataset based on day
bythoenvwide<-reshape(bythoenv, idvar = c("Tankno","Year"), timevar = "Day", direction = "wide")

#physical
bytho2122<-read.csv("2122PhysicalBythoResults.csv")
bytho2122$Barb<-as.factor(bytho2122$Barb)
bytho2122$BarbNoIntro<-paste(bytho2122$Barb,bytho2122$NoIntro)
#Add IMO threshold factor
bytho2122$IMOLiving<-bytho2122$TotalLive
bytho2122$IMOLiving[bytho2122$TotalLive>=10] <- "Greater"
bytho2122$IMOLiving[bytho2122$TotalLive<10&bytho2122$TotalLive>0] <- "Less"
bytho2122$IMOLiving[bytho2122$TotalLive<1] <- "Zero"
bytho2122$IMOTotal<-bytho2122$TotalBythoIndiv
bytho2122$IMOTotal[bytho2122$TotalBythoIndiv>=10] <- "Greater"
bytho2122$IMOTotal[bytho2122$TotalBythoIndiv<10&bytho2122$TotalBythoIndiv>0] <- "Less"
bytho2122$IMOTotal[bytho2122$TotalBythoIndiv<1] <- "Zero"

#Change order of levels
bytho2122$BarbNoIntro<- factor(bytho2122$BarbNoIntro, levels = c("NA 0", "2barb 1", "3barb 1", "2barb 5", "3barb 5", "3barb 10"))
bytho2122$BarbNoIntro<-revalue(bytho2122$BarbNoIntro, c("NA 0"="None", "2barb 1"="1, 2nd instar", "3barb 1"="1, 3rd instar",
                                                        "2barb 5"="5, 2nd instars", "3barb 5"="5, 3rd instars", 
                                                        "3barb 10"="10, 3rd instars"))
#create ordered barbnointro variable
bytho2122$BarbNoIntroOrdered<-as.numeric(revalue(bytho2122$BarbNoIntro, c("None"=0, "1, 2nd instar"=1, "1, 3rd instar"=2, "5, 2nd instars"=3,
                                                               "5, 3rd instars"=4,"10, 3rd instars"=5)))
#Create "tank and year" variable
bytho2122$tankyear<-paste(bytho2122$Tankno,bytho2122$Year)
#qPCR
bythoqPCR<-read.csv("qPCRBytho.csv")

#standardize qPCR dataset
#Remove faulty Ct value 
bythoqPCR$Ct[bythoqPCR$Ct<13] <- NA
#Use conversion factor and standard curve generated from gblock run 9.26.23
bythoqPCR$CtConverted<-bythoqPCR$Ct-bythoqPCR$Conversion
bythoqPCR$DNACopiesperuLncconv<-(10^(11-(0.26*(bythoqPCR$CtConverted))))
bythoqPCR$DNACopiesperuLnc<-(10^(11-(0.26*(bythoqPCR$Ct))))

#Standardize based on dilution factor of extracted DNA and volume filtered
bythoqPCR$DNACopiesperuLconv<-bythoqPCR$DNACopiesperuLncconv*bythoqPCR$DilutionFactor
bythoqPCR$DNACopiesperuL<-bythoqPCR$DNACopiesperuLnc*bythoqPCR$DilutionFactor
bythoqPCR$DNACopiesperuLLconv<-bythoqPCR$DNACopiesperuLconv/bythoqPCR$VolumeFilteredL
bythoqPCR$DNACopiesperuLL<-bythoqPCR$DNACopiesperuL/bythoqPCR$VolumeFilteredL

bythoqPCR$DNACopiesperuLLconv[is.na(bythoqPCR$DNACopiesperuLLconv)]<-0
bythoqPCR$DNACopiesperuLL[is.na(bythoqPCR$DNACopiesperuLL)]<-0

bythoqPCRmeansconv<-aggregate(DNACopiesperuLLconv ~ Day + Tankno + Year, data=bythoqPCR, mean)
bythoqPCRmeans<-aggregate(DNACopiesperuLL ~ Day + Tankno + Year, data=bythoqPCR, mean)
write.csv(bythoqPCRmeans,"bythoqPCRmeans.csv")
#use this to create day of first detection dataset

#Set detection limit (Dil15)
bythoqPCRmeansconv$DNACopiesperuLLconv[bythoqPCRmeansconv$DNACopiesperuLL<0.0008] <- 0
bythoqPCRmeans$DNACopiesperuLL[bythoqPCRmeans$DNACopiesperuLL<0.0008] <- 0

#Create wide dataset based on day
bythoqPCRwideconv<-reshape(bythoqPCRmeansconv, idvar = c("Tankno","Year"), timevar = "Day", direction = "wide")
bythoqPCRwide<-reshape(bythoqPCRmeans, idvar = c("Tankno","Year"), timevar = "Day", direction = "wide")
#create total eDNA conc variable
bythoqPCRwideconv$totDNAcopies<-rowSums(bythoqPCRwideconv[,3:11],na.rm=T)
bythoqPCRwide$totDNAcopies<-rowSums(bythoqPCRwide[,3:11],na.rm=T)

#merge physical,environmental, and qPCR datasets long
bythoenvqPCRlongconv<-merge(bythoenv,bythoqPCRmeansconv,all=T)
bythoenvqPCRlong<-merge(bythoenv,bythoqPCRmeans,all=T)

bytho2122qPCRlongconv<-merge(bytho2122,bythoqPCRmeansconv, all=T)
bytho2122qPCRlong<-merge(bytho2122,bythoqPCRmeans, all=T)

bytho1qelongconv<-merge(bythoenvqPCRlongconv,bytho2122)
bytho1qelong<-merge(bythoenvqPCRlong,bytho2122)

#remove negative control tanks
bytho2qelongnoctconv<-subset(bytho1qelongconv,BarbNoIntro!="None")
bytho2qelongnoct<-subset(bytho1qelong,BarbNoIntro!="None")

#remove day of intro sampling
bytho2qelongnoct2<-subset(bytho2qelongnoct,Day>2)

#merge environmental, physical and qPCR datasets wide
bytho2122qPCRwideconv<-merge(bytho2122,bythoqPCRwideconv)
bytho2122qPCRwide<-merge(bytho2122,bythoqPCRwide)

bythoenvqPCRwideconv<-merge(bythoenvwide,bythoqPCRwideconv)
bythoenvqPCRwide<-merge(bythoenvwide,bythoqPCRwide)

bytho2qewideconv<-merge(bythoenvqPCRwideconv,bytho2122)
bytho2qewide<-merge(bythoenvqPCRwide,bytho2122)

#remove negative control tanks
bytho2qewidenoctconv<-subset(bytho2qewideconv,BarbNoIntro!="None")
bytho2qewidenoct<-subset(bytho2qewide,BarbNoIntro!="None")

#upload day of first detection dataset
bythoday1st<-read.csv("Day1stDetect.csv")
#merge with physical
day1phys<-merge(bythoday1st,bytho2122)
day1physnocont<-subset(day1phys,Barb!="NA")
#merge with environmental
Day1physenv<-merge(day1phys,bythoenvwide)
#subset into 2021 and 2022
day1phys21<-subset(day1phys, Year==2021)
day1phys22<-subset(day1phys, Year==2022)

day1physnocont21<-subset(day1phys,Barb!="NA" & Year==2021)
day1physnocont22<-subset(day1phys,Barb!="NA" & Year==2021)

Day1physenv21<-subset(Day1physenv, Year==2021)
Day1physenv22<-subset(Day1physenv, Year==2022)

#####################################################
#examine negative control tanks

options(scipen = 999)
ggplot(bytho2qelongnoct, aes(x=Day, y=DNACopiesperuLL, color=as.factor(Tankno))) + 
  geom_jitter(size=2, width=.1, height=.1) +
  scale_y_continuous(trans="log10")+
  labs(x ="Day", 
       y = "eDNA Copies per µL per L", color="Tank")+
  scale_color_manual(guide="none", values=tankcolvec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

ggplot(bytho2qelongnoctconv, aes(x=Day, y=DNACopiesperuLLconv)) + 
  geom_jitter(size=2, width=.1, height=.1) +
  scale_y_continuous(trans="log10")+
  labs(x ="Day", 
       y = "DNA Copies/uL/L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

#don't use converted, because doesn't make any difference, and just adds another method to explain.
#only two hits - could be because of swf in ambient water or contamination. ignore for now

#Examine correlations between physical datasets and eDNA signal

bythocorrelations<-cor(bytho2qewidenoct[,c(3:31,33,35:75,80)],use="pairwise.complete.obs")
#correlations greater than |0.8| between morph and eDNA
#one barb live and copies on day 6 0.95 correlation (and barren)
#copies on day 3 and one barb dead 0.99 (and barren)
#copies on day 3 and total dead 0.88 (and barren)

#plot these correlations
ggplot(bytho2qewidenoct, aes(x=OneBarbLiveBarren, y=DNACopiesperuLL.6)) + 
  geom_point() +
  labs(x ="Density of Live, 1-barb individuals at draining", 
       y = "DNA Copies per µL collected Day 6")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))
#only 0's and 1's, so not a confident correlation
ggplot(bytho2qewidenoct, aes(x=OneBarbDeadBarren, y=DNACopiesperuLL.3)) + 
  geom_point() +
  labs(x ="Density of Dead, 1-barb individuals at draining", 
       y = "DNA Copies per µL collected Day 3")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))
#only one individual found at the high value, so not a confident correlation

#most of the correlations are between other eDNA amounts in other days

#Normality tests for env and dna variables
shapiro.test(bytho2qelongnoct$DNACopiesperuLL)
#non normal p<0.01
hist(bytho2qelongnoct$DNACopiesperuLL)
shapiro.test(log10(bytho2qelongnoct$DNACopiesperuLL+1))
#non normal p<0.01
hist(log10(bytho2qelongnoct$DNACopiesperuLL+1))
#better but still very skewed. May need to run analyses as both quantitative and detect vs. undetect
boxplot(log10(bytho2qelongnoct$DNACopiesperuLL+1))
#may be an outlier
boxplot.stats(bytho2qelongnoct$DNACopiesperuLL)$out
#a lot of outliers, but one is particularly bad (3399.14716488 ). Remove that point in a new variable
bytho2qelongnoct$DNACpuLLno<- replace(bytho2qelongnoct$DNACopiesperuLL, bytho2qelongnoct$DNACopiesperuLL>3000, NA)
bytho2qelongnoct2$DNACpuLLno<- replace(bytho2qelongnoct2$DNACopiesperuLL, bytho2qelongnoct2$DNACopiesperuLL>3000, NA)

boxplot.stats(bytho2qelongnoct$DNACpuLLno)$out
#fixed the worst one
boxplot.stats(log10(bytho2qelongnoct$DNACpuLLno+1))$out
#log transforming helps
bytho2qelongnoct$logDNACpuLLno<-log10(bytho2qelongnoct$DNACpuLLno+1)
bytho2qelongnoct2$logDNACpuLLno<-log10(bytho2qelongnoct2$DNACpuLLno+1)

#binary variable
bytho2qelongnoct$DNADetect<- as.numeric(bytho2qelongnoct$DNACpuLL>0)
bytho2qelongnoct2$DNADetect<- as.numeric(bytho2qelongnoct2$DNACpuLL>0)

#Day of first detection
shapiro.test(Day1physenv$FirstDay)
#not normal
shapiro.test(log10(Day1physenv$FirstDay))
#Better
hist(Day1physenv$FirstDay)
hist(log10(Day1physenv$FirstDay))
#use log transformed

#TempC
shapiro.test(bytho2qelongnoct$TempC)
#normal 0.29
#Temp
ggplot(bythoenv, aes(x=Day, y=TempC)) + 
  geom_jitter(size=2, width=.3, height=.3) +
  labs(x ="Day", 
       y = "Temperature")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Year)
#Temperature higher on starting day in 2022 than 2021
#range ~19-20.5 - a difference of 1 degree can change the establishment threshold, so keep in
cor.test(bytho2qelongnoct$TempC,bytho2qelongnoct$Day)
#not correlated
cor.test(bytho2qelongnoct$TempC,bytho2qelongnoct$DOmgL)
#cor 0.288

#CondmScm
shapiro.test(bytho2qelongnoct$CondmScm)
#not normal 
hist(bytho2qelongnoct$CondmScm)
#binomial, very different for 2 experiments
ggplot(bythoenv, aes(x=Day, y=CondmScm)) + 
  geom_jitter(size=2, width=.3, height=.3) +
  labs(x ="Day", 
       y = "Conductivity")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Year)
#Conductivity much lower in 2022 compared to 2021
#173-177 in 2021 and 162-167 in 2022

#pH
shapiro.test(bytho2qelongnoct$pH)
#normal 0.45
ggplot(bythoenv, aes(x=Day, y=pH)) + 
  geom_jitter(size=2, width=.3, height=.3) +
  labs(x ="Day", 
       y = "pH")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Year)
#pH slightly lower to start in 2022
#range 7.5 to 8.5 - not a biologically relevant range, on the other hand barnes found significant change within this range so keep
cor.test(bytho2qelongnoct$pH,bytho2qelongnoct$Day)
#pH and day highly correlated, so use Day on analyses
source("group_by_summary_stats.R")
group_by_summary_stats(bytho2qelongnoct, pH, Day, Year)

#TurbidityFNU
shapiro.test(bytho2qelongnoct$TurbidityFNU)
#not normal
hist(bytho2qelongnoct$TurbidityFNU)
#several values very high
shapiro.test(log10(bytho2qelongnoct$TurbidityFNU))
hist(log10(bytho2qelongnoct$TurbidityFNU))
#log10 transformation helps, but still not normal
bytho2qelongnoct$logTurbidityFNU<-log10(bytho2qelongnoct$TurbidityFNU)
ggplot(bythoenv, aes(x=Day, y=TurbidityFNU)) + 
  geom_jitter(size=2, width=.3, height=.3) +
  labs(x ="Day", 
       y = "Turbidity")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Year)
#Turbidity much higher in 2022

#DOmgL
shapiro.test(bytho2qelongnoct$DOmgL)
#notnormal
hist(bytho2qelongnoct$DOmgL)
#right skewed, so don't transform
ggplot(bythoenv, aes(x=Day, y=DOmgL)) + 
  geom_jitter(size=2, width=.3, height=.3) +
  labs(x ="Day", 
       y = "Dissolved Oxygen")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Year)
#Do slighltly higher on innoculation day in 2022
#range 7.4-8.6, not an ecologically relevant range for bytho or eDNA so remove
cor.test(bytho2qelongnoct$DOmgL,bytho2qelongnoct$Day)
#not correlated
summary(lm(DOmgL~as.factor(Year)*Day,data=bytho2qelongnoct))
#no clear consistent pattern

#TotalChlorugL
shapiro.test(bytho2qelongnoct$TotalChlorugL)
#not normal
hist(bytho2qelongnoct$TotalChlorugL)
#binomial distribution between two experiments
ggplot(bythoenv, aes(x=Day, y=TotalChlorugL)) + 
  geom_jitter(size=2, width=.3, height=.3) +
  labs(x ="Day", 
       y = "Chlorophyl")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Year)
#Chlorophyl a much higher in 2022
#This range could cause a change in eDNA so keep
cor.test(bytho2qelongnoct$TotalChlorugL,bytho2qelongnoct$Day)
#not correlated
cor.test(bytho2qelongnoct$TotalChlorugL,bytho2qelongnoct$Year)
#very correlated
#The correlation is a through time, but different for each year. Since correlated with day, keep day and remove variable

#examine correlations with day as a variable (i.e., long)
bythocorrelationslong<-cor(bytho2qelongnoct[,c(3:10,12,14:41,46,48:50)],use="pairwise.complete.obs")
#correlations over |0.8|
write.csv(bythocorrelationslong,"bythocorrs.csv")

#CondmScm TotalChlorugL -0.94862350
#TurbidityFNU TotalChlorugL 0.89165369
#Didn't check correlations among bytho physical parameters
#keep total chlorophyl and remove cond and turb

#First look at day of first detection and relationship to what was added
#convert NA's to 0's for those where it was never detected
day1physnocont$FirstDay0<-day1physnocont$FirstDay
day1physnocont$FirstDay0[is.na(day1physnocont$FirstDay0)] <- 0
ggplot(day1physnocont, aes(x=NoIntro, y=FirstDay0)) + 
  geom_jitter(size=2, width=.3, height=.3) +
  labs(x ="Number Introduced", 
       y = "Day first detected")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Year)

#run model without environmental parameters
lm1stdetectdemfull<- lm(log10(FirstDay) ~ (NoIntro+Barb+as.factor(Year))^2,
                          data = day1phys)
stepAIC(lm1stdetectdemfull, direction="backward")
#Year significant, best model

#Run model again with environmental parameters as well as year
lm1stdetectenvfull<- lm(log10(FirstDay) ~ NoIntro*Barb+TempC.2+TotalChlorugL.2+pH.2+DOmgL.2+as.factor(Year),
                 data = Day1physenv)
stepAIC(lm1stdetectenvfull, direction="backward")
#best model has pH and Dissolved oxygen
lm1stdetectbest<-lm(log10(FirstDay)~pH.2+DOmgL.2,data=Day1physenv)
summary(lm1stdetectbest)
#multiple rsquared 0.63
#log10(firstday)=4.65+-0.78pH+0.24DO
#pH std error 0.26 t value -3.01 p <0.01
#DO std error 0.08 t value 3.16 p <0.01

#plot pH
ggplot(Day1physenv, aes(x=pH.2, y=(FirstDay-2))) + 
  geom_point(aes(color=factor(Year))) +
  geom_smooth(method=lm, se=F, color="black")+
  scale_y_continuous(trans="log10")+
  labs(x ="pH at Inoculation", 
       y = "Day of First Detection", color="Year")+
  scale_color_manual(values=yearcolvec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))
#plot DO
ggplot(Day1physenv, aes(x=DOmgL.2, y=(FirstDay-2))) + 
  geom_point(aes(color=factor(Year))) +
  geom_smooth(method=lm, se=F, color="black", )+
  scale_y_continuous(trans="log10")+
  labs(x ="Dissolved Oxygen (mg/L) at Introduction", 
       y = "Day of First Detection", color="Year")+
  scale_color_manual(values=yearcolvec)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

#repeat linear model with nontransformed day first detect so easier to interpret (since doesn't change result)
lm1stdetectbestut<-lm(FirstDay~pH.2+DOmgL.2,data=Day1physenv)
summary(lm1stdetectbestut)
#multiple rsquared 0.54
#firstday=62.388+-11.397pH+3.62DO
#pH std error 4.57 t value -2.50 p=0.02
#DO std error 1.335 t value 2.708 p=0.01

#Run for each year separately, since different design: 2021
day1physnocont21$FirstDay0<-day1physnocont21$FirstDay
day1physnocont21$FirstDay0[is.na(day1physnocont21$FirstDay0)] <- 0

#run model without environmental parameters
lm1stdetectdemfull21<- lm(log10(FirstDay) ~ NoIntro*Barb, data = day1phys21)
stepAIC(lm1stdetectdemfull21, direction="both")
#nothing significant

#Run model again with environmental parameters as well as year
lm1stdetectenvfull21<- lm(log10(FirstDay) ~ NoIntro*Barb+TempC.2+TotalChlorugL.2+pH.2+DOmgL.2,
                        data = Day1physenv21)
stepAIC(lm1stdetectenvfull21, direction="both")
#best model has chlorophyl, pH and, DO
lm1stdetectbest21<-lm(formula = log10(FirstDay) ~ NoIntro + TotalChlorugL.2 + pH.2 + 
                        DOmgL.2, data = Day1physenv21)
summary(lm1stdetectbest21)
#overfit, so delete chlorophyl
lm1stdetectenvfull21red<- lm(log10(FirstDay) ~ NoIntro*Barb+TempC.2+pH.2+DOmgL.2,
                          data = Day1physenv21)
stepAIC(lm1stdetectenvfull21red, direction="both")
lm1stdetectbest21red<-lm(formula = log10(FirstDay) ~ Barb + TempC.2 + pH.2, data = Day1physenv21)
summary(lm1stdetectbest21red)
#overfit, so delete barb
lm1stdetectenvfull21red<- lm(log10(FirstDay) ~ NoIntro+TempC.2+pH.2+DOmgL.2,
                             data = Day1physenv21)
stepAIC(lm1stdetectenvfull21red, direction="both")
lm1stdetectbest21red<-lm(formula = log10(FirstDay) ~ NoIntro + TempC.2 + pH.2, data = Day1physenv21)
summary(lm1stdetectbest21red)
#overfit, so delete no intro
lm1stdetectenvfull21red<- lm(log10(FirstDay) ~ TempC.2+pH.2+DOmgL.2,
                             data = Day1physenv21)
stepAIC(lm1stdetectenvfull21red, direction="both")
lm1stdetectbest21red<-lm(formula = log10(FirstDay) ~ pH.2, data = Day1physenv21)
summary(lm1stdetectbest21red)
#pH significant
#p=0.05. -0.78 estimate, standard error 0.35, t value -2.18
#for interpretation, run with linear response
lm1stdetectbest21redlin<-lm(formula = FirstDay ~ pH.2, data = Day1physenv21)
summary(lm1stdetectbest21redlin)


#Try again with 2022
#Run for each year separately, since different design: 2022
day1physnocont22$FirstDay0<-day1physnocont22$FirstDay
day1physnocont22$FirstDay0[is.na(day1physnocont22$FirstDay0)] <- 0

#run model without environmental parameters
lm1stdetectdemfull22<- lm(log10(FirstDay) ~ NoIntro*Barb, data = day1phys22)
stepAIC(lm1stdetectdemfull22, direction="both")
#nothing significant

#Run model again with environmental parameters as well as year
lm1stdetectenvfull22<- lm(log10(FirstDay) ~ NoIntro*Barb+TempC.2+TotalChlorugL.2+pH.2+DOmgL.2,
                          data = Day1physenv22)
stepAIC(lm1stdetectenvfull22, direction="both")
lm1stdetectbest22<-lm(formula = log10(FirstDay) ~ NoIntro + Barb + TempC.2 + TotalChlorugL.2 + 
                        pH.2, data = Day1physenv22)
summary(lm1stdetectbest22)
#overfit, so delete no intro
lm1stdetectenvfull22red<- lm(log10(FirstDay) ~ Barb+TempC.2+TotalChlorugL.2+pH.2+DOmgL.2,
                          data = Day1physenv22)
stepAIC(lm1stdetectenvfull22red, direction="both")
lm1stdetectbest22red<-lm(formula = log10(FirstDay) ~ pH.2 + DOmgL.2, data = Day1physenv22)
summary(lm1stdetectbest22red)
#overfit, so delete pH
lm1stdetectbest22red<-lm(formula = log10(FirstDay) ~ DOmgL.2, data = Day1physenv22)
summary(lm1stdetectbest22red)
#nothing significant

#Test effect of day, demographics, environmental variables on DNA concentration

#look at correlations among all variables in this model to do variable reduction

#build linear model using stewise AIC
names(bytho2qelongnoct2)
summary(bytho2qelongnoct2)
#Creat new dataset with NA's removed
bytho2qelongnoct2nar<-subset(bytho2qelongnoct2,logDNACpuLLno!="NA" & TempC!="NA")

#don't use each individual tank/year because it's the same as all the bytho population variables
#don't use 2 barb live barren, because only 2 individuals found, so not a lot of variability
#try with counts first, then do DW
#try with divided counts by demographics rather than pooled
#first use subset dataset to use env variables
fulllmbythocountenv<-lm(logDNACpuLLno ~ (Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                       ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                       ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound+
                       as.factor(Year))*Day+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(fulllmbythocountenv, direction="both")
#AIC -147

bestlmbythocountenv<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + OneBarbLiveBarren + 
                          TwoBarbLiveGravid + ThreeBarbLiveGravid + OneBarbDeadBarren + 
                          TwoBarbDeadGravid + ThreeBarbDeadGravid + OneBarbSpine + 
                          TwoBarbSpine + ThreeBarbSpine + kink + TotalMoultsRound + 
                          as.factor(Year) + Day + TempC + Barb:NoIntro + Barb:Day + 
                          NoIntro:Day + OneBarbLiveBarren:Day + TwoBarbLiveGravid:Day + 
                          ThreeBarbLiveGravid:Day + OneBarbDeadBarren:Day + TwoBarbDeadGravid:Day + 
                          ThreeBarbDeadGravid:Day + OneBarbSpine:Day + ThreeBarbSpine:Day + 
                          kink:Day + TotalMoultsRound:Day, data = bytho2qelongnoct2nar)
summary(bestlmbythocountenv)
#r squared 0.88, adjusted r squared is 0.68, so definitely overfit. Try again with less variables
fulllmbythocountenvred<-lm(logDNACpuLLno ~ Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                                           ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                                           ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound+
                                           as.factor(Year)+Day+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(fulllmbythocountenvred, direction="both")
#AIC -108

bestlmbythocountenvred<-lm(formula = logDNACpuLLno ~ NoIntro + TwoBarbLiveGravid + ThreeBarbDeadGravid + 
                             OneBarbSpine + TwoBarbSpine + ThreeBarbSpine + TotalMoultsRound + 
                             as.factor(Year) + TempC, data = bytho2qelongnoct2nar)
summary(bestlmbythocountenvred)
#No environmental predictors so try without env

#the full model without the env
fulllmbythocount<-lm(logDNACpuLLno ~ (Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                          ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                          ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound)*Day+
                       as.factor(Year),data = bytho2qelongnoct2)
stepAIC(fulllmbythocount, direction="both")
#AIC -406

bestlmbythocount<-lm(formula = logDNACpuLLno ~ Barb + TwoBarbLiveGravid + ThreeBarbDeadGravid + 
                       OneBarbSpine + TwoBarbSpine + kink + TotalMoultsRound + Day + 
                       TwoBarbLiveGravid:Day + ThreeBarbDeadGravid:Day + OneBarbSpine:Day + 
                       TwoBarbSpine:Day + TotalMoultsRound:Day, data = bytho2qelongnoct2)
summary(bestlmbythocount)
#overfit, so try combining spines and kinks
fulllmbythocountred<-lm(logDNACpuLLno ~ (Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                                        ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                                        ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound)*Day+
                       as.factor(Year),data = bytho2qelongnoct2)
stepAIC(fulllmbythocountred, direction="both")
#AIC -407

bestlmbythocountred<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TwoBarbDeadGravid + 
                          ThreeBarbDeadGravid + TotalMoultsRound + Day + NoIntro:Day + 
                          TwoBarbDeadGravid:Day + ThreeBarbDeadGravid:Day + TotalMoultsRound:Day, 
                        data = bytho2qelongnoct2)
summary(bestlmbythocountred)
#overfit
fulllmbythocountred<-lm(logDNACpuLLno ~ Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                                           ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                                           ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+Day+
                          as.factor(Year),data = bytho2qelongnoct2)
stepAIC(fulllmbythocountred, direction="both")
#AIC -406

bestlmbythocountred<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + ThreeBarbDeadGravid + 
                          Day, data = bytho2qelongnoct2)
summary(bestlmbythocountred)
#Barb and Day significant
#Do with untransformed response to interprettability
bestlmbythocountrednt<-lm(formula = DNACpuLLno ~ Barb + NoIntro + ThreeBarbDeadGravid + 
                          Day, data = bytho2qelongnoct2)
summary(bestlmbythocountrednt)

#plot

ggplot(bytho2qelongnoct2, aes(x=Day, y=logDNACpuLLno, group=Day)) + 
  geom_boxplot() +
  labs(x ="Day", 
       y = "log(DNA Copies per µL per L)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))
ggplot(bytho2qelongnoct2, aes(x=Day, y=DNACpuLLno)) + 
  geom_jitter() +
  scale_y_continuous(trans="log10")+
  labs(x ="Day", 
       y = "log(DNA Copies per µL per L)")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))
#These effects extremely subtle

#Try again with DW
fulllmbythobm<-lm(logDNACpuLLno ~ (Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                    ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                    TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+OneBarbSpineDW+TwoBarbSpineDW+
                    ThreeBarbSpineDW+kinkDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL)*Day,
                  data = bytho2qelongnoct2)
stepAIC(fulllmbythobm, direction="both")
#AIC -149

bestlmbythobm<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + OneBarbLiveBarrenDW + 
                    TwoBarbLiveGravidDW + ThreeBarbLiveGravidDW + OneBarbDeadDWBarren + 
                    TwoBarbDeadGravidDW + ThreeBarbDeadGravidDW + OneBarbSpineDW + 
                    TwoBarbSpineDW + ThreeBarbSpineDW + kinkDW + TotalMoultsRound + 
                    as.factor(Year) + TempC + Day + Barb:StartLivingBiomass + 
                    Barb:Day + StartLivingBiomass:Day + OneBarbLiveBarrenDW:Day + 
                    TwoBarbLiveGravidDW:Day + ThreeBarbLiveGravidDW:Day + OneBarbDeadDWBarren:Day + 
                    TwoBarbDeadGravidDW:Day + ThreeBarbDeadGravidDW:Day + ThreeBarbSpineDW:Day + 
                    kinkDW:Day + TotalMoultsRound:Day + TempC:Day, data = bytho2qelongnoct2)
summary(bestlmbythobm)
#rsquared shows overfit, switch to without environmental

fulllmbythobmdem<-lm(logDNACpuLLno ~ (Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                                     ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+OneBarbSpineDW+TwoBarbSpineDW+
                                     ThreeBarbSpineDW+kinkDW+TotalMoultsRound+as.factor(Year))*Day,
                  data = bytho2qelongnoct2)
stepAIC(fulllmbythobmdem, direction="both")
#AIC=408

bestlmbythobmdem<-lm(formula = logDNACpuLLno ~ TwoBarbLiveGravidDW + ThreeBarbDeadGravidDW + 
                       TwoBarbSpineDW + kinkDW + TotalMoultsRound + Day + Barb + 
                       TwoBarbLiveGravidDW:Day + ThreeBarbDeadGravidDW:Day + TwoBarbSpineDW:Day + 
                       TotalMoultsRound:Day, data = bytho2qelongnoct2)
summary(bestlmbythobmdem)
#not overfit, 162 df residual standard error 0.18 p<0.01, 
#all extremely small estimates, so try pooling totalspineskinks

fulllmbythobmdemred<-lm(logDNACpuLLno ~ (Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                                        ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                                        TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+TotalMoultsRound+as.factor(Year))*Day,
                     data = bytho2qelongnoct2)
stepAIC(fulllmbythobmdemred, direction="both")
#AIC411
bestlmbythobmdemred<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + ThreeBarbLiveGravidDW + 
                          ThreeBarbDeadGravidDW + Day + ThreeBarbLiveGravidDW:Day, 
                        data = bytho2qelongnoct2)
summary(bestlmbythobmdemred)
#3 barblivegraviddw*Day interaction extremely small, so take out and rerun
fulllmbythobmdemred<-lm(logDNACpuLLno ~ Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                                           ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                                           TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+TotalMoultsRound+as.factor(Year)+Day,
                        data = bytho2qelongnoct2)
stepAIC(fulllmbythobmdemred, direction="both")
#AIC `406
bestlmbythobmdemred<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + ThreeBarbDeadGravidDW + 
                          Day, data = bytho2qelongnoct2)
summary(bestlmbythobmdemred)
#Day significant, positive association

#Try again with pooled demographic variables ()
fulllmbythocountenvgb<-lm(logDNACpuLLno ~ (BarbNoIntro+TotalGravidLive+TotalGravidDead+TotalBarrenLive+
                                             OneBarbDeadBarren+TotalSpinesKinks+TotalMoultsRound+
                                           as.factor(Year))*Day+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(fulllmbythocountenvgb, direction="both")
#AIC=-129

bestlmbythocountenvgb<-lm(formula = logDNACpuLLno ~ BarbNoIntro + TotalGravidLive + 
                            TotalGravidDead + TotalBarrenLive + OneBarbDeadBarren + TotalSpinesKinks + 
                            TotalMoultsRound + as.factor(Year) + Day + TempC + BarbNoIntro:Day + 
                            TotalGravidLive:Day + TotalBarrenLive:Day + OneBarbDeadBarren:Day + 
                            TotalSpinesKinks:Day + TotalMoultsRound:Day, data = bytho2qelongnoct2nar)
summary(bestlmbythocountenvgb)
#extremely overfit - take out env variables

fulllmbythocountgb<-lm(logDNACpuLLno ~ (Barb*NoIntro+TotalGravidLive+TotalGravidDead+TotalBarrenLive+
                                          OneBarbDeadBarren+TotalSpinesKinks+TotalMoultsRound+
                                             as.factor(Year))*Day,data = bytho2qelongnoct2)
stepAIC(fulllmbythocountgb, direction="both")
#AIC=-410

bestlmbythocountgb<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TotalGravidDead + 
                         TotalMoultsRound + Day + NoIntro:Day + TotalGravidDead:Day + 
                         TotalMoultsRound:Day, data = bytho2qelongnoct2)
summary(bestlmbythocountgb)
#not overfit, p<0.01, rsquared 0.16, Barb, no intro, total moults, total gravid dead*Day and totalmoults*Day significant
#Since total moults has the opposite direction estimate compared to earlier model, take out and re-run

fulllmbythocountgbred<-lm(logDNACpuLLno ~ (Barb*NoIntro+TotalGravidLive+TotalGravidDead+TotalBarrenLive+
                                          OneBarbDeadBarren+TotalSpinesKinks+
                                          as.factor(Year))*Day,data = bytho2qelongnoct2)
stepAIC(fulllmbythocountgbred, direction="both")
#AIC=-409

bestlmbythocountgbred<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TotalGravidLive + 
                            TotalGravidDead + Day + TotalGravidLive:Day + TotalGravidDead:Day, 
                          data = bytho2qelongnoct2)
summary(bestlmbythocountgbred)
#rsquared 0.15, not overfit, p<0.01
#barb and totalgravidlive*Day interaction significant

#Visualize
ggplot(bytho2qelongnoct, aes(x=TotalGravidLive, y=DNACpuLLno)) + 
  geom_point() +
  scale_y_continuous(trans="log10")+
  labs(x ="3-Barb Live Gravid Females at the end of the experiment", 
       y = "DNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(~Day)

bytho2qelongnoct$Barb<-revalue(bytho2qelongnoct$Barb, c("2barb"="2nd instar", "3barb"="3rd instar"))
ggplot(bytho2qelongnoct, aes(x=Barb, y=DNACpuLLno)) + 
  geom_boxplot() +
  scale_y_continuous(trans="log10")+
  labs(x ="Life Stage Inoculated", 
       y = "eDNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

#Pooled demographic variables and biomass
fulllmbythobmenvgb<-lm(logDNACpuLLno ~ (BarbNoIntro+GravidLiveBiomass+GravidDeadBiomass+BarrenLiveBiomass+
                                          OneBarbDeadDWBarren+as.factor(Year))*Day+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(fulllmbythobmenvgb, direction="both")
#AIC=-124.2

bestlmbythobmenvgb<-lm(formula = logDNACpuLLno ~ BarbNoIntro + TotalGravidLive + 
                            TotalGravidDead + TotalBarrenLive + OneBarbDeadBarren + TotalSpinesKinks + 
                            as.factor(Year) + Day + TempC + BarbNoIntro:Day + TotalBarrenLive:Day + 
                            OneBarbDeadBarren:Day + TotalSpinesKinks:Day, data = bytho2qelongnoct2nar)
summary(bestlmbythobmenvgb)
#extremely overfit - take out env variables

fulllmbythobmgb<-lm(logDNACpuLLno ~ (Barb*NoIntro+GravidLiveBiomass+GravidDeadBiomass+BarrenLiveBiomass+
                                          OneBarbDeadDWBarren+as.factor(Year))*Day,data = bytho2qelongnoct2)
stepAIC(fulllmbythobmgb, direction="both")
#AIC=-410

bestlmbythobmgb<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + GravidLiveBiomass + 
                      GravidDeadBiomass + Day + GravidLiveBiomass:Day + GravidDeadBiomass:Day, 
                    data = bytho2qelongnoct2)
summary(bestlmbythobmgb)
#gravidlivebiomass*Day interaction so tiny remove
fulllmbythobmgbred<-lm(logDNACpuLLno ~ Barb*NoIntro+GravidLiveBiomass+GravidDeadBiomass+BarrenLiveBiomass+
                                       OneBarbDeadDWBarren+as.factor(Year)+Day,data = bytho2qelongnoct2)
stepAIC(fulllmbythobmgbred, direction="both")
#AIC=-406

bestlmbythobmgbred<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + Day + OneBarbDeadDWBarren, 
                       data = bytho2qelongnoct2)
summary(bestlmbythobmgbred)
#Barb and Day significant

#Try with pooled counts
fulllmbythocountenvld<-lm(logDNACpuLLno ~ (Barb*NoIntro+TotalLive+TotalDead+TotalSpinesKinks+as.factor(Year))*Day+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(fulllmbythocountenvld, direction="both")
#AIC=-120

bestlmbythocountenvld<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TotalLive + TotalDead + 
                            TotalSpinesKinks + as.factor(Year) + Day + TempC + NoIntro:Day + 
                            TotalDead:Day + TotalSpinesKinks:Day, data = bytho2qelongnoct2nar)
summary(bestlmbythocountenvld)
#extremely overfit - take out env variables
fulllmbythocountld<-lm(logDNACpuLLno ~ (Barb*NoIntro+TotalLive+TotalDead+TotalSpinesKinks+as.factor(Year))*Day,data = bytho2qelongnoct2)
stepAIC(fulllmbythocountld, direction="both")
#AIC=-405

bestlmbythocountld<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TotalDead + TotalSpinesKinks + 
                         Day + NoIntro:Day + TotalSpinesKinks:Day, data = bytho2qelongnoct2)
summary(bestlmbythocountld)
#totalspineskinks*day significant, but extremely small, so take out interaction
fulllmbythocountldred<-lm(logDNACpuLLno ~ Barb*NoIntro+TotalLive+TotalDead+TotalSpinesKinks+as.factor(Year)+Day,data = bytho2qelongnoct2)
stepAIC(fulllmbythocountldred, direction="both")
#AIC=-406

bestlmbythocountldred<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TotalDead + Day, 
                          data = bytho2qelongnoct2)
summary(bestlmbythocountldred)
#Day and barb added significant

#now with biomass
fulllmbythobmenvld<-lm(logDNACpuLLno ~ (Barb*StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+as.factor(Year))*Day+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(fulllmbythobmenvld, direction="both")
#AIC=-114

bestlmbythobmenvld<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + EndLivingBiomass + 
                         EndDeadBiomass + as.factor(Year) + Day + TempC + StartLivingBiomass:Day + 
                         EndLivingBiomass:Day, data = bytho2qelongnoct2nar)
summary(bestlmbythobmenvld)
#overfit, take out env variables

fulllmbythobmld<-lm(logDNACpuLLno ~ (Barb*StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+as.factor(Year))*Day,data = bytho2qelongnoct2)
stepAIC(fulllmbythobmld, direction="both")
#AIC=-407

bestlmbythobmld<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + EndLivingBiomass + 
                      Day + EndLivingBiomass:Day, data = bytho2qelongnoct2)
summary(bestlmbythobmld)
#Extremely small estimates for interaction, so take out
fulllmbythobmldred<-lm(logDNACpuLLno ~ Barb*StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+as.factor(Year)+Day,data = bytho2qelongnoct2)
stepAIC(fulllmbythobmldred, direction="both")
#AIC=-405

bestlmbythobmldred<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + Day, 
                       data = bytho2qelongnoct2)
summary(bestlmbythobmldred)
#Start living biomass and Day significant, much lower multiple r squared then other models
ggplot(bytho2qelongnoct, aes(x=StartLivingBiomass, y=DNACpuLLno)) + 
  geom_point() +
  scale_y_continuous(trans="log10")+
  labs(x ="Starting Live Biomass", 
       y = "DNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

#Summary of models
#full demographics counts - barb and day sig, AIC=-406, rsquared = 0.1
#full demographics bm - day sig, AIC-408, rsquared 0.1
#gravid counts - barb and day*totalgravidLive sig, AIC-409, rsquared 0.15
#gravid bm - Barba nd Day significant, AIC -406, rsquared 0.1
#total counts - Barb and Day significant, AIC -406, rquared 0.1
#total bm - start living biomass and day significant, AIC-405, rsquared 0.09

#Final result - barb and day significant


###Use binary response variable now
fulllmbythocountenvbi<-lm(DNADetect ~ (Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                                           ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                                           ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound+
                                           as.factor(Year))*Day+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(fulllmbythocountenvbi, direction="both")
#AIC -66
bestlmbythocountenvbi<-lm(formula = DNADetect ~ Barb + NoIntro + OneBarbLiveBarren + 
                          TwoBarbLiveGravid + ThreeBarbLiveGravid + OneBarbDeadBarren + 
                          TwoBarbDeadGravid + ThreeBarbDeadGravid + OneBarbSpine + 
                          TwoBarbSpine + ThreeBarbSpine + kink + TotalMoultsRound + 
                          as.factor(Year) + Day + DOmgL + Barb:Day + NoIntro:Day + 
                          OneBarbLiveBarren:Day + TwoBarbLiveGravid:Day + ThreeBarbLiveGravid:Day + 
                          OneBarbDeadBarren:Day + OneBarbSpine:Day + TwoBarbSpine:Day + 
                          kink:Day + TotalMoultsRound:Day, data = bytho2qelongnoct2nar)
summary(bestlmbythocountenvbi)
#not significant, remove day interaction
redlmbythocountenvbi<-lm(DNADetect ~ Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                                         ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                                         ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound+
                                         as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(redlmbythocountenvbi, direction="both")
#AIC -63
bestlmbythocountenvbired<-lm(formula = DNADetect ~ TwoBarbLiveGravid + ThreeBarbLiveGravid + 
                               ThreeBarbDeadGravid + OneBarbSpine + TwoBarbSpine + kink + 
                               TotalMoultsRound + as.factor(Year), data = bytho2qelongnoct2nar)
summary(bestlmbythocountenvbired)
#overfit, remove environmental
redlmbythocountenvbi<-lm(DNADetect ~ Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                           ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                           ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound+
                           as.factor(Year),data = bytho2qelongnoct2)
stepAIC(redlmbythocountenvbi, direction="both")
#AIC -249
bestlmbythocountenvbired<-lm(formula = DNADetect ~ TwoBarbLiveGravid + ThreeBarbDeadGravid + 
                               OneBarbSpine + TotalMoultsRound + as.factor(Year), data = bytho2qelongnoct2)
summary(bestlmbythocountenvbired)
#explains so little variation, not a good model. nothing significant



#Try again with DW
fulllmbythobmbi<-lm(DNADetect ~ (Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                                     ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+OneBarbSpineDW+TwoBarbSpineDW+
                                     ThreeBarbSpineDW+kinkDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL)*Day,
                  data = bytho2qelongnoct2)
stepAIC(fulllmbythobmbi, direction="both")
#AIC -149
bestlmbythobmbi<-lm(formula = DNADetect ~ Barb + StartLivingBiomass + OneBarbLiveBarrenDW + 
                      TwoBarbLiveGravidDW + ThreeBarbLiveGravidDW + OneBarbDeadDWBarren + 
                      TwoBarbDeadGravidDW + ThreeBarbDeadGravidDW + OneBarbSpineDW + 
                      TwoBarbSpineDW + ThreeBarbSpineDW + kinkDW + TotalMoultsRound + 
                      as.factor(Year) + TempC + pH + Day + Barb:Day + StartLivingBiomass:Day + 
                      OneBarbLiveBarrenDW:Day + TwoBarbLiveGravidDW:Day + ThreeBarbLiveGravidDW:Day + 
                      OneBarbDeadDWBarren:Day + TwoBarbDeadGravidDW:Day + ThreeBarbDeadGravidDW:Day + 
                      OneBarbSpineDW:Day + ThreeBarbSpineDW:Day + kinkDW:Day + 
                      TotalMoultsRound:Day + TempC:Day, data = bytho2qelongnoct2)
summary(bestlmbythobmbi)
#overfit, take out interaction and environmental
redlmbythobmbi<-lm(DNADetect ~ Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                                   ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                                   TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+OneBarbSpineDW+TwoBarbSpineDW+
                                   ThreeBarbSpineDW+kinkDW+TotalMoultsRound+as.factor(Year),
                    data = bytho2qelongnoct2)
stepAIC(redlmbythobmbi, direction="both")
#AIC -149
bestlmbythobmbired<-lm(formula = DNADetect ~ TwoBarbLiveGravidDW + ThreeBarbDeadGravidDW + 
                         OneBarbSpineDW + TotalMoultsRound + as.factor(Year), data = bytho2qelongnoct2)
summary(bestlmbythobmbired)
#very very small rsquared again, nothing significant


#Try again with pooled demographic variables ()
fulllmbythocountenvgbbi<-lm(DNADetect ~ (BarbNoIntro+TotalGravidLive+TotalGravidDead+TotalBarrenLive+
                                             OneBarbDeadBarren+TotalSpinesKinks+TotalMoultsRound+
                                             as.factor(Year))*Day+TempC+DOmgL,data = bytho2qelongnoct2nar)
stepAIC(fulllmbythocountenvgbbi, direction="both")
#AIC=-80

bestlmbythocountenvgbbi<-lm(formula = DNADetect ~ BarbNoIntro + OneBarbDeadBarren + TotalMoultsRound + 
                              Day + BarbNoIntro:Day + TotalMoultsRound:Day, data = bytho2qelongnoct2nar)
summary(bestlmbythocountenvgbbi)
#remove interaction and environmental variablres
redlmbythocountgbbi<-lm(DNADetect ~ BarbNoIntro+TotalGravidLive+TotalGravidDead+TotalBarrenLive+
                                           OneBarbDeadBarren+TotalSpinesKinks+TotalMoultsRound+
                                           as.factor(Year)+Day,data = bytho2qelongnoct2)
stepAIC(redlmbythocountgbbi, direction="both")
#AIC=-245

bestlmbythocountenvgbbi<-lm(formula = DNADetect ~ TotalMoultsRound + as.factor(Year) + 
                              Day, data = bytho2qelongnoct2)
summary(bestlmbythocountenvgbbi)
#rsquared so small, nothing significant

#summary of models
#count demoraphics - nothing significant
#biomass demographics, nothing significant
#counts gravid - nothing significant






#Last day conditions

#Now look at correlations among the last day conditions
ggplot(subset(bytho2qelongnoct, Day==16), aes(x=ThreeBarbLiveGravid, y=DNACpuLLno)) + 
  geom_point() +
  scale_y_continuous(trans="log10")+
  labs(x ="End 3-Barb Live Gravid", 
       y = "End DNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

#linear models for just the last day
bytho2qelongnoctlast<-subset(bytho2qelongnoct, Day==16)

#Examine normality again
#Normality tests for env and dna variables
shapiro.test(bytho2qelongnoctlast$DNACopiesperuLL)
#non normal p<0.01
hist(bytho2qelongnoctlast$DNACopiesperuLL)
shapiro.test(log10(bytho2qelongnoctlast$DNACopiesperuLL+1))
#non normal p<0.01
hist(log10(bytho2qelongnoctlast$DNACopiesperuLL+1))
#better but still very skewed. May need to run analyses as both quantitative and detect vs. undetect
bytho2qelongnoctlast$logDNACpuLLno<-log10(bytho2qelongnoctlast$DNACpuLLno+1)

#binary variable
bytho2qelongnoctlast$DNADetect<- as.numeric(bytho2qelongnoctlast$DNACpuLL>0)

#TempC
shapiro.test(bytho2qelongnoctlast$TempC)
#normal 0.71

#CondmScm
shapiro.test(bytho2qelongnoctlast$CondmScm)
#not normal 
hist(bytho2qelongnoctlast$CondmScm)
#binomial, very different for 2 experiments

#pH
shapiro.test(bytho2qelongnoctlast$pH)
#normal 0.47

#TurbidityFNU
shapiro.test(bytho2qelongnoctlast$TurbidityFNU)
#not normal
hist(bytho2qelongnoctlast$TurbidityFNU)
#several values very high
shapiro.test(log10(bytho2qelongnoctlast$TurbidityFNU))
hist(log10(bytho2qelongnoctlast$TurbidityFNU))
#log10 transformation helps, but still not normal
bytho2qelongnoctlast$logTurbidityFNU<-log10(bytho2qelongnoctlast$TurbidityFNU)

#DOmgL
shapiro.test(bytho2qelongnoctlast$DOmgL)
#notnormal
hist(bytho2qelongnoctlast$DOmgL)
#right skewed, so don't transform

#TotalChlorugL
shapiro.test(bytho2qelongnoctlast$TotalChlorugL)
#not normal
hist(bytho2qelongnoctlast$TotalChlorugL)
#binomial distribution between two experiments

bytho16correlations<-cor(bytho2qelongnoctlast[,c(2:10,12,14:55,61:64)],use="pairwise.complete.obs")
#Conductivity, total chlorophyl, and turbidity and year correlated, keep year

#build linear model using stewise AIC
names(bytho2qelongnoctlast)
summary(bytho2qelongnoctlast)
#Creat new dataset with NA's removed
bytho2qelongnoctlastnar<-subset(bytho2qelongnoctlast,logDNACpuLLno!="NA" & TempC!="NA")

#don't use each individual tank/year because it's the same as all the bytho population variables
#try with counts first, then do DW
#try with divided counts by demographics rather than pooled
#first use subset dataset to use env variables
fulllmbythocountenv16<-lm(logDNACpuLLno ~ Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                                           ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                                           ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound+
                                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(fulllmbythocountenv16, direction="both")
#AIC -94

bestlmbythocountenv16<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TwoBarbLiveGravid + 
                          ThreeBarbLiveGravid + ThreeBarbDeadGravid + OneBarbSpine + 
                          TwoBarbSpine + ThreeBarbSpine + kink + TotalMoultsRound + 
                          as.factor(Year) + TempC + Barb:NoIntro, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenv16)
#definitely overfit. Try again without interaction
redlmbythocountenv16<-lm(logDNACpuLLno ~ Barb+NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                            ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                            ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound+
                            as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -95

bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ Barb + TwoBarbLiveGravid + ThreeBarbLiveGravid + 
                            OneBarbDeadBarren + ThreeBarbDeadGravid + OneBarbSpine + 
                            TwoBarbSpine + ThreeBarbSpine + TotalMoultsRound + as.factor(Year) + 
                            TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#definitely overfit. combining spines
redlmbythocountenv16<-lm(logDNACpuLLno ~ Barb+NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                           ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                           ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -95

bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ Barb + ThreeBarbLiveGravid + OneBarbDeadBarren + 
                               ThreeBarbDeadGravid + TotalSpinesKinks + TotalMoultsRound + 
                               as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#definitely overfit. remove OneBarbDeadBarren
redlmbythocountenv16<-lm(logDNACpuLLno ~ Barb+NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                           ThreeBarbLiveGravid+TwoBarbDeadGravid+
                           ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -95

bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ Barb + OneBarbLiveBarren + ThreeBarbLiveGravid + 
                               ThreeBarbDeadGravid + TotalSpinesKinks + TotalMoultsRound + 
                               as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#overfit, remove OneBarbLiveBarren
redlmbythocountenv16<-lm(logDNACpuLLno ~ Barb+NoIntro+TwoBarbLiveGravid+
                           ThreeBarbLiveGravid+TwoBarbDeadGravid+
                           ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -94

bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ Barb + ThreeBarbLiveGravid + ThreeBarbDeadGravid + 
                               TotalSpinesKinks + TotalMoultsRound + as.factor(Year) + TempC, 
                             data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#overfit, remove barb
redlmbythocountenv16<-lm(logDNACpuLLno ~ NoIntro+TwoBarbLiveGravid+
                           ThreeBarbLiveGravid+TwoBarbDeadGravid+
                           ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -94

bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ TwoBarbLiveGravid + ThreeBarbLiveGravid + 
                               TwoBarbDeadGravid + ThreeBarbDeadGravid + TotalSpinesKinks + 
                               TotalMoultsRound + as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#overfit, remove TwoBarbDeadGravid
redlmbythocountenv16<-lm(logDNACpuLLno ~ NoIntro+TwoBarbLiveGravid+
                           ThreeBarbLiveGravid+ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -94
bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ TwoBarbLiveGravid + ThreeBarbLiveGravid + 
                               ThreeBarbDeadGravid + TotalSpinesKinks + TotalMoultsRound + 
                               as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#remove TwoBarbLiveGravid
redlmbythocountenv16<-lm(logDNACpuLLno ~ NoIntro+ThreeBarbLiveGravid+ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -94
bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ ThreeBarbLiveGravid + ThreeBarbDeadGravid + 
                               TotalSpinesKinks + TotalMoultsRound + as.factor(Year) + TempC, 
                             data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#overfit, remove TotalMoultsRound
redlmbythocountenv16<-lm(logDNACpuLLno ~ NoIntro+ThreeBarbLiveGravid+ThreeBarbDeadGravid+TotalSpinesKinks+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -88
bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ ThreeBarbLiveGravid + ThreeBarbDeadGravid + 
                               TotalSpinesKinks + as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#overfit, remove TotalSpinesKinks
redlmbythocountenv16<-lm(logDNACpuLLno ~ NoIntro+ThreeBarbLiveGravid+ThreeBarbDeadGravid+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -83
bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ NoIntro + ThreeBarbLiveGravid + 
                               ThreeBarbDeadGravid + as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#overfit, remove number introduced
redlmbythocountenv16<-lm(logDNACpuLLno ~ ThreeBarbLiveGravid+ThreeBarbDeadGravid+
                           as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -82
bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ ThreeBarbLiveGravid + ThreeBarbDeadGravid + 
                               as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#overfit, remove ThreeBarbDeadGravid
redlmbythocountenv16<-lm(logDNACpuLLno ~ ThreeBarbLiveGravid+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -79
bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ ThreeBarbLiveGravid + as.factor(Year) + 
                               TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#overfit, remove temp
redlmbythocountenv16<-lm(logDNACpuLLno ~ ThreeBarbLiveGravid+as.factor(Year)+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16, direction="both")
#AIC -78
bestlmbythocountenvred16<-lm(formula = logDNACpuLLno ~ ThreeBarbLiveGravid + as.factor(Year), 
                             data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16)
#3barblivegravid significant
#Run again without transformation
bestlmbythocountenvred16ut<-lm(formula = DNACpuLLno ~ ThreeBarbLiveGravid + as.factor(Year), 
                             data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvred16ut)


#Now try DW
fulllmbythobm16<-lm(logDNACpuLLno ~ Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                                     ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+OneBarbSpineDW+TwoBarbSpineDW+
                                     ThreeBarbSpineDW+kinkDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                  data = bytho2qelongnoctlastnar)
stepAIC(fulllmbythobm16, direction="both")
#AIC -96

bestlmbythobm16<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + TwoBarbLiveGravidDW + 
                    ThreeBarbLiveGravidDW + ThreeBarbDeadGravidDW + OneBarbSpineDW + 
                    TwoBarbSpineDW + ThreeBarbSpineDW + kinkDW + TotalMoultsRound + 
                    as.factor(Year) + TempC + Barb:StartLivingBiomass, data = bytho2qelongnoctlastnar)
summary(bestlmbythobm16)
#rsquared shows overfit, remove interaction
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                      ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                      TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+OneBarbSpineDW+TwoBarbSpineDW+
                      ThreeBarbSpineDW+kinkDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                    data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -96
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ Barb + ThreeBarbLiveGravidDW + OneBarbDeadDWBarren + 
                         ThreeBarbDeadGravidDW + TwoBarbSpineDW + ThreeBarbSpineDW + 
                         kinkDW + TotalMoultsRound + as.factor(Year) + TempC + pH + 
                         DOmgL, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove spine and kink DWs
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                     ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -91
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + ThreeBarbLiveGravidDW + 
                         OneBarbDeadDWBarren + ThreeBarbDeadGravidDW + TotalMoultsRound + 
                         as.factor(Year) + TempC + pH + DOmgL, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove total moults round
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                     ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+as.factor(Year)+TempC+pH+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -91
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + TwoBarbLiveGravidDW + 
                         ThreeBarbLiveGravidDW + OneBarbDeadDWBarren + ThreeBarbDeadGravidDW + 
                         as.factor(Year) + TempC + pH + DOmgL, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove TwoBarbLiveGravidDW
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+StartLivingBiomass+OneBarbLiveBarrenDW+ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+as.factor(Year)+TempC+pH+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -91
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + ThreeBarbLiveGravidDW + 
                         OneBarbDeadDWBarren + ThreeBarbDeadGravidDW + as.factor(Year) + 
                         TempC + pH + DOmgL, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove OneBarbDeadDWBarren
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+StartLivingBiomass+OneBarbLiveBarrenDW+ThreeBarbLiveGravidDW+
                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+as.factor(Year)+TempC+pH+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -91
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + OneBarbLiveBarrenDW + 
                         ThreeBarbLiveGravidDW + ThreeBarbDeadGravidDW + as.factor(Year) + 
                         TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove OneBarbLiveBarrenDW
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+StartLivingBiomass+ThreeBarbLiveGravidDW+
                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+as.factor(Year)+TempC+pH+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -86
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ Barb + StartLivingBiomass + ThreeBarbLiveGravidDW + 
                         ThreeBarbDeadGravidDW + as.factor(Year) + TempC + pH, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove StartLivingBiomass
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+ThreeBarbLiveGravidDW+
                     TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+as.factor(Year)+TempC+pH+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -84
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ Barb + ThreeBarbLiveGravidDW + ThreeBarbDeadGravidDW + 
                         TempC + pH, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove ThreeBarbDeadGravidDW
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+ThreeBarbLiveGravidDW+
                     TwoBarbDeadGravidDW+as.factor(Year)+TempC+pH+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -79
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ ThreeBarbLiveGravidDW + as.factor(Year) + 
                         TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove temperature
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+ThreeBarbLiveGravidDW+
                     TwoBarbDeadGravidDW+as.factor(Year)+pH+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -78
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ Barb + ThreeBarbLiveGravidDW + pH, 
                       data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#overfit, remove pH
redlmbythobm16<-lm(logDNACpuLLno ~ Barb+ThreeBarbLiveGravidDW+
                     TwoBarbDeadGravidDW+as.factor(Year)+DOmgL,
                   data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16, direction="both")
#AIC -78
bestlmbythobmred16<-lm(formula = logDNACpuLLno ~ ThreeBarbLiveGravidDW + as.factor(Year), 
                       data = bytho2qelongnoctlastnar)
summary(bestlmbythobmred16)
#ThreeBarbLiveGravidDW significantly positive, but very small



#Try again with pooled demographic variables ()
fulllmbythocountenvgb16<-lm(logDNACpuLLno ~ Barb*NoIntro+TotalGravidLive+TotalGravidDead+TotalBarrenLive+
                                             OneBarbDeadBarren+TotalSpinesKinks+TotalMoultsRound+
                                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(fulllmbythocountenvgb16, direction="both")
#AIC=-89

bestlmbythocountenvgb16<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TotalGravidLive + 
                              TotalGravidDead + TotalBarrenLive + as.factor(Year) + TempC, 
                            data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvgb16)
#extremely overfit - , take out number introduced
redlmbythocountenvgb16<-lm(logDNACpuLLno ~ Barb+TotalGravidLive+TotalGravidDead+TotalBarrenLive+
                              OneBarbDeadBarren+TotalSpinesKinks+TotalMoultsRound+
                              as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvgb16, direction="both")
#AIC=-89

bestlmbythocountenvgbred16<-lm(formula = logDNACpuLLno ~ Barb + TotalGravidLive + TotalGravidDead + 
                              TotalBarrenLive + as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvgbred16)
#overfit, remove totalbarrenlive
redlmbythocountenvgb16<-lm(logDNACpuLLno ~ Barb+TotalGravidLive+TotalGravidDead+
                             OneBarbDeadBarren+TotalSpinesKinks+TotalMoultsRound+
                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvgb16, direction="both")
#AIC=-85

bestlmbythocountenvgbred16<-lm(formula = logDNACpuLLno ~ Barb + TotalGravidLive + TotalGravidDead + 
                                 TotalSpinesKinks + as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvgbred16)
#overfit, remove TotalSpinesKinks
redlmbythocountenvgb16<-lm(logDNACpuLLno ~ Barb+TotalGravidLive+TotalGravidDead+
                             OneBarbDeadBarren+TotalMoultsRound+
                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvgb16, direction="both")
#AIC=-85

bestlmbythocountenvgbred16<-lm(formula = logDNACpuLLno ~ Barb + TotalGravidLive + TotalGravidDead + 
                                 as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvgbred16)
#overfit, remove barb
redlmbythocountenvgb16<-lm(logDNACpuLLno ~ TotalGravidLive+TotalGravidDead+
                             OneBarbDeadBarren+TotalMoultsRound+
                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvgb16, direction="both")
#AIC=-85

bestlmbythocountenvgbred16<-lm(formula = logDNACpuLLno ~ TotalGravidLive + TotalGravidDead + 
                                 as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvgbred16)
#overfit, remove totalgraviddead
redlmbythocountenvgb16<-lm(logDNACpuLLno ~ TotalGravidLive+OneBarbDeadBarren+TotalMoultsRound+
                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvgb16, direction="both")
#AIC=-78
bestlmbythocountenvgbred16<-lm(formula = logDNACpuLLno ~ TotalGravidLive + TotalMoultsRound + 
                                 as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvgbred16)
#overfit, remove totalmoultsround
redlmbythocountenvgb16<-lm(logDNACpuLLno ~ TotalGravidLive+OneBarbDeadBarren+
                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvgb16, direction="both")
#AIC=-78
bestlmbythocountenvgbred16<-lm(formula = logDNACpuLLno ~ TotalGravidLive + as.factor(Year) + 
                                 TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvgbred16)
#overfit, remove temp
redlmbythocountenvgb16<-lm(logDNACpuLLno ~ TotalGravidLive+OneBarbDeadBarren+
                             as.factor(Year)+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvgb16, direction="both")
#AIC=-76
bestlmbythocountenvgbred16<-lm(formula = logDNACpuLLno ~ TotalGravidLive + as.factor(Year), 
                               data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvgbred16)
#TotalGravidLive significant


#Pooled demographic variables and biomass
fulllmbythobmenvgb16<-lm(logDNACpuLLno ~ Barb*NoIntro+GravidLiveBiomass+GravidDeadBiomass+BarrenLiveBiomass+
                                          OneBarbDeadDWBarren+TotalMoultsRound+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(fulllmbythobmenvgb16, direction="both")
#AIC=-93
bestlmbythobmenvgb16<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + GravidLiveBiomass + 
                           GravidDeadBiomass + BarrenLiveBiomass + as.factor(Year) + 
                           TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvgb16)
#overfit, take out barb
redlmbythobmenvgb16<-lm(logDNACpuLLno ~ NoIntro+GravidLiveBiomass+GravidDeadBiomass+BarrenLiveBiomass+
                           OneBarbDeadDWBarren+TotalMoultsRound+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvgb16, direction="both")
#AIC=-92
bestlmbythobmenvgbred16<-lm(formula = logDNACpuLLno ~ NoIntro + GravidLiveBiomass + GravidDeadBiomass + 
                              BarrenLiveBiomass + as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvgbred16)
#overfit, take out introduced
redlmbythobmenvgb16<-lm(logDNACpuLLno ~ GravidLiveBiomass+GravidDeadBiomass+BarrenLiveBiomass+
                          OneBarbDeadDWBarren+TotalMoultsRound+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvgb16, direction="both")
#AIC=-89
bestlmbythobmenvgbred16<-lm(formula = logDNACpuLLno ~ GravidLiveBiomass + GravidDeadBiomass + 
                              BarrenLiveBiomass + as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvgbred16)
#overfit, take out BarrenLiveBiomass
redlmbythobmenvgb16<-lm(logDNACpuLLno ~ GravidLiveBiomass+GravidDeadBiomass+
                          OneBarbDeadDWBarren+TotalMoultsRound+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvgb16, direction="both")
#AIC=-86
bestlmbythobmenvgbred16<-lm(formula = logDNACpuLLno ~ GravidLiveBiomass + GravidDeadBiomass + 
                              as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvgbred16)
#overfit, take out graviddeadbiomass
redlmbythobmenvgb16<-lm(logDNACpuLLno ~ GravidLiveBiomass+
                          OneBarbDeadDWBarren+TotalMoultsRound+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvgb16, direction="both")
#AIC=-82
bestlmbythobmenvgbred16<-lm(formula = logDNACpuLLno ~ GravidLiveBiomass + TotalMoultsRound + 
                              as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvgbred16)
#overfit, take out TotalMoultsRound
redlmbythobmenvgb16<-lm(logDNACpuLLno ~ GravidLiveBiomass+
                          OneBarbDeadDWBarren+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvgb16, direction="both")
#AIC=-80
bestlmbythobmenvgbred16<-lm(formula = logDNACpuLLno ~ GravidLiveBiomass + as.factor(Year) + 
                              TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvgbred16)
#overfit, take out temp
redlmbythobmenvgb16<-lm(logDNACpuLLno ~ GravidLiveBiomass+
                          OneBarbDeadDWBarren+as.factor(Year)+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvgb16, direction="both")
#AIC=-77
bestlmbythobmenvgbred16<-lm(formula = logDNACpuLLno ~ GravidLiveBiomass + as.factor(Year), 
                            data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvgbred16)
#GravidLiveBiomass significatn

#Try with pooled counts
fulllmbythocountenvld16<-lm(logDNACpuLLno ~ Barb*NoIntro+TotalLive+TotalDead+TotalSpinesKinks+TotalMoultsRound+as.factor(Year)+TempC+
                            DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(fulllmbythocountenvld16, direction="both")
#AIC=-108
bestlmbythocountenvld16<-lm(formula = logDNACpuLLno ~ Barb + NoIntro + TotalLive + TotalDead + 
                              as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvld16)
#extremely overfit - take out number introduced
redlmbythocountenvld16<-lm(logDNACpuLLno ~ Barb+TotalLive+TotalDead+TotalSpinesKinks+TotalMoultsRound+as.factor(Year)+TempC+
                              DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvld16, direction="both")
#AIC=-83
bestlmbythocountenvldred16<-lm(formula = logDNACpuLLno ~ Barb + TotalLive + TotalDead + as.factor(Year) + 
                              TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvldred16)
#overfit, take out barb
redlmbythocountenvld16<-lm(logDNACpuLLno ~ TotalLive+TotalDead+TotalSpinesKinks+TotalMoultsRound+as.factor(Year)+TempC+
                              DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvld16, direction="both")
#AIC=-82
bestlmbythocountenvldred16<-lm(formula = logDNACpuLLno ~ TotalLive + TotalDead + as.factor(Year) + 
                                 TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvldred16)
#overfit, take out total dead
redlmbythocountenvld16<-lm(logDNACpuLLno ~ TotalLive+TotalSpinesKinks+TotalMoultsRound+as.factor(Year)+TempC+
                             DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvld16, direction="both")
#AIC=-80
bestlmbythocountenvldred16<-lm(formula = logDNACpuLLno ~ TotalLive + as.factor(Year) + TempC, 
                               data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvldred16)
#overfit, take out year
redlmbythocountenvld16<-lm(logDNACpuLLno ~ TotalLive+TotalSpinesKinks+TotalMoultsRound+TempC+
                             DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenvld16, direction="both")
#AIC=-74
bestlmbythocountenvldred16<-lm(formula = logDNACpuLLno ~ TotalLive, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenvldred16)
#TotalLivesignificant



#now with biomass
fulllmbythobmenvld16<-lm(logDNACpuLLno ~ Barb*StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(fulllmbythobmenvld16, direction="both")
#AIC=-87
bestlmbythobmenvld16<-lm(formula = logDNACpuLLno ~ Barb + EndLivingBiomass + EndDeadBiomass + 
                         as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvld16)
#overfit, take out barb
redlmbythobmenvld16<-lm(logDNACpuLLno ~ StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvld16, direction="both")
#AIC=-86
bestlmbythobmenvldred16<-lm(formula = logDNACpuLLno ~ EndLivingBiomass + EndDeadBiomass + 
                              as.factor(Year) + TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvldred16)
#overfit, take out eaddeadbiomass
redlmbythobmenvld16<-lm(logDNACpuLLno ~ StartLivingBiomass+EndLivingBiomass+as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvld16, direction="both")
#AIC=-81
bestlmbythobmenvldred16<-lm(formula = logDNACpuLLno ~ EndLivingBiomass + as.factor(Year) + 
                              TempC, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvldred16)
#overfit, take out year
redlmbythobmenvld16<-lm(logDNACpuLLno ~ StartLivingBiomass+EndLivingBiomass+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobmenvld16, direction="both")
#AIC=-76
bestlmbythobmenvldred16<-lm(formula = logDNACpuLLno ~ EndLivingBiomass, data = bytho2qelongnoctlastnar)
summary(bestlmbythobmenvldred16)
#end living biomass significant but very small

###summary
#count full demographic splits - 3 barb live gravid
#biomass full demographic splits - 3 Barb Live Gravid
#count total gravid splits - total gravid live
#biomass total gravid splits - gravidlive significant
#count combined - total live significant
#biomass combined - end living biomass significant

#visualize

ggplot(subset(bytho2qelongnoct, Day==16), aes(x=ThreeBarbLiveGravid, y=DNACpuLLno)) + 
  geom_point() +
  geom_smooth(method=lm, se=F, color="black")+
  scale_y_continuous(trans="log10")+
  labs(x ="End Living, Gravid 3rd Instars", 
       y = "End eDNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

ggplot(subset(bytho2qelongnoct, Day==16), aes(x=TotalGravidLive, y=DNACpuLLno)) + 
  geom_point() +
  scale_y_continuous(trans="log10")+
  labs(x ="End Live Gravid", 
       y = "End DNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

ggplot(subset(bytho2qelongnoct, Day==16), aes(x=TotalLive, y=DNACpuLLno)) + 
  geom_point() +
  scale_y_continuous(trans="log10")+
  labs(x ="End Live", 
       y = "End DNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

#Do again with detection as response
fulllmbythocountenv16bi<-lm(DNADetect ~ Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                            ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                            ThreeBarbDeadGravid+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+TotalMoultsRound+
                            as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(fulllmbythocountenv16bi, direction="both")
#AIC -47
bestlmbythocountenv16bi<-lm(formula = DNADetect ~ TwoBarbLiveGravid + ThreeBarbLiveGravid + 
                              TwoBarbDeadGravid + ThreeBarbDeadGravid + OneBarbSpine + 
                              TwoBarbSpine + ThreeBarbSpine + kink + TotalMoultsRound + 
                              as.factor(Year), data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenv16bi)
#definitely overfit. Try again combinig spines and kinks
redlmbythocountenv16bi<-lm(DNADetect ~ Barb*NoIntro+OneBarbLiveBarren+TwoBarbLiveGravid+
                              ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                              ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+
                              as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16bi, direction="both")
#AIC -46
bestlmbythocountenv16bired<-lm(formula = DNADetect ~ Barb + TwoBarbLiveGravid + ThreeBarbDeadGravid + 
                                 TotalSpinesKinks + TotalMoultsRound + as.factor(Year) + TempC, 
                               data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenv16bired)
#overfit, remove TwoBarbLiveGravid
redlmbythocountenv16bi<-lm(DNADetect ~ Barb*NoIntro+OneBarbLiveBarren+
                             ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                             ThreeBarbDeadGravid+TotalSpinesKinks+TotalMoultsRound+
                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16bi, direction="both")
#AIC -46
bestlmbythocountenv16bired<-lm(formula = DNADetect ~ Barb + ThreeBarbDeadGravid + TotalSpinesKinks + 
                                 TotalMoultsRound, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenv16bired)
#overfit, remove totalspineskinsk
redlmbythocountenv16bi<-lm(DNADetect ~ Barb*NoIntro+OneBarbLiveBarren+
                             ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                             ThreeBarbDeadGravid+TotalMoultsRound+
                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16bi, direction="both")
#AIC -44
bestlmbythocountenv16bired<-lm(formula = DNADetect ~ Barb + NoIntro + ThreeBarbDeadGravid + 
                                 TotalMoultsRound, data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenv16bired)
#overfit, remove number introduced
redlmbythocountenv16bi<-lm(DNADetect ~ Barb+OneBarbLiveBarren+
                             ThreeBarbLiveGravid+OneBarbDeadBarren+TwoBarbDeadGravid+
                             ThreeBarbDeadGravid+TotalMoultsRound+
                             as.factor(Year)+TempC+DOmgL,data = bytho2qelongnoctlastnar)
stepAIC(redlmbythocountenv16bi, direction="both")
#AIC -42
bestlmbythocountenv16bired<-lm(formula = DNADetect ~ Barb + ThreeBarbDeadGravid + TotalMoultsRound, 
                               data = bytho2qelongnoctlastnar)
summary(bestlmbythocountenv16bired)
#nothing significant


#Now try DW
fulllmbythobm16bi<-lm(DNADetect ~ Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                      ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                      TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+OneBarbSpineDW+TwoBarbSpineDW+
                      ThreeBarbSpineDW+kinkDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                    data = bytho2qelongnoctlastnar)
stepAIC(fulllmbythobm16bi, direction="both")
#AIC -47
bestlmbythobm16bi<-lm(formula = DNADetect ~ TwoBarbLiveGravidDW + ThreeBarbLiveGravidDW + 
                        TwoBarbDeadGravidDW + ThreeBarbDeadGravidDW + OneBarbSpineDW + 
                        TwoBarbSpineDW + ThreeBarbSpineDW + kinkDW + TotalMoultsRound + 
                        as.factor(Year), data = bytho2qelongnoctlastnar)
summary(bestlmbythobm16bi)
#overfit, take olut spine and kink deadweights
redlmbythobm16bi<-lm(DNADetect ~ Barb*StartLivingBiomass+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                        ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                        TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                      data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16bi, direction="both")
#AIC -45
bestlmbythobm16bired<-lm(formula = DNADetect ~ Barb + StartLivingBiomass + ThreeBarbLiveGravidDW + 
                           ThreeBarbDeadGravidDW + TotalMoultsRound + pH, data = bytho2qelongnoctlastnar)
summary(bestlmbythobm16bired)
#overfit, remove startlivingbiomass
redlmbythobm16bi<-lm(DNADetect ~ Barb+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                       ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                       TwoBarbDeadGravidDW+ThreeBarbDeadGravidDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                     data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16bi, direction="both")
#AIC -43
bestlmbythobm16bired<-lm(formula = DNADetect ~ Barb + ThreeBarbLiveGravidDW + ThreeBarbDeadGravidDW + 
                           TotalMoultsRound + TempC + pH, data = bytho2qelongnoctlastnar)
summary(bestlmbythobm16bired)
#overfit, remove 3barbdeadgraviddw
redlmbythobm16bi<-lm(DNADetect ~ Barb+OneBarbLiveBarrenDW+TwoBarbLiveGravidDW+
                       ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                       TwoBarbDeadGravidDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                     data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16bi, direction="both")
#AIC -45
bestlmbythobm16bired<-lm(formula = DNADetect ~ Barb + TwoBarbLiveGravidDW + OneBarbDeadDWBarren + 
                           TempC + pH, data = bytho2qelongnoctlastnar)
summary(bestlmbythobm16bired)
#overfit, remove twobarblivegraviddw
redlmbythobm16bi<-lm(DNADetect ~ Barb+OneBarbLiveBarrenDW+ThreeBarbLiveGravidDW+OneBarbDeadDWBarren+
                       TwoBarbDeadGravidDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                     data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16bi, direction="both")
#AIC -43
bestlmbythobm16bired<-lm(formula = DNADetect ~ Barb + ThreeBarbLiveGravidDW + OneBarbDeadDWBarren + 
                           TotalMoultsRound + pH, data = bytho2qelongnoctlastnar)
summary(bestlmbythobm16bired)
#overfit, remove onebarbdeaddwbarren
redlmbythobm16bi<-lm(DNADetect ~ Barb+OneBarbLiveBarrenDW+ThreeBarbLiveGravidDW+
                       TwoBarbDeadGravidDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                     data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16bi, direction="both")
#AIC -42
bestlmbythobm16bired<-lm(formula = DNADetect ~ Barb + ThreeBarbLiveGravidDW + TotalMoultsRound + 
                           pH, data = bytho2qelongnoctlastnar)
summary(bestlmbythobm16bired)
#overfit, remove threebarblivegraviddw
redlmbythobm16bi<-lm(DNADetect ~ Barb+OneBarbLiveBarrenDW+
                       TwoBarbDeadGravidDW+TotalMoultsRound+as.factor(Year)+TempC+pH+DOmgL,
                     data = bytho2qelongnoctlastnar)
stepAIC(redlmbythobm16bi, direction="both")
#AIC -43
bestlmbythobm16bired<-lm(formula = DNADetect ~ Barb + pH, data = bytho2qelongnoctlastnar)
summary(bestlmbythobm16bired)
#nothing significant



#####################
#summary of results
#count demographics - nothing significant


