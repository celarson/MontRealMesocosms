#Analysis of spiny water flea eDNA in mesocosm experiments
#Datasets generated from physical enumeration and qPCR analysis

#########
#libraries
library(ggplot2)
library(plyr)
library(lme4)
library(sjPlot)


#################################
#Upload and standardize datasets

#Environmental datasets
bythoenv<-read.csv("EnvVarsBytho.csv")
#Create wide dataset based on day
bythoenvwide<-reshape(bythoenv, idvar = c("Tankno","Year"), timevar = "Day", direction = "wide")

#physical
bytho2122<-read.csv("2122PhysicalBythoResults.csv")
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
bytho2122$BarbNoIntro<-revalue(bytho2122$BarbNoIntro, c("NA 0"="None", "2barb 1"="1, 2-barb", "3barb 1"="1, 3-barb",
                                                        "2barb 5"="5, 2-barbs", "3barb 5"="5, 3-barbs", 
                                                        "3barb 10"="10, 3-barbs"))
#create ordered barbnointro variable
bytho2122$BarbNoIntroOrdered<-as.numeric(revalue(bytho2122$BarbNoIntro, c("None"=0, "1, 2-barb"=1, "1, 3-barb"=2, "5, 2-barbs"=3,
                                                               "5, 3-barbs"=4,"10, 3-barbs"=5)))
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

#####################################################
#examine negative control tanks

options(scipen = 999)
ggplot(bytho2qelongnoct, aes(x=Day, y=DNACopiesperuLL, color=as.factor(Tankno))) + 
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

bythocorrelations<-cor(bytho2qewidenoct[,c(3:31,33,35:62,67)],use="pairwise.complete.obs")
#correlations greater than |0.8| between morph and eDNA
#one barb live and copies on day 6 0.95 correlation
#copies on day 3 and one barb dead 0.99
#copies on day 3 and total dead 0.88

#plot these correlations
ggplot(bytho2qewidenoct, aes(x=OneBarbLive, y=DNACopiesperuLL.6)) + 
  geom_point() +
  labs(x ="Density of Live, 1-barb individuals at draining", 
       y = "DNA Copies per µL collected Day 6")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))
#only 0's and 1's, so not a confident correlation
ggplot(bytho2qewidenoct, aes(x=OneBarbDead, y=DNACopiesperuLL.3)) + 
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
boxplot.stats(bytho2qelongnoct$DNACpuLLno)$out
#fixed the worst one
boxplot.stats(log10(bytho2qelongnoct$DNACpuLLno+1))$out
#log transforming helps
bytho2qelongnoct$logDNACpuLLno<-log10(bytho2qelongnoct$DNACpuLLno+1)
#binary variable
bytho2qelongnoct$DNADetect<- as.numeric(bytho2qelongnoct$DNACpuLL>0)
bytho2qelongnoct2$DNADetect<- as.numeric(bytho2qelongnoct2$DNACpuLL>0)


#TempC
shapiro.test(bytho2qelongnoct$TempC)
#normal 0.29

#CondmScm
shapiro.test(bytho2qelongnoct$CondmScm)
#not normal 
hist(bytho2qelongnoct$CondmScm)
#binomial, very different for 2 experiments

#pH
shapiro.test(bytho2qelongnoct$pH)
#normal 0.45

#TurbidityFNU
shapiro.test(bytho2qelongnoct$TurbidityFNU)
#not normal
hist(bytho2qelongnoct$TurbidityFNU)
#several values very high
shapiro.test(log10(bytho2qelongnoct$TurbidityFNU))
hist(log10(bytho2qelongnoct$TurbidityFNU))
#log10 transformation helps, but still not normal
bytho2qelongnoct$logTurbidityFNU<-log10(bytho2qelongnoct$TurbidityFNU)

#DOmgL
shapiro.test(bytho2qelongnoct$DOmgL)
#notnormal
hist(bytho2qelongnoct$DOmgL)
#right skewed, so don't transform

#TotalChlorugL
shapiro.test(bytho2qelongnoct$TotalChlorugL)
#not normal
hist(bytho2qelongnoct$TotalChlorugL)
#binomial distribution between two experiments

#examine correlations with day as a variable (i.e., long)
bythocorrelationslong<-cor(bytho2qelongnoct[,c(3:10,12,14:41,46,48:50)],use="pairwise.complete.obs")
#correlations over |0.8|

#CondmScm TotalChlorugL -0.94862350
#TurbidityFNU TotalChlorugL 0.89165369
#Didn't check correlations among bytho physical parameters
#keep total chlorophyl and remove cond and turb

#build linear mixed effects model
names(bytho2qelongnoct)
lmebythoDNAcountfull<- lmer(logDNACpuLLno ~ TempC+pH+DOmgL+TotalChlorugL+BarbNoIntro+OneBarbLive+TwoBarbLive+ThreeBarbLive+
                          OneBarbDead+TwoBarbDead+ThreeBarbDead+OneBarbSpine+TwoBarbSpine+ThreeBarbSpine+kink+(1 | Day)+(1|Year),
              data = bytho2qelongnoct)
#singular, start smaller, also only 76 observations because of NA's
lmebythoDNAcount<- lmer(logDNACpuLLno ~ BarbNoIntro+(1 | Day)+(1|Year),
                            data = bytho2qelongnoct)
#still singular, so try a linear model instead
lmbythoDNAcount<- lm(logDNACpuLLno ~ BarbNoIntro+Day+Year,
                            data = bytho2qelongnoct)
summary(lmbythoDNAcount)
lmbythoDNAbiomass<- lm(logDNACpuLLno ~ StartLivingBiomass+Day+Year,
                     data = bytho2qelongnoct)
summary(lmbythoDNAbiomass)
#full, but can't include env variables because of missing values
lmbythoDNAbiomassfull<- lm(logDNACpuLLno ~ StartLivingBiomass+Day+Year+OneBarbLiveDW+TwoBarbLiveDW+ThreeBarbLiveDW+
                             OneBarbDeadDW+TwoBarbDeadDW+ThreeBarbDeadDW+OneBarbSpineDW+TwoBarbSpine+ThreeBarbSpineDW+kinkDW+TotalMoultsRound,
                       data = bytho2qelongnoct)
summary(lmbythoDNAbiomassfull)
#day significant
#end
lmbythoDNAbiomassend<- lm(logDNACpuLLno ~ StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+TotalMoultsRound+Day+Year,
                           data = bytho2qelongnoct)
summary(lmbythoDNAbiomassend)
#day and start living biomass significant
#try interactions
lmbythoDNAbiomassendint<- lm(logDNACpuLLno ~ (StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+Day+Year)^2,
                          data = bytho2qelongnoct)
summary(lmbythoDNAbiomassendint)
#end living biomass and day interaction

#try again with counts
#end
lmbythoDNAcountend<- lm(logDNACpuLLno ~ BarbNoIntro+TotalLive+TotalDead+TotalSpinesKinks+Day+Year,
                          data = bytho2qelongnoct)
summary(lmbythoDNAcountend)
#day significant
#try interactions
lmbythoDNAcountendint<- lm(logDNACpuLLno ~ (BarbNoIntro+TotalLive+TotalDead+TotalSpinesKinks+Day+Year)^2,
                        data = bytho2qelongnoct)
summary(lmbythoDNAcountendint)
#nothing, biomass is better predictor

#visualize biomass results

#day, start living biomass, end living biomass
ggplot(bytho2qelongnoct, aes(x=Day, y=DNACpuLLno, color=log10(EndLivingBiomass+1))) + 
  geom_point() +
  scale_y_continuous(trans="log10")+
  labs(x ="Day", 
       y = "DNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~BarbNoIntro)

ggplot(bytho2qelongnoct, aes(x=EndLivingBiomass, y=DNACpuLLno)) + 
  geom_point() +
  scale_y_continuous(trans="log10")+
  labs(x ="End Living Biomass", 
       y = "DNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

#Day 0 may be skewing results. Take out first day of sampling, then run again

lmbythoDNAcount2<- lm(logDNACpuLLno ~ BarbNoIntro+Day+Year,
                     data = bytho2qelongnoct2)
summary(lmbythoDNAcount2)
#daynointro 53 and 103 and day
lmbythoDNAbiomass2<- lm(logDNACpuLLno ~ StartLivingBiomass+Day+Year,
                       data = bytho2qelongnoct2)
summary(lmbythoDNAbiomass2)
#startlivingbiomass and day
#full, but can't include env variables because of missing values
lmbythoDNAbiomassfull2<- lm(logDNACpuLLno ~ StartLivingBiomass+Day+Year+OneBarbLiveDW+TwoBarbLiveDW+ThreeBarbLiveDW+
                             OneBarbDeadDW+TwoBarbDeadDW+ThreeBarbDeadDW+OneBarbSpineDW+TwoBarbSpine+ThreeBarbSpineDW+kinkDW,
                           data = bytho2qelongnoct2)
summary(lmbythoDNAbiomassfull2)
#day significant
#end
lmbythoDNAbiomassend2<- lm(logDNACpuLLno ~ StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+Day+Year,
                          data = bytho2qelongnoct2)
summary(lmbythoDNAbiomassend2)
#day and start living biomass significant
#try interactions
lmbythoDNAbiomassendint2<- lm(logDNACpuLLno ~ (StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+Day+Year)^2,
                             data = bytho2qelongnoct2)
summary(lmbythoDNAbiomassendint2)
#end living biomass and day interaction

#try again with counts
#end
lmbythoDNAcountend2<- lm(logDNACpuLLno ~ BarbNoIntro+TotalLive+TotalDead+TotalSpinesKinks+Day+Year,
                        data = bytho2qelongnoct2)
summary(lmbythoDNAcountend2)
#day significant
#try interactions
lmbythoDNAcountendint2<- lm(logDNACpuLLno ~ (BarbNoIntro+TotalLive+TotalDead+TotalSpinesKinks+Day+Year)^2,
                           data = bytho2qelongnoct2)
summary(lmbythoDNAcountendint2)
#nothing, biomass is better predictor

#visualize biomass results

#day, start living biomass, end living biomass


ggplot(bytho2qelongnoct, aes(x=StartLivingBiomass, y=DNACpuLLno)) + 
  geom_jitter() +
  scale_y_continuous(trans="log10")+
  labs(x ="Start Living Biomass", 
       y = "DNA Copies per µL per L")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

#Try again with logistic regression
#remove NA's
bytho2qelongnoct2num<-subset(bytho2qelongnoct2,DNADetect!="NA")
lmbythoDNAcount2log<- glm(DNADetect ~ BarbNoIntro+Day+Year,
                      data = bytho2qelongnoct2num,family="binomial")
summary(lmbythoDNAcount2log)
#daynointro 13 and 103
lmbythoDNAbiomass2log<- glm(DNADetect ~ StartLivingBiomass+Day+Year,
                        data = bytho2qelongnoct2num,family="binomial")
summary(lmbythoDNAbiomass2log)
#nothing sig
#full, but can't include env variables because of missing values
lmbythoDNAbiomassfull2log<- glm(DNADetect ~ StartLivingBiomass+Day+Year+OneBarbLiveDW+TwoBarbLiveDW+ThreeBarbLiveDW+
                              OneBarbDeadDW+TwoBarbDeadDW+ThreeBarbDeadDW+OneBarbSpineDW+TwoBarbSpine+ThreeBarbSpineDW+kinkDW,
                            data = bytho2qelongnoct2num, family="binomial")
summary(lmbythoDNAbiomassfull2log)
#3barb live DW, onebarbspinedw, 2 barb spine dw, model not significant
#end
lmbythoDNAbiomassend2log<- glm(DNADetect ~ StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+TotalMoultsRound+Day+Year,
                           data = bytho2qelongnoct2num, family="binomial")
summary(lmbythoDNAbiomassend2log)
#end living biomass significant, model not significant, try reducing variables
lmbythoDNAbiomassend2logred<- glm(DNADetect ~ EndLivingBiomass,
                              data = bytho2qelongnoct2num, family="binomial")
summary(lmbythoDNAbiomassend2logred)
#try interactions
lmbythoDNAbiomassendint2log<- glm(DNADetect ~ (StartLivingBiomass+EndLivingBiomass+EndDeadBiomass+Day+Year)^2,
                              data = bytho2qelongnoct2num, family="binomial")
summary(lmbythoDNAbiomassendint2log)
#no significant

#plot model
plot_model(lmbythoDNAbiomassend2logred,
           type = "pred",
           terms = "EndLivingBiomass"
) +
  labs(y = "Probability of DNA Detection")

ggplot(bytho2qelongnoct, aes(x=EndLivingBiomass, y=DNADetect)) + 
  geom_point() +
  labs(x ="End Living Biomass", 
       y = "DNA Detection")+
  geom_smooth(method="glm",family="binomial")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))

