#Analysis of spiny water flea eDNA in mesocosm experiments
#Datasets generated from physical enumeration and qPCR analysis

#########
#libraries
library(ggplot2)
library(plyr)


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
bythoqPCR$DNACopiesperuLnc<-(10^(11-(0.26*(bythoqPCR$CtConverted))))
#Standardize based on dilution factor of extracted DNA
bythoqPCR$DNACopiesperuL<-bythoqPCR$DNACopiesperuLnc*bythoqPCR$DilutionFactor
bythoqPCR$DNACopiesperuL[is.na(bythoqPCR$DNACopiesperuL)]<-0
bythoqPCRmeans<-aggregate(DNACopiesperuL ~ Day + Tankno + Year, data=bythoqPCR, mean)
#Set detection limit (Dil15)
bythoqPCRmeans$DNACopiesperuL[bythoqPCRmeans$DNACopiesperuL<0.0008] <- 0

#Create wide dataset based on day
bythoqPCRwide<-reshape(bythoqPCRmeans, idvar = c("Tankno","Year"), timevar = "Day", direction = "wide")
#create total eDNA conc variable
bythoqPCRwide$totDNAcopies<-rowSums(bythoqPCRwide[,3:11],na.rm=T)

#merge physical,environmental, and qPCR datasets long
bythoenvqPCRlong<-merge(bythoenv,bythoqPCRmeans,all=T)
bytho2122qPCRlong<-merge(bytho2122,bythoqPCRmeans, all=T)
bytho1qelong<-merge(bythoenvqPCRlong,bytho2122)
#remove negative control tanks
bytho2qelongnoct<-subset(bytho1qelong,BarbNoIntro!="None")

#merge environmental, physical and qPCR datasets wide
bytho2122qPCRwide<-merge(bytho2122,bythoqPCRwide)
bythoenvqPCRwide<-merge(bythoenvwide,bythoqPCRwide)
bytho2qewide<-merge(bythoenvqPCRwide,bytho2122)
#remove negative control tanks
bytho2qewidenoct<-subset(bytho2qewide,BarbNoIntro!="None")

#####################################################
#examine negative control tanks

options(scipen = 999)
ggplot(bytho2qelongnoct, aes(x=Day, y=DNACopiesperuL)) + 
  geom_jitter(size=2, width=.1, height=.1) +
  scale_y_continuous(trans="log10")+
  labs(x ="Day", 
       y = "DNA Copies per µL")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

#only two hits - could be because of swf in ambient water or contamination. ignore for now

#Examine correlations between physical datasets and eDNA signal

bythocorrelations<-cor(bytho2qewidenoct[,c(3:30,34:61,66)],use="pairwise.complete.obs")
#correlations greater than |0.8| between morph and eDNA
#one barb live and copies on day 6 0.95 correlation
#copies on day 3 and one barb dead 0.99
#copies on day 3 and total dead 0.88

#plot these correlations
ggplot(bytho2122qPCRwidenoct, aes(x=OneBarbLive, y=DNACopiesperuL.6)) + 
  geom_point() +
  labs(x ="Density of Live, 1-barb individuals at draining", 
       y = "DNA Copies per µL collected Day 6")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))
#only 0's and 1's, so not a confident correlation
ggplot(bytho2122qPCRwidenoct, aes(x=OneBarbDead, y=DNACopiesperuL.3)) + 
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



