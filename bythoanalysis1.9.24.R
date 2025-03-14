#bytho experiment graphs

#load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(MASS) 
library(reshape2) 
library(reshape) 

#programs
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#upload datasets
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
#DNA
bythoqPCR<-read.csv("PrelimqPCR1.9.24bytho.csv")

#visualize just morphological
melt_bytho2122<- melt(bytho2122, id = c("Tankno","DrainDay", "Year", "NoIntro","Barb","BarbNoIntro","Order")) 
melt_bytho2122$variable<-revalue(melt_bytho2122$variable, c("TotalLive"="Living",
                                                                "TotalDead"="Dead",
                                                                "TotalSpinesKinks"="Spines/Kinks"))

ggplot(subset(melt_bytho2122, variable=="Living" | 
                variable=="Dead" | variable=="Spines/Kinks"),
       aes(x=as.factor(Order), y=value, fill=variable)) + 
  geom_bar(position="stack", stat="identity") +
  labs(y = "Individuals",fill="End\nDemographics")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=20),
        axis.text.x=element_blank(),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro,scales="free_x")

ggplot(subset(melt_bytho2122, (variable=="Living" | 
                variable=="Dead" | variable=="Spines/Kinks") & BarbNoIntro!="None"),
       aes(x=as.factor(Order), y=value, fill=variable)) + 
  geom_bar(position="stack", stat="identity") +
  labs(y =bquote('SWF Individuals per'~m^3))+
  scale_fill_manual(values=c("#1b9e77", "#d95f02","#7570b3"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=20),
        axis.text.x=element_blank(),axis.text.y = element_text(size=14),
        legend.title=element_blank(),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro,scales="free_x")


meltse_bytho2122<-summarySE(subset(melt_bytho2122, variable=="Living" | 
                                     variable=="Dead" | variable=="Spines/Kinks"), 
                            measurevar=c("value"), groupvars=c("Year","NoIntro",
                                                               "Barb","BarbNoIntro", "variable"), na.rm=TRUE)
ggplot(meltse_bytho2122, aes(x=variable, y=value, color=variable)) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.4) + 
  labs(y = "Individuals",color="End\nDemographics")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_blank(),axis.title.y=element_text(size=20),
        axis.text.x=element_blank(),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)


#Standardize qPCR dataset
bythoqPCR$CtConverted<-bythoqPCR$Ct-bythoqPCR$Conversion
#max 53.38
bythoqPCR$CtConverted[is.na(bythoqPCR$CtConverted)] <- 54

bythoqPCRmeans<-aggregate(CtConverted ~ Day + Tankno + Year + DilutionFactor, data=bythoqPCR, mean)
bythoqPCRmeans$DNACopiesperuLnc<-(10^(11-(0.26*(bythoqPCRmeans$CtConverted))))
bythoqPCRmeans$DNACopiesperuL<-bythoqPCRmeans$DNACopiesperuLnc*bythoqPCRmeans$DilutionFactor
#Set detection limit
bythoqPCRmeans$DNACopiesperuL[bythoqPCRmeans$DNACopiesperuL<0.0008] <- 0
#merge datasets
bytho2122qPCR<-merge(bytho2122,bythoqPCRmeans)

options(scipen = 999)
ggplot(subset(bytho2122qPCR,BarbNoIntro!="None"), aes(x=Day, y=DNACopiesperuL)) + 
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

ggplot(subset(bytho2122qPCR,BarbNoIntro!="None"), aes(x=Day, y=DNACopiesperuL, color=IMOLiving)) + 
  geom_jitter(size=2, width=.1, height=.1) +
  scale_y_continuous(trans="log10")+
  labs(x ="Day", 
       y = "DNA Copies per µL", color="IMO\nThreshold")+
  scale_color_manual(values=c("#990033","#E69F00","#336600"))+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

ggplot(subset(bytho2122qPCR,BarbNoIntro!="None"), aes(x=Day, y=DNACopiesperuL, color=IMOTotal)) + 
  geom_jitter(size=2, width=.1, height=.1) +
  scale_y_continuous(trans="log10")+
  labs(x ="Day", 
       y = "DNACopies per µL")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)


qPCRSum<-summarySE(bythoqPCR, measurevar=c("RelDNA"), groupvars=c("Tankno","Day","Year"), na.rm=TRUE)

#merge datasets
bytho2122qPCRsum<-merge(bytho2122,qPCRSum)


ggplot(subset(bytho2122qPCRsum, BarbNoIntro!="None"), aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  labs(x ="Day", 
       y = "relative DNA Density", color="Total End\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

bytho2122qPCRsum$log10

ggplot(subset(bytho2122qPCRsum, BarbNoIntro!="None"), aes(x=Day, y=RelDNA, color=TotalBythoIndiv)) + 
  geom_point(size=2) +
  labs(x ="Day", 
       y = "relative DNA Density", color="Total End\nIndividuals")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

ggplot(subset(bytho2122qPCRsum, BarbNoIntro!="None"), aes(x=Day, y=RelDNA, color=EndDeadBiomass)) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Dead\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

ggplot(subset(bytho2122qPCRsum, BarbNoIntro!="None"), aes(x=Day, y=RelDNA, color=EndLivingBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Living\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

ggplot(subset(bytho2122qPCRsum, BarbNoIntro!="None"), aes(x=Day, y=RelDNA, color=as.factor(Tankno))) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relative DNA Density", color="Tank Code")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~BarbNoIntro)

#subset to just 2021 1 2 barb added
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="1, 2-barb" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relative DNA Density", color="End Total\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

#subset to just 2021 1 3 barb added
#total
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 1" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Total\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#Live at end
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 1" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndLivingBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Living\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#Dead at end
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 1" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndDeadBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Dead\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

#subset to just 2021 5 3 barb added
#total
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 5" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Total\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#live
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 5" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndLivingBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Living\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#dead
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 5" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndDeadBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Dead\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

#subset to just 2021 10 3 barb added
#total
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 10" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Total\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#live
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 10" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndLivingBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Living\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#dead
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 10" & Year==2021),
       aes(x=Day, y=RelDNA, color=EndDeadBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Dead\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

#subset to just 2022 1 2 barb added
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="2barb 1" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Total\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="2barb 1" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndLivingBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Living\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="2barb 1" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndDeadBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Dead\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

#subset to just 2022 1 3 barb added
#total
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 1" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Total\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#Live at end
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 1" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndLivingBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Living\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#Dead at end
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 1" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndDeadBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Dead\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

#subset to just 2022 5 2 barb added
#total
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="2barb 5" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Total\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#live
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="2barb 5" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndLivingBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Living\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#dead
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="2barb 5" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndDeadBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Dead\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

#subset to just 2022 5 3 barb added
#total
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 5" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndTotalBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Total\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#live
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 5" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndLivingBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Living\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)
#dead
ggplot(subset(bytho2122qPCRsum,BarbNoIntro=="3barb 5" & Year==2022),
       aes(x=Day, y=RelDNA, color=EndDeadBiomass)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  geom_line()+
  labs(x ="Day", 
       y = "relativeDNADensity", color="End Dead\nBiomass")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(.~Tankno)

ggplot(subset(bytho2122qPCRsum,Day==16 & BarbNoIntro!="None"),
       aes(x=EndTotalBiomass, y=RelDNA, color=BarbNoIntro)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  labs(x ="End Total Biomass", 
       y = "End relative DNA Density", color="Treatment")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~.)

ggplot(subset(bytho2122qPCRsum,Day==16 & BarbNoIntro!="None"),
       aes(x=EndDeadBiomass, y=RelDNA, color=BarbNoIntro)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  labs(x ="End Dead Biomass", 
       y = "End relative DNA Density", color="Treatment")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~.)

ggplot(subset(bytho2122qPCRsum,Day==16 & BarbNoIntro!="None"),
       aes(x=EndLivingBiomass, y=RelDNA, color=BarbNoIntro)) + 
  geom_errorbar(aes(ymin=RelDNA-se, ymax=RelDNA+se), width=.4) +  
  geom_point(size=2) +
  labs(x ="End Living Biomass", 
       y = "End relative DNA Density", color="Treatment")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=14),axis.text.y = element_text(size=14),
        legend.title=element_text(size=20),legend.text = element_text(size=16),
        strip.text.x = element_text(size = 20),strip.text.y=element_text(size = 20))+
  facet_grid(Year~.)