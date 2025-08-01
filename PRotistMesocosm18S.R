#Analysis of protist 18S metabarcoded dataset from mesocosm experiments

#Libraries
library(vegan)
library(ggplot2)
library(tidyverse)
library(ggpubr)

#Color vectors
ExpYearcolvec<-c("#1b9e77","#7570b3","#66a61e","#e7298a","#d95f02")

#Upload datasets

#Beta divesity
unifrac18S<-read.csv("PM18SWeightedUnifrac.csv", header=T, check.names=F)

#Upload phylum table
p18S<-read.csv("PM18SUMGC_asv_table_p.csv", header=T, check.names=F)
names(p18S)
#format data table so OTU name is row name
row.names(p18S)<-p18S[,1]
#Delete OTU id column now that OTU ID is rowname
p18S$'sample-id'<-NULL
#reverse so columns are samples
p18St<-data.frame(t(p18S), check.names=F)
p18St$sample.id<-row.names(p18St)

#Upload genus table
g18S<-read.csv("PM18SUMGC_asv_table_g.csv", header=T, check.names=F)
names(g18S)
#format data table so OTU name is row name
row.names(g18S)<-g18S[,1]
#Delete OTU id column now that OTU ID is rowname
g18S$'sample-id'<-NULL
#reverse so columns are samples
g18St<-data.frame(t(g18S), check.names=F)
g18St$sample.id<-row.names(g18St)

#Upload species table
s18S<-read.csv("PM18SUMGC_asv_table_s.csv", header=T, check.names=F)
names(s18S)
#format data table so OTU name is row name
row.names(s18S)<-s18S[,1]
#Delete OTU id column now that OTU ID is rowname
s18S$'sample-id'<-NULL
#reverse so columns are samples
s18St<-data.frame(t(s18S), check.names=F)
s18St$sample.id<-row.names(s18St)

#Upload ASV table
a18S<-read.csv("PM18SUMGC_asv_table.csv", header=T, check.names=F)
names(a18S)
#format data table so OTU name is row name
row.names(a18S)<-a18S[,1]
#Delete OTU id column now that OTU ID is rowname
a18S$'sample-id'<-NULL
#reverse so columns are samples
a18St<-data.frame(t(a18S), check.names=F)
a18St$sample.id<-row.names(a18St)

#Upload inoculant densities
Innoculation<-read.csv("PR_Inoculant_densities.csv")

#Anything in the genus Haematococcus represents the surrogate for experiment 1 2021 and experiment 2 2022 (H pluvialis) (Chlorophyta)
#Use species trachelomonas abrupta Experiment 2 2021 (Trachelomonas abrupta) Euglenida algae
#experiment 1 2022 experiment 3 2022 (Chrysosphaerella sp) Ochrophytina

#merge unifrac distance with innoculant metadata
#Combine distance matrix with environmental variables
uni18SInnoc<-merge(unifrac18S, Innoculation, by="sample.id")
row.names(uni18SInnoc)<-uni18SInnoc[,1]
uni18SInnoc<-uni18SInnoc[,-c(1)]
names(uni18SInnoc)
PM_18S_uni<-as.matrix(uni18SInnoc[,c(1:493)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
uni18SInnoc_env<-uni18SInnoc[,c(494:ncol(uni18SInnoc))]

#merge ASV table and innoculation densities
#find differences between the two sample ID names
p18St$sample.id[!(p18St$sample.id %in% Innoculation$sample.id)]
#DOMC1
Innoculation$sample.id[!(Innoculation$sample.id %in% p18St$sample.id)]
#[1] "D14MTP2Exp12021"  "D14MTP16Exp12021" "D7MTP18Exp12021"  "D22MTP1Exp22021"  "D22MTP7Exp22021"  "D22MTP11Exp22021" "D21MTP7Exp12022"  "D29MTP17Exp22022" "D28MTP6Exp32022"

#Merge phylum level
p18SInnoc<-merge(p18St,Innoculation,by="sample.id")

#split up by experiment
p18SInnocExpYear<- split(p18SInnoc, p18SInnoc$ExpYear)
p18SInnoc12021 <- p18SInnocExpYear[[1]]
p18SInnoc12022 <- p18SInnocExpYear[[2]]
p18SInnoc22021 <- p18SInnocExpYear[[3]]
p18SInnoc22022 <- p18SInnocExpYear[[4]]
p18SInnoc32022 <- p18SInnocExpYear[[5]]

#Merge genus level
g18SInnoc<-merge(g18St,Innoculation,by="sample.id")

#split up by experiment
g18SInnocExpYear<- split(g18SInnoc, g18SInnoc$ExpYear)
g18SInnoc12021 <- g18SInnocExpYear[[1]]
g18SInnoc12022 <- g18SInnocExpYear[[2]]
g18SInnoc22021 <- g18SInnocExpYear[[3]]
g18SInnoc22022 <- g18SInnocExpYear[[4]]
g18SInnoc32022 <- g18SInnocExpYear[[5]]

#Merge species level
s18SInnoc<-merge(s18St,Innoculation,by="sample.id")

#split up by experiment
s18SInnocExpYear<- split(s18SInnoc, s18SInnoc$ExpYear)
s18SInnoc12021 <- s18SInnocExpYear[[1]]
s18SInnoc12022 <- s18SInnocExpYear[[2]]
s18SInnoc22021 <- s18SInnocExpYear[[3]]
s18SInnoc22022 <- s18SInnocExpYear[[4]]
s18SInnoc32022 <- s18SInnocExpYear[[5]]

#Merge ASV level
a18SInnoc<-merge(a18St,Innoculation,by="sample.id")

#split up by experiment
a18SInnocExpYear<- split(a18SInnoc, a18SInnoc$ExpYear)
a18SInnoc12021 <- a18SInnocExpYear[[1]]
a18SInnoc12022 <- a18SInnocExpYear[[2]]
a18SInnoc22021 <- a18SInnocExpYear[[3]]
a18SInnoc22022 <- a18SInnocExpYear[[4]]
a18SInnoc32022 <- a18SInnocExpYear[[5]]

########################
#Analysis

#look at correlations among env variables
cor(Innoculation[,c(4,9,14:20)])
#salinity and conductivity correlated 0.83, just keep salinity
#DO and temp correlated 0.85, keep temp

range(Innoculation$TempC)

range(Innoculation$Salinity_PSU)

range(Innoculation$pH)

range(Innoculation$Turbidity_FNU)

#################Takes a long time to run #############################
#UNI-Overall permanova with unifrac distances
adonis2(as.dist(PM_18S_uni) ~ (Day+ExpYear+treatment_cellsperml+Species)^2+TempC+Salinity_PSU+pH+Turbidity_FNU+TChlor_ug.L, by="terms", data=uni18SInnoc_env,
        permutations=999)
#Day, ExpYear, Species, TempC, pH, Turbidity_FNU, DayxExpYear, DayxSpecies, ExpYearxspecies significant

#Visualize via nmds
PM_18S_NMDS<-metaMDS(as.dist(PM_18S_uni))
#stress 0.3

#Stressplot macroinvertebrate Nmds
stressplot(PM_18S_NMDS)

#NMDS plot 
PM_18S_env<-envfit(PM_18S_NMDS, uni18SInnoc_env[,c(13,15:17,19)], na.rm=T)
plot(PM_18S_NMDS, type="n")
with(PM_18S_NMDS, points(PM_18S_NMDS, display="sites", col=ExpYearcolvec[as.factor(uni18SInnoc_env$ExpYear)], pch=19))
plot(PM_18S_env, col="black", lwd=8, cex=1.5,arr.width=3)
with(PM_18S_NMDS, legend("topleft", legend=levels(as.factor(uni18SInnoc_env$ExpYear)), bty="n", col=ExpYearcolvec,
                             pch=19, pt.bg=ExpYearcolvec))
############################################################################

#Generate innoculant density eDNA datasets

#look for correlations at the species level to see if they match innoculants

#experiment 1 2021
correlations18S12021s<-data.frame(cor(s18SInnoc12021[,c(2:575,584,587)],use="pairwise.complete.obs"))
#d__Eukaryota;p__Chlorophyta;c__Chlorophyceae;o__Chlorophyceae;f__Chlorophyceae;g__Chlorophyceae;s__Haematococcus_lacustris correlates with cells by 0.91 and relative by 0.76
#rename H. Lacustris
names(s18SInnoc12021)[names(s18SInnoc12021) == "d__Eukaryota;p__Chlorophyta;c__Chlorophyceae;o__Chlorophyceae;f__Chlorophyceae;g__Chlorophyceae;s__Haematococcus_lacustris"] <- "Hpluvialis"

#plot correlation
ggplot(s18SInnoc12021, aes(x=observed_treatment_cellperml, y=Hpluvialis)) + 
  geom_point()+
  stat_cor(method="pearson")+
  labs(x ="H. pluvialis observed cells", 
       y = "H. pluvialis relative sequence abundance")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=10),
        legend.title=element_text(size=24),legend.text = element_text(size=14),
        legend.position="bottom")

#Experiment 1 2022
correlations18S12022s<-data.frame(cor(s18SInnoc12022[,c(2:575,584,587)],use="pairwise.complete.obs"))
#d__Eukaryota;p__Ochrophyta;c__Chrysophyceae;o__Synurales;f__Synurales;g__Synura;__ .70 correlation with observed cells and .80 correlation with relative
#rename Synura
names(s18SInnoc12022)[names(s18SInnoc12022) == "d__Eukaryota;p__Ochrophyta;c__Chrysophyceae;o__Synurales;f__Synurales;g__Synura;__"] <- "Synura"
#plot correlation
ggplot(s18SInnoc12022, aes(x=observed_treatment_cellperml, y=Synura)) + 
  geom_point()+
  stat_cor(method="pearson")+
  labs(x ="Chrysosphaerella observed cells", 
       y = "Synura relative sequence abundance")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=10),
        legend.title=element_text(size=24),legend.text = element_text(size=14),
        legend.position="bottom")

#Experiment 2 2021
correlations18S22021s<-data.frame(cor(s18SInnoc22021[,c(2:575,584,587)],use="pairwise.complete.obs"))
#d__Eukaryota;p__Ochrophyta;c__Chrysophyceae;o__Synurales;f__Synurales;g__Mallomonas;s__Mallomonas_caudata 0.54 correlation with observed cells and 0.65 with relative
#See how trachelomonas tracks
cor.test(s18SInnoc22021$'d__Eukaryota;p__Euglenozoa;c__Euglenida;o__Euglenophyceae;f__Euglenaceae;g__Trachelomonas;s__Trachelomonas_sp.',s18SInnoc22021$observed_treatment_cellperml)
#0.18, significant
cor.test(s18SInnoc22021$'d__Eukaryota;p__Euglenozoa;c__Euglenida;o__Euglenophyceae;f__Euglenaceae;g__Trachelomonas;s__Trachelomonas_sp.',s18SInnoc22021$relative_inoculant_density_percent)
#0.17, not significant

#Experiment 2 2022
correlations18S22022s<-data.frame(cor(s18SInnoc22022[,c(2:575,584,587)],use="pairwise.complete.obs"))
#d__Eukaryota;p__Chlorophyta;c__Chlorophyceae;o__Chlorophyceae;f__Chlorophyceae;g__Chlorophyceae;s__Haematococcus_lacustris 0.67 correlation with observed cells and 0.74 correlation with relative
#rename H. Lacustris
names(s18SInnoc22022)[names(s18SInnoc22022) == "d__Eukaryota;p__Chlorophyta;c__Chlorophyceae;o__Chlorophyceae;f__Chlorophyceae;g__Chlorophyceae;s__Haematococcus_lacustris"] <- "Hpluvialis"
#plot correlation
ggplot(s18SInnoc22022, aes(x=observed_treatment_cellperml, y=Hpluvialis)) + 
  geom_point()+
  stat_cor(method="pearson")+
  labs(x ="H. pluvialis observed cells", 
       y = "H. pluvialis relative sequence abundance")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=10),
        legend.title=element_text(size=24),legend.text = element_text(size=14),
        legend.position="bottom")

#Experiment 3 2022
correlations18S32022s<-data.frame(cor(s18SInnoc32022[,c(2:575,584,587)],use="pairwise.complete.obs"))
#d__Eukaryota;p__Ochrophyta;c__Chrysophyceae;o__Synurales;f__Synurales;g__Synura;__ 0.94 correlation with observed and 0.95 correlation with relative
#rename Synura
names(s18SInnoc32022)[names(s18SInnoc32022) == "d__Eukaryota;p__Ochrophyta;c__Chrysophyceae;o__Synurales;f__Synurales;g__Synura;__"] <- "Synura"
#plot correlation
ggplot(s18SInnoc32022, aes(x=observed_treatment_cellperml, y=Synura)) + 
  geom_point()+
  stat_cor(method="pearson")+
  labs(x ="Chrysosphaerella observed cells", 
       y = "Synura relative sequence abundance")+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=10),axis.text.y = element_text(size=10),
        legend.title=element_text(size=24),legend.text = element_text(size=14),
        legend.position="bottom")



#look for correlations at the genus level to see if they match innoculants
correlations18S22021g<-data.frame(cor(g18SInnoc22021[,c(2:331,340,343)],use="pairwise.complete.obs"))
#d__Eukaryota;p__Diatomea;c__Coscinodiscophytina;o__Coscinodiscophytina;f__Melosirids;g__Aulacoseira 0.52 observed 0.49 relative

#look for correlations for ASV level
correlations18S22021ao<-t(data.frame(cor(a18SInnoc22021[,c(68400,68403)],a18SInnoc22021[,c(2:68391)],use="pairwise.complete.obs")))
#nothing correlations more strong than 0.5

#Subset excel sheet to create surrogate invaders eDNA spreadsheet for Abby and upload
SurrogateInvaderseDNAt<-read.csv("SurrogateInvaderseDNAt.csv")
eDNASurrInnoc<-merge(SurrogateInvaderseDNAt,Innoculation,by="sample.id")
write.csv(eDNASurrInnoc, 'SurrogateInvaderseDNAmeta.csv')
