#Analysis of protist 16S metabarcoded dataset from mesocosm experiments

#Libraries
library(vegan)

#Color vectors
ExpYearcolvec<-c("#1b9e77","#7570b3","#66a61e","#e7298a","#d95f02")

#Upload datasets
#Beta divesity
unifrac16S<-read.csv("PMV3V4WeightedUnifrac.csv", header=T, check.names=F)

#Upload ASV table
ASV16S<-read.csv("PMV3V4UMGC_asv_table.csv", header=T, check.names=F)
names(ASV16S)
#format data table so OTU name is row name
row.names(ASV16S)<-ASV16S[,1]
#Delete OTU id column now that OTU ID is rowname
ASV16S$'sample-id'<-NULL
#reverse so columns are samples
ASV16St<-data.frame(t(ASV16S), check.names=F)
ASV16St$sample.id<-row.names(ASV16St)

#Upload genus table
g16S<-read.csv("PMV3V4UMGC_asv_table_g.csv", header=T, check.names=F)
names(g16S)
#format data table so OTU name is row name
row.names(g16S)<-g16S[,1]
#Delete OTU id column now that OTU ID is rowname
g16S$'sample-id'<-NULL
#reverse so columns are samples
g16St<-data.frame(t(g16S), check.names=F)
g16St$sample.id<-row.names(g16St)

#Upload inoculant densities
Innoculation<-read.csv("PR_Inoculant_densities.csv")

#merge tables and innoculation densities
#find differences between the two sample ID names
g16St$sample.id[!(g16St$sample.id %in% Innoculation$sample.id)]
#Blank11, Blank 3, Blank 5, Blank 7, Blank 9, D7Control2Exp32022, MockCommunity
Innoculation$sample.id[!(Innoculation$sample.id %in% g16St$sample.id)]
#"D0MTP7Exp12021"   "D29MTP2Exp12022"  "D7MTP18Exp12022"  "D28MTP6Exp32022"  "D0MTP8Exp32022"   "D28MTP21Exp32022"

#merve unifrac distances and sample IDs
uni16SInnoc<-merge(unifrac16S, Innoculation, by="sample.id")
row.names(uni16SInnoc)<-uni16SInnoc[,1]
uni16SInnoc<-uni16SInnoc[,-c(1)]
names(uni16SInnoc)
PM_16S_uni<-as.matrix(uni16SInnoc[,c(1:495)])
#UNI-Create overall environmental data matrix for community analysis with uni distances
uni16SInnoc_env<-uni16SInnoc[,c(496:ncol(uni16SInnoc))]

#merge genus level table and sample ID's
g16SInnoc<-merge(g16St,Innoculation,by="sample.id")
#look for correlations

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
adonis2(as.dist(PM_16S_uni) ~ (Day+ExpYear+treatment_cellsperml+Species)^2+TempC+Salinity_PSU+pH+Turbidity_FNU+TChlor_ug.L, by="terms", data=uni16SInnoc_env,
        permutations=999)
#Day, ExpYear, Species, Salinity, pH, Turbidity_FNU, TChlor_ug.L, DayxExpYear, DayxSpecies, ExpYearxspecies significant

#Visualize via nmds
PM_16S_NMDS<-metaMDS(as.dist(PM_16S_uni))
#stress 0.3

#Stressplot macroinvertebrate Nmds
stressplot(PM_16S_NMDS)

#NMDS plot 
PM_16S_env<-envfit(PM_16S_NMDS, uni16SInnoc_env[,c(13,15:17,19)], na.rm=T)
plot(PM_16S_NMDS, type="n")
with(PM_16S_NMDS, points(PM_16S_NMDS, display="sites", col=ExpYearcolvec[as.factor(uni16SInnoc_env$ExpYear)], pch=19))
plot(PM_16S_env, col="black", lwd=8, cex=1.5,arr.width=3)
with(PM_16S_NMDS, legend("topleft", legend=levels(as.factor(uni16SInnoc_env$ExpYear)), bty="n", col=ExpYearcolvec,
                         pch=19, pt.bg=ExpYearcolvec))


