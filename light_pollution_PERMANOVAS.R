############################################################
##################### phyloseq #############################
############################################################
library(phyloseq)
library(ggplot2)
library(plyr)
library(gridExtra)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/light_pollution_R/")

asv = read.table("outputs/light_pollution_ASVs_noContanoOut.raw.txt", header = T, row.names = 1)
met = read.table("Inputs/Final_metadata.txt", header = T, sep = "\t", row.names = 1)

otu.t= otu_table(asv[, 1:95], taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(met))
tax.t= tax_table(as.matrix(asv[, 97:ncol(asv)]))

phy.all= phyloseq(otu.t, tax.t,  sam.t)
#P3=c("#e36600ff", "#008000ff", "#ccccccff")

phy.t=microbiome::transform(phy.all, transform = "clr", target = "OTU", shift = 0, scale = 1)
PCA = ordinate(phy.t, method = "RDA", distance = "euclidean")
plot_ordination(phy.t,PCA, color = "Season", shape = "Site") + geom_point(size = 4, alpha = 1) + theme_bw()  + ggtitle("Seasons + Site") +  theme_classic() #+ scale_colour_manual(values=P4)
plot_ordination(phy.t,PCA, color = "Site", shape = "moon.phase") + geom_point(size = 4, alpha = 1) + theme_bw()  + ggtitle("Moon + Site") +  theme_classic() #+ scale_colour_manual(values=P4)


summer_phy=subset_samples(phy.t, Season == "Summer")
summer_PCA = ordinate(summer_phy, method = "RDA", distance = "euclidean")
su=plot_ordination(summer_phy,summer_PCA, color = "Site", shape = "moon.phase") + geom_point(size = 4, alpha = 1) + theme_bw()  + ggtitle("Summer") +  theme_classic() #+ scale_colour_manual(values=P4)

spring_phy=subset_samples(phy.t, Season == "Spring")
spring_PCA = ordinate(spring_phy, method = "RDA", distance = "euclidean")
sp=plot_ordination(spring_phy,spring_PCA, color = "Site", shape = "moon.phase") + geom_point(size = 4, alpha = 1) + theme_bw()  + ggtitle("Spring") +  theme_classic() #+ scale_colour_manual(values=P4)

winter_phy=subset_samples(phy.t, Season == "Winter")
winter_PCA = ordinate(winter_phy, method = "RDA", distance = "euclidean")
wi=plot_ordination(winter_phy,winter_PCA, color = "Site", shape = "moon.phase") + geom_point(size = 4, alpha = 1) + theme_bw()  + ggtitle("Winter") +  theme_classic() #+ scale_colour_manual(values=P4)

fall_phy=subset_samples(phy.t, Season == "Fall")
fall_PCA = ordinate(fall_phy, method = "RDA", distance = "euclidean")
fa=plot_ordination(fall_phy,fall_PCA, color = "Site", shape = "moon.phase") + geom_point(size = 4, alpha = 1) + theme_bw()  + ggtitle("Fall") +  theme_classic() #+ scale_colour_manual(values=P4)


pdf("./outputs/skeleton16S_ordination.pdf", width=8,height=3, pointsize = 12)
grid.arrange(sp,su,fa,wi, ncol=2, nrow=2)
dev.off()


#####################################################
######## Stats on community composition ############
####################################################
library(vegan)
library(pairwiseAdonis)
library(compositions)
#asv.n=as.data.frame(t(sweep(asv[, 1:95],2,colSums(asv[, 1:95]),`/`))) # relative abundances
asv.n=as.data.frame(t(apply(asv[, 1:95],2,clr))) # clr-transformed counts
asv.n$Season=met$Season[match(rownames(asv.n), rownames(met))]
asv.n$Moon=met$moon.phase[match(rownames(asv.n), rownames(met))]
asv.n$Time=met$Time[match(rownames(asv.n), rownames(met))]
asv.n$Site=met$Site[match(rownames(asv.n), rownames(met))]

##overal model
asv_adonis=adonis(asv.n[,1:18897]~ asv.n$Site + asv.n$Season + asv.n$Moon + asv.n$Time , method = "euclidean") # 0.6666667
asv_adonis_df=as.data.frame(asv_adonis[["aov.tab"]])
asv_adonis_df
write.table(asv_adonis_df, "outputs/endoliths_overal_adonis_ASVs.txt", sep = "\t", row.names = T, quote = F)

## plot
asv_adonis_df$Factor=gsub(".*\\$", "", rownames(asv_adonis_df))
asv_adonis_df$Factor=factor(asv_adonis_df$Factor, levels= c("Time", "Moon", "Season","Site" ))
ggplot() +geom_bar(aes(y = R2, x = reorder(Factor,R2)), data = asv_adonis_df[1:4,], stat="identity", fill="black") +
  labs(x="", y="PERMANOVA R2") + 
  coord_flip() + 
  theme_classic() + 
  theme(legend.position = "none") + 
  annotate("text", x = 1, y = 0.014, label = "*") + 
  annotate("text", x =2, y = 0.016, label = "**") + 
  annotate("text", x = 3, y = 0.041, label = "***") + 
  annotate("text", x =4, y = 0.085, label = "***")  

###################################
## Comparions suggested by Yaeli ##
###################################

## Comparison 1: between sites ##

adonis(asv.n[,1:18897]~ asv.n$Site, method = "euclidean" )

## Comparison 2: IUI between diel cycles ##
IUI_samples=subset(asv.n,Site == "IUI")
adonis(IUI_samples[,1:18897]~ IUI_samples$Time, method = "euclidean" )

## Comparison 3: Kisosk between diel cycles ##
Kisosk_samples=subset(asv.n,Site == "Kisosk")
adonis(Kisosk_samples[,1:18897]~ Kisosk_samples$Time, method = "euclidean" )

## Comparison 4: night IUI between moons ##
IUI_samples_night=subset(asv.n,Site == "IUI" & Time == "night" )
adonis(IUI_samples_night[,1:18897]~ IUI_samples_night$Moon, method = "euclidean")

## Comparison 5: night Kisosk between moons ##
Kisosk_samples_night=subset(asv.n,Site == "Kisosk" & Time == "night"  )
adonis(Kisosk_samples_night[,1:18897]~Kisosk_samples_night$Moon, method = "euclidean")

#####################################
## Comparions by seasons and sites ##
#####################################

#comparing diel cycles within each season and site

summer_diel=subset(asv.n,  Season == "Summer") 
summer_diel$groups=paste(summer_diel$Site, summer_diel$Time)
summer_diel_df=pairwise.adonis(summer_diel[,1:18897], summer_diel$groups,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) # 0.6666667
summer_diel_df

winter_diel=subset(asv.n,  Season == "Winter") 
winter_diel$groups=paste(winter_diel$Site, winter_diel$Time)
winter_diel_df=pairwise.adonis(winter_diel[,1:18897], winter_diel$groups,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) # 0.6666667
winter_diel_df

spring_diel=subset(asv.n,  Season == "Spring") 
spring_diel$groups=paste(spring_diel$Site, spring_diel$Time)
spring_diel_df=pairwise.adonis(spring_diel[,1:18897], spring_diel$groups,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) # 0.6666667
spring_diel_df

fall_diel=subset(asv.n,  Season == "Fall") 
fall_diel$groups=paste(fall_diel$Site, fall_diel$Time)
fall_diel_df=pairwise.adonis(fall_diel[,1:18897], fall_diel$groups,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) # 0.6666667
fall_diel_df

#comparing moon cycles within each season and site

summer_diel=subset(asv.n,  Season == "Summer") 
summer_diel$groups=paste(summer_diel$Site, summer_diel$Moon)
summer_diel_df=pairwise.adonis(summer_diel[,1:18897], summer_diel$groups,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) # 0.6666667
summer_diel_df

winter_diel=subset(asv.n,  Season == "Winter") 
winter_diel$groups=paste(winter_diel$Site, winter_diel$Moon)
winter_diel_df=pairwise.adonis(winter_diel[,1:18897], winter_diel$groups,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) # 0.6666667
winter_diel_df

spring_diel=subset(asv.n,  Season == "Spring") 
spring_diel$groups=paste(spring_diel$Site, spring_diel$Moon)
spring_diel_df=pairwise.adonis(spring_diel[,1:18897], spring_diel$groups,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) # 0.6666667
spring_diel_df

fall_diel=subset(asv.n,  Season == "Fall") 
fall_diel$groups=paste(fall_diel$Site, fall_diel$Moon)
fall_diel_df=pairwise.adonis(fall_diel[,1:18897], fall_diel$groups,  sim.method = "euclidean", p.adjust.m = "fdr", perm = 999) # 0.6666667
fall_diel_df
