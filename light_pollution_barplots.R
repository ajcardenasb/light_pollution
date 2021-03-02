library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/Coral_endoliths/16S/ASVs/")

#############################################################
#####Taxonomic profiles of the 20 most abundant families#####
#############################################################

asv = read.table("outputs/light_pollution_ASVs_noContanoOut.raw.txt", header = T, row.names = 1)
map = read.table("Inputs/Final_metadata.txt", header = T, sep = "\t", row.names = 1)

names(asv)
asv.tax.ag=aggregate(asv[, 1:95], by = list(asv[, 100]), FUN =  sum) #define sample range and group factor
#topFamilies=asv.tax.ag[order(rowSums(asv.tax.ag[, 2:ncol(asv.tax.ag)]),decreasing = TRUE),][1:20,1]
topFamilies=asv.tax.ag[order(rowSums(asv.tax.ag[,2:ncol(asv.tax.ag)]),decreasing = TRUE),][1:10,1] # top only in skeleton control samples 
fam.top=subset(asv.tax.ag, asv.tax.ag$Group.1 %in% topFamilies) 
fam.bot=subset(asv.tax.ag, !asv.tax.ag$Group.1 %in% topFamilies) 
fam.bot$Group.1=gsub(".*","zOthers", fam.bot$Group.1)
others=aggregate(fam.bot[, 2:ncol(fam.bot)], by = list(fam.bot[, 1]), FUN =  sum)
all.2 =rbind(fam.top, others)
all.l=reshape2::melt(all.2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")

## Add sample information
all.l$Site=map$Site[match(all.l$Sample, rownames(map))]
all.l$Season=map$Season[match(all.l$Sample, rownames(map))]
all.l$Time=map$Time[match(all.l$Sample, rownames(map))]
all.l$Moon=map$moon.phase[match(all.l$Sample, rownames(map))]
all.l$Sample=rownames(map)[match(all.l$Sample, rownames(map))]

ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Family), data = subset(all.l, Season == "Spring" ), stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony", main = "Spring") + scale_fill_manual(values=P10) + facet_wrap(~Site+Time+Moon, ncol=4, scales="free") + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))
ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Family), data = subset(all.l, Season == "Summer" ), stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony", main = "Spring") + scale_fill_manual(values=P10) + facet_wrap(~Site+Time+Moon, ncol=4, scales="free") + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))
ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Family), data = subset(all.l, Season == "Fall" ), stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony", main = "Spring") + scale_fill_manual(values=P10) + facet_wrap(~Site+Time+Moon, ncol=4, scales="free") + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))
ggplot() +geom_bar(aes(y = Abundance, x = Sample, fill = Family), data = subset(all.l, Season == "Winter" ), stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony", main = "Spring") + scale_fill_manual(values=P10) + facet_wrap(~Site+Time+Moon, ncol=4, scales="free") + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))


## Plot
P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
P10=c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A", "#66A61E", "#E6AB02" ,"#780116","#A6761D", "#2364aa", "#3da5d9", "#ababab")
final=all.l %>% group_by(Site, Season,Time, Moon, Sample) %>% summarise(Abundance=sum(Abundance))

pdf("./outputs/mean_barplots_Goniastrea_16S.pdf",  width = 7, height =3, pointsize = 12) 
ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Family), data = final, stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + scale_fill_manual(values=P10) + facet_grid(~Tissue) + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))
dev.off() 

# pdf("./outputs/mean_barplots_Goniastrea_16S.pdf",  width = 7, height =7, pointsize = 12) 
# ggplot() +geom_bar(aes(y = Abundance, x = Species, fill = Family), data = final, stat="identity", position = "fill") +  theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + scale_fill_manual(values=P21) + facet_grid(~Tissue) + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom') +  theme_classic() + guides(fill=guide_legend(ncol=1))
# dev.off() 
