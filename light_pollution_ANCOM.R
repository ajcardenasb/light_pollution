library(ANCOMBC)
library(phyloseq)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/light_pollution_R/")

asv = read.table("outputs/light_pollution_ASVs_noContanoOut.raw.txt", header = T, row.names = 1)
met = read.table("Inputs/Final_metadata.txt", header = T, sep = "\t", row.names = 1)

otu.t= otu_table(asv[, 1:95], taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(met))
tax.t= tax_table(as.matrix(asv[, 97:ncol(asv)]))

phy.all= phyloseq(otu.t, tax.t,  sam.t)


#################################
## Comparison 1: between sites ##
#################################

res1=ancombc(phyloseq=phy.all,formula="Site",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Site",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_df=data.frame(Beta=res1[["res"]][["beta"]], Beta_se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])
colnames(res1_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_sig=subset(res1_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_sig$Diff_more_abundant=ifelse( res1_sig$W < 0 , "IUI", "Kisosk")
res1_sig$ASV=rownames(res1_sig)
res1_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "-")[match(rownames(res1_sig),rownames(asv))]
res1_sig$Comparison="Comparison 1: between sites"
message("Number of DA ASVs: ", nrow(res1_sig), "\nNumber of DA ASVs enriched in IUI: ", nrow(subset(res1_sig, Diff_more_abundant == "IUI" )), "\nNumber of DA ASVs enriched in Kisosk: ", nrow(subset(res1_sig, Diff_more_abundant == "Kisosk" )))
#write.table(res1_sig,  "outputs/ANCOMBC_ASVs_Comp1_betweenSites.txt", sep = "\t", quote = F, row.names = T )

###########################################
## Comparison 2: IUI between diel cycles ##
###########################################
IUI_phy=subset_samples(phy.all, Site == "IUI")
res2=ancombc(phyloseq=IUI_phy,formula="Time",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Time",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res2_df=data.frame(Beta=res2[["res"]][["beta"]], Beta_se=res2[["res"]][["se"]], W=res2[["res"]][["W"]],pval=res2[["res"]][["p_val"]], qval=res2[["res"]][["q_val"]],DA=res2[["res"]][["diff_abn"]])
colnames(res2_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res2_sig=subset(res2_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res2_sig$Diff_more_abundant=ifelse( res2_sig$W < 0 , "day", "night")
res2_sig$ASV=rownames(res2_sig)
res2_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "-")[match(rownames(res2_sig),rownames(asv))]
res2_sig$Comparison="Comparison 2: IUI between diel cycle"
message("Number of DA ASVs: ", nrow(res2_sig), "\nNumber of DA ASVs enriched in day: ", nrow(subset(res2_sig, Diff_more_abundant == "day" )), "\nNumber of DA ASVs enriched in night: ", nrow(subset(res2_sig, Diff_more_abundant == "night" )))
#write.table(res2_sig,  "outputs/ANCOMBC_ASVs_Comp2_IUIbetweenDiel.txt", sep = "\t", quote = F, row.names = T )

###########################################
## Comparison 3: IUI between diel cycles ##
###########################################
Kisosk_phy=subset_samples(phy.all, Site == "Kisosk")
res3=ancombc(phyloseq=Kisosk_phy,formula="Time",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Time",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res3_df=data.frame(Beta=res3[["res"]][["beta"]], Beta_se=res3[["res"]][["se"]], W=res3[["res"]][["W"]],pval=res3[["res"]][["p_val"]], qval=res3[["res"]][["q_val"]],DA=res3[["res"]][["diff_abn"]])
colnames(res3_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res3_sig=subset(res3_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res3_sig$Diff_more_abundant=ifelse( res3_sig$W < 0 , "day", "night")
res3_sig$ASV=rownames(res3_sig)
res3_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "-")[match(rownames(res3_sig),rownames(asv))]
res3_sig$Comparison="Comparison 3: Kisosk between diel cycle"
message("Number of DA ASVs: ", nrow(res3_sig), "\nNumber of DA ASVs enriched in day: ", nrow(subset(res3_sig, Diff_more_abundant == "day" )), "\nNumber of DA ASVs enriched in night: ", nrow(subset(res3_sig, Diff_more_abundant == "night" )))
#write.table(res3_sig,  "outputs/ANCOMBC_ASVs_Comp3_KisoskbetweenDiel.txt", sep = "\t", quote = F, row.names = T )

###########################################
## Comparison 4: IUI between diel cycles ##
###########################################
IUI_night_phy=subset_samples(phy.all, Site == "IUI" & Time == "night")
res4=ancombc(phyloseq=IUI_night_phy,formula="moon.phase",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "moon.phase",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res4_df=data.frame(Beta=res4[["res"]][["beta"]], Beta_se=res4[["res"]][["se"]], W=res4[["res"]][["W"]],pval=res4[["res"]][["p_val"]], qval=res4[["res"]][["q_val"]],DA=res4[["res"]][["diff_abn"]])
colnames(res4_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res4_sig=subset(res4_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res4_sig$Diff_more_abundant=ifelse( res4_sig$W < 0 , "full", "new")
res4_sig$ASV=rownames(res4_sig)
res4_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "-")[match(rownames(res4_sig),rownames(asv))]
res4_sig$Comparison="Comparison 4: night IUI between moons"
message("Number of DA ASVs: ", nrow(res4_sig), "\nNumber of DA ASVs enriched in full: ", nrow(subset(res4_sig, Diff_more_abundant == "full" )), "\nNumber of DA ASVs enriched in new: ", nrow(subset(res4_sig, Diff_more_abundant == "new" )))
#write.table(res4_sig,  "outputs/ANCOMBC_ASVs_Comp4_nightIUIbetweenMoons.txt", sep = "\t", quote = F, row.names = T )

###########################################
## Comparison 5: IUI between diel cycles ##
###########################################
Kisosk_night_phy=subset_samples(phy.all, Site == "Kisosk" & Time == "night")
res5=ancombc(phyloseq=Kisosk_phy,formula="moon.phase",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "moon.phase",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res5_df=data.frame(Beta=res5[["res"]][["beta"]], Beta_se=res5[["res"]][["se"]], W=res5[["res"]][["W"]],pval=res5[["res"]][["p_val"]], qval=res5[["res"]][["q_val"]],DA=res5[["res"]][["diff_abn"]])
colnames(res5_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res5_sig=subset(res5_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res5_sig$Diff_more_abundant=ifelse( res5_sig$W < 0 , "full", "new")
res5_sig$ASV=rownames(res5_sig)
res5_sig$Taxa=paste(asv$Phylum, asv$Family, asv$Genus,sep = "-")[match(rownames(res5_sig),rownames(asv))]
res5_sig$Comparison="Comparison 5: night Kisosk between moons"
message("Number of DA ASVs: ", nrow(res5_sig), "\nNumber of DA ASVs enriched in full: ", nrow(subset(res5_sig, Diff_more_abundant == "full" )), "\nNumber of DA ASVs enriched in new: ", nrow(subset(res5_sig, Diff_more_abundant == "new" )))
#write.table(res5_sig,  "outputs/ANCOMBC_ASVs_Comp5_nightKisosketweenMoons.txt", sep = "\t", quote = F, row.names = T )

ANCOMresults=rbind(res1_sig,res2_sig,res3_sig,res4_sig,res5_sig)
write.table(ANCOMresults,  "outputs/ANCOMBC_ASVs_results.txt", sep = "\t", quote = F, row.names = T )

library(reshape2)
library(ggplot2)

ANCOMresults_plot=ANCOMresults %>% group_by(Diff_more_abundant, Comparison) %>% tally()
ANCOMresults_plot$n=ifelse(ANCOMresults_plot$Diff_more_abundant %in% c("IUI", "day", "full"), ANCOMresults_plot$n*-1, ANCOMresults_plot$n*1)

pdf("./outputs/ANCOMBC_DA_barplots.pdf", width=6,height=4, pointsize = 12)
ggplot(data=ANCOMresults_plot, aes(x=Comparison, y=n)) + geom_bar(stat="identity", position = "dodge")  + geom_text(aes(label=n), vjust=0.5, color="white", position = position_dodge(0.8), size=3) +  theme_classic() + theme(axis.text.x=element_text(angle=90,hjust=1)) 
dev.off()

#data analysis
ANCOMresults$Family = asv$Family[match(ANCOMresults$ASV,rownames(asv))]
ANCOMresults$Class = asv$Class[match(ANCOMresults$ASV,rownames(asv))]
temp1=ANCOMresults %>% group_by(Diff_more_abundant, Comparison, Family, Class) %>% tally()
