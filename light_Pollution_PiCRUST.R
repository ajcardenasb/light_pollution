setwd("~/Documents/Bioinformatics_scripts/R_scripts/light_pollution_R/")

asv = read.table("outputs/light_pollution_ASVs_noContanoOut.raw.txt", header = T, row.names = 1)
met = read.table("Inputs/Final_metadata.txt", header = T, sep = "\t", row.names = 1)


############################## prepare dataset for picrust ##############################
library(seqinr)
#create fasta
write.fasta(as.list(asv$Sequence), rownames(asv), "outputs/light_pollution_ASV.fasta", open = "w", nbchar = 60, as.string = FALSE)
 
#export ASV table with no tax and seq data, ASVs hat did not pass the first filter on PICRUSt need to be removed
discarted=c('ASV18860',	'ASV13779',	'ASV18182',	'ASV18290',	'ASV18408',	'ASV18838',	'ASV18469',	'ASV18544',	'ASV17984',	'ASV11742',	'ASV18168',	'ASV18658',	'ASV4806',	'ASV18640',	'ASV4392',	'ASV16589',	'ASV18807',	'ASV18465',	'ASV18782',	'ASV16565',	'ASV18749',	'ASV18143',	'ASV17873',	'ASV18468',	'ASV10597',	'ASV15053',	'ASV18855',	'ASV18144',	'ASV17756',	'ASV17763',	'ASV18696',	'ASV18170',	'ASV5179',	'ASV18784',	'ASV18033',	'ASV17110',	'ASV15378',	'ASV17764',	'ASV18055',	'ASV18101',	'ASV18642',	'ASV18596',	'ASV17435',	'ASV14055',	'ASV16851',	'ASV13908',	'ASV18646',	'ASV10579',	'ASV18427',	'ASV18454',	'ASV18190',	'ASV18318',	'ASV9748',	'ASV3107',	'ASV1182',	'ASV18464',	'ASV18396',	'ASV16848',	'ASV18301',	'ASV18645',	'ASV15957',	'ASV8863',	'ASV12685',	'ASV18806',	'ASV18815',	'ASV16663',	'ASV7381',	'ASV6178',	'ASV13778',	'ASV8365',	'ASV15559',	'ASV18460',	'ASV16342',	'ASV18227',	'ASV18399')
asv2=asv[!rownames(asv) %in% discarted, 1:95] # need to be 18822
write.table(asv2, "outputs/light_pollution_ASV_picrust.txt", quote = F, row.names = T, sep = "\t")
# need to mannually add a name for the ASV ID column

####################################################
######### Analyze PICRUSt data #####################
####################################################

marker_NSTI=read.table("outputs/weighted_nsti.tsv",  header = T)
marker_NSTI$NSTI=(1-marker_NSTI$weighted_NSTI)*100 # convert to percentage
marker_NSTI_s1=subset(marker_NSTI, !NSTI < 0) #remove negative nsti 

pdf("outputs/light_pollution_picrust_NSTI.pdf", width = 7, height = 7, pointsize = 12)
ggplot(marker_NSTI_s1, aes(x=sample, y=NSTI)) + geom_bar(stat='identity',  width=1)  +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  geom_hline(yintercept = 86.5, linetype="dashed", color = "red", size=1) +
  geom_hline(yintercept = 94.5, linetype="dashed", color = "blue", size=1) + 
  labs(y="NSTI", x = "Samples") + coord_flip()
dev.off() 

round(mean(marker_NSTI_s1$NSTI),2)

##################
### ANCOM KEGG ###
##################

pred=read.table("outputs/pred_metagenome_unstrat.tsv", header = T, row.names = 1, sep = "\t")
kegg_meta=read.table("~/Documents/Bioinformatics_scripts/KEGG_files/MostUpdated_hierarchy_ko00001", sep = "\t", header = T, fill = T)#kegg=read.table("~/scripts/hierarchy_ko00001", sep = "\t", header = T, fill = T)
pred$Pathways=kegg_meta$MD[match(rownames(pred), kegg_meta$KO)]
pred_pw=aggregate(pred[,1:95], by=list(pred$Pathways), sum)
rownames(pred_pw)=pred_pw$Group.1
pred_pw=pred_pw[,-1]

met = read.table("Inputs/Final_metadata.txt", header = T, sep = "\t", row.names = 1)
met$Site=gsub("IUI", "Non-urban", met$Site)
met$Site=gsub("Kisosk", "Urban", met$Site)


#kegg_cfix=read.table("Inputs/Carbon_fix_metadata.txt", sep = "\t", header = T, fill = T)#kegg=read.table("~/scripts/hierarchy_ko00001", sep = "\t", header = T, fill = T)
#kegg_nfix=read.table("Inputs/Nitrogen_fix_metadata.txt", sep = "\t", header = T, fill = T)#kegg=read.table("~/scripts/hierarchy_ko00001", sep = "\t", header = T, fill = T)

otu.t= otu_table(pred_pw, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(met))

phy.all= phyloseq(otu.t,  sam.t)


res1=ancombc(phyloseq=phy.all,formula="Site",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Site",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_df=data.frame(Beta=res1[["res"]][["beta"]], Beta_se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])
colnames(res1_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_sig=subset(res1_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_sig$Diff_more_abundant=ifelse( res1_sig$W < 0 , "Non rban", "Urban")
res1_sig$Pathway=kegg_meta$module[match(rownames(res1_sig), kegg_meta$MD)]
res1_sig$L3=kegg_meta$L3[match(rownames(res1_sig), kegg_meta$MD)]
subset(res1_sig, rownames(res1_sig) %in% c("M00173", "M00374", "M00375" ,"M00376", "M00377" ,"M00579", "M00620", "M00175") )
#write.table(res1_sig,  "outputs/ANCOMBC_ASVs_Comp1_betweenSites.txt", sep = "\t", quote = F, row.names = T )


### Plotting
library(compositions)

clr=as.data.frame(t(apply(pred_pw,2,clr)))
clr$Site=met$Site[match(rownames(clr), rownames(met))]



clr_l=reshape2::melt(clr)
clr_l$pathway=kegg_meta$module[match(clr_l$variable,kegg_meta$MD )]
clr_l$label=paste(clr_l$pathway," (",clr_l$variable, ")", sep = "" )

clr_sub=subset(clr_l, variable %in% c("M00175","M00161", "M00163", "M00206", "M00222", "M00212"))
  
p1=ggplot(subset(clr_l, variable == "M00175"), aes(x=Site, y=value, fill=Site)) + geom_boxplot()  + theme_classic() + 
  annotate("text", x=1.5, y=3, label="q-val = 8.79 e-05") + 
  scale_fill_manual(values=c("#8EB7CE","#EE6363" )) + labs(title = "Nitrogen fixation", x="", y = "clr-trasformed counts")

p2=ggplot(subset(clr_l, variable == "M00161"), aes(x=Site, y=value, fill=Site)) + geom_boxplot()  + theme_classic() + 
  annotate("text", x=1.5, y=3.5, label="q-val = 6.25e-06") + 
  scale_fill_manual(values=c("#8EB7CE","#EE6363" )) + labs(title = "Photosystem I", x="", y = "clr-trasformed counts")

p3=ggplot(subset(clr_l, variable == "M00163"), aes(x=Site, y=value, fill=Site)) + geom_boxplot()  + theme_classic() + 
  annotate("text", x=1.5, y=3.2, label="q-val = 6.86e-06") + 
  scale_fill_manual(values=c("#8EB7CE","#EE6363" )) + labs(title = "Photosystem II", x="", y = "clr-trasformed counts")

p4=ggplot(subset(clr_l, variable == "M00206"), aes(x=Site, y=value, fill=Site)) + geom_boxplot()  + theme_classic() + 
  annotate("text", x=1.5, y=1, label="q-val = 0.02") + 
  scale_fill_manual(values=c("#8EB7CE","#EE6363" )) + labs(title = "Cellobiose transport system", x="", y = "clr-trasformed counts")

p5=ggplot(subset(clr_l, variable == "M00222"), aes(x=Site, y=value, fill=Site)) + geom_boxplot()  + theme_classic() + 
  annotate("text", x=1.5, y=5.2, label="q-val = 1.03e-02") + 
  scale_fill_manual(values=c("#8EB7CE","#EE6363" )) + labs(title = "Phosphate transport system", x="", y = "clr-trasformed counts")

p6=ggplot(subset(clr_l, variable == "M00212"), aes(x=Site, y=value, fill=Site)) + geom_boxplot()  + theme_classic() + 
  annotate("text", x=1.5, y=4.5, label="q-val = 2.32e-11") + 
  scale_fill_manual(values=c("#8EB7CE","#EE6363" )) + labs(title = "Ribose transport system", x="", y = "clr-trasformed counts")

pdf("outputs/light_pollution_picrust.pdf", width = 12, height = 7, pointsize = 12)
(p1+p2+p3)/(p4+p5+p6)
dev.off() 

