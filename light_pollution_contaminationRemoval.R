###################################################################
#####Identifying and removing contaminant ASVs normalized data#####
###################################################################

setwd("~/Documents/Bioinformatics_scripts/R_scripts/light_pollution_R/")

#read in asv table generated in mothur
asv = read.table("Inputs/light_pollution_ASV_table.txt", header = T)
met = read.table("Inputs/Final_metadata.txt", header = T, sep = "\t", row.names = 1)
names(asv)
tax=asv[,c(98:105)]
#Remove samples with low read number
asv.n=apply(asv[,c(1:96)], 2, as.numeric)
asv.o=asv.n[, colSums(asv.n) > 1000]
dim(asv.o)
rownames(asv.o)=rownames(asv)
message(ncol(asv.o)," samples with > 1000 reads were retained out of ", ncol(asv.n), " total samples")

#Identify and removing contaminant ASVs raw data
asv.r=as.data.frame(sweep(asv.o,2,colSums(asv.o),`/`))
asv.r=transform(asv.r,  Sum = rowSums(asv.r[,2:ncol(asv.r)]))
names(asv.r)
asv.r=transform(asv.r,  SumNegs = asv.r[,94]) # define negative controls here
names(asv.r)
asv.r=transform(asv.r,  contaFactor=(asv.r$SumNegs/asv.r$Sum)*100)
rownames(asv.r)=rownames(asv)
Conta=subset(asv.r, asv.r$contaFactor > 10)
Conta$Family=asv$Family[match(rownames(Conta), rownames(asv))]
message("Number of total ASVs: ", nrow(asv))
message("Number of identified contaminant ASVs removed from the analysis: ", length(rownames(Conta)), "\n", Conta$Family[1],"\n", Conta$Family[2],"\n", Conta$Family[3],"\n", Conta$Family[4],"\n", Conta$Family[5])

# Export normalized and raw ASV tables
asv.noConta=subset(asv.o, !rownames(asv.o) %in% rownames(Conta))[,-94]
colnames(asv.noConta)
dim(asv.noConta)
asv.noConta.f=merge(asv.noConta, tax, by="row.names")
write.table(asv.noConta.f, "./outputs/light_pollution_ASVs_noContanoOut.raw.txt",  quote = FALSE, row.names=F, sep = "\t") #define sample range
message("Number of ASVs used in the analysis: ", length(rownames(asv.noConta.f)))

