if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChAMP")
#install via terminal
git clone https://github.com/YuanTian1991/ChAMP.git
sudo R CMD INSTALL ChAMP
library(ChAMP)

setwd("/media/ankitv/Archivio2/ankit/Array21/ICF_hg38")
dir()
#myLoad <- champ.filter(beta=myImport$beta, arraytype = "EPIC", filterDetP=TRUE, detPcut=0.01)
#data(EPIC.manifest.hg38)
#LOAD FILE (method=ChAMP) PER SWAN NORMALIZATION###########
myLoad_2<-champ.load(directory = getwd(),
                     method="ChAMP",
                     methValue="B",
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.01,
                     filterBeads=TRUE,
                     beadCutoff=0.05,
                     filterNoCG=TRUE,
                     filterSNPs=TRUE,
                     population=NULL,
                     filterMultiHit=TRUE,
                     filterXY=TRUE,
                     force=FALSE,
                     arraytype="EPIC")
#write.table(myLoad_2$beta,"myLoad_2_beta.txt",row.names=TRUE,quote=FALSE)
summary(myLoad_2)
head(myLoad_2$beta)
head(myLoad_2$pd)
dim(myLoad_2$beta)
myLoad_3 <- data.frame(myLoad_2$beta[order(rownames(myLoad_2$beta)),])
head(myLoad_3,2)
dim(myLoad_3)
#QC
QC.GUI(myLoad_2$beta,
       pheno=myLoad_2$pd$Sample_Group, 
       arraytype="EPIC")
champ.QC(beta = myLoad_2$beta,
         pheno=myLoad_2$pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_myLoad/")
head(myLoad_2$beta)
DATA_myLoad2 <- cbind.data.frame(rownames(myLoad_2$beta), myLoad_2$beta)
head(DATA_myLoad2)
colnames(DATA_myLoad2) <- c("TargetID",colnames(myLoad_2$beta))
head(DATA_myLoad2)
dim(DATA_myLoad2)
DATA_myLoad2sort <- DATA_myLoad2[order(DATA_myLoad2$TargetID),]
head(DATA_myLoad2sort)

# NORMALIZE DATASET (ChAMP load)
myNorm2<- champ.norm(myLoad_2$beta,
                    method="BMIQ",
                    arraytype="EPIC")
dim(myNorm2)
head(myNorm2)
write.table(myNorm2,"myNorm2.txt",row.names=TRUE,quote=FALSE,append = F,sep = "\t")

#myNorm2 <- read.table("myNorm2.txt", header = TRUE)
myNorm3 <- data.frame(myNorm2[order(rownames(myNorm2)),])
head(myNorm3)
#write.table(myNorm3,"myNorm3.txt",row.names=TRUE,quote=FALSE,append = F,sep = "\t")

#other features not used
#QUALITY CONTROL-TIME CONSUMING BE-CARE
QC.GUI(myNorm2,
       pheno=myLoad_2$pd$Sample_Group,
       arraytype="EPIC")
champ.QC(beta = myNorm2,
         pheno=myLoad_2$pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=TRUE,
         dendrogram=TRUE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages_myNorm2/")
# Do SVD check on data set
champ.SVD(myNorm2,
          rgSet=NULL,
          pd=myLoad_2$pd)

#save as Batch_effects_SVD.svg
#######################################################
##################### No Batch correction #############
#######################################################

#write.table(myCombat2,"myCombat_2.txt",row.names=TRUE,quote=FALSE)
# If Batch detected, run champ.runCombat() here.
#MERGE myNorm2 to annotations of Illumina Manifest
#before Add TargetID to first header of myNorm2.txt
#reload modified file
#DATA_myNorm_2<-read.delim("myNorm2.txt", sep=" ")
DATA_myNorm_2 <- cbind.data.frame(rownames(myNorm2), myNorm2)
head(DATA_myNorm_2)
colnames(DATA_myNorm_2) <- c("TargetID",colnames(myNorm2))
head(DATA_myNorm_2)
dim(DATA_myNorm_2)
#LOAD MANIFEST hg19 B4 https://bioconductor.org/packages/release/data/annotation/manuals/IlluminaHumanMethylationEPICanno.ilm10b4.hg19/man/IlluminaHumanMethylationEPICanno.ilm10b4.hg19.pdf. Rearranged to csv
MANIFEST <- read.csv("/media/ankitv/Archivio1/testing_array/Loca/test/MethylationEPIC_v-1-0_B4_ext.csv")

#https://bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
library(liftOver)
#Setup data
head(MANIFEST)
MANIFEST_cords <- paste0(MANIFEST$CHR, "%", MANIFEST$MAPINFO,"%", MANIFEST$MAPINFO,"%", MANIFEST$Name)
cur19 <- MANIFEST[,c(12,13,13,15,1)]
head(cur19)
cur19["chr"] <- paste0("chr", cur19$CHR)
cur19 <- cur19[,c(6,2:5)]
write.table(cur19,"cur19.txt",row.names=FALSE,quote=FALSE, col.names = F, sep="\t", append = F)
grep cg cur19.txt | sort -k1,1 -k2,2n >  cur19_cg.txt
cur19_cg <- read.table("cur19_cg.txt", header = F, stringsAsFactors = F)
head(cur19_cg)
colnames(cur19_cg) <- c("chr", "start", "end", "strand", "ProbeID")
cur19_cg <- cur19_cg[,c(1,2,3,5)]
curgR <- makeGRangesFromDataFrame(cur19_cg, keep.extra.columns=TRUE)
#Resource: The chain file for hg19 to hg38 transformation
library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch = import.chain(path)
ch

str(ch[[1]])

#Action: liftOver
seqlevelsStyle(curgR) = "UCSC"  # necessary
cur38 = liftOver(curgR, ch)
class(cur38)

cur38 = unlist(cur38)
genome(cur38) = "hg38"
cur38

length(curgR)-length(cur38)

curhg38 <- data.frame(seqnames=seqnames(cur38),
                 starts=start(cur38),
                 ends=end(cur38),
                 names=elementMetadata(cur38)$ProbeID)
head(curhg38)
write.table(curhg38, file="curhg38.bed", append = F ,quote=F, sep="\t", row.names=F, col.names=F)
colnames(curhg38) <- c("chrhg38", "starthg38", "endhg38", "IlmnID")
colnames(MANIFEST)
dim(MANIFEST)
head(MANIFEST)
curhg38Manifest <- merge(curhg38,MANIFEST,by.x="IlmnID",by.y="IlmnID")
dim(curhg38Manifest)
colnames(curhg38Manifest)
head(curhg38Manifest)
write.table(curhg38Manifest, "curhg38Manifest.txt", quote=F, sep="\t", append = F, row.names = F)
curhg38Manifest_cords <- paste0(curhg38Manifest$chrhg38, "%", curhg38Manifest$starthg38, "%", curhg38Manifest$endhg38, "%", curhg38Manifest$Name)
head(curhg38Manifest_cords)
#To gain confidence over conversion of hg19 to hg38, I downloaded a precompiled file from https://zwdzwd.github.io/InfiniumAnnotation: Basic%20hg38%20annotation%20with%20suggested%20overall%20masking%20(EPIC
#Gunzip and read tsv file
EPIC.hg38.manifest <- read.table(file='/media/ankitv/Archivio1/testing_array/Loca/test/annoepic/EPIC.hg38.manifest.tsv', sep = '\t', header = TRUE, fill = TRUE)
head(EPIC.hg38.manifest)
#This file has 1-based coordinate format so I just add +1 to start  and consider it as both start and end
EPIC.hg38.manifest["Chr"] <- EPIC.hg38.manifest$CpG_chrm
EPIC.hg38.manifest["Start"] <- EPIC.hg38.manifest$CpG_beg +1
EPIC.hg38.manifest["End"] <- EPIC.hg38.manifest$CpG_beg +1

EPIC.hg38.manifest_cords <- paste0(EPIC.hg38.manifest$Chr,"%", EPIC.hg38.manifest$Start, "%", EPIC.hg38.manifest$End, "%", EPIC.hg38.manifest$probeID)

#Check overlap
#hg19 and hg38 converted 1_vennn
library(gplots) 
venn(list(MANIFEST_cords, curhg38Manifest_cords))
#hg38 epic and hg38 converted 2_venn
venn(list(EPIC.hg38.manifest_cords, curhg38Manifest_cords))
#hg38 epic and hg38 converted 3_venn
venn(list(EPIC.hg38.manifest_cords, MANIFEST_cords))

#Since 859497 probes are common between epic_hg38 and liftover curhg38, I used curhg38 for further analysi.
myNorm2 <- read.table("myNorm2.txt", header = T)
head(myNorm2)
tmyNorm2 = t(myNorm2)
#head(tmyNorm2)
dim(tmyNorm2)
#rownames(tmyNorm2)
tmyNorm2 = data.frame(tmyNorm2)
#head(tmyNorm2)
#write.table(tmyNorm2 , "tmyNorm2dedup.txt", sep="\t", quote = FALSE, append = FALSE)
#tmyNorm2["Color"] <-  c("cNormalEpic16","bBWSEpic16","bBWSEpic16","cNormalEpic16","bBWSEpic16","cNormalEpic16","cNormalEpic16","bBWSEpic16","cNormalEpic16","cNormalEpic16","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aaControlLoca","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16","aControlEpic16")
tmyNorm2["Color"] <-  c("corrClone","corrClone","corrClone","corrClone","Control","Case","Case","Control")
dim(tmyNorm2)
dfx <-tmyNorm2[c(1:729682)]
#head(dfx)
PC<-prcomp(dfx, center = TRUE, scale. = TRUE)
#head(PC)
PCi<-data.frame(PC$x,Color=tmyNorm2$Color)
percentageAll <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentageAll <- paste( colnames(PCi), "(", paste( as.character(percentageAll), "%", ")", sep="") )
library(ggplot2)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p <- ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentageAll[1]) + ylab(percentageAll[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("#FFB6B3","black","#BDE7BD","grey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  scale_shape_manual(values=c(0,1,2))
p <- p+theme_classic()
#p + xlim(-50,50)+ ylim(-50,50)
p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA-tmyNorm2.svg", width=11, height=8.6, units="cm", dpi=96)
svg(filename="boxplot_myNorm2.svg", width=10, height=5, pointsize=12)
boxplot(myNorm2, main="Normalized beta values", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("blue","blue","blue","blue","blue","blue","blue","blue","#FFB6B3","#FFB6B3","#FFB6B3","#FFB6B3","#FFB6B3","#FFB6B3","blue","blue","blue","blue","blue","blue"))
dev.off()


library(ggpubr)
library(ggplot2)
library(plyr)

head(myNorm3)



myNorm3re <- as.matrix(myNorm3)
head(myNorm3re)
head(myNorm3re,2)
dim(myNorm3re)
myiNorm3reprobed <- cbind.data.frame(rownames(myNorm3re),
                                     myNorm3re[,1:8])

head(myiNorm3reprobed)
colnames(myiNorm3reprobed) <- c("TargetID",colnames(myNorm2))
MERGE_myiNorm3re <- merge(myiNorm3reprobed,curhg38Manifest,by.x="TargetID",by.y="IlmnID")
head(MERGE_myiNorm3re)
MERGE_myiNorm3re_pos <- MERGE_myiNorm3re[,c(10,11,12,1:9)]
#check for dimensions
dim(MERGE_myiNorm3re_pos) #All cg annotated with coordinates
head(MERGE_myiNorm3re_pos)
tail(MERGE_myiNorm3re_pos)
write.table(MERGE_myiNorm3re_pos, "MERGE_myiNorm3re_pos.txt", col.names = F, quote = F, row.names = F)

head(myNorm2)
myNorm2sort <- myNorm2[order(rownames(myNorm2)),]
#Density distribution curve
myNorm2 <- data.frame(myNorm2)
head(myNorm2)
myNorm2stack <- stack(myNorm2)
head(myNorm2stack)
ggplot(myNorm2stack) + geom_density(aes(x = values, color = ind)) + scale_colour_manual(values = c("#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","darkgrey","#FFB6B3","#FFB6B3","darkgrey"))
# Density_myNorm2stack.svg







awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_pos.txt | sort -k1,1 -k2,2n > MERGE_myiNorm3re_pos_chr.txt
bedtools intersect -wa -wb -a MERGE_myiNorm3re_pos_chr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed > MERGE_myiNorm3re_human_ICR.txt
awk '{print $16"\t"$17"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_human_ICR.txt > MERGE_myiNorm3re_human_ICR.rearranged.txt

library(ggpubr)
library(ggplot2)
library(plyr)

#PCA with  Human ICR
inormdata1 <- read.table("MERGE_myiNorm3re_human_ICR.rearranged.txt", header = F)
colnames(inormdata1) <- c("DMR", "DMRType","chr","start","end","TargetID", colnames(myNorm2))
rownames(inormdata1)
colnames(inormdata1)
head(inormdata1)
dim(inormdata1)
Count_inormdataICR1 <- count(inormdata1, "DMR")
head(Count_inormdataICR1)
dim(Count_inormdataICR1)
scomaggregate1 = aggregate(inormdata1[,7:14],by=list(inormdata1$DMR), mean)
head(scomaggregate1, 2)
#Plot PCA by group color and labelling
inormdataICR1=scomaggregate1
rownames(inormdataICR1)
inormdataICR1[,1]
rownames(inormdataICR1)=inormdataICR1[,1]
rownames(inormdataICR1)
colnames(inormdataICR1)
inormdataICR1 = inormdataICR1[,-1]
head(inormdataICR1)
dim(inormdataICR1)
inormdataICR1_counted <- cbind(inormdataICR1, Count_inormdataICR1)
head(inormdataICR1_counted)
dim(inormdataICR1_counted)
inormdataICR1_counted1 <- inormdataICR1_counted[which(inormdataICR1_counted$freq >= 3),]
head(inormdataICR1_counted1)
dim(inormdataICR1_counted1)
write.table(inormdataICR1_counted1, "inormdataICR1_counted1.filtmin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)


inormdataICR1_counted2 <- inormdataICR1_counted1[,1:8]
head(inormdataICR1_counted2)
summary(inormdataICR1_counted2)
#df <- as.inormdata.frame(inormdataICR1)
inormdfICR <- inormdataICR1_counted2
head(inormdfICR)
dim(inormdfICR)
inormdfICR = data.frame(inormdfICR)
#write.table(inormdfICR , "inormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
tinormdfICR = t(inormdfICR)
tinormdfICR = data.frame(tinormdfICR)
#write.table(tinormdfICR , "tinormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfICR)
tinormdfICR["Color"] <-  c("corrClone","corrClone","corrClone","corrClone","Control","Case","Case","Control")
#head(tinormdfICR)
#write.table(tinormdfICR , "tinormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfICR)
inormdfx <-tinormdfICR[c(1:43)]
PcomC<-prcomp(inormdfx, center = TRUE, scale. = TRUE)
PcomCi<-data.frame(PcomC$x,Color=tinormdfICR$Color)
percentagecomICR <- round(PcomC$sdev^2 / sum(PcomC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomICR <- paste( colnames(PcomCi), "(", paste( as.character(percentagecomICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pcom1 <-ggplot(PcomCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomICR[1]) + ylab(percentagecomICR[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("#FFB6B3","darkgrey","#BDE7BD"))+
  scale_shape_manual(values=c(0,1,2))
pcom1 <- pcom1 +theme_classic()
pcom1 + xlim(-20,20)+ ylim(-20,20)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_tinormdfICR_CntrlIndiv.svg", width=10*1.25, height=8*1.25, units="cm", dpi=96)


################################### Heatmap_indiv human ICR ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(inormdataICR1_counted1)
dim(inormdataICR1_counted1)
sum(inormdataICR1_counted1$freq) #714-705=9 CpGs  removed while min 3 filtering

#Heatmap_indiv with 12 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
inormdataICR2 = as.matrix(inormdataICR1_counted2)
head(inormdataICR2)
dim(inormdataICR2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
inormdataICR3 <- inormdataICR2
head(inormdataICR3)
write.table(inormdataICR3, "inormdataICR3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
#Control Samples
#Heatmap_indiv.2(inormdataICR2, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
#Let it cluster
svg(filename="Heatmap_com_indiv_pimethhed.avg_human_ICR.rearranged_SVG.svg", width=5, height=10, pointsize=12)
#Heatmap_indiv.2(inormdataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(inormdataICR3, col = colfunc, Colv = "NA",dendrogram = c("none"), trace = "none", keysize=0.8, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

svg(filename="Heatmap_com_indiv_pimethhed.avg_human_ICR.rearranged_SVG_clusterboth.svg", width=15, height=15, pointsize=12)
#Heatmap_indiv.2(inormdataICR3, trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
heatmap.2(inormdataICR3, col = colfunc, dendrogram = c("both"), trace = "none", key=TRUE, density.info=c("none"), scale = "none", sepwidth=c(0.001, 0.001), cexRow=0.8, font=3, cexCol = 1, margins =c(5,8), breaks = seq(0,1, length.out = 100))
dev.off()

#Box and Violin Plot human ICR 
rownames(inormdata1)
colnames(inormdata1)
head(inormdata1)
dim(inormdata1)
Count_inormdataICR1 <- count(inormdata1, "DMR")
head(Count_inormdataICR1)
inormagregateicr1 = aggregate(inormdata1[,7:14],by=list(inormdata1$DMR), mean)
head(inormagregateicr1, 2)
inormdataICR1=inormagregateicr1
rownames(inormdataICR1)
inormdataICR1[,1]
rownames(inormdataICR1)=inormdataICR1[,1]
rownames(inormdataICR1)
colnames(inormdataICR1)
inormdataICR1 = inormdataICR1[,-1]
head(inormdataICR1)
dim(inormdataICR1)
write.table(inormdataICR1, "inormaggregatedicr_inormdataICR1.txt", sep="\t", quote = FALSE, append = FALSE)
inormdataICR1_counted <- cbind(inormdataICR1, Count_inormdataICR1)
head(inormdataICR1_counted)
dim(inormdataICR1_counted)
inormdataICR1_counted1 <- inormdataICR1_counted[which(inormdataICR1_counted$freq >= 3),]
head(inormdataICR1_counted1)
dim(inormdataICR1_counted1)
write.table(inormdataICR1_counted1, "inormdataICR1_counted1.filtmin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
inormdataICR1_counted2 <- inormdataICR1_counted1[,1:8]
head(inormdataICR1_counted2)
VionormIndvCR <- data.frame(inormdataICR1_counted2[,c(5,8,6,7,1:4)])
head(VionormIndvCR)
write.table(VionormIndvCR, "VionormIndvCR.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_VionormIndvICR.svg", width=10, height=5, pointsize=12)
boxplot(VionormIndvCR, main="Average methylation at human ICRs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","navy","navy","navy","navy"))
dev.off()
dim(VionormIndvCR)
VionormIndvCR <- data.frame(VionormIndvCR)
VionormIndvCR1 <- VionormIndvCR
head(VionormIndvCR1,1)
VionormIndvCR2 <- stack(VionormIndvCR1)
head(VionormIndvCR2)
colnames(VionormIndvCR2) <- c("Methylation", "inormdatasets")
head(VionormIndvCR2)
ggplot(VionormIndvCR2, aes(x=inormdatasets, y=Methylation, color=inormdatasets, fill=inormdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_bw() +scale_color_manual(values=c("black","black","black","black","black","black","black","black"))+scale_fill_manual(values=c("grey","grey","#FFB6B3","#FFB6B3","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD"))
ggsave("Violin_plot_VionormIndvCR2_pimethhed_humanICR.svg", width=15*1.25, height=8*1.25, units="cm", dpi=96)
ggsave("Violin_plot_VionormIndvCR2_pimethhed_humanICR.png", width=15*1.25, height=8*1.25, units="cm", dpi=96)


#Box and Violin Plot Global methylation 
head(myNorm3re)
dim(myNorm3re)
inormcomdataGlob <- myNorm3re
rownames(inormcomdataGlob)
inormcomdataGlob[,1]
inormcomdataGlob1 <- data.frame(inormcomdataGlob)
head(inormcomdataGlob1)
str(inormcomdataGlob1)
VionormGlobIndv <- inormcomdataGlob1
head(VionormGlobIndv)

svg(filename="boxplot_myNorm3re.svg", width=10, height=5, pointsize=12)
boxplot(VionormGlobIndv, main="Global average methylation", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("navy","navy","navy","navy","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))
dev.off()
head(VionormGlobIndv)
VionormGlobIndv1 <- VionormGlobIndv[,c(5,8,6,7,1:4)]
head(VionormGlobIndv1)
dim(VionormGlobIndv1)
VionormGlobIndv2 <- stack(VionormGlobIndv1)
colnames(VionormGlobIndv2) <- c("Methylation", "inormdatasets")
head(VionormGlobIndv2)
ggplot(VionormGlobIndv2, aes(x=inormdatasets, y=Methylation, color=inormdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","#FFB6B3","#FFB6B3","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","navy","navy","navy","navy"))
ggsave("Violin_plot_VionormGlobIndv2_pimethhed_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)


#Heatmap of control normalized Log other samples merged
head(inormdataICR1_counted2)
dim(inormdataICR1_counted2)
#rearrange columns
inormdataICR1_counted2ratio <- as.matrix(inormdataICR1_counted2[,c(5,8,6,7,1:4)])

head(inormdataICR1_counted2ratio)
dim(inormdataICR1_counted2ratio)
winormdataICR1_counted2ratioavg <- data.frame(cbind((rowMeans(inormdataICR1_counted2ratio[,1:2])),
                                                    (inormdataICR1_counted2ratio[,1:8])))

head(winormdataICR1_counted2ratioavg)
colnames(winormdataICR1_counted2ratioavg) <- c("Allcontrol",colnames(inormdataICR1_counted2ratio))
head(winormdataICR1_counted2ratioavg)
winormdataICR1_counted2ratioavg1 <- winormdataICR1_counted2ratioavg
winormdataICR1_counted2ratioavg2 <- as.matrix(winormdataICR1_counted2ratioavg1)
head(winormdataICR1_counted2ratioavg2)
dim(winormdataICR1_counted2ratioavg2)
winormdataICR1_counted2ratioavg3 <- data.frame(cbind((winormdataICR1_counted2ratioavg2[,1]/winormdataICR1_counted2ratioavg2[,1]),
                                                     (winormdataICR1_counted2ratioavg2[,2]/winormdataICR1_counted2ratioavg2[,1]),
                                                     (winormdataICR1_counted2ratioavg2[,3]/winormdataICR1_counted2ratioavg2[,1]),
                                                     (winormdataICR1_counted2ratioavg2[,4]/winormdataICR1_counted2ratioavg2[,1]),
                                                     (winormdataICR1_counted2ratioavg2[,5]/winormdataICR1_counted2ratioavg2[,1]),
                                                     (winormdataICR1_counted2ratioavg2[,6]/winormdataICR1_counted2ratioavg2[,1]),
                                                     (winormdataICR1_counted2ratioavg2[,7]/winormdataICR1_counted2ratioavg2[,1]),
                                                     (winormdataICR1_counted2ratioavg2[,8]/winormdataICR1_counted2ratioavg2[,1]),
                                                     (winormdataICR1_counted2ratioavg2[,9]/winormdataICR1_counted2ratioavg2[,1])))


head(winormdataICR1_counted2ratioavg3)
colnames(winormdataICR1_counted2ratioavg3) <- colnames(winormdataICR1_counted2ratioavg2)
head(winormdataICR1_counted2ratioavg3)
winormdataICR1_counted2ratioavg3 <- as.matrix(winormdataICR1_counted2ratioavg3)
winormdataICR1_counted2ratioavg3Log <- log2(winormdataICR1_counted2ratioavg3)
head(winormdataICR1_counted2ratioavg3Log)
write.table(winormdataICR1_counted2ratioavg3Log , "Heatmap_Loca_Padi6_normcontrolwinormdataICR1_counted2ratioavg2Log.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)

winormdataICR1_counted2ratioavg3Logminuscontrol <- winormdataICR1_counted2ratioavg3Log[,c(2:9)]
head(winormdataICR1_counted2ratioavg3Logminuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(winormdataICR1_counted2ratioavg3Logminuscontrol))
row.names(my_sample_col2) <- colnames(winormdataICR1_counted2ratioavg3Logminuscontrol)
my_colour2 = list(Annotations = c("D250"= "darkgrey","UN"= "darkgrey","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(winormdataICR1_counted2ratioavg3Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_winormdataICR1_counted2ratioavg.normcontrol.svg

minormdataICR1_counted2ratioavg3 <- data.frame(cbind((winormdataICR1_counted2ratioavg2[,1:9]-winormdataICR1_counted2ratioavg2[,1])))


head(minormdataICR1_counted2ratioavg3)
minormdataICR1_counted2ratioavg3 <- as.matrix(minormdataICR1_counted2ratioavg3)

minormdataICR1_counted2ratioavg3minuscontrol <- minormdataICR1_counted2ratioavg3[,c(2:9)]
head(minormdataICR1_counted2ratioavg3minuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(minormdataICR1_counted2ratioavg3minuscontrol))
row.names(my_sample_col2) <- colnames(minormdataICR1_counted2ratioavg3minuscontrol)
my_colour2 = list(Annotations = c("D250"= "darkgrey","UN"= "darkgrey","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(minormdataICR1_counted2ratioavg3minuscontrol,
         color = colorRampPalette(c("navy", "white", "#D47400"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "black",
         cutree_cols = 3)

#save as pHeatmap_minormdataICR1_counted2ratioavg.normcontrol.svg

minormdataICR1_counted2ratioavg3Medcontrol <- minormdataICR1_counted2ratioavg3[,-c(2:3)]
head(minormdataICR1_counted2ratioavg3Medcontrol)
library(pheatmap)
my_sample_col6 <- data.frame(Annotations= colnames(minormdataICR1_counted2ratioavg3Medcontrol))
row.names(my_sample_col6) <- colnames(minormdataICR1_counted2ratioavg3Medcontrol)
my_colour6 = list(Annotations = c("Allcontrol"= "darkgrey","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

pheatmap(minormdataICR1_counted2ratioavg3Medcontrol,
         color = colorRampPalette(c("navy", "white", "#D47400"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour6,
         fontsize = 8,
         annotation_col = my_sample_col6,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 1)

#save as pHeatmap_minormdataICR1_counted2ratioavg3Medcontrol.svg
MERGE_myiNorm3re_pos_chr <- read.table("MERGE_myiNorm3re_pos_chr.txt", header = F, stringsAsFactors = F)
head(MERGE_myiNorm3re_pos_chr)
colnames(MERGE_myiNorm3re_pos_chr) <- c("chr","start","end","probeID",colnames(myNorm2))
rownames(MERGE_myiNorm3re_pos_chr) <- paste0(MERGE_myiNorm3re_pos_chr$chr, 
                                             "%",
                                             MERGE_myiNorm3re_pos_chr$start,
                                             "%",
                                             MERGE_myiNorm3re_pos_chr$end,
                                             "%",
                                             MERGE_myiNorm3re_pos_chr$probeID)

MERGE_myiNorm3re_pos_chrPre <- MERGE_myiNorm3re_pos_chr[,-(1:4)]
head(MERGE_myiNorm3re_pos_chrPre)
MERGE_myiNorm3re_pos_chrPre <- MERGE_myiNorm3re_pos_chrPre[,c(5,8,6,7,1:4)]
head(MERGE_myiNorm3re_pos_chrPre)
dim(MERGE_myiNorm3re_pos_chrPre)
MERGE_myiNorm3re_pos_chrPre <- as.matrix(MERGE_myiNorm3re_pos_chrPre)
MERGE_myiNorm3re_pos_chravg <- data.frame(cbind((rowMeans(MERGE_myiNorm3re_pos_chrPre[,1:2])),
                                        (MERGE_myiNorm3re_pos_chrPre[,1:8])))

head(MERGE_myiNorm3re_pos_chravg)
colnames(MERGE_myiNorm3re_pos_chravg) <- c("Allcontrol",colnames(MERGE_myiNorm3re_pos_chrPre))
dim(MERGE_myiNorm3re_pos_chravg)
#Subtraction, Control -Control
MERGE_myiNorm3re_pos_chravg["D250_UN"] <- MERGE_myiNorm3re_pos_chravg$D250 - MERGE_myiNorm3re_pos_chravg$UN

#Subtraction, Case -Control
MERGE_myiNorm3re_pos_chravg["PG_Con"] <- MERGE_myiNorm3re_pos_chravg$PG - MERGE_myiNorm3re_pos_chravg$Allcontrol
MERGE_myiNorm3re_pos_chravg["PR_Con"] <- MERGE_myiNorm3re_pos_chravg$PR - MERGE_myiNorm3re_pos_chravg$Allcontrol
MERGE_myiNorm3re_pos_chravg["C7_Con"] <- MERGE_myiNorm3re_pos_chravg$C7 - MERGE_myiNorm3re_pos_chravg$Allcontrol
MERGE_myiNorm3re_pos_chravg["C13_Con"] <- MERGE_myiNorm3re_pos_chravg$C13 - MERGE_myiNorm3re_pos_chravg$Allcontrol
MERGE_myiNorm3re_pos_chravg["C35_Con"] <- MERGE_myiNorm3re_pos_chravg$C35 - MERGE_myiNorm3re_pos_chravg$Allcontrol
MERGE_myiNorm3re_pos_chravg["C50_Con"] <- MERGE_myiNorm3re_pos_chravg$C50 - MERGE_myiNorm3re_pos_chravg$Allcontrol
MERGE_myiNorm3re_pos_chravg["PG_D250"] <- MERGE_myiNorm3re_pos_chravg$PG - MERGE_myiNorm3re_pos_chravg$D250
MERGE_myiNorm3re_pos_chravg["PR_D250"] <- MERGE_myiNorm3re_pos_chravg$PR - MERGE_myiNorm3re_pos_chravg$D250
MERGE_myiNorm3re_pos_chravg["C7_D250"] <- MERGE_myiNorm3re_pos_chravg$C7 - MERGE_myiNorm3re_pos_chravg$D250
MERGE_myiNorm3re_pos_chravg["C13_D250"] <- MERGE_myiNorm3re_pos_chravg$C13 - MERGE_myiNorm3re_pos_chravg$D250
MERGE_myiNorm3re_pos_chravg["C35_D250"] <- MERGE_myiNorm3re_pos_chravg$C35 - MERGE_myiNorm3re_pos_chravg$D250
MERGE_myiNorm3re_pos_chravg["C50_D250"] <- MERGE_myiNorm3re_pos_chravg$C50 - MERGE_myiNorm3re_pos_chravg$D250
MERGE_myiNorm3re_pos_chravg["PG_UN"] <- MERGE_myiNorm3re_pos_chravg$PG - MERGE_myiNorm3re_pos_chravg$UN
MERGE_myiNorm3re_pos_chravg["PR_UN"] <- MERGE_myiNorm3re_pos_chravg$PR - MERGE_myiNorm3re_pos_chravg$UN
MERGE_myiNorm3re_pos_chravg["C7_UN"] <- MERGE_myiNorm3re_pos_chravg$C7 - MERGE_myiNorm3re_pos_chravg$UN
MERGE_myiNorm3re_pos_chravg["C13_UN"] <- MERGE_myiNorm3re_pos_chravg$C13 - MERGE_myiNorm3re_pos_chravg$UN
MERGE_myiNorm3re_pos_chravg["C35_UN"] <- MERGE_myiNorm3re_pos_chravg$C35 - MERGE_myiNorm3re_pos_chravg$UN
MERGE_myiNorm3re_pos_chravg["C50_UN"] <- MERGE_myiNorm3re_pos_chravg$C50 - MERGE_myiNorm3re_pos_chravg$UN
MERGE_myiNorm3re_pos_chravg["C7_PR"] <- MERGE_myiNorm3re_pos_chravg$C7 - MERGE_myiNorm3re_pos_chravg$PR
MERGE_myiNorm3re_pos_chravg["C13_PG"] <- MERGE_myiNorm3re_pos_chravg$C13 - MERGE_myiNorm3re_pos_chravg$PG
MERGE_myiNorm3re_pos_chravg["C35_PR"] <- MERGE_myiNorm3re_pos_chravg$C35 - MERGE_myiNorm3re_pos_chravg$PR
MERGE_myiNorm3re_pos_chravg["C50_PG"] <- MERGE_myiNorm3re_pos_chravg$C50 - MERGE_myiNorm3re_pos_chravg$PG

head(MERGE_myiNorm3re_pos_chravg)
write.table(MERGE_myiNorm3re_pos_chravg, "MERGE_myiNorm3re_pos_chravg.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)

library(splitstackshape)
MERGE_myiNorm3re_pos_chravg1 <- MERGE_myiNorm3re_pos_chravg
MERGE_myiNorm3re_pos_chravg1["row"] <- rownames(MERGE_myiNorm3re_pos_chravg1)
MERGE_myiNorm3re_pos_chravg1 <- cSplit(MERGE_myiNorm3re_pos_chravg1, "row", "%")
head(MERGE_myiNorm3re_pos_chravg1)
dim(MERGE_myiNorm3re_pos_chravg1)
MERGE_myiNorm3re_pos_chravg1 <- MERGE_myiNorm3re_pos_chravg1[,c(33:36,1:32)]
head(MERGE_myiNorm3re_pos_chravg1)
write.table(MERGE_myiNorm3re_pos_chravg1, "MERGE_myiNorm3re_pos_chravg1.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)




#Regions where two controls are close, first select sample specific diff meth probes thanprobes where two controls are close
pMERGE_myiNorm3re_pos_chravg2 <- data.frame(MERGE_myiNorm3re_pos_chravg)
dim(pMERGE_myiNorm3re_pos_chravg2)
head(pMERGE_myiNorm3re_pos_chravg2)
#PG diff
pMERGE_myiNorm3re_pos_chravg_PG <- pMERGE_myiNorm3re_pos_chravg2[which(pMERGE_myiNorm3re_pos_chravg2$PG_D250 > 0.2 | 
                                                                          pMERGE_myiNorm3re_pos_chravg2$PG_D250 < -0.2),]

miMERGE_myiNorm3re_pos_chravg_PG <- pMERGE_myiNorm3re_pos_chravg_PG[which(pMERGE_myiNorm3re_pos_chravg_PG$PG_UN > 0.2 | 
                                                                            pMERGE_myiNorm3re_pos_chravg_PG$PG_UN < -0.2),]

miMERGE_myiNorm3re_pos_chravg_PG <- miMERGE_myiNorm3re_pos_chravg_PG[which(miMERGE_myiNorm3re_pos_chravg_PG$D250_UN < 0.1 & 
                                                                             miMERGE_myiNorm3re_pos_chravg_PG$D250_UN > -0.1),]

head(miMERGE_myiNorm3re_pos_chravg_PG)
dim(miMERGE_myiNorm3re_pos_chravg_PG)
miMERGE_myiNorm3re_pos_chravg_PG["row"] <- rownames(miMERGE_myiNorm3re_pos_chravg_PG)
library(splitstackshape)
miMERGE_myiNorm3re_chrpos_avg_PG <- cSplit(miMERGE_myiNorm3re_pos_chravg_PG, "row", "%")
head(miMERGE_myiNorm3re_chrpos_avg_PG)
dim(miMERGE_myiNorm3re_chrpos_avg_PG)
miMERGE_myiNorm3re_chrpos_avg_PG <- miMERGE_myiNorm3re_chrpos_avg_PG[,c(33:36,1:32)]
write.table(miMERGE_myiNorm3re_chrpos_avg_PG, "miMERGE_myiNorm3re_chrpos_avg_PG.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myiNorm3re_chrpos_avg_PG$row_4, "miMERGE_myiNorm3re_chrpos_avg_PG.id", sep = "\t", quote = F, append = F, row.names = T, col.names = T)

#PR diff
pMERGE_myiNorm3re_pos_chravg_PR <- pMERGE_myiNorm3re_pos_chravg2[which(pMERGE_myiNorm3re_pos_chravg2$PR_D250 > 0.2 | 
                                                                         pMERGE_myiNorm3re_pos_chravg2$PR_D250 < -0.2),]

miMERGE_myiNorm3re_pos_chravg_PR <- pMERGE_myiNorm3re_pos_chravg_PR[which(pMERGE_myiNorm3re_pos_chravg_PR$PR_UN > 0.2 | 
                                                                            pMERGE_myiNorm3re_pos_chravg_PR$PR_UN < -0.2),]

miMERGE_myiNorm3re_pos_chravg_PR <- miMERGE_myiNorm3re_pos_chravg_PR[which(miMERGE_myiNorm3re_pos_chravg_PR$D250_UN < 0.1 & 
                                                                             miMERGE_myiNorm3re_pos_chravg_PR$D250_UN > -0.1),]

head(miMERGE_myiNorm3re_pos_chravg_PR)
dim(miMERGE_myiNorm3re_pos_chravg_PR)
miMERGE_myiNorm3re_pos_chravg_PR["row"] <- rownames(miMERGE_myiNorm3re_pos_chravg_PR)
library(splitstackshape)
miMERGE_myiNorm3re_chrpos_avg_PR <- cSplit(miMERGE_myiNorm3re_pos_chravg_PR, "row", "%")
head(miMERGE_myiNorm3re_chrpos_avg_PR)
dim(miMERGE_myiNorm3re_chrpos_avg_PR)
miMERGE_myiNorm3re_chrpos_avg_PR <- miMERGE_myiNorm3re_chrpos_avg_PR[,c(33:36,1:32)]
write.table(miMERGE_myiNorm3re_chrpos_avg_PR, "miMERGE_myiNorm3re_chrpos_avg_PR.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myiNorm3re_chrpos_avg_PR$row_4, "miMERGE_myiNorm3re_chrpos_avg_PR.id", sep = "\t", quote = F, append = F, row.names = T, col.names = T)

# Extract Imprinted CpGs
head(inormdata1)
dim(inormdata1)

#All Imprinted region CpGs
write.table(data.frame(inormdata1$TargetID), "inormdata1_TargetID.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
#fgrep -f inormdata1_TargetID.txt MERGE_myiNorm3re_pos_chravg1.txt -w > MERGE_myiNorm3re_pos_chravg1_imp.txt

#PG specific imprinted region CpGs
Manually make diff_imp_loci_pG.txt list by visualizing heatmap
diff_imp_loci_pG <- read.table("diff_imp_loci_pG.txt", header = F, stringsAsFactors = F)
diff_imp_loci_pG_CpGs <-  merge(inormdata1, diff_imp_loci_pG, by.x = "DMR", by.y="V1")
dim(diff_imp_loci_pG_CpGs)
pG_specific_impCpGs <- data.frame(diff_imp_loci_pG_CpGs$TargetID)
colnames(pG_specific_impCpGs) <- c("TargetID")
dim(pG_specific_impCpGs)
length(unique(pG_specific_impCpGs$TargetID))
dim(miMERGE_myiNorm3re_chrpos_avg_PG)
miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs <-  merge(miMERGE_myiNorm3re_chrpos_avg_PG, pG_specific_impCpGs, by.x = "row_4", by.y="TargetID")
head(miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs)
dim(miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs)
cat("Out of",dim(miMERGE_myiNorm3re_chrpos_avg_PG)[1]  ," DMPs in matrixPG and ",dim(pG_specific_impCpGs)[1]," in PG_specific_CpGs only ",dim(miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs)[1]," overlapped. Rest all removed while filtering at closeness of D250 and UN and difference of +/- 0.2 among two Controls v/s PG")


#PR specific imprinted region CpGs
Manually make diff_imp_loci_pR.txt list by visualizing heatmap
diff_imp_loci_pR <- read.table("diff_imp_loci_pR.txt", header = F, stringsAsFactors = F)
diff_imp_loci_pR_CpGs <-  merge(inormdata1, diff_imp_loci_pR, by.x = "DMR", by.y="V1")
dim(diff_imp_loci_pR_CpGs)
pR_specific_impCpGs <- data.frame(diff_imp_loci_pR_CpGs$TargetID)
colnames(pR_specific_impCpGs) <- c("TargetID")
dim(pR_specific_impCpGs)
length(unique(pR_specific_impCpGs$TargetID))
dim(miMERGE_myiNorm3re_chrpos_avg_PR)
miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs <-  merge(miMERGE_myiNorm3re_chrpos_avg_PR, pR_specific_impCpGs, by.x = "row_4", by.y="TargetID")
head(miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs)
dim(miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs)
cat("Out of",dim(miMERGE_myiNorm3re_chrpos_avg_PR)[1]  ," DMPs in matrixPR and ",dim(pR_specific_impCpGs)[1]," in PR_specific_CpGs only ",dim(miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs)[1]," overlapped. Rest all removed while filtering at closeness of D250 and UN and difference of +/- 0.2 among two Controls v/s PR")



# Extract PCDH CpGs
bedtools intersect -wa -wb  -a MERGE_myiNorm3re_pos_chr.txt -b PCDHcluster.txt > ipcdhdata1.txt
ipcdhdata1 <- read.table("ipcdhdata1.txt", header = F, stringsAsFactors = F)
dim(ipcdhdata1)
length(unique(ipcdhdata1$V4))
#All PCDH region CpGs
write.table(data.frame(ipcdhdata1$V4), "ipcdhdata1_TargetID.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
#fgrep -f ipcdhdata1_TargetID.txt MERGE_myiNorm3re_pos_chravg1.txt -w > MERGE_myiNorm3re_pos_chravg1_pcdh.txt


ipcdhdata1CpGs <- data.frame(ipcdhdata1$V4)
colnames(ipcdhdata1CpGs) <- c("TargetID")
head(miMERGE_myiNorm3re_chrpos_avg_PG)
head(ipcdhdata1CpGs)
#PG specific PCDH region CpGs
miMERGE_myiNorm3re_chrpos_avg_PG_pG_pcdhCpGs <-  merge(miMERGE_myiNorm3re_chrpos_avg_PG, ipcdhdata1CpGs, by.x = "row_4", by.y="TargetID")
head(miMERGE_myiNorm3re_chrpos_avg_PG_pG_pcdhCpGs)
dim(miMERGE_myiNorm3re_chrpos_avg_PG_pG_pcdhCpGs)
cat("Out of",dim(miMERGE_myiNorm3re_chrpos_avg_PG)[1]  ," DMPs in matrixPG and ",dim(ipcdhdata1CpGs)[1]," in PG_specific_CpGs only ",dim(miMERGE_myiNorm3re_chrpos_avg_PG_pG_pcdhCpGs)[1]," overlapped. Rest all removed while filtering at closeness of D250 and UN and difference of +/- 0.2 among two Controls v/s PG")


#PR specific PCDH region CpGs
miMERGE_myiNorm3re_chrpos_avg_PR_pR_pcdhCpGs <-  merge(miMERGE_myiNorm3re_chrpos_avg_PR, ipcdhdata1CpGs, by.x = "row_4", by.y="TargetID")
head(miMERGE_myiNorm3re_chrpos_avg_PR_pR_pcdhCpGs)
dim(miMERGE_myiNorm3re_chrpos_avg_PR_pR_pcdhCpGs)
cat("Out of",dim(miMERGE_myiNorm3re_chrpos_avg_PR)[1]  ," DMPs in matrixPR and ",dim(ipcdhdata1CpGs)[1]," in PR_specific_CpGs only ",dim(miMERGE_myiNorm3re_chrpos_avg_PR_pR_pcdhCpGs)[1]," overlapped. Rest all removed while filtering at closeness of D250 and UN and difference of +/- 0.2 among two Controls v/s PR")



#Mean methylation loss and recovery for all diff meth among cases and corrected clones
#----->All_PG
dim(miMERGE_myiNorm3re_pos_chravg_PG)
head(miMERGE_myiNorm3re_pos_chravg_PG)
library(dplyr)
miMean_myiNorm3re_pos_chravg_PG <- select(miMERGE_myiNorm3re_pos_chravg_PG, PG_Con, C13_Con, C50_Con)
head(miMean_myiNorm3re_pos_chravg_PG)
stmiMean_myiNorm3re_pos_chravg_PG <- stack(as.matrix(miMean_myiNorm3re_pos_chravg_PG))
head(stmiMean_myiNorm3re_pos_chravg_PG)
colnames(stmiMean_myiNorm3re_pos_chravg_PG) <- c("region", "sample", "meandelta_meth")
stmiMean_myiNorm3re_pos_chravg_PG <- data.frame(stmiMean_myiNorm3re_pos_chravg_PG)
ggplot(stmiMean_myiNorm3re_pos_chravg_PG, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#FFB6B3", "#BDE7BD", "#BDE7BD"))+
  theme_bw()+ ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myiNorm3re_pos_chravg_PG.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myiNorm3re_pos_chravg_PG.png", width=10, height=8, units="cm", dpi=96)

#----->All_PR
dim(miMERGE_myiNorm3re_pos_chravg_PR)
head(miMERGE_myiNorm3re_pos_chravg_PR)
miMean_myiNorm3re_pos_chravg_PR <- select(miMERGE_myiNorm3re_pos_chravg_PR, PR_Con, C7_Con, C35_Con)
head(miMean_myiNorm3re_pos_chravg_PR)
stmiMean_myiNorm3re_pos_chravg_PR <- stack(as.matrix(miMean_myiNorm3re_pos_chravg_PR))
head(stmiMean_myiNorm3re_pos_chravg_PR)
colnames(stmiMean_myiNorm3re_pos_chravg_PR) <- c("region", "sample", "meandelta_meth")
stmiMean_myiNorm3re_pos_chravg_PR <- data.frame(stmiMean_myiNorm3re_pos_chravg_PR)
ggplot(stmiMean_myiNorm3re_pos_chravg_PR, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#FFB6B3", "#BDE7BD", "#BDE7BD"))+
  theme_bw() + ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myiNorm3re_pos_chravg_PR.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myiNorm3re_pos_chravg_PR.png", width=10, height=8, units="cm", dpi=96)


#
#For_Imprinted_regions_take only hypomethylated dmrs as hypermethylated are only minorly up
#------>Imprinted_PG
miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs <- data.frame(miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs)
rownames(miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs) <- paste0(miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs$row_1,
                                                                "%",
                                                                miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs$row_2,
                                                                "%",
                                                                miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs$row_3,
                                                                "%",
                                                                miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs$row_4,
                                                                "%",
                                                                "Imprinted")

head(miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs)
miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG <- select(miMERGE_myiNorm3re_chrpos_avg_PG_pG_impCpGs, PG_Con, C13_Con, C50_Con)
head(miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG)
stmiMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG <- stack(as.matrix(miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG))
head(stmiMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG)
colnames(stmiMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG) <- c("region", "sample", "meandelta_meth")
stmiMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG <- data.frame(stmiMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG)

ggplot(stmiMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#FFB6B3", "#BDE7BD", "#BDE7BD"))+
  theme_bw()+ ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG.png", width=10, height=8, units="cm", dpi=96)

#------>Imprinted_PR
miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs <- data.frame(miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs)
rownames(miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs) <- paste0(miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs$row_1,
                                                                "%",
                                                                miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs$row_2,
                                                                "%",
                                                                miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs$row_3,
                                                                "%",
                                                                miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs$row_4,
                                                                "%",
                                                                "Imprinted")

head(miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs)
miMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR <- select(miMERGE_myiNorm3re_chrpos_avg_PR_pR_impCpGs, PR_Con, C7_Con, C35_Con)
head(miMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR)
stmiMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR <- stack(as.matrix(miMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR))
head(stmiMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR)
colnames(stmiMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR) <- c("region", "sample", "meandelta_meth")
stmiMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR <- data.frame(stmiMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR)

ggplot(stmiMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#FFB6B3", "#BDE7BD", "#BDE7BD"))+
  theme_bw()+ ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR.png", width=10, height=8, units="cm", dpi=96)



#------>PCDH_PG
miMERGE_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs <- data.frame(miMERGE_myiNorm3re_chrpos_avg_PG_pG_pcdhCpGs)
rownames(miMERGE_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs) <- paste0(miMERGE_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs$row_1,
                                                                 "%",
                                                                 miMERGE_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs$row_2,
                                                                 "%",
                                                                 miMERGE_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs$row_3,
                                                                 "%",
                                                                 miMERGE_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs$row_4,
                                                                 "%",
                                                                 "PCDH")

head(miMERGE_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs)
miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG <- select(miMERGE_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs, PG_Con, C13_Con, C50_Con)
head(miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG)
stmiMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG <- stack(as.matrix(miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG))
head(stmiMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG)
colnames(stmiMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG) <- c("region", "sample", "meandelta_meth")
stmiMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG <- data.frame(stmiMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG)

ggplot(stmiMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#FFB6B3", "#BDE7BD", "#BDE7BD"))+
  theme_bw()+ ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG.png", width=10, height=8, units="cm", dpi=96)

#------>PCDH_PR
miMERGE_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs <- data.frame(miMERGE_myiNorm3re_chrpos_avg_PR_pR_pcdhCpGs)
rownames(miMERGE_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs) <- paste0(miMERGE_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs$row_1,
                                                                 "%",
                                                                 miMERGE_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs$row_2,
                                                                 "%",
                                                                 miMERGE_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs$row_3,
                                                                 "%",
                                                                 miMERGE_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs$row_4,
                                                                 "%",
                                                                 "PCDH")

head(miMERGE_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs)
miMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR <- select(miMERGE_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs, PR_Con, C7_Con, C35_Con)
head(miMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR)
stmiMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR <- stack(as.matrix(miMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR))
head(stmiMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR)
colnames(stmiMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR) <- c("region", "sample", "meandelta_meth")
stmiMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR <- data.frame(stmiMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR)

ggplot(stmiMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#FFB6B3", "#BDE7BD", "#BDE7BD"))+
  theme_bw()+ ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR.png", width=10, height=8, units="cm", dpi=96)



#Histogram
mirowMean_myiNorm3re_pos_chravg_PG <- data.frame(colMeans(miMean_myiNorm3re_pos_chravg_PG))
mirowMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG <- data.frame(colMeans(miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG))
mirowMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG <- data.frame(colMeans(miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG))
mirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG <- cbind.data.frame(mirowMean_myiNorm3re_pos_chravg_PG,
                                                                         mirowMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG,
                                                                         mirowMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG)

colnames(mirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG) <- c("All_filtered_probes", "diffImprinted", "PCDH")
stmirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG <- stack(as.matrix(mirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG))
stmirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG["samples"] <- c("All_Filt_pG", "All_Filt_zC13","All_Filt_zC50",
                                                                       "Imprinted_pG","Imprinted_zC13","Imprinted_zC50",
                                                                       "PCDH_pG","PCDH_zC13","PCDH_zC50")

head(stmirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG)
stmirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG <- data.frame(stmirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG)
ggplot(stmirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG, aes(x=samples,y=value)) +
  geom_col(aes(fill=samples, color=samples), position = "dodge") +ylim(-1,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = c(rep("darkgrey",3),rep("orange",3),rep("blue",3)))+
  scale_color_manual(values = rep(c("#FFB6B3", "#BDE7BD", "#BDE7BD"), 3))+
  ggtitle("Absolute methylation loss per pG and respective corrected clones")

ggplot(stmirowMean_myiNorm3re_chrpos_avg_PG_pG_impPCDHCpGs_PG, aes(x=samples,y=value)) +
  geom_col(aes(fill=samples, color=samples), position = "dodge") +ylim(-1,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = c(rep("darkgrey",3),rep("orange",3),rep("blue",3)))+
  scale_color_manual(values = rep(c("#FFB6B3", "#BDE7BD", "#BDE7BD"), 3))+
  ggtitle("Absolute methylation loss per pG and respective corrected clones")



#Boxplot
#pG
miMean_myiNorm3re_pos_chravg_PG_all <- miMean_myiNorm3re_pos_chravg_PG
colnames(miMean_myiNorm3re_pos_chravg_PG_all) <- paste0(colnames(miMean_myiNorm3re_pos_chravg_PG_all), "_All")
miMean_myiNorm3re_pos_chravg_PG_all <- stack(as.matrix(miMean_myiNorm3re_pos_chravg_PG_all))
miMean_myiNorm3re_pos_chravg_PG_all <- data.frame(miMean_myiNorm3re_pos_chravg_PG_all)
miMean_myiNorm3re_pos_chravg_PG_all["Color"] <- "All"

miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined <- miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG
colnames(miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined) <- paste0(colnames(miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined), "_diffImprinted")
miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined <- stack(as.matrix(miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined))
miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined <- data.frame(miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined)
miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined["Color"] <- "diffImprinted"

miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh <- miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG
colnames(miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh) <- paste0(colnames(miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh), "_PCDH")
miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh <- stack(as.matrix(miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh))
miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh <- data.frame(miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh)
miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh["Color"] <- "PCDH"
miMean_myiNorm3re_pos_chravg_PG_imp_PCDH <- rbind.data.frame(miMean_myiNorm3re_pos_chravg_PG_all,
                                                             miMean_myiNorm3re_chrpos_avg_PG_pG_impCpGs_PG_imprined,
                                                             miMean_myiNorm3re_chrpos_avg_PG_pG_PCDHCpGs_PG_pcdh)

head(miMean_myiNorm3re_pos_chravg_PG_imp_PCDH)
dim(miMean_myiNorm3re_pos_chravg_PG_imp_PCDH)
ggplot(data=miMean_myiNorm3re_pos_chravg_PG_imp_PCDH, aes(x=col, y=value, fill=Color)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+ coord_flip()+
  theme_minimal()
head(miMean_myiNorm3re_pos_chravg_PG_imp_PCDH)
ggplot(data=miMean_myiNorm3re_pos_chravg_PG_imp_PCDH, aes(x=col, y=value)) +
  geom_boxplot(aes(fill=col),position=position_dodge())+
  ylim(-1,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = rep(c("#FFB6B3", "#BDE7BD", "#BDE7BD"), 3))+
  facet_wrap(~Color, scale="free") +geom_hline(yintercept=0, linetype="dashed", color = "blue")

ggsave("miMean_myiNorm3re_pos_chravg_PG_imp_PCDH.png", width=17*1.25, height=15*1.25, units="cm", dpi=96)


#pR

miMean_myiNorm3re_pos_chravg_PR_all <- miMean_myiNorm3re_pos_chravg_PR
colnames(miMean_myiNorm3re_pos_chravg_PR_all) <- paste0(colnames(miMean_myiNorm3re_pos_chravg_PR_all), "_All")
miMean_myiNorm3re_pos_chravg_PR_all <- stack(as.matrix(miMean_myiNorm3re_pos_chravg_PR_all))
miMean_myiNorm3re_pos_chravg_PR_all <- data.frame(miMean_myiNorm3re_pos_chravg_PR_all)
miMean_myiNorm3re_pos_chravg_PR_all["Color"] <- "All"

miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined <- miMean_myiNorm3re_chrpos_avg_PR_pR_impCpGs_PR
colnames(miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined) <- paste0(colnames(miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined), "_diffImprinted")
miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined <- stack(as.matrix(miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined))
miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined <- data.frame(miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined)
miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined["Color"] <- "diffImprinted"

miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh <- miMean_myiNorm3re_chrpos_avg_PR_pR_PCDHCpGs_PR
colnames(miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh) <- paste0(colnames(miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh), "_PCDH")
miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh <- stack(as.matrix(miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh))
miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh <- data.frame(miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh)
miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh["Color"] <- "PCDH"
miMean_myiNorm3re_pos_chravg_PR_imp_PCDH <- rbind.data.frame(miMean_myiNorm3re_pos_chravg_PR_all,
                                                             miMean_myiNorm3re_chrpos_avg_PR_PR_impCpGs_PR_imprined,
                                                             miMean_myiNorm3re_chrpos_avg_PR_PR_PCDHCpGs_PR_pcdh)

head(miMean_myiNorm3re_pos_chravg_PR_imp_PCDH)
dim(miMean_myiNorm3re_pos_chravg_PR_imp_PCDH)
ggplot(data=miMean_myiNorm3re_pos_chravg_PR_imp_PCDH, aes(x=col, y=value)) +
  geom_boxplot(aes(fill=col),position=position_dodge())+
  ylim(-1,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = rep(c("#FFB6B3", "#BDE7BD", "#BDE7BD"), 3))+
  facet_wrap(~Color, scale="free") +geom_hline(yintercept=0, linetype="dashed", color = "blue")

ggsave("miMean_myiNorm3re_pos_chravg_PR_imp_PCDH.png", width=17*1.25, height=15*1.25, units="cm", dpi=96)


#All DMR
bedtools intersect -wa -wb -a MERGE_myiNorm3re_pos_chr.txt -b /home/ankitv/ref_av/hg38_DMRs/All_known_DMRs_hg38.txt > MERGE_myiNorm3re_human_DMR.txt
awk '{print $16"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_human_DMR.txt > MERGE_myiNorm3re_human_DMR.rearranged.txt

#PCA with  Human ICR
inormsub1 <- read.table("MERGE_myiNorm3re_human_DMR.rearranged.txt", header = F)
colnames(inormsub1) <- c("DMR","chr","start","end","TargetID", colnames(myNorm2))
rownames(inormsub1)
colnames(inormsub1)
head(inormsub1)
dim(inormsub1)
detach("package:dplyr")
Count_inormsubDMR1 <- count(inormsub1, "DMR")
head(Count_inormsubDMR1)
dim(Count_inormsubDMR1)
sdmrcomaggregate1 = aggregate(inormsub1[,6:13],by=list(inormsub1$DMR), mean)
head(sdmrcomaggregate1, 2)
#Plot PCA by group color and labelling
inormsubDMR1=sdmrcomaggregate1
rownames(inormsubDMR1)
inormsubDMR1[,1]
rownames(inormsubDMR1)=inormsubDMR1[,1]
rownames(inormsubDMR1)
colnames(inormsubDMR1)
inormsubDMR1 = inormsubDMR1[,-1]
head(inormsubDMR1)
dim(inormsubDMR1)
inormsubDMR1_counted <- cbind(inormsubDMR1, Count_inormsubDMR1)
head(inormsubDMR1_counted)
tail(inormsubDMR1_counted)
dim(inormsubDMR1_counted)
inormsubDMR1_counted1 <- inormsubDMR1_counted[which(inormsubDMR1_counted$freq >= 3),]
head(inormsubDMR1_counted1)
dim(inormsubDMR1_counted1)
write.table(inormsubDMR1_counted1, "inormsubDMR1_counted1.filtmin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)


inormsubDMR1_counted2 <- inormsubDMR1_counted1[,1:8]
head(inormsubDMR1_counted2)
summary(inormsubDMR1_counted2)

#Heatmap of control normalized Log other samples merged
head(inormsubDMR1_counted2)
dim(inormsubDMR1_counted2)
#rearrange columns
inormsubDMR1_counted2ratio <- as.matrix(inormsubDMR1_counted2[,c(5,8,6,7,1:4)])

head(inormsubDMR1_counted2ratio)
dim(inormsubDMR1_counted2ratio)
winormsubDMR1_counted2ratioavg <- data.frame(cbind((rowMeans(inormsubDMR1_counted2ratio[,1:2])),
                                                   (inormsubDMR1_counted2ratio[,1:8])))

head(winormsubDMR1_counted2ratioavg)
colnames(winormsubDMR1_counted2ratioavg) <- c("Allcontrol",colnames(inormsubDMR1_counted2ratio))
head(winormsubDMR1_counted2ratioavg)
winormsubDMR1_counted2ratioavg1 <- winormsubDMR1_counted2ratioavg
winormsubDMR1_counted2ratioavg2 <- as.matrix(winormsubDMR1_counted2ratioavg1)
head(winormsubDMR1_counted2ratioavg2)
dim(winormsubDMR1_counted2ratioavg2)

minormsubDMR1_counted2ratioavg3 <- data.frame(cbind((winormsubDMR1_counted2ratioavg2[,1:9]-winormsubDMR1_counted2ratioavg2[,1])))


head(minormsubDMR1_counted2ratioavg3)
minormsubDMR1_counted2ratioavg3 <- as.matrix(minormsubDMR1_counted2ratioavg3)

minormsubDMR1_counted2ratioavg3minuscontrol <- minormsubDMR1_counted2ratioavg3[,c(2:9)]
head(minormsubDMR1_counted2ratioavg3minuscontrol)
library(pheatmap)
my_sample_col7 <- data.frame(Annotations= colnames(minormsubDMR1_counted2ratioavg3minuscontrol))
row.names(my_sample_col7) <- colnames(minormsubDMR1_counted2ratioavg3minuscontrol)
my_colour7 = list(Annotations = c("D250"= "grey","UN"= "grey","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(minormsubDMR1_counted2ratioavg3minuscontrol,
         color = colorRampPalette(c("navy", "white", "#D47400"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour7,
         fontsize = 8,
         annotation_col = my_sample_col7,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "black",
         cutree_cols = 2)

#save as pHeatmap_minormsubDMR1_counted2ratioavg.normcontrol.png



#UCSC
#Prepare files for UCSC:
head(MERGE_myiNorm3re_pos_chr)
dim(MERGE_myiNorm3re_pos_chr)

write.table(MERGE_myiNorm3re_pos_chr[,c(1:3,5)], "ucsc.C7.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiNorm3re_pos_chr[,c(1:3,6)], "ucsc.C13.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiNorm3re_pos_chr[,c(1:3,7)], "ucsc.C35.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiNorm3re_pos_chr[,c(1:3,8)], "ucsc.C50.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiNorm3re_pos_chr[,c(1:3,9)], "ucsc.D250.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiNorm3re_pos_chr[,c(1:3,10)], "ucsc.PG.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiNorm3re_pos_chr[,c(1:3,11)], "ucsc.PR.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiNorm3re_pos_chr[,c(1:3,12)], "ucsc.UN.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)

sed '1s/^/track type=bedGraph name="C7_beta_values" description="C7_beta_values" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=51,255,51 \n/' ucsc.C7.hg38.chr.bed > ucsc.C7.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="C13_beta_values" description="C13_beta_values" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=51,255,51 \n/' ucsc.C13.hg38.chr.bed > ucsc.C13.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="C35_beta_values" description="C35_beta_values" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=51,255,51 \n/' ucsc.C35.hg38.chr.bed > ucsc.C35.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="C50_beta_values" description="C50_beta_values" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=51,255,51 \n/' ucsc.C50.hg38.chr.bed > ucsc.C50.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="D250_beta_values" description="D250_beta_values" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc.D250.hg38.chr.bed > ucsc.D250.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="PG_beta_values" description="PG_beta_values" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=255,0,0 \n/' ucsc.PG.hg38.chr.bed > ucsc.PG.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="PR_beta_values" description="PR_beta_values" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=255,0,0 \n/' ucsc.PR.hg38.chr.bed > ucsc.PR.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="UN_beta_values" description="UN_beta_values" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc.UN.hg38.chr.bed > ucsc.UN.hg38.chr.bedGraph


#Peek into genome wide view
bedtools makewindows -g /home/ankitv/ref_av/hg38/hg38.chrom.sizes -w 1000 > hg38_1000.bed


#Bin1Kb 

bedtools intersect -wa -wb -a MERGE_myiNorm3re_pos_chr.txt -b hg38_1000.bed > MERGE_myiNorm3re_human_Bin1Kb.txt
awk '{print $13"%"$14"%"$15"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_human_Bin1Kb.txt > MERGE_myiNorm3re_human_Bin1Kb.rearranged.txt

awk '{print $13"\t"$14"\t"$15"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_human_Bin1Kb.txt | sort -k1,1 -k2,2n > MERGE_myiNorm3re_human_Bin1Kb.rearranged_prechr.txt

#PCA with  Human Bin1Kb
inormdataBin1Kb <- read.table("MERGE_myiNorm3re_human_Bin1Kb.rearranged.txt", header = F, stringsAsFactors = F)
colnames(inormdataBin1Kb) <- c("Bin1Kb","chr","start","end","TargetID", colnames(myNorm2))
rownames(inormdataBin1Kb)
colnames(inormdataBin1Kb)
head(inormdataBin1Kb)
dim(inormdataBin1Kb)
Count_inormdataBin1Kb1 <- count(inormdataBin1Kb, "Bin1Kb")
head(Count_inormdataBin1Kb1)
dim(Count_inormdataBin1Kb1)
snormaggregateBin1Kb = aggregate(inormdataBin1Kb[,6:13],by=list(inormdataBin1Kb$Bin1Kb), mean)
head(snormaggregateBin1Kb, 2)
#Plot PCA by group color and labelling
inormdataBin1Kb1=snormaggregateBin1Kb
rownames(inormdataBin1Kb1)
inormdataBin1Kb1[,1]
rownames(inormdataBin1Kb1)=inormdataBin1Kb1[,1]
rownames(inormdataBin1Kb1)
colnames(inormdataBin1Kb1)
inormdataBin1Kb1 = inormdataBin1Kb1[,-1]
head(inormdataBin1Kb1)
dim(inormdataBin1Kb1)
inormdataBin1Kb1_counted <- cbind(inormdataBin1Kb1, Count_inormdataBin1Kb1)
head(inormdataBin1Kb1_counted)
dim(inormdataBin1Kb1_counted)
tail(inormdataBin1Kb1_counted)
inormdataBin1Kb1_counted1 <- inormdataBin1Kb1_counted[which(inormdataBin1Kb1_counted$freq >= 5),]
head(inormdataBin1Kb1_counted1)
dim(inormdataBin1Kb1_counted1)
write.table(inormdataBin1Kb1_counted1, "inormdataBin1Kb1_counted1.filtmin.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

inormdataBin1Kb1_counted1chr <- cSplit(inormdataBin1Kb1_counted1, "Bin1Kb","%")
head(inormdataBin1Kb1_counted1chr)
inormdataBin1Kb1_counted1chr <- inormdataBin1Kb1_counted1chr[,c(10:12,1:9)]
write.table(inormdataBin1Kb1_counted1chr, "inormdataBin1Kb1_counted1filtminchr.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
sort -k1,1 -k2,2n inormdataBin1Kb1_counted1filtminchr.txt > inormdataBin1Kb1_counted1filtminchr.sort.txt
bedtools closest -a inormdataBin1Kb1_counted1filtminchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > inormdataBin1Kb1_counted1filtminchr_genehg38.txt

inormdataBin1Kb1_counted2 <- inormdataBin1Kb1_counted1[,1:8]
head(inormdataBin1Kb1_counted2)
summary(inormdataBin1Kb1_counted2)
#df <- as.inormdata.frame(inormdataBin1Kb1)
inormdfBin1Kb <- inormdataBin1Kb1_counted2
head(inormdfBin1Kb)
dim(inormdfBin1Kb)
inormdfBin1Kb = data.frame(inormdfBin1Kb)
#write.table(inormdfBin1Kb , "inormdfBin1Kbdedup.txt", sep="\t", quote = FALSE, append = FALSE)
tinormdfBin1Kb = t(inormdfBin1Kb)
tinormdfBin1Kb = data.frame(tinormdfBin1Kb)
#write.table(tinormdfBin1Kb , "tinormdfBin1Kbdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfBin1Kb)
tinormdfBin1Kb["Color"] <-  c("corrClone","corrClone","corrClone","corrClone","Control","Case","Case","Control")
#head(tinormdfBin1Kb)
#write.table(tinormdfBin1Kb , "tinormdfBin1Kbdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfBin1Kb)
inormdfx <-tinormdfBin1Kb[c(1:(length(tinormdfBin1Kb)-1))]
PBin1KbC<-prcomp(inormdfx, center = TRUE, scale. = TRUE)
PBin1KbCi<-data.frame(PBin1KbC$x,Color=tinormdfBin1Kb$Color)
percentagecomBin1Kb <- round(PBin1KbC$sdev^2 / sum(PBin1KbC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomBin1Kb <- paste( colnames(PBin1KbCi), "(", paste( as.character(percentagecomBin1Kb), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pBin1Kb1 <-ggplot(PBin1KbCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomBin1Kb[1]) + ylab(percentagecomBin1Kb[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("#FFB6B3","darkgrey","#BDE7BD"))+
  scale_shape_manual(values=c(0,1,2))
pBin1Kb1 <- pBin1Kb1 +theme_classic()
pBin1Kb1 + xlim(-500,500)+ ylim(-500,500)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_tinormdfBin1Kb_CntrlIndiv.svg", width=10*1.25, height=8*1.25, units="cm", dpi=96)


################################### Heatmap_indiv human Bin1Kb ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(inormdataBin1Kb1_counted1)
dim(inormdataBin1Kb1_counted1)
sum(inormdataBin1Kb1_counted1$freq)

#Heatmap_indiv with 2 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
inormdataBin1Kb2 = as.matrix(inormdataBin1Kb1_counted2)
head(inormdataBin1Kb2)
dim(inormdataBin1Kb2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
inormdataBin1Kb3 <- inormdataBin1Kb2
head(inormdataBin1Kb3)
write.table(inormdataBin1Kb3, "inormdataBin1Kb3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

#Box and Violin Plot human Bin1Kb using 12 Controls from MethHed
rownames(inormdataBin1Kb)
colnames(inormdataBin1Kb)
head(inormdataBin1Kb)
dim(inormdataBin1Kb)
Count_inormdataBin1Kb1 <- count(inormdataBin1Kb, "Bin1Kb")
head(Count_inormdataBin1Kb1)
inormagregateBin1Kb1 = aggregate(inormdataBin1Kb[,6:13],by=list(inormdataBin1Kb$Bin1Kb), mean)
head(inormagregateBin1Kb1, 2)
inormdataBin1Kb1=inormagregateBin1Kb1
rownames(inormdataBin1Kb1)
inormdataBin1Kb1[,1]
rownames(inormdataBin1Kb1)=inormdataBin1Kb1[,1]
rownames(inormdataBin1Kb1)
colnames(inormdataBin1Kb1)
inormdataBin1Kb1 = inormdataBin1Kb1[,-1]
head(inormdataBin1Kb1)
dim(inormdataBin1Kb1)
write.table(inormdataBin1Kb1, "inormaggregatedBin1Kb_inormdataBin1Kb1.txt", sep="\t", quote = FALSE, append = FALSE)
inormdataBin1Kb1_counted <- cbind(inormdataBin1Kb1, Count_inormdataBin1Kb1)
head(inormdataBin1Kb1_counted)
dim(inormdataBin1Kb1_counted)
write.table(inormdataBin1Kb1_counted, "inormaggregatedBin1Kb_counted.txt", sep="\t", quote = FALSE, append = FALSE)

inormdataBin1Kb1_counted1 <- inormdataBin1Kb1_counted[which(inormdataBin1Kb1_counted$freq >= 5),]
head(inormdataBin1Kb1_counted1)
dim(inormdataBin1Kb1_counted1)
inormdataBin1Kb1_counted2 <- inormdataBin1Kb1_counted1[,1:8]
head(inormdataBin1Kb1_counted2)
VionormIndvBin1Kb <- data.frame(inormdataBin1Kb1_counted2[,c(5,8,6,7,1:4)])
head(VionormIndvBin1Kb)
write.table(VionormIndvBin1Kb, "VionormIndvBin1Kb.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_VionormIndvBin1Kb.svg", width=10, height=5, pointsize=12)
boxplot(VionormIndvBin1Kb, main="Average methylation at human Bin1Kbs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","navy","navy","navy","navy"))
dev.off()
dim(VionormIndvBin1Kb)
VionormIndvBin1Kb <- data.frame(VionormIndvBin1Kb)
VionormIndvBin1Kb1 <- VionormIndvBin1Kb
head(VionormIndvBin1Kb1,1)
VionormIndvBin1Kb2 <- stack(VionormIndvBin1Kb1)
head(VionormIndvBin1Kb2)
colnames(VionormIndvBin1Kb2) <- c("Methylation", "inormdatasets")
head(VionormIndvBin1Kb2)
ggplot(VionormIndvBin1Kb2, aes(x=inormdatasets, y=Methylation, color=inormdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","#FFB6B3","#FFB6B3","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","navy","navy","navy"))
ggsave("Violin_plot_VionormIndvBin1Kb2_pimethhed_humanBin1Kb.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)



#Heatmap of control normalized Log other samples merged
head(inormdataBin1Kb1_counted2)
dim(inormdataBin1Kb1_counted2)
inormdataBin1Kb1_counted2ratio <- as.matrix(inormdataBin1Kb1_counted2[,c(5,8,6,7,1:4)])

head(inormdataBin1Kb1_counted2ratio)
dim(inormdataBin1Kb1_counted2ratio)
winormdataBin1Kb1_counted2ratioavg <- data.frame(cbind((rowMeans(inormdataBin1Kb1_counted2ratio[,1:2])),
                                                          (inormdataBin1Kb1_counted2ratio[,1:8])))

head(winormdataBin1Kb1_counted2ratioavg)
colnames(winormdataBin1Kb1_counted2ratioavg) <- c("Allcontrol",colnames(inormdataBin1Kb1_counted2ratio))
head(winormdataBin1Kb1_counted2ratioavg)
winormdataBin1Kb1_counted2ratioavg1 <- winormdataBin1Kb1_counted2ratioavg
winormdataBin1Kb1_counted2ratioavg2 <- as.matrix(winormdataBin1Kb1_counted2ratioavg1)
head(winormdataBin1Kb1_counted2ratioavg2)
dim(winormdataBin1Kb1_counted2ratioavg2)
winormdataBin1Kb1_counted2ratioavg3 <- data.frame(cbind((winormdataBin1Kb1_counted2ratioavg2[,1]/winormdataBin1Kb1_counted2ratioavg2[,1]),
                                                           (winormdataBin1Kb1_counted2ratioavg2[,2]/winormdataBin1Kb1_counted2ratioavg2[,1]),
                                                           (winormdataBin1Kb1_counted2ratioavg2[,3]/winormdataBin1Kb1_counted2ratioavg2[,1]),
                                                           (winormdataBin1Kb1_counted2ratioavg2[,4]/winormdataBin1Kb1_counted2ratioavg2[,1]),
                                                           (winormdataBin1Kb1_counted2ratioavg2[,5]/winormdataBin1Kb1_counted2ratioavg2[,1]),
                                                           (winormdataBin1Kb1_counted2ratioavg2[,6]/winormdataBin1Kb1_counted2ratioavg2[,1]),
                                                           (winormdataBin1Kb1_counted2ratioavg2[,7]/winormdataBin1Kb1_counted2ratioavg2[,1]),
                                                           (winormdataBin1Kb1_counted2ratioavg2[,8]/winormdataBin1Kb1_counted2ratioavg2[,1]),
                                                           (winormdataBin1Kb1_counted2ratioavg2[,9]/winormdataBin1Kb1_counted2ratioavg2[,1])))


head(winormdataBin1Kb1_counted2ratioavg3)
colnames(winormdataBin1Kb1_counted2ratioavg3) <- colnames(winormdataBin1Kb1_counted2ratioavg2)
head(winormdataBin1Kb1_counted2ratioavg3)
winormdataBin1Kb1_counted2ratioavg3 <- as.matrix(winormdataBin1Kb1_counted2ratioavg3)
winormdataBin1Kb1_counted2ratioavg3Log <- log2(winormdataBin1Kb1_counted2ratioavg3)
head(winormdataBin1Kb1_counted2ratioavg3Log)
write.table(winormdataBin1Kb1_counted2ratioavg3Log , "Heatmap_ICF_normcontrolwinormdataBin1Kb1_counted2ratioavg2Log.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)

winormdataBin1Kb1_counted2ratioavg3Logminuscontrol <- winormdataBin1Kb1_counted2ratioavg3Log[,c(2:9)]
head(winormdataBin1Kb1_counted2ratioavg3Logminuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(winormdataBin1Kb1_counted2ratioavg3Logminuscontrol))
row.names(my_sample_col2) <- colnames(winormdataBin1Kb1_counted2ratioavg3Logminuscontrol)
my_colour2 = list(Annotations = c("D250"= "darkgrey","UN"= "darkgrey","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(winormdataBin1Kb1_counted2ratioavg3Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_winormdataBin1Kb1_counted2ratioavg.normcontrol.svg

minormdataBin1Kb1_counted2ratioavg3 <- data.frame(cbind((winormdataBin1Kb1_counted2ratioavg2[,1:9]-winormdataBin1Kb1_counted2ratioavg2[,1])))


head(minormdataBin1Kb1_counted2ratioavg3)
minormdataBin1Kb1_counted2ratioavg3 <- as.matrix(minormdataBin1Kb1_counted2ratioavg3)

minormdataBin1Kb1_counted2ratioavg3minuscontrol <- minormdataBin1Kb1_counted2ratioavg3[,c(2:9)]
head(minormdataBin1Kb1_counted2ratioavg3minuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(minormdataBin1Kb1_counted2ratioavg3minuscontrol))
row.names(my_sample_col2) <- colnames(minormdataBin1Kb1_counted2ratioavg3minuscontrol)
my_colour2 = list(Annotations = c("D250"= "darkgrey","UN"= "darkgrey","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(minormdataBin1Kb1_counted2ratioavg3minuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataBin1Kb1_counted2ratioavg.normcontrol.svg

pheatmap(minormdataBin1Kb1_counted2ratioavg3minuscontrol[,1:4],
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataBin1Kb1_counted2ratioavg.normcontrol_casecntrl.svg

#Differentially methylated probes specifc to cases
pinormdataBin1Kb1_counted2ratioavg2 <- data.frame(winormdataBin1Kb1_counted2ratioavg2)
head(pinormdataBin1Kb1_counted2ratioavg2)
#Subtraction, Control -Control
pinormdataBin1Kb1_counted2ratioavg2["D250_UN"] <- pinormdataBin1Kb1_counted2ratioavg2$D250 - pinormdataBin1Kb1_counted2ratioavg2$UN

#Subtraction, Case -Control
pinormdataBin1Kb1_counted2ratioavg2["PG_Con"] <- pinormdataBin1Kb1_counted2ratioavg2$PG - pinormdataBin1Kb1_counted2ratioavg2$Allcontrol
pinormdataBin1Kb1_counted2ratioavg2["PR_Con"] <- pinormdataBin1Kb1_counted2ratioavg2$PR - pinormdataBin1Kb1_counted2ratioavg2$Allcontrol
pinormdataBin1Kb1_counted2ratioavg2["C7_Con"] <- pinormdataBin1Kb1_counted2ratioavg2$C7 - pinormdataBin1Kb1_counted2ratioavg2$Allcontrol
pinormdataBin1Kb1_counted2ratioavg2["C13_Con"] <- pinormdataBin1Kb1_counted2ratioavg2$C13 - pinormdataBin1Kb1_counted2ratioavg2$Allcontrol
pinormdataBin1Kb1_counted2ratioavg2["C35_Con"] <- pinormdataBin1Kb1_counted2ratioavg2$C35 - pinormdataBin1Kb1_counted2ratioavg2$Allcontrol
pinormdataBin1Kb1_counted2ratioavg2["C50_Con"] <- pinormdataBin1Kb1_counted2ratioavg2$C50 - pinormdataBin1Kb1_counted2ratioavg2$Allcontrol
pinormdataBin1Kb1_counted2ratioavg2["PG_D250"] <- pinormdataBin1Kb1_counted2ratioavg2$PG - pinormdataBin1Kb1_counted2ratioavg2$D250
pinormdataBin1Kb1_counted2ratioavg2["PR_D250"] <- pinormdataBin1Kb1_counted2ratioavg2$PR - pinormdataBin1Kb1_counted2ratioavg2$D250
pinormdataBin1Kb1_counted2ratioavg2["C7_D250"] <- pinormdataBin1Kb1_counted2ratioavg2$C7 - pinormdataBin1Kb1_counted2ratioavg2$D250
pinormdataBin1Kb1_counted2ratioavg2["C13_D250"] <- pinormdataBin1Kb1_counted2ratioavg2$C13 - pinormdataBin1Kb1_counted2ratioavg2$D250
pinormdataBin1Kb1_counted2ratioavg2["C35_D250"] <- pinormdataBin1Kb1_counted2ratioavg2$C35 - pinormdataBin1Kb1_counted2ratioavg2$D250
pinormdataBin1Kb1_counted2ratioavg2["C50_D250"] <- pinormdataBin1Kb1_counted2ratioavg2$C50 - pinormdataBin1Kb1_counted2ratioavg2$D250
pinormdataBin1Kb1_counted2ratioavg2["PG_UN"] <- pinormdataBin1Kb1_counted2ratioavg2$PG - pinormdataBin1Kb1_counted2ratioavg2$UN
pinormdataBin1Kb1_counted2ratioavg2["PR_UN"] <- pinormdataBin1Kb1_counted2ratioavg2$PR - pinormdataBin1Kb1_counted2ratioavg2$UN
pinormdataBin1Kb1_counted2ratioavg2["C7_UN"] <- pinormdataBin1Kb1_counted2ratioavg2$C7 - pinormdataBin1Kb1_counted2ratioavg2$UN
pinormdataBin1Kb1_counted2ratioavg2["C13_UN"] <- pinormdataBin1Kb1_counted2ratioavg2$C13 - pinormdataBin1Kb1_counted2ratioavg2$UN
pinormdataBin1Kb1_counted2ratioavg2["C35_UN"] <- pinormdataBin1Kb1_counted2ratioavg2$C35 - pinormdataBin1Kb1_counted2ratioavg2$UN
pinormdataBin1Kb1_counted2ratioavg2["C50_UN"] <- pinormdataBin1Kb1_counted2ratioavg2$C50 - pinormdataBin1Kb1_counted2ratioavg2$UN
pinormdataBin1Kb1_counted2ratioavg2["C7_PR"] <- pinormdataBin1Kb1_counted2ratioavg2$C7 - pinormdataBin1Kb1_counted2ratioavg2$PR
pinormdataBin1Kb1_counted2ratioavg2["C13_PG"] <- pinormdataBin1Kb1_counted2ratioavg2$C13 - pinormdataBin1Kb1_counted2ratioavg2$PG
pinormdataBin1Kb1_counted2ratioavg2["C35_PR"] <- pinormdataBin1Kb1_counted2ratioavg2$C35 - pinormdataBin1Kb1_counted2ratioavg2$PR
pinormdataBin1Kb1_counted2ratioavg2["C50_PG"] <- pinormdataBin1Kb1_counted2ratioavg2$C50 - pinormdataBin1Kb1_counted2ratioavg2$PG

#Region where two controls are close
pinormdataBin1Kb1_counted2ratioavg2 <- pinormdataBin1Kb1_counted2ratioavg2[which(pinormdataBin1Kb1_counted2ratioavg2$D250_UN < 0.1 & 
                                                                                   pinormdataBin1Kb1_counted2ratioavg2$D250_UN > -0.1),]

dim(pinormdataBin1Kb1_counted2ratioavg2)
#PG
pinormdataBin1Kb1_counted2ratioavg2_PG <- pinormdataBin1Kb1_counted2ratioavg2[which(pinormdataBin1Kb1_counted2ratioavg2$PG_D250 > 0.2 | 
                                                                                      pinormdataBin1Kb1_counted2ratioavg2$PG_D250 < -0.2),]


mpinormdataBin1Kb1_counted2ratioavg2_PG <- pinormdataBin1Kb1_counted2ratioavg2_PG[which(pinormdataBin1Kb1_counted2ratioavg2_PG$PG_UN > 0.2 | 
                                                                                          pinormdataBin1Kb1_counted2ratioavg2_PG$PG_UN < -0.2),]

head(mpinormdataBin1Kb1_counted2ratioavg2_PG)
dim(mpinormdataBin1Kb1_counted2ratioavg2_PG)

write.table(rownames(mpinormdataBin1Kb1_counted2ratioavg2_PG) , "mpinormdataBin1Kb1_counted2ratioavg2_PGchr.regions.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = F)

mpinormdataBin1Kb1_counted2ratioavg2_PGchr <- data.frame(mpinormdataBin1Kb1_counted2ratioavg2_PG)
mpinormdataBin1Kb1_counted2ratioavg2_PGchr["pos"] <- rownames(mpinormdataBin1Kb1_counted2ratioavg2_PGchr)
mpinormdataBin1Kb1_counted2ratioavg2_PGchr <- cSplit(mpinormdataBin1Kb1_counted2ratioavg2_PGchr, "pos", "%")
head(mpinormdataBin1Kb1_counted2ratioavg2_PGchr)
dim(mpinormdataBin1Kb1_counted2ratioavg2_PGchr)
mpinormdataBin1Kb1_counted2ratioavg2_PGchr <- mpinormdataBin1Kb1_counted2ratioavg2_PGchr[,c(33:35,1:2)]
write.table(mpinormdataBin1Kb1_counted2ratioavg2_PGchr , "mpinormdataBin1Kb1_counted2ratioavg2_PGchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = T)


#get beta-matrix on hypo and hyper DMR PG
fgrep -f mpinormdataBin1Kb1_counted2ratioavg2_PGchr.regions.txt inormaggregatedBin1Kb_counted.txt  -w > inormaggregatedBin1Kb_countedBin1kb_PG.txt
inormaggregatedBin1Kb_countedBin1kb_PG <- read.table("inormaggregatedBin1Kb_countedBin1kb_PG.txt", header = F, stringsAsFactors = F, row.names = 1)
head(inormaggregatedBin1Kb_countedBin1kb_PG)
colnames(inormaggregatedBin1Kb_countedBin1kb_PG) <- colnames(inormdataBin1Kb1_counted)
inormaggregatedBin1Kb_countedBin1kb_PGchr <- data.frame(inormaggregatedBin1Kb_countedBin1kb_PG)
inormaggregatedBin1Kb_countedBin1kb_PGchr["pos"] <- rownames(inormaggregatedBin1Kb_countedBin1kb_PGchr)
inormaggregatedBin1Kb_countedBin1kb_PGchr <- cSplit(inormaggregatedBin1Kb_countedBin1kb_PGchr, "pos", "%")
head(inormaggregatedBin1Kb_countedBin1kb_PGchr)
dim(inormaggregatedBin1Kb_countedBin1kb_PGchr)
inormaggregatedBin1Kb_countedBin1kb_PGchr <- inormaggregatedBin1Kb_countedBin1kb_PGchr[,c(11:13,9:10,1:8)]
write.table(inormaggregatedBin1Kb_countedBin1kb_PGchr , "inormaggregatedBin1Kb_countedBin1kb_PGchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = T)
rownames(inormaggregatedBin1Kb_countedBin1kb_PGchr)  <- inormaggregatedBin1Kb_countedBin1kb_PGchr$Bin1Kb
pheatmap(inormaggregatedBin1Kb_countedBin1kb_PGchr[,6:13])
sort -k1,1 -k2,2n inormaggregatedBin1Kb_countedBin1kb_PGchr.txt | grep chr > inormaggregatedBin1Kb_countedBin1kb_PGchr.sort.txt
sort -k1,1 -k2,2n mpinormdataBin1Kb1_counted2ratioavg2_PGchr.txt | grep chr > mpinormdataBin1Kb1_counted2ratioavg2_PGchr.sort.txt

bedtools closest -wa -wb -a mpinormdataBin1Kb1_counted2ratioavg2_PGchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > mpinormdataBin1Kb1_counted2ratioavg2_PGchr_genehg38.txt
bedtools closest -wa -wb -a inormaggregatedBin1Kb_countedBin1kb_PGchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d  > inormaggregatedBin1Kb_countedBin1kb_PGchr_genehg38.txt




#PR
pinormdataBin1Kb1_counted2ratioavg2_PR <- pinormdataBin1Kb1_counted2ratioavg2[which(pinormdataBin1Kb1_counted2ratioavg2$PR_D250 > 0.2 | 
                                                                                      pinormdataBin1Kb1_counted2ratioavg2$PR_D250 < -0.2),]


mpinormdataBin1Kb1_counted2ratioavg2_PR <- pinormdataBin1Kb1_counted2ratioavg2_PR[which(pinormdataBin1Kb1_counted2ratioavg2_PR$PR_UN > 0.2 | 
                                                                                          pinormdataBin1Kb1_counted2ratioavg2_PR$PR_UN < -0.2),]

head(mpinormdataBin1Kb1_counted2ratioavg2_PR)
dim(mpinormdataBin1Kb1_counted2ratioavg2_PR)

write.table(rownames(mpinormdataBin1Kb1_counted2ratioavg2_PR) , "mpinormdataBin1Kb1_counted2ratioavg2_PRchr.regions.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = F)

mpinormdataBin1Kb1_counted2ratioavg2_PRchr <- data.frame(mpinormdataBin1Kb1_counted2ratioavg2_PR)
mpinormdataBin1Kb1_counted2ratioavg2_PRchr["pos"] <- rownames(mpinormdataBin1Kb1_counted2ratioavg2_PRchr)
mpinormdataBin1Kb1_counted2ratioavg2_PRchr <- cSplit(mpinormdataBin1Kb1_counted2ratioavg2_PRchr, "pos", "%")
head(mpinormdataBin1Kb1_counted2ratioavg2_PRchr)
dim(mpinormdataBin1Kb1_counted2ratioavg2_PRchr)
mpinormdataBin1Kb1_counted2ratioavg2_PRchr <- mpinormdataBin1Kb1_counted2ratioavg2_PRchr[,c(33:35,1:32)]
write.table(mpinormdataBin1Kb1_counted2ratioavg2_PRchr , "mpinormdataBin1Kb1_counted2ratioavg2_PRchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = T)


#get beta-matrix on hypo and hyper DMR PR
fgrep -f mpinormdataBin1Kb1_counted2ratioavg2_PRchr.regions.txt inormaggregatedBin1Kb_counted.txt  -w > inormaggregatedBin1Kb_countedBin1kb_PR.txt
inormaggregatedBin1Kb_countedBin1kb_PR <- read.table("inormaggregatedBin1Kb_countedBin1kb_PR.txt", header = F, stringsAsFactors = F, row.names = 1)
head(inormaggregatedBin1Kb_countedBin1kb_PR)
colnames(inormaggregatedBin1Kb_countedBin1kb_PR) <- colnames(inormdataBin1Kb1_counted)
inormaggregatedBin1Kb_countedBin1kb_PRchr <- data.frame(inormaggregatedBin1Kb_countedBin1kb_PR)
inormaggregatedBin1Kb_countedBin1kb_PRchr["pos"] <- rownames(inormaggregatedBin1Kb_countedBin1kb_PRchr)
inormaggregatedBin1Kb_countedBin1kb_PRchr <- cSplit(inormaggregatedBin1Kb_countedBin1kb_PRchr, "pos", "%")
head(inormaggregatedBin1Kb_countedBin1kb_PRchr)
dim(inormaggregatedBin1Kb_countedBin1kb_PRchr)
inormaggregatedBin1Kb_countedBin1kb_PRchr <- inormaggregatedBin1Kb_countedBin1kb_PRchr[,c(11:13,9:10,1:8)]
write.table(inormaggregatedBin1Kb_countedBin1kb_PRchr , "inormaggregatedBin1Kb_countedBin1kb_PRchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = T)
rownames(inormaggregatedBin1Kb_countedBin1kb_PRchr)  <- inormaggregatedBin1Kb_countedBin1kb_PRchr$Bin1Kb
pheatmap(inormaggregatedBin1Kb_countedBin1kb_PRchr[,6:13])
sort -k1,1 -k2,2n inormaggregatedBin1Kb_countedBin1kb_PRchr.txt | grep chr > inormaggregatedBin1Kb_countedBin1kb_PRchr.sort.txt
sort -k1,1 -k2,2n mpinormdataBin1Kb1_counted2ratioavg2_PRchr.txt | grep chr > mpinormdataBin1Kb1_counted2ratioavg2_PRchr.sort.txt

bedtools closest -wa -wb -a mpinormdataBin1Kb1_counted2ratioavg2_PRchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > mpinormdataBin1Kb1_counted2ratioavg2_PRchr_genehg38.txt
bedtools closest -wa -wb -a inormaggregatedBin1Kb_countedBin1kb_PRchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d  > inormaggregatedBin1Kb_countedBin1kb_PRchr_genehg38.txt




#Bin 500 bp

#Peek into genome wide view
bedtools makewindows -g /home/ankitv/ref_av/hg38/hg38.chrom.sizes -w 500 > hg38_500.bed

#Bin500bp 
bedtools intersect -wa -wb -a MERGE_myiNorm3re_pos_chr.txt -b hg38_500.bed > MERGE_myiNorm3re_human_Bin500bp.txt
awk '{print $13"%"$14"%"$15"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_human_Bin500bp.txt > MERGE_myiNorm3re_human_Bin500bp.rearranged.txt

awk '{print $13"\t"$14"\t"$15"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_human_Bin500bp.txt | sort -k1,1 -k2,2n > MERGE_myiNorm3re_human_Bin500bp.rearranged_prechr.txt

#PCA with  Human Bin500bp
inormdataBin500bp <- read.table("MERGE_myiNorm3re_human_Bin500bp.rearranged.txt", header = F, stringsAsFactors = F)
colnames(inormdataBin500bp) <- c("Bin500bp","chr","start","end","TargetID", colnames(myNorm2))
rownames(inormdataBin500bp)
colnames(inormdataBin500bp)
head(inormdataBin500bp)
dim(inormdataBin500bp)
Count_inormdataBin500bp1 <- count(inormdataBin500bp, "Bin500bp")
head(Count_inormdataBin500bp1)
dim(Count_inormdataBin500bp1)
snormaggregateBin500bp = aggregate(inormdataBin500bp[,6:13],by=list(inormdataBin500bp$Bin500bp), mean)
head(snormaggregateBin500bp, 2)
#Plot PCA by group color and labelling
inormdataBin500bp1=snormaggregateBin500bp
rownames(inormdataBin500bp1)
inormdataBin500bp1[,1]
rownames(inormdataBin500bp1)=inormdataBin500bp1[,1]
rownames(inormdataBin500bp1)
colnames(inormdataBin500bp1)
inormdataBin500bp1 = inormdataBin500bp1[,-1]
head(inormdataBin500bp1)
dim(inormdataBin500bp1)
inormdataBin500bp1_counted <- cbind(inormdataBin500bp1, Count_inormdataBin500bp1)
head(inormdataBin500bp1_counted)
dim(inormdataBin500bp1_counted)
tail(inormdataBin500bp1_counted)
inormdataBin500bp1_counted1 <- inormdataBin500bp1_counted[which(inormdataBin500bp1_counted$freq >= 3),]
head(inormdataBin500bp1_counted1)
dim(inormdataBin500bp1_counted1)
write.table(inormdataBin500bp1_counted1, "inormdataBin500bp1_counted1.filtmin.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

inormdataBin500bp1_counted1chr <- cSplit(inormdataBin500bp1_counted1, "Bin500bp","%")
head(inormdataBin500bp1_counted1chr)
inormdataBin500bp1_counted1chr <- inormdataBin500bp1_counted1chr[,c(10:12,1:9)]
write.table(inormdataBin500bp1_counted1chr, "inormdataBin500bp1_counted1filtminchr.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
sort -k1,1 -k2,2n inormdataBin500bp1_counted1filtminchr.txt > inormdataBin500bp1_counted1filtminchr.sort.txt
bedtools closest -a inormdataBin500bp1_counted1filtminchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > inormdataBin500bp1_counted1filtminchr_genehg38.txt

inormdataBin500bp1_counted2 <- inormdataBin500bp1_counted1[,1:8]
head(inormdataBin500bp1_counted2)
summary(inormdataBin500bp1_counted2)
#df <- as.inormdata.frame(inormdataBin500bp1)
inormdfBin500bp <- inormdataBin500bp1_counted2
head(inormdfBin500bp)
dim(inormdfBin500bp)
inormdfBin500bp = data.frame(inormdfBin500bp)
#write.table(inormdfBin500bp , "inormdfBin500bpdedup.txt", sep="\t", quote = FALSE, append = FALSE)
tinormdfBin500bp = t(inormdfBin500bp)
tinormdfBin500bp = data.frame(tinormdfBin500bp)
#write.table(tinormdfBin500bp , "tinormdfBin500bpdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfBin500bp)
tinormdfBin500bp["Color"] <-  c("corrClone","corrClone","corrClone","corrClone","Control","Case","Case","Control")
#head(tinormdfBin500bp)
#write.table(tinormdfBin500bp , "tinormdfBin500bpdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfBin500bp)
inormdfx <-tinormdfBin500bp[c(1:(length(tinormdfBin500bp)-1))]
PBin500bpC<-prcomp(inormdfx, center = TRUE, scale. = TRUE)
PBin500bpCi<-data.frame(PBin500bpC$x,Color=tinormdfBin500bp$Color)
percentagecomBin500bp <- round(PBin500bpC$sdev^2 / sum(PBin500bpC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomBin500bp <- paste( colnames(PBin500bpCi), "(", paste( as.character(percentagecomBin500bp), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pBin500bp1 <-ggplot(PBin500bpCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomBin500bp[1]) + ylab(percentagecomBin500bp[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("#FFB6B3","darkgrey","#BDE7BD"))+
  scale_shape_manual(values=c(0,1,2,0))
pBin500bp1 <- pBin500bp1 +theme_classic()
pBin500bp1 + xlim(-500,500)+ ylim(-500,500)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_tinormdfBin500bp_CntrlIndiv.svg", width=10*1.25, height=8*1.25, units="cm", dpi=96)


################################### Heatmap_indiv human Bin500bp ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(inormdataBin500bp1_counted1)
dim(inormdataBin500bp1_counted1)
sum(inormdataBin500bp1_counted1$freq)

#Heatmap_indiv with 2 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
inormdataBin500bp2 = as.matrix(inormdataBin500bp1_counted2)
head(inormdataBin500bp2)
dim(inormdataBin500bp2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
inormdataBin500bp3 <- inormdataBin500bp2
head(inormdataBin500bp3)
write.table(inormdataBin500bp3, "inormdataBin500bp3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

#Box and Violin Plot human Bin500bp using 12 Controls from MethHed
rownames(inormdataBin500bp)
colnames(inormdataBin500bp)
head(inormdataBin500bp)
dim(inormdataBin500bp)
Count_inormdataBin500bp1 <- count(inormdataBin500bp, "Bin500bp")
head(Count_inormdataBin500bp1)
inormagregateBin500bp1 = aggregate(inormdataBin500bp[,6:13],by=list(inormdataBin500bp$Bin500bp), mean)
head(inormagregateBin500bp1, 2)
inormdataBin500bp1=inormagregateBin500bp1
rownames(inormdataBin500bp1)
inormdataBin500bp1[,1]
rownames(inormdataBin500bp1)=inormdataBin500bp1[,1]
rownames(inormdataBin500bp1)
colnames(inormdataBin500bp1)
inormdataBin500bp1 = inormdataBin500bp1[,-1]
head(inormdataBin500bp1)
dim(inormdataBin500bp1)
write.table(inormdataBin500bp1, "inormaggregatedBin500bp_inormdataBin500bp1.txt", sep="\t", quote = FALSE, append = FALSE)
inormdataBin500bp1_counted <- cbind(inormdataBin500bp1, Count_inormdataBin500bp1)
head(inormdataBin500bp1_counted)
dim(inormdataBin500bp1_counted)
write.table(inormdataBin500bp1_counted, "inormaggregatedBin500bp_counted.txt", sep="\t", quote = FALSE, append = FALSE)

inormdataBin500bp1_counted1 <- inormdataBin500bp1_counted[which(inormdataBin500bp1_counted$freq >= 3),]
head(inormdataBin500bp1_counted1)
dim(inormdataBin500bp1_counted1)
inormdataBin500bp1_counted2 <- inormdataBin500bp1_counted1[,1:8]
head(inormdataBin500bp1_counted2)
VionormIndvBin500bp <- data.frame(inormdataBin500bp1_counted2[,c(5,8,6,7,1:4)])
head(VionormIndvBin500bp)
write.table(VionormIndvBin500bp, "VionormIndvBin500bp.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_VionormIndvBin500bp.svg", width=10, height=5, pointsize=12)
boxplot(VionormIndvBin500bp, main="Average methylation at human Bin500bps", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","navy","navy","navy","navy"))
dev.off()
dim(VionormIndvBin500bp)
VionormIndvBin500bp <- data.frame(VionormIndvBin500bp)
VionormIndvBin500bp1 <- VionormIndvBin500bp
head(VionormIndvBin500bp1,1)
VionormIndvBin500bp2 <- stack(VionormIndvBin500bp1)
head(VionormIndvBin500bp2)
colnames(VionormIndvBin500bp2) <- c("Methylation", "inormdatasets")
head(VionormIndvBin500bp2)
ggplot(VionormIndvBin500bp2, aes(x=inormdatasets, y=Methylation, color=inormdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","#FFB6B3","#FFB6B3","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","navy","navy","navy"))
ggsave("Violin_plot_VionormIndvBin500bp2_pimethhed_humanBin500bp.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)



#Heatmap of control normalized Log other samples merged
head(inormdataBin500bp1_counted2)
dim(inormdataBin500bp1_counted2)
inormdataBin500bp1_counted2ratio <- as.matrix(inormdataBin500bp1_counted2[,c(5,8,6,7,1:4)])

head(inormdataBin500bp1_counted2ratio)
dim(inormdataBin500bp1_counted2ratio)
winormdataBin500bp1_counted2ratioavg <- data.frame(cbind((rowMeans(inormdataBin500bp1_counted2ratio[,1:2])),
                                                         (inormdataBin500bp1_counted2ratio[,1:8])))

head(winormdataBin500bp1_counted2ratioavg)
colnames(winormdataBin500bp1_counted2ratioavg) <- c("Allcontrol",colnames(inormdataBin500bp1_counted2ratio))
head(winormdataBin500bp1_counted2ratioavg)
winormdataBin500bp1_counted2ratioavg1 <- winormdataBin500bp1_counted2ratioavg
winormdataBin500bp1_counted2ratioavg2 <- as.matrix(winormdataBin500bp1_counted2ratioavg1)
head(winormdataBin500bp1_counted2ratioavg2)
dim(winormdataBin500bp1_counted2ratioavg2)
winormdataBin500bp1_counted2ratioavg3 <- data.frame(cbind((winormdataBin500bp1_counted2ratioavg2[,1]/winormdataBin500bp1_counted2ratioavg2[,1]),
                                                          (winormdataBin500bp1_counted2ratioavg2[,2]/winormdataBin500bp1_counted2ratioavg2[,1]),
                                                          (winormdataBin500bp1_counted2ratioavg2[,3]/winormdataBin500bp1_counted2ratioavg2[,1]),
                                                          (winormdataBin500bp1_counted2ratioavg2[,4]/winormdataBin500bp1_counted2ratioavg2[,1]),
                                                          (winormdataBin500bp1_counted2ratioavg2[,5]/winormdataBin500bp1_counted2ratioavg2[,1]),
                                                          (winormdataBin500bp1_counted2ratioavg2[,6]/winormdataBin500bp1_counted2ratioavg2[,1]),
                                                          (winormdataBin500bp1_counted2ratioavg2[,7]/winormdataBin500bp1_counted2ratioavg2[,1]),
                                                          (winormdataBin500bp1_counted2ratioavg2[,8]/winormdataBin500bp1_counted2ratioavg2[,1]),
                                                          (winormdataBin500bp1_counted2ratioavg2[,9]/winormdataBin500bp1_counted2ratioavg2[,1])))


head(winormdataBin500bp1_counted2ratioavg3)
colnames(winormdataBin500bp1_counted2ratioavg3) <- colnames(winormdataBin500bp1_counted2ratioavg2)
head(winormdataBin500bp1_counted2ratioavg3)
winormdataBin500bp1_counted2ratioavg3 <- as.matrix(winormdataBin500bp1_counted2ratioavg3)
winormdataBin500bp1_counted2ratioavg3Log <- log2(winormdataBin500bp1_counted2ratioavg3)
head(winormdataBin500bp1_counted2ratioavg3Log)
write.table(winormdataBin500bp1_counted2ratioavg3Log , "Heatmap_ICF_normcontrolwinormdataBin500bp1_counted2ratioavg2Log.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)

winormdataBin500bp1_counted2ratioavg3Logminuscontrol <- winormdataBin500bp1_counted2ratioavg3Log[,c(2:9)]
head(winormdataBin500bp1_counted2ratioavg3Logminuscontrol)

minormdataBin500bp1_counted2ratioavg3 <- data.frame(cbind((winormdataBin500bp1_counted2ratioavg2[,1:9]-winormdataBin500bp1_counted2ratioavg2[,1])))


head(minormdataBin500bp1_counted2ratioavg3)
minormdataBin500bp1_counted2ratioavg3 <- as.matrix(minormdataBin500bp1_counted2ratioavg3)

minormdataBin500bp1_counted2ratioavg3minuscontrol <- minormdataBin500bp1_counted2ratioavg3[,c(2:9)]
head(minormdataBin500bp1_counted2ratioavg3minuscontrol)

#Differentially methylated probes specifc to cases
pinormdataBin500bp1_counted2ratioavg2 <- data.frame(winormdataBin500bp1_counted2ratioavg2)
head(pinormdataBin500bp1_counted2ratioavg2)
#Subtraction, Control -Control
pinormdataBin500bp1_counted2ratioavg2["D250_UN"] <- pinormdataBin500bp1_counted2ratioavg2$D250 - pinormdataBin500bp1_counted2ratioavg2$UN

#Subtraction, Case -Control
pinormdataBin500bp1_counted2ratioavg2["PG_Con"] <- pinormdataBin500bp1_counted2ratioavg2$PG - pinormdataBin500bp1_counted2ratioavg2$Allcontrol
pinormdataBin500bp1_counted2ratioavg2["PR_Con"] <- pinormdataBin500bp1_counted2ratioavg2$PR - pinormdataBin500bp1_counted2ratioavg2$Allcontrol
pinormdataBin500bp1_counted2ratioavg2["C7_Con"] <- pinormdataBin500bp1_counted2ratioavg2$C7 - pinormdataBin500bp1_counted2ratioavg2$Allcontrol
pinormdataBin500bp1_counted2ratioavg2["C13_Con"] <- pinormdataBin500bp1_counted2ratioavg2$C13 - pinormdataBin500bp1_counted2ratioavg2$Allcontrol
pinormdataBin500bp1_counted2ratioavg2["C35_Con"] <- pinormdataBin500bp1_counted2ratioavg2$C35 - pinormdataBin500bp1_counted2ratioavg2$Allcontrol
pinormdataBin500bp1_counted2ratioavg2["C50_Con"] <- pinormdataBin500bp1_counted2ratioavg2$C50 - pinormdataBin500bp1_counted2ratioavg2$Allcontrol
pinormdataBin500bp1_counted2ratioavg2["PG_D250"] <- pinormdataBin500bp1_counted2ratioavg2$PG - pinormdataBin500bp1_counted2ratioavg2$D250
pinormdataBin500bp1_counted2ratioavg2["PR_D250"] <- pinormdataBin500bp1_counted2ratioavg2$PR - pinormdataBin500bp1_counted2ratioavg2$D250
pinormdataBin500bp1_counted2ratioavg2["C7_D250"] <- pinormdataBin500bp1_counted2ratioavg2$C7 - pinormdataBin500bp1_counted2ratioavg2$D250
pinormdataBin500bp1_counted2ratioavg2["C13_D250"] <- pinormdataBin500bp1_counted2ratioavg2$C13 - pinormdataBin500bp1_counted2ratioavg2$D250
pinormdataBin500bp1_counted2ratioavg2["C35_D250"] <- pinormdataBin500bp1_counted2ratioavg2$C35 - pinormdataBin500bp1_counted2ratioavg2$D250
pinormdataBin500bp1_counted2ratioavg2["C50_D250"] <- pinormdataBin500bp1_counted2ratioavg2$C50 - pinormdataBin500bp1_counted2ratioavg2$D250
pinormdataBin500bp1_counted2ratioavg2["PG_UN"] <- pinormdataBin500bp1_counted2ratioavg2$PG - pinormdataBin500bp1_counted2ratioavg2$UN
pinormdataBin500bp1_counted2ratioavg2["PR_UN"] <- pinormdataBin500bp1_counted2ratioavg2$PR - pinormdataBin500bp1_counted2ratioavg2$UN
pinormdataBin500bp1_counted2ratioavg2["C7_UN"] <- pinormdataBin500bp1_counted2ratioavg2$C7 - pinormdataBin500bp1_counted2ratioavg2$UN
pinormdataBin500bp1_counted2ratioavg2["C13_UN"] <- pinormdataBin500bp1_counted2ratioavg2$C13 - pinormdataBin500bp1_counted2ratioavg2$UN
pinormdataBin500bp1_counted2ratioavg2["C35_UN"] <- pinormdataBin500bp1_counted2ratioavg2$C35 - pinormdataBin500bp1_counted2ratioavg2$UN
pinormdataBin500bp1_counted2ratioavg2["C50_UN"] <- pinormdataBin500bp1_counted2ratioavg2$C50 - pinormdataBin500bp1_counted2ratioavg2$UN
pinormdataBin500bp1_counted2ratioavg2["C7_PR"] <- pinormdataBin500bp1_counted2ratioavg2$C7 - pinormdataBin500bp1_counted2ratioavg2$PR
pinormdataBin500bp1_counted2ratioavg2["C13_PG"] <- pinormdataBin500bp1_counted2ratioavg2$C13 - pinormdataBin500bp1_counted2ratioavg2$PG
pinormdataBin500bp1_counted2ratioavg2["C35_PR"] <- pinormdataBin500bp1_counted2ratioavg2$C35 - pinormdataBin500bp1_counted2ratioavg2$PR
pinormdataBin500bp1_counted2ratioavg2["C50_PG"] <- pinormdataBin500bp1_counted2ratioavg2$C50 - pinormdataBin500bp1_counted2ratioavg2$PG

#Region where two controls are close
pinormdataBin500bp1_counted2ratioavg2 <- pinormdataBin500bp1_counted2ratioavg2[which(pinormdataBin500bp1_counted2ratioavg2$D250_UN < 0.1 & 
                                                                                       pinormdataBin500bp1_counted2ratioavg2$D250_UN > -0.1),]

dim(pinormdataBin500bp1_counted2ratioavg2)
#PG
pinormdataBin500bp1_counted2ratioavg2_PG <- pinormdataBin500bp1_counted2ratioavg2[which(pinormdataBin500bp1_counted2ratioavg2$PG_D250 > 0.2 | 
                                                                                          pinormdataBin500bp1_counted2ratioavg2$PG_D250 < -0.2),]


mpinormdataBin500bp1_counted2ratioavg2_PG <- pinormdataBin500bp1_counted2ratioavg2_PG[which(pinormdataBin500bp1_counted2ratioavg2_PG$PG_UN > 0.2 | 
                                                                                              pinormdataBin500bp1_counted2ratioavg2_PG$PG_UN < -0.2),]

head(mpinormdataBin500bp1_counted2ratioavg2_PG)
dim(mpinormdataBin500bp1_counted2ratioavg2_PG)

write.table(rownames(mpinormdataBin500bp1_counted2ratioavg2_PG) , "mpinormdataBin500bp1_counted2ratioavg2_PGchr.regions.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = F)

mpinormdataBin500bp1_counted2ratioavg2_PGchr <- data.frame(mpinormdataBin500bp1_counted2ratioavg2_PG)
mpinormdataBin500bp1_counted2ratioavg2_PGchr["pos"] <- rownames(mpinormdataBin500bp1_counted2ratioavg2_PGchr)
mpinormdataBin500bp1_counted2ratioavg2_PGchr <- cSplit(mpinormdataBin500bp1_counted2ratioavg2_PGchr, "pos", "%")
head(mpinormdataBin500bp1_counted2ratioavg2_PGchr)
dim(mpinormdataBin500bp1_counted2ratioavg2_PGchr)
mpinormdataBin500bp1_counted2ratioavg2_PGchr <- mpinormdataBin500bp1_counted2ratioavg2_PGchr[,c(33:35,1:32)]
write.table(mpinormdataBin500bp1_counted2ratioavg2_PGchr , "mpinormdataBin500bp1_counted2ratioavg2_PGchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = T)


#get beta-matrix on hypo and hyper DMR PG
fgrep -f mpinormdataBin500bp1_counted2ratioavg2_PGchr.regions.txt inormaggregatedBin500bp_counted.txt  -w > inormaggregatedBin500bp_countedBin500bp_PG.txt
inormaggregatedBin500bp_countedBin500bp_PG <- read.table("inormaggregatedBin500bp_countedBin500bp_PG.txt", header = F, stringsAsFactors = F, row.names = 1)
head(inormaggregatedBin500bp_countedBin500bp_PG)
dim(inormaggregatedBin500bp_countedBin500bp_PG)
colnames(inormaggregatedBin500bp_countedBin500bp_PG) <- colnames(inormdataBin500bp1_counted)
inormaggregatedBin500bp_countedBin500bp_PGchr <- data.frame(inormaggregatedBin500bp_countedBin500bp_PG)
inormaggregatedBin500bp_countedBin500bp_PGchr["pos"] <- rownames(inormaggregatedBin500bp_countedBin500bp_PGchr)
inormaggregatedBin500bp_countedBin500bp_PGchr <- cSplit(inormaggregatedBin500bp_countedBin500bp_PGchr, "pos", "%")
head(inormaggregatedBin500bp_countedBin500bp_PGchr)
dim(inormaggregatedBin500bp_countedBin500bp_PGchr)
inormaggregatedBin500bp_countedBin500bp_PGchr <- inormaggregatedBin500bp_countedBin500bp_PGchr[,c(11:13,9:10,1:8)]
write.table(inormaggregatedBin500bp_countedBin500bp_PGchr , "inormaggregatedBin500bp_countedBin500bp_PGchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = T)
rownames(inormaggregatedBin500bp_countedBin500bp_PGchr)  <- inormaggregatedBin500bp_countedBin500bp_PGchr$Bin500bp
pheatmap(inormaggregatedBin500bp_countedBin500bp_PGchr[,6:13])
sort -k1,1 -k2,2n inormaggregatedBin500bp_countedBin500bp_PGchr.txt | grep chr > inormaggregatedBin500bp_countedBin500bp_PGchr.sort.txt
sort -k1,1 -k2,2n mpinormdataBin500bp1_counted2ratioavg2_PGchr.txt | grep chr > mpinormdataBin500bp1_counted2ratioavg2_PGchr.sort.txt

Homo_sapiens.GRCh38.104.chr.gtf  <- rtracklayer::import("/home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.gtf")
head(Homo_sapiens.GRCh38.104.chr.gtf)
Homo_sapiens.GRCh38.104.chr.gtf =as.data.frame(Homo_sapiens.GRCh38.104.chr.gtf) 
Homo_sapiens.GRCh38.104.chr.gene.gtf <- Homo_sapiens.GRCh38.104.chr.gtf[which(Homo_sapiens.GRCh38.104.chr.gtf$type == 'gene'),]
head(Homo_sapiens.GRCh38.104.chr.gene.gtf)
dim(Homo_sapiens.GRCh38.104.chr.gene.gtf)
Homo_sapiens.GRCh38.104.chr.exon.gtf <- Homo_sapiens.GRCh38.104.chr.gtf[which(Homo_sapiens.GRCh38.104.chr.gtf$type == 'exon'),]
head(Homo_sapiens.GRCh38.104.chr.exon.gtf)
dim(Homo_sapiens.GRCh38.104.chr.exon.gtf)
Homo_sapiens.GRCh38.104.chr.five_prime_utr.gtf <- Homo_sapiens.GRCh38.104.chr.gtf[which(Homo_sapiens.GRCh38.104.chr.gtf$type == 'five_prime_utr'),]
head(Homo_sapiens.GRCh38.104.chr.five_prime_utr.gtf)
dim(Homo_sapiens.GRCh38.104.chr.five_prime_utr.gtf)
Homo_sapiens.GRCh38.104.chr.three_prime_utr.gtf <- Homo_sapiens.GRCh38.104.chr.gtf[which(Homo_sapiens.GRCh38.104.chr.gtf$type == 'three_prime_utr'),]
head(Homo_sapiens.GRCh38.104.chr.three_prime_utr.gtf)
dim(Homo_sapiens.GRCh38.104.chr.three_prime_utr.gtf)
Homo_sapiens.GRCh38.104.chr.CDS.gtf <- Homo_sapiens.GRCh38.104.chr.gtf[which(Homo_sapiens.GRCh38.104.chr.gtf$type == 'CDS'),]
head(Homo_sapiens.GRCh38.104.chr.CDS.gtf)
dim(Homo_sapiens.GRCh38.104.chr.CDS.gtf)

homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff <- read.delim("/home/ankitv/ref_av/gencodes/human/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff", header=F, stringsAsFactors = F)
head(homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff)
Homo_sapiens.GRCh38.104.chr.promoter.gtf <- homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff[which(homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff$V3 == 'promoter'),]
head(Homo_sapiens.GRCh38.104.chr.promoter.gtf)
dim(Homo_sapiens.GRCh38.104.chr.promoter.gtf)
Homo_sapiens.GRCh38.104.chr.promoter.gtf <- Homo_sapiens.GRCh38.104.chr.promoter.gtf[,c(1,4,5)]
write.table(Homo_sapiens.GRCh38.104.chr.promoter.gtf, "/home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.gtf", sep="\t", quote = F, append = F, row.names = F, col.names = F)
awk '{print "chr"$1"\t"$2"\t"$3}' Homo_sapiens.GRCh38.104.chr.promoter.gtf | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.104.chr.promoter.chr.gtf
#Promoter region overlapped with genes and closest to < 1 Kb
bedtools closest -wa -wb -a Homo_sapiens.GRCh38.104.chr.promoter.chr.gtf -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d | awk '{if($10 < 1000) print $0}' > Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf


awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.gtf | grep chrM -v | sort -k1,1 -k2,2n > /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf

bedtools closest -wa -wb -a mpinormdataBin500bp1_counted2ratioavg2_PGchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38.txt
bedtools closest -wa -wb -a inormaggregatedBin500bp_countedBin500bp_PGchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d  > inormaggregatedBin500bp_countedBin500bp_PGchr_genehg38.txt


bedtools closest -wa -wb -a mpinormdataBin500bp1_counted2ratioavg2_PGchr.sort.txt -b /home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf -d > mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38.txt
bedtools closest -wa -wb -a inormaggregatedBin500bp_countedBin500bp_PGchr.sort.txt -b /home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf -d  > inormaggregatedBin500bp_countedBin500bp_PGchr_promoterhg38.txt


#PR
pinormdataBin500bp1_counted2ratioavg2_PR <- pinormdataBin500bp1_counted2ratioavg2[which(pinormdataBin500bp1_counted2ratioavg2$PR_D250 > 0.2 | 
                                                                                          pinormdataBin500bp1_counted2ratioavg2$PR_D250 < -0.2),]


mpinormdataBin500bp1_counted2ratioavg2_PR <- pinormdataBin500bp1_counted2ratioavg2_PR[which(pinormdataBin500bp1_counted2ratioavg2_PR$PR_UN > 0.2 | 
                                                                                              pinormdataBin500bp1_counted2ratioavg2_PR$PR_UN < -0.2),]

head(mpinormdataBin500bp1_counted2ratioavg2_PR)
dim(mpinormdataBin500bp1_counted2ratioavg2_PR)

write.table(rownames(mpinormdataBin500bp1_counted2ratioavg2_PR) , "mpinormdataBin500bp1_counted2ratioavg2_PRchr.regions.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = F)

mpinormdataBin500bp1_counted2ratioavg2_PRchr <- data.frame(mpinormdataBin500bp1_counted2ratioavg2_PR)
mpinormdataBin500bp1_counted2ratioavg2_PRchr["pos"] <- rownames(mpinormdataBin500bp1_counted2ratioavg2_PRchr)
mpinormdataBin500bp1_counted2ratioavg2_PRchr <- cSplit(mpinormdataBin500bp1_counted2ratioavg2_PRchr, "pos", "%")
head(mpinormdataBin500bp1_counted2ratioavg2_PRchr)
dim(mpinormdataBin500bp1_counted2ratioavg2_PRchr)
mpinormdataBin500bp1_counted2ratioavg2_PRchr <- mpinormdataBin500bp1_counted2ratioavg2_PRchr[,c(33:35,1:32)]
write.table(mpinormdataBin500bp1_counted2ratioavg2_PRchr , "mpinormdataBin500bp1_counted2ratioavg2_PRchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = T)


#get beta-matrix on hypo and hyper DMR PR
fgrep -f mpinormdataBin500bp1_counted2ratioavg2_PRchr.regions.txt inormaggregatedBin500bp_counted.txt  -w > inormaggregatedBin500bp_countedBin500bp_PR.txt
inormaggregatedBin500bp_countedBin500bp_PR <- read.table("inormaggregatedBin500bp_countedBin500bp_PR.txt", header = F, stringsAsFactors = F, row.names = 1)
head(inormaggregatedBin500bp_countedBin500bp_PR)
colnames(inormaggregatedBin500bp_countedBin500bp_PR) <- colnames(inormdataBin500bp1_counted)
inormaggregatedBin500bp_countedBin500bp_PRchr <- data.frame(inormaggregatedBin500bp_countedBin500bp_PR)
inormaggregatedBin500bp_countedBin500bp_PRchr["pos"] <- rownames(inormaggregatedBin500bp_countedBin500bp_PRchr)
inormaggregatedBin500bp_countedBin500bp_PRchr <- cSplit(inormaggregatedBin500bp_countedBin500bp_PRchr, "pos", "%")
head(inormaggregatedBin500bp_countedBin500bp_PRchr)
dim(inormaggregatedBin500bp_countedBin500bp_PRchr)
inormaggregatedBin500bp_countedBin500bp_PRchr <- inormaggregatedBin500bp_countedBin500bp_PRchr[,c(11:13,9:10,1:8)]
write.table(inormaggregatedBin500bp_countedBin500bp_PRchr , "inormaggregatedBin500bp_countedBin500bp_PRchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = T)
rownames(inormaggregatedBin500bp_countedBin500bp_PRchr)  <- inormaggregatedBin500bp_countedBin500bp_PRchr$Bin500bp
pheatmap(inormaggregatedBin500bp_countedBin500bp_PRchr[,6:13])
sort -k1,1 -k2,2n inormaggregatedBin500bp_countedBin500bp_PRchr.txt | grep chr > inormaggregatedBin500bp_countedBin500bp_PRchr.sort.txt
sort -k1,1 -k2,2n mpinormdataBin500bp1_counted2ratioavg2_PRchr.txt | grep chr > mpinormdataBin500bp1_counted2ratioavg2_PRchr.sort.txt

bedtools closest -wa -wb -a mpinormdataBin500bp1_counted2ratioavg2_PRchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38.txt
bedtools closest -wa -wb -a inormaggregatedBin500bp_countedBin500bp_PRchr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d  > inormaggregatedBin500bp_countedBin500bp_PRchr_genehg38.txt


bedtools closest -wa -wb -a mpinormdataBin500bp1_counted2ratioavg2_PRchr.sort.txt -b /home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf -d > mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38.txt
bedtools closest -wa -wb -a inormaggregatedBin500bp_countedBin500bp_PRchr.sort.txt -b /home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf -d  > inormaggregatedBin500bp_countedBin500bp_PRchr_promoterhg38.txt


dim(miMERGE_myiNorm3re_pos_chravg_PG)
head(miMERGE_myiNorm3re_pos_chravg_PG)
dim(miMERGE_myiNorm3re_pos_chravg_PR)
head(miMERGE_myiNorm3re_pos_chravg_PR)
#Scatter Plot
head(myNorm2)
myNorm2ratio <- as.matrix(myNorm2[,c(5,8,6,7,1:4)])

head(myNorm2ratio)
dim(myNorm2ratio)
wmyNorm2ratioavg <- data.frame(cbind((rowMeans(myNorm2ratio[,1:2])),
                                     (myNorm2ratio[,1:8])))

head(wmyNorm2ratioavg)
colnames(wmyNorm2ratioavg) <- c("Allcontrol",colnames(myNorm2ratio))
head(wmyNorm2ratioavg)
wmyNorm2ratioavg <- as.matrix(wmyNorm2ratioavg)
head(wmyNorm2ratioavg)
dim(wmyNorm2ratioavg)
wmyNorm2ratioavg <- data.frame(wmyNorm2ratioavg)

head(inormdata1)
dim(inormdata1)
inormdata2 <- inormdata1[,6:14]
head(inormdata2)
rownames(inormdata2) <- inormdata2$TargetID
inormdata2 <- inormdata2[,-1]
inormdata2ratio <- as.matrix(inormdata2[,c(5,8,6,7,1:4)])

head(inormdata2ratio)
dim(inormdata2ratio)
winormdata2ratioavg <- data.frame(cbind((rowMeans(inormdata2ratio[,1:2])),
                                        (inormdata2ratio[,1:8])))

head(winormdata2ratioavg)
colnames(winormdata2ratioavg) <- c("Allcontrol",colnames(inormdata2ratio))
head(winormdata2ratioavg)
winormdata2ratioavg <- as.matrix(winormdata2ratioavg)
head(winormdata2ratioavg)
dim(winormdata2ratioavg)
winormdata2ratioavg <- data.frame(winormdata2ratioavg)

with(wmyNorm2ratioavg, plot(Allcontrol, PG, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#FFB6B3", bg="#FFB6B3",ylab="pG (Beta-value)", xlab="Control (Beta-value)", bty = 'n'))
with(subset(wmyNorm2ratioavg), points(Allcontrol, C13, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') #reddish
with(subset(wmyNorm2ratioavg), points(Allcontrol, C50, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') #reddish


with(subset(winormdata2ratioavg), points(Allcontrol, PG, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="magenta", bg="magenta"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdata2ratioavg), points(Allcontrol, C13, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdata2ratioavg), points(Allcontrol, C50, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish

lines(x = c(0,1), y = c(0,1), col = "#464647", lty = 1,lwd=3)
#save as wmyNorm2ratioavg_PG.png


with(wmyNorm2ratioavg, plot(Allcontrol, PR, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#FFB6B3", bg="#FFB6B3",ylab="PR (Beta-value)", xlab="Control (Beta-value)", bty = 'n'))
with(subset(wmyNorm2ratioavg), points(Allcontrol, C7, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') #reddish
with(subset(wmyNorm2ratioavg), points(Allcontrol, C35, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') #reddish

with(subset(winormdata2ratioavg), points(Allcontrol, PR, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="magenta", bg="magenta"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdata2ratioavg), points(Allcontrol, C7, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdata2ratioavg), points(Allcontrol, C35, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish

lines(x = c(0,1), y = c(0,1), col = "#464647", lty = 1,lwd=3)
#save as wmyNorm2ratioavg_PR.png

#Scatter with only diff meth specific probes
#Also add diff meth case specific probes
dim(miMERGE_myiNorm3re_chrpos_avg_PG)

with(miMERGE_myiNorm3re_chrpos_avg_PG, plot(Allcontrol, PG, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#FFB6B3", bg="#FFB6B3",ylab="pG (Beta-value)", xlab="Control (Beta-value)", bty = 'n'))
with(subset(miMERGE_myiNorm3re_chrpos_avg_PG), points(Allcontrol, C13, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') 
with(subset(miMERGE_myiNorm3re_chrpos_avg_PG), points(Allcontrol, C50, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') 


with(subset(winormdata2ratioavg), points(Allcontrol, PG, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="magenta", bg="magenta"),ylab="", xlab="", bty = 'n') 
with(subset(winormdata2ratioavg), points(Allcontrol, C13, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') 
with(subset(winormdata2ratioavg), points(Allcontrol, C50, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') 

lines(x = c(0,1), y = c(0,1), col = "#464647", lty = 1,lwd=3)
#save as miMERGE_myiNorm3re_chrpos_avg_PG_PG.png

dim(miMERGE_myiNorm3re_chrpos_avg_PR)

with(miMERGE_myiNorm3re_chrpos_avg_PR, plot(Allcontrol, PR, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#FFB6B3", bg="#FFB6B3",ylab="PR (Beta-value)", xlab="Control (Beta-value)", bty = 'n'))
with(subset(miMERGE_myiNorm3re_chrpos_avg_PR), points(Allcontrol, C7, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') 
with(subset(miMERGE_myiNorm3re_chrpos_avg_PR), points(Allcontrol, C35, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') 

with(subset(winormdata2ratioavg), points(Allcontrol, PR, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="magenta", bg="magenta"),ylab="", xlab="", bty = 'n') 
with(subset(winormdata2ratioavg), points(Allcontrol, C7, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') 
with(subset(winormdata2ratioavg), points(Allcontrol, C35, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') 

lines(x = c(0,1), y = c(0,1), col = "#464647", lty = 1,lwd=3)
#save as miMERGE_myiNorm3re_chrpos_avg_PR_PR.png

#Functional analysis_gene overlapped probes
mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38 <- read.table("mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38.txt", header =F)
#Gene distance <1000 bp
mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38 <- mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38[which(mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38$V42 < 1000),]
head(mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38)
dim(mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38)
mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38.symbol <- unique(mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38$V41)
write.table(data.frame(mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38.symbol), "mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38.symbol", quote=F, append=F, row.names = F, col.names = F)
mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38 <- read.table("mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38.txt", header =F)
mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38 <- mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38[which(mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38$V42 < 1000),]
head(mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38)
dim(mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38)
mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38.symbol <- unique(mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38$V41)
write.table(data.frame(mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38.symbol), "mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38.symbol", quote=F, append=F, row.names = F, col.names = F)

library(gprofiler2)
gost.diff_meth_500bp.PG <- gost(query=as.character(mpinormdataBin500bp1_counted2ratioavg2_PGchr_genehg38.symbol), organism="hsapiens")
gostplot(gost.diff_meth_500bp.PG)
gost.diff_meth_500bp.PR <- gost(query=as.character(mpinormdataBin500bp1_counted2ratioavg2_PRchr_genehg38.symbol), organism="hsapiens")
gostplot(gost.diff_meth_500bp.PR)
#Get excel
gost.diff_meth_500bp.PG.res.sorted <- gost.diff_meth_500bp.PG$result[order(gost.diff_meth_500bp.PG$result$p_value),]
head(gost.diff_meth_500bp.PG.res.sorted)
write.table(as.matrix(gost.diff_meth_500bp.PG.res.sorted),"gost.diff_meth_500bp.PG.res.sorted.txt",quote = F, append = F, sep = "\t")
write_xlsx((gost.diff_meth_500bp.PG.res.sorted),"gost.diff_meth_500bp.PG.res.sorted.xlsx")

gost.diff_meth_500bp.PR.res.sorted <- gost.diff_meth_500bp.PR$result[order(gost.diff_meth_500bp.PR$result$p_value),]
head(gost.diff_meth_500bp.PR.res.sorted)
write.table(as.matrix(gost.diff_meth_500bp.PR.res.sorted),"gost.diff_meth_500bp.PR.res.sorted.txt",quote = F, append = F, sep = "\t")
write_xlsx((gost.diff_meth_500bp.PR.res.sorted),"gost.diff_meth_500bp.PR.res.sorted.xlsx")

############----------------------
#Note: p_value output of gost function is Hypergeometric p-value after correction for multiple testing, see https://biit.cs.ut.ee/gprofiler_beta/page/apis#gost_query_results
gost.diff_meth_500bp.PG.res.sorted <- gost.diff_meth_500bp.PG$result[order(gost.diff_meth_500bp.PG$result$p_value),]
head(gost.diff_meth_500bp.PG.res.sorted)
#gost.diff_meth_500bp.PG.res.sorted <- gost.diff_meth_500bp.PG$result[order(gost.diff_meth_500bp.PG$result$p_value),]
head(gost.diff_meth_500bp.PG.res.sorted)
#write.table(as.matrix(gost.diff_meth_500bp.PG.res.sorted),"gost.diff_meth_500bp.PG.res.sorted.txt",quote = F, append = F, sep = "\t")
#GO:BP
gost.diff_meth_500bp.PG.res.sorted_gobp <- gost.diff_meth_500bp.PG.res.sorted[which(gost.diff_meth_500bp.PG.res.sorted$source == "GO:BP"),]
head(gost.diff_meth_500bp.PG.res.sorted_gobp)
dim(gost.diff_meth_500bp.PG.res.sorted_gobp)
gost.diff_meth_500bp.PG.res.sorted_gobp_bar <- gost.diff_meth_500bp.PG.res.sorted_gobp[,c(11,3)]
head(gost.diff_meth_500bp.PG.res.sorted_gobp_bar)
gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top <- head(gost.diff_meth_500bp.PG.res.sorted_gobp_bar,10)
gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top$term_name)
gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top$p_value)
gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top <- data.frame(gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#FFB6B3", 
          fill = "#FFB6B3" ,
          sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          sort.val = "asc",
          rotate = TRUE, 
          position = position_dodge(),
          ggtheme = theme_bw())

ggsave("gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top.svg", width=18, height=8, units="cm", dpi=96)
ggsave("gost.diff_meth_500bp.PG.res.sorted_gobp_bar_top.png", width=18, height=8, units="cm", dpi=96)


#Sort by p-value
gost.diff_meth_500bp.PR.res.sorted <- gost.diff_meth_500bp.PR$result[order(gost.diff_meth_500bp.PR$result$p_value),]
head(gost.diff_meth_500bp.PR.res.sorted)
#write.table(as.matrix(gost.diff_meth_500bp.PR.res.sorted),"gost.diff_meth_500bp.PR.res.sorted.txt",quote = F, append = F, sep = "\t")

#gost.diff_meth_500bp.PR.res.sorted <- gost.diff_meth_500bp.PR$result[order(gost.diff_meth_500bp.PR$result$p_value),]
#head(gost.diff_meth_500bp.PR.res.sorted)
#gost.diff_meth_500bp.PR.res.sorted <- gost.diff_meth_500bp.PR$result[order(gost.diff_meth_500bp.PR$result$p_value),]
#head(gost.diff_meth_500bp.PR.res.sorted)

#GO:BP
gost.diff_meth_500bp.PR.res.sorted_gobp <- gost.diff_meth_500bp.PR.res.sorted[which(gost.diff_meth_500bp.PR.res.sorted$source == "GO:BP"),]
head(gost.diff_meth_500bp.PR.res.sorted_gobp)
dim(gost.diff_meth_500bp.PR.res.sorted_gobp)
gost.diff_meth_500bp.PR.res.sorted_gobp_bar <- gost.diff_meth_500bp.PR.res.sorted_gobp[,c(11,3)]
head(gost.diff_meth_500bp.PR.res.sorted_gobp_bar)
gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top <- head(gost.diff_meth_500bp.PR.res.sorted_gobp_bar,10)
gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top$term_name)
gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top$p_value)
gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top <- data.frame(gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#FFB6B3", 
          fill = "#FFB6B3" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top.svg", width=18, height=8, units="cm", dpi=96)
ggsave("gost.diff_meth_500bp.PR.res.sorted_gobp_bar_top.png", width=18, height=8, units="cm", dpi=96)

bedtools closest -wa -wb -a mpinormdataBin500bp1_counted2ratioavg2_PGchr.sort.txt -b /home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf -d > mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38.txt
bedtools closest -wa -wb -a inormaggregatedBin500bp_countedBin500bp_PGchr.sort.txt -b /home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf -d  > inormaggregatedBin500bp_countedBin500bp_PGchr_promoterhg38.txt

bedtools closest -wa -wb -a mpinormdataBin500bp1_counted2ratioavg2_PRchr.sort.txt -b /home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf -d > mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38.txt
bedtools closest -wa -wb -a inormaggregatedBin500bp_countedBin500bp_PRchr.sort.txt -b /home/ankitv/ref_av/gencodes/human/Homo_sapiens.GRCh38.104.chr.promoter.chr.closestgene.gtf -d  > inormaggregatedBin500bp_countedBin500bp_PRchr_promoterhg38.txt


#Functional analysis_gene overlapped probes
mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38 <- read.table("mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38.txt", header =F)
head(mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38,1)
#Gene distance <1000 bp
mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38 <- mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38[which(mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38$V46 < 1000),]
head(mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38)
dim(mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38)
mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38.symbol <- unique(mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38$V44)
write.table(data.frame(mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38.symbol), "mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38.symbol", quote=F, append=F, row.names = F, col.names = F)
mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38 <- read.table("mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38.txt", header =F)
mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38 <- mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38[which(mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38$V46 < 1000),]
head(mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38)
dim(mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38)
mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38.symbol <- unique(mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38$V44)
write.table(data.frame(mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38.symbol), "mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38.symbol", quote=F, append=F, row.names = F, col.names = F)

library(gprofiler2)
gost.diff_meth_promoter_500bp.PG <- gost(query=as.character(mpinormdataBin500bp1_counted2ratioavg2_PGchr_promoterhg38.symbol), organism="hsapiens")
gostplot(gost.diff_meth_promoter_500bp.PG)
gost.diff_meth_promoter_500bp.PR <- gost(query=as.character(mpinormdataBin500bp1_counted2ratioavg2_PRchr_promoterhg38.symbol), organism="hsapiens")
gostplot(gost.diff_meth_promoter_500bp.PR)
#Get excel
gost.diff_meth_promoter_500bp.PG.res.sorted <- gost.diff_meth_promoter_500bp.PG$result[order(gost.diff_meth_promoter_500bp.PG$result$p_value),]
head(gost.diff_meth_promoter_500bp.PG.res.sorted)
write.table(as.matrix(gost.diff_meth_promoter_500bp.PG.res.sorted),"gost.diff_meth_promoter_500bp.PG.res.sorted.txt",quote = F, append = F, sep = "\t")
write_xlsx((gost.diff_meth_promoter_500bp.PG.res.sorted),"gost.diff_meth_promoter_500bp.PG.res.sorted.xlsx")

gost.diff_meth_promoter_500bp.PR.res.sorted <- gost.diff_meth_promoter_500bp.PR$result[order(gost.diff_meth_promoter_500bp.PR$result$p_value),]
head(gost.diff_meth_promoter_500bp.PR.res.sorted)
write.table(as.matrix(gost.diff_meth_promoter_500bp.PR.res.sorted),"gost.diff_meth_promoter_500bp.PR.res.sorted.txt",quote = F, append = F, sep = "\t")
write_xlsx((gost.diff_meth_promoter_500bp.PR.res.sorted),"gost.diff_meth_promoter_500bp.PR.res.sorted.xlsx")

############----------------------
#Note: p_value output of gost function is Hypergeometric p-value after correction for multiple testing, see https://biit.cs.ut.ee/gprofiler_beta/page/apis#gost_query_results
gost.diff_meth_promoter_500bp.PG.res.sorted <- gost.diff_meth_promoter_500bp.PG$result[order(gost.diff_meth_promoter_500bp.PG$result$p_value),]
head(gost.diff_meth_promoter_500bp.PG.res.sorted)
#gost.diff_meth_promoter_500bp.PG.res.sorted <- gost.diff_meth_promoter_500bp.PG$result[order(gost.diff_meth_promoter_500bp.PG$result$p_value),]
head(gost.diff_meth_promoter_500bp.PG.res.sorted)
#write.table(as.matrix(gost.diff_meth_promoter_500bp.PG.res.sorted),"gost.diff_meth_promoter_500bp.PG.res.sorted.txt",quote = F, append = F, sep = "\t")
#GO:BP
gost.diff_meth_promoter_500bp.PG.res.sorted_gobp <- gost.diff_meth_promoter_500bp.PG.res.sorted[which(gost.diff_meth_promoter_500bp.PG.res.sorted$source == "GO:BP"),]
head(gost.diff_meth_promoter_500bp.PG.res.sorted_gobp)
dim(gost.diff_meth_promoter_500bp.PG.res.sorted_gobp)
gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar <- gost.diff_meth_promoter_500bp.PG.res.sorted_gobp[,c(11,3)]
head(gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar)
gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top <- head(gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar,10)
gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top$term_name)
gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top$p_value)
gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top <- data.frame(gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#FFB6B3", 
          fill = "#FFB6B3" ,
          sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          sort.val = "asc",
          rotate = TRUE, 
          position = position_dodge(),
          ggtheme = theme_bw())

ggsave("gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top.svg", width=18, height=8, units="cm", dpi=96)
ggsave("gost.diff_meth_promoter_500bp.PG.res.sorted_gobp_bar_top.png", width=18, height=8, units="cm", dpi=96)


#Sort by p-value
gost.diff_meth_promoter_500bp.PR.res.sorted <- gost.diff_meth_promoter_500bp.PR$result[order(gost.diff_meth_promoter_500bp.PR$result$p_value),]
head(gost.diff_meth_promoter_500bp.PR.res.sorted)
#write.table(as.matrix(gost.diff_meth_promoter_500bp.PR.res.sorted),"gost.diff_meth_promoter_500bp.PR.res.sorted.txt",quote = F, append = F, sep = "\t")

#gost.diff_meth_promoter_500bp.PR.res.sorted <- gost.diff_meth_promoter_500bp.PR$result[order(gost.diff_meth_promoter_500bp.PR$result$p_value),]
#head(gost.diff_meth_promoter_500bp.PR.res.sorted)
#gost.diff_meth_promoter_500bp.PR.res.sorted <- gost.diff_meth_promoter_500bp.PR$result[order(gost.diff_meth_promoter_500bp.PR$result$p_value),]
#head(gost.diff_meth_promoter_500bp.PR.res.sorted)

#GO:BP
gost.diff_meth_promoter_500bp.PR.res.sorted_gobp <- gost.diff_meth_promoter_500bp.PR.res.sorted[which(gost.diff_meth_promoter_500bp.PR.res.sorted$source == "GO:BP"),]
head(gost.diff_meth_promoter_500bp.PR.res.sorted_gobp)
dim(gost.diff_meth_promoter_500bp.PR.res.sorted_gobp)
gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar <- gost.diff_meth_promoter_500bp.PR.res.sorted_gobp[,c(11,3)]
head(gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar)
gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top <- head(gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar,10)
gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top$term_name)
gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top$p_value)
gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top <- data.frame(gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#FFB6B3", 
          fill = "#FFB6B3" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top.svg", width=18, height=8, units="cm", dpi=96)
ggsave("gost.diff_meth_promoter_500bp.PR.res.sorted_gobp_bar_top.png", width=18, height=8, units="cm", dpi=96)


#ComplexHeatmap
#histone_counts_at_1kb_bins.txt provided by Varsha
library(ComplexHeatmap)
library(circlize)
mat <- read.table("histone_counts_at_1kb_bins.txt", header = T, row.names = 1)
depth_readcount<-c(48236923,42403316,42578879,24398414,19214480,58102285,46219210,39960630,22816468,23150450,53584323,43762985,37122298,22264590,30049140,42170864,42986895,39732568,26010420,
                   32343816,51699787,39567924,45383552,25073548,25427094,57820998,38224806,44385864,24406289,24037939,53306867,41599398,44406006,41599398,44406006)

head(mat,4)
mat <- as.matrix(mat)
matx <- cbind.data.frame(rownames(mat),
                         (mat[,1]/depth_readcount[1]) * 1000000,
                         (mat[,2]/depth_readcount[2]) * 1000000,
                         (mat[,3]/depth_readcount[3]) * 1000000,
                         (mat[,4]/depth_readcount[4]) * 1000000,
                         (mat[,5]/depth_readcount[5]) * 1000000,
                         (mat[,6]/depth_readcount[6]) * 1000000,
                         (mat[,7]/depth_readcount[7]) * 1000000,
                         (mat[,8]/depth_readcount[8]) * 1000000,
                         (mat[,9]/depth_readcount[9]) * 1000000,
                         (mat[,10]/depth_readcount[10]) * 1000000,
                         (mat[,11]/depth_readcount[11]) * 1000000,
                         (mat[,12]/depth_readcount[12]) * 1000000,
                         (mat[,13]/depth_readcount[13]) * 1000000,
                         (mat[,14]/depth_readcount[14]) * 1000000,
                         (mat[,15]/depth_readcount[15]) * 1000000,
                         (mat[,16]/depth_readcount[16]) * 1000000,
                         (mat[,17]/depth_readcount[17]) * 1000000,
                         (mat[,18]/depth_readcount[18]) * 1000000,
                         (mat[,19]/depth_readcount[19]) * 1000000,
                         (mat[,20]/depth_readcount[20]) * 1000000,
                         (mat[,21]/depth_readcount[21]) * 1000000,
                         (mat[,22]/depth_readcount[22]) * 1000000,
                         (mat[,23]/depth_readcount[23]) * 1000000,
                         (mat[,24]/depth_readcount[24]) * 1000000,
                         (mat[,25]/depth_readcount[25]) * 1000000,
                         (mat[,26]/depth_readcount[26]) * 1000000,
                         (mat[,27]/depth_readcount[27]) * 1000000,
                         (mat[,28]/depth_readcount[28]) * 1000000,
                         (mat[,29]/depth_readcount[29]) * 1000000,
                         (mat[,30]/depth_readcount[30]) * 1000000,
                         (mat[,31]/depth_readcount[31]) * 1000000,
                         (mat[,32]/depth_readcount[32]) * 1000000,
                         (mat[,33]/depth_readcount[33]) * 1000000,
                         (mat[,34]/depth_readcount[34]) * 1000000,
                         (mat[,35]/depth_readcount[35]) * 1000000)




head(matx)
colnames(matx) <- c("row",colnames(mat))


matxNorm <- cbind.data.frame(rownames(matx),
                             (matx$WT1_K36_1-matx$WT1_in),
                             (matx$WT1_K36_2-matx$WT1_in),
                             (matx$WT1_K4_1-matx$WT1_in),
                             (matx$WT1_K4_2-matx$WT1_in),
                             (matx$pR_K36_1-matx$pR_in),
                             (matx$pR_K36_2-matx$pR_in),
                             (matx$pR_K4_1-matx$pR_in),
                             (matx$pR_K4_2-matx$pR_in),
                             (matx$c7_K36_1-matx$c7_in),
                             (matx$c7_K36_2-matx$c7_in),
                             (matx$c7_K4_1-matx$c7_in),
                             (matx$c7_K4_2-matx$c7_in),
                             (matx$c35_K36_1-matx$c35_in),
                             (matx$c35_K36_2-matx$c35_in),
                             (matx$c35_K4_1-matx$c35_in),
                             (matx$c35_K4_2-matx$c35_in),
                             (matx$pG_K36_1-matx$pG_in),
                             (matx$pG_K36_2-matx$pG_in),
                             (matx$pG_K4_1-matx$pG_in),
                             (matx$pG_K4_2-matx$pG_in),
                             (matx$c13_K36_1-matx$c13_in),
                             (matx$c13_K36_2-matx$c13_in),
                             (matx$c13_K4_1-matx$c13_in),
                             (matx$c13_K4_2-matx$c13_in),
                             (matx$c50_K36_1-matx$c50_in),
                             (matx$c50_K36_2-matx$c50_in),
                             (matx$c50_K4_1-matx$c50_in),
                             (matx$c50_K4_2-matx$c50_in))
head(matxNorm)
dim(matxNorm)

colnames(matxNorm) <- c("row","WT1_K36_1","WT1_K36_2","WT1_K4_1","WT1_K4_2","pR_K36_1","pR_K36_2","pR_K4_1","pR_K4_2","c7_K36_1","c7_K36_2","c7_K4_1","c7_K4_2","c35_K36_1","c35_K36_2","c35_K4_1","c35_K4_2","pG_K36_1","pG_K36_2","pG_K4_1","pG_K4_2","c13_K36_1","c13_K36_2","c13_K4_1","c13_K4_2","c50_K36_1","c50_K36_2","c50_K4_1","c50_K4_2")
rownames(matxNorm) <- matxNorm$row

inormaggregatedBin1Kb_countedBin1kbchr_hg38 <- read.table("inormaggregatedBin1Kb_countedBin1kbchr.txt", header = T)
dim(inormaggregatedBin1Kb_countedBin1kbchr_hg38)

minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38 <- read.table("minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr.txt", header = T)
head(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38)

dim(matxNorm)
head(inormaggregatedBin1Kb_countedBin1kbchr_hg38)
head(matxNorm)
minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx <- cbind.data.frame(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38, matxNorm)
head(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx)
dim(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx)
matxNorm1 <- as.matrix(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx[,c(seq(4,11,by =1))])
matxNorm2 <- as.matrix(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx[,c(seq(13,40,by =4))])
matxNorm3 <- as.matrix(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx[,c(seq(14,40,by =4))])
matxNorm4 <- as.matrix(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx[,c(seq(15,40,by =4))])
matxNorm5 <- as.matrix(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx[,c(seq(16,40,by =4))])

library(circlize)
array_col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht1 = Heatmap(matxNorm1, name = "Array_Beta_dmbins", 
              column_order =  colnames(matxNorm1), 
              show_row_names = T,
              col = array_col_fun,
              column_title = "Genomic DMR: Array",
              row_names_gp = gpar(fontsize=6),
              column_labels = colnames(matxNorm1),
              row_names_side = "left",
              column_names_gp = gpar(fontsize=12),
              use_raster = T, raster_quality = 8)
K36_col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "orange"))

ht2 = Heatmap(matxNorm2, name = "K36_1_dmbins", 
              column_order =  colnames(matxNorm2), 
              show_row_names = T,
              col = K36_col_fun,
              column_title = "K36_1",
              row_names_gp = gpar(fontsize=6),
              column_labels = colnames(matxNorm2),
              row_names_side = "left",
              column_names_gp = gpar(fontsize=12),
              use_raster = T, raster_quality = 8)

K36_col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "orange"))

ht3 = Heatmap(matxNorm3, name = "K36_2_dmbins", 
              column_order =  colnames(matxNorm3), 
              show_row_names = T,
              col = K36_col_fun,
              column_title = "K36_2",
              row_names_gp = gpar(fontsize=6),
              column_labels = colnames(matxNorm3),
              row_names_side = "left",
              column_names_gp = gpar(fontsize=12),
              use_raster = T, raster_quality = 8)


K4_col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "orange"))

ht4 = Heatmap(matxNorm4, name = "K4_1_dmbins", 
              column_order =  colnames(matxNorm4), 
              show_row_names = T,
              col = K4_col_fun,
              column_title = "K4_1",
              row_names_gp = gpar(fontsize=6),
              column_labels = colnames(matxNorm4),
              row_names_side = "left",
              column_names_gp = gpar(fontsize=12),
              use_raster = T, raster_quality = 8)

K4_col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "orange"))

ht5 = Heatmap(matxNorm5, name = "K4_2_dmbins", 
              column_order =  colnames(matxNorm5), 
              show_row_names = T,
              col = K4_col_fun,
              column_title = "K4_2",
              row_names_gp = gpar(fontsize=6),
              column_labels = colnames(matxNorm5),
              row_names_side = "left",
              column_names_gp = gpar(fontsize=12),
              use_raster = T, raster_quality = 8)



ht1+ ht2+ht3+ ht4+ ht5



head(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx)
dim(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx)
rownames(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx) <- paste0(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx$pos_1, "_",minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx$pos_2, "_", minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx$pos_3)
minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1 <- as.matrix(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx[,c(4:11,13:40)])
minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1 <- minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1[,c(1:8,9,13,17,21,25,29,33,10,14,18,22,26,30,34,11,15,19,23,27,31,35,12,16,20,24,28,32,36)]
dim(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1)
head(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1)
miBin1Kb1_c2rdiff3chr_hg38_matx1 <- data.frame(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1)
ht6 = Heatmap(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1, name = "dmbins", 
              column_order =  colnames(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1), 
              show_row_names = T,
              col = array_col_fun,
              column_title = "Array with Histone",
              row_names_gp = gpar(fontsize=6),
              column_labels = colnames(minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1),
              row_names_side = "left",
              column_names_gp = gpar(fontsize=12),
              use_raster = F, 
              raster_quality = 20,
              cluster_columns = TRUE, 
              cluster_rows = TRUE, 
              show_row_dend = TRUE,
              row_km = 2, 
              column_split  = c(rep("A",8),rep("B",7),rep("C",7),rep("D",7),rep("E",7) ))

draw(ht6)
#minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1.pdf
#minormdataBin1Kb1_counted2ratiodiff3minuscontrolchr_hg38_matx1.svg
Monk_human_ICR_hg38 <- read.table("/home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed", header=F, stringsAsFactors = F)
colnames(Monk_human_ICR_hg38) <- c("chr","start","end","DMR", "DMRtype")
minormdataICR1_counted2ratioavg3 <- read.table("minormdataICR1_counted2ratioavg3.txt", header = T, stringsAsFactors = F)
minormdataICR1_counted2ratioavg4 <- minormdataICR1_counted2ratioavg3
minormdataICR1_counted2ratioavg4["DMR"] <- rownames(minormdataICR1_counted2ratioavg4)
minormdataICR1_counted2ratioavgDMR <- merge(minormdataICR1_counted2ratioavg4, Monk_human_ICR_hg38, by="DMR")
minormdataICR1_counted2ratioavgDMR["row"] <- minormdataICR1_counted2ratioavgDMR$DMRtype
Monk_ICR_dnmt3b <- read.table("Monk_ICR_dnmt3b_IP_counts.txt", header = T, stringsAsFactors = F)
Monk_ICR_dnmt3bx <- Monk_ICR_dnmt3b[,c(22:28)]
colnames(Monk_ICR_dnmt3bx) <- paste0(c(colnames(Monk_ICR_dnmt3bx)),"_DNMT3b")
Monk_ICR_dnmt3bx["row"] <- rownames(Monk_ICR_dnmt3bx)
Monk_ICR_h3k36me3 <- read.table("/home/ankitv/Downloads/Monk_ICR_h3k36me3_IP_counts.txt")
Monk_ICR_h3k36me3x <- Monk_ICR_h3k36me3[,c(22:28)]
colnames(Monk_ICR_h3k36me3x) <- paste0(c(colnames(Monk_ICR_h3k36me3x)),"_K36")
Monk_ICR_h3k36me3x["row"] <- rownames(Monk_ICR_h3k36me3x)
Monk_ICR_h3k4me3 <- read.table("/home/ankitv/Downloads/Monk_ICR_h3k4me3_IP_counts.txt")
Monk_ICR_h3k4me3x <- Monk_ICR_h3k4me3[,c(22:28)]
colnames(Monk_ICR_h3k4me3x) <- paste0(c(colnames(Monk_ICR_h3k4me3x)),"_K4")
Monk_ICR_h3k4me3x["row"] <- rownames(Monk_ICR_h3k4me3x)


minormdataICR1_counted2ratioavgDMR_dnmt3b <- merge(minormdataICR1_counted2ratioavgDMR, Monk_ICR_dnmt3bx, by="row")
head(minormdataICR1_counted2ratioavgDMR_dnmt3b)
minormdataICR1_counted2ratioavgDMR_dnmt3b_k36 <- merge(minormdataICR1_counted2ratioavgDMR_dnmt3b, Monk_ICR_h3k36me3x, by="row")
minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4 <- merge(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36, Monk_ICR_h3k4me3x, by="row")
head(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4,1)
dim(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4)
rownames(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4)  <- minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4$row
head(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4,1)
write.table(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4,"minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4.txt",row.names=F,quote=FALSE,append = F,sep = "\t")

inormdataICR1_counted1 <- read.table("inormdataICR1_counted1.filtmin3.txt", header = T, stringsAsFactors = F)
inormdataICR1_counted4 <- inormdataICR1_counted1
inormdataICR1_counted4["DMR"] <- rownames(inormdataICR1_counted4)
inormdataICR1_countedDMR <- merge(inormdataICR1_counted4, Monk_human_ICR_hg38, by="DMR")
inormdataICR1_countedDMR["row"] <- inormdataICR1_countedDMR$DMRtype


inormdataICR1_countedDMR_dnmt3b <- merge(inormdataICR1_countedDMR, Monk_ICR_dnmt3bx, by="row")
head(inormdataICR1_countedDMR_dnmt3b)
inormdataICR1_countedDMR_dnmt3b_k36 <- merge(inormdataICR1_countedDMR_dnmt3b, Monk_ICR_h3k36me3x, by="row")
inormdataICR1_countedDMR_dnmt3b_k36_k4 <- merge(inormdataICR1_countedDMR_dnmt3b_k36, Monk_ICR_h3k4me3x, by="row")
head(inormdataICR1_countedDMR_dnmt3b_k36_k4,1)
dim(inormdataICR1_countedDMR_dnmt3b_k36_k4)
rownames(inormdataICR1_countedDMR_dnmt3b_k36_k4)  <- inormdataICR1_countedDMR_dnmt3b_k36_k4$row
head(inormdataICR1_countedDMR_dnmt3b_k36_k4,1)
write.table(inormdataICR1_countedDMR_dnmt3b_k36_k4,"inormdataICR1_countedDMR_dnmt3b_k36_k4.txt",row.names=F,quote=FALSE,append = F,sep = "\t")

library(ComplexHeatmap)
library(circlize)
minorm_ICR_only <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4[,c(4,5,7,8,10,6,9,11)])
inorm_ICR_only <- as.matrix(inormdataICR1_countedDMR_dnmt3b_k36_k4[,c(7,10,8,9,3:6)])
dnmt3b_only <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4[,c(seq(16,22,by =1))])
k36_only <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4[,c(seq(23,29,by =1))])
k4_only <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4[,c(seq(30,36,by =1))])

array_col_fun = colorRamp2(c(-0.5, 0, 0.5), c("navy", "white", "darkred"))
normarray_col_fun = colorRamp2(c(0, 0.5, 1), c("navy", "white", "darkred"))

minorm_ICRonlyheat = Heatmap(minorm_ICR_only, name = "Array_Beta", 
                             column_order =  colnames(minorm_ICR_only), 
                             show_row_names = T,
                             col = array_col_fun,
                             column_title = "Genomic DMR: Array",
                             row_names_gp = gpar(fontsize=6),
                             column_labels = colnames(minorm_ICR_only),
                             row_names_side = "left",
                             column_names_gp = gpar(fontsize=12),
                             use_raster = T, 
                             raster_quality = 8,
                             row_km = 3,
                             border =FALSE,
                             rect_gp = gpar(col = "black", lwd = 2))
DNMT3b_col_fun = colorRamp2(c(0, 1, 2), c("green", "white", "orange"))
#DNMT3b_col_fun = colorRamp2(c(0,  1), c("white", "orange"))
#DNMT3b_col_fun = colorRamp2(c(0, 1, 2), c("white","darkgreen", "orange"))

dnmt3b_onlyheat = Heatmap(dnmt3b_only, name = "DNMT3b_IP", 
                          column_order =  colnames(dnmt3b_only), 
                          show_row_names = T,
                          col = DNMT3b_col_fun,
                          column_title = "DNMT3b_IP",
                          row_names_gp = gpar(fontsize=6),
                          column_labels = colnames(dnmt3b_only),
                          row_names_side = "left",
                          column_names_gp = gpar(fontsize=12),
                          use_raster = T, 
                          raster_quality = 8,
                          row_km = 3,
                          border =TRUE,
                          rect_gp = gpar(col = "black", lwd = 2))

K36_col_fun = colorRamp2(c(0, 1, 2), c("darkgrey", "white", "#DEC20B"))
#K36_col_fun = colorRamp2(c(0, 1), c( "white", "darkgrey"))
#K36_col_fun = colorRamp2(c(0, 1, 2), c("white","darkgrey",  "#DEC20B"))

k36_onlyheat = Heatmap(k36_only, name = "H3K36me3", 
                       column_order =  colnames(k36_only), 
                       show_row_names = T,
                       col = K36_col_fun,
                       column_title = "K36_1",
                       row_names_gp = gpar(fontsize=6),
                       column_labels = colnames(k36_only),
                       row_names_side = "left",
                       column_names_gp = gpar(fontsize=12),
                       use_raster = T, 
                       raster_quality = 8,
                       row_km = 3,
                       border =TRUE,
                       rect_gp = gpar(col = "black", lwd = 2))


K4_col_fun = colorRamp2(c(0, 3, 6), c("#84BDA9", "white", "#BD9A7A"))
#K4_col_fun = colorRamp2(c(0, 3, 6), c("white", "#84BDA9", "#BD9A7A"))
#K4_col_fun = colorRamp2(c(0,  4), c("white", "#BD9A7A"))

k4_onlyheat = Heatmap(k4_only, name = "H3K4me3", 
                      column_order =  colnames(k4_only), 
                      show_row_names = T,
                      col = K4_col_fun,
                      column_title = "K4_1",
                      row_names_gp = gpar(fontsize=6),
                      column_labels = colnames(k4_only),
                      row_names_side = "left",
                      column_names_gp = gpar(fontsize=12),
                      use_raster = T, 
                      raster_quality = 8,
                      row_km = 3,
                      border =TRUE,
                      rect_gp = gpar(col = "black", lwd = 2))


minorm_ICRonlyheat + dnmt3b_onlyheat + k36_onlyheat + k4_onlyheat
#minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41_col.pdf
#minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41_col.svg

minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41 <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4[,c(3:11,16:36)])
dim(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41)
head(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41,1)
dim(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41)
dnmt3b_k36_k4_ICR6 = Heatmap(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41, name = "dmbins", 
                             column_order =  colnames(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41), 
                             show_row_names = T,
                             col = array_col_fun,
                             column_title = "Array with Histone",
                             row_names_gp = gpar(fontsize=6),
                             column_labels = colnames(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41),
                             row_names_side = "left",
                             column_names_gp = gpar(fontsize=12),
                             use_raster = F, 
                             raster_quality = 20,
                             cluster_columns = TRUE, 
                             cluster_rows = TRUE, 
                             show_row_dend = TRUE,
                             row_km = 2, 
                             column_split  = c(rep("A",9),rep("B",7),rep("C",7),rep("D",7)))

draw(dnmt3b_k36_k4_ICR6)
#minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41.pdf
#minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k41.svg


#Extract only subset ICRs
fgrep -f diff_imp_loci_pG_pR.txt minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4.txt > minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR.txt
minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR <- read.table("minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR.txt", header=F, stringsAsFactors = F)
colnames(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR) <- colnames(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4)
rownames(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR) <-  minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR$row
head(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR)
dim(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR)
fgrep -f diff_imp_loci_pG_pR.txt inormdataICR1_countedDMR_dnmt3b_k36_k4.txt > inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR.txt
inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR <- read.table("inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR.txt", header=F, stringsAsFactors = F)
colnames(inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR) <- colnames(inormdataICR1_countedDMR_dnmt3b_k36_k4)
rownames(inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR) <-  inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR$row
head(inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR)
dim(inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR)


library(ComplexHeatmap)
library(circlize)
minorm_ICR_subsetonly <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR[,c(4,5,7,8,10,6,9,11)])
inorm_ICR_subsetonly <- as.matrix(inormdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipGpR[,c(7,10,9,3,5,8,4,6)])
dnmt3b_subsetonly <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR[,c(seq(16,22,by =1))])
k36_subsetonly <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR[,c(seq(23,29,by =1))])
k4_subsetonly <- as.matrix(minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR[,c(seq(30,36,by =1))])


minorm_ICRsubsetonlyheat = Heatmap(minorm_ICR_subsetonly, name = "Array_Beta", 
                                   column_order =  colnames(minorm_ICR_subsetonly), 
                                   show_row_names = T,
                                   col = array_col_fun,
                                   column_title = "Genomic DMR: Array",
                                   row_names_gp = gpar(fontsize=6),
                                   column_labels = colnames(minorm_ICR_subsetonly),
                                   row_names_side = "left",
                                   column_names_gp = gpar(fontsize=12),
                                   use_raster = T, 
                                   raster_quality = 8, 
                                   row_km = 3,
                                   border = TRUE,
                                   rect_gp = gpar(col = "black", lwd = 2))

inorm_ICRsubsetonlyheat = Heatmap(inorm_ICR_subsetonly, name = "Array_Beta", 
                                  column_order =  colnames(inorm_ICR_subsetonly), 
                                  show_row_names = T,
                                  col = normarray_col_fun,
                                  column_title = "Genomic DMR: Array",
                                  row_names_gp = gpar(fontsize=6),
                                  column_labels = colnames(inorm_ICR_subsetonly),
                                  row_names_side = "left",
                                  column_names_gp = gpar(fontsize=12),
                                  use_raster = T, 
                                  raster_quality = 8, 
                                  row_km = 3,
                                  border = TRUE,
                                  rect_gp = gpar(col = "black", lwd = 2))

dnmt3b_subsetonlyheat = Heatmap(dnmt3b_subsetonly, name = "DNMT3b_IP", 
                                column_order =  colnames(dnmt3b_subsetonly), 
                                show_row_names = T,
                                col = DNMT3b_col_fun,
                                column_title = "DNMT3b_IP",
                                row_names_gp = gpar(fontsize=6),
                                column_labels = colnames(dnmt3b_subsetonly),
                                row_names_side = "left",
                                column_names_gp = gpar(fontsize=12),
                                use_raster = T, 
                                raster_quality = 8, 
                                row_km = 3,
                                border = TRUE,
                                rect_gp = gpar(col = "black", lwd = 2))


k36_subsetonlyheat = Heatmap(k36_subsetonly, name = "H3K36me3", 
                             column_order =  colnames(k36_subsetonly), 
                             show_row_names = T,
                             col = K36_col_fun,
                             column_title = "K36_1",
                             row_names_gp = gpar(fontsize=6),
                             column_labels = colnames(k36_subsetonly),
                             row_names_side = "left",
                             column_names_gp = gpar(fontsize=12),
                             use_raster = T, 
                             raster_quality = 8, 
                             row_km = 3,
                             border = TRUE,
                             rect_gp = gpar(col = "black", lwd = 2))



k4_subsetonlyheat = Heatmap(k4_subsetonly, name = "H3K4me3", 
                            column_order =  colnames(k4_subsetonly), 
                            show_row_names = T,
                            col = K4_col_fun,
                            column_title = "K4_1",
                            row_names_gp = gpar(fontsize=6),
                            column_labels = colnames(k4_subsetonly),
                            row_names_side = "left",
                            column_names_gp = gpar(fontsize=12),
                            use_raster = T, 
                            raster_quality = 8, 
                            row_km = 3,
                            border = TRUE,
                            rect_gp = gpar(col = "black", lwd = 2))


minorm_ICRsubsetonlyheat + dnmt3b_subsetonlyheat + k36_subsetonlyheat + k4_subsetonlyheat
#minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR1_colsubsetonly.pdf
#minormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR1_colsubsetonly.svg

inorm_ICRsubsetonlyheat + dnmt3b_subsetonlyheat + k36_subsetonlyheat + k4_subsetonlyheat
#inormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR1_colsubsetonly.pdf
#inormdataICR1_counted2ratioavgDMR_dnmt3b_k36_k4_subset_diffimplocipGpR1_colsubsetonly.svg


#Pie Chart
head(curhg38Manifest)
curhg38Manifest_features <- curhg38Manifest[,c(1,2,3,4,23,30)]
head(curhg38Manifest_features)
write.table(curhg38Manifest_features, "curhg38Manifest_features.txt", quote=F, sep="\t", append = F, row.names = F)
#Pie Chart
awk '{print $5}' miMERGE_myiNorm3re_chrpos_avg_PG.txt | sort -k1,1 -u > miMERGE_myiNorm3re_chrpos_avg_PG.probesID
awk '{print $5}' miMERGE_myiNorm3re_chrpos_avg_PR.txt | sort -k1,1 -u > miMERGE_myiNorm3re_chrpos_avg_PR.probesID

fgrep -f miMERGE_myiNorm3re_chrpos_avg_PG.probesID curhg38Manifest_features.txt -w > miMERGE_myiNorm3re_chrpos_avg_PG.features.txt
fgrep -f miMERGE_myiNorm3re_chrpos_avg_PR.probesID curhg38Manifest_features.txt -w > miMERGE_myiNorm3re_chrpos_avg_PR.features.txt

#Possible features Relation_to_UCSC_CpG_Island
unique(curhg38Manifest_features$Relation_to_UCSC_CpG_Island)
count_curhg38Manifest_features <- count(curhg38Manifest_features,"Relation_to_UCSC_CpG_Island")

miMERGE_myiNorm3re_chrpos_avg_PG.features <- read.delim("miMERGE_myiNorm3re_chrpos_avg_PG.features.txt", header = F)
miMERGE_myiNorm3re_chrpos_avg_PG.features <- data.frame(miMERGE_myiNorm3re_chrpos_avg_PG.features)
unique(miMERGE_myiNorm3re_chrpos_avg_PG.features$V5)

count_miMERGE_myiNorm3re_chrpos_avg_PG.features <- count(miMERGE_myiNorm3re_chrpos_avg_PG.features,"V5")
count_miMERGE_myiNorm3re_chrpos_avg_PG.features["col"] <- c("aOthers", "Island","N_Shelf","N_Shore","S_Shelf","S_Shore")
rownames(count_miMERGE_myiNorm3re_chrpos_avg_PG.features) <- count_miMERGE_myiNorm3re_chrpos_avg_PG.features$col
count_miMERGE_myiNorm3re_chrpos_avg_PG.features <- count_miMERGE_myiNorm3re_chrpos_avg_PG.features[,-1]
library(ggplot2)
# Barplot-Pie
pie_miMERGE_myiNorm3re_chrpos_avg_PG.features <- ggplot(count_miMERGE_myiNorm3re_chrpos_avg_PG.features, aes(x="", y=freq, fill=col))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)
pie_miMERGE_myiNorm3re_chrpos_avg_PG.features + scale_fill_brewer(palette="Reds")+
  theme_minimal()
ggsave("pie_miMERGE_myiNorm3re_chrpos_avg_PG.features.svg", width=12, height=10, units="cm", dpi=96)

miMERGE_myiNorm3re_chrpos_avg_PR.features <- read.delim("miMERGE_myiNorm3re_chrpos_avg_PR.features.txt", header = F)
miMERGE_myiNorm3re_chrpos_avg_PR.features <- data.frame(miMERGE_myiNorm3re_chrpos_avg_PR.features)
unique(miMERGE_myiNorm3re_chrpos_avg_PR.features$V5)

count_miMERGE_myiNorm3re_chrpos_avg_PR.features <- count(miMERGE_myiNorm3re_chrpos_avg_PR.features,"V5")
count_miMERGE_myiNorm3re_chrpos_avg_PR.features["col"] <- c("aOthers", "Island","N_Shelf","N_Shore","S_Shelf","S_Shore")
rownames(count_miMERGE_myiNorm3re_chrpos_avg_PR.features) <- count_miMERGE_myiNorm3re_chrpos_avg_PR.features$col
count_miMERGE_myiNorm3re_chrpos_avg_PR.features <- count_miMERGE_myiNorm3re_chrpos_avg_PR.features[,-1]
library(ggplot2)
# Barplot-Pie
pie_miMERGE_myiNorm3re_chrpos_avg_PR.features <- ggplot(count_miMERGE_myiNorm3re_chrpos_avg_PR.features, aes(x="", y=freq, fill=col))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)
pie_miMERGE_myiNorm3re_chrpos_avg_PR.features + scale_fill_brewer(palette="Reds")+
  theme_minimal()
ggsave("pie_miMERGE_myiNorm3re_chrpos_avg_PR.features.svg", width=12, height=10, units="cm", dpi=96)


#Regulatory Features Group
unique(curhg38Manifest_features$Regulatory_Feature_Group)
count_curhg38Manifest_regulatory_features <- count(curhg38Manifest_features,"Regulatory_Feature_Group")

miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features <- data.frame(miMERGE_myiNorm3re_chrpos_avg_PG.features)
unique(miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features$V6)

count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features <- count(miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features,"V6")
count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features["col"] <- c("aOthers",as.character(count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features$V6)[2:length(count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features$V6)])
rownames(count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features) <- count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features$col
count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features <- count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features[,-1]
library(ggplot2)
# Barplot-Pie
pie_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features <- ggplot(count_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features, aes(x="", y=freq, fill=col))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)
pie_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features + scale_fill_brewer(palette="Reds")+
  theme_minimal()
ggsave("pie_miMERGE_myiNorm3re_chrpos_avg_PG.regulatory_features.svg", width=20, height=18, units="cm", dpi=96)

miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features <- data.frame(miMERGE_myiNorm3re_chrpos_avg_PR.features)
unique(miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features$V6)

count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features <- count(miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features,"V6")
count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features["col"] <- c("aOthers",as.character(count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features$V6)[2:length(count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features$V6)])
rownames(count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features) <- count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features$col
count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features <- count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features[,-1]
library(ggplot2)
# Barplot-Pie
pie_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features <- ggplot(count_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features, aes(x="", y=freq, fill=col))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)
pie_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features + scale_fill_brewer(palette="Reds")+
  theme_minimal()
ggsave("pie_miMERGE_myiNorm3re_chrpos_avg_PR.regulatory_features.svg", width=20, height=18, units="cm", dpi=96)


#Scatter Plot
head(inormdataBin500bp1_counted1)
dim(inormdataBin500bp1_counted1)
inormdataBin500bp1_counted1ratio <- as.matrix(inormdataBin500bp1_counted1[,c(5,8,6,7,1:4)])

head(inormdataBin500bp1_counted1ratio)
dim(inormdataBin500bp1_counted1ratio)
winormdataBin500bp1_counted1ratioavg <- data.frame(cbind((rowMeans(inormdataBin500bp1_counted1ratio[,1:2])),
                                                         (inormdataBin500bp1_counted1ratio[,1:8])))

head(winormdataBin500bp1_counted1ratioavg)
colnames(winormdataBin500bp1_counted1ratioavg) <- c("Allcontrol",colnames(inormdataBin500bp1_counted1ratio))
head(winormdataBin500bp1_counted1ratioavg)
winormdataBin500bp1_counted1ratioavg <- as.matrix(winormdataBin500bp1_counted1ratioavg)
head(winormdataBin500bp1_counted1ratioavg)
dim(winormdataBin500bp1_counted1ratioavg)
winormdataBin500bp1_counted1ratioavg <- data.frame(winormdataBin500bp1_counted1ratioavg)

head(inormdataICR1_counted2)
dim(inormdataICR1_counted2)
inormdataICR1_counted3 <- inormdataICR1_counted2[,c(5,8,6,7,1:4)]
inormdataICR1_counted3ratio <- as.matrix(inormdataICR1_counted3)

head(inormdataICR1_counted3ratio)
dim(inormdataICR1_counted3ratio)
winormdataICR1_counted3ratioavg <- data.frame(cbind((rowMeans(inormdataICR1_counted3ratio[,1:2])),
                                                    (inormdataICR1_counted3ratio[,1:8])))

head(winormdataICR1_counted3ratioavg)
colnames(winormdataICR1_counted3ratioavg) <- c("Allcontrol",colnames(inormdataICR1_counted3ratio))
head(winormdataICR1_counted3ratioavg)
winormdataICR1_counted3ratioavg <- as.matrix(winormdataICR1_counted3ratioavg)
head(winormdataICR1_counted3ratioavg)
dim(winormdataICR1_counted3ratioavg)
winormdataICR1_counted3ratioavg <- data.frame(winormdataICR1_counted3ratioavg)

#500bp
with(winormdataBin500bp1_counted1ratioavg, plot(Allcontrol, PG, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#FFB6B3", bg="#FFB6B3",ylab="pG (Beta-value)", xlab="Control (Beta-value)", bty = 'n'))
with(subset(winormdataBin500bp1_counted1ratioavg), points(Allcontrol, C13, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdataBin500bp1_counted1ratioavg), points(Allcontrol, C50, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') #reddish


with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, PG, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="magenta", bg="magenta"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, C13, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, C50, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish

lines(x = c(0,1), y = c(0,1), col = "#464647", lty = 1,lwd=2)
#save as winormdataBin500bp1_counted1ratioavg_PG.png


with(winormdataBin500bp1_counted1ratioavg, plot(Allcontrol, PR, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#FFB6B3", bg="#FFB6B3",ylab="PR (Beta-value)", xlab="Control (Beta-value)", bty = 'n'))
with(subset(winormdataBin500bp1_counted1ratioavg), points(Allcontrol, C7, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdataBin500bp1_counted1ratioavg), points(Allcontrol, C35, pch=21, cex = 0.2, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') #reddish

with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, PR, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="magenta", bg="magenta"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, C7, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish
with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, C35, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish

lines(x = c(0,1), y = c(0,1), col = "#464647", lty = 1,lwd=2)
#save as winormdataBin500bp1_counted1ratioavg_PR.png

#Scatter with only diff meth specific probes
#Also add diff meth case specific probes
dim(inormaggregatedBin500bp_countedBin500bp_PG)
head(inormaggregatedBin500bp_countedBin500bp_PG)
miMERGE_myiBin500bp_chrpos_avg_PG <- inormaggregatedBin500bp_countedBin500bp_PG

miMERGE_myiBin500bp_chrpos_avg_PG <- inormaggregatedBin500bp_countedBin500bp_PG[,c(5,8,6,7,1:4)]
miMERGE_myiBin500bp_chrpos_avg_PGratio <- as.matrix(miMERGE_myiBin500bp_chrpos_avg_PG)
head(miMERGE_myiBin500bp_chrpos_avg_PGratio)
dim(miMERGE_myiBin500bp_chrpos_avg_PGratio)
wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg <- data.frame(cbind((rowMeans(miMERGE_myiBin500bp_chrpos_avg_PGratio[,1:2])),
                                                               (miMERGE_myiBin500bp_chrpos_avg_PGratio[,1:8])))

head(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg)
colnames(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg) <- c("Allcontrol",colnames(miMERGE_myiBin500bp_chrpos_avg_PGratio))
head(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg)
wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg <- as.matrix(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg)
head(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg)
dim(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg)
wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg <- data.frame(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg)
with(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg, plot(Allcontrol, PG, pch=21, cex = 0.3, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#FFB6B3", bg="#FFB6B3",ylab="pG (Beta-value)", xlab="Control (Beta-value)", bty = 'n'))
with(subset(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg), points(Allcontrol, C13, pch=21, cex = 0.3, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') 
with(subset(wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg), points(Allcontrol, C50, pch=21, cex = 0.3, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') 


with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, PG, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="magenta", bg="magenta"),ylab="", xlab="", bty = 'n') 
with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, C13, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') 
with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, C50, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') 

lines(x = c(0,1), y = c(0,1), col = "#464647", lty = 1,lwd=2)
#save as wmiMERGE_myiBin500bp_chrpos_avg_PGratioavg.png

dim(inormaggregatedBin500bp_countedBin500bp_PR)
head(inormaggregatedBin500bp_countedBin500bp_PR)
miMERGE_myiBin500bp_chrpos_avg_PR <- inormaggregatedBin500bp_countedBin500bp_PR

miMERGE_myiBin500bp_chrpos_avg_PR <- inormaggregatedBin500bp_countedBin500bp_PR[,c(5,8,6,7,1:4)]
miMERGE_myiBin500bp_chrpos_avg_PRratio <- as.matrix(miMERGE_myiBin500bp_chrpos_avg_PR)
head(miMERGE_myiBin500bp_chrpos_avg_PRratio)
dim(miMERGE_myiBin500bp_chrpos_avg_PRratio)
wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg <- data.frame(cbind((rowMeans(miMERGE_myiBin500bp_chrpos_avg_PRratio[,1:2])),
                                                               (miMERGE_myiBin500bp_chrpos_avg_PRratio[,1:8])))

head(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg)
colnames(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg) <- c("Allcontrol",colnames(miMERGE_myiBin500bp_chrpos_avg_PRratio))
head(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg)
wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg <- as.matrix(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg)
head(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg)
dim(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg)
wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg <- data.frame(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg)
with(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg, plot(Allcontrol, PR, pch=21, cex = 0.3, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#FFB6B3", bg="#FFB6B3",ylab="PR (Beta-value)", xlab="Control (Beta-value)", bty = 'n'))
with(subset(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg), points(Allcontrol, C13, pch=21, cex = 0.3, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') 
with(subset(wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg), points(Allcontrol, C50, pch=21, cex = 0.3, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="#BDE7BD", bg="#BDE7BD"),ylab="", xlab="", bty = 'n') 


with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, PR, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="magenta", bg="magenta"),ylab="", xlab="", bty = 'n') 
with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, C13, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') 
with(subset(winormdataICR1_counted3ratioavg), points(Allcontrol, C50, pch=21, cex = 0.4, xlim=c(0,1), ylim=c(0,1), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') 

lines(x = c(0,1), y = c(0,1), col = "#464647", lty = 1,lwd=2)
#save as wmiMERGE_myiBin500bp_chrpos_avg_PRratioavg.png

#################  END OF  ANALYSIS ####################














#MeBin 500 bp

bedtools intersect -wa -wb -a MERGE_myiNorm3re_pos_chr.txt -b UN_tiles_with_meth_level_80perc.txt > MERGE_myiNorm3re_human_MeBin.txt
awk '{print $13"%"$14"%"$15"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_human_MeBin.txt > MERGE_myiNorm3re_human_MeBin.rearranged.txt


#PCA with  Human MeBin
inormdataMeBin <- read.table("MERGE_myiNorm3re_human_MeBin.rearranged.txt", header = F, stringsAsFactors = F)
colnames(inormdataMeBin) <- c("MeBin","chr","start","end","TargetID", colnames(myNorm2))
rownames(inormdataMeBin)
colnames(inormdataMeBin)
head(inormdataMeBin)
dim(inormdataMeBin)
Count_inormdataMeBin1 <- count(inormdataMeBin, "MeBin")
head(Count_inormdataMeBin1)
dim(Count_inormdataMeBin1)
snormaggregateMeBin = aggregate(inormdataMeBin[,6:13],by=list(inormdataMeBin$MeBin), mean)
head(snormaggregateMeBin, 2)
#Plot PCA by group color and labelling
inormdataMeBin1=snormaggregateMeBin
rownames(inormdataMeBin1)
inormdataMeBin1[,1]
rownames(inormdataMeBin1)=inormdataMeBin1[,1]
rownames(inormdataMeBin1)
colnames(inormdataMeBin1)
inormdataMeBin1 = inormdataMeBin1[,-1]
head(inormdataMeBin1)
dim(inormdataMeBin1)
inormdataMeBin1_counted <- cbind(inormdataMeBin1, Count_inormdataMeBin1)
head(inormdataMeBin1_counted)
dim(inormdataMeBin1_counted)
tail(inormdataMeBin1_counted)
inormdataMeBin1_counted1 <- inormdataMeBin1_counted[which(inormdataMeBin1_counted$freq >= 5),]
head(inormdataMeBin1_counted1)
dim(inormdataMeBin1_counted1)
write.table(inormdataMeBin1_counted1, "inormdataMeBin1_counted1.filtmin.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)


inormdataMeBin1_counted2 <- inormdataMeBin1_counted1[,1:8]
head(inormdataMeBin1_counted2)
summary(inormdataMeBin1_counted2)
#df <- as.inormdata.frame(inormdataMeBin1)
inormdfMeBin <- inormdataMeBin1_counted2
head(inormdfMeBin)
dim(inormdfMeBin)
inormdfMeBin = data.frame(inormdfMeBin)
#write.table(inormdfMeBin , "inormdfMeBindedup.txt", sep="\t", quote = FALSE, append = FALSE)
tinormdfMeBin = t(inormdfMeBin)
tinormdfMeBin = data.frame(tinormdfMeBin)
#write.table(tinormdfMeBin , "tinormdfMeBindedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfMeBin)
tinormdfMeBin["Color"] <-  c("corrClone","corrClone","corrClone","corrClone","Control","Case","Case","Control")
#head(tinormdfMeBin)
#write.table(tinormdfMeBin , "tinormdfMeBindedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfMeBin)
inormdfx <-tinormdfMeBin[c(1:(length(tinormdfMeBin)-1))]
PMeBinC<-prcomp(inormdfx, center = TRUE, scale. = TRUE)
PMeBinCi<-data.frame(PMeBinC$x,Color=tinormdfMeBin$Color)
percentagecomMeBin <- round(PMeBinC$sdev^2 / sum(PMeBinC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomMeBin <- paste( colnames(PMeBinCi), "(", paste( as.character(percentagecomMeBin), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pMeBin1 <-ggplot(PMeBinCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomMeBin[1]) + ylab(percentagecomMeBin[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("#FFB6B3","black","#BDE7BD","navy","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","navy","navy","pink","navy","pink","pink","navy","pink","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  scale_shape_manual(values=c(0,1,2,0,1,1,1,1,25,25,25,25,25,25,25,22,19,19,22,21,22,22,21,22,24,24,24,24,24,24,23,23,23,23))
pMeBin1 <- pMeBin1 +theme_classic()
pMeBin1 + xlim(-500,500)+ ylim(-500,500)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_tinormdfMeBin_CntrlIndiv.svg", width=10*1.25, height=8*1.25, units="cm", dpi=96)


################################### Heatmap_indiv human MeBin ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(inormdataMeBin1_counted1)
dim(inormdataMeBin1_counted1)
sum(inormdataMeBin1_counted1$freq)

#Heatmap_indiv with 2 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
inormdataMeBin2 = as.matrix(inormdataMeBin1_counted2)
head(inormdataMeBin2)
dim(inormdataMeBin2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
inormdataMeBin3 <- inormdataMeBin2
head(inormdataMeBin3)
write.table(inormdataMeBin3, "inormdataMeBin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
#Box and Violin Plot human MeBin using 12 Controls from MethHed
rownames(inormdataMeBin)
colnames(inormdataMeBin)
head(inormdataMeBin)
dim(inormdataMeBin)
Count_inormdataMeBin1 <- count(inormdataMeBin, "MeBin")
head(Count_inormdataMeBin1)
inormagregateMeBin1 = aggregate(inormdataMeBin[,6:13],by=list(inormdataMeBin$MeBin), mean)
head(inormagregateMeBin1, 2)
inormdataMeBin1=inormagregateMeBin1
rownames(inormdataMeBin1)
inormdataMeBin1[,1]
rownames(inormdataMeBin1)=inormdataMeBin1[,1]
rownames(inormdataMeBin1)
colnames(inormdataMeBin1)
inormdataMeBin1 = inormdataMeBin1[,-1]
head(inormdataMeBin1)
dim(inormdataMeBin1)
write.table(inormdataMeBin1, "inormaggregatedMeBin_inormdataMeBin1.txt", sep="\t", quote = FALSE, append = FALSE)
inormdataMeBin1_counted <- cbind(inormdataMeBin1, Count_inormdataMeBin1)
head(inormdataMeBin1_counted)
dim(inormdataMeBin1_counted)
inormdataMeBin1_counted1 <- inormdataMeBin1_counted[which(inormdataMeBin1_counted$freq >= 5),]
head(inormdataMeBin1_counted1)
dim(inormdataMeBin1_counted1)
inormdataMeBin1_counted2 <- inormdataMeBin1_counted1[,1:8]
head(inormdataMeBin1_counted2)
VionormIndvMeBin <- data.frame(inormdataMeBin1_counted2[,c(5,8,6,7,1:4)])
head(VionormIndvMeBin)
write.table(VionormIndvMeBin, "VionormIndvMeBin.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_VionormIndvMeBin.svg", width=10, height=5, pointsize=12)
boxplot(VionormIndvMeBin, main="Average methylation at human MeBins", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","navy","navy","navy","navy"))
dev.off()
dim(VionormIndvMeBin)
VionormIndvMeBin <- data.frame(VionormIndvMeBin)
VionormIndvMeBin1 <- VionormIndvMeBin
head(VionormIndvMeBin1,1)
VionormIndvMeBin2 <- stack(VionormIndvMeBin1)
head(VionormIndvMeBin2)
colnames(VionormIndvMeBin2) <- c("Methylation", "inormdatasets")
head(VionormIndvMeBin2)
ggplot(VionormIndvMeBin2, aes(x=inormdatasets, y=Methylation, color=inormdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","#FFB6B3","#FFB6B3","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","navy","navy","navy"))
ggsave("Violin_plot_VionormIndvMeBin2_pimethhed_humanMeBin.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)



#Heatmap of control normalized Log other samples merged
head(inormdataMeBin1_counted2)
dim(inormdataMeBin1_counted2)
inormdataMeBin1_counted2ratio <- as.matrix(inormdataMeBin1_counted2[,c(5,8,6,7,1:4)])

head(inormdataMeBin1_counted2ratio)
dim(inormdataMeBin1_counted2ratio)
winormdataMeBin1_counted2ratioavg <- data.frame(cbind((rowMeans(inormdataMeBin1_counted2ratio[,1:2])),
                                                        (inormdataMeBin1_counted2ratio[,1:8])))

head(winormdataMeBin1_counted2ratioavg)
colnames(winormdataMeBin1_counted2ratioavg) <- c("Allcontrol",colnames(inormdataMeBin1_counted2ratio))
head(winormdataMeBin1_counted2ratioavg)
winormdataMeBin1_counted2ratioavg1 <- winormdataMeBin1_counted2ratioavg
winormdataMeBin1_counted2ratioavg2 <- as.matrix(winormdataMeBin1_counted2ratioavg1)
head(winormdataMeBin1_counted2ratioavg2)
dim(winormdataMeBin1_counted2ratioavg2)
winormdataMeBin1_counted2ratioavg3 <- data.frame(cbind((winormdataMeBin1_counted2ratioavg2[,1]/winormdataMeBin1_counted2ratioavg2[,1]),
                                                         (winormdataMeBin1_counted2ratioavg2[,2]/winormdataMeBin1_counted2ratioavg2[,1]),
                                                         (winormdataMeBin1_counted2ratioavg2[,3]/winormdataMeBin1_counted2ratioavg2[,1]),
                                                         (winormdataMeBin1_counted2ratioavg2[,4]/winormdataMeBin1_counted2ratioavg2[,1]),
                                                         (winormdataMeBin1_counted2ratioavg2[,5]/winormdataMeBin1_counted2ratioavg2[,1]),
                                                         (winormdataMeBin1_counted2ratioavg2[,6]/winormdataMeBin1_counted2ratioavg2[,1]),
                                                         (winormdataMeBin1_counted2ratioavg2[,7]/winormdataMeBin1_counted2ratioavg2[,1]),
                                                         (winormdataMeBin1_counted2ratioavg2[,8]/winormdataMeBin1_counted2ratioavg2[,1]),
                                                         (winormdataMeBin1_counted2ratioavg2[,9]/winormdataMeBin1_counted2ratioavg2[,1])))


head(winormdataMeBin1_counted2ratioavg3)
colnames(winormdataMeBin1_counted2ratioavg3) <- colnames(winormdataMeBin1_counted2ratioavg2)
head(winormdataMeBin1_counted2ratioavg3)
winormdataMeBin1_counted2ratioavg3 <- as.matrix(winormdataMeBin1_counted2ratioavg3)
winormdataMeBin1_counted2ratioavg3Log <- log2(winormdataMeBin1_counted2ratioavg3)
head(winormdataMeBin1_counted2ratioavg3Log)
write.table(winormdataMeBin1_counted2ratioavg3Log , "Heatmap_ICF_normcontrolwinormdataMeBin1_counted2ratioavg2Log.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)

winormdataMeBin1_counted2ratioavg3Logminuscontrol <- winormdataMeBin1_counted2ratioavg3Log[,c(2:9)]
head(winormdataMeBin1_counted2ratioavg3Logminuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(winormdataMeBin1_counted2ratioavg3Logminuscontrol))
row.names(my_sample_col2) <- colnames(winormdataMeBin1_counted2ratioavg3Logminuscontrol)
my_colour2 = list(Annotations = c("D250"= "black","UN"= "black","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(winormdataMeBin1_counted2ratioavg3Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_winormdataMeBin1_counted2ratioavg.normcontrol.svg

minormdataMeBin1_counted2ratioavg3 <- data.frame(cbind((winormdataMeBin1_counted2ratioavg2[,1:9]-winormdataMeBin1_counted2ratioavg2[,1])))


head(minormdataMeBin1_counted2ratioavg3)
minormdataMeBin1_counted2ratioavg3 <- as.matrix(minormdataMeBin1_counted2ratioavg3)

minormdataMeBin1_counted2ratioavg3minuscontrol <- minormdataMeBin1_counted2ratioavg3[,c(2:9)]
head(minormdataMeBin1_counted2ratioavg3minuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(minormdataMeBin1_counted2ratioavg3minuscontrol))
row.names(my_sample_col2) <- colnames(minormdataMeBin1_counted2ratioavg3minuscontrol)
my_colour2 = list(Annotations = c("D250"= "black","UN"= "black","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(minormdataMeBin1_counted2ratioavg3minuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataMeBin1_counted2ratioavg.normcontrol.svg

pheatmap(minormdataMeBin1_counted2ratioavg3minuscontrol[,1:4],
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataMeBin1_counted2ratioavg.normcontrol_casecntrl.svg

#Regions where two controls are close
pinormdataMeBin1_counted2ratioavg2 <- data.frame(winormdataMeBin1_counted2ratioavg2)
pinormdataMeBin1_counted2ratioavg2["diffcntrl"] <- pinormdataMeBin1_counted2ratioavg2[,2] - pinormdataMeBin1_counted2ratioavg2[,1]
head(pinormdataMeBin1_counted2ratioavg2)
dim(pinormdataMeBin1_counted2ratioavg2)
write.table(pinormdataMeBin1_counted2ratioavg2,"pinormdataMeBin1_counted2ratioavg2.txt", sep = "\t", append = F, quote=F)
pinormdataMeBin1_counted2ratioavg3 <- pinormdataMeBin1_counted2ratioavg2[which(pinormdataMeBin1_counted2ratioavg2$diffcntrl < 0.1 & 
                                                                                 pinormdataMeBin1_counted2ratioavg2$diffcntrl > -0.1),]

head(pinormdataMeBin1_counted2ratioavg3)
dim(pinormdataMeBin1_counted2ratioavg3)
minormdataMeBin1_counted2ratiodiff3 <- data.frame(cbind((pinormdataMeBin1_counted2ratioavg3[,1:9]-pinormdataMeBin1_counted2ratioavg3[,1])))


head(minormdataMeBin1_counted2ratiodiff3)

minormdataMeBin1_counted2ratiodiff3 <- minormdataMeBin1_counted2ratiodiff3[which(minormdataMeBin1_counted2ratiodiff3$PG > 0.2 | 
                                                                                   minormdataMeBin1_counted2ratiodiff3$PG < -0.2),]

minormdataMeBin1_counted2ratiodiff3 <- minormdataMeBin1_counted2ratiodiff3[which(minormdataMeBin1_counted2ratiodiff3$PR > 0.2 | 
                                                                                   minormdataMeBin1_counted2ratiodiff3$PR < -0.2),]

head(minormdataMeBin1_counted2ratiodiff3)
dim(minormdataMeBin1_counted2ratiodiff3)
minormdataMeBin1_counted2ratiodiff3 <- as.matrix(minormdataMeBin1_counted2ratiodiff3)

minormdataMeBin1_counted2ratiodiff3minuscontrol <- minormdataMeBin1_counted2ratiodiff3[,c(2:9)]
head(minormdataMeBin1_counted2ratiodiff3minuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(minormdataMeBin1_counted2ratiodiff3minuscontrol))
row.names(my_sample_col2) <- colnames(minormdataMeBin1_counted2ratiodiff3minuscontrol)
my_colour2 = list(Annotations = c("D250"= "black","UN"= "black","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(minormdataMeBin1_counted2ratiodiff3minuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataMeBin1_counted2ratiodiff3.normcontrol.svg

pheatmap(minormdataMeBin1_counted2ratiodiff3minuscontrol[,1:4],
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataMeBin1_counted2ratiodiff3.normcontrol_casecntrl.svg

library(splitstackshape)
minormdataMeBin1_counted2ratiodiff3minuscontrolchr <- data.frame(minormdataMeBin1_counted2ratiodiff3minuscontrol)
minormdataMeBin1_counted2ratiodiff3minuscontrolchr["pos"] <- rownames(minormdataMeBin1_counted2ratiodiff3minuscontrolchr)
minormdataMeBin1_counted2ratiodiff3minuscontrolchr <- cSplit(minormdataMeBin1_counted2ratiodiff3minuscontrolchr, "pos", "%")
head(minormdataMeBin1_counted2ratiodiff3minuscontrolchr)
dim(minormdataMeBin1_counted2ratiodiff3minuscontrolchr)
minormdataMeBin1_counted2ratiodiff3minuscontrolchr <- minormdataMeBin1_counted2ratiodiff3minuscontrolchr[,c(9:11,1:8)]
write.table(minormdataMeBin1_counted2ratiodiff3minuscontrolchr , "minormdataMeBin1_counted2ratiodiff3minuscontrolchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = F)


bedtools intersect -wa -wb -a minormdataMeBin1_counted2ratiodiff3minuscontrolchr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed

#CGIs_MeBin 500 bp

bedtools intersect -wa -wb -a MERGE_myiNorm3re_pos_chr.txt -b CpGi_islands_with_UN_meth_80percent.txt > MERGE_myiNorm3re_human_CGIs_MeBin.txt
awk '{print $13"%"$14"%"$15"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' MERGE_myiNorm3re_human_CGIs_MeBin.txt > MERGE_myiNorm3re_human_CGIs_MeBin.rearranged.txt


#PCA with  Human CGIs_MeBin
inormdataCGIs_MeBin <- read.table("MERGE_myiNorm3re_human_CGIs_MeBin.rearranged.txt", header = F, stringsAsFactors = F)
colnames(inormdataCGIs_MeBin) <- c("CGIs_MeBin","chr","start","end","TargetID", colnames(myNorm2))
rownames(inormdataCGIs_MeBin)
colnames(inormdataCGIs_MeBin)
head(inormdataCGIs_MeBin)
dim(inormdataCGIs_MeBin)
Count_inormdataCGIs_MeBin1 <- count(inormdataCGIs_MeBin, "CGIs_MeBin")
head(Count_inormdataCGIs_MeBin1)
dim(Count_inormdataCGIs_MeBin1)
snormaggregateCGIs_MeBin = aggregate(inormdataCGIs_MeBin[,6:13],by=list(inormdataCGIs_MeBin$CGIs_MeBin), mean)
head(snormaggregateCGIs_MeBin, 2)
#Plot PCA by group color and labelling
inormdataCGIs_MeBin1=snormaggregateCGIs_MeBin
rownames(inormdataCGIs_MeBin1)
inormdataCGIs_MeBin1[,1]
rownames(inormdataCGIs_MeBin1)=inormdataCGIs_MeBin1[,1]
rownames(inormdataCGIs_MeBin1)
colnames(inormdataCGIs_MeBin1)
inormdataCGIs_MeBin1 = inormdataCGIs_MeBin1[,-1]
head(inormdataCGIs_MeBin1)
dim(inormdataCGIs_MeBin1)
inormdataCGIs_MeBin1_counted <- cbind(inormdataCGIs_MeBin1, Count_inormdataCGIs_MeBin1)
head(inormdataCGIs_MeBin1_counted)
dim(inormdataCGIs_MeBin1_counted)
tail(inormdataCGIs_MeBin1_counted)
inormdataCGIs_MeBin1_counted1 <- inormdataCGIs_MeBin1_counted[which(inormdataCGIs_MeBin1_counted$freq >= 5),]
head(inormdataCGIs_MeBin1_counted1)
dim(inormdataCGIs_MeBin1_counted1)
write.table(inormdataCGIs_MeBin1_counted1, "inormdataCGIs_MeBin1_counted1.filtmin.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)


inormdataCGIs_MeBin1_counted2 <- inormdataCGIs_MeBin1_counted1[,1:8]
head(inormdataCGIs_MeBin1_counted2)
summary(inormdataCGIs_MeBin1_counted2)
#df <- as.inormdata.frame(inormdataCGIs_MeBin1)
inormdfCGIs_MeBin <- inormdataCGIs_MeBin1_counted2
head(inormdfCGIs_MeBin)
dim(inormdfCGIs_MeBin)
inormdfCGIs_MeBin = data.frame(inormdfCGIs_MeBin)
#write.table(inormdfCGIs_MeBin , "inormdfCGIs_MeBindedup.txt", sep="\t", quote = FALSE, append = FALSE)
tinormdfCGIs_MeBin = t(inormdfCGIs_MeBin)
tinormdfCGIs_MeBin = data.frame(tinormdfCGIs_MeBin)
#write.table(tinormdfCGIs_MeBin , "tinormdfCGIs_MeBindedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfCGIs_MeBin)
tinormdfCGIs_MeBin["Color"] <-  c("corrClone","corrClone","corrClone","corrClone","Control","Case","Case","Control")
#head(tinormdfCGIs_MeBin)
#write.table(tinormdfCGIs_MeBin , "tinormdfCGIs_MeBindedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tinormdfCGIs_MeBin)
inormdfx <-tinormdfCGIs_MeBin[c(1:(length(tinormdfCGIs_MeBin)-1))]
PCGIs_MeBinC<-prcomp(inormdfx, center = TRUE, scale. = TRUE)
PCGIs_MeBinCi<-data.frame(PCGIs_MeBinC$x,Color=tinormdfCGIs_MeBin$Color)
percentagecomCGIs_MeBin <- round(PCGIs_MeBinC$sdev^2 / sum(PCGIs_MeBinC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomCGIs_MeBin <- paste( colnames(PCGIs_MeBinCi), "(", paste( as.character(percentagecomCGIs_MeBin), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pCGIs_MeBin1 <-ggplot(PCGIs_MeBinCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomCGIs_MeBin[1]) + ylab(percentagecomCGIs_MeBin[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("#FFB6B3","black","#BDE7BD","navy","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","pink","navy","navy","pink","navy","pink","pink","navy","pink","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  scale_shape_manual(values=c(0,1,2,0,1,1,1,1,25,25,25,25,25,25,25,22,19,19,22,21,22,22,21,22,24,24,24,24,24,24,23,23,23,23))
pCGIs_MeBin1 <- pCGIs_MeBin1 +theme_classic()
pCGIs_MeBin1 + xlim(-20,20)+ ylim(-20,20)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_tinormdfCGIs_MeBin_CntrlIndiv.svg", width=10*1.25, height=8*1.25, units="cm", dpi=96)


################################### Heatmap_indiv human CGIs_MeBin ######################################
head(inormdataCGIs_MeBin1_counted1)
dim(inormdataCGIs_MeBin1_counted1)
sum(inormdataCGIs_MeBin1_counted1$freq)

#Heatmap_indiv with 2 Controls Avg Filtered by minimum 3 CpGs
library(ggplot2)
library(gplots)
inormdataCGIs_MeBin2 = as.matrix(inormdataCGIs_MeBin1_counted2)
head(inormdataCGIs_MeBin2)
dim(inormdataCGIs_MeBin2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
inormdataCGIs_MeBin3 <- inormdataCGIs_MeBin2
head(inormdataCGIs_MeBin3)
write.table(inormdataCGIs_MeBin3, "inormdataCGIs_MeBin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
#Box and Violin Plot human CGIs_MeBin using 12 Controls from MethHed
rownames(inormdataCGIs_MeBin)
colnames(inormdataCGIs_MeBin)
head(inormdataCGIs_MeBin)
dim(inormdataCGIs_MeBin)
Count_inormdataCGIs_MeBin1 <- count(inormdataCGIs_MeBin, "CGIs_MeBin")
head(Count_inormdataCGIs_MeBin1)
inormagregateCGIs_MeBin1 = aggregate(inormdataCGIs_MeBin[,6:13],by=list(inormdataCGIs_MeBin$CGIs_MeBin), mean)
head(inormagregateCGIs_MeBin1, 2)
inormdataCGIs_MeBin1=inormagregateCGIs_MeBin1
rownames(inormdataCGIs_MeBin1)
inormdataCGIs_MeBin1[,1]
rownames(inormdataCGIs_MeBin1)=inormdataCGIs_MeBin1[,1]
rownames(inormdataCGIs_MeBin1)
colnames(inormdataCGIs_MeBin1)
inormdataCGIs_MeBin1 = inormdataCGIs_MeBin1[,-1]
head(inormdataCGIs_MeBin1)
dim(inormdataCGIs_MeBin1)
write.table(inormdataCGIs_MeBin1, "inormaggregatedCGIs_MeBin_inormdataCGIs_MeBin1.txt", sep="\t", quote = FALSE, append = FALSE)
inormdataCGIs_MeBin1_counted <- cbind(inormdataCGIs_MeBin1, Count_inormdataCGIs_MeBin1)
head(inormdataCGIs_MeBin1_counted)
dim(inormdataCGIs_MeBin1_counted)
inormdataCGIs_MeBin1_counted1 <- inormdataCGIs_MeBin1_counted[which(inormdataCGIs_MeBin1_counted$freq >= 5),]
head(inormdataCGIs_MeBin1_counted1)
dim(inormdataCGIs_MeBin1_counted1)
inormdataCGIs_MeBin1_counted2 <- inormdataCGIs_MeBin1_counted1[,1:8]
head(inormdataCGIs_MeBin1_counted2)
VionormIndvCGIs_MeBin <- data.frame(inormdataCGIs_MeBin1_counted2[,c(5,8,6,7,1:4)])
head(VionormIndvCGIs_MeBin)
write.table(VionormIndvCGIs_MeBin, "VionormIndvCGIs_MeBin.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_VionormIndvCGIs_MeBin.svg", width=10, height=5, pointsize=12)
boxplot(VionormIndvCGIs_MeBin, main="Average methylation at human CGIs_MeBins", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","navy","navy","navy","navy"))
dev.off()
dim(VionormIndvCGIs_MeBin)
VionormIndvCGIs_MeBin <- data.frame(VionormIndvCGIs_MeBin)
VionormIndvCGIs_MeBin1 <- VionormIndvCGIs_MeBin
head(VionormIndvCGIs_MeBin1,1)
VionormIndvCGIs_MeBin2 <- stack(VionormIndvCGIs_MeBin1)
head(VionormIndvCGIs_MeBin2)
colnames(VionormIndvCGIs_MeBin2) <- c("Methylation", "inormdatasets")
head(VionormIndvCGIs_MeBin2)
ggplot(VionormIndvCGIs_MeBin2, aes(x=inormdatasets, y=Methylation, color=inormdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("darkgrey","darkgrey","#FFB6B3","#FFB6B3","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","#BDE7BD","navy","navy","navy"))
ggsave("Violin_plot_VionormIndvCGIs_MeBin2_pimethhed_humanCGIs_MeBin.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)



#Heatmap of control normalized Log other samples merged
head(inormdataCGIs_MeBin1_counted2)
dim(inormdataCGIs_MeBin1_counted2)
inormdataCGIs_MeBin1_counted2ratio <- as.matrix(inormdataCGIs_MeBin1_counted2[,c(5,8,6,7,1:4)])

head(inormdataCGIs_MeBin1_counted2ratio)
dim(inormdataCGIs_MeBin1_counted2ratio)
winormdataCGIs_MeBin1_counted2ratioavg <- data.frame(cbind((rowMeans(inormdataCGIs_MeBin1_counted2ratio[,1:2])),
                                                           (inormdataCGIs_MeBin1_counted2ratio[,1:8])))

head(winormdataCGIs_MeBin1_counted2ratioavg)
colnames(winormdataCGIs_MeBin1_counted2ratioavg) <- c("Allcontrol",colnames(inormdataCGIs_MeBin1_counted2ratio))
head(winormdataCGIs_MeBin1_counted2ratioavg)
winormdataCGIs_MeBin1_counted2ratioavg1 <- winormdataCGIs_MeBin1_counted2ratioavg
winormdataCGIs_MeBin1_counted2ratioavg2 <- as.matrix(winormdataCGIs_MeBin1_counted2ratioavg1)
head(winormdataCGIs_MeBin1_counted2ratioavg2)
dim(winormdataCGIs_MeBin1_counted2ratioavg2)
winormdataCGIs_MeBin1_counted2ratioavg3 <- data.frame(cbind((winormdataCGIs_MeBin1_counted2ratioavg2[,1]/winormdataCGIs_MeBin1_counted2ratioavg2[,1]),
                                                            (winormdataCGIs_MeBin1_counted2ratioavg2[,2]/winormdataCGIs_MeBin1_counted2ratioavg2[,1]),
                                                            (winormdataCGIs_MeBin1_counted2ratioavg2[,3]/winormdataCGIs_MeBin1_counted2ratioavg2[,1]),
                                                            (winormdataCGIs_MeBin1_counted2ratioavg2[,4]/winormdataCGIs_MeBin1_counted2ratioavg2[,1]),
                                                            (winormdataCGIs_MeBin1_counted2ratioavg2[,5]/winormdataCGIs_MeBin1_counted2ratioavg2[,1]),
                                                            (winormdataCGIs_MeBin1_counted2ratioavg2[,6]/winormdataCGIs_MeBin1_counted2ratioavg2[,1]),
                                                            (winormdataCGIs_MeBin1_counted2ratioavg2[,7]/winormdataCGIs_MeBin1_counted2ratioavg2[,1]),
                                                            (winormdataCGIs_MeBin1_counted2ratioavg2[,8]/winormdataCGIs_MeBin1_counted2ratioavg2[,1]),
                                                            (winormdataCGIs_MeBin1_counted2ratioavg2[,9]/winormdataCGIs_MeBin1_counted2ratioavg2[,1])))


head(winormdataCGIs_MeBin1_counted2ratioavg3)
colnames(winormdataCGIs_MeBin1_counted2ratioavg3) <- colnames(winormdataCGIs_MeBin1_counted2ratioavg2)
head(winormdataCGIs_MeBin1_counted2ratioavg3)
winormdataCGIs_MeBin1_counted2ratioavg3 <- as.matrix(winormdataCGIs_MeBin1_counted2ratioavg3)
winormdataCGIs_MeBin1_counted2ratioavg3Log <- log2(winormdataCGIs_MeBin1_counted2ratioavg3)
head(winormdataCGIs_MeBin1_counted2ratioavg3Log)
write.table(winormdataCGIs_MeBin1_counted2ratioavg3Log , "Heatmap_ICF_normcontrolwinormdataCGIs_MeBin1_counted2ratioavg2Log.txt", sep="\t", quote = FALSE, append = FALSE,row.names = T)

winormdataCGIs_MeBin1_counted2ratioavg3Logminuscontrol <- winormdataCGIs_MeBin1_counted2ratioavg3Log[,c(2:9)]
head(winormdataCGIs_MeBin1_counted2ratioavg3Logminuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(winormdataCGIs_MeBin1_counted2ratioavg3Logminuscontrol))
row.names(my_sample_col2) <- colnames(winormdataCGIs_MeBin1_counted2ratioavg3Logminuscontrol)
my_colour2 = list(Annotations = c("D250"= "black","UN"= "black","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(winormdataCGIs_MeBin1_counted2ratioavg3Logminuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_winormdataCGIs_MeBin1_counted2ratioavg.normcontrol.svg

minormdataCGIs_MeBin1_counted2ratioavg3 <- data.frame(cbind((winormdataCGIs_MeBin1_counted2ratioavg2[,1:9]-winormdataCGIs_MeBin1_counted2ratioavg2[,1])))


head(minormdataCGIs_MeBin1_counted2ratioavg3)
minormdataCGIs_MeBin1_counted2ratioavg3 <- as.matrix(minormdataCGIs_MeBin1_counted2ratioavg3)

minormdataCGIs_MeBin1_counted2ratioavg3minuscontrol <- minormdataCGIs_MeBin1_counted2ratioavg3[,c(2:9)]
head(minormdataCGIs_MeBin1_counted2ratioavg3minuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(minormdataCGIs_MeBin1_counted2ratioavg3minuscontrol))
row.names(my_sample_col2) <- colnames(minormdataCGIs_MeBin1_counted2ratioavg3minuscontrol)
my_colour2 = list(Annotations = c("D250"= "black","UN"= "black","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(minormdataCGIs_MeBin1_counted2ratioavg3minuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataCGIs_MeBin1_counted2ratioavg.normcontrol.svg

pheatmap(minormdataCGIs_MeBin1_counted2ratioavg3minuscontrol[,1:4],
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataCGIs_MeBin1_counted2ratioavg.normcontrol_casecntrl.svg

#Regions where two controls are close
pinormdataCGIs_MeBin1_counted2ratioavg2 <- data.frame(winormdataCGIs_MeBin1_counted2ratioavg2)
pinormdataCGIs_MeBin1_counted2ratioavg2["diffcntrl"] <- pinormdataCGIs_MeBin1_counted2ratioavg2[,2] - pinormdataCGIs_MeBin1_counted2ratioavg2[,1]
head(pinormdataCGIs_MeBin1_counted2ratioavg2)
dim(pinormdataCGIs_MeBin1_counted2ratioavg2)
write.table(pinormdataCGIs_MeBin1_counted2ratioavg2,"pinormdataCGIs_MeBin1_counted2ratioavg2.txt", sep = "\t", append = F, quote=F)
pinormdataCGIs_MeBin1_counted2ratioavg3 <- pinormdataCGIs_MeBin1_counted2ratioavg2[which(pinormdataCGIs_MeBin1_counted2ratioavg2$diffcntrl < 0.1 & 
                                                                                           pinormdataCGIs_MeBin1_counted2ratioavg2$diffcntrl > -0.1),]

head(pinormdataCGIs_MeBin1_counted2ratioavg3)
dim(pinormdataCGIs_MeBin1_counted2ratioavg3)
minormdataCGIs_MeBin1_counted2ratiodiff3 <- data.frame(cbind((pinormdataCGIs_MeBin1_counted2ratioavg3[,1:9]-pinormdataCGIs_MeBin1_counted2ratioavg3[,1])))


head(minormdataCGIs_MeBin1_counted2ratiodiff3)

minormdataCGIs_MeBin1_counted2ratiodiff3 <- minormdataCGIs_MeBin1_counted2ratiodiff3[which(minormdataCGIs_MeBin1_counted2ratiodiff3$PG > 0.2 | 
                                                                                             minormdataCGIs_MeBin1_counted2ratiodiff3$PG < -0.2),]

minormdataCGIs_MeBin1_counted2ratiodiff3 <- minormdataCGIs_MeBin1_counted2ratiodiff3[which(minormdataCGIs_MeBin1_counted2ratiodiff3$PR > 0.2 | 
                                                                                             minormdataCGIs_MeBin1_counted2ratiodiff3$PR < -0.2),]

head(minormdataCGIs_MeBin1_counted2ratiodiff3)
dim(minormdataCGIs_MeBin1_counted2ratiodiff3)
minormdataCGIs_MeBin1_counted2ratiodiff3 <- as.matrix(minormdataCGIs_MeBin1_counted2ratiodiff3)

minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrol <- minormdataCGIs_MeBin1_counted2ratiodiff3[,c(2:9)]
head(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrol)
library(pheatmap)
my_sample_col2 <- data.frame(Annotations= colnames(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrol))
row.names(my_sample_col2) <- colnames(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrol)
my_colour2 = list(Annotations = c("D250"= "black","UN"= "black","PG"= "#FFB6B3","PR"= "#FFB6B3","C7"= "#BDE7BD","C13"= "#BDE7BD","C35"= "#BDE7BD","C50"= "#BDE7BD"))

breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrol,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataCGIs_MeBin1_counted2ratiodiff3.normcontrol.svg

pheatmap(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrol[,1:4],
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 2)

#save as pHeatmap_minormdataCGIs_MeBin1_counted2ratiodiff3.normcontrol_casecntrl.svg

library(splitstackshape)
minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr <- data.frame(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrol)
minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr["pos"] <- rownames(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr)
minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr <- cSplit(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr, "pos", "%")
head(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr)
dim(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr)
minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr <- minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr[,c(9:11,1:8)]
write.table(minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr , "minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr.txt", sep="\t", quote = FALSE, append = FALSE,row.names = F, col.names = F)


bedtools intersect -wa -wb -a minormdataCGIs_MeBin1_counted2ratiodiff3minuscontrolchr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed
