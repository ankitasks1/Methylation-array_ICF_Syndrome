print("Methylation Array Data: ICF Syndrome iPSCs")
print("Array type : EPIC/450K")
print("Human cell line : iPSCs")
print("Data type : merged array processing")

library(ChAMP)
library(splitstackshape)
library(ggpubr)
library(ggplot2)
library(plyr)
library(gplots)
library(pheatmap)
library(annotatr)
setwd("/media/ankitv/Archivio2/ankit/Array21/Controls")
dir()

#LOAD FILE (method=ChAMP)  ###########
myLoad_2 <-champ.load(directory = getwd(),
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
                     arraytype="450K")

head(myLoad_2$beta)
#QUALITY CONTROL-TIME CONSUMING BE-CARE
QC.GUI(as.matrix(myLoad_2$beta),
       pheno=myLoad_2$pd$Sample_Group,
       arraytype="450K")
summary(myLoad_2$beta)
myNorm4<- champ.norm(myLoad_2$beta,
                     method="BMIQ",
                     arraytype="450K")


dim(myNorm4)
head(myNorm4)

#QUALITY CONTROL-TIME CONSUMING BE-CARE
QC.GUI(as.matrix(myNorm4),
       pheno=myLoad_2$pd$Sample_Group,
       arraytype="450K")

# Do SVD check on data set
champ.SVD(myNorm4,
          rgSet=NULL,
          pd=myLoad_2$pd)


#Batch_effect_detected.svg
summary(myNorm4)
myNorm4 <- data.frame(myNorm4)
myNorm4["TargetID"] <- rownames(myNorm4)

head(myNorm4,3)
head(myLoad_2$beta,3)
myNorm4stack <- stack((myNorm4))
head(myNorm4stack)

myNorm4stack <- data.frame(myNorm4stack)
#https://homepage.divms.uiowa.edu/~luke/classes/STAT4580/histdens.html
#ggplot(myNorm4stack) + geom_density(aes(x = values, color = ind)) + scale_colour_manual(values =c("red","red","red","blue","blue","blue","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))
#save as myNorm4stack.svg





#Load ICF data myNorm
myNorm2Array21 <- read.table("/media/ankitv/Archivio2/ankit/Array21/ICF_hg38/myNorm2.txt", header=T, stringsAsFactors = F)
head(myNorm2Array21)
myNorm2Array21 <- data.frame(myNorm2Array21)
myNorm2Array21["TargetID"] <- rownames(myNorm2Array21)
#Take only shared CpGs among public controls and ICF array data
myNorm2Array21_pub <- merge(myNorm2Array21, myNorm4 , by ="TargetID")
rownames(myNorm2Array21_pub) <- myNorm2Array21_pub$TargetID
myNorm2Array21_pub <- myNorm2Array21_pub[,-1]
head(myNorm2Array21_pub,1)
adjPd <-  read.table("adjPd.txt", header =T, stringsAsFactors = F)

sample_info1 <- adjPd
champ.SVD(myNorm2Array21_pub,
          rgSet=NULL,
          pd=sample_info1)
#save as Batch_effects_SVD_myNorm2Array21_pub.svg
head(myNorm2Array21_pub)
dim(myNorm2Array21_pub)
boxplot(myNorm2Array21_pub, col = c("green","green","green","green","green","green","green","green","blue","blue","red","red","red","orange","orange","orange","orange"))
#boxplot_myNorm2Array21_pub.svg


#Density distribution curve

myNorm2Array21_pubstack <- stack((myNorm2Array21_pub))
head(myNorm2Array21_pubstack)

myNorm2Array21_pubstack <- data.frame(myNorm2Array21_pubstack)
#https://homepage.divms.uiowa.edu/~luke/classes/STAT4580/histdens.html
ggplot(myNorm2Array21_pubstack) + geom_density(aes(x = values, color = ind)) + scale_colour_manual(values =c("green","green","green","green","green","green","green","green","red","red","blue","blue","blue","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))
#save as myNorm2Array21_pubstack.png

colnames(myNorm2Array21_pub)
head(myNorm2Array21_pub)
dim(myNorm2Array21_pub)
myCombat1 <- champ.runCombat(beta=as.matrix(myNorm2Array21_pub),
                             pd=adjPd,
                             variablename="Sample_Group",
                             batchname=c("Batch"),
                             logitTrans=TRUE)



myCombat2 <- myCombat1
boxplot(myCombat2, col = c("green","green","green","green","green","green","green","green","blue","blue","red","red","red","orange","orange","orange","orange"))
#boxplot_myCombat2.png
#myCombat2 <- myCombat2[,c(1:4,7,5,6,8,9:18)]
#QUALITY CONTROL-TIME CONSUMING BE-CARE
QC.GUI(as.matrix(myCombat2),
       pheno=sample_info1,
       arraytype="450K")
champ.SVD(myCombat2,
          rgSet=NULL,
          pd=sample_info1)
#save as Batch_effects_SVD_corrected.svg
#Import MANIFEST
MANIFEST <- read.delim("curhg38Manifest_features.txt", header=T, stringsAsFactors = F)
colnames(MANIFEST)
MANIFEST <- data.frame(MANIFEST)
dim(MANIFEST)
head(MANIFEST)

dim(myCombat2)
head(myCombat2,1)

head(myCombat2)
#Density distribution curve
myCombat2 <- data.frame(myCombat2)
head(myCombat2)
myCombat2stack <- stack((myCombat2))
head(myCombat2stack)

myCombat2stack <- data.frame(myCombat2stack)
#https://homepage.divms.uiowa.edu/~luke/classes/STAT4580/histdens.html
ggplot(myCombat2stack) + geom_density(aes(x = values, color = ind)) + scale_colour_manual(values =c("green","green","green","green","green","green","green","green","red","red","blue","blue","blue","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))
#save as myCombat2stack.png

write.table(myCombat2,"myCombat2.txt",row.names=TRUE,quote=FALSE)
myCombat2 <- read.table("myCombat2.txt",header=TRUE, stringsAsFactors = F)

print('PCA for all probes')
head(myCombat2)
tmyCombat2 = t(myCombat2)
#head(tmyCombat2)
dim(tmyCombat2)
#rownames(tmyCombat2)
tmyCombat2 = data.frame(tmyCombat2)
#head(tmyCombat2)
tmyCombat2["Color"] <-  c("corrClone","corrClone","corrClone","corrClone","Control","Case","Case","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control")
dim(tmyCombat2)
dfxcom <-tmyCombat2[c(1:(length(tmyCombat2)-1))]
#head(dfxcom)
PCcom<-prcomp(dfxcom, center = TRUE, scale. = TRUE)
#head(PC)
PCcomi<-data.frame(PCcom$x,Color=tmyCombat2$Color)
percentagecomAll <- round(PCcom$sdev^2 / sum(PCcom$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomAll <- paste( colnames(PCcomi), "(", paste( as.character(percentagecomAll), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pcom <- ggplot(PCcomi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomAll[1]) + ylab(percentagecomAll[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("red","darkgrey","green","green","darkgrey","red","red","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"))+
  scale_shape_manual(values=c(0,1,2))
pcom <- pcom+theme_classic()
#p + xlim(-50,50)+ ylim(-50,50)
pcom + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA-tmyCombat2.svg", width=11, height=8.6, units="cm", dpi=96)


myCombat3 <- data.frame(myCombat2[order(rownames(myCombat2)),])
head(myCombat3)
#Take median of controls
head(myCombat3)
dim(myCombat3)
myCombat3re <- as.matrix(myCombat3)
head(myCombat3re)


dim(myCombat3re)
#*********************************   Individual samples ************************************
head(myCombat3re,2)
dim(myCombat3re)


myiCombat3reprobed <- cbind.data.frame(rownames(myCombat3re),
                                       myCombat3re[,1:length(colnames(myCombat3re))])

head(myiCombat3reprobed,1)

colnames(myiCombat3reprobed) <- c("TargetID", colnames(myCombat3re))

MERGE_myiCombat3re <- merge(myiCombat3reprobed,MANIFEST,by.x="TargetID",by.y="IlmnID")
head(MERGE_myiCombat3re)
MERGE_myiCombat3re_pos <- MERGE_myiCombat3re[,c(as.numeric(which(colnames(MERGE_myiCombat3re)=="chrhg38")),as.numeric(which(colnames(MERGE_myiCombat3re)=="starthg38")),as.numeric(which(colnames(MERGE_myiCombat3re)=="endhg38")),1:length(myiCombat3reprobed))]
#check for dimensions
dim(MERGE_myiCombat3re_pos)
head(MERGE_myiCombat3re_pos,1)
tail(MERGE_myiCombat3re_pos)
#write.table(MERGE_myiCombat3re_pos, "MERGE_myiCombat3re_pos.txt", col.names = F, quote = F, row.names = F)
head(MERGE_myiCombat3re_pos)
colnames(MERGE_myiCombat3re_pos) <- c("chr","start","end","probeID",colnames(myCombat3re))
dim(MERGE_myiCombat3re_pos)
head(MERGE_myiCombat3re_pos,1)
MERGE_myiCombat3re_pos_chr <- MERGE_myiCombat3re_pos[order(MERGE_myiCombat3re_pos$chr, MERGE_myiCombat3re_pos$start),]
head(MERGE_myiCombat3re_pos_chr)

write.table(MERGE_myiCombat3re_pos_chr, "MERGE_myiCombat3re_pos_chr.txt", col.names = F, quote = F, row.names = F, sep = "\t")


print('Take Median  of Controls')
head(myCombat3re,2)
dim(myCombat3re)
myCombat3pre <- myCombat3re[,c(5,8,9:17,6,7,1:4)]
head(myCombat3pre)
myCombat3reavg <- data.frame(cbind((rowMedians(myCombat3pre[,1:11])),
                                   (myCombat3pre[,1:17])))
head(myCombat3reavg)


colnames(myCombat3reavg) <- c("Allcontrol", colnames(myCombat3pre))
myCombat3reavg["TargetID"] <- rownames(myCombat3reavg)
MERGE_myiCombat3reavg <- merge(myCombat3reavg,MANIFEST,by.x="TargetID",by.y="IlmnID")
head(MERGE_myiCombat3reavg)
dim(MERGE_myiCombat3reavg)

cat("Note: ", length(myCombat2$iPSC_control9) - length(MERGE_myiCombat3reavg$TargetID), " probes could not be assigned to MANIFEST as it is removed while conversion to hg38")

MERGE_myiCombat3reavg_pos <- MERGE_myiCombat3reavg[,c(as.numeric(which(colnames(MERGE_myiCombat3reavg)=="chrhg38")),as.numeric(which(colnames(MERGE_myiCombat3reavg)=="starthg38")),as.numeric(which(colnames(MERGE_myiCombat3reavg)=="endhg38")),1:length(myCombat3reavg))]
#check for dimensions
dim(MERGE_myiCombat3reavg_pos)
head(MERGE_myiCombat3reavg_pos,1)
tail(MERGE_myiCombat3reavg_pos)
write.table(MERGE_myiCombat3reavg_pos, "MERGE_myiCombat3reavg_pos.txt", col.names = F, quote = F, row.names = F)
#MERGE_myiCombat3reavg_pos <- read.table("MERGE_myiCombat3reavg_pos.txt")
head(MERGE_myiCombat3reavg_pos,1)
dim(MERGE_myiCombat3reavg_pos)
head(MERGE_myiCombat3reavg_pos,1)
MERGE_myiCombat3reavg_pos_chr <- MERGE_myiCombat3reavg_pos[order(MERGE_myiCombat3reavg_pos$chrhg38, MERGE_myiCombat3reavg_pos$starthg38),]
head(MERGE_myiCombat3reavg_pos_chr)
dim(MERGE_myiCombat3reavg_pos_chr)
write.table(MERGE_myiCombat3reavg_pos_chr, "MERGE_myiCombat3reavg_pos_chr.txt", append = F,col.names = F, quote = F, row.names = F, sep = "\t")



#Box and Violin Plot Global methylation 
head(myCombat3reavg)
dim(myCombat3reavg)
ipubcomcomdataGlob <- myCombat3reavg
rownames(ipubcomcomdataGlob)
ipubcomcomdataGlob[,1]
ipubcomcomdataGlob1 <- data.frame(ipubcomcomdataGlob)
head(ipubcomcomdataGlob1)
str(ipubcomcomdataGlob1)
ViopubcomGlobIndv <- ipubcomcomdataGlob1
head(ViopubcomGlobIndv)

dim(ViopubcomGlobIndv)
ViopubcomGlobIndv1 <- ViopubcomGlobIndv[,-c(2:12,19)]
head(ViopubcomGlobIndv1)
dim(ViopubcomGlobIndv1)
ViopubcomGlobIndv2 <- stack(ViopubcomGlobIndv1[,c(1,2,5,7,3,4,6)])
colnames(ViopubcomGlobIndv2) <- c("Methylation", "ipubcomdatasets")
head(ViopubcomGlobIndv2)
ggplot(ViopubcomGlobIndv2, aes(x=ipubcomdatasets, y=Methylation, color=ipubcomdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_fill_manual(values=c("white","white","white","white","white","white","white"))+scale_color_manual(values=c("darkgrey","#EC7063","#48C9B0","#52BE80","#922B21","#148F77","#1E8449"))
ggsave("Violin_plot_ViopubcomGlobIndv2_pimethhed_Globalmeth.svg", width=12*1.25, height=10*1.25, units="cm", dpi=96)



#Find genome-wide differentially methylated CpGs
dim(MERGE_myiCombat3re_pos_chr)
head(MERGE_myiCombat3re_pos_chr)

MERGE_myiCombat3re_pos_chr_Betaval <- MERGE_myiCombat3re_pos_chr[,4:21]
rownames(MERGE_myiCombat3re_pos_chr_Betaval) <- MERGE_myiCombat3re_pos_chr_Betaval$probeID
MERGE_myiCombat3re_pos_chr_Betaval <- MERGE_myiCombat3re_pos_chr_Betaval[,-1]
head(MERGE_myiCombat3re_pos_chr_Betaval)
#Convert Beta to M-values
MERGE_myiCombat3re_pos_chr_Mval <- lumi::beta2m(MERGE_myiCombat3re_pos_chr_Betaval)
head(MERGE_myiCombat3re_pos_chr_Mval)
summary(MERGE_myiCombat3re_pos_chr_Mval)
colnames(MERGE_myiCombat3re_pos_chr_Mval) <- paste0(colnames(MERGE_myiCombat3re_pos_chr_Mval), "_Mval")

MERGE_myiCombat3re_pos_chr_Mvalre <- MERGE_myiCombat3re_pos_chr_Mval[,c(5,8,9:17,6,7,1:4)]
head(MERGE_myiCombat3re_pos_chr_Mvalre,3)

MERGE_myiCombat3re_pos_chr_Mvalre <- as.matrix(MERGE_myiCombat3re_pos_chr_Mvalre)
MERGE_myiCombat3re_pos_chr_Mvalavg <- data.frame(cbind(
  (rowMeans(MERGE_myiCombat3re_pos_chr_Mvalre[,1:11])),
  3 * (rowSds(MERGE_myiCombat3re_pos_chr_Mvalre[,1:11])),
  (MERGE_myiCombat3re_pos_chr_Mvalre[,1:17])
  ))

#Check if the functions rowMeans and rowSds works correctly
testrowfval <- data.frame(t(MERGE_myiCombat3re_pos_chr_Mvalre[1:2,1:11]))
colnames(testrowfval) <- c("val1","val2")
mean(c(testrowfval$val1))
testrowfvalsd <- sd(c(testrowfval$val1)) 
testrowfvalsd * 3
mean(c(testrowfval$val2))
testrowfvalsd2 <- sd(c(testrowfval$val2)) 
testrowfvalsd2 * 3

#So both function match with manual check

head(MERGE_myiCombat3re_pos_chr_Mvalavg,2)
colnames(MERGE_myiCombat3re_pos_chr_Mvalavg) <- c("Meancontrol_Mval","Mval_3Sdcontrol",colnames(MERGE_myiCombat3re_pos_chr_Mvalre))
dim(MERGE_myiCombat3re_pos_chr_Mvalavg)
head(MERGE_myiCombat3re_pos_chr_Mvalavg)

#Subtraction, Control -Control
MERGE_myiCombat3re_pos_chr_Mvalavg["Con_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["D250_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$D250_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["UN_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$UN_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control1_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control1_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control2_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control2_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control3_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control3_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control4_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control4_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control5_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control5_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control6_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control6_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control7_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control7_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control8_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control8_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["iPSC_control9_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$iPSC_control9_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)

#Subtraction, Case -Control
MERGE_myiCombat3re_pos_chr_Mvalavg["PG_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$PG_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["PR_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$PR_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["C7_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$C7_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["C13_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$C13_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["C35_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$C35_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["C50_Con_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$C50_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["C7_PR_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$C7_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$PR_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["C13_PG_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$C13_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$PG_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["C35_PR_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$C35_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$PR_Mval)
MERGE_myiCombat3re_pos_chr_Mvalavg["C50_PG_Mval"] <- abs(MERGE_myiCombat3re_pos_chr_Mvalavg$C50_Mval - MERGE_myiCombat3re_pos_chr_Mvalavg$PG_Mval)

MERGE_myiCombat3re_pos_chr_Mvalavg$probeID <- rownames(MERGE_myiCombat3re_pos_chr_Mvalavg)
head(MERGE_myiCombat3re_pos_chr_Mvalavg)
dim(MERGE_myiCombat3re_pos_chr_Mvalavg)

#Process beta value 
MERGE_myipubCombat3re_pos_chr <- MERGE_myiCombat3re_pos_chr
head(MERGE_myipubCombat3re_pos_chr)
dim(MERGE_myipubCombat3re_pos_chr) 

rownames(MERGE_myipubCombat3re_pos_chr) <- paste0(MERGE_myipubCombat3re_pos_chr$chr, 
                                                  "%",
                                                  MERGE_myipubCombat3re_pos_chr$start,
                                                  "%",
                                                  MERGE_myipubCombat3re_pos_chr$end,
                                                  "%",
                                                  MERGE_myipubCombat3re_pos_chr$probeID)

MERGE_myipubCombat3re_pos_chrPre <- MERGE_myipubCombat3re_pos_chr[,-(1:4)]
head(MERGE_myipubCombat3re_pos_chrPre)
dim(MERGE_myipubCombat3re_pos_chrPre)
MERGE_myipubCombat3re_pos_chrPre <- MERGE_myipubCombat3re_pos_chrPre[,c(5,8,9:17,6,7,1:4)]
head(MERGE_myipubCombat3re_pos_chrPre)
dim(MERGE_myipubCombat3re_pos_chrPre)
MERGE_myipubCombat3re_pos_chrPre <- as.matrix(MERGE_myipubCombat3re_pos_chrPre)
MERGE_myipubCombat3re_pos_chravg <- data.frame(cbind((rowMedians(MERGE_myipubCombat3re_pos_chrPre[,1:11])),
                                                     (MERGE_myipubCombat3re_pos_chrPre[,1:17])))

head(MERGE_myipubCombat3re_pos_chravg)
colnames(MERGE_myipubCombat3re_pos_chravg) <- c("Allcontrol",colnames(MERGE_myipubCombat3re_pos_chrPre))
dim(MERGE_myipubCombat3re_pos_chravg)


#Subtraction, Control -Control
MERGE_myipubCombat3re_pos_chravg["Con_Con"] <- MERGE_myipubCombat3re_pos_chravg$Allcontrol - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["D250_Con"] <- MERGE_myipubCombat3re_pos_chravg$D250 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["UN_Con"] <- MERGE_myipubCombat3re_pos_chravg$UN - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control1_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control1 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control2_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control2 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control3_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control3 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control4_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control4 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control5_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control5 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control6_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control6 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control7_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control7 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control8_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control8 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["iPSC_control9_Con"] <- MERGE_myipubCombat3re_pos_chravg$iPSC_control9 - MERGE_myipubCombat3re_pos_chravg$Allcontrol

head(MERGE_myipubCombat3re_pos_chravg)
#Subtraction, Case -Control
MERGE_myipubCombat3re_pos_chravg["PG_Con"] <- MERGE_myipubCombat3re_pos_chravg$PG - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["PR_Con"] <- MERGE_myipubCombat3re_pos_chravg$PR - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["C7_Con"] <- MERGE_myipubCombat3re_pos_chravg$C7 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["C13_Con"] <- MERGE_myipubCombat3re_pos_chravg$C13 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["C35_Con"] <- MERGE_myipubCombat3re_pos_chravg$C35 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["C50_Con"] <- MERGE_myipubCombat3re_pos_chravg$C50 - MERGE_myipubCombat3re_pos_chravg$Allcontrol
MERGE_myipubCombat3re_pos_chravg["C7_PR"] <- MERGE_myipubCombat3re_pos_chravg$C7 - MERGE_myipubCombat3re_pos_chravg$PR
MERGE_myipubCombat3re_pos_chravg["C13_PG"] <- MERGE_myipubCombat3re_pos_chravg$C13 - MERGE_myipubCombat3re_pos_chravg$PG
MERGE_myipubCombat3re_pos_chravg["C35_PR"] <- MERGE_myipubCombat3re_pos_chravg$C35 - MERGE_myipubCombat3re_pos_chravg$PR
MERGE_myipubCombat3re_pos_chravg["C50_PG"] <- MERGE_myipubCombat3re_pos_chravg$C50 - MERGE_myipubCombat3re_pos_chravg$PG

head(MERGE_myipubCombat3re_pos_chravg)
dim(MERGE_myipubCombat3re_pos_chravg)
write.table(MERGE_myipubCombat3re_pos_chravg, "MERGE_myipubCombat3re_pos_chravg.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)
MERGE_myipubCombat3re_pos_chravg1 <- MERGE_myipubCombat3re_pos_chravg
MERGE_myipubCombat3re_pos_chravg1["row"] <- rownames(MERGE_myipubCombat3re_pos_chravg1)
MERGE_myipubCombat3re_pos_chravg1 <- cSplit(MERGE_myipubCombat3re_pos_chravg1, "row", "%")
head(MERGE_myipubCombat3re_pos_chravg1)
dim(MERGE_myipubCombat3re_pos_chravg1)
MERGE_myipubCombat3re_pos_chravg1 <- MERGE_myipubCombat3re_pos_chravg1[,c(41:44,1:40)]
head(MERGE_myipubCombat3re_pos_chravg1)
write.table(MERGE_myipubCombat3re_pos_chravg1, "MERGE_myipubCombat3re_pos_chravg1.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(MERGE_myipubCombat3re_pos_chravg1, "MERGE_myipubCombat3re_pos_chravg2.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)


#Sample specific diff meth CpGs
pMERGE_myipubCombat3re_pos_chravg2 <- merge(MERGE_myipubCombat3re_pos_chravg1, MERGE_myiCombat3re_pos_chr_Mvalavg, by.x="row_4",by.y="probeID", all.x=T) 


dim(pMERGE_myipubCombat3re_pos_chravg2)
head(pMERGE_myipubCombat3re_pos_chravg2)

#Apply the formula for  abs(meanControlMval - sampleMval) >= 3SD

pMERGE_myipubCombat3re_pos_chravg2_PGdev <- pMERGE_myipubCombat3re_pos_chravg2[which(pMERGE_myipubCombat3re_pos_chravg2$PG_Con_Mval >= pMERGE_myipubCombat3re_pos_chravg2$Mval_3Sdcontrol),]
dim(pMERGE_myipubCombat3re_pos_chravg2_PGdev)

#Now filter by beta value difference
#PG specific diff meth CpGs +/- 0.2
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg2_PGdev[which(pMERGE_myipubCombat3re_pos_chravg2_PGdev$PG_Con > 0.2),]
#Should not be Hyper in any controls 
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$D250_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$UN_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control1_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control2_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control3_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control4_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control5_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control6_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control7_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control8_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper$iPSC_control9_Con < 0.2),]
dim(pMERGE_myipubCombat3re_pos_chravg_PG_Hyper)
miMERGE_myipubCombat3re_pos_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PG_Hyper
head(miMERGE_myipubCombat3re_pos_chravg_PG_Hyper)
dim(miMERGE_myipubCombat3re_pos_chravg_PG_Hyper)
miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper <- miMERGE_myipubCombat3re_pos_chravg_PG_Hyper[,c(2:4,1,5:85)]
head(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper)
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper, "miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper$row_4, "miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper.id", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper, "miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper_chr.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)


pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg2_PGdev[which(pMERGE_myipubCombat3re_pos_chravg2_PGdev$PG_Con < -0.2),]
#Should not be Hypo in any controls 
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$D250_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$UN_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control1_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control2_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control3_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control4_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control5_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control6_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control7_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control8_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo$iPSC_control9_Con > -0.2),]
dim(pMERGE_myipubCombat3re_pos_chravg_PG_Hypo)


miMERGE_myipubCombat3re_pos_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PG_Hypo
head(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo)
dim(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo)
miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo <- miMERGE_myipubCombat3re_pos_chravg_PG_Hypo[,c(2:4,1,5:85)]
head(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo, "miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo$row_4, "miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo.id", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo, "miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)

#Apply the formula for  abs(meanControlMval - sampleMval) >= 3SD

pMERGE_myipubCombat3re_pos_chravg2_PRdev <- pMERGE_myipubCombat3re_pos_chravg2[which(pMERGE_myipubCombat3re_pos_chravg2$PR_Con_Mval >= pMERGE_myipubCombat3re_pos_chravg2$Mval_3Sdcontrol),]
dim(pMERGE_myipubCombat3re_pos_chravg2_PRdev)

#Now filter by difference of beta values
#PR specific diff meth CpGs +/- 0.2
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg2_PRdev[which(pMERGE_myipubCombat3re_pos_chravg2_PRdev$PR_Con > 0.2),]
#Should not be Hyper in any controls 
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$D250_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$UN_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control1_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control2_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control3_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control4_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control5_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control6_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control7_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control8_Con < 0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper$iPSC_control9_Con < 0.2),]
dim(pMERGE_myipubCombat3re_pos_chravg_PR_Hyper)
miMERGE_myipubCombat3re_pos_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos_chravg_PR_Hyper
head(miMERGE_myipubCombat3re_pos_chravg_PR_Hyper)
dim(miMERGE_myipubCombat3re_pos_chravg_PR_Hyper)
miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper <- miMERGE_myipubCombat3re_pos_chravg_PR_Hyper[,c(2:4,1,5:85)]
head(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper, "miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper$row_4, "miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper.id", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper, "miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper_chr.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)


pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg2_PRdev[which(pMERGE_myipubCombat3re_pos_chravg2_PRdev$PR_Con < -0.2),]
#Should not be Hypo in any controls 
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$D250_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$UN_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control1_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control2_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control3_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control4_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control5_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control6_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control7_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control8_Con > -0.2),]
pMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo$iPSC_control9_Con > -0.2),]
dim(pMERGE_myipubCombat3re_pos_chravg_PR_Hypo)


miMERGE_myipubCombat3re_pos_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos_chravg_PR_Hypo
head(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo)
dim(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo)
miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo <- miMERGE_myipubCombat3re_pos_chravg_PR_Hypo[,c(2:4,1,5:85)]
head(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo, "miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo$row_4, "miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo.id", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo, "miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)

miMERGE_myipubCombat3re_chrpos_avg_PG <- rbind.data.frame(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo, miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper)
miMERGE_myipubCombat3re_chrpos_avg_PR <- rbind.data.frame(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo, miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper)
dim(miMERGE_myipubCombat3re_chrpos_avg_PG)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR)
#miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt")
#miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper_chr.txt")
#miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt")
#miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper_chr.txt")
write.table(miMERGE_myipubCombat3re_chrpos_avg_PG, "miMERGE_myipubCombat3re_chrpos_avg_PG.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(miMERGE_myipubCombat3re_chrpos_avg_PR, "miMERGE_myipubCombat3re_chrpos_avg_PR.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
#Just in case same name order required but the objects assigned below are same
miMERGE_myipubCombat3re_pos_chravg_PG = miMERGE_myipubCombat3re_chrpos_avg_PG
miMERGE_myipubCombat3re_pos_chravg_PR = miMERGE_myipubCombat3re_chrpos_avg_PR

detach("package:dplyr")
#Peek into genome wide view
system("bedtools makewindows -g /home/ankitv/ref_av/hg38/hg38.chrom.sizes -w 500 > hg38_500.bed")

print('--------------------  Identify differentially methylated 500 bp bins ---------------------------')
system("bedtools intersect -wa -wb -a MERGE_myiCombat3re_pos_chr.txt -b hg38_500.bed > MERGE_myiCombat3re_pos500bp_chr.txt")
MERGE_myiCombat3re_pos500bp_chr <- read.table("MERGE_myiCombat3re_pos500bp_chr.txt", header = F, stringsAsFactors = F)
head(MERGE_myiCombat3re_pos500bp_chr)
colnames(MERGE_myiCombat3re_pos500bp_chr) <- c(colnames(MERGE_myiCombat3re_pos_chr), "Binchr","Binstart","Binend")

dim(MERGE_myiCombat3re_pos500bp_chr)
head(MERGE_myiCombat3re_pos500bp_chr)

MERGE_myiCombat3re_pos500bp_chr["Bin500bp"] <- paste0(MERGE_myiCombat3re_pos500bp_chr$Binchr, "%",
                                                      MERGE_myiCombat3re_pos500bp_chr$Binstart, "%",
                                                      MERGE_myiCombat3re_pos500bp_chr$Binend)

head(MERGE_myiCombat3re_pos500bp_chr[,c(25,5:21)])
dim(MERGE_myiCombat3re_pos500bp_chr)
MERGE_myiCombat3re_pos500bp_chr_rearranged <- MERGE_myiCombat3re_pos500bp_chr[,c(25,5:21)]
head(MERGE_myiCombat3re_pos500bp_chr_rearranged)
ipubcomdataBin500bp_filt <- MERGE_myiCombat3re_pos500bp_chr_rearranged
rownames(ipubcomdataBin500bp_filt)
colnames(ipubcomdataBin500bp_filt)
head(ipubcomdataBin500bp_filt)
dim(ipubcomdataBin500bp_filt)
Count_ipubcomdataBin500bp_filt1 <- count(ipubcomdataBin500bp_filt, "Bin500bp")
head(Count_ipubcomdataBin500bp_filt1)
dim(Count_ipubcomdataBin500bp_filt1)
Total_CpGs_ineach_500bp <- Count_ipubcomdataBin500bp_filt1
sfiltpubcomaggregateBin500bp = aggregate(ipubcomdataBin500bp_filt,by=list(ipubcomdataBin500bp_filt$Bin500bp), mean)
head(sfiltpubcomaggregateBin500bp, 2)
ipubcomdataBin500bp_filt1=sfiltpubcomaggregateBin500bp
rownames(ipubcomdataBin500bp_filt1)
ipubcomdataBin500bp_filt1[,1]
rownames(ipubcomdataBin500bp_filt1)=ipubcomdataBin500bp_filt1[,1]
rownames(ipubcomdataBin500bp_filt1)
colnames(ipubcomdataBin500bp_filt1)
ipubcomdataBin500bp_filt1 = ipubcomdataBin500bp_filt1[,-c(1:2)]
head(ipubcomdataBin500bp_filt1)
dim(ipubcomdataBin500bp_filt1)
ipubcomdataBin500bp_filt1_counted <- cbind(ipubcomdataBin500bp_filt1, Count_ipubcomdataBin500bp_filt1)
head(ipubcomdataBin500bp_filt1_counted)
dim(ipubcomdataBin500bp_filt1_counted)
tail(ipubcomdataBin500bp_filt1_counted)
ipubcomdataBin500bp_filt1_counted1 <- ipubcomdataBin500bp_filt1_counted[which(ipubcomdataBin500bp_filt1_counted$freq >= 3),]
head(ipubcomdataBin500bp_filt1_counted1)
dim(ipubcomdataBin500bp_filt1_counted1)
write.table(ipubcomdataBin500bp_filt1_counted1, "ipubcomdataBin500bp_filt1_counted1.filtmin.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE, col.names = T)
#ipubcomdataBin500bp_filt1_counted1 <- read.table("ipubcomdataBin500bp_filt1_counted1.filtmin.txt", header = T, stringsAsFactors = F)

MERGE_myiCombat3re_pos500bp_chr_Betaval <- ipubcomdataBin500bp_filt1_counted1[,1:17]
head(MERGE_myiCombat3re_pos500bp_chr_Betaval)
#Convert Beta to M-values
MERGE_myiCombat3re_pos500bp_chr_Mval <- lumi::beta2m(MERGE_myiCombat3re_pos500bp_chr_Betaval)
head(MERGE_myiCombat3re_pos500bp_chr_Mval)
summary(MERGE_myiCombat3re_pos500bp_chr_Mval)
colnames(MERGE_myiCombat3re_pos500bp_chr_Mval) <- paste0(colnames(MERGE_myiCombat3re_pos500bp_chr_Mval), "_Mval")

MERGE_myiCombat3re_pos500bp_chr_Mvalre <- MERGE_myiCombat3re_pos500bp_chr_Mval[,c(5,8,9:17,6,7,1:4)]
head(MERGE_myiCombat3re_pos500bp_chr_Mvalre,3)

MERGE_myiCombat3re_pos500bp_chr_Mvalre <- as.matrix(MERGE_myiCombat3re_pos500bp_chr_Mvalre)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg <- data.frame(cbind(
  (rowMeans(MERGE_myiCombat3re_pos500bp_chr_Mvalre[,1:11])),
  3 * (rowSds(MERGE_myiCombat3re_pos500bp_chr_Mvalre[,1:11])),
  (MERGE_myiCombat3re_pos500bp_chr_Mvalre[,1:17])
))

#Check if the functions rowMeans and rowSds works correctly
test500rowfval <- data.frame(t(MERGE_myiCombat3re_pos500bp_chr_Mvalre[1:2,1:11]))
colnames(test500rowfval) <- c("val1","val2")
mean(c(test500rowfval$val1))
test500rowfvalsd <- sd(c(test500rowfval$val1)) 
test500rowfvalsd * 3
mean(c(test500rowfval$val2))
test500rowfvalsd2 <- sd(c(test500rowfval$val2)) 
test500rowfvalsd2 * 3

#So both function match with manual check

head(MERGE_myiCombat3re_pos500bp_chr_Mvalavg,2)
colnames(MERGE_myiCombat3re_pos500bp_chr_Mvalavg) <- c("Meancontrol_Mval","Mval_3Sdcontrol",colnames(MERGE_myiCombat3re_pos500bp_chr_Mvalre))
dim(MERGE_myiCombat3re_pos500bp_chr_Mvalavg)
head(MERGE_myiCombat3re_pos500bp_chr_Mvalavg)

#Subtraction, Control -Control
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["Con_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["D250_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$D250_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["UN_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$UN_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control1_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control1_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control2_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control2_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control3_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control3_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control4_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control4_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control5_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control5_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control6_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control6_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control7_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control7_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control8_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control8_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["iPSC_control9_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$iPSC_control9_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)

#Subtraction, Case -Control
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["PG_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$PG_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["PR_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$PR_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["C7_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$C7_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["C13_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$C13_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["C35_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$C35_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["C50_Con_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$C50_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$Meancontrol_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["C7_PR_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$C7_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$PR_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["C13_PG_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$C13_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$PG_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["C35_PR_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$C35_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$PR_Mval)
MERGE_myiCombat3re_pos500bp_chr_Mvalavg["C50_PG_Mval"] <- abs(MERGE_myiCombat3re_pos500bp_chr_Mvalavg$C50_Mval - MERGE_myiCombat3re_pos500bp_chr_Mvalavg$PG_Mval)

MERGE_myiCombat3re_pos500bp_chr_Mvalavg$probeID <- rownames(MERGE_myiCombat3re_pos500bp_chr_Mvalavg)
head(MERGE_myiCombat3re_pos500bp_chr_Mvalavg)
dim(MERGE_myiCombat3re_pos500bp_chr_Mvalavg)

#Process beta value 
MERGE_myipubCombat3re_pos500bp_chr <- MERGE_myiCombat3re_pos500bp_chr_Betaval
head(MERGE_myipubCombat3re_pos500bp_chr)
dim(MERGE_myipubCombat3re_pos500bp_chr) 

MERGE_myipubCombat3re_pos500bp_chrPre <- MERGE_myipubCombat3re_pos500bp_chr[,c(5,8,9:17,6,7,1:4)]
head(MERGE_myipubCombat3re_pos500bp_chrPre)
dim(MERGE_myipubCombat3re_pos500bp_chrPre)
MERGE_myipubCombat3re_pos500bp_chrPre <- as.matrix(MERGE_myipubCombat3re_pos500bp_chrPre)
MERGE_myipubCombat3re_pos500bp_chravg <- data.frame(cbind((rowMedians(MERGE_myipubCombat3re_pos500bp_chrPre[,1:11])),
                                                          (MERGE_myipubCombat3re_pos500bp_chrPre[,1:17])))

head(MERGE_myipubCombat3re_pos500bp_chravg)
colnames(MERGE_myipubCombat3re_pos500bp_chravg) <- c("Allcontrol",colnames(MERGE_myipubCombat3re_pos500bp_chrPre))
dim(MERGE_myipubCombat3re_pos500bp_chravg)


#Subtraction, Control -Control
MERGE_myipubCombat3re_pos500bp_chravg["Con_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["D250_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$D250 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["UN_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$UN - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control1_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control1 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control2_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control2 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control3_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control3 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control4_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control4 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control5_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control5 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control6_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control6 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control7_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control7 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control8_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control8 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["iPSC_control9_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$iPSC_control9 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol

head(MERGE_myipubCombat3re_pos500bp_chravg)
#Subtraction, Case -Control
MERGE_myipubCombat3re_pos500bp_chravg["PG_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$PG - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["PR_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$PR - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["C7_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$C7 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["C13_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$C13 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["C35_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$C35 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["C50_Con"] <- MERGE_myipubCombat3re_pos500bp_chravg$C50 - MERGE_myipubCombat3re_pos500bp_chravg$Allcontrol
MERGE_myipubCombat3re_pos500bp_chravg["C7_PR"] <- MERGE_myipubCombat3re_pos500bp_chravg$C7 - MERGE_myipubCombat3re_pos500bp_chravg$PR
MERGE_myipubCombat3re_pos500bp_chravg["C13_PG"] <- MERGE_myipubCombat3re_pos500bp_chravg$C13 - MERGE_myipubCombat3re_pos500bp_chravg$PG
MERGE_myipubCombat3re_pos500bp_chravg["C35_PR"] <- MERGE_myipubCombat3re_pos500bp_chravg$C35 - MERGE_myipubCombat3re_pos500bp_chravg$PR
MERGE_myipubCombat3re_pos500bp_chravg["C50_PG"] <- MERGE_myipubCombat3re_pos500bp_chravg$C50 - MERGE_myipubCombat3re_pos500bp_chravg$PG

head(MERGE_myipubCombat3re_pos500bp_chravg)
dim(MERGE_myipubCombat3re_pos500bp_chravg)
write.table(MERGE_myipubCombat3re_pos500bp_chravg, "MERGE_myipubCombat3re_pos500bp_chravg.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)
MERGE_myipubCombat3re_pos500bp_chravg1 <- MERGE_myipubCombat3re_pos500bp_chravg
MERGE_myipubCombat3re_pos500bp_chravg1["row"] <- rownames(MERGE_myipubCombat3re_pos500bp_chravg1)


#Sample specific diff meth CpGs
pMERGE_myipubCombat3re_pos500bp_chravg2 <- merge(MERGE_myipubCombat3re_pos500bp_chravg1, MERGE_myiCombat3re_pos500bp_chr_Mvalavg, by.x="row",by.y="probeID", all.x=T) 


dim(pMERGE_myipubCombat3re_pos500bp_chravg2)
head(pMERGE_myipubCombat3re_pos500bp_chravg2)

#Apply the formula for  abs(meanControlMval - sampleMval) >= 3SD

pMERGE_myipubCombat3re_pos500bp_chravg2_PGdev <- pMERGE_myipubCombat3re_pos500bp_chravg2[which(pMERGE_myipubCombat3re_pos500bp_chravg2$PG_Con_Mval >= pMERGE_myipubCombat3re_pos500bp_chravg2$Mval_3Sdcontrol),]
dim(pMERGE_myipubCombat3re_pos500bp_chravg2_PGdev)

#Now filter by beta value difference
#PG specific diff meth CpGs +/- 0.2
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg2_PGdev[which(pMERGE_myipubCombat3re_pos500bp_chravg2_PGdev$PG_Con > 0.2),]
#Should not be Hyper in any controls 
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$D250_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$UN_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control1_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control2_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control3_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control4_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control5_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control6_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control7_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control8_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper$iPSC_control9_Con < 0.2),]
dim(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper)
miMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper1 <- cSplit(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper, "row", "%")
head(miMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper1)
dim(miMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper1)
miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper <- miMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper1[,c(82:84,1:81)]
head(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper)
dim(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper, "miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper, "miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper_chr.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)


pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg2_PGdev[which(pMERGE_myipubCombat3re_pos500bp_chravg2_PGdev$PG_Con < -0.2),]
#Should not be Hypo in any controls 
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$D250_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$UN_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control1_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control2_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control3_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control4_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control5_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control6_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control7_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control8_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo$iPSC_control9_Con > -0.2),]
dim(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo)


miMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo1 <- cSplit(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo, "row", "%")
head(miMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo1)
dim(miMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo1)
miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo <- miMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo1[,c(82:84,1:81)]
head(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo)
dim(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo, "miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo, "miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo_chr.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)

#Apply the formula for  abs(meanControlMval - sampleMval) >= 3SD

pMERGE_myipubCombat3re_pos500bp_chravg2_PRdev <- pMERGE_myipubCombat3re_pos500bp_chravg2[which(pMERGE_myipubCombat3re_pos500bp_chravg2$PR_Con_Mval >= pMERGE_myipubCombat3re_pos500bp_chravg2$Mval_3Sdcontrol),]
dim(pMERGE_myipubCombat3re_pos500bp_chravg2_PRdev)

#Now filter by difference of beta values
#PR specific diff meth CpGs +/- 0.2
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg2_PRdev[which(pMERGE_myipubCombat3re_pos500bp_chravg2_PRdev$PR_Con > 0.2),]
#Should not be Hyper in any controls 
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$D250_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$UN_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control1_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control2_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control3_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control4_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control5_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control6_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control7_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control8_Con < 0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper$iPSC_control9_Con < 0.2),]
dim(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper)
miMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper1 <- cSplit(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper, "row", "%")
head(miMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper1)
dim(miMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper1)
miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper <- miMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper1[,c(82:84,1:81)]
head(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper)
dim(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper, "miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper, "miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper_chr.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)


pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg2_PRdev[which(pMERGE_myipubCombat3re_pos500bp_chravg2_PRdev$PR_Con < -0.2),]
#Should not be Hypo in any controls 
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$D250_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$UN_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control1_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control2_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control3_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control4_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control5_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control6_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control7_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control8_Con > -0.2),]
pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo <- pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo[which(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo$iPSC_control9_Con > -0.2),]
dim(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo)


miMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo1 <- cSplit(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo, "row", "%")
head(miMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo1)
dim(miMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo1)
miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo <- miMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo1[,c(82:84,1:81)]
head(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo)
dim(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo, "miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo.txt", sep = "\t", quote = F, append = F, row.names = T, col.names = T)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo, "miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo_chr.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)

miMERGE_myipubCombat3re_chrpos500bp_avg_PG <- rbind.data.frame(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo, miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper)
miMERGE_myipubCombat3re_chrpos500bp_avg_PR <- rbind.data.frame(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo, miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper)
dim(miMERGE_myipubCombat3re_chrpos500bp_avg_PG)
dim(miMERGE_myipubCombat3re_chrpos500bp_avg_PR)
#miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo_chr.txt")
#miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper_chr.txt")
#miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo_chr.txt")
#miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper_chr.txt")
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PG, "miMERGE_myipubCombat3re_chrpos500bp_avg_PG.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
write.table(miMERGE_myipubCombat3re_chrpos500bp_avg_PR, "miMERGE_myipubCombat3re_chrpos500bp_avg_PR.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)

#Bin500bp_Hyper PG
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper_chr.txt -b miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper_chr.txt > MERGE_myipubCombat3re_PG_Bin500bp_Hyper.txt")
MERGE_myipubCombat3re_PG_Bin500bp_Hyper <- read.table("MERGE_myipubCombat3re_PG_Bin500bp_Hyper.txt", header = F, stringsAsFactors = F)
head(MERGE_myipubCombat3re_PG_Bin500bp_Hyper)
colnames(MERGE_myipubCombat3re_PG_Bin500bp_Hyper) <- c(colnames(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hyper), colnames(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper))
MERGE_myipubCombat3re_PG_Bin500bp_Hyper["Bin500bp_Hyper"] <- paste0(MERGE_myipubCombat3re_PG_Bin500bp_Hyper$row_1, "%",
                                                                    MERGE_myipubCombat3re_PG_Bin500bp_Hyper$row_2, "%",
                                                                    MERGE_myipubCombat3re_PG_Bin500bp_Hyper$row_3)

Count_MERGE_myipubCombat3re_PG_Bin500bp_Hyper <- count(MERGE_myipubCombat3re_PG_Bin500bp_Hyper, "Bin500bp_Hyper")
head(Count_MERGE_myipubCombat3re_PG_Bin500bp_Hyper)
dim(Count_MERGE_myipubCombat3re_PG_Bin500bp_Hyper)

Count_ipubcomdataBin500bp_Hyper_PG_1_filtmin2 <- Count_MERGE_myipubCombat3re_PG_Bin500bp_Hyper[which(Count_MERGE_myipubCombat3re_PG_Bin500bp_Hyper$freq >=2),]
head(Count_ipubcomdataBin500bp_Hyper_PG_1_filtmin2)
dim(Count_ipubcomdataBin500bp_Hyper_PG_1_filtmin2)

dim(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper)

ipubcomdataBin500bp_Hyper_PG_1_filtmin2<- merge(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hyper, Count_ipubcomdataBin500bp_Hyper_PG_1_filtmin2, by.x="row",by.y = "Bin500bp_Hyper", all.y = TRUE ) 
head(ipubcomdataBin500bp_Hyper_PG_1_filtmin2)
dim(ipubcomdataBin500bp_Hyper_PG_1_filtmin2)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_freq <- ddply(MERGE_myipubCombat3re_PG_Bin500bp_Hyper, .(Bin500bp_Hyper), summarize, Probes = toString(row_4))
dim(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_freq)
head(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_freq)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_freq$Probes <- gsub(" ", "", ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_freq$Probes)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2re <- merge(ipubcomdataBin500bp_Hyper_PG_1_filtmin2, ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_freq, by.x="row", by.y="Bin500bp_Hyper", all.x=TRUE)
rownames(ipubcomdataBin500bp_Hyper_PG_1_filtmin2re) <- ipubcomdataBin500bp_Hyper_PG_1_filtmin2re$row
ipubcomdataBin500bp_Hyper_PG_1_filtmin2re1 <- merge(ipubcomdataBin500bp_Hyper_PG_1_filtmin2re, Total_CpGs_ineach_500bp,by.x="row", by.y="Bin500bp", all.x=TRUE)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2re1 <- cSplit(ipubcomdataBin500bp_Hyper_PG_1_filtmin2re1, "row", "%")
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr <- ipubcomdataBin500bp_Hyper_PG_1_filtmin2re1[,c(85:87,1:81,84,82,83)]
head(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr)
dim(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr)
write.table(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr, "ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
system("sort -k1,1 -k2,2n ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr.txt > ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr.sort.txt")
system("bedtools closest -wa -wb -a ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38.txt")
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38 <- read.table("ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38.txt", header = F, stringsAsFactors = F)
head(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38)
colnames(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38) <- c(colnames(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr), "GeneChr", "GeneStart", "GeneEnd", "Type","ENSID","Gene","Distance")
head(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38)
dim(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38 <- data.frame(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38["Bin500bp_Hyper"] <- paste0(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38$row_1,"%", ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38$row_2, "%",ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38$row_3)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_genedis <- ddply(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_genehg38, .(Bin500bp_Hyper), summarize, Gene = toString(Gene), Distance = toString(Distance))
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_genedis$Gene <- gsub(" ", "", ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_genedis$Gene)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_genedis$Distance <- gsub(" ", "", ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_genedis$Distance)
dim(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_genedis)

ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_bin <- data.frame(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_bin["Bin500bp_Hyper"] <- paste0(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_bin$row_1,"%", ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_bin$row_2, "%",ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_bin$row_3)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis <- merge(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_bin,ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_colapprob_genedis,by.x="Bin500bp_Hyper", by.y="Bin500bp_Hyper", all.x=TRUE)
head(ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis)
dim(ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis)
ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis <- ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis[,c(2:4,6:22,86:90)]
ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis["Meth_level"] <- "Hyper"
write.table(ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis, "ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = T)

#Bin500bp_Hyper PR
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper_chr.txt -b miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper_chr.txt > MERGE_myipubCombat3re_PR_Bin500bp_Hyper.txt")
MERGE_myipubCombat3re_PR_Bin500bp_Hyper <- read.table("MERGE_myipubCombat3re_PR_Bin500bp_Hyper.txt", header = F, stringsAsFactors = F)
head(MERGE_myipubCombat3re_PR_Bin500bp_Hyper)
colnames(MERGE_myipubCombat3re_PR_Bin500bp_Hyper) <- c(colnames(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hyper), colnames(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper))
MERGE_myipubCombat3re_PR_Bin500bp_Hyper["Bin500bp_Hyper"] <- paste0(MERGE_myipubCombat3re_PR_Bin500bp_Hyper$row_1, "%",
                                                                    MERGE_myipubCombat3re_PR_Bin500bp_Hyper$row_2, "%",
                                                                    MERGE_myipubCombat3re_PR_Bin500bp_Hyper$row_3)

Count_MERGE_myipubCombat3re_PR_Bin500bp_Hyper <- count(MERGE_myipubCombat3re_PR_Bin500bp_Hyper, "Bin500bp_Hyper")
head(Count_MERGE_myipubCombat3re_PR_Bin500bp_Hyper)
dim(Count_MERGE_myipubCombat3re_PR_Bin500bp_Hyper)

Count_ipubcomdataBin500bp_Hyper_PR_1_filtmin2 <- Count_MERGE_myipubCombat3re_PR_Bin500bp_Hyper[which(Count_MERGE_myipubCombat3re_PR_Bin500bp_Hyper$freq >=2),]
head(Count_ipubcomdataBin500bp_Hyper_PR_1_filtmin2)
dim(Count_ipubcomdataBin500bp_Hyper_PR_1_filtmin2)

dim(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper)

ipubcomdataBin500bp_Hyper_PR_1_filtmin2<- merge(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hyper, Count_ipubcomdataBin500bp_Hyper_PR_1_filtmin2, by.x="row",by.y = "Bin500bp_Hyper", all.y = TRUE ) 
head(ipubcomdataBin500bp_Hyper_PR_1_filtmin2)
dim(ipubcomdataBin500bp_Hyper_PR_1_filtmin2)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_freq <- ddply(MERGE_myipubCombat3re_PR_Bin500bp_Hyper, .(Bin500bp_Hyper), summarize, Probes = toString(row_4))
dim(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_freq)
head(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_freq)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_freq$Probes <- gsub(" ", "", ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_freq$Probes)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2re <- merge(ipubcomdataBin500bp_Hyper_PR_1_filtmin2, ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_freq, by.x="row", by.y="Bin500bp_Hyper", all.x=TRUE)
rownames(ipubcomdataBin500bp_Hyper_PR_1_filtmin2re) <- ipubcomdataBin500bp_Hyper_PR_1_filtmin2re$row
ipubcomdataBin500bp_Hyper_PR_1_filtmin2re1 <- merge(ipubcomdataBin500bp_Hyper_PR_1_filtmin2re, Total_CpGs_ineach_500bp,by.x="row", by.y="Bin500bp", all.x=TRUE)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2re1 <- cSplit(ipubcomdataBin500bp_Hyper_PR_1_filtmin2re1, "row", "%")
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr <- ipubcomdataBin500bp_Hyper_PR_1_filtmin2re1[,c(85:87,1:81,84,82,83)]
head(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr)
dim(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr)
write.table(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr, "ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
system("sort -k1,1 -k2,2n ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr.txt > ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr.sort.txt")
system("bedtools closest -wa -wb -a ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38.txt")
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38 <- read.table("ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38.txt", header = F, stringsAsFactors = F)
head(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38)
colnames(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38) <- c(colnames(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr), "GeneChr", "GeneStart", "GeneEnd", "Type","ENSID","Gene","Distance")
head(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38)
dim(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38 <- data.frame(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38["Bin500bp_Hyper"] <- paste0(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38$row_1,"%", ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38$row_2, "%",ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38$row_3)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_genedis <- ddply(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_genehg38, .(Bin500bp_Hyper), summarize, Gene = toString(Gene), Distance = toString(Distance))
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_genedis$Gene <- gsub(" ", "", ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_genedis$Gene)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_genedis$Distance <- gsub(" ", "", ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_genedis$Distance)
dim(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_genedis)

ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_bin <- data.frame(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_bin["Bin500bp_Hyper"] <- paste0(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_bin$row_1,"%", ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_bin$row_2, "%",ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_bin$row_3)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis <- merge(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_bin,ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_colapprob_genedis,by.x="Bin500bp_Hyper", by.y="Bin500bp_Hyper", all.x=TRUE)
head(ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis)
dim(ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis)
ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis <- ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis[,c(2:4,6:22,86:90)]
ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis["Meth_level"] <- "Hyper"
write.table(ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis, "ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = T)


#Bin500bp_Hypo PG
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo_chr.txt -b miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt > MERGE_myipubCombat3re_PG_Bin500bp_Hypo.txt")
MERGE_myipubCombat3re_PG_Bin500bp_Hypo <- read.table("MERGE_myipubCombat3re_PG_Bin500bp_Hypo.txt", header = F, stringsAsFactors = F)
head(MERGE_myipubCombat3re_PG_Bin500bp_Hypo)
colnames(MERGE_myipubCombat3re_PG_Bin500bp_Hypo) <- c(colnames(miMERGE_myipubCombat3re_chrpos500bp_avg_PG_Hypo), colnames(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo))
MERGE_myipubCombat3re_PG_Bin500bp_Hypo["Bin500bp_Hypo"] <- paste0(MERGE_myipubCombat3re_PG_Bin500bp_Hypo$row_1, "%",
                                                                  MERGE_myipubCombat3re_PG_Bin500bp_Hypo$row_2, "%",
                                                                  MERGE_myipubCombat3re_PG_Bin500bp_Hypo$row_3)

Count_MERGE_myipubCombat3re_PG_Bin500bp_Hypo <- count(MERGE_myipubCombat3re_PG_Bin500bp_Hypo, "Bin500bp_Hypo")
head(Count_MERGE_myipubCombat3re_PG_Bin500bp_Hypo)
dim(Count_MERGE_myipubCombat3re_PG_Bin500bp_Hypo)

Count_ipubcomdataBin500bp_Hypo_PG_1_filtmin2 <- Count_MERGE_myipubCombat3re_PG_Bin500bp_Hypo[which(Count_MERGE_myipubCombat3re_PG_Bin500bp_Hypo$freq >=2),]
head(Count_ipubcomdataBin500bp_Hypo_PG_1_filtmin2)
dim(Count_ipubcomdataBin500bp_Hypo_PG_1_filtmin2)

dim(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo)

ipubcomdataBin500bp_Hypo_PG_1_filtmin2<- merge(pMERGE_myipubCombat3re_pos500bp_chravg_PG_Hypo, Count_ipubcomdataBin500bp_Hypo_PG_1_filtmin2, by.x="row",by.y = "Bin500bp_Hypo", all.y = TRUE ) 
head(ipubcomdataBin500bp_Hypo_PG_1_filtmin2)
dim(ipubcomdataBin500bp_Hypo_PG_1_filtmin2)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_freq <- ddply(MERGE_myipubCombat3re_PG_Bin500bp_Hypo, .(Bin500bp_Hypo), summarize, Probes = toString(row_4))
dim(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_freq)
head(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_freq)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_freq$Probes <- gsub(" ", "", ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_freq$Probes)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2re <- merge(ipubcomdataBin500bp_Hypo_PG_1_filtmin2, ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_freq, by.x="row", by.y="Bin500bp_Hypo", all.x=TRUE)
rownames(ipubcomdataBin500bp_Hypo_PG_1_filtmin2re) <- ipubcomdataBin500bp_Hypo_PG_1_filtmin2re$row
ipubcomdataBin500bp_Hypo_PG_1_filtmin2re1 <- merge(ipubcomdataBin500bp_Hypo_PG_1_filtmin2re, Total_CpGs_ineach_500bp,by.x="row", by.y="Bin500bp", all.x=TRUE)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2re1 <- cSplit(ipubcomdataBin500bp_Hypo_PG_1_filtmin2re1, "row", "%")
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr <- ipubcomdataBin500bp_Hypo_PG_1_filtmin2re1[,c(85:87,1:81,84,82,83)]
head(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr)
dim(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr)
write.table(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr, "ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
system("sort -k1,1 -k2,2n ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr.txt > ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr.sort.txt")
system("bedtools closest -wa -wb -a ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38.txt")

ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38 <- read.table("ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38.txt", header = F, stringsAsFactors = F)
head(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38)
colnames(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38) <- c(colnames(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr), "GeneChr", "GeneStart", "GeneEnd", "Type","ENSID","Gene","Distance")
head(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38)
dim(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38 <- data.frame(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38["Bin500bp_Hypo"] <- paste0(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38$row_1,"%", ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38$row_2, "%",ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38$row_3)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_genedis <- ddply(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_genehg38, .(Bin500bp_Hypo), summarize, Gene = toString(Gene), Distance = toString(Distance))
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_genedis$Gene <- gsub(" ", "", ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_genedis$Gene)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_genedis$Distance <- gsub(" ", "", ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_genedis$Distance)
dim(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_genedis)

ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_bin <- data.frame(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_bin["Bin500bp_Hypo"] <- paste0(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_bin$row_1,"%", ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_bin$row_2, "%",ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_bin$row_3)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis <- merge(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_bin,ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_colapprob_genedis,by.x="Bin500bp_Hypo", by.y="Bin500bp_Hypo", all.x=TRUE)
head(ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis)
dim(ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis)
ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis <- ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis[,c(2:4,6:22,86:90)]
ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis["Meth_level"] <- "Hypo"
write.table(ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis, "ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = T)

#Bin500bp_Hypo PR
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo_chr.txt -b miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt > MERGE_myipubCombat3re_PR_Bin500bp_Hypo.txt")
MERGE_myipubCombat3re_PR_Bin500bp_Hypo <- read.table("MERGE_myipubCombat3re_PR_Bin500bp_Hypo.txt", header = F, stringsAsFactors = F)
head(MERGE_myipubCombat3re_PR_Bin500bp_Hypo)
colnames(MERGE_myipubCombat3re_PR_Bin500bp_Hypo) <- c(colnames(miMERGE_myipubCombat3re_chrpos500bp_avg_PR_Hypo), colnames(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo))
MERGE_myipubCombat3re_PR_Bin500bp_Hypo["Bin500bp_Hypo"] <- paste0(MERGE_myipubCombat3re_PR_Bin500bp_Hypo$row_1, "%",
                                                                  MERGE_myipubCombat3re_PR_Bin500bp_Hypo$row_2, "%",
                                                                  MERGE_myipubCombat3re_PR_Bin500bp_Hypo$row_3)

Count_MERGE_myipubCombat3re_PR_Bin500bp_Hypo <- count(MERGE_myipubCombat3re_PR_Bin500bp_Hypo, "Bin500bp_Hypo")
head(Count_MERGE_myipubCombat3re_PR_Bin500bp_Hypo)
dim(Count_MERGE_myipubCombat3re_PR_Bin500bp_Hypo)

Count_ipubcomdataBin500bp_Hypo_PR_1_filtmin2 <- Count_MERGE_myipubCombat3re_PR_Bin500bp_Hypo[which(Count_MERGE_myipubCombat3re_PR_Bin500bp_Hypo$freq >=2),]
head(Count_ipubcomdataBin500bp_Hypo_PR_1_filtmin2)
dim(Count_ipubcomdataBin500bp_Hypo_PR_1_filtmin2)

dim(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo)

ipubcomdataBin500bp_Hypo_PR_1_filtmin2<- merge(pMERGE_myipubCombat3re_pos500bp_chravg_PR_Hypo, Count_ipubcomdataBin500bp_Hypo_PR_1_filtmin2, by.x="row",by.y = "Bin500bp_Hypo", all.y = TRUE ) 
head(ipubcomdataBin500bp_Hypo_PR_1_filtmin2)
dim(ipubcomdataBin500bp_Hypo_PR_1_filtmin2)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_freq <- ddply(MERGE_myipubCombat3re_PR_Bin500bp_Hypo, .(Bin500bp_Hypo), summarize, Probes = toString(row_4))
dim(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_freq)
head(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_freq)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_freq$Probes <- gsub(" ", "", ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_freq$Probes)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2re <- merge(ipubcomdataBin500bp_Hypo_PR_1_filtmin2, ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_freq, by.x="row", by.y="Bin500bp_Hypo", all.x=TRUE)
rownames(ipubcomdataBin500bp_Hypo_PR_1_filtmin2re) <- ipubcomdataBin500bp_Hypo_PR_1_filtmin2re$row
ipubcomdataBin500bp_Hypo_PR_1_filtmin2re1 <- merge(ipubcomdataBin500bp_Hypo_PR_1_filtmin2re, Total_CpGs_ineach_500bp,by.x="row", by.y="Bin500bp", all.x=TRUE)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2re1 <- cSplit(ipubcomdataBin500bp_Hypo_PR_1_filtmin2re1, "row", "%")
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr <- ipubcomdataBin500bp_Hypo_PR_1_filtmin2re1[,c(85:87,1:81,84,82,83)]
head(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr)
dim(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr)
write.table(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr, "ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
system("sort -k1,1 -k2,2n ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr.txt > ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr.sort.txt")
system("bedtools closest -wa -wb -a ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38.txt")
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38 <- read.table("ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38.txt", header = F, stringsAsFactors = F)
head(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38)
colnames(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38) <- c(colnames(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr), "GeneChr", "GeneStart", "GeneEnd", "Type","ENSID","Gene","Distance")
head(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38)
dim(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38 <- data.frame(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38["Bin500bp_Hypo"] <- paste0(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38$row_1,"%", ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38$row_2, "%",ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38$row_3)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_genedis <- ddply(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_genehg38, .(Bin500bp_Hypo), summarize, Gene = toString(Gene), Distance = toString(Distance))
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_genedis$Gene <- gsub(" ", "", ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_genedis$Gene)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_genedis$Distance <- gsub(" ", "", ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_genedis$Distance)
dim(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_genedis)

ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_bin <- data.frame(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_bin["Bin500bp_Hypo"] <- paste0(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_bin$row_1,"%", ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_bin$row_2, "%",ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_bin$row_3)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis <- merge(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_bin,ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_colapprob_genedis,by.x="Bin500bp_Hypo", by.y="Bin500bp_Hypo", all.x=TRUE)
head(ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis)
dim(ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis)
ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis <- ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis[,c(2:4,6:22,86:90)]
ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis["Meth_level"] <- "Hypo"
write.table(ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis, "ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = T)


dim(ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr)#921
dim(ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr)#1264

Re_ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_aggregate <- ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr[,c(1:21)]
Re_ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_aggregate <- ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr[,c(1:21)]

dim(ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr)#19
dim(ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr)#21

Re_ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_aggregate <- ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr[,c(1:21)]
Re_ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_aggregate <- ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr[,c(1:21)]
#Combine Hypermethylated and Hypomethylated bins
Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate <- rbind.data.frame(Re_ipubcomdataBin500bp_Hypo_PG_1_filtmin2chr_aggregate, Re_ipubcomdataBin500bp_Hyper_PG_1_filtmin2chr_aggregate)
Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate <- rbind.data.frame(Re_ipubcomdataBin500bp_Hypo_PR_1_filtmin2chr_aggregate, Re_ipubcomdataBin500bp_Hyper_PR_1_filtmin2chr_aggregate)

Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate <- data.frame(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate)
Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate["Bin500bp"] <- paste0(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate$row_1,
                                                                                   "%",
                                                                                   Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate$row_2,
                                                                                   "%",
                                                                                   Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate$row_3)

dim(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate)
Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate <- Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate[,c(22,4:21)]
head(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate,1)
write.table(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate, "Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)

Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate <- data.frame(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate)
Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate["Bin500bp"] <- paste0(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate$row_1,
                                                                                   "%",
                                                                                   Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate$row_2,
                                                                                   "%",
                                                                                   Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate$row_3)

dim(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate)
Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate <- Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate[,c(22,4:21)]
head(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate,1)
write.table(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate, "Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)


#Prepare Supplementary Table
Supplementary_table_PG_Hypo_Hyper <- rbind.data.frame(ipubcomdataBin500bp_Hypo_PG_1_filtmin2_featured_genedis, ipubcomdataBin500bp_Hyper_PG_1_filtmin2_featured_genedis)
dim(Supplementary_table_PG_Hypo_Hyper)
write.table(Supplementary_table_PG_Hypo_Hyper, "Supplementary_table_PG_Hypo_Hyper.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)

Supplementary_table_PR_Hypo_Hyper <- rbind.data.frame(ipubcomdataBin500bp_Hypo_PR_1_filtmin2_featured_genedis, ipubcomdataBin500bp_Hyper_PR_1_filtmin2_featured_genedis)
dim(Supplementary_table_PR_Hypo_Hyper)
write.table(Supplementary_table_PR_Hypo_Hyper, "Supplementary_table_PR_Hypo_Hyper.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)

#All bins
ipubcomdataBin500bp_filt1_counted1avg <- MERGE_myipubCombat3re_pos500bp_chravg[,1:18]
dim(ipubcomdataBin500bp_filt1_counted1avg)
head(ipubcomdataBin500bp_filt1_counted1avg)
ipubcomdataBin500bp_filt1_counted1avg["Bin500bp"] <- rownames(ipubcomdataBin500bp_filt1_counted1avg)
ipubcomdataBin500bp_filt1_counted1avg <- ipubcomdataBin500bp_filt1_counted1avg[,c(19,1:18)]
#----------Scatter Plot with affected imprinted regions belong to 500bp bins being highlighted---------#
#500bp
#PG, C13, C50
with(ipubcomdataBin500bp_filt1_counted1avg, plot(Allcontrol, PG, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="grey", bg="grey",ylab="pG (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
with(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate, points(Allcontrol, PG, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#EC7063", bg="#EC7063",ylab="pG (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
#with(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate_DMR, points(Allcontrol, PG, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 2,col="black", bg="#EC7063",ylab="pG (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
abline(lm(PG ~ Allcontrol, data = ipubcomdataBin500bp_filt1_counted1avg), col = "darkgrey", lty = 2,lwd=2)
abline(lm(PG ~ Allcontrol, data = Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate), col = "#EC7063", lty = 2,lwd=4)
#save as Scatter_Bin500bp_Re_Hypo_Hyper_PG.svg
 
with(ipubcomdataBin500bp_filt1_counted1avg, plot(Allcontrol, C13, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="grey", bg="grey",ylab="C13 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
with(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate, points(Allcontrol, C13, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#48C9B0", bg="#48C9B0",ylab="C13 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
#with(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate_DMR, points(Allcontrol, C13, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 2,col="black", bg="#48C9B0",ylab="C13 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
abline(lm(C13 ~ Allcontrol, data = ipubcomdataBin500bp_filt1_counted1avg), col = "darkgrey", lty = 2,lwd=2)
abline(lm(C13 ~ Allcontrol, data = Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate), col = "#48C9B0", lty = 2,lwd=4)
#save as Scatter_Bin500bp_Re_Hypo_Hyper_C13.svg
with(ipubcomdataBin500bp_filt1_counted1avg, plot(Allcontrol, C50, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="grey", bg="grey",ylab="C50 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
with(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate, points(Allcontrol, C50, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#52BE80", bg="#52BE80",ylab="C50 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
#with(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate_DMR, points(Allcontrol, C50, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 2,col="black", bg="#52BE80",ylab="C50 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
abline(lm(C50 ~ Allcontrol, data = ipubcomdataBin500bp_filt1_counted1avg), col = "darkgrey", lty = 2,lwd=2)
abline(lm(C50 ~ Allcontrol, data = Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate), col = "#52BE80", lty = 2,lwd=4)
#save as Scatter_Bin500bp_Re_Hypo_Hyper_C50.svg
#Use ggplot2
Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate.id <- data.frame(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate[,1])
colnames(Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate.id) <- "Bin500bp"
Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate.id["Color"] <- "diff"
ipubcomdataBin500bp_filt1_counted1_rePG <- merge(ipubcomdataBin500bp_filt1_counted1avg, Re_ipubcomdataBin500bp_Hypo_Hyper_PG_1_filtmin2chr_aggregate.id, by.x = "Bin500bp", by.y = "Bin500bp", all.x = T)
ipubcomdataBin500bp_filt1_counted1_rePG$Color[is.na(ipubcomdataBin500bp_filt1_counted1_rePG$Color)] <- "All"

ggplot(ipubcomdataBin500bp_filt1_counted1_rePG, aes(x=Allcontrol, y=PG, fill=Color, color=Color)) +
  geom_point(size=1, shape=21, alpha = 0.5)+theme_classic()+
  scale_fill_manual(values = c("grey","#EC7063"))+
  scale_color_manual(values = c("grey","#EC7063"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linetype="dashed") + ylim(c(0,1))
ggsave("Scatter_ipubcomdataBin500bp_filt1_counted1_rePG.svg", width=5, height=4, units="in", dpi=96)


ggplot(ipubcomdataBin500bp_filt1_counted1_rePG, aes(x=Allcontrol, y=C13, fill=Color, color=Color)) +
  geom_point(size=1, shape=21, alpha = 0.5)+theme_classic()+
  scale_fill_manual(values = c("grey","#48C9B0"))+
  scale_color_manual(values = c("grey","#48C9B0"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linetype="dashed") + ylim(c(0,1))
ggsave("Scatter_ipubcomdataBin500bp_filt1_counted1_reC13.svg", width=5, height=4, units="in", dpi=96)

ggplot(ipubcomdataBin500bp_filt1_counted1_rePG, aes(x=Allcontrol, y=C50, fill=Color, color=Color)) +
  geom_point(size=1, shape=21, alpha = 0.5)+theme_classic()+
  scale_fill_manual(values = c("grey","#52BE80"))+
  scale_color_manual(values = c("grey","#52BE80"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linetype="dashed") + ylim(c(0,1))
ggsave("Scatter_ipubcomdataBin500bp_filt1_counted1_reCC50.svg", width=5, height=4, units="in", dpi=96)


#PR, C7, C35

with(ipubcomdataBin500bp_filt1_counted1avg, plot(Allcontrol, PR, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="grey", bg="grey",ylab="pR (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
with(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate, points(Allcontrol, PR, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#922B21", bg="#922B21",ylab="PR (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
#with(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate_DMR, points(Allcontrol, PR, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 2,col="black", bg="#922B21",ylab="PR (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
abline(lm(PR ~ Allcontrol, data = ipubcomdataBin500bp_filt1_counted1avg), col = "darkgrey", lty = 2,lwd=2)
abline(lm(PR ~ Allcontrol, data = Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate), col = "#922B21", lty = 2,lwd=4)
#save as Scatter_Bin500bp_Hypo_Hyper_PR.svg

with(ipubcomdataBin500bp_filt1_counted1avg, plot(Allcontrol, C7, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="grey", bg="grey",ylab="C7 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
with(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate, points(Allcontrol, C7, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#148F77", bg="#148F77",ylab="C7 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
#with(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate_DMR, points(Allcontrol, C7, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 2,col="black", bg="#148F77",ylab="C7 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
abline(lm(C7 ~ Allcontrol, data = ipubcomdataBin500bp_filt1_counted1avg), col = "darkgrey", lty = 2,lwd=2)
abline(lm(C7 ~ Allcontrol, data = Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate), col = "#148F77", lty = 2,lwd=4)
#save as Scatter_Bin500bp_Hypo_Hyper_C7.svg

with(ipubcomdataBin500bp_filt1_counted1avg, plot(Allcontrol, C35, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="grey", bg="grey",ylab="C35 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
with(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate, points(Allcontrol, C35, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 0.3,col="#1E8449", bg="#1E8449",ylab="C35 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
#with(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate_DMR, points(Allcontrol, C35, pch=21, cex = 0.8, main="Scatter Plot", xlim=c(0,1), ylim=c(0,1),lwd = 2,col="black", bg="#1E8449",ylab="C35 (Beta-value)", xlab="Median Controls (Beta-value)", bty = 'n'))
abline(lm(C35 ~ Allcontrol, data = ipubcomdataBin500bp_filt1_counted1avg), col = "darkgrey", lty = 2,lwd=2)
abline(lm(C35 ~ Allcontrol, data = Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate), col = "#1E8449", lty = 2,lwd=4)
#save as Scatter_Bin500bp_Hypo_Hyper_C35.svg


#Use ggplot2
Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate.id <- data.frame(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate[,1])
colnames(Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate.id) <- "Bin500bp"
Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate.id["Color"] <- "diff"
ipubcomdataBin500bp_filt1_counted1_rePR <- merge(ipubcomdataBin500bp_filt1_counted1avg, Re_ipubcomdataBin500bp_Hypo_Hyper_PR_1_filtmin2chr_aggregate.id, by.x = "Bin500bp", by.y = "Bin500bp", all.x = T)
ipubcomdataBin500bp_filt1_counted1_rePR$Color[is.na(ipubcomdataBin500bp_filt1_counted1_rePR$Color)] <- "All"

ggplot(ipubcomdataBin500bp_filt1_counted1_rePR, aes(x=Allcontrol, y=PR, fill=Color, color=Color)) +
  geom_point(size=1, shape=21, alpha = 0.5)+theme_classic()+
  scale_fill_manual(values = c("grey","#922B21"))+
  scale_color_manual(values = c("grey","#922B21"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linetype="dashed") + ylim(c(0,1))
ggsave("Scatter_ipubcomdataBin500bp_filt1_counted1_rePR.svg", width=5, height=4, units="in", dpi=96)


ggplot(ipubcomdataBin500bp_filt1_counted1_rePR, aes(x=Allcontrol, y=C7, fill=Color, color=Color)) +
  geom_point(size=1, shape=21, alpha = 0.5)+theme_classic()+
  scale_fill_manual(values = c("grey","#148F77"))+
  scale_color_manual(values = c("grey","#148F77"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linetype="dashed") + ylim(c(0,1))
ggsave("Scatter_ipubcomdataBin500bp_filt1_counted1_reC7.svg", width=5, height=4, units="in", dpi=96)

ggplot(ipubcomdataBin500bp_filt1_counted1_rePR, aes(x=Allcontrol, y=C35, fill=Color, color=Color)) +
  geom_point(size=1, shape=21, alpha = 0.5)+theme_classic()+
  scale_fill_manual(values = c("grey","#1E8449"))+
  scale_color_manual(values = c("grey","#1E8449"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linetype="dashed") + ylim(c(0,1))
ggsave("Scatter_ipubcomdataBin500bp_filt1_counted1_reCC35.svg", width=5, height=4, units="in", dpi=96)




print('--------------------  Bin 500 bp :Single dm CpG_Hypos overlap bins way---------------------------')
detach("package:dplyr")
#Peek into genome wide view
system("bedtools makewindows -g /home/ankitv/ref_av/hg38/hg38.chrom.sizes -w 500 > hg38_500.bed")

#Bin500bp PG_Hypo
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt -b hg38_500.bed > MERGE_myipubCombat3re_PG_Hypo_Bin500bp.txt")
MERGE_myipubCombat3re_PG_Hypo_Bin500bp <- read.table("MERGE_myipubCombat3re_PG_Hypo_Bin500bp.txt", header = F, stringsAsFactors = F)
head(MERGE_myipubCombat3re_PG_Hypo_Bin500bp)
colnames(MERGE_myipubCombat3re_PG_Hypo_Bin500bp) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo),"Binchr","Binstart","Binend")
MERGE_myipubCombat3re_PG_Hypo_Bin500bp["Bin500bp"] <- paste0(MERGE_myipubCombat3re_PG_Hypo_Bin500bp$Binchr, "%",
                                                             MERGE_myipubCombat3re_PG_Hypo_Bin500bp$Binstart, "%",
                                                             MERGE_myipubCombat3re_PG_Hypo_Bin500bp$Binend)

head(MERGE_myipubCombat3re_PG_Hypo_Bin500bp)
dim(MERGE_myipubCombat3re_PG_Hypo_Bin500bp)
MERGE_myipubCombat3re_PG_Hypo_Bin500bp_rearranged <- MERGE_myipubCombat3re_PG_Hypo_Bin500bp[,c(89,1:88)]
head(MERGE_myipubCombat3re_PG_Hypo_Bin500bp_rearranged)
ipubcomdataBin500bp_PG_Hypo <- MERGE_myipubCombat3re_PG_Hypo_Bin500bp_rearranged
rownames(ipubcomdataBin500bp_PG_Hypo)
colnames(ipubcomdataBin500bp_PG_Hypo)
head(ipubcomdataBin500bp_PG_Hypo)
dim(ipubcomdataBin500bp_PG_Hypo)
Count_ipubcomdataBin500bp_PG_Hypo_1 <- count(ipubcomdataBin500bp_PG_Hypo, "Bin500bp")
head(Count_ipubcomdataBin500bp_PG_Hypo_1)
dim(Count_ipubcomdataBin500bp_PG_Hypo_1)

Count_ipubcomdataBin500bp_PG_Hypo_1_filtmin3 <- Count_ipubcomdataBin500bp_PG_Hypo_1[which(Count_ipubcomdataBin500bp_PG_Hypo_1$freq >=3),]
head(Count_ipubcomdataBin500bp_PG_Hypo_1_filtmin3)
dim(Count_ipubcomdataBin500bp_PG_Hypo_1_filtmin3)

ipubcomdataBin500bp_PG_Hypo_1_filtmin3 <- merge(ipubcomdataBin500bp_PG_Hypo, Count_ipubcomdataBin500bp_PG_Hypo_1_filtmin3, by = "Bin500bp", all.y = TRUE ) 
head(ipubcomdataBin500bp_PG_Hypo_1_filtmin3)
dim(ipubcomdataBin500bp_PG_Hypo_1_filtmin3)
ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr <- ipubcomdataBin500bp_PG_Hypo_1_filtmin3[,c(87:90,1:86)]
head(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr)
write.table(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr, "ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
system("sort -k1,1 -k2,2n ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr.txt > ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr.sort.txt")
system("bedtools closest -wa -wb -a ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38.txt")
#Collapse genenames and Prepare Table
ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38 <- read.table("ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38.txt", header = F, stringsAsFactors = F)
head(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38)
colnames(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38) <- c(colnames(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr), "GeneChr", "GeneStart", "GeneEnd", "Type","ENSID","Gene","Distance")
head(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38)
ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38dis <- ddply(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38, .(Bin500bp), summarize, Gene = toString(Gene), Dis = toString(Distance))
dim(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38dis)
head(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38dis)
ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_named <- merge(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr, 
                                                         ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_genehg38dis,
                                                         by = "Bin500bp")

ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_named$Genes <- gsub(" ", "", ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_named$Gene)
ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_named$Distance <- gsub(" ", "", ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_named$Dis)
head(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_named)
write.table(ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_named, "ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_named.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)


system("bedtools merge -i ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr.sort.txt -d 2 | awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' | sort -k1,1  -k2,2n > ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_clustered.txt")
system("bedtools closest -a ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_clustered.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene.txt")

ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene <- read.table("ipubcomdataBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene.txt", header = F, stringsAsFactors = F)
ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene["row"] <- paste0(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene$V1,
                                                                   "%",
                                                                   ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene$V2,
                                                                   "%",
                                                                   ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene$V3,
                                                                   "%",
                                                                   ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene$V4)

head(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene)
ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedis <- ddply(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgene, .(row), summarize, Gene = toString(V10), Dis = toString(V11))
dim(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedis)
head(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedis)
ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedis$Genes <- gsub(" ", "", ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedis$Gene)
ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedis$Distance <- gsub(" ", "", ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedis$Dis)
ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr <- cSplit(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedis, "row", "%")
head(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr)
ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr <- ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr[,c(5:8,3,4)]
ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr <- ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr[order(-ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr$row_4),]
head(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr,10)
write.table(ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr, "ipubcBin500bp_PG_Hypo_1_filtmin3chr_clusteredgenedischr.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)



#Bin500bp PR
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt -b hg38_500.bed > MERGE_myipubCombat3re_PR_Hypo_Bin500bp.txt")
MERGE_myipubCombat3re_PR_Hypo_Bin500bp <- read.table("MERGE_myipubCombat3re_PR_Hypo_Bin500bp.txt", header = F, stringsAsFactors = F)
head(MERGE_myipubCombat3re_PR_Hypo_Bin500bp)
colnames(MERGE_myipubCombat3re_PR_Hypo_Bin500bp) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo),"Binchr","Binstart","Binend")
MERGE_myipubCombat3re_PR_Hypo_Bin500bp["Bin500bp"] <- paste0(MERGE_myipubCombat3re_PR_Hypo_Bin500bp$Binchr, "%",
                                                             MERGE_myipubCombat3re_PR_Hypo_Bin500bp$Binstart, "%",
                                                             MERGE_myipubCombat3re_PR_Hypo_Bin500bp$Binend)

head(MERGE_myipubCombat3re_PR_Hypo_Bin500bp)
dim(MERGE_myipubCombat3re_PR_Hypo_Bin500bp)
MERGE_myipubCombat3re_PR_Hypo_Bin500bp_rearranged <- MERGE_myipubCombat3re_PR_Hypo_Bin500bp[,c(89,1:88)]
head(MERGE_myipubCombat3re_PR_Hypo_Bin500bp_rearranged)
ipubcomdataBin500bp_PR_Hypo <- MERGE_myipubCombat3re_PR_Hypo_Bin500bp_rearranged
rownames(ipubcomdataBin500bp_PR_Hypo)
colnames(ipubcomdataBin500bp_PR_Hypo)
head(ipubcomdataBin500bp_PR_Hypo)
dim(ipubcomdataBin500bp_PR_Hypo)
Count_ipubcomdataBin500bp_PR_Hypo_1 <- count(ipubcomdataBin500bp_PR_Hypo, "Bin500bp")
head(Count_ipubcomdataBin500bp_PR_Hypo_1)
dim(Count_ipubcomdataBin500bp_PR_Hypo_1)

Count_ipubcomdataBin500bp_PR_Hypo_1_filtmin3 <- Count_ipubcomdataBin500bp_PR_Hypo_1[which(Count_ipubcomdataBin500bp_PR_Hypo_1$freq >=3),]
head(Count_ipubcomdataBin500bp_PR_Hypo_1_filtmin3)
dim(Count_ipubcomdataBin500bp_PR_Hypo_1_filtmin3)

ipubcomdataBin500bp_PR_Hypo_1_filtmin3 <- merge(ipubcomdataBin500bp_PR_Hypo, Count_ipubcomdataBin500bp_PR_Hypo_1_filtmin3, by = "Bin500bp", all.y = TRUE ) 
head(ipubcomdataBin500bp_PR_Hypo_1_filtmin3)
dim(ipubcomdataBin500bp_PR_Hypo_1_filtmin3)
ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr <- ipubcomdataBin500bp_PR_Hypo_1_filtmin3[,c(87:90,1:86)]
head(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr)
write.table(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr, "ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = F, col.names = F)
system("sort -k1,1 -k2,2n ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr.txt > ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr.sort.txt")
system("bedtools closest -wa -wb -a ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38.txt")
#Collapse genenames and Prepare Table
ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38 <- read.table("ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38.txt", header = F, stringsAsFactors = F)
head(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38)
colnames(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38) <- c(colnames(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr), "GeneChr", "GeneStart", "GeneEnd", "Type","ENSID","Gene","Distance")
head(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38)
ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38dis <- ddply(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38, .(Bin500bp), summarize, Gene = toString(Gene), Dis = toString(Distance))
dim(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38dis)
head(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38dis)
ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_named <- merge(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr, 
                                                         ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_genehg38dis,
                                                         by = "Bin500bp")

ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_named$Genes <- gsub(" ", "", ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_named$Gene)
ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_named$Distance <- gsub(" ", "", ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_named$Dis)
head(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_named)
write.table(ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_named, "ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_named.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)


system("bedtools merge -i ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr.sort.txt -d 2 | awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' | sort -k1,1  -k2,2n > ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_clustered.txt")
system("bedtools closest -a ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_clustered.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene.txt")

ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene <- read.table("ipubcomdataBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene.txt", header = F, stringsAsFactors = F)
ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene["row"] <- paste0(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene$V1,
                                                                   "%",
                                                                   ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene$V2,
                                                                   "%",
                                                                   ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene$V3,
                                                                   "%",
                                                                   ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene$V4)

head(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene)
ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedis <- ddply(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgene, .(row), summarize, Gene = toString(V10), Dis = toString(V11))
dim(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedis)
head(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedis)
ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedis$Genes <- gsub(" ", "", ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedis$Gene)
ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedis$Distance <- gsub(" ", "", ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedis$Dis)
ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr <- cSplit(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedis, "row", "%")
head(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr)
ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr <- ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr[,c(5:8,3,4)]
ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr <- ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr[order(-ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr$row_4),]
head(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr,10)
write.table(ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr, "ipubcBin500bp_PR_Hypo_1_filtmin3chr_clusteredgenedischr.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)

#Note with closest all bins are assigned with genes
system("bedtools merge -i ipubcomaggregatedBin500bp_countedBin500bp_PGchr_genehg38.txt -d 2 | awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' | sort -k1,1  -k2,2n > ipubcomaggregatedBin500bp_countedBin500bp_PGchr_genehg38_clustered.txt")
system("bedtools closest -a ipubcomaggregatedBin500bp_countedBin500bp_PGchr_genehg38_clustered.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomaggregatedBin500bp_countedBin500bp_PGchr_genehg38_clusteredgene.txt")

ipubcBin500bp_PGhg38_clusteredgene <- read.table("ipubcomaggregatedBin500bp_countedBin500bp_PGchr_genehg38_clusteredgene.txt", header = F, stringsAsFactors = F)
ipubcBin500bp_PGhg38_clusteredgene["row"] <- paste0(ipubcBin500bp_PGhg38_clusteredgene$V1,
                                                    "%",
                                                    ipubcBin500bp_PGhg38_clusteredgene$V2,
                                                    "%",
                                                    ipubcBin500bp_PGhg38_clusteredgene$V3,
                                                    "%",
                                                    ipubcBin500bp_PGhg38_clusteredgene$V4)

head(ipubcBin500bp_PGhg38_clusteredgene)
library(plyr)
ipubcBin500bp_PGhg38_clusteredgenedis <- ddply(ipubcBin500bp_PGhg38_clusteredgene, .(row), summarize, Gene = toString(V10), Dis = toString(V11))
dim(ipubcBin500bp_PGhg38_clusteredgenedis)
head(ipubcBin500bp_PGhg38_clusteredgenedis)
ipubcBin500bp_PGhg38_clusteredgenedis$Genes <- gsub(" ", "", ipubcBin500bp_PGhg38_clusteredgenedis$Gene)
ipubcBin500bp_PGhg38_clusteredgenedis$Distance <- gsub(" ", "", ipubcBin500bp_PGhg38_clusteredgenedis$Dis)
library(splitstackshape)
ipubcBin500bp_PGhg38_clusteredgenedischr <- cSplit(ipubcBin500bp_PGhg38_clusteredgenedis, "row", "%")
head(ipubcBin500bp_PGhg38_clusteredgenedischr)
ipubcBin500bp_PGhg38_clusteredgenedischr <- ipubcBin500bp_PGhg38_clusteredgenedischr[,c(5:8,3,4)]
ipubcBin500bp_PGhg38_clusteredgenedischr <- ipubcBin500bp_PGhg38_clusteredgenedischr[order(-ipubcBin500bp_PGhg38_clusteredgenedischr$row_4),]
head(ipubcBin500bp_PGhg38_clusteredgenedischr,20)
write.table(ipubcBin500bp_PGhg38_clusteredgenedischr, "ipubcBin500bp_PGhg38_clusteredgenedischr.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)


#xxxxx
#Note with closest all bins are assigned with genes
system("bedtools merge -i ipubcomaggregatedBin500bp_countedBin500bp_PRchr_genehg38.txt -d 2 | awk '{print $1"\t"$2"\t"$3"\t"$3-$2}' | sort -k1,1  -k2,2n > ipubcomaggregatedBin500bp_countedBin500bp_PRchr_genehg38_clustered.txt")
system("bedtools closest -a ipubcomaggregatedBin500bp_countedBin500bp_PRchr_genehg38_clustered.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomaggregatedBin500bp_countedBin500bp_PRchr_genehg38_clusteredgene.txt")

ipubcBin500bp_PRhg38_clusteredgene <- read.table("ipubcomaggregatedBin500bp_countedBin500bp_PRchr_genehg38_clusteredgene.txt", header = F, stringsAsFactors = F)
ipubcBin500bp_PRhg38_clusteredgene["row"] <- paste0(ipubcBin500bp_PRhg38_clusteredgene$V1,
                                                    "%",
                                                    ipubcBin500bp_PRhg38_clusteredgene$V2,
                                                    "%",
                                                    ipubcBin500bp_PRhg38_clusteredgene$V3,
                                                    "%",
                                                    ipubcBin500bp_PRhg38_clusteredgene$V4)

head(ipubcBin500bp_PRhg38_clusteredgene)
ipubcBin500bp_PRhg38_clusteredgenedis <- ddply(ipubcBin500bp_PRhg38_clusteredgene, .(row), summarize, Gene = toString(V10), Dis = toString(V11))
dim(ipubcBin500bp_PRhg38_clusteredgenedis)
head(ipubcBin500bp_PRhg38_clusteredgenedis)
ipubcBin500bp_PRhg38_clusteredgenedis$Genes <- gsub(" ", "", ipubcBin500bp_PRhg38_clusteredgenedis$Gene)
ipubcBin500bp_PRhg38_clusteredgenedis$Distance <- gsub(" ", "", ipubcBin500bp_PRhg38_clusteredgenedis$Dis)
ipubcBin500bp_PRhg38_clusteredgenedischr <- cSplit(ipubcBin500bp_PRhg38_clusteredgenedis, "row", "%")
head(ipubcBin500bp_PRhg38_clusteredgenedischr)
ipubcBin500bp_PRhg38_clusteredgenedischr <- ipubcBin500bp_PRhg38_clusteredgenedischr[,c(5:8,3,4)]
ipubcBin500bp_PRhg38_clusteredgenedischr <- ipubcBin500bp_PRhg38_clusteredgenedischr[order(-ipubcBin500bp_PRhg38_clusteredgenedischr$row_4),]
head(ipubcBin500bp_PRhg38_clusteredgenedischr,20)
write.table(ipubcBin500bp_PRhg38_clusteredgenedischr, "ipubcBin500bp_PRhg38_clusteredgenedischr.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)


#Extract CpGs inside perturbed cluster
system("grep TNXB ipubcBin500bp_PGhg38_clusteredgenedischr.txt > ipubcBin500bp_PGhg38_clusteredgenedischr_TNXB.txt")
system("bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b ipubcBin500bp_PGhg38_clusteredgenedischr_TNXB.txt > MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXB.txt")
MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXB <- read.table("MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXB.txt")
colnames(MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXB) <- c(colnames(MERGE_myiCombat3reavg_pos_chr), "chr","start","end","span","gene","distance")
head(MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXB)
MipubcB5CpG_PGhg38_clust_TNXB <- MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXB


MipubcB5CpG_PGhg38_clust_TNXB["Row"] <- paste0(MipubcB5CpG_PGhg38_clust_TNXB$chrhg38,"%",MipubcB5CpG_PGhg38_clust_TNXB$starthg38,"%",MipubcB5CpG_PGhg38_clust_TNXB$endhg38,"%",MipubcB5CpG_PGhg38_clust_TNXB$TargetID)
rownames(MipubcB5CpG_PGhg38_clust_TNXB) <- MipubcB5CpG_PGhg38_clust_TNXB$Row
head(MipubcB5CpG_PGhg38_clust_TNXB)
dim(MipubcB5CpG_PGhg38_clust_TNXB)
MipubcB5CpG_PGhg38_clust_TNXB <- MipubcB5CpG_PGhg38_clust_TNXB[,-c(1:5,23:29)]
head(MipubcB5CpG_PGhg38_clust_TNXB)
MipubcB5CpG_PGhg38_clust_TNXB <- as.matrix(MipubcB5CpG_PGhg38_clust_TNXB)
MipubcB5CpG_PGhg38_clust_TNXBst <- stack(MipubcB5CpG_PGhg38_clust_TNXB)
head(MipubcB5CpG_PGhg38_clust_TNXBst)
MipubcB5CpG_PGhg38_clust_TNXBst <- data.frame(MipubcB5CpG_PGhg38_clust_TNXBst)
MipubcB5CpG_PGhg38_clust_TNXBst <- MipubcB5CpG_PGhg38_clust_TNXBst[,c(2,1,3)]
colnames(MipubcB5CpG_PGhg38_clust_TNXBst) <- c("Group", "Csites", "Value")
head(MipubcB5CpG_PGhg38_clust_TNXBst)


MipubcB5CpG_PGhg38_clust_TNXBst.sep <- cSplit(MipubcB5CpG_PGhg38_clust_TNXBst, "Csites", "%")
ggplot(MipubcB5CpG_PGhg38_clust_TNXBst.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(32095500,	32097500), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20)+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_MipubcB5CpG_PGhg38_clust_TNXBst.sep.png", width=25, height=5, units="cm", dpi=96)


system("grep TNXB ipubcBin500bp_PGhg38_clusteredgenedischr.txt | grep 2000 > ipubcBin500bp_PGhg38_clusteredgenedischr_TNXBcontspan.txt")
system("bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b ipubcBin500bp_PGhg38_clusteredgenedischr_TNXBcontspan.txt > MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXBcontspan.txt")
MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXBcontspan <- read.table("MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXBcontspan.txt")
colnames(MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXBcontspan) <- c(colnames(MERGE_myiCombat3reavg_pos_chr), "chr","start","end","span","gene","distance")
head(MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXBcontspan)
MipubcB5CpG_PGhg38_clust_TNXBcontspan <- MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_TNXBcontspan


MipubcB5CpG_PGhg38_clust_TNXBcontspan["Row"] <- paste0(MipubcB5CpG_PGhg38_clust_TNXBcontspan$chrhg38,"%",MipubcB5CpG_PGhg38_clust_TNXBcontspan$starthg38,"%",MipubcB5CpG_PGhg38_clust_TNXBcontspan$endhg38,"%",MipubcB5CpG_PGhg38_clust_TNXBcontspan$TargetID)
rownames(MipubcB5CpG_PGhg38_clust_TNXBcontspan) <- MipubcB5CpG_PGhg38_clust_TNXBcontspan$Row
head(MipubcB5CpG_PGhg38_clust_TNXBcontspan)
dim(MipubcB5CpG_PGhg38_clust_TNXBcontspan)
MipubcB5CpG_PGhg38_clust_TNXBcontspan <- MipubcB5CpG_PGhg38_clust_TNXBcontspan[,-c(1:5,23:29)]
head(MipubcB5CpG_PGhg38_clust_TNXBcontspan)
MipubcB5CpG_PGhg38_clust_TNXBcontspan <- as.matrix(MipubcB5CpG_PGhg38_clust_TNXBcontspan)
MipubcB5CpG_PGhg38_clust_TNXBcontspanst <- stack(MipubcB5CpG_PGhg38_clust_TNXBcontspan)
head(MipubcB5CpG_PGhg38_clust_TNXBcontspanst)
MipubcB5CpG_PGhg38_clust_TNXBcontspanst <- data.frame(MipubcB5CpG_PGhg38_clust_TNXBcontspanst)
MipubcB5CpG_PGhg38_clust_TNXBcontspanst <- MipubcB5CpG_PGhg38_clust_TNXBcontspanst[,c(2,1,3)]
colnames(MipubcB5CpG_PGhg38_clust_TNXBcontspanst) <- c("Group", "Csites", "Value")
head(MipubcB5CpG_PGhg38_clust_TNXBcontspanst)


MipubcB5CpG_PGhg38_clust_TNXBcontspanst.sep <- cSplit(MipubcB5CpG_PGhg38_clust_TNXBcontspanst, "Csites", "%")
ggplot(MipubcB5CpG_PGhg38_clust_TNXBcontspanst.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(32095500,	32097500), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20)+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_MipubcB5CpG_PGhg38_clust_TNXBcontspanst.sep.png", width=25, height=5, units="cm", dpi=96)

#No need to do for pR because full matrix is aleady covered (both pG and pR and respective correction) and exactly same TNXB region is present in pG and pR

#DIP2C
system("grep DIP2C ipubcBin500bp_PGhg38_clusteredgenedischr.txt | grep 485 > ipubcBin500bp_PGhg38_clusteredgenedischr_DIP2C.txt")
system("bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b ipubcBin500bp_PGhg38_clusteredgenedischr_DIP2C.txt > MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_DIP2C.txt")
MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_DIP2C <- read.table("MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_DIP2C.txt")
colnames(MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_DIP2C) <- c(colnames(MERGE_myiCombat3reavg_pos_chr), "chr","start","end","span","gene","distance")
head(MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_DIP2C)
MipubcB5CpG_PGhg38_clust_DIP2C <- MERGE_ipubcBin500bp_PGhg38_clusteredgenedischr_DIP2C


MipubcB5CpG_PGhg38_clust_DIP2C["Row"] <- paste0(MipubcB5CpG_PGhg38_clust_DIP2C$chrhg38,"%",MipubcB5CpG_PGhg38_clust_DIP2C$starthg38,"%",MipubcB5CpG_PGhg38_clust_DIP2C$endhg38,"%",MipubcB5CpG_PGhg38_clust_DIP2C$TargetID)
rownames(MipubcB5CpG_PGhg38_clust_DIP2C) <- MipubcB5CpG_PGhg38_clust_DIP2C$Row
head(MipubcB5CpG_PGhg38_clust_DIP2C)
dim(MipubcB5CpG_PGhg38_clust_DIP2C)
MipubcB5CpG_PGhg38_clust_DIP2C <- MipubcB5CpG_PGhg38_clust_DIP2C[,-c(1:5,23:29)]
head(MipubcB5CpG_PGhg38_clust_DIP2C)
MipubcB5CpG_PGhg38_clust_DIP2C <- as.matrix(MipubcB5CpG_PGhg38_clust_DIP2C)
MipubcB5CpG_PGhg38_clust_DIP2Cst <- stack(MipubcB5CpG_PGhg38_clust_DIP2C)
head(MipubcB5CpG_PGhg38_clust_DIP2Cst)
MipubcB5CpG_PGhg38_clust_DIP2Cst <- data.frame(MipubcB5CpG_PGhg38_clust_DIP2Cst)
MipubcB5CpG_PGhg38_clust_DIP2Cst <- MipubcB5CpG_PGhg38_clust_DIP2Cst[,c(2,1,3)]
colnames(MipubcB5CpG_PGhg38_clust_DIP2Cst) <- c("Group", "Csites", "Value")
head(MipubcB5CpG_PGhg38_clust_DIP2Cst)


MipubcB5CpG_PGhg38_clust_DIP2Cst.sep <- cSplit(MipubcB5CpG_PGhg38_clust_DIP2Cst, "Csites", "%")
ggplot(MipubcB5CpG_PGhg38_clust_DIP2Cst.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(484500,	485500), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20)+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_MipubcB5CpG_PGhg38_clust_DIP2Cst.sep.png", width=25, height=5, units="cm", dpi=96)




#DIP2C
system("grep DIP2C ipubcBin500bp_PRhg38_clusteredgenedischr.txt | grep 484 > ipubcBin500bp_PRhg38_clusteredgenedischr_DIP2C.txt")
system("bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b ipubcBin500bp_PRhg38_clusteredgenedischr_DIP2C.txt > MERGE_ipubcBin500bp_PRhg38_clusteredgenedischr_DIP2C.txt")
MERGE_ipubcBin500bp_PRhg38_clusteredgenedischr_DIP2C <- read.table("MERGE_ipubcBin500bp_PRhg38_clusteredgenedischr_DIP2C.txt")
colnames(MERGE_ipubcBin500bp_PRhg38_clusteredgenedischr_DIP2C) <- c(colnames(MERGE_myiCombat3reavg_pos_chr), "chr","start","end","span","gene","distance")
head(MERGE_ipubcBin500bp_PRhg38_clusteredgenedischr_DIP2C)
MipubcB5CpG_PRhg38_clust_DIP2C <- MERGE_ipubcBin500bp_PRhg38_clusteredgenedischr_DIP2C


MipubcB5CpG_PRhg38_clust_DIP2C["Row"] <- paste0(MipubcB5CpG_PRhg38_clust_DIP2C$chrhg38,"%",MipubcB5CpG_PRhg38_clust_DIP2C$starthg38,"%",MipubcB5CpG_PRhg38_clust_DIP2C$endhg38,"%",MipubcB5CpG_PRhg38_clust_DIP2C$TargetID)
rownames(MipubcB5CpG_PRhg38_clust_DIP2C) <- MipubcB5CpG_PRhg38_clust_DIP2C$Row
head(MipubcB5CpG_PRhg38_clust_DIP2C)
dim(MipubcB5CpG_PRhg38_clust_DIP2C)
MipubcB5CpG_PRhg38_clust_DIP2C <- MipubcB5CpG_PRhg38_clust_DIP2C[,-c(1:5,23:29)]
head(MipubcB5CpG_PRhg38_clust_DIP2C)
MipubcB5CpG_PRhg38_clust_DIP2C <- as.matrix(MipubcB5CpG_PRhg38_clust_DIP2C)
MipubcB5CpG_PRhg38_clust_DIP2Cst <- stack(MipubcB5CpG_PRhg38_clust_DIP2C)
head(MipubcB5CpG_PRhg38_clust_DIP2Cst)
MipubcB5CpG_PRhg38_clust_DIP2Cst <- data.frame(MipubcB5CpG_PRhg38_clust_DIP2Cst)
MipubcB5CpG_PRhg38_clust_DIP2Cst <- MipubcB5CpG_PRhg38_clust_DIP2Cst[,c(2,1,3)]
colnames(MipubcB5CpG_PRhg38_clust_DIP2Cst) <- c("Group", "Csites", "Value")
head(MipubcB5CpG_PRhg38_clust_DIP2Cst)


MipubcB5CpG_PRhg38_clust_DIP2Cst.sep <- cSplit(MipubcB5CpG_PRhg38_clust_DIP2Cst, "Csites", "%")
ggplot(MipubcB5CpG_PRhg38_clust_DIP2Cst.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(484500,	486500), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20)+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_MipubcB5CpG_PRhg38_clust_DIP2Cst.sep.png", width=25, height=5, units="cm", dpi=96)


#Only PR specific is ok as it span hypo regions for both pG and pR


#--------------------------------------

print('-------------Assign probes to Human ICRs: hg38-------------------')
system("bedtools intersect -wa -wb -a MERGE_myiCombat3re_pos_chr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed > MERGE_myiCombat3re_human_ICR.txt")

MERGE_myiCombat3re_human_ICR <- read.table("MERGE_myiCombat3re_human_ICR.txt", header = F)
head(MERGE_myiCombat3re_human_ICR)
colnames(MERGE_myiCombat3re_human_ICR) <- c(colnames(MERGE_myiCombat3re_pos_chr),"chr","start","end","DMR", "DMRType")
#PCA with  Human ICR
icomdata1 <- MERGE_myiCombat3re_human_ICR[,c((ncol(MERGE_myiCombat3re_human_ICR)-1),ncol(MERGE_myiCombat3re_human_ICR), 1:(ncol(MERGE_myiCombat3re_human_ICR)-5))]
rownames(icomdata1)
colnames(icomdata1)
head(icomdata1,1)
dim(icomdata1)
Count_icomdataICR1 <- count(icomdata1, "DMR")
head(Count_icomdataICR1)
dim(Count_icomdataICR1)
scomaggregate1 = aggregate(icomdata1[,7:ncol(icomdata1)],by=list(icomdata1$DMR), mean)
head(scomaggregate1, 2)
#Plot PCA by group color and labelling
icomdataICR1=scomaggregate1
rownames(icomdataICR1)
icomdataICR1[,1]
rownames(icomdataICR1)=icomdataICR1[,1]
rownames(icomdataICR1)
colnames(icomdataICR1)
icomdataICR1 = icomdataICR1[,-1]
head(icomdataICR1)
dim(icomdataICR1)
icomdataICR1_counted <- cbind(icomdataICR1, Count_icomdataICR1)
head(icomdataICR1_counted)
dim(icomdataICR1_counted)
icomdataICR1_counted1 <- icomdataICR1_counted[which(icomdataICR1_counted$freq >= 3),]
head(icomdataICR1_counted1)
dim(icomdataICR1_counted1)
write.table(icomdataICR1_counted1, "icomdataICR1_counted1.filtmin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

icomdataICR1_counted2 <- icomdataICR1_counted1[,1:(ncol(icomdataICR1_counted1)-2)]
head(icomdataICR1_counted2)
summary(icomdataICR1_counted2)

print('PCA for Human ICR')

icomdfICR <- icomdataICR1_counted2
head(icomdfICR)
dim(icomdfICR)
icomdfICR = data.frame(icomdfICR)
#write.table(inormdfICR , "inormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
ticomdfICR = t(icomdfICR)
ticomdfICR = data.frame(ticomdfICR)
#write.table(tinormdfICR , "tinormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(ticomdfICR)
ticomdfICR["Color"] <- adjPd$Sample_Group
#write.table(tinormdfICR , "tinormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(ticomdfICR)
icomdfx <-ticomdfICR[c(1:nrow(icomdfICR))]
PicomC<-prcomp(icomdfx, center = TRUE, scale. = TRUE)
PicomCi<-data.frame(PicomC$x,Color=ticomdfICR$Color)
percentagecomICR <- round(PicomC$sdev^2 / sum(PicomC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagecomICR <- paste( colnames(PicomCi), "(", paste( as.character(percentagecomICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pcom1 <-ggplot(PicomCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagecomICR[1]) + ylab(percentagecomICR[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("red","darkgrey","green"))+
  scale_shape_manual(values= c(0,1,2))
pcom1 <- pcom1 +theme_classic()
pcom1 + xlim(-20,20)+ ylim(-20,20)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_ticomdfICR_CntrlIndiv.svg", width=15*1.25, height=10*1.25, units="cm", dpi=96)

################################### Heatmap_indiv human ICR ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(icomdataICR1_counted1)
dim(icomdataICR1_counted1)
sum(icomdataICR1_counted1$freq)

#Heatmap_indiv with 11 Controls median Filtered by minimum 3 CpGs

icomdataICR2 = as.matrix(icomdataICR1_counted2)
head(icomdataICR2)
dim(icomdataICR2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
icomdataICR3 <- icomdataICR2
head(icomdataICR3)
write.table(icomdataICR3, "Heatmap_indiv_icomdataICR3_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
head(icomdataICR1_counted2)

#Heatmap of other samples merged
head(icomdataICR1_counted2)
dim(icomdataICR1_counted2)
newarray_iPSC_Control_ICR <- icomdataICR1_counted2
my_sample_col2 <- data.frame(Annotations= c("Corrclone","Corrclone","Corrclone", "Corrclone","Control","Case","Case","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control"))
row.names(my_sample_col2) <- colnames(newarray_iPSC_Control_ICR)
my_colour2 = list(Annotations = c("Control"= "darkgrey","Case"= "red", "Corrclone"="green"))

breaksList3 = seq(0, 1, by = 0.01)

pheatmap(newarray_iPSC_Control_ICR,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList3)),
         breaks = breaksList3,
         annotation_colors = my_colour2,
         fontsize = 8,
         annotation_col = my_sample_col2,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 3)
#save as pheatmap_icomdataICR1_counted2ratio.svg
head(newarray_iPSC_Control_ICR,1)
newarray_iPSC_Control_ICR_counted2ration <- as.matrix(newarray_iPSC_Control_ICR[,c(5,8,9:17,1:4,6,7)])

head(newarray_iPSC_Control_ICR_counted2ration)
dim(newarray_iPSC_Control_ICR_counted2ration)
wnewarray_iPSC_Control_ICR_counted2rationavg <- data.frame(cbind((rowMedians(newarray_iPSC_Control_ICR_counted2ration[,1:11])),
                                                    (newarray_iPSC_Control_ICR_counted2ration[,1:17])))

head(wnewarray_iPSC_Control_ICR_counted2rationavg)
colnames(wnewarray_iPSC_Control_ICR_counted2rationavg) <- c("Allcontrol",colnames(newarray_iPSC_Control_ICR_counted2ration))
head(wnewarray_iPSC_Control_ICR_counted2rationavg)
wnewarray_iPSC_Control_ICR_counted2rationavg1 <- wnewarray_iPSC_Control_ICR_counted2rationavg
wnewarray_iPSC_Control_ICR_counted2rationavg2 <- as.matrix(wnewarray_iPSC_Control_ICR_counted2rationavg1)
head(wnewarray_iPSC_Control_ICR_counted2rationavg2)
dim(wnewarray_iPSC_Control_ICR_counted2rationavg2)
mnewarray_iPSC_Control_ICR_counted2rationavg3 <- data.frame(cbind((wnewarray_iPSC_Control_ICR_counted2rationavg2[,1:18]-wnewarray_iPSC_Control_ICR_counted2rationavg2[,1])))


head(mnewarray_iPSC_Control_ICR_counted2rationavg3)
mnewarray_iPSC_Control_ICR_counted2rationavg3 <- as.matrix(mnewarray_iPSC_Control_ICR_counted2rationavg3)
dim(mnewarray_iPSC_Control_ICR_counted2rationavg3)
mnewarray_iPSC_Control_ICR_counted2rationavg3minuscontrol <- mnewarray_iPSC_Control_ICR_counted2rationavg3[,c(2:18)]
head(mnewarray_iPSC_Control_ICR_counted2rationavg3minuscontrol,1)
my_sample_col5 <- data.frame(Annotations= c("Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Control","Corrclone7","Corrclone13","Corrclone35","Corrclone50","CasepG","CasepR"))
row.names(my_sample_col5) <- colnames(mnewarray_iPSC_Control_ICR_counted2rationavg3minuscontrol)
my_colour5 = list(Annotations = c("Control"= "darkgrey","CasepG"= "#EC7063","CasepR"= "#922B21", "Corrclone13"="#48C9B0", "Corrclone7"="#148F77", "Corrclone50"="#52BE80", "Corrclone35"="#1E8449"))

#Remove VTRNA2
breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(mnewarray_iPSC_Control_ICR_counted2rationavg3minuscontrol[c(1:36,38:42),],
         color = colorRampPalette(c("navy", "white", "#D47400"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour5,
         fontsize = 8,
         annotation_col = my_sample_col5,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cutree_cols = 3)

#save as pHeatmap_mnewarray_iPSC_Control_ICR_counted2rationavg.normcontrol.svg
#save as pHeatmap_mnewarray_iPSC_Control_ICR_counted2rationavg.normcontrol.png



############ Other plots with full data comine #############
head(icomdataICR1_counted2,1)
ViopubcomIndvCR <- data.frame(icomdataICR1_counted2[,c(5,8,9:17,6,7,1:4)])
head(ViopubcomIndvCR)
write.table(ViopubcomIndvCR, "ViopubcomIndvCR.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
svg(filename="boxplot_ViopubcomIndvICR.svg", width=18, height=5, pointsize=12)
boxplot(ViopubcomIndvCR, main="Average methylation at human ICRs", xlab="Individuals", ylab="Avg Methylation", ylim= c(0,1), col= c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","red","red","green","green","green","green"))
dev.off()
dim(ViopubcomIndvCR)
ViopubcomIndvCR <- data.frame(ViopubcomIndvCR)
ViopubcomIndvCR1 <- ViopubcomIndvCR[,c(1:11,12,15,17,13:14,16)]
head(ViopubcomIndvCR1,1)
ViopubcomIndvCR2 <- stack(ViopubcomIndvCR1)
head(ViopubcomIndvCR2)
colnames(ViopubcomIndvCR2) <- c("Methylation", "ipubcomdatasets")
head(ViopubcomIndvCR2)
ggplot(ViopubcomIndvCR2, aes(x=ipubcomdatasets, y=Methylation, color=ipubcomdatasets, fill=ipubcomdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() +scale_color_manual(values=c("grey","grey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#48C9B0","#52BE80","#922B21","#148F77","#1E8449"))+
  scale_fill_manual(values=c("white","white","white","white","white","white","white","white","white","white","white","white","white","white","white","white","white"))
ggsave("Violin_plot_ViopubcomIndvCR2_pimethhed_humanICR.svg", width=25*1.25, height=8*1.25, units="cm", dpi=96)
ggsave("Violin_plot_ViopubcomIndvCR2_pimethhed_humanICR.png", width=25*1.25, height=8*1.25, units="cm", dpi=96)


print('Do the same with median controls')
##########
#No need to take filtered CpGs files as I want to evaluate all CpGs present inside Imprinted regions or PCDH
#And compare it to all CpGs 
system("bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed > MERGE_myiCombat3reavg_human_ICR.txt")

#MERGE_myiCombat3reavg_pos_chr <- read.table("MERGE_myiCombat3reavg_pos_chr.txt", header = F)
MERGE_myiCombat3reavg_human_ICR <- read.table("MERGE_myiCombat3reavg_human_ICR.txt", header = F)
head(MERGE_myiCombat3reavg_human_ICR)
dim(MERGE_myiCombat3reavg_human_ICR)
colnames(MERGE_myiCombat3reavg_human_ICR) <- c(colnames(MERGE_myiCombat3reavg_pos_chr),"chr","start","end","DMR", "DMRType")
#PCA with  Human ICR
iAvgcomdata1 <- MERGE_myiCombat3reavg_human_ICR[,c((ncol(MERGE_myiCombat3reavg_human_ICR)-1),ncol(MERGE_myiCombat3reavg_human_ICR), 1:(ncol(MERGE_myiCombat3reavg_human_ICR)-5))]
rownames(iAvgcomdata1)
colnames(iAvgcomdata1)
head(iAvgcomdata1,1)
dim(iAvgcomdata1)
Count_iAvgcomdataICR1 <- count(iAvgcomdata1, "DMR")
head(Count_iAvgcomdataICR1)
dim(Count_iAvgcomdataICR1)
sgicrcomaggregate1 = aggregate(iAvgcomdata1[,7:ncol(iAvgcomdata1)],by=list(iAvgcomdata1$DMR), mean)
head(sgicrcomaggregate1, 2)
#Plot PCA by group color and labelling
iAvgcomdataICR1=sgicrcomaggregate1
rownames(iAvgcomdataICR1)
iAvgcomdataICR1[,1]
rownames(iAvgcomdataICR1)=iAvgcomdataICR1[,1]
rownames(iAvgcomdataICR1)
colnames(iAvgcomdataICR1)
iAvgcomdataICR1 = iAvgcomdataICR1[,-1]
head(iAvgcomdataICR1)
dim(iAvgcomdataICR1)
iAvgcomdataICR1_counted <- cbind(iAvgcomdataICR1, Count_iAvgcomdataICR1)
head(iAvgcomdataICR1_counted)
dim(iAvgcomdataICR1_counted)
iAvgcomdataICR1_counted1 <- iAvgcomdataICR1_counted[which(iAvgcomdataICR1_counted$freq >= 3),]
head(iAvgcomdataICR1_counted1)
dim(iAvgcomdataICR1_counted1)
write.table(iAvgcomdataICR1_counted1, "iAvgcomdataICR1_counted1.filtmin3.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)

iAvgcomdataICR1_counted2 <- iAvgcomdataICR1_counted1[,1:(ncol(iAvgcomdataICR1_counted1)-2)]
iAvgcomdataICR1_counted3 <- iAvgcomdataICR1_counted2[,-c(2:12)]
head(iAvgcomdataICR1_counted3)
summary(iAvgcomdataICR1_counted3)
#df <- as.iAvgcomdata.frame(iAvgcomdataICR1)
iAvgcomdfICR <- iAvgcomdataICR1_counted3
head(iAvgcomdfICR)
dim(iAvgcomdfICR)
iAvgcomdfICR = data.frame(iAvgcomdfICR)
#write.table(inormdfICR , "inormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
tiAvgcomdfICR = t(iAvgcomdfICR)
tiAvgcomdfICR = data.frame(tiAvgcomdfICR)
#write.table(tinormdfICR , "tinormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tiAvgcomdfICR)

adjPdavg <-  read.table("adjPdavg.txt", header =T, stringsAsFactors = F)

sample_infoavg <- adjPdavg


tiAvgcomdfICR["Color"] <- adjPdavg$Sample_Group
#write.table(tinormdfICR , "tinormdfICRdedup.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tiAvgcomdfICR)
iAvgcomdfx <-tiAvgcomdfICR[c(1:nrow(iAvgcomdfICR))]
GicomC<-prcomp(iAvgcomdfx, center = TRUE, scale. = TRUE)
GicomCi<-data.frame(GicomC$x,Color=tiAvgcomdfICR$Color)
percentagegomICR <- round(GicomC$sdev^2 / sum(GicomC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentagegomICR <- paste( colnames(GicomCi), "(", paste( as.character(percentagegomICR), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

gcom1 <-ggplot(GicomCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentagegomICR[1]) + ylab(percentagegomICR[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("darkgrey","red","green"))+
  scale_shape_manual(values= c(0,1,2))
gcom1 <- gcom1 +theme_classic()
gcom1 + xlim(-20,20)+ ylim(-20,20)
#p + xlim(-1000,1000)+ ylim(-1000,1000)
ggsave("PCA_tiAvgcomdfICR_CntrlIndiv.svg", width=15*1.25, height=10*1.25, units="cm", dpi=96)

################################### Heatmap_indiv human ICR ######################################
#Count number of probes present in the final list: i.e. Inside 43 DMRs
head(iAvgcomdataICR1_counted1)
dim(iAvgcomdataICR1_counted1)
sum(iAvgcomdataICR1_counted1$freq)

#Heatmap_indiv with 11 Controls Median Filtered by minimum 3 CpGs

iAvgcomdataICR2 = as.matrix(iAvgcomdataICR1_counted2)
head(iAvgcomdataICR2)
dim(iAvgcomdataICR2)
colfunc <- colorRampPalette(c("navy","white", "darkred"))
iAvgcomdataICR3 <- iAvgcomdataICR2
head(iAvgcomdataICR3)
write.table(iAvgcomdataICR3, "Heatmap_indiv_iAvgcomdataICR3_sorted_counted.txt" , quote = FALSE, append = FALSE, sep = "\t", row.names = TRUE)
head(iAvgcomdataICR1_counted2)

newavgrray_iPSC_Control_ICR <- iAvgcomdataICR1_counted3
my_sample_col9 <- data.frame(Annotations= c("Control","Case","Case","Corrclone","Corrclone","Corrclone", "Corrclone"))
row.names(my_sample_col9) <- colnames(newavgrray_iPSC_Control_ICR)
my_colour9 = list(Annotations = c("Control"= "darkgrey","Case"= "red", "Corrclone"="green"))

breaksList9 = seq(0, 1, by = 0.01)

pheatmap(newavgrray_iPSC_Control_ICR,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList9)),
         breaks = breaksList9,
         annotation_colors = my_colour9,
         fontsize = 8,
         annotation_col = my_sample_col9,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 3)
#save as pheatmap_iAvgcomdataICR1_counted2ratio.svg
head(newavgrray_iPSC_Control_ICR,1)

wnewavgrray_iPSC_Control_ICR_counted2rationavg1 <- newavgrray_iPSC_Control_ICR
wnewavgrray_iPSC_Control_ICR_counted2rationavg2 <- as.matrix(wnewavgrray_iPSC_Control_ICR_counted2rationavg1)
head(wnewavgrray_iPSC_Control_ICR_counted2rationavg2)
dim(wnewavgrray_iPSC_Control_ICR_counted2rationavg2)
mnewavgrray_iPSC_Control_ICR_counted2rationavg3 <- data.frame(cbind((wnewavgrray_iPSC_Control_ICR_counted2rationavg2[,1:7]-wnewavgrray_iPSC_Control_ICR_counted2rationavg2[,1])))


head(mnewavgrray_iPSC_Control_ICR_counted2rationavg3)
mnewavgrray_iPSC_Control_ICR_counted2rationavg3 <- as.matrix(mnewavgrray_iPSC_Control_ICR_counted2rationavg3)
dim(mnewavgrray_iPSC_Control_ICR_counted2rationavg3)
mnewavgrray_iPSC_Control_ICR_counted2rationavg3minuscontrol <- mnewavgrray_iPSC_Control_ICR_counted2rationavg3[,c(1:7)]
head(mnewavgrray_iPSC_Control_ICR_counted2rationavg3minuscontrol,1)
my_sample_col11 <- data.frame(Annotations= c("Allontrol","Case","Case","Corrclone","Corrclone","Corrclone", "Corrclone"))
row.names(my_sample_col11) <- colnames(mnewavgrray_iPSC_Control_ICR_counted2rationavg3minuscontrol)
my_colour11 = list(Annotations = c("Allontrol"= "darkgrey","Case"= "#FFB6B3", "Corrclone"="#BDE7BD"))


breaksListx = seq(-0.5, 0.5, by = 0.01)
pheatmap(mnewavgrray_iPSC_Control_ICR_counted2rationavg3minuscontrol,
         color = colorRampPalette(c("navy", "white", "#D47400"))(length(breaksListx)),
         breaks = breaksListx,
         annotation_colors = my_colour11,
         fontsize = 8,
         annotation_col = my_sample_col11,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cutree_cols = 3)

#save as pHeatmap_mnewavgrray_iPSC_Control_ICR_counted2rationavg.normcontrol.svg
#save as pHeatmap_mnewavgrray_iPSC_Control_ICR_counted2rationavg.normcontrol.png
head(iAvgcomdataICR1_counted3)
VioavgpubcomIndvCR <- data.frame(iAvgcomdataICR1_counted3)
VioavgpubcomIndvCR1 <- VioavgpubcomIndvCR[,c(1,2,5,7,3,4,6)]
head(VioavgpubcomIndvCR1,1)
VioavgpubcomIndvCR2 <- stack(VioavgpubcomIndvCR1)
head(VioavgpubcomIndvCR2)
colnames(VioavgpubcomIndvCR2) <- c("Methylation", "ipubcomdatasets")
head(VioavgpubcomIndvCR2)
ggplot(VioavgpubcomIndvCR2, aes(x=ipubcomdatasets, y=Methylation, color=ipubcomdatasets, fill=ipubcomdatasets)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_bw() +scale_fill_manual(values=c("white","white","white","white","white","white","white"))+scale_color_manual(values=c("darkgrey","#EC7063","#48C9B0","#52BE80","#922B21","#148F77","#1E8449"))
ggsave("Violin_plot_VioavgpubcomIndvCR2_pimethhed_humanICR.svg", width=15*1.25, height=8*1.25, units="cm", dpi=96)
ggsave("Violin_plot_VioavgpubcomIndvCR2_pimethhed_humanICR.png", width=15*1.25, height=8*1.25, units="cm", dpi=96)





# Extract Imprinted CpGs
head(iAvgcomdata1)
dim(iAvgcomdata1)

#All Imprinted region CpGs
write.table(data.frame(iAvgcomdata1$TargetID), "iAvgcomdata1_TargetID.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
#system("fsystem("grep -f iAvgcomdata1_TargetID.txt MERGE_myipubCombat3re_pos_chravg1.txt -w > MERGE_myipubCombat3re_pos_chravg1_imp.txt")
iAvgcomdata1_TargetID <- data.frame(iAvgcomdata1$TargetID)
colnames(iAvgcomdata1_TargetID) <- "TargetID"
MERGE_myipubCombat3re_pos_chravg1_imp <-  merge(iAvgcomdata1_TargetID, MERGE_myipubCombat3re_pos_chravg1, by.x = "TargetID", by.y="row_4")
dim(MERGE_myipubCombat3re_pos_chravg1_imp)
#----------Identification of ICF-iPSCs specific CpGs (both hypo and hyper) in ICRs-----------#
dim(miMERGE_myipubCombat3re_chrpos_avg_PG)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR)
miMERGE_myipubCombat3re_chrpos_avg_PG_ICRs <-  merge(iAvgcomdata1, miMERGE_myipubCombat3re_chrpos_avg_PG, by.x = "TargetID", by.y="row_4")
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_ICRs)
miMERGE_myipubCombat3re_chrpos_avg_PR_ICRs <-  merge(iAvgcomdata1, miMERGE_myipubCombat3re_chrpos_avg_PR, by.x = "TargetID", by.y="row_4")
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_ICRs)


#Convert beta values to M-values and then apply Mean +/- 3SD test
#Later compare DMRs list obtained from individual dmCpGs and from average beta for ICRs
#Find genome-wide differentially methylated CpGs
dim(iAvgcomdataICR1_counted1)
head(iAvgcomdataICR1_counted1)

iAvgcomdataICR1_counted1_Betaval <- iAvgcomdataICR1_counted1[,2:18]
head(iAvgcomdataICR1_counted1_Betaval)
#Convert Beta to M-values
iAvgcomdataICR1_counted1_Mval <- lumi::beta2m(iAvgcomdataICR1_counted1_Betaval)
head(iAvgcomdataICR1_counted1_Mval)
summary(iAvgcomdataICR1_counted1_Mval)
colnames(iAvgcomdataICR1_counted1_Mval) <- paste0(colnames(iAvgcomdataICR1_counted1_Mval), "_Mval")

iAvgcomdataICR1_counted1_Mvalre <- iAvgcomdataICR1_counted1_Mval
head(iAvgcomdataICR1_counted1_Mvalre,3)

iAvgcomdataICR1_counted1_Mvalre <- as.matrix(iAvgcomdataICR1_counted1_Mvalre)
iAvgcomdataICR1_counted1_Mvalavg <- data.frame(cbind(
  (rowMeans(iAvgcomdataICR1_counted1_Mvalre[,1:11])),
  3 * (rowSds(iAvgcomdataICR1_counted1_Mvalre[,1:11])),
  (iAvgcomdataICR1_counted1_Mvalre[,1:17])
))

#Check if the functions rowMeans and rowSds works correctly
icrtestrowfval <- data.frame(t(iAvgcomdataICR1_counted1_Mvalre[1:2,1:11]))
colnames(icrtestrowfval) <- c("val1","val2")
mean(c(icrtestrowfval$val1))
icrtestrowfvalsd <- sd(c(icrtestrowfval$val1)) 
icrtestrowfvalsd * 3
mean(c(icrtestrowfval$val2))
icrtestrowfvalsd2 <- sd(c(icrtestrowfval$val2)) 
icrtestrowfvalsd2 * 3

#So both function match with manual check

head(iAvgcomdataICR1_counted1_Mvalavg,2)
colnames(iAvgcomdataICR1_counted1_Mvalavg) <- c("Meancontrol_Mval","Mval_3Sdcontrol",colnames(iAvgcomdataICR1_counted1_Mvalre))
dim(iAvgcomdataICR1_counted1_Mvalavg)
head(iAvgcomdataICR1_counted1_Mvalavg)

#Subtraction, Control -Control
iAvgcomdataICR1_counted1_Mvalavg["Con_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["D250_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$D250_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["UN_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$UN_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control1_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control1_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control2_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control2_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control3_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control3_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control4_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control4_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control5_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control5_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control6_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control6_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control7_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control7_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control8_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control8_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["iPSC_control9_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$iPSC_control9_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)

#Subtraction, Case -Control
iAvgcomdataICR1_counted1_Mvalavg["PG_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$PG_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["PR_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$PR_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["C7_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$C7_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["C13_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$C13_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["C35_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$C35_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["C50_Con_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$C50_Mval - iAvgcomdataICR1_counted1_Mvalavg$Meancontrol_Mval)
iAvgcomdataICR1_counted1_Mvalavg["C7_PR_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$C7_Mval - iAvgcomdataICR1_counted1_Mvalavg$PR_Mval)
iAvgcomdataICR1_counted1_Mvalavg["C13_PG_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$C13_Mval - iAvgcomdataICR1_counted1_Mvalavg$PG_Mval)
iAvgcomdataICR1_counted1_Mvalavg["C35_PR_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$C35_Mval - iAvgcomdataICR1_counted1_Mvalavg$PR_Mval)
iAvgcomdataICR1_counted1_Mvalavg["C50_PG_Mval"] <- abs(iAvgcomdataICR1_counted1_Mvalavg$C50_Mval - iAvgcomdataICR1_counted1_Mvalavg$PG_Mval)

iAvgcomdataICR1_counted1_Mvalavg$DMR <- rownames(iAvgcomdataICR1_counted1_Mvalavg)
head(iAvgcomdataICR1_counted1_Mvalavg)
dim(iAvgcomdataICR1_counted1_Mvalavg)

#Process beta value 
iAvgcomdataICR1_counted1_specificity <- iAvgcomdataICR1_counted1
head(iAvgcomdataICR1_counted1_specificity)
#Subtraction, Control - Control oe Sample - Control
iAvgcomdataICR1_counted1_specificity["Con_Con"] <- iAvgcomdataICR1_counted1_specificity$Allcontrol - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["D250_Con"] <- iAvgcomdataICR1_counted1_specificity$D250 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["UN_Con"] <- iAvgcomdataICR1_counted1_specificity$UN - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control1_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control1 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control2_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control2 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control3_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control3 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control4_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control4 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control5_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control5 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control6_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control6 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control7_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control7 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control8_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control8 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["iPSC_control9_Con"] <- iAvgcomdataICR1_counted1_specificity$iPSC_control9 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["PG_Con"] <- iAvgcomdataICR1_counted1_specificity$PG - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["PR_Con"] <- iAvgcomdataICR1_counted1_specificity$PR - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["C7_Con"] <- iAvgcomdataICR1_counted1_specificity$C7 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["C13_Con"] <- iAvgcomdataICR1_counted1_specificity$C13 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["C35_Con"] <- iAvgcomdataICR1_counted1_specificity$C35 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["C50_Con"] <- iAvgcomdataICR1_counted1_specificity$C50 - iAvgcomdataICR1_counted1_specificity$Allcontrol
iAvgcomdataICR1_counted1_specificity["C7_PR"] <- iAvgcomdataICR1_counted1_specificity$C7 - iAvgcomdataICR1_counted1_specificity$PR
iAvgcomdataICR1_counted1_specificity["C13_PG"] <- iAvgcomdataICR1_counted1_specificity$C13 - iAvgcomdataICR1_counted1_specificity$PG
iAvgcomdataICR1_counted1_specificity["C35_PR"] <- iAvgcomdataICR1_counted1_specificity$C35 - iAvgcomdataICR1_counted1_specificity$PR
iAvgcomdataICR1_counted1_specificity["C50_PG"] <- iAvgcomdataICR1_counted1_specificity$C50 - iAvgcomdataICR1_counted1_specificity$PG

write.table(iAvgcomdataICR1_counted1_specificity, "iAvgcomdataICR1_counted1_specificity.txt", sep="\t", quote = F, append = F, row.names = F)
head(iAvgcomdataICR1_counted1_specificity)

#Sample specific diff meth CpGs
piAvgcomdataICR1_counted1_specificity2 <- merge(iAvgcomdataICR1_counted1_specificity, iAvgcomdataICR1_counted1_Mvalavg, by.x="DMR",by.y="DMR", all.x=T) 


dim(piAvgcomdataICR1_counted1_specificity2)
head(piAvgcomdataICR1_counted1_specificity2)

#Apply the formula for  abs(meanControlMval - sampleMval) >= 3SD

piAvgcomdataICR1_counted1_specificity2_PGdev <- piAvgcomdataICR1_counted1_specificity2[which(piAvgcomdataICR1_counted1_specificity2$PG_Con_Mval >= piAvgcomdataICR1_counted1_specificity2$Mval_3Sdcontrol),]
dim(piAvgcomdataICR1_counted1_specificity2_PGdev)

#Now filter by beta value difference

#PG specific diff meth CpGs +/- 0.2
piAvgcomdataICR1_counted1_specificity_PG <- piAvgcomdataICR1_counted1_specificity2_PGdev[which(piAvgcomdataICR1_counted1_specificity2_PGdev$PG_Con > 0.20 | 
                                                                                                 piAvgcomdataICR1_counted1_specificity2_PGdev$PG_Con < -0.20),]



miiAvgcomdataICR1_counted1_specificity_PG <- piAvgcomdataICR1_counted1_specificity_PG[which(piAvgcomdataICR1_counted1_specificity_PG$D250_Con < 0.20 & 
                                                                                              piAvgcomdataICR1_counted1_specificity_PG$D250_Con > -0.20),]


miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$UN_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$UN_Con > -0.20),]


miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control1_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control1_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control2_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control2_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control3_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control3_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control4_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control4_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control5_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control5_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control6_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control6_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control7_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control7_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control8_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control8_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PG <- miiAvgcomdataICR1_counted1_specificity_PG[which(miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control9_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PG$iPSC_control9_Con > -0.20),]

head(miiAvgcomdataICR1_counted1_specificity_PG)
dim(miiAvgcomdataICR1_counted1_specificity_PG)
unique(miiAvgcomdataICR1_counted1_specificity_PG$DMR)

write.table(miiAvgcomdataICR1_counted1_specificity_PG$DMR, "miiAvgcomdataICR1_counted1_specificity_PG.dmr", sep ="\t", row.names = F, append = F, quote=F, col.names = F)
write.table(miiAvgcomdataICR1_counted1_specificity_PG, "miiAvgcomdataICR1_counted1_specificity_PG.txt", sep ="\t", row.names = T, append = F, quote=F, col.names = T)


#Apply the formula for  abs(meanControlMval - sampleMval) >= 3SD

piAvgcomdataICR1_counted1_specificity2_PRdev <- piAvgcomdataICR1_counted1_specificity2[which(piAvgcomdataICR1_counted1_specificity2$PR_Con_Mval >= piAvgcomdataICR1_counted1_specificity2$Mval_3Sdcontrol),]
dim(piAvgcomdataICR1_counted1_specificity2_PRdev)
#PR specific diff meth CpGs +/- 0.2
piAvgcomdataICR1_counted1_specificity_PR <- piAvgcomdataICR1_counted1_specificity2_PRdev[which(piAvgcomdataICR1_counted1_specificity2_PRdev$PR_Con > 0.20 | 
                                                                                                 piAvgcomdataICR1_counted1_specificity2_PRdev$PR_Con < -0.20),]



miiAvgcomdataICR1_counted1_specificity_PR <- piAvgcomdataICR1_counted1_specificity_PR[which(piAvgcomdataICR1_counted1_specificity_PR$D250_Con < 0.20 & 
                                                                                              piAvgcomdataICR1_counted1_specificity_PR$D250_Con > -0.20),]


miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$UN_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$UN_Con > -0.20),]


miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control1_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control1_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control2_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control2_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control3_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control3_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control4_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control4_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control5_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control5_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control6_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control6_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control7_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control7_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control8_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control8_Con > -0.20),]

miiAvgcomdataICR1_counted1_specificity_PR <- miiAvgcomdataICR1_counted1_specificity_PR[which(miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control9_Con < 0.20 & 
                                                                                               miiAvgcomdataICR1_counted1_specificity_PR$iPSC_control9_Con > -0.20),]

head(miiAvgcomdataICR1_counted1_specificity_PR)
dim(miiAvgcomdataICR1_counted1_specificity_PR)
unique(miiAvgcomdataICR1_counted1_specificity_PR$DMR)
write.table(miiAvgcomdataICR1_counted1_specificity_PR$DMR, "miiAvgcomdataICR1_counted1_specificity_PR.dmr", sep ="\t", row.names = F, append = F, quote=F, col.names = F)
write.table(miiAvgcomdataICR1_counted1_specificity_PR, "miiAvgcomdataICR1_counted1_specificity_PR.txt", sep ="\t", row.names = T, append = F, quote=F, col.names = T)

#Compare single CpGs
unique(miMERGE_myipubCombat3re_chrpos_avg_PG_ICRs$DMR)
unique(miMERGE_myipubCombat3re_chrpos_avg_PR_ICRs$DMR)

#Overlap differentially methylated DMR with differentially methylated CpGs
system("wc -l miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt")
system("wc -l miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper_chr.txt")
system("wc -l miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt")
system("wc -l miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper_chr.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed | awk '{print $90}' | sort -k1,1 -u > miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onICR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper_chr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed | awk '{print $90}' | sort -k1,1 -u > miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper_chr_onICR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed | awk '{print $90}' | sort -k1,1 -u > miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onICR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper_chr.txt -b /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed | awk '{print $90}' | sort -k1,1 -u > miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper_chr_onICR.txt")
#These extra DMRs can be put in line chart
#PG specific imprinted region CpGs
#Keep ZNF597_tss out as it is hypermethylated
diff_imp_loci_pG <- as.character(unique(miiAvgcomdataICR1_counted1_specificity_PG$DMR))
write.table(diff_imp_loci_pG,"diff_imp_loci_pG.txt", sep = "\t",append = F, quote = F, row.names = F, col.names = F)
diff_imp_loci_pG_Hypo <- diff_imp_loci_pG[!diff_imp_loci_pG %in% c("ZNF597_tss")]
diff_imp_loci_pG_Hypo <- data.frame(diff_imp_loci_pG_Hypo)
write.table(diff_imp_loci_pG_Hypo,"diff_imp_loci_pG_Hypo.txt", sep = "\t",append = F, quote = F, row.names = F, col.names = F)

diff_imp_loci_pG_Hypo <- read.table("diff_imp_loci_pG_Hypo.txt")
colnames(diff_imp_loci_pG_Hypo) <- "diff_imp_loci_pG"
#Get CpGs present inside dmrs which hypometh in pG
diff_imp_loci_pG_Hypo_CpGs <-  merge(iAvgcomdata1, diff_imp_loci_pG_Hypo, by.x = "DMR", by.y="diff_imp_loci_pG")
dim(diff_imp_loci_pG_Hypo_CpGs) #104
pG_specific_imp_Hypo_CpGs <- data.frame(diff_imp_loci_pG_Hypo_CpGs$TargetID)
colnames(pG_specific_imp_Hypo_CpGs) <- c("TargetID")
dim(pG_specific_imp_Hypo_CpGs)
length(unique(pG_specific_imp_Hypo_CpGs$TargetID))
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)
miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs <-  merge(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo, pG_specific_imp_Hypo_CpGs, by.x = "row_4", by.y="TargetID")
head(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs)
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs) 
cat("Out of Total ",dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)[1]  ," DMPs (Hypometh) in matrixPG and ",dim(pG_specific_imp_Hypo_CpGs)[1]," in PG_specific__Hypo_CpGs only ",dim(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs)[1]," overlapped.")
cat("So ",dim(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs)[1]," Hypo_CpGs are affected with |B|>0.2 and ", length(pG_specific_imp_Hypo_CpGs$TargetID)- length(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$row_4), " are not affected with this criteria")

#PR specific imprinted region CpGs
#Keep ZNF597_tss out as it is hypermethylated
diff_imp_loci_pR <- as.character(unique(miiAvgcomdataICR1_counted1_specificity_PR$DMR))
write.table(diff_imp_loci_pR,"diff_imp_loci_pR.txt", sep = "\t",append = F, quote = F, row.names = F, col.names = F)
diff_imp_loci_pR_Hypo <- diff_imp_loci_pR[!diff_imp_loci_pR %in% c("ZNF597_tss")]
diff_imp_loci_pR_Hypo <- data.frame(diff_imp_loci_pR_Hypo)
write.table(diff_imp_loci_pR_Hypo,"diff_imp_loci_pR_Hypo.txt", sep = "\t",append = F, quote = F, row.names = F, col.names = F)

diff_imp_loci_pR_Hypo <- read.table("diff_imp_loci_pR_Hypo.txt")
colnames(diff_imp_loci_pR_Hypo) <- "diff_imp_loci_pR"
#Get CpGs present inside dmrs which affected in pR
diff_imp_loci_pR_Hypo_CpGs <-  merge(iAvgcomdata1, diff_imp_loci_pR_Hypo, by.x = "DMR", by.y="diff_imp_loci_pR")
dim(diff_imp_loci_pR_Hypo_CpGs)
pR_specific_imp_Hypo_CpGs <- data.frame(diff_imp_loci_pR_Hypo_CpGs$TargetID)
colnames(pR_specific_imp_Hypo_CpGs) <- c("TargetID")
dim(pR_specific_imp_Hypo_CpGs)
length(unique(pR_specific_imp_Hypo_CpGs$TargetID))
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)
miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs <-  merge(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo, pR_specific_imp_Hypo_CpGs, by.x = "row_4", by.y="TargetID")
head(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs)
cat("Out of",dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)[1]  ," DMPs in matrixPR and ",dim(pR_specific_imp_Hypo_CpGs)[1]," in PR_specific__Hypo_CpGs only ",dim(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs)[1]," overlapped. ")
cat("So ",dim(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs)[1]," Hypo_CpGs are affected with |B|>0.2 and ", length(pR_specific_imp_Hypo_CpGs$TargetID)- length(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$row_4), " are not affected with this criteria")


MERGE_myipubCombat3re_pos_chravg1_imp <-  merge(iAvgcomdata1_TargetID, MERGE_myipubCombat3re_pos_chravg1, by.x = "TargetID", by.y="row_4")
dim(MERGE_myipubCombat3re_pos_chravg1_imp)#Total CpGs covered in ICRs 654
# Extract PCDH CpGs
system("bedtools intersect -wa -wb  -a MERGE_myiCombat3reavg_pos_chr.txt -b PCDHcluster.txt > ipcdhdata1.txt")
ipcdhdata1 <- read.table("ipcdhdata1.txt", header = F, stringsAsFactors = F)
dim(ipcdhdata1)#620 CpGs
length(unique(ipcdhdata1$V4))
#All PCDH region CpGs
write.table(data.frame(ipcdhdata1$V4), "ipcdhdata1_TargetID.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
#system("grep -f ipcdhdata1_TargetID.txt MERGE_myipubCombat3re_pos_chravg1.txt -w > MERGE_myipubCombat3re_pos_chravg1_pcdh.txt")
ipcdhdata1CpGs <- data.frame(ipcdhdata1$V4)
colnames(ipcdhdata1CpGs) <- c("TargetID")
dim(MERGE_myipubCombat3re_pos_chravg1)
#Get any CpGs covered in PCDH cluster
MERGE_myipubCombat3re_pos_chravg1_pcdh <-  merge(ipcdhdata1CpGs, MERGE_myipubCombat3re_pos_chravg1, by.x = "TargetID", by.y="row_4")
dim(MERGE_myipubCombat3re_pos_chravg1_pcdh)

#Get hypomethylated CpGs covered in PCDH cluster
dim(ipcdhdata1CpGs)

#PG specific PCDH region Hypo CpGs
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)
miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs <-  merge(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo, ipcdhdata1CpGs, by.x = "row_4", by.y="TargetID")
head(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs)
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs)
cat("Out of",dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)[1]  ," DMPs in matrixPG and ",dim(ipcdhdata1CpGs)[1]," in PG_specific_CpGs only ",dim(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs)[1]," overlapped. ")


#PR specific PCDH region Hypo CpGs
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)
miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs <-  merge(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo, ipcdhdata1CpGs, by.x = "row_4", by.y="TargetID")
head(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs)
cat("Out of",dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)[1]  ," DMPs in matrixPR and ",dim(ipcdhdata1CpGs)[1]," in PR_specific_CpGs only ",dim(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs)[1]," overlapped. ")

print("#---Mean methylation loss and recovery for all hypomethylated among cases and corrected clones--#")

#----->All_PG
dim(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo)
head(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo)
library(dplyr)
miMean_myipubCombat3re_pos_chravg_PG_Hypo <- select(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo, PG_Con, C13_Con, C50_Con)
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo)
stmiMean_myipubCombat3re_pos_chravg_PG_Hypo <- stack(as.matrix(miMean_myipubCombat3re_pos_chravg_PG_Hypo))
head(stmiMean_myipubCombat3re_pos_chravg_PG_Hypo)
colnames(stmiMean_myipubCombat3re_pos_chravg_PG_Hypo) <- c("region", "sample", "meandelta_meth")
stmiMean_myipubCombat3re_pos_chravg_PG_Hypo <- data.frame(stmiMean_myipubCombat3re_pos_chravg_PG_Hypo)
ggplot(stmiMean_myipubCombat3re_pos_chravg_PG_Hypo, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#EC7063","#48C9B0","#52BE80"))+
  theme_bw()+ ylim(-1,1) + geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myipubCombat3re_pos_chravg_PG_Hypo.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myipubCombat3re_pos_chravg_PG_Hypo.png", width=10, height=8, units="cm", dpi=96)

#----->All_PR
dim(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo)
head(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo)
miMean_myipubCombat3re_pos_chravg_PR_Hypo <- select(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo, PR_Con, C7_Con, C35_Con)
head(miMean_myipubCombat3re_pos_chravg_PR_Hypo)
stmiMean_myipubCombat3re_pos_chravg_PR_Hypo <- stack(as.matrix(miMean_myipubCombat3re_pos_chravg_PR_Hypo))
head(stmiMean_myipubCombat3re_pos_chravg_PR_Hypo)
colnames(stmiMean_myipubCombat3re_pos_chravg_PR_Hypo) <- c("region", "sample", "meandelta_meth")
stmiMean_myipubCombat3re_pos_chravg_PR_Hypo <- data.frame(stmiMean_myipubCombat3re_pos_chravg_PR_Hypo)
ggplot(stmiMean_myipubCombat3re_pos_chravg_PR_Hypo, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#922B21","#148F77","#1E8449"))+
  theme_bw() + ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myipubCombat3re_pos_chravg_PR_Hypo.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myipubCombat3re_pos_chravg_PR_Hypo.png", width=10, height=8, units="cm", dpi=96)


#
#For_Imprinted_regions_take only hypomethylated dmrs as hypermethylated are only minorly up
#------>Imprinted_PG Hypomethylated
miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs <- data.frame(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs)
rownames(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs) <- paste0(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$row_1,
                                                                           "%",
                                                                           miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$row_2,
                                                                           "%",
                                                                           miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$row_3,
                                                                           "%",
                                                                           miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$row_4,
                                                                           "%",
                                                                           "Imprinted")

head(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs)
miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG <- select(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs, PG_Con, C13_Con, C50_Con)
head(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG)
stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG <- stack(as.matrix(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG))
head(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG)
colnames(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG) <- c("region", "sample", "meandelta_meth")
stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG <- data.frame(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG)

ggplot(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#EC7063","#48C9B0","#52BE80"))+
  theme_bw()+ ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG.png", width=10, height=8, units="cm", dpi=96)

#------>Imprinted_PR Hypomethylated
miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs <- data.frame(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs)
rownames(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs) <- paste0(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$row_1,
                                                                           "%",
                                                                           miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$row_2,
                                                                           "%",
                                                                           miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$row_3,
                                                                           "%",
                                                                           miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$row_4,
                                                                           "%",
                                                                           "Imprinted")

head(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs)
miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR <- select(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs, PR_Con, C7_Con, C35_Con)
head(miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR)
stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR <- stack(as.matrix(miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR))
head(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR)
colnames(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR) <- c("region", "sample", "meandelta_meth")
stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR <- data.frame(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR)

ggplot(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#922B21","#148F77","#1E8449"))+
  theme_bw()+ ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR.png", width=10, height=8, units="cm", dpi=96)



#------>PCDH_PG Hypomethylated
miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs <- data.frame(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs)
rownames(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs) <- paste0(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$row_1,
                                                                            "%",
                                                                            miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$row_2,
                                                                            "%",
                                                                            miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$row_3,
                                                                            "%",
                                                                            miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$row_4,
                                                                            "%",
                                                                            "PCDH")

head(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs)
miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG <- select(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs, PG_Con, C13_Con, C50_Con)
head(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG)
stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG <- stack(as.matrix(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG))
head(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG)
colnames(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG) <- c("region", "sample", "meandelta_meth")
stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG <- data.frame(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG)

ggplot(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#EC7063","#48C9B0","#52BE80"))+
  theme_bw()+ ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG.png", width=10, height=8, units="cm", dpi=96)

#------>PCDH_PR Hypomethylated
miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs <- data.frame(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs)
rownames(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs) <- paste0(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$row_1,
                                                                            "%",
                                                                            miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$row_2,
                                                                            "%",
                                                                            miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$row_3,
                                                                            "%",
                                                                            miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$row_4,
                                                                            "%",
                                                                            "PCDH")

head(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs)
miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR <- select(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs, PR_Con, C7_Con, C35_Con)
head(miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR)
stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR <- stack(as.matrix(miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR))
head(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR)
colnames(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR) <- c("region", "sample", "meandelta_meth")
stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR <- data.frame(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR)

ggplot(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR, aes(x=sample, y=meandelta_meth, fill=sample)) + 
  geom_boxplot(position=position_dodge(1)) + scale_fill_manual(values=c("#922B21","#148F77","#1E8449"))+
  theme_bw()+ ylim(-1,1)+ geom_hline(yintercept=0, linetype="dashed",  color = "blue", size=0.5)
ggsave("stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR.svg", width=10, height=8, units="cm", dpi=96)
ggsave("stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR.png", width=10, height=8, units="cm", dpi=96)



#----- Combined Boxplot facet_wrap --------#
#The boxplot represents CpGs diff meth by +/- 0.2 for each feature

#pG-All
miMean_myipubCombat3re_pos_chravg_PG_Hypo_all <- miMean_myipubCombat3re_pos_chravg_PG_Hypo
colnames(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all) <- paste0(colnames(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all), "_All")
miMean_myipubCombat3re_pos_chravg_PG_Hypo_all <- stack(as.matrix(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all))
miMean_myipubCombat3re_pos_chravg_PG_Hypo_all <- data.frame(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all)
miMean_myipubCombat3re_pos_chravg_PG_Hypo_all["Color"] <- "All"
#pG-Imprinted
miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined <- miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG
colnames(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined) <- paste0(colnames(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined), "_diffImprinted")
miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined <- stack(as.matrix(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined))
miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined <- data.frame(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined)
miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined["Color"] <- "diffImprinted"
#pG-PCDH
miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh <- miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG
colnames(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh) <- paste0(colnames(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh), "_PCDH")
miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh <- stack(as.matrix(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh))
miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh <- data.frame(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh)
miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh["Color"] <- "PCDH"
miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH <- rbind.data.frame(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all,
                                                                       miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined,
                                                                       miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh)

head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH)
dim(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH)
ggplot(data=miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH, aes(x=col, y=value, fill=Color)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+ coord_flip()+
  theme_minimal()
head(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH)
ggplot(data=miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH, aes(x=col, y=value)) +
  geom_boxplot(aes(col=col,fill=col),position=position_dodge())+
  ylim(-1,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = rep(c("white","white","white"), 3))+
  scale_color_manual(values = rep(c("#EC7063","#48C9B0","#52BE80"), 3))+
  facet_wrap(~Color, scale="free") +geom_hline(yintercept=0, linetype="dashed", color = "blue")

ggsave("miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH.png", width=17*1.25, height=15*1.25, units="cm", dpi=96)
ggsave("miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH.svg", width=17*1.25, height=15*1.25, units="cm", dpi=96)
write.table(miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH, "miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH.txt", row.names = F, quote = F, append = F, sep="\t")


#pR-All
miMean_myipubCombat3re_pos_chravg_PR_Hypo_all <- miMean_myipubCombat3re_pos_chravg_PR_Hypo
colnames(miMean_myipubCombat3re_pos_chravg_PR_Hypo_all) <- paste0(colnames(miMean_myipubCombat3re_pos_chravg_PR_Hypo_all), "_All")
miMean_myipubCombat3re_pos_chravg_PR_Hypo_all <- stack(as.matrix(miMean_myipubCombat3re_pos_chravg_PR_Hypo_all))
miMean_myipubCombat3re_pos_chravg_PR_Hypo_all <- data.frame(miMean_myipubCombat3re_pos_chravg_PR_Hypo_all)
miMean_myipubCombat3re_pos_chravg_PR_Hypo_all["Color"] <- "All"
#pR-Imprinted
miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined <- miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR
colnames(miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined) <- paste0(colnames(miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined), "_diffImprinted")
miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined <- stack(as.matrix(miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined))
miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined <- data.frame(miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined)
miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined["Color"] <- "diffImprinted"
#pR-PCDH
miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh <- miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR
colnames(miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh) <- paste0(colnames(miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh), "_PCDH")
miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh <- stack(as.matrix(miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh))
miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh <- data.frame(miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh)
miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh["Color"] <- "PCDH"
miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH <- rbind.data.frame(miMean_myipubCombat3re_pos_chravg_PR_Hypo_all,
                                                                       miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined,
                                                                       miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh)

head(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH)
dim(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH)
ggplot(data=miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH, aes(x=col, y=value)) +
  geom_boxplot(aes(col=col,fill=col),position=position_dodge())+
  ylim(-1,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = rep(c("white","white","white"), 3))+
  scale_color_manual(values = rep(c("#922B21","#148F77","#1E8449"), 3))+
  facet_wrap(~Color, scale="free") +geom_hline(yintercept=0, linetype="dashed", color = "blue")

ggsave("miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH.png", width=17*1.25, height=15*1.25, units="cm", dpi=96)
ggsave("miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH.svg", width=17*1.25, height=15*1.25, units="cm", dpi=96)

write.table(miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH, "miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH.txt", row.names = F, quote = F, append = F, sep="\t")



#Statistical test
head(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo)
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$Allcontrol, miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$PG, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$Allcontrol, miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$C13, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$Allcontrol, miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$C50, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$PG, miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$C13, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$PG, miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$C50, paired = TRUE, alternative = "two.sided")


head(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo)
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$Allcontrol, miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$PR, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$Allcontrol, miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$C7, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$Allcontrol, miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$C35, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$PR, miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$C7, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$PR, miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$C35, paired = TRUE, alternative = "two.sided")

head(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs)
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$PG, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$C13, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$C50, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$PG, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$C13, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$PG, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$C50, paired = TRUE, alternative = "two.sided")

head(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs)
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$PR, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$C7, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$C35, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$PR, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$C7, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$PR, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$C35, paired = TRUE, alternative = "two.sided")

head(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs)
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$PG, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$C13, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$C50, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$PG, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$C13, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$PG, miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$C50, paired = TRUE, alternative = "two.sided")

head(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs)
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$PR, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$C7, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$Allcontrol, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$C35, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$PR, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$C7, paired = TRUE, alternative = "two.sided")
wilcox.test(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$PR, miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$C35, paired = TRUE, alternative = "two.sided")

median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$Allcontrol)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$PG))
median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$Allcontrol)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$C13))
median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$Allcontrol)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$C50))
median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$PG)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$C13))
median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$PG)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo$C50))


median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$Allcontrol)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$PR))
median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$Allcontrol)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$C7))
median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$Allcontrol)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$C35))
median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$PR)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$C7))
median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$PR)) - median(c(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo$C35))


median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$PG))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$C13))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$C50))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$PG)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$C13))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$PG)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs$C50))

median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$PR))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$C7))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$C35))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$PR)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$C7))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$PR)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs$C35))

median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$PG))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$C13))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$C50))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$PG)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$C13))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$PG)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs$C50))

median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$PR))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$C7))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$Allcontrol)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$C35))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$PR)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$C7))
median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$PR)) - median(c(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs$C35))


print("#---Boxplot Mean methylation loss and recovery for all hypomethylated among cases and corrected clones--#")


#----->All_PG
dim(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo)
head(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo)
library(dplyr)
boxMean_myipubCombat3re_pos_chravg_PG_Hypo <- select(miMERGE_myipubCombat3re_pos_chravg_PG_Hypo, Allcontrol, PG, C13, C50)
head(boxMean_myipubCombat3re_pos_chravg_PG_Hypo)
stboxMean_myipubCombat3re_pos_chravg_PG_Hypo <- stack(as.matrix(boxMean_myipubCombat3re_pos_chravg_PG_Hypo))
head(stboxMean_myipubCombat3re_pos_chravg_PG_Hypo)
colnames(stboxMean_myipubCombat3re_pos_chravg_PG_Hypo) <- c("region", "sample", "meandelta_meth")
stboxMean_myipubCombat3re_pos_chravg_PG_Hypo <- data.frame(stboxMean_myipubCombat3re_pos_chravg_PG_Hypo)

#----->All_PR
dim(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo)
head(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo)
boxMean_myipubCombat3re_pos_chravg_PR_Hypo <- select(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo, Allcontrol, PR, C7, C35)
head(boxMean_myipubCombat3re_pos_chravg_PR_Hypo)
stboxMean_myipubCombat3re_pos_chravg_PR_Hypo <- stack(as.matrix(boxMean_myipubCombat3re_pos_chravg_PR_Hypo))
head(stboxMean_myipubCombat3re_pos_chravg_PR_Hypo)
colnames(stboxMean_myipubCombat3re_pos_chravg_PR_Hypo) <- c("region", "sample", "meandelta_meth")
stboxMean_myipubCombat3re_pos_chravg_PR_Hypo <- data.frame(stboxMean_myipubCombat3re_pos_chravg_PR_Hypo)


#
#For_Imprinted_regions_take only hypomethylated dmrs as hypermethylated are only minorly up
#------>Imprinted_PG Hypomethylated

head(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs)
boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG <- select(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs, Allcontrol, PG, C13, C50)
head(boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG)
stboxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG <- stack(as.matrix(boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG))
head(stboxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG)
colnames(stboxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG) <- c("region", "sample", "meandelta_meth")
stboxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG <- data.frame(stboxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG)

#------>Imprinted_PR Hypomethylated
head(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs)
boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR <- select(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs,Allcontrol, PR, C7, C35)
head(boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR)
stboxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR <- stack(as.matrix(boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR))
head(stboxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR)
colnames(stboxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR) <- c("region", "sample", "meandelta_meth")
stboxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR <- data.frame(stboxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR)



#------>PCDH_PG Hypomethylated
head(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs)
boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG <- select(miMERGE_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs, Allcontrol, PG, C13, C50)
head(boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG)
stboxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG <- stack(as.matrix(boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG))
head(stboxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG)
colnames(stboxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG) <- c("region", "sample", "meandelta_meth")
stboxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG <- data.frame(stboxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG)


#------>PCDH_PR Hypomethylated
head(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs)
boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR <- select(miMERGE_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs,Allcontrol, PR, C7, C35)
head(boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR)
stboxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR <- stack(as.matrix(boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR))
head(stboxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR)
colnames(stboxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR) <- c("region", "sample", "meandelta_meth")
stboxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR <- data.frame(stboxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR)


#----- Combined Boxplot facet_wrap --------#
#The boxplot represents CpGs diff meth by +/- 0.2 for each feature

#pG-All
boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all <- boxMean_myipubCombat3re_pos_chravg_PG_Hypo
colnames(boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all) <- paste0(colnames(boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all), "_All")
boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all <- stack(as.matrix(boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all))
boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all <- data.frame(boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all)
boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all["Color"] <- "All"
#pG-Imprinted
boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined <- boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG
colnames(boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined) <- paste0(colnames(boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined), "_diffImprinted")
boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined <- stack(as.matrix(boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined))
boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined <- data.frame(boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined)
boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined["Color"] <- "diffImprinted"
#pG-PCDH
boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh <- boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG
colnames(boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh) <- paste0(colnames(boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh), "_PCDH")
boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh <- stack(as.matrix(boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh))
boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh <- data.frame(boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh)
boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh["Color"] <- "PCDH"
boxMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH <- rbind.data.frame(boxMean_myipubCombat3re_pos_chravg_PG_Hypo_all,
                                                                        boxMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined,
                                                                        boxMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh)

head(boxMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH)
dim(boxMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH)
ggplot(data=boxMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH, aes(x=col, y=value)) +
  geom_boxplot(aes(col=col,fill=col),position=position_dodge())+
  ylim(0,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = rep(c("white","white","white","white"), 3))+
  scale_color_manual(values = rep(c("darkgrey", "#EC7063","#48C9B0","#52BE80"), 3))+
  facet_wrap(~Color, scale="free")

ggsave("boxMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH.png", width=17*1.25, height=10*1.25, units="cm", dpi=96)
ggsave("boxMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH.svg", width=17*1.25, height=10*1.25, units="cm", dpi=96)
write.table(boxMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH, "boxMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH.txt", row.names = F, quote = F, append = F, sep="\t")

#pR-All
boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all <- boxMean_myipubCombat3re_pos_chravg_PR_Hypo
colnames(boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all) <- paste0(colnames(boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all), "_All")
boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all <- stack(as.matrix(boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all))
boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all <- data.frame(boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all)
boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all["Color"] <- "All"
#pR-Imprinted
boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined <- boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR
colnames(boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined) <- paste0(colnames(boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined), "_diffImprinted")
boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined <- stack(as.matrix(boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined))
boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined <- data.frame(boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined)
boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined["Color"] <- "diffImprinted"
#pR-PCDH
boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh <- boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR
colnames(boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh) <- paste0(colnames(boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh), "_PCDH")
boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh <- stack(as.matrix(boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh))
boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh <- data.frame(boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh)
boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh["Color"] <- "PCDH"
boxMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH <- rbind.data.frame(boxMean_myipubCombat3re_pos_chravg_PR_Hypo_all,
                                                                        boxMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR_imprined,
                                                                        boxMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR_pcdh)

head(boxMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH)
dim(boxMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH)
ggplot(data=boxMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH, aes(x=col, y=value)) +
  geom_boxplot(aes(col=col,fill=col),position=position_dodge())+
  ylim(0,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = rep(c("white","white","white","white"), 3))+
  scale_color_manual(values = rep(c("darkgrey","#922B21","#148F77","#1E8449"), 3))+
  facet_wrap(~Color, scale="free") 
ggsave("boxMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH.png", width=17*1.25, height=10*1.25, units="cm", dpi=96)
ggsave("boxMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH.svg", width=17*1.25, height=10*1.25, units="cm", dpi=96)

write.table(boxMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH, "boxMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH.txt", row.names = F, quote = F, append = F, sep="\t")


#-------------------- Bar plot con Line trend  ----------------------#
Monk_human_ICR_hg38 <- read.table("/home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed", header=F, stringsAsFactors = F)
colnames(Monk_human_ICR_hg38) <- c("chr","start","end","DMR", "DMRtype")
Monk_ICR_dnmt3b <- read.table("Monk_ICR_dnmt3b_IP_counts.txt", header = T, stringsAsFactors = F)
Monk_ICR_dnmt3bx <- Monk_ICR_dnmt3b[,c(22:28)]
colnames(Monk_ICR_dnmt3bx) <- paste0(c(colnames(Monk_ICR_dnmt3bx)),"_DNMT3b")
Monk_ICR_dnmt3bx["row"] <- rownames(Monk_ICR_dnmt3bx)
Monk_ICR_h3k36me3 <- read.table("Monk_ICR_h3k36me3_IP_counts.txt")
Monk_ICR_h3k36me3x <- Monk_ICR_h3k36me3[,c(22:28)]
colnames(Monk_ICR_h3k36me3x) <- paste0(c(colnames(Monk_ICR_h3k36me3x)),"_K36")
Monk_ICR_h3k36me3x["row"] <- rownames(Monk_ICR_h3k36me3x)
Monk_ICR_h3k4me3 <- read.table("Monk_ICR_h3k4me3_IP_counts.txt")
Monk_ICR_h3k4me3x <- Monk_ICR_h3k4me3[,c(22:28)]
colnames(Monk_ICR_h3k4me3x) <- paste0(c(colnames(Monk_ICR_h3k4me3x)),"_K4")
Monk_ICR_h3k4me3x["row"] <- rownames(Monk_ICR_h3k4me3x)


iAvgcomdataICR1_counted1 <- read.table("iAvgcomdataICR1_counted1.filtmin3.txt", header = T, stringsAsFactors = F)
iAvgcomdataICR1_counted4 <- iAvgcomdataICR1_counted1
iAvgcomdataICR1_counted4["DMR"] <- rownames(iAvgcomdataICR1_counted4)
iAvgcomdataICR1_countedDMR <- merge(iAvgcomdataICR1_counted4, Monk_human_ICR_hg38, by="DMR")
iAvgcomdataICR1_countedDMR["row"] <- iAvgcomdataICR1_countedDMR$DMRtype


icomdataICR1_countedDMR_dnmt3b <- merge(iAvgcomdataICR1_countedDMR, Monk_ICR_dnmt3bx, by="row")
head(icomdataICR1_countedDMR_dnmt3b)
icomdataICR1_countedDMR_dnmt3b_k36 <- merge(icomdataICR1_countedDMR_dnmt3b, Monk_ICR_h3k36me3x, by="row")
icomdataICR1_countedDMR_dnmt3b_k36_k4 <- merge(icomdataICR1_countedDMR_dnmt3b_k36, Monk_ICR_h3k4me3x, by="row")
head(icomdataICR1_countedDMR_dnmt3b_k36_k4,1)
dim(icomdataICR1_countedDMR_dnmt3b_k36_k4)
rownames(icomdataICR1_countedDMR_dnmt3b_k36_k4)  <- icomdataICR1_countedDMR_dnmt3b_k36_k4$row
head(icomdataICR1_countedDMR_dnmt3b_k36_k4,1)
write.table(icomdataICR1_countedDMR_dnmt3b_k36_k4,"icomdataICR1_countedDMR_dnmt3b_k36_k4.txt",row.names=F,quote=FALSE,append = F,sep = "\t")

system("grep -f diff_imp_loci_pG_ZNF597_tss.txt icomdataICR1_countedDMR_dnmt3b_k36_k4.txt > icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG.txt")
system("grep -f diff_imp_loci_pR_ZNF597_tss.txt icomdataICR1_countedDMR_dnmt3b_k36_k4.txt > icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR.txt")

icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG <- read.table("icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG.txt", header=F, stringsAsFactors = F)
colnames(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG) <- colnames(icomdataICR1_countedDMR_dnmt3b_k36_k4)
rownames(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG) <-  icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG$row
head(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG)
dim(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG)


icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR <- read.table("icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR.txt", header=F, stringsAsFactors = F)
colnames(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR) <- colnames(icomdataICR1_countedDMR_dnmt3b_k36_k4)
rownames(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR) <-  icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR$row
head(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR)
dim(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR)


library(ComplexHeatmap)
library(circlize)
icom_ICR_subsetonlypG <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG[,c(3,15,18,20)])
dnmt3b_subsetonlypG <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG[,c(26,30,31,32)])
k36_subsetonlypG <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG[,c(33,37,38,39)])
k4_subsetonlypG <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipG[,c(40,44,45,46)])

icom_ICR_subsetonlypR <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR[,c(3,16,17,19)])
dnmt3b_subsetonlypR <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR[,c(26,27,28,29)])
k36_subsetonlypR <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR[,c(33,34,35,36)])
k4_subsetonlypR <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4_subset_diffimplocipR[,c(40,41,42,43)])


# Barplot PR
sticom_ICR_subsetonly_PR <- data.frame(stack(as.matrix(icom_ICR_subsetonlypR[,c(1,2:4)])))

sticom_ICR_subsetonly_PR %>%
  ggplot() + 
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) + theme_classic()+
  scale_color_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  scale_fill_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))
ggsave("bar_sticom_ICR_subsetonly_PR.svg", width=15, height=2, units="in", dpi=96)
ggsave("bar_sticom_ICR_subsetonly_PR.png", width=15, height=2, units="in", dpi=96)

# Barplot
stdnmt3b_subsetonly_PR <- data.frame(stack(as.matrix(dnmt3b_subsetonlypR[,c(1,2:4)])))
stdnmt3b_subsetonly_PR %>%
  ggplot() + scale_y_continuous(breaks = seq(0, 3, by = 0.5))+
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +theme_classic()+
  scale_color_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  scale_fill_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))

ggsave("bar_stdnmt3b_subsetonly_PR.svg", width=15, height=2, units="in", dpi=96)
ggsave("bar_stdnmt3b_subsetonly_PR.png", width=15, height=2, units="in", dpi=96)

# Barplot
stk36_subsetonly_PR <- data.frame(stack(as.matrix(k36_subsetonlypR[,c(1,2:4)])))

stk36_subsetonly_PR %>%
  ggplot() +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +theme_classic()+
  scale_color_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  scale_fill_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))

ggsave("bar_stk36_subsetonly_PR.svg", width=15, height=2, units="in", dpi=96)
ggsave("bar_stk36_subsetonly_PR.png", width=15, height=2, units="in", dpi=96)

# Barplot
stk4_subsetonly_PR <- data.frame(stack(as.matrix(k4_subsetonlypR[,c(1,2:4)])))

stk4_subsetonly_PR %>%
  ggplot() +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +theme_classic()+
  scale_color_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  scale_fill_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))

ggsave("bar_stk4_subsetonly_PR.svg", width=15, height=2, units="in", dpi=96)
ggsave("bar_stk4_subsetonly_PR.png", width=15, height=2, units="in", dpi=96)



#pR-icr
icom_ICR_subsetonlypR_icr <- icom_ICR_subsetonlypR[,c(1,2:4)]
colnames(icom_ICR_subsetonlypR_icr) <- paste0(colnames(icom_ICR_subsetonlypR_icr), "_icr")
icom_ICR_subsetonlypR_icr <- stack(as.matrix(icom_ICR_subsetonlypR_icr))
icom_ICR_subsetonlypR_icr <- data.frame(icom_ICR_subsetonlypR_icr)
icom_ICR_subsetonlypR_icr["Color"] <- "aicr"
#pR-dnmt3b
dnmt3b_subsetonlypR_dnmt3b <- dnmt3b_subsetonlypR[,c(1,2:4)]
colnames(dnmt3b_subsetonlypR_dnmt3b) <- paste0(colnames(dnmt3b_subsetonlypR_dnmt3b), "_dnmt3b")
dnmt3b_subsetonlypR_dnmt3b <- stack(as.matrix(dnmt3b_subsetonlypR_dnmt3b))
dnmt3b_subsetonlypR_dnmt3b <- data.frame(dnmt3b_subsetonlypR_dnmt3b)
dnmt3b_subsetonlypR_dnmt3b["Color"] <- "dnmt3b"
#pR-k4
k4_subsetonlypR_k4 <- k4_subsetonlypR[,c(1,2:4)]
colnames(k4_subsetonlypR_k4) <- paste0(colnames(k4_subsetonlypR_k4), "_k4")
k4_subsetonlypR_k4 <- stack(as.matrix(k4_subsetonlypR_k4))
k4_subsetonlypR_k4 <- data.frame(k4_subsetonlypR_k4)
k4_subsetonlypR_k4["Color"] <- "k4"
features_icr_dnmt3b_k4_pR <- rbind.data.frame(icom_ICR_subsetonlypR_icr,
                                              dnmt3b_subsetonlypR_dnmt3b,
                                              k4_subsetonlypR_k4)

head(features_icr_dnmt3b_k4_pR)
dim(features_icr_dnmt3b_k4_pR)
summary(features_icr_dnmt3b_k4_pR)
ggplot(data=features_icr_dnmt3b_k4_pR, aes(x=row, y=value)) +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"),panel.spacing = unit(1, "lines"))+
  scale_fill_manual(values = rep(c("white","white","white","white"), 9))+ggtitle("pR")+
  scale_color_manual(values = rep(c("darkgrey","#922B21","#148F77","#1E8449"), 9))+
  facet_grid(Color~. ,scales = "free") +
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))
ggsave("features_icr_dnmt3b_k4_pR.svg", width=8.04, height=5.41, units="in", dpi=96)

# Barplot PG
sticom_ICR_subsetonly_PG <- data.frame(stack(as.matrix(icom_ICR_subsetonlypG[,c(1,2:4)])))

sticom_ICR_subsetonly_PG %>%
  ggplot() +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +theme_classic()+
  scale_color_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  scale_fill_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))
ggsave("bar_sticom_ICR_subsetonly_PG.svg", width=15, height=2, units="in", dpi=96)
ggsave("bar_sticom_ICR_subsetonly_PG.png", width=15, height=2, units="in", dpi=96)

# Barplot
stdnmt3b_subsetonly_PG <- data.frame(stack(as.matrix(dnmt3b_subsetonlypG[,c(1,2:4)])))
stdnmt3b_subsetonly_PG %>%
  ggplot() + scale_y_continuous(breaks = seq(0, 3, by = 0.5))+
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +theme_classic()+
  scale_color_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  scale_fill_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))

ggsave("bar_stdnmt3b_subsetonly_PG.svg", width=15, height=2, units="in", dpi=96)
ggsave("bar_stdnmt3b_subsetonly_PG.png", width=15, height=2, units="in", dpi=96)

# Barplot
stk36_subsetonly_PG <- data.frame(stack(as.matrix(k36_subsetonlypG[,c(1,2:4)])))

stk36_subsetonly_PG %>%
  ggplot() +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +theme_classic()+
  scale_color_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  scale_fill_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))

ggsave("bar_stk36_subsetonly_PG.svg", width=15, height=2, units="in", dpi=96)
ggsave("bar_stk36_subsetonly_PG.png", width=15, height=2, units="in", dpi=96)

# Barplot
stk4_subsetonly_PG <- data.frame(stack(as.matrix(k4_subsetonlypG[,c(1,2:4)])))

stk4_subsetonly_PG %>%
  ggplot() +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +theme_classic()+
  scale_color_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  scale_fill_manual(values=c("darkgrey","#FFB6B3","#BDE7BD","#8be28b"))+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))

ggsave("bar_stk4_subsetonly_PG.svg", width=15, height=2, units="in", dpi=96)
ggsave("bar_stk4_subsetonly_PG.png", width=15, height=2, units="in", dpi=96)



#pG-icr
icom_ICR_subsetonlypG_icr <- icom_ICR_subsetonlypG[,c(1,2:4)]
colnames(icom_ICR_subsetonlypG_icr) <- paste0(colnames(icom_ICR_subsetonlypG_icr), "_icr")
icom_ICR_subsetonlypG_icr <- stack(as.matrix(icom_ICR_subsetonlypG_icr))
icom_ICR_subsetonlypG_icr <- data.frame(icom_ICR_subsetonlypG_icr)
icom_ICR_subsetonlypG_icr["Color"] <- "aicr"
#pG-dnmt3b
dnmt3b_subsetonlypG_dnmt3b <- dnmt3b_subsetonlypG[,c(1,2:4)]
colnames(dnmt3b_subsetonlypG_dnmt3b) <- paste0(colnames(dnmt3b_subsetonlypG_dnmt3b), "_dnmt3b")
dnmt3b_subsetonlypG_dnmt3b <- stack(as.matrix(dnmt3b_subsetonlypG_dnmt3b))
dnmt3b_subsetonlypG_dnmt3b <- data.frame(dnmt3b_subsetonlypG_dnmt3b)
dnmt3b_subsetonlypG_dnmt3b["Color"] <- "dnmt3b"
#pG-k4
k4_subsetonlypG_k4 <- k4_subsetonlypG[,c(1,2:4)]
colnames(k4_subsetonlypG_k4) <- paste0(colnames(k4_subsetonlypG_k4), "_k4")
k4_subsetonlypG_k4 <- stack(as.matrix(k4_subsetonlypG_k4))
k4_subsetonlypG_k4 <- data.frame(k4_subsetonlypG_k4)
k4_subsetonlypG_k4["Color"] <- "k4"
features_icr_dnmt3b_k4_pG <- rbind.data.frame(icom_ICR_subsetonlypG_icr,
                                              dnmt3b_subsetonlypG_dnmt3b,
                                              k4_subsetonlypG_k4)

head(features_icr_dnmt3b_k4_pG)
dim(features_icr_dnmt3b_k4_pG)
summary(features_icr_dnmt3b_k4_pG)
ggplot(data=features_icr_dnmt3b_k4_pG, aes(x=row, y=value)) +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"),panel.spacing = unit(1, "lines"))+
  scale_fill_manual(values = rep(c("white","white","white","white"), 9))+ggtitle("pG")+
  scale_color_manual(values = rep(c("darkgrey","#EC7063","#48C9B0","#52BE80"), 9))+
  facet_grid(Color~. ,scales = "free") +geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))

ggsave("features_icr_dnmt3b_k4_pG.svg", width=8.04, height=5.41, units="in", dpi=96)

#All ICR profile
icom_ICR_allpG <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4[,c(3,15,18,20)])
dnmt3b_allpG <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4[,c(26,30,31,32)])
k36_allpG <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4[,c(33,37,38,39)])
k4_allpG <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4[,c(40,44,45,46)])

icom_ICR_allpR <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4[,c(3,16,17,19)])
dnmt3b_allpR <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4[,c(26,27,28,29)])
k36_allpR <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4[,c(33,34,35,36)])
k4_allpR <- as.matrix(icomdataICR1_countedDMR_dnmt3b_k36_k4[,c(40,41,42,43)])



#pR-icr
icom_ICR_allpR_icr <- icom_ICR_allpR[,c(1,2:4)]
colnames(icom_ICR_allpR_icr) <- paste0(colnames(icom_ICR_allpR_icr), "_icr")
icom_ICR_allpR_icr <- stack(as.matrix(icom_ICR_allpR_icr))
icom_ICR_allpR_icr <- data.frame(icom_ICR_allpR_icr)
icom_ICR_allpR_icr["Color"] <- "aicr"
#pR-dnmt3b
dnmt3b_allpR_dnmt3b <- dnmt3b_allpR[,c(1,2:4)]
colnames(dnmt3b_allpR_dnmt3b) <- paste0(colnames(dnmt3b_allpR_dnmt3b), "_dnmt3b")
dnmt3b_allpR_dnmt3b <- stack(as.matrix(dnmt3b_allpR_dnmt3b))
dnmt3b_allpR_dnmt3b <- data.frame(dnmt3b_allpR_dnmt3b)
dnmt3b_allpR_dnmt3b["Color"] <- "dnmt3b"
#pR-k4
k4_allpR_k4 <- k4_allpR[,c(1,2:4)]
colnames(k4_allpR_k4) <- paste0(colnames(k4_allpR_k4), "_k4")
k4_allpR_k4 <- stack(as.matrix(k4_allpR_k4))
k4_allpR_k4 <- data.frame(k4_allpR_k4)
k4_allpR_k4["Color"] <- "k4"
featuresall_icr_dnmt3b_k4_pR <- rbind.data.frame(icom_ICR_allpR_icr,
                                                 dnmt3b_allpR_dnmt3b,
                                                 k4_allpR_k4)

head(featuresall_icr_dnmt3b_k4_pR)
dim(featuresall_icr_dnmt3b_k4_pR)
summary(featuresall_icr_dnmt3b_k4_pR)
ggplot(data=featuresall_icr_dnmt3b_k4_pR, aes(x=row, y=value)) +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"),panel.spacing = unit(1, "lines"))+
  scale_fill_manual(values = rep(c("white","white","white","white"), 9))+ggtitle("pR")+
  scale_color_manual(values = rep(c("darkgrey","#922B21","#148F77","#1E8449"), 9))+
  facet_grid(Color~. ,scales = "free") +
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))
ggsave("featuresall_icr_dnmt3b_k4_pR.svg", width=40.04, height=5.41, units="in", dpi=96)
ggsave("featuresall_icr_dnmt3b_k4_pR.png", width=20.04, height=5.41, units="in", dpi=396)



#pG-icr
icom_ICR_allpG_icr <- icom_ICR_allpG[,c(1,2:4)]
colnames(icom_ICR_allpG_icr) <- paste0(colnames(icom_ICR_allpG_icr), "_icr")
icom_ICR_allpG_icr <- stack(as.matrix(icom_ICR_allpG_icr))
icom_ICR_allpG_icr <- data.frame(icom_ICR_allpG_icr)
icom_ICR_allpG_icr["Color"] <- "aicr"
#pG-dnmt3b
dnmt3b_allpG_dnmt3b <- dnmt3b_allpG[,c(1,2:4)]
colnames(dnmt3b_allpG_dnmt3b) <- paste0(colnames(dnmt3b_allpG_dnmt3b), "_dnmt3b")
dnmt3b_allpG_dnmt3b <- stack(as.matrix(dnmt3b_allpG_dnmt3b))
dnmt3b_allpG_dnmt3b <- data.frame(dnmt3b_allpG_dnmt3b)
dnmt3b_allpG_dnmt3b["Color"] <- "dnmt3b"
#pG-k4
k4_allpG_k4 <- k4_allpG[,c(1,2:4)]
colnames(k4_allpG_k4) <- paste0(colnames(k4_allpG_k4), "_k4")
k4_allpG_k4 <- stack(as.matrix(k4_allpG_k4))
k4_allpG_k4 <- data.frame(k4_allpG_k4)
k4_allpG_k4["Color"] <- "k4"
featuresall_icr_dnmt3b_k4_pG <- rbind.data.frame(icom_ICR_allpG_icr,
                                                 dnmt3b_allpG_dnmt3b,
                                                 k4_allpG_k4)

head(featuresall_icr_dnmt3b_k4_pG)
dim(featuresall_icr_dnmt3b_k4_pG)
summary(featuresall_icr_dnmt3b_k4_pG)
ggplot(data=featuresall_icr_dnmt3b_k4_pG, aes(x=row, y=value)) +
  geom_bar(aes(x=row, y=value, size=0.1,fill =col,col=col),stat = "identity",size =0.8, width = 0.4, position = position_dodge(0.8)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"),panel.spacing = unit(1, "lines"))+
  scale_fill_manual(values = rep(c("white","white","white","white"), 9))+ggtitle("pG")+
  scale_color_manual(values = rep(c("darkgrey","#EC7063","#48C9B0","#52BE80"), 9))+
  facet_grid(Color~. ,scales = "free") +geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))

ggsave("featuresall_icr_dnmt3b_k4_pG.svg", width=40.04, height=5.41, units="in", dpi=96)
ggsave("featuresall_icr_dnmt3b_k4_pG.png", width=20.04, height=5.41, units="in", dpi=396)

#Stats
corArray_DNMT3B_pR <- rbind.data.frame(cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][1,]), c(dnmt3b_subsetonlypR[,c(1,2:4)][1,])),
                                       cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][2,]), c(dnmt3b_subsetonlypR[,c(1,2:4)][2,])),
                                       cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][3,]), c(dnmt3b_subsetonlypR[,c(1,2:4)][3,])),
                                       cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][4,]), c(dnmt3b_subsetonlypR[,c(1,2:4)][4,])),
                                       cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][5,]), c(dnmt3b_subsetonlypR[,c(1,2:4)][5,])),
                                       cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][6,]), c(dnmt3b_subsetonlypR[,c(1,2:4)][6,])),
                                       cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][7,]), c(dnmt3b_subsetonlypR[,c(1,2:4)][7,])))

corArray_K36me_pR <- rbind.data.frame(cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][1,]), c(k36_subsetonlypR[,c(1,2:4)][1,])),
                                      cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][2,]), c(k36_subsetonlypR[,c(1,2:4)][2,])),
                                      cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][3,]), c(k36_subsetonlypR[,c(1,2:4)][3,])),
                                      cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][4,]), c(k36_subsetonlypR[,c(1,2:4)][4,])),
                                      cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][5,]), c(k36_subsetonlypR[,c(1,2:4)][5,])),
                                      cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][6,]), c(k36_subsetonlypR[,c(1,2:4)][6,])),
                                      cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][7,]), c(k36_subsetonlypR[,c(1,2:4)][7,])))

corArray_K4me_pR <- rbind.data.frame(cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][1,]), c(k4_subsetonlypR[,c(1,2:4)][1,])),
                                     cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][2,]), c(k4_subsetonlypR[,c(1,2:4)][2,])),
                                     cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][3,]), c(k4_subsetonlypR[,c(1,2:4)][3,])),
                                     cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][4,]), c(k4_subsetonlypR[,c(1,2:4)][4,])),
                                     cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][5,]), c(k4_subsetonlypR[,c(1,2:4)][5,])),
                                     cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][6,]), c(k4_subsetonlypR[,c(1,2:4)][6,])),
                                     cor(c(icom_ICR_subsetonlypR[,c(1,2:4)][7,]), c(k4_subsetonlypR[,c(1,2:4)][7,])))

corArray_DNMT3B_K36me_K4me_pR <- cbind.data.frame(corArray_DNMT3B_pR,corArray_K36me_pR,corArray_K4me_pR)
colnames(corArray_DNMT3B_K36me_K4me_pR) <- c("corArray_DNMT3B_pR","corArray_K36me_pR","corArray_K4me_pR")
rownames(corArray_DNMT3B_K36me_K4me_pR) <- rownames(icom_ICR_subsetonlypR)
pheatmap(corArray_DNMT3B_K36me_K4me_pR)

breaksListd = seq(-1, 1, by = 0.001)
pheatmap(as.matrix(corArray_DNMT3B_K36me_K4me_pR),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListd)),
         breaks = breaksListd,
         fontsize = 8,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white",
         gaps_row = c(1:6),
         gaps_col = c(1:3),
         cellwidth = 6)
#export as scale_corArray_DNMT3B_K36me_K4me_pR.svg width 354, height 241, Then zoom into area ma closeup and screen capture shift+stamp and paste in fig


#Stats
corArray_DNMT3B_pG <- rbind.data.frame(cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][1,]), c(dnmt3b_subsetonlypG[,c(1,2:4)][1,])),
                                       cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][2,]), c(dnmt3b_subsetonlypG[,c(1,2:4)][2,])),
                                       cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][3,]), c(dnmt3b_subsetonlypG[,c(1,2:4)][3,])),
                                       cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][4,]), c(dnmt3b_subsetonlypG[,c(1,2:4)][4,])),
                                       cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][5,]), c(dnmt3b_subsetonlypG[,c(1,2:4)][5,])),
                                       cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][6,]), c(dnmt3b_subsetonlypG[,c(1,2:4)][6,])),
                                       cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][7,]), c(dnmt3b_subsetonlypG[,c(1,2:4)][7,])),
                                       cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][8,]), c(dnmt3b_subsetonlypG[,c(1,2:4)][8,])))

corArray_K36me_pG <- rbind.data.frame(cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][1,]), c(k36_subsetonlypG[,c(1,2:4)][1,])),
                                      cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][2,]), c(k36_subsetonlypG[,c(1,2:4)][2,])),
                                      cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][3,]), c(k36_subsetonlypG[,c(1,2:4)][3,])),
                                      cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][4,]), c(k36_subsetonlypG[,c(1,2:4)][4,])),
                                      cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][5,]), c(k36_subsetonlypG[,c(1,2:4)][5,])),
                                      cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][6,]), c(k36_subsetonlypG[,c(1,2:4)][6,])),
                                      cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][7,]), c(k36_subsetonlypG[,c(1,2:4)][7,])),
                                      cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][8,]), c(k36_subsetonlypG[,c(1,2:4)][8,])))

corArray_K4me_pG <- rbind.data.frame(cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][1,]), c(k4_subsetonlypG[,c(1,2:4)][1,])),
                                     cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][2,]), c(k4_subsetonlypG[,c(1,2:4)][2,])),
                                     cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][3,]), c(k4_subsetonlypG[,c(1,2:4)][3,])),
                                     cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][4,]), c(k4_subsetonlypG[,c(1,2:4)][4,])),
                                     cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][5,]), c(k4_subsetonlypG[,c(1,2:4)][5,])),
                                     cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][6,]), c(k4_subsetonlypG[,c(1,2:4)][6,])),
                                     cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][7,]), c(k4_subsetonlypG[,c(1,2:4)][7,])),
                                     cor(c(icom_ICR_subsetonlypG[,c(1,2:4)][8,]), c(k4_subsetonlypG[,c(1,2:4)][8,])))

corArray_DNMT3B_K36me_K4me_pG <- cbind.data.frame(corArray_DNMT3B_pG,corArray_K36me_pG,corArray_K4me_pG)
colnames(corArray_DNMT3B_K36me_K4me_pG) <- c("corArray_DNMT3B_pG","corArray_K36me_pG","corArray_K4me_pG")
rownames(corArray_DNMT3B_K36me_K4me_pG) <- rownames(icom_ICR_subsetonlypG)
pheatmap(corArray_DNMT3B_K36me_K4me_pG)

breaksListd = seq(-1, 1, by = 0.001)
pheatmap(as.matrix(corArray_DNMT3B_K36me_K4me_pG),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksListd)),
         breaks = breaksListd,
         fontsize = 8,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white",
         gaps_row = c(1:7),
         gaps_col = c(1:3),
         cellwidth = 6)
#export as scale_corArray_DNMT3B_K36me_K4me_pG.svg width 354, height 256, Then zoom into area ma closeup and screen capture shift+stamp and paste in fig


#EXPORT TABLE OF CORRELATION
head(corArray_DNMT3B_K36me_K4me_pR)
corArray_DNMT3B_K36me_K4me_pR_tab <- data.frame(stack(as.matrix(corArray_DNMT3B_K36me_K4me_pR)))
head(corArray_DNMT3B_K36me_K4me_pG)
corArray_DNMT3B_K36me_K4me_pG_tab <- data.frame(stack(as.matrix(corArray_DNMT3B_K36me_K4me_pG)))

corArray_DNMT3B_K36me_K4me_pR_pG_tab <- rbind.data.frame(corArray_DNMT3B_K36me_K4me_pR_tab,
                                                         corArray_DNMT3B_K36me_K4me_pG_tab)

##Suppl. Table
write.table(corArray_DNMT3B_K36me_K4me_pR_pG_tab, "corArray_DNMT3B_K36me_K4me_pR_pG_tab.txt", sep = "\t", quote = F, row.names = F)
#ucsc_partial
#Prepare files for ucsc_partial:
head(MERGE_myiCombat3reavg_pos_chr)
dim(MERGE_myiCombat3reavg_pos_chr)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,5)], "ucsc_partial.Allcontrol.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,6)], "ucsc_partial.D250.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,7)], "ucsc_partial.UN.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,8)], "ucsc_partial.iPSC_control1.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,9)], "ucsc_partial.iPSC_control2.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,10)], "ucsc_partial.iPSC_control3.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,11)], "ucsc_partial.iPSC_control4.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,12)], "ucsc_partial.iPSC_control5.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,13)], "ucsc_partial.iPSC_control6.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,14)], "ucsc_partial.iPSC_control7.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,15)], "ucsc_partial.iPSC_control8.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,16)], "ucsc_partial.iPSC_control9.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,17)], "ucsc_partial.PG.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,18)], "ucsc_partial.PR.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,19)], "ucsc_partial.C7.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,20)], "ucsc_partial.C13.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,21)], "ucsc_partial.C35.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)
write.table(MERGE_myiCombat3reavg_pos_chr[,c(1:3,22)], "ucsc_partial.C50.hg38.chr.bed",sep="\t", quote = F,append = F,row.names = F, col.names = F)



sed '1s/^/track type=bedGraph name="D250_beta_values_partial" description="D250_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.D250.hg38.chr.bed > ucsc_partial.D250.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="UN_beta_values_partial" description="UN_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.UN.hg38.chr.bed > ucsc_partial.UN.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control1_beta_values_partial" description="iPSC_control1_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control1.hg38.chr.bed > ucsc_partial.iPSC_control1.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control2_beta_values_partial" description="iPSC_control2_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control2.hg38.chr.bed > ucsc_partial.iPSC_control2.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control3_beta_values_partial" description="iPSC_control3_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control3.hg38.chr.bed > ucsc_partial.iPSC_control3.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control4_beta_values_partial" description="iPSC_control4_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control4.hg38.chr.bed > ucsc_partial.iPSC_control4.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control5_beta_values_partial" description="iPSC_control5_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control5.hg38.chr.bed > ucsc_partial.iPSC_control5.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control6_beta_values_partial" description="iPSC_control6_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control6.hg38.chr.bed > ucsc_partial.iPSC_control6.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control7_beta_values_partial" description="iPSC_control7_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control7.hg38.chr.bed > ucsc_partial.iPSC_control7.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control8_beta_values_partial" description="iPSC_control8_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control8.hg38.chr.bed > ucsc_partial.iPSC_control8.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="iPSC_control9_beta_values_partial" description="iPSC_control9_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=0,0,0 \n/' ucsc_partial.iPSC_control9.hg38.chr.bed > ucsc_partial.iPSC_control9.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="PG_beta_values_partial" description="PG_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=255,0,0 \n/' ucsc_partial.PG.hg38.chr.bed > ucsc_partial.PG.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="PR_beta_values_partial" description="PR_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=255,0,0 \n/' ucsc_partial.PR.hg38.chr.bed > ucsc_partial.PR.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="C7_beta_values_partial" description="C7_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=51,255,51 \n/' ucsc_partial.C7.hg38.chr.bed > ucsc_partial.C7.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="C13_beta_values_partial" description="C13_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=51,255,51 \n/' ucsc_partial.C13.hg38.chr.bed > ucsc_partial.C13.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="C35_beta_values_partial" description="C35_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=51,255,51 \n/' ucsc_partial.C35.hg38.chr.bed > ucsc_partial.C35.hg38.chr.bedGraph
sed '1s/^/track type=bedGraph name="C50_beta_values_partial" description="C50_beta_values_partial" visibility=full viewLimits=0:1 autoScale=off  maxHeightPixels=50:33:11  color=51,255,51 \n/' ucsc_partial.C50.hg38.chr.bed > ucsc_partial.C50.hg38.chr.bedGraph

#For custom bed tracks
grep TNXB mpipubcomdataBin500bp1_counted2ratioavg2_PGchr_genehg38.txt | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"pG","0",".",$2,$3,"236,112,99"}'
grep TNXB mpipubcomdataBin500bp1_counted2ratioavg2_PRchr_genehg38.txt | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"pR","0",".",$2,$3,"146,43,33"}'
grep PCDH mpipubcomdataBin500bp1_counted2ratioavg2_PGchr_genehg38.txt | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"pG","0",".",$2,$3,"236,112,99"}'
grep PCDH mpipubcomdataBin500bp1_counted2ratioavg2_PRchr_genehg38.txt | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"pR","0",".",$2,$3,"146,43,33"}'





#Assign diffmeth CpGs to genomic elements
dim(miMERGE_myipubCombat3re_chrpos_avg_PG)
miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper <- miMERGE_myipubCombat3re_chrpos_avg_PG[which(miMERGE_myipubCombat3re_chrpos_avg_PG$PG_Con > 0.2),]
head(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper)
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper)
miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo <- miMERGE_myipubCombat3re_chrpos_avg_PG[which(miMERGE_myipubCombat3re_chrpos_avg_PG$PG_Con < -0.2),]
head(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)

dim(miMERGE_myipubCombat3re_chrpos_avg_PR)
miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper <- miMERGE_myipubCombat3re_chrpos_avg_PR[which(miMERGE_myipubCombat3re_chrpos_avg_PR$PR_Con > 0.2),]
head(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper)
miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo <- miMERGE_myipubCombat3re_chrpos_avg_PR[which(miMERGE_myipubCombat3re_chrpos_avg_PR$PR_Con < -0.2),]
head(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)


#Annotation to regions
library(annotatr)
annots_basic_gelements = c('hg38_cpgs','hg38_basicgenes', 'hg38_genes_intergenic','hg38_genes_cds')

annotations_b_hg38_el = build_annotations(genome = 'hg38', annotations = annots_basic_gelements)


annotations_b_hg38_eldf <- data.frame(annotations_b_hg38_el)
annotations_b_hg38_eldf <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$start > 0),]
annotations_b_hg38_eldf <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$end > 0),]

head(annotations_b_hg38_eldf)
dim(annotations_b_hg38_eldf)
write.table(data.frame(annotations_b_hg38_eldf), "annotations_b_hg38_eldf.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

unique(annotations_b_hg38_eldf$type)

annotations_b_hg38_eldf_promoter <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_genes_promoters"),]
write.table(data.frame(annotations_b_hg38_eldf_promoter), "annotations_b_hg38_eldf_promoter.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_CGI <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_cpg_islands"),]
write.table(data.frame(annotations_b_hg38_eldf_CGI), "annotations_b_hg38_eldf_CGI.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_CGshores <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_cpg_shores"),]
write.table(data.frame(annotations_b_hg38_eldf_CGshores), "annotations_b_hg38_eldf_CGshores.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_CGshelves <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_cpg_shelves"),]
write.table(data.frame(annotations_b_hg38_eldf_CGshelves), "annotations_b_hg38_eldf_CGshelves.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_CGinter <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_cpg_inter"),]
write.table(data.frame(annotations_b_hg38_eldf_CGinter), "annotations_b_hg38_eldf_CGinter.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_cds <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_genes_cds"),]
write.table(data.frame(annotations_b_hg38_eldf_cds), "annotations_b_hg38_eldf_cds.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_5UTR <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_genes_5UTRs"),]
write.table(data.frame(annotations_b_hg38_eldf_5UTR), "annotations_b_hg38_eldf_5UTR.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_3UTR <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_genes_3UTRs"),]
write.table(data.frame(annotations_b_hg38_eldf_3UTR), "annotations_b_hg38_eldf_3UTR.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_intergenic <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_genes_intergenic"),]
write.table(data.frame(annotations_b_hg38_eldf_intergenic), "annotations_b_hg38_eldf_intergenic.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_exon <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_genes_exons"),]
write.table(data.frame(annotations_b_hg38_eldf_exon), "annotations_b_hg38_eldf_exon.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_eldf_intron <- annotations_b_hg38_eldf[which(annotations_b_hg38_eldf$type == "hg38_genes_introns"),]
write.table(data.frame(annotations_b_hg38_eldf_intron), "annotations_b_hg38_eldf_intron.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)


#Repeats
rmsk <- read.table("rmsk.txt", header = F)
head(rmsk)
rmsk <- rmsk[,c(6:8,10,12)]
rmsk
colnames(rmsk) <- c("chr", "start","end", "strand", "repeats")

system("cat miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper_chr.txt > miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt")
system("cat miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper_chr.txt > miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt")

miMERGE_myipubCombat3re_chrpos_avg_PG <- read.table("miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt")
miMERGE_myipubCombat3re_chrpos_avg_PR <- read.table("miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt")
dim(miMERGE_myipubCombat3re_chrpos_avg_PG)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR)
colnames(miMERGE_myipubCombat3re_chrpos_avg_PG) <- colnames(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)
colnames(miMERGE_myipubCombat3re_chrpos_avg_PR) <- colnames(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)

#PG
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_promoter.txt > miMERGE_myipubCombat3re_avg_PG_promoter.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_CGI.txt > miMERGE_myipubCombat3re_avg_PG_CGI.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_CGshores.txt > miMERGE_myipubCombat3re_avg_PG_CGshores.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_CGshelves.txt > miMERGE_myipubCombat3re_avg_PG_CGshelves.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_CGinter.txt > miMERGE_myipubCombat3re_avg_PG_CGinter.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_cds.txt > miMERGE_myipubCombat3re_avg_PG_cds.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_5UTR.txt > miMERGE_myipubCombat3re_avg_PG_5UTR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_3UTR.txt > miMERGE_myipubCombat3re_avg_PG_3UTR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_intergenic.txt > miMERGE_myipubCombat3re_avg_PG_intergenic.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_exon.txt > miMERGE_myipubCombat3re_avg_PG_exon.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_eldf_intron.txt > miMERGE_myipubCombat3re_avg_PG_intron.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_rmsk_LINE.txt > miMERGE_myipubCombat3re_avg_PG_LINE.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_rmsk_LTR.txt > miMERGE_myipubCombat3re_avg_PG_LTR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_rmsk_SINE.txt > miMERGE_myipubCombat3re_avg_PG_SINE.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_chr.txt -b annotations_b_hg38_rmsk_Satellite.txt > miMERGE_myipubCombat3re_avg_PG_Satellite.txt")

miMERGE_myipubCombat3re_avg_PG_promoter <- read.table("miMERGE_myipubCombat3re_avg_PG_promoter.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_promoter) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_promoter)
miMERGE_myipubCombat3re_avg_PG_promoter <- miMERGE_myipubCombat3re_avg_PG_promoter[!duplicated(miMERGE_myipubCombat3re_avg_PG_promoter$row_4),]
miMERGE_myipubCombat3re_avg_PG_promoter_Hypo <- miMERGE_myipubCombat3re_avg_PG_promoter[which(miMERGE_myipubCombat3re_avg_PG_promoter$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_promoter_Hyper <- miMERGE_myipubCombat3re_avg_PG_promoter[which(miMERGE_myipubCombat3re_avg_PG_promoter$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_promoter_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_promoter_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_promoter_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_promoter_Hypo)
mean_myipubCombat3re_avg_PG_promoter_Hypo["promoter"] <- "promoter"
colnames(mean_myipubCombat3re_avg_PG_promoter_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_promoter_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_promoter_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_promoter_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_promoter_Hyper)
mean_myipubCombat3re_avg_PG_promoter_Hyper["promoter"] <- "promoter"
colnames(mean_myipubCombat3re_avg_PG_promoter_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PG_CGI <- read.table("miMERGE_myipubCombat3re_avg_PG_CGI.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_CGI) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_CGI)
miMERGE_myipubCombat3re_avg_PG_CGI <- miMERGE_myipubCombat3re_avg_PG_CGI[!duplicated(miMERGE_myipubCombat3re_avg_PG_CGI$row_4),]
miMERGE_myipubCombat3re_avg_PG_CGI_Hypo <- miMERGE_myipubCombat3re_avg_PG_CGI[which(miMERGE_myipubCombat3re_avg_PG_CGI$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_CGI_Hyper <- miMERGE_myipubCombat3re_avg_PG_CGI[which(miMERGE_myipubCombat3re_avg_PG_CGI$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_CGI_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_CGI_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_CGI_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_CGI_Hypo)
mean_myipubCombat3re_avg_PG_CGI_Hypo["CGI"] <- "CGI"
colnames(mean_myipubCombat3re_avg_PG_CGI_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_CGI_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_CGI_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_CGI_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_CGI_Hyper)
mean_myipubCombat3re_avg_PG_CGI_Hyper["CGI"] <- "CGI"
colnames(mean_myipubCombat3re_avg_PG_CGI_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PG_CGshores <- read.table("miMERGE_myipubCombat3re_avg_PG_CGshores.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_CGshores) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_CGshores)
miMERGE_myipubCombat3re_avg_PG_CGshores <- miMERGE_myipubCombat3re_avg_PG_CGshores[!duplicated(miMERGE_myipubCombat3re_avg_PG_CGshores$row_4),]
miMERGE_myipubCombat3re_avg_PG_CGshores_Hypo <- miMERGE_myipubCombat3re_avg_PG_CGshores[which(miMERGE_myipubCombat3re_avg_PG_CGshores$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_CGshores_Hyper <- miMERGE_myipubCombat3re_avg_PG_CGshores[which(miMERGE_myipubCombat3re_avg_PG_CGshores$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_CGshores_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_CGshores_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_CGshores_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_CGshores_Hypo)
mean_myipubCombat3re_avg_PG_CGshores_Hypo["CGshores"] <- "CGshores"
colnames(mean_myipubCombat3re_avg_PG_CGshores_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_CGshores_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_CGshores_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_CGshores_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_CGshores_Hyper)
mean_myipubCombat3re_avg_PG_CGshores_Hyper["CGshores"] <- "CGshores"
colnames(mean_myipubCombat3re_avg_PG_CGshores_Hyper) <- c("Values","Sample","Feature")


miMERGE_myipubCombat3re_avg_PG_CGshelves <- read.table("miMERGE_myipubCombat3re_avg_PG_CGshelves.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_CGshelves) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_CGshelves)
miMERGE_myipubCombat3re_avg_PG_CGshelves <- miMERGE_myipubCombat3re_avg_PG_CGshelves[!duplicated(miMERGE_myipubCombat3re_avg_PG_CGshelves$row_4),]
miMERGE_myipubCombat3re_avg_PG_CGshelves_Hypo <- miMERGE_myipubCombat3re_avg_PG_CGshelves[which(miMERGE_myipubCombat3re_avg_PG_CGshelves$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_CGshelves_Hyper <- miMERGE_myipubCombat3re_avg_PG_CGshelves[which(miMERGE_myipubCombat3re_avg_PG_CGshelves$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_CGshelves_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_CGshelves_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_CGshelves_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_CGshelves_Hypo)
mean_myipubCombat3re_avg_PG_CGshelves_Hypo["CGshelves"] <- "CGshelves"
colnames(mean_myipubCombat3re_avg_PG_CGshelves_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_CGshelves_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_CGshelves_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_CGshelves_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_CGshelves_Hyper)
mean_myipubCombat3re_avg_PG_CGshelves_Hyper["CGshelves"] <- "CGshelves"
colnames(mean_myipubCombat3re_avg_PG_CGshelves_Hyper) <- c("Values","Sample","Feature")



miMERGE_myipubCombat3re_avg_PG_CGinter <- read.table("miMERGE_myipubCombat3re_avg_PG_CGinter.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_CGinter) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_CGinter)
miMERGE_myipubCombat3re_avg_PG_CGinter <- miMERGE_myipubCombat3re_avg_PG_CGinter[!duplicated(miMERGE_myipubCombat3re_avg_PG_CGinter$row_4),]
miMERGE_myipubCombat3re_avg_PG_CGinter_Hypo <- miMERGE_myipubCombat3re_avg_PG_CGinter[which(miMERGE_myipubCombat3re_avg_PG_CGinter$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_CGinter_Hyper <- miMERGE_myipubCombat3re_avg_PG_CGinter[which(miMERGE_myipubCombat3re_avg_PG_CGinter$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_CGinter_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_CGinter_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_CGinter_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_CGinter_Hypo)
mean_myipubCombat3re_avg_PG_CGinter_Hypo["CGinter"] <- "CGinter"
colnames(mean_myipubCombat3re_avg_PG_CGinter_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_CGinter_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_CGinter_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_CGinter_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_CGinter_Hyper)
mean_myipubCombat3re_avg_PG_CGinter_Hyper["CGinter"] <- "CGinter"
colnames(mean_myipubCombat3re_avg_PG_CGinter_Hyper) <- c("Values","Sample","Feature")



miMERGE_myipubCombat3re_avg_PG_cds <- read.table("miMERGE_myipubCombat3re_avg_PG_cds.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_cds) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_cds)
miMERGE_myipubCombat3re_avg_PG_cds <- miMERGE_myipubCombat3re_avg_PG_cds[!duplicated(miMERGE_myipubCombat3re_avg_PG_cds$row_4),]
miMERGE_myipubCombat3re_avg_PG_cds_Hypo <- miMERGE_myipubCombat3re_avg_PG_cds[which(miMERGE_myipubCombat3re_avg_PG_cds$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_cds_Hyper <- miMERGE_myipubCombat3re_avg_PG_cds[which(miMERGE_myipubCombat3re_avg_PG_cds$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_cds_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_cds_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_cds_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_cds_Hypo)
mean_myipubCombat3re_avg_PG_cds_Hypo["cds"] <- "cds"
colnames(mean_myipubCombat3re_avg_PG_cds_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_cds_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_cds_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_cds_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_cds_Hyper)
mean_myipubCombat3re_avg_PG_cds_Hyper["cds"] <- "cds"
colnames(mean_myipubCombat3re_avg_PG_cds_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PG_5UTR <- read.table("miMERGE_myipubCombat3re_avg_PG_5UTR.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_5UTR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_5UTR)
miMERGE_myipubCombat3re_avg_PG_5UTR <- miMERGE_myipubCombat3re_avg_PG_5UTR[!duplicated(miMERGE_myipubCombat3re_avg_PG_5UTR$row_4),]
miMERGE_myipubCombat3re_avg_PG_5UTR_Hypo <- miMERGE_myipubCombat3re_avg_PG_5UTR[which(miMERGE_myipubCombat3re_avg_PG_5UTR$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_5UTR_Hyper <- miMERGE_myipubCombat3re_avg_PG_5UTR[which(miMERGE_myipubCombat3re_avg_PG_5UTR$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_5UTR_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_5UTR_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_5UTR_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_5UTR_Hypo)
mean_myipubCombat3re_avg_PG_5UTR_Hypo["5UTR"] <- "5UTR"
colnames(mean_myipubCombat3re_avg_PG_5UTR_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_5UTR_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_5UTR_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_5UTR_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_5UTR_Hyper)
mean_myipubCombat3re_avg_PG_5UTR_Hyper["5UTR"] <- "5UTR"
colnames(mean_myipubCombat3re_avg_PG_5UTR_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PG_3UTR <- read.table("miMERGE_myipubCombat3re_avg_PG_3UTR.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_3UTR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_3UTR)
miMERGE_myipubCombat3re_avg_PG_3UTR <- miMERGE_myipubCombat3re_avg_PG_3UTR[!duplicated(miMERGE_myipubCombat3re_avg_PG_3UTR$row_4),]
miMERGE_myipubCombat3re_avg_PG_3UTR_Hypo <- miMERGE_myipubCombat3re_avg_PG_3UTR[which(miMERGE_myipubCombat3re_avg_PG_3UTR$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_3UTR_Hyper <- miMERGE_myipubCombat3re_avg_PG_3UTR[which(miMERGE_myipubCombat3re_avg_PG_3UTR$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_3UTR_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_3UTR_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_3UTR_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_3UTR_Hypo)
mean_myipubCombat3re_avg_PG_3UTR_Hypo["3UTR"] <- "3UTR"
colnames(mean_myipubCombat3re_avg_PG_3UTR_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_3UTR_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_3UTR_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_3UTR_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_3UTR_Hyper)
mean_myipubCombat3re_avg_PG_3UTR_Hyper["3UTR"] <- "3UTR"
colnames(mean_myipubCombat3re_avg_PG_3UTR_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PG_intergenic <- read.table("miMERGE_myipubCombat3re_avg_PG_intergenic.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_intergenic) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_intergenic)
miMERGE_myipubCombat3re_avg_PG_intergenic <- miMERGE_myipubCombat3re_avg_PG_intergenic[!duplicated(miMERGE_myipubCombat3re_avg_PG_intergenic$row_4),]
miMERGE_myipubCombat3re_avg_PG_intergenic_Hypo <- miMERGE_myipubCombat3re_avg_PG_intergenic[which(miMERGE_myipubCombat3re_avg_PG_intergenic$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_intergenic_Hyper <- miMERGE_myipubCombat3re_avg_PG_intergenic[which(miMERGE_myipubCombat3re_avg_PG_intergenic$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_intergenic_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_intergenic_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_intergenic_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_intergenic_Hypo)
mean_myipubCombat3re_avg_PG_intergenic_Hypo["intergenic"] <- "intergenic"
colnames(mean_myipubCombat3re_avg_PG_intergenic_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_intergenic_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_intergenic_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_intergenic_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_intergenic_Hyper)
mean_myipubCombat3re_avg_PG_intergenic_Hyper["intergenic"] <- "intergenic"
colnames(mean_myipubCombat3re_avg_PG_intergenic_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PG_exon <- read.table("miMERGE_myipubCombat3re_avg_PG_exon.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_exon) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_exon)
miMERGE_myipubCombat3re_avg_PG_exon <- miMERGE_myipubCombat3re_avg_PG_exon[!duplicated(miMERGE_myipubCombat3re_avg_PG_exon$row_4),]
miMERGE_myipubCombat3re_avg_PG_exon_Hypo <- miMERGE_myipubCombat3re_avg_PG_exon[which(miMERGE_myipubCombat3re_avg_PG_exon$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_exon_Hyper <- miMERGE_myipubCombat3re_avg_PG_exon[which(miMERGE_myipubCombat3re_avg_PG_exon$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_exon_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_exon_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_exon_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_exon_Hypo)
mean_myipubCombat3re_avg_PG_exon_Hypo["exon"] <- "exon"
colnames(mean_myipubCombat3re_avg_PG_exon_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_exon_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_exon_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_exon_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_exon_Hyper)
mean_myipubCombat3re_avg_PG_exon_Hyper["exon"] <- "exon"
colnames(mean_myipubCombat3re_avg_PG_exon_Hyper) <- c("Values","Sample","Feature")



miMERGE_myipubCombat3re_avg_PG_intron <- read.table("miMERGE_myipubCombat3re_avg_PG_intron.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_intron) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PG_intron)
miMERGE_myipubCombat3re_avg_PG_intron <- miMERGE_myipubCombat3re_avg_PG_intron[!duplicated(miMERGE_myipubCombat3re_avg_PG_intron$row_4),]
miMERGE_myipubCombat3re_avg_PG_intron_Hypo <- miMERGE_myipubCombat3re_avg_PG_intron[which(miMERGE_myipubCombat3re_avg_PG_intron$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_intron_Hyper <- miMERGE_myipubCombat3re_avg_PG_intron[which(miMERGE_myipubCombat3re_avg_PG_intron$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_intron_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_intron_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_intron_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_intron_Hypo)
mean_myipubCombat3re_avg_PG_intron_Hypo["intron"] <- "intron"
colnames(mean_myipubCombat3re_avg_PG_intron_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_intron_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_intron_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_intron_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_intron_Hyper)
mean_myipubCombat3re_avg_PG_intron_Hyper["intron"] <- "intron"
colnames(mean_myipubCombat3re_avg_PG_intron_Hyper) <- c("Values","Sample","Feature")


annotations_b_hg38_rmsk_LINE <- rmsk[which(rmsk$repeats == "LINE"),]
write.table(data.frame(annotations_b_hg38_rmsk_LINE), "annotations_b_hg38_rmsk_LINE.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_rmsk_LTR <- rmsk[which(rmsk$repeats == "LTR"),]
write.table(data.frame(annotations_b_hg38_rmsk_LTR), "annotations_b_hg38_rmsk_LTR.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_rmsk_SINE <- rmsk[which(rmsk$repeats == "SINE"),]
write.table(data.frame(annotations_b_hg38_rmsk_SINE), "annotations_b_hg38_rmsk_SINE.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)
annotations_b_hg38_rmsk_Satellite <- rmsk[which(rmsk$repeats == "Satellite"),]
write.table(data.frame(annotations_b_hg38_rmsk_Satellite), "annotations_b_hg38_rmsk_Satellite.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)


miMERGE_myipubCombat3re_avg_PG_LINE <- read.table("miMERGE_myipubCombat3re_avg_PG_LINE.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_LINE) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(rmsk))
dim(miMERGE_myipubCombat3re_avg_PG_LINE)
miMERGE_myipubCombat3re_avg_PG_LINE <- miMERGE_myipubCombat3re_avg_PG_LINE[!duplicated(miMERGE_myipubCombat3re_avg_PG_LINE$row_4),]
miMERGE_myipubCombat3re_avg_PG_LINE_Hypo <- miMERGE_myipubCombat3re_avg_PG_LINE[which(miMERGE_myipubCombat3re_avg_PG_LINE$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_LINE_Hyper <- miMERGE_myipubCombat3re_avg_PG_LINE[which(miMERGE_myipubCombat3re_avg_PG_LINE$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_LINE_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_LINE_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_LINE_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_LINE_Hypo)
mean_myipubCombat3re_avg_PG_LINE_Hypo["LINE"] <- "LINE"
colnames(mean_myipubCombat3re_avg_PG_LINE_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_LINE_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_LINE_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_LINE_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_LINE_Hyper)
mean_myipubCombat3re_avg_PG_LINE_Hyper["LINE"] <- "LINE"
colnames(mean_myipubCombat3re_avg_PG_LINE_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PG_LTR <- read.table("miMERGE_myipubCombat3re_avg_PG_LTR.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_LTR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(rmsk))
dim(miMERGE_myipubCombat3re_avg_PG_LTR)
miMERGE_myipubCombat3re_avg_PG_LTR <- miMERGE_myipubCombat3re_avg_PG_LTR[!duplicated(miMERGE_myipubCombat3re_avg_PG_LTR$row_4),]
miMERGE_myipubCombat3re_avg_PG_LTR_Hypo <- miMERGE_myipubCombat3re_avg_PG_LTR[which(miMERGE_myipubCombat3re_avg_PG_LTR$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_LTR_Hyper <- miMERGE_myipubCombat3re_avg_PG_LTR[which(miMERGE_myipubCombat3re_avg_PG_LTR$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_LTR_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_LTR_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_LTR_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_LTR_Hypo)
mean_myipubCombat3re_avg_PG_LTR_Hypo["LTR"] <- "LTR"
colnames(mean_myipubCombat3re_avg_PG_LTR_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_LTR_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_LTR_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_LTR_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_LTR_Hyper)
mean_myipubCombat3re_avg_PG_LTR_Hyper["LTR"] <- "LTR"
colnames(mean_myipubCombat3re_avg_PG_LTR_Hyper) <- c("Values","Sample","Feature")



miMERGE_myipubCombat3re_avg_PG_SINE <- read.table("miMERGE_myipubCombat3re_avg_PG_SINE.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_SINE) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(rmsk))
dim(miMERGE_myipubCombat3re_avg_PG_SINE)
miMERGE_myipubCombat3re_avg_PG_SINE <- miMERGE_myipubCombat3re_avg_PG_SINE[!duplicated(miMERGE_myipubCombat3re_avg_PG_SINE$row_4),]
miMERGE_myipubCombat3re_avg_PG_SINE_Hypo <- miMERGE_myipubCombat3re_avg_PG_SINE[which(miMERGE_myipubCombat3re_avg_PG_SINE$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_SINE_Hyper <- miMERGE_myipubCombat3re_avg_PG_SINE[which(miMERGE_myipubCombat3re_avg_PG_SINE$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_SINE_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_SINE_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_SINE_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_SINE_Hypo)
mean_myipubCombat3re_avg_PG_SINE_Hypo["SINE"] <- "SINE"
colnames(mean_myipubCombat3re_avg_PG_SINE_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_SINE_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_SINE_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_SINE_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_SINE_Hyper)
mean_myipubCombat3re_avg_PG_SINE_Hyper["SINE"] <- "SINE"
colnames(mean_myipubCombat3re_avg_PG_SINE_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PG_Satellite <- read.table("miMERGE_myipubCombat3re_avg_PG_Satellite.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PG_Satellite) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG), colnames(rmsk))
dim(miMERGE_myipubCombat3re_avg_PG_Satellite)
miMERGE_myipubCombat3re_avg_PG_Satellite <- miMERGE_myipubCombat3re_avg_PG_Satellite[!duplicated(miMERGE_myipubCombat3re_avg_PG_Satellite$row_4),]
miMERGE_myipubCombat3re_avg_PG_Satellite_Hypo <- miMERGE_myipubCombat3re_avg_PG_Satellite[which(miMERGE_myipubCombat3re_avg_PG_Satellite$PG_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PG_Satellite_Hyper <- miMERGE_myipubCombat3re_avg_PG_Satellite[which(miMERGE_myipubCombat3re_avg_PG_Satellite$PG_Con > 0.2),]
mean_myipubCombat3re_avg_PG_Satellite_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_Satellite_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PG_Satellite_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_Satellite_Hypo)
mean_myipubCombat3re_avg_PG_Satellite_Hypo["Satellite"] <- "Satellite"
colnames(mean_myipubCombat3re_avg_PG_Satellite_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PG_Satellite_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PG_Satellite_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PG_Satellite_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PG_Satellite_Hyper)
mean_myipubCombat3re_avg_PG_Satellite_Hyper["Satellite"] <- "Satellite"
colnames(mean_myipubCombat3re_avg_PG_Satellite_Hyper) <- c("Values","Sample","Feature")


mean_myipubCombat3re_avg_PG_gen_el_Hypo <- rbind.data.frame(mean_myipubCombat3re_avg_PG_promoter_Hypo,
                                                            mean_myipubCombat3re_avg_PG_CGI_Hypo,
                                                            mean_myipubCombat3re_avg_PG_CGshores_Hypo,
                                                            mean_myipubCombat3re_avg_PG_CGshelves_Hypo,
                                                            mean_myipubCombat3re_avg_PG_CGinter_Hypo,
                                                            mean_myipubCombat3re_avg_PG_cds_Hypo,
                                                            mean_myipubCombat3re_avg_PG_5UTR_Hypo,
                                                            mean_myipubCombat3re_avg_PG_3UTR_Hypo,
                                                            mean_myipubCombat3re_avg_PG_intergenic_Hypo,
                                                            mean_myipubCombat3re_avg_PG_exon_Hypo,
                                                            mean_myipubCombat3re_avg_PG_intron_Hypo,
                                                            mean_myipubCombat3re_avg_PG_LINE_Hypo,
                                                            mean_myipubCombat3re_avg_PG_LTR_Hypo,
                                                            mean_myipubCombat3re_avg_PG_SINE_Hypo,
                                                            mean_myipubCombat3re_avg_PG_Satellite_Hypo)

library(tidyr)
tabmean_myipubCombat3re_avg_PG_gen_el_Hypo <- spread(mean_myipubCombat3re_avg_PG_gen_el_Hypo, Sample, Values)
rownames(tabmean_myipubCombat3re_avg_PG_gen_el_Hypo) <- tabmean_myipubCombat3re_avg_PG_gen_el_Hypo$Feature
tabmean_myipubCombat3re_avg_PG_gen_el_Hypo <- tabmean_myipubCombat3re_avg_PG_gen_el_Hypo[,c(1,2,17,18,3,5,6,4)]
rownames(tabmean_myipubCombat3re_avg_PG_gen_el_Hypo) <- paste(c("b","c","g","h","i","j","k","d","e","f","l","m","a","n","o"),"_",tabmean_myipubCombat3re_avg_PG_gen_el_Hypo$Feature, sep="")
tabmean_myipubCombat3re_avg_PG_gen_el_Hypo <- tabmean_myipubCombat3re_avg_PG_gen_el_Hypo[,-1]
tabmean_myipubCombat3re_avg_PG_gen_el_Hypo <- as.matrix(tabmean_myipubCombat3re_avg_PG_gen_el_Hypo)
tabmean_myipubCombat3re_avg_PG_gen_el_Hypo <- tabmean_myipubCombat3re_avg_PG_gen_el_Hypo[order(rownames(tabmean_myipubCombat3re_avg_PG_gen_el_Hypo)),]

pheatmap(tabmean_myipubCombat3re_avg_PG_gen_el_Hypo,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 3)

sticom_gen_el_PG_Hypo <- data.frame(stack(as.matrix(tabmean_myipubCombat3re_avg_PG_gen_el_Hypo)))

sticom_gen_el_PG_Hypo %>%
  ggplot() +
  geom_point(aes(x=row, y=value,colour = col),size=3,alpha=1,stat = "identity", position = position_dodge2(0.8))+theme_classic()+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))+
  scale_color_manual(values=c("darkgrey","#EC7063","#922B21","#48C9B0","#52BE80","#148F77","#1E8449"))+ylim(0,1)+
  scale_fill_manual(values=c("white","white","white","white","white","white","white"))
ggsave("point_sticom_gen_el_PG_Hypo.svg", width=15, height=2, units="in", dpi=96)
ggsave("point_sticom_gen_el_PG_Hypo.png", width=15, height=2, units="in", dpi=96)





mean_myipubCombat3re_avg_PG_gen_el_Hyper <- rbind.data.frame(mean_myipubCombat3re_avg_PG_promoter_Hyper,
                                                             mean_myipubCombat3re_avg_PG_CGI_Hyper,
                                                             mean_myipubCombat3re_avg_PG_CGshores_Hyper,
                                                             mean_myipubCombat3re_avg_PG_CGshelves_Hyper,
                                                             mean_myipubCombat3re_avg_PG_CGinter_Hyper,
                                                             mean_myipubCombat3re_avg_PG_cds_Hyper,
                                                             mean_myipubCombat3re_avg_PG_5UTR_Hyper,
                                                             mean_myipubCombat3re_avg_PG_3UTR_Hyper,
                                                             mean_myipubCombat3re_avg_PG_intergenic_Hyper,
                                                             mean_myipubCombat3re_avg_PG_exon_Hyper,
                                                             mean_myipubCombat3re_avg_PG_intron_Hyper,
                                                             mean_myipubCombat3re_avg_PG_LINE_Hyper,
                                                             mean_myipubCombat3re_avg_PG_LTR_Hyper,
                                                             mean_myipubCombat3re_avg_PG_SINE_Hyper)

tabmean_myipubCombat3re_avg_PG_gen_el_Hyper <- spread(mean_myipubCombat3re_avg_PG_gen_el_Hyper, Sample, Values)
rownames(tabmean_myipubCombat3re_avg_PG_gen_el_Hyper) <- tabmean_myipubCombat3re_avg_PG_gen_el_Hyper$Feature
tabmean_myipubCombat3re_avg_PG_gen_el_Hyper <- tabmean_myipubCombat3re_avg_PG_gen_el_Hyper[,c(1,2,17,18,3,5,6,4)]
rownames(tabmean_myipubCombat3re_avg_PG_gen_el_Hyper) <- paste(c("b","c","g","h","i","j","k","d","e","f","l","m","a","n"),"_",tabmean_myipubCombat3re_avg_PG_gen_el_Hyper$Feature, sep="")
tabmean_myipubCombat3re_avg_PG_gen_el_Hyper <- tabmean_myipubCombat3re_avg_PG_gen_el_Hyper[,-1]
tabmean_myipubCombat3re_avg_PG_gen_el_Hyper <- as.matrix(tabmean_myipubCombat3re_avg_PG_gen_el_Hyper)
tabmean_myipubCombat3re_avg_PG_gen_el_Hyper <- tabmean_myipubCombat3re_avg_PG_gen_el_Hyper[order(rownames(tabmean_myipubCombat3re_avg_PG_gen_el_Hyper)),]

pheatmap(tabmean_myipubCombat3re_avg_PG_gen_el_Hyper,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 3)

sticom_gen_el_PG_Hyper <- data.frame(stack(as.matrix(tabmean_myipubCombat3re_avg_PG_gen_el_Hyper)))

sticom_gen_el_PG_Hyper %>%
  ggplot() +
  geom_point(aes(x=row, y=value,colour = col),size=3,alpha=1,stat = "identity", position = position_dodge2(0.8))+theme_classic()+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))+
  scale_color_manual(values=c("darkgrey","#EC7063","#922B21","#48C9B0","#52BE80","#148F77","#1E8449"))+ylim(0,1)+
  scale_fill_manual(values=c("white","white","white","white","white","white","white"))
ggsave("point_sticom_gen_el_PG_Hyper.svg", width=15, height=2, units="in", dpi=96)
ggsave("point_sticom_gen_el_PG_Hyper.png", width=15, height=2, units="in", dpi=96)



#PR
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_promoter.txt > miMERGE_myipubCombat3re_avg_PR_promoter.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_CGI.txt > miMERGE_myipubCombat3re_avg_PR_CGI.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_CGshores.txt > miMERGE_myipubCombat3re_avg_PR_CGshores.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_CGshelves.txt > miMERGE_myipubCombat3re_avg_PR_CGshelves.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_CGinter.txt > miMERGE_myipubCombat3re_avg_PR_CGinter.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_cds.txt > miMERGE_myipubCombat3re_avg_PR_cds.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_5UTR.txt > miMERGE_myipubCombat3re_avg_PR_5UTR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_3UTR.txt > miMERGE_myipubCombat3re_avg_PR_3UTR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_intergenic.txt > miMERGE_myipubCombat3re_avg_PR_intergenic.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_exon.txt > miMERGE_myipubCombat3re_avg_PR_exon.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_eldf_intron.txt > miMERGE_myipubCombat3re_avg_PR_intron.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_rmsk_LINE.txt > miMERGE_myipubCombat3re_avg_PR_LINE.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_rmsk_LTR.txt > miMERGE_myipubCombat3re_avg_PR_LTR.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_rmsk_SINE.txt > miMERGE_myipubCombat3re_avg_PR_SINE.txt")
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_chr.txt -b annotations_b_hg38_rmsk_Satellite.txt > miMERGE_myipubCombat3re_avg_PR_Satellite.txt")

miMERGE_myipubCombat3re_avg_PR_promoter <- read.table("miMERGE_myipubCombat3re_avg_PR_promoter.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_promoter) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_promoter)
miMERGE_myipubCombat3re_avg_PR_promoter <- miMERGE_myipubCombat3re_avg_PR_promoter[!duplicated(miMERGE_myipubCombat3re_avg_PR_promoter$row_4),]
miMERGE_myipubCombat3re_avg_PR_promoter_Hypo <- miMERGE_myipubCombat3re_avg_PR_promoter[which(miMERGE_myipubCombat3re_avg_PR_promoter$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_promoter_Hyper <- miMERGE_myipubCombat3re_avg_PR_promoter[which(miMERGE_myipubCombat3re_avg_PR_promoter$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_promoter_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_promoter_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_promoter_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_promoter_Hypo)
mean_myipubCombat3re_avg_PR_promoter_Hypo["promoter"] <- "promoter"
colnames(mean_myipubCombat3re_avg_PR_promoter_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_promoter_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_promoter_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_promoter_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_promoter_Hyper)
mean_myipubCombat3re_avg_PR_promoter_Hyper["promoter"] <- "promoter"
colnames(mean_myipubCombat3re_avg_PR_promoter_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_CGI <- read.table("miMERGE_myipubCombat3re_avg_PR_CGI.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_CGI) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_CGI)
miMERGE_myipubCombat3re_avg_PR_CGI <- miMERGE_myipubCombat3re_avg_PR_CGI[!duplicated(miMERGE_myipubCombat3re_avg_PR_CGI$row_4),]
miMERGE_myipubCombat3re_avg_PR_CGI_Hypo <- miMERGE_myipubCombat3re_avg_PR_CGI[which(miMERGE_myipubCombat3re_avg_PR_CGI$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_CGI_Hyper <- miMERGE_myipubCombat3re_avg_PR_CGI[which(miMERGE_myipubCombat3re_avg_PR_CGI$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_CGI_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_CGI_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_CGI_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_CGI_Hypo)
mean_myipubCombat3re_avg_PR_CGI_Hypo["CGI"] <- "CGI"
colnames(mean_myipubCombat3re_avg_PR_CGI_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_CGI_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_CGI_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_CGI_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_CGI_Hyper)
mean_myipubCombat3re_avg_PR_CGI_Hyper["CGI"] <- "CGI"
colnames(mean_myipubCombat3re_avg_PR_CGI_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_CGshores <- read.table("miMERGE_myipubCombat3re_avg_PR_CGshores.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_CGshores) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_CGshores)
miMERGE_myipubCombat3re_avg_PR_CGshores <- miMERGE_myipubCombat3re_avg_PR_CGshores[!duplicated(miMERGE_myipubCombat3re_avg_PR_CGshores$row_4),]
miMERGE_myipubCombat3re_avg_PR_CGshores_Hypo <- miMERGE_myipubCombat3re_avg_PR_CGshores[which(miMERGE_myipubCombat3re_avg_PR_CGshores$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_CGshores_Hyper <- miMERGE_myipubCombat3re_avg_PR_CGshores[which(miMERGE_myipubCombat3re_avg_PR_CGshores$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_CGshores_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_CGshores_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_CGshores_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_CGshores_Hypo)
mean_myipubCombat3re_avg_PR_CGshores_Hypo["CGshores"] <- "CGshores"
colnames(mean_myipubCombat3re_avg_PR_CGshores_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_CGshores_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_CGshores_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_CGshores_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_CGshores_Hyper)
mean_myipubCombat3re_avg_PR_CGshores_Hyper["CGshores"] <- "CGshores"
colnames(mean_myipubCombat3re_avg_PR_CGshores_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_CGshelves <- read.table("miMERGE_myipubCombat3re_avg_PR_CGshelves.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_CGshelves) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_CGshelves)
miMERGE_myipubCombat3re_avg_PR_CGshelves <- miMERGE_myipubCombat3re_avg_PR_CGshelves[!duplicated(miMERGE_myipubCombat3re_avg_PR_CGshelves$row_4),]
miMERGE_myipubCombat3re_avg_PR_CGshelves_Hypo <- miMERGE_myipubCombat3re_avg_PR_CGshelves[which(miMERGE_myipubCombat3re_avg_PR_CGshelves$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_CGshelves_Hyper <- miMERGE_myipubCombat3re_avg_PR_CGshelves[which(miMERGE_myipubCombat3re_avg_PR_CGshelves$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_CGshelves_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_CGshelves_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_CGshelves_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_CGshelves_Hypo)
mean_myipubCombat3re_avg_PR_CGshelves_Hypo["CGshelves"] <- "CGshelves"
colnames(mean_myipubCombat3re_avg_PR_CGshelves_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_CGshelves_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_CGshelves_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_CGshelves_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_CGshelves_Hyper)
mean_myipubCombat3re_avg_PR_CGshelves_Hyper["CGshelves"] <- "CGshelves"
colnames(mean_myipubCombat3re_avg_PR_CGshelves_Hyper) <- c("Values","Sample","Feature")


miMERGE_myipubCombat3re_avg_PR_CGinter <- read.table("miMERGE_myipubCombat3re_avg_PR_CGinter.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_CGinter) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_CGinter)
miMERGE_myipubCombat3re_avg_PR_CGinter <- miMERGE_myipubCombat3re_avg_PR_CGinter[!duplicated(miMERGE_myipubCombat3re_avg_PR_CGinter$row_4),]
miMERGE_myipubCombat3re_avg_PR_CGinter_Hypo <- miMERGE_myipubCombat3re_avg_PR_CGinter[which(miMERGE_myipubCombat3re_avg_PR_CGinter$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_CGinter_Hyper <- miMERGE_myipubCombat3re_avg_PR_CGinter[which(miMERGE_myipubCombat3re_avg_PR_CGinter$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_CGinter_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_CGinter_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_CGinter_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_CGinter_Hypo)
mean_myipubCombat3re_avg_PR_CGinter_Hypo["CGinter"] <- "CGinter"
colnames(mean_myipubCombat3re_avg_PR_CGinter_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_CGinter_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_CGinter_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_CGinter_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_CGinter_Hyper)
mean_myipubCombat3re_avg_PR_CGinter_Hyper["CGinter"] <- "CGinter"
colnames(mean_myipubCombat3re_avg_PR_CGinter_Hyper) <- c("Values","Sample","Feature")


miMERGE_myipubCombat3re_avg_PR_cds <- read.table("miMERGE_myipubCombat3re_avg_PR_cds.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_cds) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_cds)
miMERGE_myipubCombat3re_avg_PR_cds <- miMERGE_myipubCombat3re_avg_PR_cds[!duplicated(miMERGE_myipubCombat3re_avg_PR_cds$row_4),]
miMERGE_myipubCombat3re_avg_PR_cds_Hypo <- miMERGE_myipubCombat3re_avg_PR_cds[which(miMERGE_myipubCombat3re_avg_PR_cds$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_cds_Hyper <- miMERGE_myipubCombat3re_avg_PR_cds[which(miMERGE_myipubCombat3re_avg_PR_cds$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_cds_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_cds_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_cds_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_cds_Hypo)
mean_myipubCombat3re_avg_PR_cds_Hypo["cds"] <- "cds"
colnames(mean_myipubCombat3re_avg_PR_cds_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_cds_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_cds_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_cds_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_cds_Hyper)
mean_myipubCombat3re_avg_PR_cds_Hyper["cds"] <- "cds"
colnames(mean_myipubCombat3re_avg_PR_cds_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_5UTR <- read.table("miMERGE_myipubCombat3re_avg_PR_5UTR.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_5UTR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_5UTR)
miMERGE_myipubCombat3re_avg_PR_5UTR <- miMERGE_myipubCombat3re_avg_PR_5UTR[!duplicated(miMERGE_myipubCombat3re_avg_PR_5UTR$row_4),]
miMERGE_myipubCombat3re_avg_PR_5UTR_Hypo <- miMERGE_myipubCombat3re_avg_PR_5UTR[which(miMERGE_myipubCombat3re_avg_PR_5UTR$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_5UTR_Hyper <- miMERGE_myipubCombat3re_avg_PR_5UTR[which(miMERGE_myipubCombat3re_avg_PR_5UTR$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_5UTR_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_5UTR_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_5UTR_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_5UTR_Hypo)
mean_myipubCombat3re_avg_PR_5UTR_Hypo["5UTR"] <- "5UTR"
colnames(mean_myipubCombat3re_avg_PR_5UTR_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_5UTR_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_5UTR_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_5UTR_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_5UTR_Hyper)
mean_myipubCombat3re_avg_PR_5UTR_Hyper["5UTR"] <- "5UTR"
colnames(mean_myipubCombat3re_avg_PR_5UTR_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_3UTR <- read.table("miMERGE_myipubCombat3re_avg_PR_3UTR.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_3UTR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_3UTR)
miMERGE_myipubCombat3re_avg_PR_3UTR <- miMERGE_myipubCombat3re_avg_PR_3UTR[!duplicated(miMERGE_myipubCombat3re_avg_PR_3UTR$row_4),]
miMERGE_myipubCombat3re_avg_PR_3UTR_Hypo <- miMERGE_myipubCombat3re_avg_PR_3UTR[which(miMERGE_myipubCombat3re_avg_PR_3UTR$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_3UTR_Hyper <- miMERGE_myipubCombat3re_avg_PR_3UTR[which(miMERGE_myipubCombat3re_avg_PR_3UTR$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_3UTR_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_3UTR_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_3UTR_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_3UTR_Hypo)
mean_myipubCombat3re_avg_PR_3UTR_Hypo["3UTR"] <- "3UTR"
colnames(mean_myipubCombat3re_avg_PR_3UTR_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_3UTR_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_3UTR_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_3UTR_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_3UTR_Hyper)
mean_myipubCombat3re_avg_PR_3UTR_Hyper["3UTR"] <- "3UTR"
colnames(mean_myipubCombat3re_avg_PR_3UTR_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_intergenic <- read.table("miMERGE_myipubCombat3re_avg_PR_intergenic.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_intergenic) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_intergenic)
miMERGE_myipubCombat3re_avg_PR_intergenic <- miMERGE_myipubCombat3re_avg_PR_intergenic[!duplicated(miMERGE_myipubCombat3re_avg_PR_intergenic$row_4),]
miMERGE_myipubCombat3re_avg_PR_intergenic_Hypo <- miMERGE_myipubCombat3re_avg_PR_intergenic[which(miMERGE_myipubCombat3re_avg_PR_intergenic$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_intergenic_Hyper <- miMERGE_myipubCombat3re_avg_PR_intergenic[which(miMERGE_myipubCombat3re_avg_PR_intergenic$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_intergenic_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_intergenic_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_intergenic_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_intergenic_Hypo)
mean_myipubCombat3re_avg_PR_intergenic_Hypo["intergenic"] <- "intergenic"
colnames(mean_myipubCombat3re_avg_PR_intergenic_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_intergenic_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_intergenic_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_intergenic_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_intergenic_Hyper)
mean_myipubCombat3re_avg_PR_intergenic_Hyper["intergenic"] <- "intergenic"
colnames(mean_myipubCombat3re_avg_PR_intergenic_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_exon <- read.table("miMERGE_myipubCombat3re_avg_PR_exon.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_exon) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_exon)
miMERGE_myipubCombat3re_avg_PR_exon <- miMERGE_myipubCombat3re_avg_PR_exon[!duplicated(miMERGE_myipubCombat3re_avg_PR_exon$row_4),]
miMERGE_myipubCombat3re_avg_PR_exon_Hypo <- miMERGE_myipubCombat3re_avg_PR_exon[which(miMERGE_myipubCombat3re_avg_PR_exon$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_exon_Hyper <- miMERGE_myipubCombat3re_avg_PR_exon[which(miMERGE_myipubCombat3re_avg_PR_exon$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_exon_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_exon_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_exon_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_exon_Hypo)
mean_myipubCombat3re_avg_PR_exon_Hypo["exon"] <- "exon"
colnames(mean_myipubCombat3re_avg_PR_exon_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_exon_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_exon_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_exon_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_exon_Hyper)
mean_myipubCombat3re_avg_PR_exon_Hyper["exon"] <- "exon"
colnames(mean_myipubCombat3re_avg_PR_exon_Hyper) <- c("Values","Sample","Feature")


miMERGE_myipubCombat3re_avg_PR_intron <- read.table("miMERGE_myipubCombat3re_avg_PR_intron.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_intron) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(annotations_b_hg38_eldf))
dim(miMERGE_myipubCombat3re_avg_PR_intron)
miMERGE_myipubCombat3re_avg_PR_intron <- miMERGE_myipubCombat3re_avg_PR_intron[!duplicated(miMERGE_myipubCombat3re_avg_PR_intron$row_4),]
miMERGE_myipubCombat3re_avg_PR_intron_Hypo <- miMERGE_myipubCombat3re_avg_PR_intron[which(miMERGE_myipubCombat3re_avg_PR_intron$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_intron_Hyper <- miMERGE_myipubCombat3re_avg_PR_intron[which(miMERGE_myipubCombat3re_avg_PR_intron$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_intron_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_intron_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_intron_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_intron_Hypo)
mean_myipubCombat3re_avg_PR_intron_Hypo["intron"] <- "intron"
colnames(mean_myipubCombat3re_avg_PR_intron_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_intron_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_intron_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_intron_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_intron_Hyper)
mean_myipubCombat3re_avg_PR_intron_Hyper["intron"] <- "intron"
colnames(mean_myipubCombat3re_avg_PR_intron_Hyper) <- c("Values","Sample","Feature")


miMERGE_myipubCombat3re_avg_PR_LINE <- read.table("miMERGE_myipubCombat3re_avg_PR_LINE.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_LINE) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(rmsk))
dim(miMERGE_myipubCombat3re_avg_PR_LINE)
miMERGE_myipubCombat3re_avg_PR_LINE <- miMERGE_myipubCombat3re_avg_PR_LINE[!duplicated(miMERGE_myipubCombat3re_avg_PR_LINE$row_4),]
miMERGE_myipubCombat3re_avg_PR_LINE_Hypo <- miMERGE_myipubCombat3re_avg_PR_LINE[which(miMERGE_myipubCombat3re_avg_PR_LINE$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_LINE_Hyper <- miMERGE_myipubCombat3re_avg_PR_LINE[which(miMERGE_myipubCombat3re_avg_PR_LINE$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_LINE_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_LINE_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_LINE_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_LINE_Hypo)
mean_myipubCombat3re_avg_PR_LINE_Hypo["LINE"] <- "LINE"
colnames(mean_myipubCombat3re_avg_PR_LINE_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_LINE_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_LINE_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_LINE_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_LINE_Hyper)
mean_myipubCombat3re_avg_PR_LINE_Hyper["LINE"] <- "LINE"
colnames(mean_myipubCombat3re_avg_PR_LINE_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_LTR <- read.table("miMERGE_myipubCombat3re_avg_PR_LTR.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_LTR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(rmsk))
dim(miMERGE_myipubCombat3re_avg_PR_LTR)
miMERGE_myipubCombat3re_avg_PR_LTR <- miMERGE_myipubCombat3re_avg_PR_LTR[!duplicated(miMERGE_myipubCombat3re_avg_PR_LTR$row_4),]
miMERGE_myipubCombat3re_avg_PR_LTR_Hypo <- miMERGE_myipubCombat3re_avg_PR_LTR[which(miMERGE_myipubCombat3re_avg_PR_LTR$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_LTR_Hyper <- miMERGE_myipubCombat3re_avg_PR_LTR[which(miMERGE_myipubCombat3re_avg_PR_LTR$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_LTR_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_LTR_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_LTR_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_LTR_Hypo)
mean_myipubCombat3re_avg_PR_LTR_Hypo["LTR"] <- "LTR"
colnames(mean_myipubCombat3re_avg_PR_LTR_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_LTR_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_LTR_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_LTR_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_LTR_Hyper)
mean_myipubCombat3re_avg_PR_LTR_Hyper["LTR"] <- "LTR"
colnames(mean_myipubCombat3re_avg_PR_LTR_Hyper) <- c("Values","Sample","Feature")


miMERGE_myipubCombat3re_avg_PR_SINE <- read.table("miMERGE_myipubCombat3re_avg_PR_SINE.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_SINE) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(rmsk))
dim(miMERGE_myipubCombat3re_avg_PR_SINE)
miMERGE_myipubCombat3re_avg_PR_SINE <- miMERGE_myipubCombat3re_avg_PR_SINE[!duplicated(miMERGE_myipubCombat3re_avg_PR_SINE$row_4),]
miMERGE_myipubCombat3re_avg_PR_SINE_Hypo <- miMERGE_myipubCombat3re_avg_PR_SINE[which(miMERGE_myipubCombat3re_avg_PR_SINE$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_SINE_Hyper <- miMERGE_myipubCombat3re_avg_PR_SINE[which(miMERGE_myipubCombat3re_avg_PR_SINE$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_SINE_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_SINE_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_SINE_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_SINE_Hypo)
mean_myipubCombat3re_avg_PR_SINE_Hypo["SINE"] <- "SINE"
colnames(mean_myipubCombat3re_avg_PR_SINE_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_SINE_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_SINE_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_SINE_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_SINE_Hyper)
mean_myipubCombat3re_avg_PR_SINE_Hyper["SINE"] <- "SINE"
colnames(mean_myipubCombat3re_avg_PR_SINE_Hyper) <- c("Values","Sample","Feature")

miMERGE_myipubCombat3re_avg_PR_Satellite <- read.table("miMERGE_myipubCombat3re_avg_PR_Satellite.txt", header = F)
colnames(miMERGE_myipubCombat3re_avg_PR_Satellite) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR), colnames(rmsk))
dim(miMERGE_myipubCombat3re_avg_PR_Satellite)
miMERGE_myipubCombat3re_avg_PR_Satellite <- miMERGE_myipubCombat3re_avg_PR_Satellite[!duplicated(miMERGE_myipubCombat3re_avg_PR_Satellite$row_4),]
miMERGE_myipubCombat3re_avg_PR_Satellite_Hypo <- miMERGE_myipubCombat3re_avg_PR_Satellite[which(miMERGE_myipubCombat3re_avg_PR_Satellite$PR_Con < -0.2),]
miMERGE_myipubCombat3re_avg_PR_Satellite_Hyper <- miMERGE_myipubCombat3re_avg_PR_Satellite[which(miMERGE_myipubCombat3re_avg_PR_Satellite$PR_Con > 0.2),]
mean_myipubCombat3re_avg_PR_Satellite_Hypo <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_Satellite_Hypo[,c(5:22)]))
mean_myipubCombat3re_avg_PR_Satellite_Hypo["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_Satellite_Hypo)
mean_myipubCombat3re_avg_PR_Satellite_Hypo["Satellite"] <- "Satellite"
colnames(mean_myipubCombat3re_avg_PR_Satellite_Hypo) <- c("Values","Sample","Feature")
mean_myipubCombat3re_avg_PR_Satellite_Hyper <- data.frame(colMeans(miMERGE_myipubCombat3re_avg_PR_Satellite_Hyper[,c(5:22)]))
mean_myipubCombat3re_avg_PR_Satellite_Hyper["Sample"] <- rownames(mean_myipubCombat3re_avg_PR_Satellite_Hyper)
mean_myipubCombat3re_avg_PR_Satellite_Hyper["Satellite"] <- "Satellite"
colnames(mean_myipubCombat3re_avg_PR_Satellite_Hyper) <- c("Values","Sample","Feature")

mean_myipubCombat3re_avg_PR_gen_el_Hypo <- rbind.data.frame(mean_myipubCombat3re_avg_PR_promoter_Hypo,
                                                            mean_myipubCombat3re_avg_PR_CGI_Hypo,
                                                            mean_myipubCombat3re_avg_PR_CGshores_Hypo,
                                                            mean_myipubCombat3re_avg_PR_CGshelves_Hypo,
                                                            mean_myipubCombat3re_avg_PR_CGinter_Hypo,
                                                            mean_myipubCombat3re_avg_PR_cds_Hypo,
                                                            mean_myipubCombat3re_avg_PR_5UTR_Hypo,
                                                            mean_myipubCombat3re_avg_PR_3UTR_Hypo,
                                                            mean_myipubCombat3re_avg_PR_intergenic_Hypo,
                                                            mean_myipubCombat3re_avg_PR_exon_Hypo,
                                                            mean_myipubCombat3re_avg_PR_intron_Hypo,
                                                            mean_myipubCombat3re_avg_PR_LINE_Hypo,
                                                            mean_myipubCombat3re_avg_PR_LTR_Hypo,
                                                            mean_myipubCombat3re_avg_PR_SINE_Hypo,
                                                            mean_myipubCombat3re_avg_PR_Satellite_Hypo)

library(tidyr)
tabmean_myipubCombat3re_avg_PR_gen_el_Hypo <- spread(mean_myipubCombat3re_avg_PR_gen_el_Hypo, Sample, Values)
rownames(tabmean_myipubCombat3re_avg_PR_gen_el_Hypo) <- tabmean_myipubCombat3re_avg_PR_gen_el_Hypo$Feature
tabmean_myipubCombat3re_avg_PR_gen_el_Hypo <- tabmean_myipubCombat3re_avg_PR_gen_el_Hypo[,c(1,2,17,18,3,5,6,4)]
rownames(tabmean_myipubCombat3re_avg_PR_gen_el_Hypo) <- paste(c("b","c","g","h","i","j","k","d","e","f","l","m","a","o","n"),"_",tabmean_myipubCombat3re_avg_PR_gen_el_Hypo$Feature, sep="")
tabmean_myipubCombat3re_avg_PR_gen_el_Hypo <- tabmean_myipubCombat3re_avg_PR_gen_el_Hypo[,-1]
tabmean_myipubCombat3re_avg_PR_gen_el_Hypo <- as.matrix(tabmean_myipubCombat3re_avg_PR_gen_el_Hypo)
tabmean_myipubCombat3re_avg_PR_gen_el_Hypo <- tabmean_myipubCombat3re_avg_PR_gen_el_Hypo[order(rownames(tabmean_myipubCombat3re_avg_PR_gen_el_Hypo)),]

pheatmap(tabmean_myipubCombat3re_avg_PR_gen_el_Hypo,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 3)

sticom_gen_el_PR_Hypo <- data.frame(stack(as.matrix(tabmean_myipubCombat3re_avg_PR_gen_el_Hypo)))

sticom_gen_el_PR_Hypo %>%
  ggplot() +
  geom_point(aes(x=row, y=value,colour = col),size=3,alpha=1,stat = "identity", position = position_dodge2(0.8))+theme_classic()+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))+
  scale_color_manual(values=c("darkgrey","#EC7063","#922B21","#48C9B0","#52BE80","#148F77","#1E8449"))+ylim(0,1)+
  scale_fill_manual(values=c("white","white","white","white","white","white","white"))
ggsave("point_sticom_gen_el_PR_Hypo.svg", width=15, height=2, units="in", dpi=96)
ggsave("point_sticom_gen_el_PR_Hypo.png", width=15, height=2, units="in", dpi=96)


mean_myipubCombat3re_avg_PR_gen_el_Hyper <- rbind.data.frame(mean_myipubCombat3re_avg_PR_promoter_Hyper,
                                                             mean_myipubCombat3re_avg_PR_CGI_Hyper,
                                                             mean_myipubCombat3re_avg_PR_CGshores_Hyper,
                                                             mean_myipubCombat3re_avg_PR_CGshelves_Hyper,
                                                             mean_myipubCombat3re_avg_PR_CGinter_Hyper,
                                                             mean_myipubCombat3re_avg_PR_cds_Hyper,
                                                             mean_myipubCombat3re_avg_PR_5UTR_Hyper,
                                                             mean_myipubCombat3re_avg_PR_3UTR_Hyper,
                                                             mean_myipubCombat3re_avg_PR_intergenic_Hyper,
                                                             mean_myipubCombat3re_avg_PR_exon_Hyper,
                                                             mean_myipubCombat3re_avg_PR_intron_Hyper,
                                                             mean_myipubCombat3re_avg_PR_LINE_Hyper,
                                                             mean_myipubCombat3re_avg_PR_LTR_Hyper,
                                                             mean_myipubCombat3re_avg_PR_SINE_Hyper,
                                                             mean_myipubCombat3re_avg_PR_Satellite_Hyper)

tabmean_myipubCombat3re_avg_PR_gen_el_Hyper <- spread(mean_myipubCombat3re_avg_PR_gen_el_Hyper, Sample, Values)
rownames(tabmean_myipubCombat3re_avg_PR_gen_el_Hyper) <- tabmean_myipubCombat3re_avg_PR_gen_el_Hyper$Feature
tabmean_myipubCombat3re_avg_PR_gen_el_Hyper <- tabmean_myipubCombat3re_avg_PR_gen_el_Hyper[,c(1,2,17,18,3,5,6,4)]
rownames(tabmean_myipubCombat3re_avg_PR_gen_el_Hyper) <- paste(c("b","c","g","h","i","j","k","d","e","f","l","m","a","o","n"),"_",tabmean_myipubCombat3re_avg_PR_gen_el_Hyper$Feature, sep="")
tabmean_myipubCombat3re_avg_PR_gen_el_Hyper <- tabmean_myipubCombat3re_avg_PR_gen_el_Hyper[,-1]
tabmean_myipubCombat3re_avg_PR_gen_el_Hyper <- as.matrix(tabmean_myipubCombat3re_avg_PR_gen_el_Hyper)
tabmean_myipubCombat3re_avg_PR_gen_el_Hyper <- tabmean_myipubCombat3re_avg_PR_gen_el_Hyper[order(rownames(tabmean_myipubCombat3re_avg_PR_gen_el_Hyper)),]

pheatmap(tabmean_myipubCombat3re_avg_PR_gen_el_Hyper,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "grey60",
         cutree_cols = 3)

sticom_gen_el_PR_Hyper <- data.frame(stack(as.matrix(tabmean_myipubCombat3re_avg_PR_gen_el_Hyper)))

sticom_gen_el_PR_Hyper %>%
  ggplot() +
  geom_point(aes(x=row, y=value,colour = col),size=3,alpha=1,stat = "identity", position = position_dodge2(0.8))+theme_classic()+
  geom_line(aes(x=row, y=value),stat = "identity", position = position_dodge2(0.8))+
  scale_color_manual(values=c("darkgrey","#EC7063","#922B21","#48C9B0","#52BE80","#148F77","#1E8449"))+ylim(0,1)+
  scale_fill_manual(values=c("white","white","white","white","white","white","white"))
ggsave("point_sticom_gen_el_PR_Hyper.svg", width=15, height=2, units="in", dpi=96)
ggsave("point_sticom_gen_el_PR_Hyper.png", width=15, height=2, units="in", dpi=96)


#Which genomic elements is the most affected?
#PG,C13,C50
#Hyper
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper)
pMicom_in_PG_Hyper <- miMERGE_myipubCombat3re_chrpos_avg_PG_Hyper
write.table(pMicom_in_PG_Hyper, "pMicom_in_PG_Hyper.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
dim(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo)
pMicom_in_PG_Hypo <- miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo
write.table(pMicom_in_PG_Hypo, "pMicom_in_PG_Hypo.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)

#For recovery in corrected clones probes meth tend to be higher than Hypo ICF cutoff (> -0.2 no upper limit) or 
#less than Hyper ICF cutoff (< 0.2 no lower limit)
dim(pMicom_in_PG_Hyper)#pMicom_in_PG_Hyper.txt")
pMicom_in_PG_Hyper_C13rec <- pMicom_in_PG_Hyper[which(pMicom_in_PG_Hyper$C13_Con < 0.2),] 
pMicom_in_PG_Hyper_C13rec_C50rec <- pMicom_in_PG_Hyper_C13rec[which(pMicom_in_PG_Hyper_C13rec$C50_Con < 0.2),] 
dim(pMicom_in_PG_Hyper_C13rec_C50rec)
write.table(pMicom_in_PG_Hyper_C13rec_C50rec, "pMicom_in_PG_Hyper_C13rec_C50rec.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
dim(pMicom_in_PG_Hypo)#pMicom_in_PG_Hypo.txt")
pMicom_in_PG_Hypo_C13rec <- pMicom_in_PG_Hypo[which(pMicom_in_PG_Hypo$C13_Con > -0.2),] 
pMicom_in_PG_Hypo_C13rec_C50rec <- pMicom_in_PG_Hypo_C13rec[which(pMicom_in_PG_Hypo_C13rec$C50_Con > -0.2),] 
dim(pMicom_in_PG_Hypo_C13rec_C50rec)
write.table(pMicom_in_PG_Hypo_C13rec_C50rec, "pMicom_in_PG_Hypo_C13rec_C50rec.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)


#PR,C7,C35
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper)
pMicom_in_PR_Hyper <- miMERGE_myipubCombat3re_chrpos_avg_PR_Hyper
write.table(pMicom_in_PR_Hyper, "pMicom_in_PR_Hyper.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
dim(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo)
pMicom_in_PR_Hypo <- miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo
write.table(pMicom_in_PR_Hypo, "pMicom_in_PR_Hypo.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)



dim(pMicom_in_PR_Hyper)#pMicom_in_PR_Hyper.txt")
pMicom_in_PR_Hyper_C7rec <- pMicom_in_PR_Hyper[which(pMicom_in_PR_Hyper$C7_Con < 0.2),] 
pMicom_in_PR_Hyper_C7rec_C35rec <- pMicom_in_PR_Hyper_C7rec[which(pMicom_in_PR_Hyper_C7rec$C35_Con < 0.2),] 
dim(pMicom_in_PR_Hyper_C7rec_C35rec)
write.table(pMicom_in_PR_Hyper_C7rec_C35rec, "pMicom_in_PR_Hyper_C7rec_C35rec.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
dim(pMicom_in_PR_Hypo)#pMicom_in_PR_Hypo.txt")
pMicom_in_PR_Hypo_C7rec <- pMicom_in_PR_Hypo[which(pMicom_in_PR_Hypo$C7_Con > -0.2),] 
pMicom_in_PR_Hypo_C7rec_C35rec <- pMicom_in_PR_Hypo_C7rec[which(pMicom_in_PR_Hypo_C7rec$C35_Con > -0.2),] 
dim(pMicom_in_PR_Hypo_C7rec_C35rec)
write.table(pMicom_in_PR_Hypo_C7rec_C35rec, "pMicom_in_PR_Hypo_C7rec_C35rec.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)

awk -F' '  '{print $6"\t"$8}' methylation_recovery_count.sh | awk  -F'_'  '{print $3"_"$4"_"$5"_"$8"_"$10}' > methylation_recovery_colnames.txt
./methylation_recovery_count.sh > methylation_recovery_values.txt

system("paste methylation_recovery_colnames.txt methylation_recovery_values.txt > methylation_recovery_colnames_values.txt")

#manully remove unwanted text and save as methylation_recovery_colnames_values_re.txt")
#Apply a python script to assign shape range
methylation_recovery_colnames_values_re <- read.table("methylation_recovery_colnames_values_re_shaped.txt")


methylation_recovery_colnames_values_rearrange <- cbind.data.frame(methylation_recovery_colnames_values_re[c(1:15, 31:45,61:75,91:105),],
                                                                   methylation_recovery_colnames_values_re[c(16:30,46:60,76:90,106:120),])


colnames(methylation_recovery_colnames_values_rearrange) <- c("ICFs","TotalDMCpGs","ShapeT","ICFsCorr","recoveredCpGs","ShapeR")
methylation_recovery_colnames_values_rearrange["NonrecovGs"] <- methylation_recovery_colnames_values_rearrange$TotalDMCpGs - methylation_recovery_colnames_values_rearrange$recoveredCpGs
head(methylation_recovery_colnames_values_rearrange)
methylation_recovery_colnames_values_rearrange["Prop_recov"] <- methylation_recovery_colnames_values_rearrange$recoveredCpGs / methylation_recovery_colnames_values_rearrange$TotalDMCpGs
methylation_recovery_colnames_values_rearrange["Prop_Nonrecov"] <- methylation_recovery_colnames_values_rearrange$NonrecovGs / methylation_recovery_colnames_values_rearrange$TotalDMCpGs
dim(methylation_recovery_colnames_values_rearrange)
rownames(methylation_recovery_colnames_values_rearrange) <- methylation_recovery_colnames_values_rearrange$ICFs
#PG_Hyper
methylation_recovery_colnames_values_rearrange_PG_Hyper <- methylation_recovery_colnames_values_rearrange[1:15,c(8:9,2:3,5:6)]
methylation_recovery_colnames_values_rearrange_PG_Hyper["FeatureRe"] <- rownames(methylation_recovery_colnames_values_rearrange_PG_Hyper)
tabmean_PG_gen_el_Hyper <- spread(mean_myipubCombat3re_avg_PG_gen_el_Hyper, Sample, Values)
tabmean_PG_gen_el_Hyper <- tabmean_PG_gen_el_Hyper[,c(1,2,17,18,3,5,6,4)]
tabmean_PG_gen_el_Hyper["C13_C50"] <- (tabmean_PG_gen_el_Hyper$C13 + tabmean_PG_gen_el_Hyper$C50)/2
tabmean_PG_gen_el_Hyper["FeatureRe"] <- paste(c("PG_Hyper"),"_",tabmean_PG_gen_el_Hyper$Feature, sep="")
methylation_recovery_colnames_values_merge_PG_Hyper <- merge(tabmean_PG_gen_el_Hyper, methylation_recovery_colnames_values_rearrange_PG_Hyper, by="FeatureRe")
rownames(methylation_recovery_colnames_values_merge_PG_Hyper) <- methylation_recovery_colnames_values_merge_PG_Hyper$FeatureRe

stmethylation_recovery_colnames_values_merge_PG_Hyper_meth <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PG_Hyper[,c(4,10)])))
stmethylation_recovery_colnames_values_merge_PG_Hyper_rec <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PG_Hyper[,c(12,11)])))
stmethylation_recovery_colnames_values_merge_PG_Hyper_size <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PG_Hyper[,c(14,14)])))
stmethylation_recovery_colnames_values_merge_PG_Hyper_rec_meth <- cbind.data.frame(stmethylation_recovery_colnames_values_merge_PG_Hyper_rec,
                                                                                   stmethylation_recovery_colnames_values_merge_PG_Hyper_meth,
                                                                                   stmethylation_recovery_colnames_values_merge_PG_Hyper_size)

colnames(stmethylation_recovery_colnames_values_merge_PG_Hyper_rec_meth) <- c("Feature","tagPropRec","valPropRec","FeatureRe","tagSample","valmeth","Featuresi","Shapetype","Size")
ggplot(stmethylation_recovery_colnames_values_merge_PG_Hyper_rec_meth, aes(x="", y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth)) +
  geom_bar(width = 1, stat = "identity",size=1.3) +
  coord_polar("y", start=0) + facet_wrap(~ Feature, ncol = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
#  scale_fill_gradientn(colours = c("navy", "white", "#D47400"),limits=c(0,1))+
  scale_color_manual(values = c("darkred","darkgreen"))

#No meth
ggplot(stmethylation_recovery_colnames_values_merge_PG_Hyper_rec_meth, aes(x=Size/2, y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth, width=Size)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start=0) +facet_wrap(~ Feature, ncol = 18)+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
  scale_fill_gradientn(colours = c("white", "white", "white"))+
  scale_color_manual(values = c("darkred","darkgreen"))
ggsave("stmethylation_recovery_colnames_values_merge_PG_Hyper_rec_meth.svg", width=47*1.25, height=7*1.25, units="cm", dpi=96)

#PR_Hyper
methylation_recovery_colnames_values_rearrange_PR_Hyper <- methylation_recovery_colnames_values_rearrange[31:45,c(8:9,2:3,5:6)]
methylation_recovery_colnames_values_rearrange_PR_Hyper["FeatureRe"] <- rownames(methylation_recovery_colnames_values_rearrange_PR_Hyper)
tabmean_PR_gen_el_Hyper <- spread(mean_myipubCombat3re_avg_PR_gen_el_Hyper, Sample, Values)
tabmean_PR_gen_el_Hyper <- tabmean_PR_gen_el_Hyper[,c(1,2,17,18,3,5,6,4)]
tabmean_PR_gen_el_Hyper["C7_C35"] <- (tabmean_PR_gen_el_Hyper$C7 + tabmean_PR_gen_el_Hyper$C35)/2
tabmean_PR_gen_el_Hyper["FeatureRe"] <- paste(c("PR_Hyper"),"_",tabmean_PR_gen_el_Hyper$Feature, sep="")
methylation_recovery_colnames_values_merge_PR_Hyper <- merge(tabmean_PR_gen_el_Hyper, methylation_recovery_colnames_values_rearrange_PR_Hyper, by="FeatureRe")
rownames(methylation_recovery_colnames_values_merge_PR_Hyper) <- methylation_recovery_colnames_values_merge_PR_Hyper$FeatureRe
stmethylation_recovery_colnames_values_merge_PR_Hyper_meth <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PR_Hyper[,c(5,10)])))
stmethylation_recovery_colnames_values_merge_PR_Hyper_rec <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PR_Hyper[,c(12,11)])))
stmethylation_recovery_colnames_values_merge_PR_Hyper_size <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PR_Hyper[,c(14,14)])))

stmethylation_recovery_colnames_values_merge_PR_Hyper_rec_meth <- cbind.data.frame(stmethylation_recovery_colnames_values_merge_PR_Hyper_rec,
                                                                                   stmethylation_recovery_colnames_values_merge_PR_Hyper_meth,
                                                                                   stmethylation_recovery_colnames_values_merge_PR_Hyper_size)

colnames(stmethylation_recovery_colnames_values_merge_PR_Hyper_rec_meth) <- c("Feature","tagPropRec","valPropRec","FeatureRe","tagSample","valmeth","Featuresi","Shapetype","Size")
ggplot(stmethylation_recovery_colnames_values_merge_PR_Hyper_rec_meth, aes(x="", y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth)) +
  geom_bar(width = 1, stat = "identity",size=1.3) +
  coord_polar("y", start=0) + facet_wrap(~ Feature, ncol = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
  scale_fill_gradientn(colours = c("navy", "white", "#D47400"),limits=c(0,1))+
  scale_color_manual(values = c("darkred","darkgreen"))

#No meth
ggplot(stmethylation_recovery_colnames_values_merge_PR_Hyper_rec_meth, aes(x=Size/2, y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth,width = Size)) +
  geom_bar(stat = "identity",size=1.3) +
  coord_polar("y", start=0) + facet_wrap(~ Feature, ncol = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
  scale_fill_gradientn(colours = c("white", "white", "white"),limits=c(0,1))+
  scale_color_manual(values = c("darkred","darkgreen"))

ggsave("stmethylation_recovery_colnames_values_merge_PR_Hyper_rec_meth.svg", width=47*1.25, height=7*1.25, units="cm", dpi=96)

#PG_Hypo
methylation_recovery_colnames_values_rearrange_PG_Hypo <- methylation_recovery_colnames_values_rearrange[16:30,c(8:9,2:3,5:6)]
methylation_recovery_colnames_values_rearrange_PG_Hypo["FeatureRe"] <- rownames(methylation_recovery_colnames_values_rearrange_PG_Hypo)
tabmean_PG_gen_el_Hypo <- spread(mean_myipubCombat3re_avg_PG_gen_el_Hypo, Sample, Values)
tabmean_PG_gen_el_Hypo <- tabmean_PG_gen_el_Hypo[,c(1,2,17,18,3,5,6,4)]
tabmean_PG_gen_el_Hypo["C13_C50"] <- (tabmean_PG_gen_el_Hypo$C13 + tabmean_PG_gen_el_Hypo$C50)/2
tabmean_PG_gen_el_Hypo["FeatureRe"] <- paste(c("PG_Hypo"),"_",tabmean_PG_gen_el_Hypo$Feature, sep="")
methylation_recovery_colnames_values_merge_PG_Hypo <- merge(tabmean_PG_gen_el_Hypo, methylation_recovery_colnames_values_rearrange_PG_Hypo, by="FeatureRe")
rownames(methylation_recovery_colnames_values_merge_PG_Hypo) <- methylation_recovery_colnames_values_merge_PG_Hypo$FeatureRe
stmethylation_recovery_colnames_values_merge_PG_Hypo_meth <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PG_Hypo[,c(4,10)])))
stmethylation_recovery_colnames_values_merge_PG_Hypo_rec <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PG_Hypo[,c(12,11)])))
stmethylation_recovery_colnames_values_merge_PG_Hypo_size <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PG_Hypo[,c(14,14)])))

stmethylation_recovery_colnames_values_merge_PG_Hypo_rec_meth <- cbind.data.frame(stmethylation_recovery_colnames_values_merge_PG_Hypo_rec,
                                                                                  stmethylation_recovery_colnames_values_merge_PG_Hypo_meth,
                                                                                  stmethylation_recovery_colnames_values_merge_PG_Hypo_size)

colnames(stmethylation_recovery_colnames_values_merge_PG_Hypo_rec_meth) <- c("Feature","tagPropRec","valPropRec","FeatureRe","tagSample","valmeth","Featuresi","Shapetype","Size")
ggplot(stmethylation_recovery_colnames_values_merge_PG_Hypo_rec_meth, aes(x="", y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth)) +
  geom_bar(width = 1, stat = "identity",size=1.3) +
  coord_polar("y", start=0) + facet_wrap(~ Feature, ncol = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
  scale_fill_gradientn(colours = c("navy", "white", "#D47400"),limits=c(0,1))+
  scale_color_manual(values = c("darkred","darkgreen"))

#No meth
ggplot(stmethylation_recovery_colnames_values_merge_PG_Hypo_rec_meth, aes(x=Size/2, y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth, width = Size)) +
  geom_bar(stat = "identity",size=1.3) +
  coord_polar("y", start=0) + facet_wrap(~ Feature, ncol = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
  scale_fill_gradientn(colours = c("white", "white", "white"),limits=c(0,1))+
  scale_color_manual(values = c("darkred","darkgreen"))


ggsave("stmethylation_recovery_colnames_values_merge_PG_Hypo_rec_meth.svg", width=47*1.25, height=7*1.25, units="cm", dpi=96)

#PR_Hypo
methylation_recovery_colnames_values_rearrange_PR_Hypo <- methylation_recovery_colnames_values_rearrange[46:60,c(8:9,2:3,5:6)]
methylation_recovery_colnames_values_rearrange_PR_Hypo["FeatureRe"] <- rownames(methylation_recovery_colnames_values_rearrange_PR_Hypo)
tabmean_PR_gen_el_Hypo <- spread(mean_myipubCombat3re_avg_PR_gen_el_Hypo, Sample, Values)
tabmean_PR_gen_el_Hypo <- tabmean_PR_gen_el_Hypo[,c(1,2,17,18,3,5,6,4)]
tabmean_PR_gen_el_Hypo["C7_C35"] <- (tabmean_PR_gen_el_Hypo$C7 + tabmean_PR_gen_el_Hypo$C35)/2
tabmean_PR_gen_el_Hypo["FeatureRe"] <- paste(c("PR_Hypo"),"_",tabmean_PR_gen_el_Hypo$Feature, sep="")
methylation_recovery_colnames_values_merge_PR_Hypo <- merge(tabmean_PR_gen_el_Hypo, methylation_recovery_colnames_values_rearrange_PR_Hypo, by="FeatureRe")
rownames(methylation_recovery_colnames_values_merge_PR_Hypo) <- methylation_recovery_colnames_values_merge_PR_Hypo$FeatureRe
stmethylation_recovery_colnames_values_merge_PR_Hypo_meth <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PR_Hypo[,c(5,10)])))
stmethylation_recovery_colnames_values_merge_PR_Hypo_rec <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PR_Hypo[,c(12,11)])))
stmethylation_recovery_colnames_values_merge_PR_Hypo_size <- data.frame(stack(as.matrix(methylation_recovery_colnames_values_merge_PR_Hypo[,c(14,14)])))

stmethylation_recovery_colnames_values_merge_PR_Hypo_rec_meth <- cbind.data.frame(stmethylation_recovery_colnames_values_merge_PR_Hypo_rec,
                                                                                  stmethylation_recovery_colnames_values_merge_PR_Hypo_meth,
                                                                                  stmethylation_recovery_colnames_values_merge_PR_Hypo_size)

colnames(stmethylation_recovery_colnames_values_merge_PR_Hypo_rec_meth) <- c("Feature","tagPropRec","valPropRec","FeatureRe","tagSample","valmeth","Featuresi","Shapetype","Size")
ggplot(stmethylation_recovery_colnames_values_merge_PR_Hypo_rec_meth, aes(x="", y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth)) +
  geom_bar(width = 1, stat = "identity",size=1.3) +
  coord_polar("y", start=0) + facet_wrap(~ Feature, ncol = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
  scale_fill_gradientn(colours = c("navy", "white", "#D47400"),limits=c(0,1))+
  scale_color_manual(values = c("darkred","darkgreen"))

#No meth
ggplot(stmethylation_recovery_colnames_values_merge_PR_Hypo_rec_meth, aes(x=Size/2, y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth, width = Size)) +
  geom_bar(stat = "identity",size=1.3) +
  coord_polar("y", start=0) + facet_wrap(~ Feature, ncol = 18) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
  scale_fill_gradientn(colours = c("white", "white", "white"),limits=c(0,1))+
  scale_color_manual(values = c("darkred","darkgreen"))


ggsave("stmethylation_recovery_colnames_values_merge_PR_Hypo_rec_meth.svg", width=47*1.25, height=7*1.25, units="cm", dpi=96)

#Combine_recovery proportion data
methylation_recovery_colnames_values_merge_diff_meth_PG <- rbind.data.frame(methylation_recovery_colnames_values_merge_PG_Hyper,
                                                                            methylation_recovery_colnames_values_merge_PG_Hypo)
methylation_recovery_colnames_values_merge_diff_meth_PR <- rbind.data.frame(methylation_recovery_colnames_values_merge_PR_Hyper,
                                                                            methylation_recovery_colnames_values_merge_PR_Hypo)
write.table(methylation_recovery_colnames_values_merge_diff_meth_PG, "methylation_recovery_colnames_values_merge_diff_meth_PG.txt", row.names = F, sep = "\t", quote = F, append = F)
write.table(methylation_recovery_colnames_values_merge_diff_meth_PR, "methylation_recovery_colnames_values_merge_diff_meth_PR.txt", row.names = F, sep = "\t", quote = F, append = F)

stmethylation_recovery_colnames_values_merge_diff_meth_rec_meth <- rbind.data.frame(stmethylation_recovery_colnames_values_merge_PR_Hypo_rec_meth,
                                                                                    stmethylation_recovery_colnames_values_merge_PG_Hypo_rec_meth,
                                                                                    stmethylation_recovery_colnames_values_merge_PR_Hyper_rec_meth,
                                                                                    stmethylation_recovery_colnames_values_merge_PG_Hyper_rec_meth)

write.table(stmethylation_recovery_colnames_values_merge_diff_meth_rec_meth, "stmethylation_recovery_colnames_values_merge_diff_meth_rec_meth.txt", row.names = F, sep = "\t", quote = F, append = F)


ggplot(stmethylation_recovery_colnames_values_merge_diff_meth_rec_meth, aes(x=Size/2, y=valPropRec, group=tagSample, color=tagPropRec, fill=valmeth, width = Size)) +
  geom_bar(stat = "identity",size=1.3) +
  coord_polar("y", start=0) + facet_wrap(~ Feature, nrow = 4) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.background = element_rect(fill = NA))+
  scale_fill_gradientn(colours = c("white", "white", "white"),limits=c(0,1))+
  scale_color_manual(values = c("darkred","darkgreen"))


ggsave("stmethylation_recovery_colnames_values_merge_diff_meth_rec_meth.svg", width=30*1.25, height=20*1.25, units="cm", dpi=96)


#Statistical analysis
dim(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all)
summary(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all[which(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "PG_Con_All"),])
summary(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all[which(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "C13_Con_All"),])
summary(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all[which(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "C50_Con_All"),])

summary(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined$col == "PG_Con_diffImprinted"),])
summary(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined$col == "C13_Con_diffImprinted"),])
summary(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined$col == "C50_Con_diffImprinted"),])

summary(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh$col == "PG_Con_PCDH"),])
summary(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh$col == "C13_Con_PCDH"),])
summary(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh$col == "C50_Con_PCDH"),])



library(dplyr)
Median_miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH <- miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH%>%
  group_by(col)%>% 
  summarise(Median=median(value))

Median_miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH <- data.frame(Median_miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH)


Median_miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH <- miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH%>%
  group_by(col)%>% 
  summarise(Median=median(value))

Median_miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH <- data.frame(Median_miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH)

Median_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH <- rbind.data.frame(Median_miMean_myipubCombat3re_pos_chravg_PG_Hypo_imp_PCDH,
                                                                              Median_miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH)

tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH <- data.frame(t(Median_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH))
colnames(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH) <- as.character(Median_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$col)
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH <- tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH[-1,]
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PG_C13_all"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PG_Con_All)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C13_Con_All))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PG_C50_all"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PG_Con_All)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C50_Con_All))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PG_C13_imprinted"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PG_Con_diffImprinted)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C13_Con_diffImprinted))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PG_C50_imprinted"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PG_Con_diffImprinted)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C50_Con_diffImprinted))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PG_C13_pcdh"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PG_Con_PCDH)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C13_Con_PCDH))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PG_C50_pcdh"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PG_Con_PCDH)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C50_Con_PCDH))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PR_C7_all"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PR_Con_All)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C7_Con_All))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PR_C35_all"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PR_Con_All)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C35_Con_All))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PR_C7_imprinted"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PR_Con_diffImprinted)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C7_Con_diffImprinted))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PR_C35_imprinted"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PR_Con_diffImprinted)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C35_Con_diffImprinted))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PR_C7_pcdh"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PR_Con_PCDH)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C7_Con_PCDH))
tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH["abs_PR_C35_pcdh"] <- abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$PR_Con_PCDH)) - abs(as.numeric(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH$C35_Con_PCDH))

write.table(tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH, "tMedian_miMean_myipubCombat3re_pos_chravg_PG_PR_Hypo_imp_PCDH.txt", sep = "\t", append = F, quote = F )


#One way Anova with repeated measures then Posthoc
library(rstatix)

stmiMean_myipubCombat3re_pos_chravg_PG_Hypo %>%
  group_by(sample) %>%
  identify_outliers(meandelta_meth)
stmiMean_myipubCombat3re_pos_chravg_PG_Hypo %>%
  group_by(sample) %>%
  shapiro_test(meandelta_meth)

ggqqplot(stmiMean_myipubCombat3re_pos_chravg_PG_Hypo, "meandelta_meth", facet.by = "sample", color = "sample")

#stmiMean_myipubCombat3re_pos_chravg_PG_Hypo
stmiMean_myipubCombat3re_pos_chravg_PG_Hypo.aov <- anova_test(data = stmiMean_myipubCombat3re_pos_chravg_PG_Hypo,
                                                              dv = meandelta_meth,
                                                              wid = region,
                                                              within = sample)
get_anova_table(stmiMean_myipubCombat3re_pos_chravg_PG_Hypo.aov)

stmiMean_myipubCombat3re_pos_chravg_PG_Hypo.pwc <- stmiMean_myipubCombat3re_pos_chravg_PG_Hypo %>%
  pairwise_t_test(
    meandelta_meth ~ sample, paired = TRUE,
    p.adjust.method = "bonferroni"
  )


#stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG
stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG.aov <- anova_test(data = stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG,
                                                                             dv = meandelta_meth,
                                                                             wid = region,
                                                                             within = sample)
get_anova_table(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG.aov)

stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG.pwc <- stmiMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG %>%
  pairwise_t_test(
    meandelta_meth ~ sample, paired = TRUE,
    p.adjust.method = "bonferroni"
  )

#stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG

stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG.aov <- anova_test(data = stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG,
                                                                              dv = meandelta_meth,
                                                                              wid = region,
                                                                              within = sample)
get_anova_table(stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG.aov)

stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG.pwc <- stmiMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG %>%
  pairwise_t_test(
    meandelta_meth ~ sample, paired = TRUE,
    p.adjust.method = "bonferroni"
  )


#stmiMean_myipubCombat3re_pos_chravg_PR_Hypo
stmiMean_myipubCombat3re_pos_chravg_PR_Hypo.aov <- anova_test(data = stmiMean_myipubCombat3re_pos_chravg_PR_Hypo,
                                                              dv = meandelta_meth,
                                                              wid = region,
                                                              within = sample)
get_anova_table(stmiMean_myipubCombat3re_pos_chravg_PR_Hypo.aov)

stmiMean_myipubCombat3re_pos_chravg_PR_Hypo.pwc <- stmiMean_myipubCombat3re_pos_chravg_PR_Hypo %>%
  pairwise_t_test(
    meandelta_meth ~ sample, paired = TRUE,
    p.adjust.method = "bonferroni"
  )


#stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR
stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR.aov <- anova_test(data = stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR,
                                                                             dv = meandelta_meth,
                                                                             wid = region,
                                                                             within = sample)
get_anova_table(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR.aov)

stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR.pwc <- stmiMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR %>%
  pairwise_t_test(
    meandelta_meth ~ sample, paired = TRUE,
    p.adjust.method = "bonferroni"
  )

#stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR

stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR.aov <- anova_test(data = stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR,
                                                                              dv = meandelta_meth,
                                                                              wid = region,
                                                                              within = sample)
get_anova_table(stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR.aov)

stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR.pwc <- stmiMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR %>%
  pairwise_t_test(
    meandelta_meth ~ sample, paired = TRUE,
    p.adjust.method = "bonferroni"
  )




dim(miMean_myipubCombat3re_pos_chravg_PG_Hypo)
dim(miMERGE_myipubCombat3re_pos_chravg_PR_Hypo)

dim(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG)
dim(miMean_myipubCombat3re_chrpos_avg_PR_pR_imp_Hypo_CpGs_PR)

dim(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG)
dim(miMean_myipubCombat3re_chrpos_avg_PR_pR_pcdh_Hypo_CpGs_PR)


#RNASEq
H19_log2RPKM <- data.frame(c((2.810938079 +	2.889594382)/2,	6.704997892,	5.292277061,	5.98803185,	5.4150453,	7.056901935))
rownames(H19_log2RPKM) <- c("aControl",	"bpR",	"c7",	"dpG",	"ec13",	"fc50")

colnames(H19_log2RPKM) <- c("Log2RPKM")
tH19_log2RPKM <- data.frame(t(H19_log2RPKM))
t.test(tH19_log2RPKM$aControl, tH19_log2RPKM$bpR)

summary(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all[which(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "PG_Con_All" | miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "C13_Con_All"),])


wilcox.test(value ~ col, data = (miMean_myipubCombat3re_pos_chravg_PG_Hypo_all[which(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "PG_Con_All" | miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "C13_Con_All"),]),
            exact = FALSE)
wilcox.test(value ~ col, data = (miMean_myipubCombat3re_pos_chravg_PG_Hypo_all[which(miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "PG_Con_All" | miMean_myipubCombat3re_pos_chravg_PG_Hypo_all$col == "C50_Con_All"),]),
            exact = FALSE)

wilcox.test(value ~ col, data = (miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined$col == "PG_Con_diffImprinted" | miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined$col == "C13_Con_diffImprinted"),]),
            exact = FALSE)
wilcox.test(value ~ col, data = (miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined$col == "PG_Con_diffImprinted" | miMean_myipubCombat3re_chrpos_avg_PG_pG_imp_Hypo_CpGs_PG_imprined$col == "C50_Con_diffImprinted"),]),
            exact = FALSE)
wilcox.test(value ~ col, data = (miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh$col == "PG_Con_PCDH" | miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh$col == "C13_Con_PCDH"),]),
            exact = FALSE)
wilcox.test(value ~ col, data = (miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh[which(miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh$col == "PG_Con_PCDH" | miMean_myipubCombat3re_chrpos_avg_PG_pG_pcdh_Hypo_CpGs_PG_pcdh$col == "C50_Con_PCDH"),]),
            exact = FALSE)

ggboxplot(stmiMean_myipubCombat3re_pos_chravg_PG_Hypo, x = "sample", y = "meandelta_meth",
          color = "sample", palette = "jco")+
  stat_compare_means(comparisons = list( c("PG_Con", "C13_Con"), c("C50_Con", "C13_Con"), c("PG_Con", "C50_Con") ))+ # Add pairwise comparisons p-value
  stat_compare_means(method = "t.test",label.y = 1)     # Add global p-value

stat_compare_means(method = "anova", label.y = 40)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0.5") 

ggplot(data=miMean_myipubCombat3re_pos_chravg_PR_Hypo_imp_PCDH, aes(x=col, y=value)) +
  geom_boxplot(aes(col=col,fill=col),position=position_dodge())+
  ylim(-1,1)+theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(fill="white", color = "black"))+
  scale_fill_manual(values = rep(c("white","white","white"), 3))+
  scale_color_manual(values = rep(c("#922B21","#148F77","#1E8449"), 3))+
  facet_wrap(~Color, scale="free") +geom_hline(yintercept=0, linetype="dashed", color = "blue")
 
#DNMT3B peaks
ICR_overlap_DNMT3B_peaks <- read.table("ICR_overlap_DNMT3B_peaks.txt", header = T, stringsAsFactors = F)
ICR_overlap_DNMT3B_peaks_UN <- ICR_overlap_DNMT3B_peaks[which(ICR_overlap_DNMT3B_peaks$IP_replicate == "UN_IP1" | 
                                                                 ICR_overlap_DNMT3B_peaks$IP_replicate == "UN_IP2"),]
ICR_overlap_DNMT3B_peaks_PG <- ICR_overlap_DNMT3B_peaks[which(ICR_overlap_DNMT3B_peaks$IP_replicate == "pG_IP1" | 
                                                                ICR_overlap_DNMT3B_peaks$IP_replicate == "pG_IP2"),]
ICR_overlap_DNMT3B_peaks_PR <- ICR_overlap_DNMT3B_peaks[which(ICR_overlap_DNMT3B_peaks$IP_replicate == "pR_IP1" | 
                                                                ICR_overlap_DNMT3B_peaks$IP_replicate == "pR_IP2"),]
ICR_overlap_DNMT3B_peaks_C7 <- ICR_overlap_DNMT3B_peaks[which(ICR_overlap_DNMT3B_peaks$IP_replicate == "cR7_IP1" | 
                                                                ICR_overlap_DNMT3B_peaks$IP_replicate == "cR7_IP2"),]
ICR_overlap_DNMT3B_peaks_C13 <- ICR_overlap_DNMT3B_peaks[which(ICR_overlap_DNMT3B_peaks$IP_replicate == "cG13_IP1" | 
                                                                ICR_overlap_DNMT3B_peaks$IP_replicate == "cG13_IP2"),]
ICR_overlap_DNMT3B_peaks_C35 <- ICR_overlap_DNMT3B_peaks[which(ICR_overlap_DNMT3B_peaks$IP_replicate == "cR35_IP1" | 
                                                                ICR_overlap_DNMT3B_peaks$IP_replicate == "cR35_IP2"),]
ICR_overlap_DNMT3B_peaks_C50 <- ICR_overlap_DNMT3B_peaks[which(ICR_overlap_DNMT3B_peaks$IP_replicate == "cG50_IP1" | 
                                                                 ICR_overlap_DNMT3B_peaks$IP_replicate == "cG50_IP2"),]


ICR_overlap_DNMT3B_iPSC_peaks <- ICR_overlap_DNMT3B_peaks[,c(7,9)]
colnames(ICR_overlap_DNMT3B_iPSC_peaks) <- c("ICR", "Condition")
#ZFP57 hESC
bedtools intersect -wa -wb -a GSE115387_hES_57_vs_hES_TI_peaks.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed | awk '{print $9"\t""ZFP57_hESC"}' | sort -k1,1 -u > ZFP57_hESC_hg19_GSE115387.bed

ICR_overlap_ZFP57_hESC_peaks <- read.table("ZFP57_hESC_hg19_GSE115387.bed")
colnames(ICR_overlap_ZFP57_hESC_peaks) <- c("ICR", "Condition")

#ZFP57_HEK293 exo ChIP
system("cat GSM2466450_ZFP57_peaks_processed_score_signal_exo.bed GSM2466451_ZFP57-rep2_peaks_processed_score_signal_exo.bed > ZFP57_HEK293T_GRCh37.bed")

bedtools intersect -wa -wb -a ZFP57_HEK293T_GRCh37.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed | awk '{print $10"\t""ZFP57_HEK293T"}' | sort -k1,1 -u > ZFP57_HEK293T_GRCh37_GSM2466450_1.bed

ICR_overlap_ZFP57_HEK293T_peaks <- read.table("ZFP57_HEK293T_GRCh37_GSM2466450_1.bed")
colnames(ICR_overlap_ZFP57_HEK293T_peaks) <- c("ICR", "Condition")

combinedICR_overlap_DNMT3B_ZFP57_peaks <- rbind.data.frame(ICR_overlap_DNMT3B_iPSC_peaks, ICR_overlap_ZFP57_hESC_peaks, ICR_overlap_ZFP57_HEK293T_peaks)
table(combinedICR_overlap_DNMT3B_ZFP57_peaks$ICR, combinedICR_overlap_DNMT3B_ZFP57_peaks$Condition)


combinedICR_overlap_ZFP57_peaks <- rbind.data.frame(ICR_overlap_ZFP57_hESC_peaks, ICR_overlap_ZFP57_HEK293T_peaks)
tab_combinedICR_overlap_ZFP57_peaks <- table(combinedICR_overlap_ZFP57_peaks$ICR, combinedICR_overlap_ZFP57_peaks$Condition)
tab_combinedICR_overlap_ZFP57_peaks <- data.frame(unclass(tab_combinedICR_overlap_ZFP57_peaks))
tab_combinedICR_overlap_ZFP57_peaks$DMR <- rownames(tab_combinedICR_overlap_ZFP57_peaks)


#Heatmap Dnmt3b all ICR
Monk_ICR_dnmt3bx_ICRname <- merge(Monk_human_ICR_hg38, Monk_ICR_dnmt3bx, by.x = "DMRtype", by.y = "row")
Monk_ICR_dnmt3bx_ICRname_zfp57 <- merge(tab_combinedICR_overlap_ZFP57_peaks, Monk_ICR_dnmt3bx_ICRname, by.x = "DMR", by.y = "DMR", all.y=T)
Monk_ICR_dnmt3bx_ICRname_zfp57$ZFP57_HEK293T[is.na(Monk_ICR_dnmt3bx_ICRname_zfp57$ZFP57_HEK293T)] <- 0 
Monk_ICR_dnmt3bx_ICRname_zfp57$ZFP57_hESC[is.na(Monk_ICR_dnmt3bx_ICRname_zfp57$ZFP57_hESC)] <- 0 
rownames(Monk_ICR_dnmt3bx_ICRname_zfp57) <- Monk_ICR_dnmt3bx_ICRname_zfp57$DMR
Monk_ICR_dnmt3bx_ICRname_zfp57 <- Monk_ICR_dnmt3bx_ICRname_zfp57[,c(2,3,8,9,12,10:11,13:14)]
breaksListd2 = seq(0, 2, by = 0.001)
pheatmap(as.matrix(Monk_ICR_dnmt3bx_ICRname_zfp57[,3:9]),
         color = colorRampPalette(c("grey", "white", "#0000ce"))(length(breaksListd2)),
         breaks = breaksListd2,
         fontsize = 8,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white",
         cellwidth = 8,
         treeheight_row = 0,
         treeheight_col = 0)
#save as Monk_ICR_dnmt3bx_ICRname.svg
pheatmap(as.matrix(Monk_ICR_dnmt3bx_ICRname_zfp57[,1:2]),
         color = colorRampPalette(c("white", "black", "black"))(length(breaksListd2)),
         breaks = breaksListd2,
         fontsize = 8,
         cluster_rows = F,
         cluster_cols = F,
         border_color = "white",
         cellwidth = 8,
         treeheight_row = 0,
         treeheight_col = 0)
#save as Monk_ICR_ICRname_zfp57.svg


#Methylation
diff_imp_loci_pG_meth <- data.frame(diff_imp_loci_pG)
diff_imp_loci_pG_meth$Condition <- "pG_affectedICR"
colnames(diff_imp_loci_pG_meth) <- c("ICR", "Condition")

diff_imp_loci_pR_meth <- data.frame(diff_imp_loci_pR)
diff_imp_loci_pR_meth$Condition <- "pR_affectedICR"
colnames(diff_imp_loci_pR_meth) <- c("ICR", "Condition")

combinedICR_overlap_DNMT3B_ZFP57_peaks_meth <- rbind.data.frame(diff_imp_loci_pG_meth, diff_imp_loci_pR_meth, ICR_overlap_ZFP57_hESC_peaks, ICR_overlap_ZFP57_HEK293T_peaks, ICR_overlap_DNMT3B_iPSC_peaks)


tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth <- table(combinedICR_overlap_DNMT3B_ZFP57_peaks_meth$ICR, combinedICR_overlap_DNMT3B_ZFP57_peaks_meth$Condition)
dim(tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth)
head(tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth)

tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth <- tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth[,c(9,12,17,18,15:16,10:11,13,1:4,14,5:8)]
gplots::heatmap.2(as.matrix(tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth1),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none' , sepcolor="black", sepwidth=c(6,6), col = colorRampPalette(c('white', 'black'))(12), key = FALSE)

write.table(tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth ,"tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth.txt", sep = "\t", quote = F, append = F)

tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth1 <- read.table("tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth_re.txt", header= TRUE, stringsAsFactors = F, row.names = 1)
tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth1 <- as.matrix(tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth1)

pheatmap(tabcombinedICR_overlap_DNMT3B_ZFP57_peaks_meth1,
         fontsize = 8,
         cluster_col=F,
         cluster_row=T,
         border_color = "grey60")

system("bedtools intersect -wa -wb -a GSE57989_wt_kap1_peaks.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed | awk '{print $10}' | sort -k1,1 -u")


head(ipubcomdataBin500bp_filt1_counted1)
dim(ipubcomdataBin500bp_filt1_counted1)

ipubcomdataBin500bp_filt1_counted1chr <- ipubcomdataBin500bp_filt1_counted1
ipubcomdataBin500bp_filt1_counted1chr <- cSplit(ipubcomdataBin500bp_filt1_counted1chr, "Bin500bp", "%")
head(ipubcomdataBin500bp_filt1_counted1chr)
dim(ipubcomdataBin500bp_filt1_counted1chr)
ipubcomdataBin500bp_filt1_counted1chr <- ipubcomdataBin500bp_filt1_counted1chr[,c(20:22,1:19)]
write.table(ipubcomdataBin500bp_filt1_counted1chr, "ipubcomdataBin500bp_filt1_counted1chr.txt", sep ="\t", quote=F, append = F, row.names = F, col.names = F)
system("sort -k1,1 -k2,2n ipubcomdataBin500bp_filt1_counted1chr.txt > ipubcomdataBin500bp_filt1_counted1chr.sort.txt")
system("bedtools closest -wa -wb -a ipubcomdataBin500bp_filt1_counted1chr.sort.txt -b /home/ankitv/ref_av/hg38/gencode.v29.annotation.gene.sort.gtf -d > ipubcomdataBin500bp_filt1_counted1chr_genehg38.txt")
#Collapse genenames and Prepare Table
ipubcomdataBin500bp_filt1_counted1chr_genehg38 <- read.table("ipubcomdataBin500bp_filt1_counted1chr_genehg38.txt", header = F, stringsAsFactors = F)
head(ipubcomdataBin500bp_filt1_counted1chr_genehg38)
colnames(ipubcomdataBin500bp_filt1_counted1chr_genehg38) <- c(colnames(ipubcomdataBin500bp_filt1_counted1chr), "GeneChr", "GeneStart", "GeneEnd", "Type","ENSID","Gene","Distance")
head(ipubcomdataBin500bp_filt1_counted1chr_genehg38)
ipubcomdataBin500bp_filt1_counted1chr_genehg38["Bin500bp"] <- paste0(ipubcomdataBin500bp_filt1_counted1chr_genehg38$Bin500bp_1, "%" ,ipubcomdataBin500bp_filt1_counted1chr_genehg38$Bin500bp_2, "%" ,ipubcomdataBin500bp_filt1_counted1chr_genehg38$Bin500bp_3)
ipubcomdataBin500bp_filt1_counted1chr_genehg38dis <- ddply(ipubcomdataBin500bp_filt1_counted1chr_genehg38, .(Bin500bp), summarize, Gene = toString(Gene), Dis = toString(Distance))
dim(ipubcomdataBin500bp_filt1_counted1chr_genehg38dis)
head(ipubcomdataBin500bp_filt1_counted1chr_genehg38dis)
ipubcomdataBin500bp_filt1_counted1chr_named <- merge(ipubcomdataBin500bp_filt1_counted1, 
                                                     ipubcomdataBin500bp_filt1_counted1chr_genehg38dis,
                                                     by = "Bin500bp")

ipubcomdataBin500bp_filt1_counted1chr_named$Genes <- gsub(" ", "", ipubcomdataBin500bp_filt1_counted1chr_named$Gene)
ipubcomdataBin500bp_filt1_counted1chr_named$Distance <- gsub(" ", "", ipubcomdataBin500bp_filt1_counted1chr_named$Dis)
head(ipubcomdataBin500bp_filt1_counted1chr_named)
ipubcomdataBin500bp_filt1_counted1chr_namedre <- cSplit(ipubcomdataBin500bp_filt1_counted1chr_named, "Bin500bp", "%")
head(ipubcomdataBin500bp_filt1_counted1chr_namedre)
dim(ipubcomdataBin500bp_filt1_counted1chr_namedre)
ipubcomdataBin500bp_filt1_counted1chr_namedre <- ipubcomdataBin500bp_filt1_counted1chr_namedre[,c(24:26, 1:19, 22:23)]
ipubcomdataBin500bp_filt1_counted1chr_namedre <- data.frame(ipubcomdataBin500bp_filt1_counted1chr_namedre)
ipubcomdataBin500bp_filt1_counted1chr_namedre["CGdensity"] <- ipubcomdataBin500bp_filt1_counted1chr_namedre$freq / 500
ipubcomdataBin500bp_filt1_counted1chr_namedre <- ipubcomdataBin500bp_filt1_counted1chr_namedre[order(-ipubcomdataBin500bp_filt1_counted1chr_namedre$CGdensity),]
write.table(ipubcomdataBin500bp_filt1_counted1chr_namedre, "ipubcomdataBin500bp_filt1_counted1chr_namedre.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)

dim(ipubcomdataBin500bp_filt1_counted1chr_named)
dim(ipubcomdataBin500bp_Hyper_PG_1_filtmin3chr)
dim(ipubcomdataBin500bp_Hyper_PR_1_filtmin3chr)
dim(ipubcomdataBin500bp_Hypo_PG_1_filtmin3chr)
dim(ipubcomdataBin500bp_Hypo_PR_1_filtmin3chr)

ipubcomdataBin500bp_Hyper_PG_1_filtmin3chr_id <- unique(ipubcomdataBin500bp_Hyper_PG_1_filtmin3chr$Bin500bp_Hyper)
ipubcomdataBin500bp_Hyper_PR_1_filtmin3chr_id <- unique(ipubcomdataBin500bp_Hyper_PR_1_filtmin3chr$Bin500bp_Hyper)
ipubcomdataBin500bp_Hypo_PG_1_filtmin3chr_id <- unique(ipubcomdataBin500bp_Hypo_PG_1_filtmin3chr$Bin500bp_Hypo)
ipubcomdataBin500bp_Hypo_PR_1_filtmin3chr_id <- unique(ipubcomdataBin500bp_Hypo_PR_1_filtmin3chr$Bin500bp_Hypo)

ipubcomdataBin500bp_Hypo_PG_1_filtmin3chr_id <- data.frame(ipubcomdataBin500bp_Hypo_PG_1_filtmin3chr_id)
ipubcomdataBin500bp_Hypo_PR_1_filtmin3chr_id <- data.frame(ipubcomdataBin500bp_Hypo_PR_1_filtmin3chr_id)
colnames(ipubcomdataBin500bp_Hypo_PG_1_filtmin3chr_id) <- c("Bin500bp")
colnames(ipubcomdataBin500bp_Hypo_PR_1_filtmin3chr_id) <- c("Bin500bp")

head(ipubcomdataBin500bp_filt1_counted1chr_named)

ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG <- merge(ipubcomdataBin500bp_filt1_counted1chr_named, ipubcomdataBin500bp_Hypo_PG_1_filtmin3chr_id, by="Bin500bp", all.y =TRUE)
dim(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG)
head(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG)
#Some are NA and need to be removed. Because 
# ipubcomdataBin500bp_filt1_counted1chr_named is all possible bins with >= 3 CpGs
# ipubcomdataBin500bp_Hypo_PG_1_filtmin3chr_id is bins with >=3 CpGs hypomethylated
# If any bin in ipubcomdataBin500bp_filt1_counted1chr_named has less < 2 hypomethylated CpG it will be removed.
ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG <- na.omit(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG)
dim(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG)
head(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG)
ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG["CGdensity"] <- ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG$freq / 500
ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG <- ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG[order(-ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG$CGdensity),]
write.table(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG, "ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)
write.xlsx(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG, "ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PG.xlsx")

ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR <- merge(ipubcomdataBin500bp_filt1_counted1chr_named, ipubcomdataBin500bp_Hypo_PR_1_filtmin3chr_id, by="Bin500bp", all.y =TRUE)
dim(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR)
head(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR)
#Some are NA and need to be removed. Because 
# ipubcomdataBin500bp_filt1_counted1chr_named is all possible bins with >= 3 CpGs
# ipubcomdataBin500bp_Hypo_PR_1_filtmin3chr_id is bins with >=3 CpGs hypomethylated
# If any bin in ipubcomdataBin500bp_filt1_counted1chr_named has less < 2 hypomethylated CpG it will be removed.
ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR <- na.omit(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR)
dim(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR)
head(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR)
ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR["CGdensity"] <- ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR$freq / 500
ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR <- ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR[order(-ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR$CGdensity),]
write.table(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR, "ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR.txt", sep="\t", append = F,quote=F, col.names=T, row.names = F)
write.xlsx(ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR,"ipubcomdataBin500bp_filt1_counted1chr_named_Hypo_PR.xlsx")


system("bedtools intersect -wa -wb -a GSM1003585_hg19_wgEncodeBroadHistoneH1hescH3k09me3StdPk.broadPeak  -b /home/ankitv/ref_av/human_ICR_hg19.bed")
system("bedtools intersect -wa -wb -a GSE99346_WT_DNMT3B-Input_peaks.bed  -b /home/ankitv/ref_av/human_ICR_hg19.bed")
system("bedtools intersect -wa -wb -a GSE115387_hES_57_vs_hES_TI_peaks.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed")
system("bedtools intersect -wa -wb -a GSE57989_wt_kap1_peaks.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed")


system("bedtools intersect -a GSE99346_WT_DNMT3B-Input_peaks.bed -b GSE115387_hES_57_vs_hES_TI_peaks.bed > DNMT3B_ZFP57_hESC.bed")
system("bedtools intersect -a DNMT3B_ZFP57_hESC.bed -b GSE57989_wt_kap1_peaks.bed > DNMT3B_ZFP57_KAP1_hESC.bed")
system("bedtools intersect -a DNMT3B_ZFP57_KAP1_hESC.bed -b GSM1003585_hg19_wgEncodeBroadHistoneH1hescH3k09me3StdPk.broadPeak > DNMT3B_ZFP57_KAP1_H3K9me3_hESC.bed")

system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_ICR_hg38.bed > miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexICR.txt")
miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexICR <- read.table("miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexICR.txt", header = F, stringsAsFactors = F)
colnames(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexICR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo), "chrf","startf","endf")
head(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexICR,50)
pheatmap(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexICR[,c(6:22)])

system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_nonICR_hg38.bed > miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexnonICR.txt")
miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexnonICR <- read.table("miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexnonICR.txt", header = F, stringsAsFactors = F)
colnames(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexnonICR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo), "chrf","startf","endf")
head(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexnonICR,50)
pheatmap(miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr_onsharedcomplexnonICR[,c(6:22)])


system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_ICR_hg38.bed > miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexICR.txt")
miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexICR <- read.table("miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexICR.txt", header = F, stringsAsFactors = F)
colnames(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexICR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo), "chrf","startf","endf")
head(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexICR,50)
pheatmap(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexICR[,c(6:22)])

system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_nonICR_hg38.bed > miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexnonICR.txt")
miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexnonICR <- read.table("miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexnonICR.txt", header = F, stringsAsFactors = F)
colnames(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexnonICR) <- c(colnames(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo), "chrf","startf","endf")
head(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexnonICR,50)
pheatmap(miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr_onsharedcomplexnonICR[,c(6:22)])

#DNMT3B_ZFP57_KAP1_H3K9me3_hESC.bed")
system("bedtools intersect -wa -wb -a DNMT3B_ZFP57_KAP1_H3K9me3_hESC.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed | awk '{print $1"\t"$2"\t"$3}' > DNMT3B_ZFP57_KAP1_H3K9me3_hESC_ICR_hg19.bed")
system("bedtools intersect -wa -wb -a DNMT3B_ZFP57_KAP1_H3K9me3_hESC.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed -v | awk '{print $1"\t"$2"\t"$3}' > DNMT3B_ZFP57_KAP1_H3K9me3_hESC_nonICR_hg19.bed")
#Liftover UCSC
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_ICR_hg38.bed | wc -l") 
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_nonICR_hg38.bed | wc -l") 
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_ICR_hg38.bed | wc -l") 
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_nonICR_hg38.bed | wc -l") 

system("bedtools intersect -wa -wb -a MERGE_myiCombat3re_pos_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_ICR_hg38.bed | wc -l") 
system("bedtools intersect -wa -wb -a MERGE_myiCombat3re_pos_chr.txt -b DNMT3B_ZFP57_KAP1_H3K9me3_hESC_nonICR_hg38.bed | wc -l") 



#DNMT3B_ZFP57_KAP1_hESC.bed")
system("bedtools intersect -wa -wb -a DNMT3B_ZFP57_KAP1_hESC.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed | awk '{print $1"\t"$2"\t"$3}' > DNMT3B_ZFP57_KAP1_hESC_ICR_hg19.bed")
system("bedtools intersect -wa -wb -a DNMT3B_ZFP57_KAP1_hESC.bed -b /home/ankitv/ref_av/human_ICR_hg19.bed -v | awk '{print $1"\t"$2"\t"$3}' > DNMT3B_ZFP57_KAP1_hESC_nonICR_hg19.bed")
#Liftover UCSC
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_hESC_ICR_hg38.bed | wc -l") #23
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_hESC_ICR_hg38.bed | wc -l") #22
system("bedtools intersect -wa -wb -a MERGE_myiCombat3re_pos_chr.txt -b DNMT3B_ZFP57_KAP1_hESC_ICR_hg38.bed | wc -l") #34
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PG_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_hESC_nonICR_hg38.bed | wc -l") #4
system("bedtools intersect -wa -wb -a miMERGE_myipubCombat3re_chrpos_avg_PR_Hypo_chr.txt -b DNMT3B_ZFP57_KAP1_hESC_nonICR_hg38.bed | wc -l") #5
system("bedtools intersect -wa -wb -a MERGE_myiCombat3re_pos_chr.txt -b DNMT3B_ZFP57_KAP1_hESC_nonICR_hg38.bed | wc -l") #44

pie(c(23,34), labels = round(100*c(23,34)/sum(c(23,34)), 1), col=c("#FAF8DF","#DDF2FD"))
pie(c(22,34), labels = round(100*c(22,34)/sum(c(22,34)), 1), col=c("#FAF8DF","#DDF2FD"))

pie(c(4,44), labels = round(100*c(4,44)/sum(c(4,44)), 1), col=c("#FAF8DF","#DDF2FD"))
pie(c(5,44), labels = round(100*c(5,44)/sum(c(5,44)), 1), col=c("#FAF8DF","#DDF2FD"))
print('#---------------------------------   END  OF   ANALSIS   ----------------------------------------#')

