print("Methylation Array Data: ICF Syndrome iPSCs")
print("Array type : EPIC/450K")
print("Human cell line : iPSCs")
print("Data type : Line charts")

library(ChAMP)
library(splitstackshape)
library(ggpubr)
library(ggplot2)
library(plyr)
library(gplots)
library(pheatmap)
setwd("/Users/ankitverma/Documents/Archivio2/Controls")

dim(myCombat3reavg)
head(myCombat3reavg)
############################## Flank Line chart ############################
awk '{print $1"\t"$2-1000"\t"$3+1000"\t"$4"\t"$5}' /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed  > Monk_human_ICR_hg38_flankplusminus1000.bed

bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b Monk_human_ICR_hg38_flankplusminus1000.bed > MERGE_myicomCGbatCG3re_human_ICR_flanked1000.txt
MERGE_myicomCGbatCG3re_human_ICR_flanked1000 <- read.table("MERGE_myicomCGbatCG3re_human_ICR_flanked1000.txt", header = F)
head(MERGE_myicomCGbatCG3re_human_ICR_flanked1000)
dim(MERGE_myicomCGbatCG3re_human_ICR_flanked1000)
head(MERGE_myicomCGbatCG3re_human_ICR_flanked1000)
colnames(MERGE_myicomCGbatCG3re_human_ICR_flanked1000) <- c("chrc","startc","endc","id", colnames(MERGE_myiCombat3reavg_pos_chr[,5:22]),"chr","start","end","DMR", "DMRType")
fi1000comCGdata1 <- MERGE_myicomCGbatCG3re_human_ICR_flanked1000[,c((ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked1000)-1),ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked1000), 1:(ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked1000)-5))]
rownames(fi1000comCGdata1)
colnames(fi1000comCGdata1)
head(fi1000comCGdata1)
dim(fi1000comCGdata1)
fi1000comCGdata2 <- fi1000comCGdata1
head(fi1000comCGdata2)
dim(fi1000comCGdata2)
fi1000comCGdata2_line <- data.frame(fi1000comCGdata2[,-c(2)])
head(fi1000comCGdata2_line)
fi1000comCGdata2_line["Row"] <- paste0(fi1000comCGdata2_line$DMR,"%",fi1000comCGdata2_line$chrc,"%",fi1000comCGdata2_line$startc,"%",fi1000comCGdata2_line$endc,"%",fi1000comCGdata2_line$id)
fi1000comCGdata2_line2 <- fi1000comCGdata2_line
rownames(fi1000comCGdata2_line2) <- fi1000comCGdata2_line2$Row
head(fi1000comCGdata2_line2)
dim(fi1000comCGdata2_line2)
fi1000comCGdata2_line2 <- fi1000comCGdata2_line2[,-c(1:6,24)]
head(fi1000comCGdata2_line2)
fi1000comCGdata2_line2 <- as.matrix(fi1000comCGdata2_line2)
fi1000comCGdata2_line2st <- stack(fi1000comCGdata2_line2)
head(fi1000comCGdata2_line2st)
fi1000comCGdata2_line2st <- data.frame(fi1000comCGdata2_line2st)
fi1000comCGdata2_line2st <- fi1000comCGdata2_line2st[,c(2,1,3)]
colnames(fi1000comCGdata2_line2st) <- c("Group", "Csites", "Value")
head(fi1000comCGdata2_line2st)


#PPEL
fi1000comCGdata2_line2st_PPEL <- fi1000comCGdata2_line2st[grep("PPEL", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_PPEL.sep <- cSplit(fi1000comCGdata2_line2st_PPEL, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_PPEL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(39558954, 39559868), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_PPEL.sep.svg", width=25, height=5, units="cm", dpi=96)
#DIRAS3_TSS
fi1000comCGdata2_line2st_DIRAS3_TSS <- fi1000comCGdata2_line2st[grep("DIRAS3_TSS", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_DIRAS3_TSS.sep <- cSplit(fi1000comCGdata2_line2st_DIRAS3_TSS, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_DIRAS3_TSS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(68049750, 68051862), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_DIRAS3_TSS.sep.svg", width=25, height=5, units="cm", dpi=96)
#DIRAS3_ex2
fi1000comCGdata2_line2st_DIRAS3_ex2 <- fi1000comCGdata2_line2st[grep("DIRAS3_ex2", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_DIRAS3_ex2.sep <- cSplit(fi1000comCGdata2_line2st_DIRAS3_ex2, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_DIRAS3_ex2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(68046822, 68047803), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_DIRAS3_ex2.sep.svg", width=25, height=5, units="cm", dpi=96)
#GPR1_AS
fi1000comCGdata2_line2st_GPR1_AS <- fi1000comCGdata2_line2st[grep("GPR1-AS", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_GPR1_AS.sep <- cSplit(fi1000comCGdata2_line2st_GPR1_AS, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_GPR1_AS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(206202243, 206204721), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_GPR1_AS.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZDBF2_GPR1
fi1000comCGdata2_line2st_ZDBF2_GPR1 <- fi1000comCGdata2_line2st[grep("ZDBF2/GPR1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_ZDBF2_GPR1.sep <- cSplit(fi1000comCGdata2_line2st_ZDBF2_GPR1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_ZDBF2_GPR1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(206249859, 206271820), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_ZDBF2_GPR1.sep.svg", width=25, height=5, units="cm", dpi=96)
#NAP1L5
fi1000comCGdata2_line2st_NAP1L5 <- fi1000comCGdata2_line2st[grep("NAP1L5", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_NAP1L5.sep <- cSplit(fi1000comCGdata2_line2st_NAP1L5, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_NAP1L5.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(88697033, 88698086), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_NAP1L5.sep.svg", width=25, height=5, units="cm", dpi=96)
#VTRNA2_1
fi1000comCGdata2_line2st_VTRNA2_1 <- fi1000comCGdata2_line2st[grep("VTRNA2-1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_VTRNA2_1.sep <- cSplit(fi1000comCGdata2_line2st_VTRNA2_1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_VTRNA2_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(136079113, 136080956), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_VTRNA2_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#FAM50B
fi1000comCGdata2_line2st_FAM50B <- fi1000comCGdata2_line2st[grep("FAM50B", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_FAM50B.sep <- cSplit(fi1000comCGdata2_line2st_FAM50B, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_FAM50B.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3848848, 3850125), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_FAM50B.sep.svg", width=25, height=5, units="cm", dpi=96)
#PLAGL1
fi1000comCGdata2_line2st_PLAGL1 <- fi1000comCGdata2_line2st[grep("PLAGL1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_PLAGL1.sep <- cSplit(fi1000comCGdata2_line2st_PLAGL1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_PLAGL1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(144006941, 144008751), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_PLAGL1.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2R
fi1000comCGdata2_line2st_IGF2R <- fi1000comCGdata2_line2st[grep("IGF2R", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_IGF2R.sep <- cSplit(fi1000comCGdata2_line2st_IGF2R, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_IGF2R.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(160005526, 160006529), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_IGF2R.sep.svg", width=25, height=5, units="cm", dpi=96)
#VDR27
fi1000comCGdata2_line2st_VDR27 <- fi1000comCGdata2_line2st[grep("VDR27", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_VDR27.sep <- cSplit(fi1000comCGdata2_line2st_VDR27, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_VDR27.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(169654408, 169655522), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_VDR27.sep.svg", width=25, height=5, units="cm", dpi=96)
#GRB10
fi1000comCGdata2_line2st_GRB10 <- fi1000comCGdata2_line2st[grep("GRB10", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_GRB10.sep <- cSplit(fi1000comCGdata2_line2st_GRB10, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_GRB10.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(50781029, 50783615), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_GRB10.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG10
fi1000comCGdata2_line2st_PEG10 <- fi1000comCGdata2_line2st[grep("PEG10", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_PEG10.sep <- cSplit(fi1000comCGdata2_line2st_PEG10, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_PEG10.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(94656225, 94658648), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_PEG10.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEST
fi1000comCGdata2_line2st_MEST <- fi1000comCGdata2_line2st[grep("MEST", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_MEST.sep <- cSplit(fi1000comCGdata2_line2st_MEST, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_MEST.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(130490281, 130494547), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_MEST.sep.svg", width=25, height=5, units="cm", dpi=96)
#SVOPL
fi1000comCGdata2_line2st_SVOPL <- fi1000comCGdata2_line2st[grep("SVOPL", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_SVOPL.sep <- cSplit(fi1000comCGdata2_line2st_SVOPL, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_SVOPL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(138663373, 138664324), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_SVOPL.sep.svg", width=25, height=5, units="cm", dpi=96)
#HTR5A
fi1000comCGdata2_line2st_HTR5A <- fi1000comCGdata2_line2st[grep("HTR5A", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_HTR5A.sep <- cSplit(fi1000comCGdata2_line2st_HTR5A, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_HTR5A.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(155071009, 155071672), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_HTR5A.sep.png", width=25, height=5, units="cm", dpi=96)
#ERILN2
fi1000comCGdata2_line2st_ERILN2 <- fi1000comCGdata2_line2st[grep("ERILN2", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_ERILN2.sep <- cSplit(fi1000comCGdata2_line2st_ERILN2, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_ERILN2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(37747474, 37748570), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_ERILN2.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG13
fi1000comCGdata2_line2st_PEG13 <- fi1000comCGdata2_line2st[grep("PEG13", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_PEG13.sep <- cSplit(fi1000comCGdata2_line2st_PEG13, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_PEG13.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(140098048, 140100982), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_PEG13.sep.svg", width=25, height=5, units="cm", dpi=96)
#FANCC
fi1000comCGdata2_line2st_FANCC <- fi1000comCGdata2_line2st[grep("FANCC", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_FANCC.sep <- cSplit(fi1000comCGdata2_line2st_FANCC, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_FANCC.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(95313118, 95313462), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_FANCC.sep.svg", width=25, height=5, units="cm", dpi=96)
#INPP5F
fi1000comCGdata2_line2st_INPP5F <- fi1000comCGdata2_line2st[grep("INPP5F", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_INPP5F.sep <- cSplit(fi1000comCGdata2_line2st_INPP5F, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_INPP5F.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(119818534, 119819215), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_INPP5F.sep.svg", width=25, height=5, units="cm", dpi=96)
#H19_IGF2
fi1000comCGdata2_line2st_H19_IGF2 <- fi1000comCGdata2_line2st[grep("H19/IGF2", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_H19_IGF2.sep <- cSplit(fi1000comCGdata2_line2st_H19_IGF2, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_H19_IGF2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(1997582, 2003510), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_H19_IGF2.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2_ex9
fi1000comCGdata2_line2st_IGF2_ex9 <- fi1000comCGdata2_line2st[grep("IGF2_ex9", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_IGF2_ex9.sep <- cSplit(fi1000comCGdata2_line2st_IGF2_ex9, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_IGF2_ex9.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2132761, 2133882), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_IGF2_ex9.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2_altTSS
fi1000comCGdata2_line2st_IGF2_altTSS <- fi1000comCGdata2_line2st[grep("IGF2_altTSS", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_IGF2_altTSS.sep <- cSplit(fi1000comCGdata2_line2st_IGF2_altTSS, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_IGF2_altTSS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2147103, 2148538), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_IGF2_altTSS.sep.svg", width=25, height=5, units="cm", dpi=96)
#KCNQ1OT1
fi1000comCGdata2_line2st_KCNQ1OT1 <- fi1000comCGdata2_line2st[grep("KCNQ1OT1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_KCNQ1OT1.sep <- cSplit(fi1000comCGdata2_line2st_KCNQ1OT1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_KCNQ1OT1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2698718, 2701029), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_KCNQ1OT1.sep.svg", width=25, height=5, units="cm", dpi=96)
#RB1
fi1000comCGdata2_line2st_RB1 <- fi1000comCGdata2_line2st[grep("RB1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_RB1.sep <- cSplit(fi1000comCGdata2_line2st_RB1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_RB1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(48318205, 48321627), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_RB1.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG3_DLK1
fi1000comCGdata2_line2st_MEG3_DLK1 <- fi1000comCGdata2_line2st[grep("MEG3/DLK1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_MEG3_DLK1.sep <- cSplit(fi1000comCGdata2_line2st_MEG3_DLK1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_MEG3_DLK1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100809090, 100811721), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_MEG3_DLK1.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG3
fi1000comCGdata2_line2st_MEG3 <- fi1000comCGdata2_line2st[grep("MEG3", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_MEG3.sep <- cSplit(fi1000comCGdata2_line2st_MEG3, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_MEG3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100824187, 100827641), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_MEG3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG8
fi1000comCGdata2_line2st_MEG8 <- fi1000comCGdata2_line2st[grep("MEG8", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_MEG8.sep <- cSplit(fi1000comCGdata2_line2st_MEG8, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_MEG8.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100904404, 100905082), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_MEG8.sep.svg", width=25, height=5, units="cm", dpi=96)
#MKRN3
fi1000comCGdata2_line2st_MKRN3 <- fi1000comCGdata2_line2st[grep("MKRN3", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_MKRN3.sep <- cSplit(fi1000comCGdata2_line2st_MKRN3, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_MKRN3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23561939, 23567348), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_MKRN3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MAGEL2
fi1000comCGdata2_line2st_MAGEL2 <- fi1000comCGdata2_line2st[grep("MAGEL2", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_MAGEL2.sep <- cSplit(fi1000comCGdata2_line2st_MAGEL2, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_MAGEL2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23647278, 23648882), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_MAGEL2.sep.svg", width=25, height=5, units="cm", dpi=96)
#NDN
fi1000comCGdata2_line2st_NDN <- fi1000comCGdata2_line2st[grep("NDN", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_NDN.sep <- cSplit(fi1000comCGdata2_line2st_NDN, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_NDN.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23686304, 23687612), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_NDN.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_alt
fi1000comCGdata2_line2st_SNRPN_alt <- fi1000comCGdata2_line2st[grep("SNRPN_alt", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_SNRPN_alt.sep <- cSplit(fi1000comCGdata2_line2st_SNRPN_alt, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_SNRPN_alt.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24823417, 24824334), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_SNRPN_alt.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_int1
fi1000comCGdata2_line2st_SNRPN_int1 <- fi1000comCGdata2_line2st[grep("SNRPN_int1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_SNRPN_int1.sep <- cSplit(fi1000comCGdata2_line2st_SNRPN_int1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_SNRPN_int1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24847861, 24848682), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_SNRPN_int1.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_int2
fi1000comCGdata2_line2st_SNRPN_int2 <- fi1000comCGdata2_line2st[grep("SNRPN_int2", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_SNRPN_int2.sep <- cSplit(fi1000comCGdata2_line2st_SNRPN_int2, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_SNRPN_int2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24877880, 24878758), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_SNRPN_int2.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNURF
fi1000comCGdata2_line2st_SNURF <- fi1000comCGdata2_line2st[grep("SNURF", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_SNURF.sep <- cSplit(fi1000comCGdata2_line2st_SNURF, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_SNURF.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24954857, 24956829), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_SNURF.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF1R
fi1000comCGdata2_line2st_IGF1R <- fi1000comCGdata2_line2st[grep("IGF1R", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_IGF1R.sep <- cSplit(fi1000comCGdata2_line2st_IGF1R, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_IGF1R.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(98865267, 98866421), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_IGF1R.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF597_1
fi1000comCGdata2_line2st_ZNF597_1 <- fi1000comCGdata2_line2st[grep("ZNF597_1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_ZNF597_1.sep <- cSplit(fi1000comCGdata2_line2st_ZNF597_1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_ZNF597_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3431801, 3432388), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_ZNF597_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF597_tss
fi1000comCGdata2_line2st_ZNF597_tss <- fi1000comCGdata2_line2st[grep("ZNF597_tss", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_ZNF597_tss.sep <- cSplit(fi1000comCGdata2_line2st_ZNF597_tss, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_ZNF597_tss.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3442828, 3444463), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_ZNF597_tss.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF331_1
fi1000comCGdata2_line2st_ZNF331_1 <- fi1000comCGdata2_line2st[grep("ZNF331_1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_ZNF331_1.sep <- cSplit(fi1000comCGdata2_line2st_ZNF331_1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_ZNF331_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(53537256, 53538958), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_ZNF331_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF331_2
fi1000comCGdata2_line2st_ZNF331_2 <- fi1000comCGdata2_line2st[grep("ZNF331_2", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_ZNF331_2.sep <- cSplit(fi1000comCGdata2_line2st_ZNF331_2, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_ZNF331_2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(53553832, 53555171), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_ZNF331_2.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG3
fi1000comCGdata2_line2st_PEG3 <- fi1000comCGdata2_line2st[grep("PEG3", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_PEG3.sep <- cSplit(fi1000comCGdata2_line2st_PEG3, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_PEG3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(56837125, 56841903), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_PEG3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MCTS2P
fi1000comCGdata2_line2st_MCTS2P <- fi1000comCGdata2_line2st[grep("MCTS2P", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_MCTS2P.sep <- cSplit(fi1000comCGdata2_line2st_MCTS2P, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_MCTS2P.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(31546860, 31548130), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_MCTS2P.sep.svg", width=25, height=5, units="cm", dpi=96)
#NNAT
fi1000comCGdata2_line2st_NNAT <- fi1000comCGdata2_line2st[grep("NNAT", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_NNAT.sep <- cSplit(fi1000comCGdata2_line2st_NNAT, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_NNAT.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(37520202, 37522126), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_NNAT.sep.svg", width=25, height=5, units="cm", dpi=96)
#L3MBTL1
fi1000comCGdata2_line2st_L3MBTL1 <- fi1000comCGdata2_line2st[grep("L3MBTL1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_L3MBTL1.sep <- cSplit(fi1000comCGdata2_line2st_L3MBTL1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_L3MBTL1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(43513725, 43515400), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_L3MBTL1.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_NESP
fi1000comCGdata2_line2st_GNAS_NESP <- fi1000comCGdata2_line2st[grep("GNAS-NESP", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_GNAS_NESP.sep <- cSplit(fi1000comCGdata2_line2st_GNAS_NESP, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_GNAS_NESP.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58838984, 58843557), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_GNAS_NESP.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_AS1
fi1000comCGdata2_line2st_GNAS_AS1 <- fi1000comCGdata2_line2st[grep("GNAS-AS1", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_GNAS_AS1.sep <- cSplit(fi1000comCGdata2_line2st_GNAS_AS1, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_GNAS_AS1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58850594, 58852978), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_GNAS_AS1.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_XL
fi1000comCGdata2_line2st_GNAS_XL <- fi1000comCGdata2_line2st[grep("GNAS-XL", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_GNAS_XL.sep <- cSplit(fi1000comCGdata2_line2st_GNAS_XL, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_GNAS_XL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58853850, 58856408), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_GNAS_XL.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNASA_B
fi1000comCGdata2_line2st_GNASA_B <- fi1000comCGdata2_line2st[grep("GNASA/B", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_GNASA_B.sep <- cSplit(fi1000comCGdata2_line2st_GNASA_B, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_GNASA_B.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58888210, 58890146), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_GNASA_B.sep.svg", width=25, height=5, units="cm", dpi=96)
#WRB
fi1000comCGdata2_line2st_WRB <- fi1000comCGdata2_line2st[grep("WRB", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_WRB.sep <- cSplit(fi1000comCGdata2_line2st_WRB, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_WRB.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(39385584, 39386350), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_WRB.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNU13
fi1000comCGdata2_line2st_SNU13 <- fi1000comCGdata2_line2st[grep("SNU13", fi1000comCGdata2_line2st$Csites),]
fi1000comCGdata2_line2st_SNU13.sep <- cSplit(fi1000comCGdata2_line2st_SNU13, "Csites", "%")
ggplot(fi1000comCGdata2_line2st_SNU13.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(41681770, 41682869), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi1000comCGdata2_line2st_SNU13.sep.svg", width=25, height=5, units="cm", dpi=96)





############################## Flank Line chart ############################

awk '{print $1"\t"$2-2000"\t"$3+2000"\t"$4"\t"$5}' /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed  > Monk_human_ICR_hg38_flankplusminus2000.bed

bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b Monk_human_ICR_hg38_flankplusminus2000.bed > MERGE_myicomCGbatCG3re_human_ICR_flanked2000.txt
MERGE_myicomCGbatCG3re_human_ICR_flanked2000 <- read.table("MERGE_myicomCGbatCG3re_human_ICR_flanked2000.txt", header = F)
head(MERGE_myicomCGbatCG3re_human_ICR_flanked2000)
dim(MERGE_myicomCGbatCG3re_human_ICR_flanked2000)
head(MERGE_myicomCGbatCG3re_human_ICR_flanked2000)
colnames(MERGE_myicomCGbatCG3re_human_ICR_flanked2000) <- c("chrc","startc","endc","id", colnames(MERGE_myiCombat3reavg_pos_chr[,5:22]),"chr","start","end","DMR", "DMRType")
fi2000comCGdata1 <- MERGE_myicomCGbatCG3re_human_ICR_flanked2000[,c((ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked2000)-1),ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked2000), 1:(ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked2000)-5))]
rownames(fi2000comCGdata1)
colnames(fi2000comCGdata1)
head(fi2000comCGdata1)
dim(fi2000comCGdata1)
fi2000comCGdata2 <- fi2000comCGdata1
head(fi2000comCGdata2)
dim(fi2000comCGdata2)
fi2000comCGdata2_line <- data.frame(fi2000comCGdata2[,-c(2)])
head(fi2000comCGdata2_line)
fi2000comCGdata2_line["Row"] <- paste0(fi2000comCGdata2_line$DMR,"%",fi2000comCGdata2_line$chrc,"%",fi2000comCGdata2_line$startc,"%",fi2000comCGdata2_line$endc,"%",fi2000comCGdata2_line$id)
fi2000comCGdata2_line2 <- fi2000comCGdata2_line
rownames(fi2000comCGdata2_line2) <- fi2000comCGdata2_line2$Row
head(fi2000comCGdata2_line2)
fi2000comCGdata2_line2 <- fi2000comCGdata2_line2[,-c(1:6,24)]
head(fi2000comCGdata2_line2)
fi2000comCGdata2_line2 <- as.matrix(fi2000comCGdata2_line2)
fi2000comCGdata2_line2st <- stack(fi2000comCGdata2_line2)
head(fi2000comCGdata2_line2st)
fi2000comCGdata2_line2st <- data.frame(fi2000comCGdata2_line2st)
fi2000comCGdata2_line2st <- fi2000comCGdata2_line2st[,c(2,1,3)]
colnames(fi2000comCGdata2_line2st) <- c("Group", "Csites", "Value")
head(fi2000comCGdata2_line2st)


#PPEL
fi2000comCGdata2_line2st_PPEL <- fi2000comCGdata2_line2st[grep("PPEL", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_PPEL.sep <- cSplit(fi2000comCGdata2_line2st_PPEL, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_PPEL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(39558954, 39559868), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_PPEL.sep.svg", width=25, height=5, units="cm", dpi=96)
#DIRAS3_TSS
fi2000comCGdata2_line2st_DIRAS3_TSS <- fi2000comCGdata2_line2st[grep("DIRAS3_TSS", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_DIRAS3_TSS.sep <- cSplit(fi2000comCGdata2_line2st_DIRAS3_TSS, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_DIRAS3_TSS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(68049750, 68051862), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_DIRAS3_TSS.sep.svg", width=25, height=5, units="cm", dpi=96)
#DIRAS3_ex2
fi2000comCGdata2_line2st_DIRAS3_ex2 <- fi2000comCGdata2_line2st[grep("DIRAS3_ex2", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_DIRAS3_ex2.sep <- cSplit(fi2000comCGdata2_line2st_DIRAS3_ex2, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_DIRAS3_ex2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(68046822, 68047803), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_DIRAS3_ex2.sep.svg", width=25, height=5, units="cm", dpi=96)
#GPR1_AS
fi2000comCGdata2_line2st_GPR1_AS <- fi2000comCGdata2_line2st[grep("GPR1-AS", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_GPR1_AS.sep <- cSplit(fi2000comCGdata2_line2st_GPR1_AS, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_GPR1_AS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(206202243, 206204721), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_GPR1_AS.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZDBF2_GPR1
fi2000comCGdata2_line2st_ZDBF2_GPR1 <- fi2000comCGdata2_line2st[grep("ZDBF2/GPR1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_ZDBF2_GPR1.sep <- cSplit(fi2000comCGdata2_line2st_ZDBF2_GPR1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_ZDBF2_GPR1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(206249859, 206271820), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_ZDBF2_GPR1.sep.svg", width=25, height=5, units="cm", dpi=96)
#NAP1L5
fi2000comCGdata2_line2st_NAP1L5 <- fi2000comCGdata2_line2st[grep("NAP1L5", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_NAP1L5.sep <- cSplit(fi2000comCGdata2_line2st_NAP1L5, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_NAP1L5.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(88697033, 88698086), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_NAP1L5.sep.svg", width=25, height=5, units="cm", dpi=96)
#VTRNA2_1
fi2000comCGdata2_line2st_VTRNA2_1 <- fi2000comCGdata2_line2st[grep("VTRNA2-1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_VTRNA2_1.sep <- cSplit(fi2000comCGdata2_line2st_VTRNA2_1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_VTRNA2_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(136079113, 136080956), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_VTRNA2_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#FAM50B
fi2000comCGdata2_line2st_FAM50B <- fi2000comCGdata2_line2st[grep("FAM50B", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_FAM50B.sep <- cSplit(fi2000comCGdata2_line2st_FAM50B, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_FAM50B.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3848848, 3850125), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_FAM50B.sep.svg", width=25, height=5, units="cm", dpi=96)
#PLAGL1
fi2000comCGdata2_line2st_PLAGL1 <- fi2000comCGdata2_line2st[grep("PLAGL1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_PLAGL1.sep <- cSplit(fi2000comCGdata2_line2st_PLAGL1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_PLAGL1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(144006941, 144008751), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_PLAGL1.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2R
fi2000comCGdata2_line2st_IGF2R <- fi2000comCGdata2_line2st[grep("IGF2R", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_IGF2R.sep <- cSplit(fi2000comCGdata2_line2st_IGF2R, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_IGF2R.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(160005526, 160006529), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_IGF2R.sep.svg", width=25, height=5, units="cm", dpi=96)
#VDR27
fi2000comCGdata2_line2st_VDR27 <- fi2000comCGdata2_line2st[grep("VDR27", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_VDR27.sep <- cSplit(fi2000comCGdata2_line2st_VDR27, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_VDR27.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(169654408, 169655522), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_VDR27.sep.svg", width=25, height=5, units="cm", dpi=96)
#GRB10
fi2000comCGdata2_line2st_GRB10 <- fi2000comCGdata2_line2st[grep("GRB10", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_GRB10.sep <- cSplit(fi2000comCGdata2_line2st_GRB10, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_GRB10.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(50781029, 50783615), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_GRB10.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG10
fi2000comCGdata2_line2st_PEG10 <- fi2000comCGdata2_line2st[grep("PEG10", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_PEG10.sep <- cSplit(fi2000comCGdata2_line2st_PEG10, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_PEG10.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(94656225, 94658648), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_PEG10.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEST
fi2000comCGdata2_line2st_MEST <- fi2000comCGdata2_line2st[grep("MEST", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_MEST.sep <- cSplit(fi2000comCGdata2_line2st_MEST, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_MEST.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(130490281, 130494547), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_MEST.sep.svg", width=25, height=5, units="cm", dpi=96)
#SVOPL
fi2000comCGdata2_line2st_SVOPL <- fi2000comCGdata2_line2st[grep("SVOPL", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_SVOPL.sep <- cSplit(fi2000comCGdata2_line2st_SVOPL, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_SVOPL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(138663373, 138664324), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_SVOPL.sep.svg", width=25, height=5, units="cm", dpi=96)
#HTR5A
fi2000comCGdata2_line2st_HTR5A <- fi2000comCGdata2_line2st[grep("HTR5A", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_HTR5A.sep <- cSplit(fi2000comCGdata2_line2st_HTR5A, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_HTR5A.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(155071009, 155071672), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_HTR5A.sep.svg", width=25, height=5, units="cm", dpi=96)
#ERILN2
fi2000comCGdata2_line2st_ERILN2 <- fi2000comCGdata2_line2st[grep("ERILN2", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_ERILN2.sep <- cSplit(fi2000comCGdata2_line2st_ERILN2, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_ERILN2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(37747474, 37748570), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_ERILN2.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG13
fi2000comCGdata2_line2st_PEG13 <- fi2000comCGdata2_line2st[grep("PEG13", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_PEG13.sep <- cSplit(fi2000comCGdata2_line2st_PEG13, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_PEG13.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(140098048, 140100982), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_PEG13.sep.svg", width=25, height=5, units="cm", dpi=96)
#FANCC
fi2000comCGdata2_line2st_FANCC <- fi2000comCGdata2_line2st[grep("FANCC", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_FANCC.sep <- cSplit(fi2000comCGdata2_line2st_FANCC, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_FANCC.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(95313118, 95313462), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_FANCC.sep.svg", width=25, height=5, units="cm", dpi=96)
#INPP5F
fi2000comCGdata2_line2st_INPP5F <- fi2000comCGdata2_line2st[grep("INPP5F", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_INPP5F.sep <- cSplit(fi2000comCGdata2_line2st_INPP5F, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_INPP5F.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(119818534, 119819215), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_INPP5F.sep.svg", width=25, height=5, units="cm", dpi=96)
#H19_IGF2
fi2000comCGdata2_line2st_H19_IGF2 <- fi2000comCGdata2_line2st[grep("H19/IGF2", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_H19_IGF2.sep <- cSplit(fi2000comCGdata2_line2st_H19_IGF2, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_H19_IGF2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(1997582, 2003510), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_H19_IGF2.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2_ex9
fi2000comCGdata2_line2st_IGF2_ex9 <- fi2000comCGdata2_line2st[grep("IGF2_ex9", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_IGF2_ex9.sep <- cSplit(fi2000comCGdata2_line2st_IGF2_ex9, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_IGF2_ex9.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2132761, 2133882), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_IGF2_ex9.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2_altTSS
fi2000comCGdata2_line2st_IGF2_altTSS <- fi2000comCGdata2_line2st[grep("IGF2_altTSS", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_IGF2_altTSS.sep <- cSplit(fi2000comCGdata2_line2st_IGF2_altTSS, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_IGF2_altTSS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2147103, 2148538), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_IGF2_altTSS.sep.svg", width=25, height=5, units="cm", dpi=96)
#KCNQ1OT1
fi2000comCGdata2_line2st_KCNQ1OT1 <- fi2000comCGdata2_line2st[grep("KCNQ1OT1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_KCNQ1OT1.sep <- cSplit(fi2000comCGdata2_line2st_KCNQ1OT1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_KCNQ1OT1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2698718, 2701029), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_KCNQ1OT1.sep.svg", width=25, height=5, units="cm", dpi=96)
#RB1
fi2000comCGdata2_line2st_RB1 <- fi2000comCGdata2_line2st[grep("RB1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_RB1.sep <- cSplit(fi2000comCGdata2_line2st_RB1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_RB1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(48318205, 48321627), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_RB1.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG3_DLK1
fi2000comCGdata2_line2st_MEG3_DLK1 <- fi2000comCGdata2_line2st[grep("MEG3/DLK1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_MEG3_DLK1.sep <- cSplit(fi2000comCGdata2_line2st_MEG3_DLK1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_MEG3_DLK1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100809090, 100811721), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_MEG3_DLK1.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG3
fi2000comCGdata2_line2st_MEG3 <- fi2000comCGdata2_line2st[grep("MEG3", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_MEG3.sep <- cSplit(fi2000comCGdata2_line2st_MEG3, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_MEG3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100824187, 100827641), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_MEG3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG8
fi2000comCGdata2_line2st_MEG8 <- fi2000comCGdata2_line2st[grep("MEG8", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_MEG8.sep <- cSplit(fi2000comCGdata2_line2st_MEG8, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_MEG8.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100904404, 100905082), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_MEG8.sep.svg", width=25, height=5, units="cm", dpi=96)
#MKRN3
fi2000comCGdata2_line2st_MKRN3 <- fi2000comCGdata2_line2st[grep("MKRN3", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_MKRN3.sep <- cSplit(fi2000comCGdata2_line2st_MKRN3, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_MKRN3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23561939, 23567348), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_MKRN3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MAGEL2
fi2000comCGdata2_line2st_MAGEL2 <- fi2000comCGdata2_line2st[grep("MAGEL2", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_MAGEL2.sep <- cSplit(fi2000comCGdata2_line2st_MAGEL2, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_MAGEL2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23647278, 23648882), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_MAGEL2.sep.svg", width=25, height=5, units="cm", dpi=96)
#NDN
fi2000comCGdata2_line2st_NDN <- fi2000comCGdata2_line2st[grep("NDN", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_NDN.sep <- cSplit(fi2000comCGdata2_line2st_NDN, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_NDN.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23686304, 23687612), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_NDN.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_alt
fi2000comCGdata2_line2st_SNRPN_alt <- fi2000comCGdata2_line2st[grep("SNRPN_alt", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_SNRPN_alt.sep <- cSplit(fi2000comCGdata2_line2st_SNRPN_alt, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_SNRPN_alt.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24823417, 24824334), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_SNRPN_alt.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_int1
fi2000comCGdata2_line2st_SNRPN_int1 <- fi2000comCGdata2_line2st[grep("SNRPN_int1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_SNRPN_int1.sep <- cSplit(fi2000comCGdata2_line2st_SNRPN_int1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_SNRPN_int1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24847861, 24848682), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_SNRPN_int1.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_int2
fi2000comCGdata2_line2st_SNRPN_int2 <- fi2000comCGdata2_line2st[grep("SNRPN_int2", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_SNRPN_int2.sep <- cSplit(fi2000comCGdata2_line2st_SNRPN_int2, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_SNRPN_int2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24877880, 24878758), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_SNRPN_int2.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNURF
fi2000comCGdata2_line2st_SNURF <- fi2000comCGdata2_line2st[grep("SNURF", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_SNURF.sep <- cSplit(fi2000comCGdata2_line2st_SNURF, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_SNURF.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24954857, 24956829), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_SNURF.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF1R
fi2000comCGdata2_line2st_IGF1R <- fi2000comCGdata2_line2st[grep("IGF1R", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_IGF1R.sep <- cSplit(fi2000comCGdata2_line2st_IGF1R, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_IGF1R.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(98865267, 98866421), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_IGF1R.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF597_1
fi2000comCGdata2_line2st_ZNF597_1 <- fi2000comCGdata2_line2st[grep("ZNF597_1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_ZNF597_1.sep <- cSplit(fi2000comCGdata2_line2st_ZNF597_1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_ZNF597_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3431801, 3432388), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_ZNF597_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF597_tss
fi2000comCGdata2_line2st_ZNF597_tss <- fi2000comCGdata2_line2st[grep("ZNF597_tss", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_ZNF597_tss.sep <- cSplit(fi2000comCGdata2_line2st_ZNF597_tss, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_ZNF597_tss.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3442828, 3444463), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_ZNF597_tss.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF331_1
fi2000comCGdata2_line2st_ZNF331_1 <- fi2000comCGdata2_line2st[grep("ZNF331_1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_ZNF331_1.sep <- cSplit(fi2000comCGdata2_line2st_ZNF331_1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_ZNF331_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(53537256, 53538958), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_ZNF331_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF331_2
fi2000comCGdata2_line2st_ZNF331_2 <- fi2000comCGdata2_line2st[grep("ZNF331_2", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_ZNF331_2.sep <- cSplit(fi2000comCGdata2_line2st_ZNF331_2, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_ZNF331_2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(53553832, 53555171), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_ZNF331_2.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG3
fi2000comCGdata2_line2st_PEG3 <- fi2000comCGdata2_line2st[grep("PEG3", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_PEG3.sep <- cSplit(fi2000comCGdata2_line2st_PEG3, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_PEG3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(56837125, 56841903), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_PEG3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MCTS2P
fi2000comCGdata2_line2st_MCTS2P <- fi2000comCGdata2_line2st[grep("MCTS2P", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_MCTS2P.sep <- cSplit(fi2000comCGdata2_line2st_MCTS2P, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_MCTS2P.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(31546860, 31548130), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_MCTS2P.sep.svg", width=25, height=5, units="cm", dpi=96)
#NNAT
fi2000comCGdata2_line2st_NNAT <- fi2000comCGdata2_line2st[grep("NNAT", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_NNAT.sep <- cSplit(fi2000comCGdata2_line2st_NNAT, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_NNAT.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(37520202, 37522126), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_NNAT.sep.svg", width=25, height=5, units="cm", dpi=96)
#L3MBTL1
fi2000comCGdata2_line2st_L3MBTL1 <- fi2000comCGdata2_line2st[grep("L3MBTL1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_L3MBTL1.sep <- cSplit(fi2000comCGdata2_line2st_L3MBTL1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_L3MBTL1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(43513725, 43515400), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_L3MBTL1.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_NESP
fi2000comCGdata2_line2st_GNAS_NESP <- fi2000comCGdata2_line2st[grep("GNAS-NESP", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_GNAS_NESP.sep <- cSplit(fi2000comCGdata2_line2st_GNAS_NESP, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_GNAS_NESP.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58838984, 58843557), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_GNAS_NESP.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_AS1
fi2000comCGdata2_line2st_GNAS_AS1 <- fi2000comCGdata2_line2st[grep("GNAS-AS1", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_GNAS_AS1.sep <- cSplit(fi2000comCGdata2_line2st_GNAS_AS1, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_GNAS_AS1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58850594, 58852978), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_GNAS_AS1.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_XL
fi2000comCGdata2_line2st_GNAS_XL <- fi2000comCGdata2_line2st[grep("GNAS-XL", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_GNAS_XL.sep <- cSplit(fi2000comCGdata2_line2st_GNAS_XL, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_GNAS_XL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58853850, 58856408), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_GNAS_XL.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNASA_B
fi2000comCGdata2_line2st_GNASA_B <- fi2000comCGdata2_line2st[grep("GNASA/B", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_GNASA_B.sep <- cSplit(fi2000comCGdata2_line2st_GNASA_B, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_GNASA_B.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58888210, 58890146), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_GNASA_B.sep.svg", width=25, height=5, units="cm", dpi=96)
#WRB
fi2000comCGdata2_line2st_WRB <- fi2000comCGdata2_line2st[grep("WRB", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_WRB.sep <- cSplit(fi2000comCGdata2_line2st_WRB, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_WRB.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(39385584, 39386350), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_WRB.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNU13
fi2000comCGdata2_line2st_SNU13 <- fi2000comCGdata2_line2st[grep("SNU13", fi2000comCGdata2_line2st$Csites),]
fi2000comCGdata2_line2st_SNU13.sep <- cSplit(fi2000comCGdata2_line2st_SNU13, "Csites", "%")
ggplot(fi2000comCGdata2_line2st_SNU13.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(41681770, 41682869), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi2000comCGdata2_line2st_SNU13.sep.svg", width=25, height=5, units="cm", dpi=96)






############################## Flank Line chart ############################

awk '{print $1"\t"$2-3000"\t"$3+3000"\t"$4"\t"$5}' /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed  > Monk_human_ICR_hg38_flankplusminus3000.bed

bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b Monk_human_ICR_hg38_flankplusminus3000.bed > MERGE_myicomCGbatCG3re_human_ICR_flanked3000.txt
MERGE_myicomCGbatCG3re_human_ICR_flanked3000 <- read.table("MERGE_myicomCGbatCG3re_human_ICR_flanked3000.txt", header = F)
head(MERGE_myicomCGbatCG3re_human_ICR_flanked3000)
dim(MERGE_myicomCGbatCG3re_human_ICR_flanked3000)
head(MERGE_myicomCGbatCG3re_human_ICR_flanked3000)
colnames(MERGE_myicomCGbatCG3re_human_ICR_flanked3000) <- c("chrc","startc","endc","id", colnames(MERGE_myiCombat3reavg_pos_chr[,5:22]),"chr","start","end","DMR", "DMRType")
fi3000comCGdata1 <- MERGE_myicomCGbatCG3re_human_ICR_flanked3000[,c((ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked3000)-1),ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked3000), 1:(ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked3000)-5))]
rownames(fi3000comCGdata1)
colnames(fi3000comCGdata1)
head(fi3000comCGdata1)
dim(fi3000comCGdata1)
Count_fi3000comCGdataICR1 <- count(fi3000comCGdata1, "DMR")
head(Count_fi3000comCGdataICR1)
dim(Count_fi3000comCGdataICR1)
#Convert meth percent to beta like 
fi3000comCGdata2 <- fi3000comCGdata1
head(fi3000comCGdata2)
dim(fi3000comCGdata2)
fi3000comCGdata2_line <- data.frame(fi3000comCGdata2[,-c(2)])
head(fi3000comCGdata2_line)
fi3000comCGdata2_line["Row"] <- paste0(fi3000comCGdata2_line$DMR,"%",fi3000comCGdata2_line$chrc,"%",fi3000comCGdata2_line$startc,"%",fi3000comCGdata2_line$endc,"%",fi3000comCGdata2_line$id)
fi3000comCGdata2_line2 <- fi3000comCGdata2_line
rownames(fi3000comCGdata2_line2) <- fi3000comCGdata2_line2$Row
head(fi3000comCGdata2_line2)
fi3000comCGdata2_line2 <- fi3000comCGdata2_line2[,-c(1:6,24)]
head(fi3000comCGdata2_line2)
fi3000comCGdata2_line2 <- as.matrix(fi3000comCGdata2_line2)
fi3000comCGdata2_line2st <- stack(fi3000comCGdata2_line2)
head(fi3000comCGdata2_line2st)
fi3000comCGdata2_line2st <- data.frame(fi3000comCGdata2_line2st)
fi3000comCGdata2_line2st <- fi3000comCGdata2_line2st[,c(2,1,3)]
colnames(fi3000comCGdata2_line2st) <- c("Group", "Csites", "Value")
head(fi3000comCGdata2_line2st)


#PPEL
fi3000comCGdata2_line2st_PPEL <- fi3000comCGdata2_line2st[grep("PPEL", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_PPEL.sep <- cSplit(fi3000comCGdata2_line2st_PPEL, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_PPEL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(39558954, 39559868), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_PPEL.sep.svg", width=25, height=5, units="cm", dpi=96)
#DIRAS3_TSS
fi3000comCGdata2_line2st_DIRAS3_TSS <- fi3000comCGdata2_line2st[grep("DIRAS3_TSS", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_DIRAS3_TSS.sep <- cSplit(fi3000comCGdata2_line2st_DIRAS3_TSS, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_DIRAS3_TSS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(68049750, 68051862), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_DIRAS3_TSS.sep.svg", width=25, height=5, units="cm", dpi=96)
#DIRAS3_ex2
fi3000comCGdata2_line2st_DIRAS3_ex2 <- fi3000comCGdata2_line2st[grep("DIRAS3_ex2", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_DIRAS3_ex2.sep <- cSplit(fi3000comCGdata2_line2st_DIRAS3_ex2, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_DIRAS3_ex2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(68046822, 68047803), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_DIRAS3_ex2.sep.svg", width=25, height=5, units="cm", dpi=96)
#GPR1_AS
fi3000comCGdata2_line2st_GPR1_AS <- fi3000comCGdata2_line2st[grep("GPR1-AS", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_GPR1_AS.sep <- cSplit(fi3000comCGdata2_line2st_GPR1_AS, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_GPR1_AS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(206202243, 206204721), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_GPR1_AS.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZDBF2_GPR1
fi3000comCGdata2_line2st_ZDBF2_GPR1 <- fi3000comCGdata2_line2st[grep("ZDBF2/GPR1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_ZDBF2_GPR1.sep <- cSplit(fi3000comCGdata2_line2st_ZDBF2_GPR1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_ZDBF2_GPR1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(206249859, 206271820), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_ZDBF2_GPR1.sep.svg", width=25, height=5, units="cm", dpi=96)
#NAP1L5
fi3000comCGdata2_line2st_NAP1L5 <- fi3000comCGdata2_line2st[grep("NAP1L5", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_NAP1L5.sep <- cSplit(fi3000comCGdata2_line2st_NAP1L5, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_NAP1L5.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(88697033, 88698086), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_NAP1L5.sep.svg", width=25, height=5, units="cm", dpi=96)
#VTRNA2_1
fi3000comCGdata2_line2st_VTRNA2_1 <- fi3000comCGdata2_line2st[grep("VTRNA2-1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_VTRNA2_1.sep <- cSplit(fi3000comCGdata2_line2st_VTRNA2_1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_VTRNA2_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(136079113, 136080956), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_VTRNA2_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#FAM50B
fi3000comCGdata2_line2st_FAM50B <- fi3000comCGdata2_line2st[grep("FAM50B", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_FAM50B.sep <- cSplit(fi3000comCGdata2_line2st_FAM50B, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_FAM50B.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3848848, 3850125), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_FAM50B.sep.svg", width=25, height=5, units="cm", dpi=96)
#PLAGL1
fi3000comCGdata2_line2st_PLAGL1 <- fi3000comCGdata2_line2st[grep("PLAGL1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_PLAGL1.sep <- cSplit(fi3000comCGdata2_line2st_PLAGL1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_PLAGL1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(144006941, 144008751), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_PLAGL1.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2R
fi3000comCGdata2_line2st_IGF2R <- fi3000comCGdata2_line2st[grep("IGF2R", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_IGF2R.sep <- cSplit(fi3000comCGdata2_line2st_IGF2R, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_IGF2R.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(160005526, 160006529), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_IGF2R.sep.svg", width=25, height=5, units="cm", dpi=96)
#VDR27
fi3000comCGdata2_line2st_VDR27 <- fi3000comCGdata2_line2st[grep("VDR27", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_VDR27.sep <- cSplit(fi3000comCGdata2_line2st_VDR27, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_VDR27.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(169654408, 169655522), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_VDR27.sep.svg", width=25, height=5, units="cm", dpi=96)
#GRB10
fi3000comCGdata2_line2st_GRB10 <- fi3000comCGdata2_line2st[grep("GRB10", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_GRB10.sep <- cSplit(fi3000comCGdata2_line2st_GRB10, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_GRB10.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(50781029, 50783615), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_GRB10.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG10
fi3000comCGdata2_line2st_PEG10 <- fi3000comCGdata2_line2st[grep("PEG10", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_PEG10.sep <- cSplit(fi3000comCGdata2_line2st_PEG10, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_PEG10.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(94656225, 94658648), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_PEG10.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEST
fi3000comCGdata2_line2st_MEST <- fi3000comCGdata2_line2st[grep("MEST", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_MEST.sep <- cSplit(fi3000comCGdata2_line2st_MEST, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_MEST.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(130490281, 130494547), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_MEST.sep.svg", width=25, height=5, units="cm", dpi=96)
#SVOPL
fi3000comCGdata2_line2st_SVOPL <- fi3000comCGdata2_line2st[grep("SVOPL", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_SVOPL.sep <- cSplit(fi3000comCGdata2_line2st_SVOPL, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_SVOPL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(138663373, 138664324), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_SVOPL.sep.svg", width=25, height=5, units="cm", dpi=96)
#HTR5A
fi3000comCGdata2_line2st_HTR5A <- fi3000comCGdata2_line2st[grep("HTR5A", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_HTR5A.sep <- cSplit(fi3000comCGdata2_line2st_HTR5A, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_HTR5A.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(155071009, 155071672), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_HTR5A.sep.svg", width=25, height=5, units="cm", dpi=96)
#ERILN2
fi3000comCGdata2_line2st_ERILN2 <- fi3000comCGdata2_line2st[grep("ERILN2", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_ERILN2.sep <- cSplit(fi3000comCGdata2_line2st_ERILN2, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_ERILN2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(37747474, 37748570), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_ERILN2.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG13
fi3000comCGdata2_line2st_PEG13 <- fi3000comCGdata2_line2st[grep("PEG13", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_PEG13.sep <- cSplit(fi3000comCGdata2_line2st_PEG13, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_PEG13.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(140098048, 140100982), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_PEG13.sep.svg", width=25, height=5, units="cm", dpi=96)
#FANCC
fi3000comCGdata2_line2st_FANCC <- fi3000comCGdata2_line2st[grep("FANCC", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_FANCC.sep <- cSplit(fi3000comCGdata2_line2st_FANCC, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_FANCC.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(95313118, 95313462), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_FANCC.sep.svg", width=25, height=5, units="cm", dpi=96)
#INPP5F
fi3000comCGdata2_line2st_INPP5F <- fi3000comCGdata2_line2st[grep("INPP5F", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_INPP5F.sep <- cSplit(fi3000comCGdata2_line2st_INPP5F, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_INPP5F.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(119818534, 119819215), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_INPP5F.sep.svg", width=25, height=5, units="cm", dpi=96)
#H19_IGF2
fi3000comCGdata2_line2st_H19_IGF2 <- fi3000comCGdata2_line2st[grep("H19/IGF2", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_H19_IGF2.sep <- cSplit(fi3000comCGdata2_line2st_H19_IGF2, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_H19_IGF2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(1997582, 2003510), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_H19_IGF2.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2_ex9
fi3000comCGdata2_line2st_IGF2_ex9 <- fi3000comCGdata2_line2st[grep("IGF2_ex9", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_IGF2_ex9.sep <- cSplit(fi3000comCGdata2_line2st_IGF2_ex9, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_IGF2_ex9.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2132761, 2133882), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_IGF2_ex9.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2_altTSS
fi3000comCGdata2_line2st_IGF2_altTSS <- fi3000comCGdata2_line2st[grep("IGF2_altTSS", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_IGF2_altTSS.sep <- cSplit(fi3000comCGdata2_line2st_IGF2_altTSS, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_IGF2_altTSS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2147103, 2148538), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_IGF2_altTSS.sep.svg", width=25, height=5, units="cm", dpi=96)
#KCNQ1OT1
fi3000comCGdata2_line2st_KCNQ1OT1 <- fi3000comCGdata2_line2st[grep("KCNQ1OT1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_KCNQ1OT1.sep <- cSplit(fi3000comCGdata2_line2st_KCNQ1OT1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_KCNQ1OT1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2698718, 2701029), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_KCNQ1OT1.sep.svg", width=25, height=5, units="cm", dpi=96)
#RB1
fi3000comCGdata2_line2st_RB1 <- fi3000comCGdata2_line2st[grep("RB1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_RB1.sep <- cSplit(fi3000comCGdata2_line2st_RB1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_RB1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(48318205, 48321627), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_RB1.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG3_DLK1
fi3000comCGdata2_line2st_MEG3_DLK1 <- fi3000comCGdata2_line2st[grep("MEG3/DLK1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_MEG3_DLK1.sep <- cSplit(fi3000comCGdata2_line2st_MEG3_DLK1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_MEG3_DLK1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100809090, 100811721), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_MEG3_DLK1.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG3
fi3000comCGdata2_line2st_MEG3 <- fi3000comCGdata2_line2st[grep("MEG3", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_MEG3.sep <- cSplit(fi3000comCGdata2_line2st_MEG3, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_MEG3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100824187, 100827641), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_MEG3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG8
fi3000comCGdata2_line2st_MEG8 <- fi3000comCGdata2_line2st[grep("MEG8", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_MEG8.sep <- cSplit(fi3000comCGdata2_line2st_MEG8, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_MEG8.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100904404, 100905082), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_MEG8.sep.svg", width=25, height=5, units="cm", dpi=96)
#MKRN3
fi3000comCGdata2_line2st_MKRN3 <- fi3000comCGdata2_line2st[grep("MKRN3", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_MKRN3.sep <- cSplit(fi3000comCGdata2_line2st_MKRN3, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_MKRN3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23561939, 23567348), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_MKRN3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MAGEL2
fi3000comCGdata2_line2st_MAGEL2 <- fi3000comCGdata2_line2st[grep("MAGEL2", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_MAGEL2.sep <- cSplit(fi3000comCGdata2_line2st_MAGEL2, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_MAGEL2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23647278, 23648882), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_MAGEL2.sep.svg", width=25, height=5, units="cm", dpi=96)
#NDN
fi3000comCGdata2_line2st_NDN <- fi3000comCGdata2_line2st[grep("NDN", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_NDN.sep <- cSplit(fi3000comCGdata2_line2st_NDN, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_NDN.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23686304, 23687612), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_NDN.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_alt
fi3000comCGdata2_line2st_SNRPN_alt <- fi3000comCGdata2_line2st[grep("SNRPN_alt", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_SNRPN_alt.sep <- cSplit(fi3000comCGdata2_line2st_SNRPN_alt, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_SNRPN_alt.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24823417, 24824334), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_SNRPN_alt.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_int1
fi3000comCGdata2_line2st_SNRPN_int1 <- fi3000comCGdata2_line2st[grep("SNRPN_int1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_SNRPN_int1.sep <- cSplit(fi3000comCGdata2_line2st_SNRPN_int1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_SNRPN_int1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24847861, 24848682), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_SNRPN_int1.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_int2
fi3000comCGdata2_line2st_SNRPN_int2 <- fi3000comCGdata2_line2st[grep("SNRPN_int2", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_SNRPN_int2.sep <- cSplit(fi3000comCGdata2_line2st_SNRPN_int2, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_SNRPN_int2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24877880, 24878758), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_SNRPN_int2.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNURF
fi3000comCGdata2_line2st_SNURF <- fi3000comCGdata2_line2st[grep("SNURF", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_SNURF.sep <- cSplit(fi3000comCGdata2_line2st_SNURF, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_SNURF.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24954857, 24956829), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_SNURF.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF1R
fi3000comCGdata2_line2st_IGF1R <- fi3000comCGdata2_line2st[grep("IGF1R", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_IGF1R.sep <- cSplit(fi3000comCGdata2_line2st_IGF1R, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_IGF1R.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(98865267, 98866421), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_IGF1R.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF597_1
fi3000comCGdata2_line2st_ZNF597_1 <- fi3000comCGdata2_line2st[grep("ZNF597_1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_ZNF597_1.sep <- cSplit(fi3000comCGdata2_line2st_ZNF597_1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_ZNF597_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3431801, 3432388), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_ZNF597_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF597_tss
fi3000comCGdata2_line2st_ZNF597_tss <- fi3000comCGdata2_line2st[grep("ZNF597_tss", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_ZNF597_tss.sep <- cSplit(fi3000comCGdata2_line2st_ZNF597_tss, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_ZNF597_tss.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3442828, 3444463), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_ZNF597_tss.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF331_1
fi3000comCGdata2_line2st_ZNF331_1 <- fi3000comCGdata2_line2st[grep("ZNF331_1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_ZNF331_1.sep <- cSplit(fi3000comCGdata2_line2st_ZNF331_1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_ZNF331_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(53537256, 53538958), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_ZNF331_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF331_2
fi3000comCGdata2_line2st_ZNF331_2 <- fi3000comCGdata2_line2st[grep("ZNF331_2", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_ZNF331_2.sep <- cSplit(fi3000comCGdata2_line2st_ZNF331_2, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_ZNF331_2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(53553832, 53555171), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_ZNF331_2.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG3
fi3000comCGdata2_line2st_PEG3 <- fi3000comCGdata2_line2st[grep("PEG3", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_PEG3.sep <- cSplit(fi3000comCGdata2_line2st_PEG3, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_PEG3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(56837125, 56841903), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_PEG3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MCTS2P
fi3000comCGdata2_line2st_MCTS2P <- fi3000comCGdata2_line2st[grep("MCTS2P", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_MCTS2P.sep <- cSplit(fi3000comCGdata2_line2st_MCTS2P, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_MCTS2P.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(31546860, 31548130), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_MCTS2P.sep.svg", width=25, height=5, units="cm", dpi=96)
#NNAT
fi3000comCGdata2_line2st_NNAT <- fi3000comCGdata2_line2st[grep("NNAT", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_NNAT.sep <- cSplit(fi3000comCGdata2_line2st_NNAT, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_NNAT.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(37520202, 37522126), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_NNAT.sep.svg", width=25, height=5, units="cm", dpi=96)
#L3MBTL1
fi3000comCGdata2_line2st_L3MBTL1 <- fi3000comCGdata2_line2st[grep("L3MBTL1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_L3MBTL1.sep <- cSplit(fi3000comCGdata2_line2st_L3MBTL1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_L3MBTL1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(43513725, 43515400), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_L3MBTL1.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_NESP
fi3000comCGdata2_line2st_GNAS_NESP <- fi3000comCGdata2_line2st[grep("GNAS-NESP", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_GNAS_NESP.sep <- cSplit(fi3000comCGdata2_line2st_GNAS_NESP, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_GNAS_NESP.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58838984, 58843557), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_GNAS_NESP.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_AS1
fi3000comCGdata2_line2st_GNAS_AS1 <- fi3000comCGdata2_line2st[grep("GNAS-AS1", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_GNAS_AS1.sep <- cSplit(fi3000comCGdata2_line2st_GNAS_AS1, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_GNAS_AS1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58850594, 58852978), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_GNAS_AS1.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_XL
fi3000comCGdata2_line2st_GNAS_XL <- fi3000comCGdata2_line2st[grep("GNAS-XL", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_GNAS_XL.sep <- cSplit(fi3000comCGdata2_line2st_GNAS_XL, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_GNAS_XL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58853850, 58856408), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_GNAS_XL.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNASA_B
fi3000comCGdata2_line2st_GNASA_B <- fi3000comCGdata2_line2st[grep("GNASA/B", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_GNASA_B.sep <- cSplit(fi3000comCGdata2_line2st_GNASA_B, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_GNASA_B.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58888210, 58890146), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_GNASA_B.sep.svg", width=25, height=5, units="cm", dpi=96)
#WRB
fi3000comCGdata2_line2st_WRB <- fi3000comCGdata2_line2st[grep("WRB", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_WRB.sep <- cSplit(fi3000comCGdata2_line2st_WRB, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_WRB.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(39385584, 39386350), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_WRB.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNU13
fi3000comCGdata2_line2st_SNU13 <- fi3000comCGdata2_line2st[grep("SNU13", fi3000comCGdata2_line2st$Csites),]
fi3000comCGdata2_line2st_SNU13.sep <- cSplit(fi3000comCGdata2_line2st_SNU13, "Csites", "%")
ggplot(fi3000comCGdata2_line2st_SNU13.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(41681770, 41682869), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi3000comCGdata2_line2st_SNU13.sep.svg", width=25, height=5, units="cm", dpi=96)




############################## Flank Line chart ############################

awk '{print $1"\t"$2-5000"\t"$3+5000"\t"$4"\t"$5}' /home/ankitv/ref_av/hg38_DMRs/Monk_human_ICR_hg38.bed  > Monk_human_ICR_hg38_flankplusminus5000.bed

bedtools intersect -wa -wb -a MERGE_myiCombat3reavg_pos_chr.txt -b Monk_human_ICR_hg38_flankplusminus5000.bed > MERGE_myicomCGbatCG3re_human_ICR_flanked5000.txt
MERGE_myicomCGbatCG3re_human_ICR_flanked5000 <- read.table("MERGE_myicomCGbatCG3re_human_ICR_flanked5000.txt", header = F)
head(MERGE_myicomCGbatCG3re_human_ICR_flanked5000)
dim(MERGE_myicomCGbatCG3re_human_ICR_flanked5000)
head(MERGE_myicomCGbatCG3re_human_ICR_flanked5000)
colnames(MERGE_myicomCGbatCG3re_human_ICR_flanked5000) <- c("chrc","startc","endc","id", colnames(MERGE_myiCombat3reavg_pos_chr[,5:22]),"chr","start","end","DMR", "DMRType")
fi5000comCGdata1 <- MERGE_myicomCGbatCG3re_human_ICR_flanked5000[,c((ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked5000)-1),ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked5000), 1:(ncol(MERGE_myicomCGbatCG3re_human_ICR_flanked5000)-5))]
rownames(fi5000comCGdata1)
colnames(fi5000comCGdata1)
head(fi5000comCGdata1)
dim(fi5000comCGdata1)
fi5000comCGdata2 <- fi5000comCGdata1
head(fi5000comCGdata2)
dim(fi5000comCGdata2)
fi5000comCGdata2_line <- data.frame(fi5000comCGdata2[,-c(2)])
head(fi5000comCGdata2_line)
fi5000comCGdata2_line["Row"] <- paste0(fi5000comCGdata2_line$DMR,"%",fi5000comCGdata2_line$chrc,"%",fi5000comCGdata2_line$startc,"%",fi5000comCGdata2_line$endc,"%",fi5000comCGdata2_line$id)
fi5000comCGdata2_line2 <- fi5000comCGdata2_line
rownames(fi5000comCGdata2_line2) <- fi5000comCGdata2_line2$Row
head(fi5000comCGdata2_line2)
fi5000comCGdata2_line2 <- fi5000comCGdata2_line2[,-c(1:6,24)]
head(fi5000comCGdata2_line2)
fi5000comCGdata2_line2 <- as.matrix(fi5000comCGdata2_line2)
fi5000comCGdata2_line2st <- stack(fi5000comCGdata2_line2)
head(fi5000comCGdata2_line2st)
fi5000comCGdata2_line2st <- data.frame(fi5000comCGdata2_line2st)
fi5000comCGdata2_line2st <- fi5000comCGdata2_line2st[,c(2,1,3)]
colnames(fi5000comCGdata2_line2st) <- c("Group", "Csites", "Value")
head(fi5000comCGdata2_line2st)


#PPEL
fi5000comCGdata2_line2st_PPEL <- fi5000comCGdata2_line2st[grep("PPEL", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_PPEL.sep <- cSplit(fi5000comCGdata2_line2st_PPEL, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_PPEL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(39558954, 39559868), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_PPEL.sep.svg", width=25, height=5, units="cm", dpi=96)
#DIRAS3_TSS
fi5000comCGdata2_line2st_DIRAS3_TSS <- fi5000comCGdata2_line2st[grep("DIRAS3_TSS", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_DIRAS3_TSS.sep <- cSplit(fi5000comCGdata2_line2st_DIRAS3_TSS, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_DIRAS3_TSS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(68049750, 68051862), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_DIRAS3_TSS.sep.svg", width=25, height=5, units="cm", dpi=96)
#DIRAS3_ex2
fi5000comCGdata2_line2st_DIRAS3_ex2 <- fi5000comCGdata2_line2st[grep("DIRAS3_ex2", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_DIRAS3_ex2.sep <- cSplit(fi5000comCGdata2_line2st_DIRAS3_ex2, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_DIRAS3_ex2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(68046822, 68047803), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_DIRAS3_ex2.sep.svg", width=25, height=5, units="cm", dpi=96)
#GPR1_AS
fi5000comCGdata2_line2st_GPR1_AS <- fi5000comCGdata2_line2st[grep("GPR1-AS", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_GPR1_AS.sep <- cSplit(fi5000comCGdata2_line2st_GPR1_AS, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_GPR1_AS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(206202243, 206204721), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_GPR1_AS.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZDBF2_GPR1
fi5000comCGdata2_line2st_ZDBF2_GPR1 <- fi5000comCGdata2_line2st[grep("ZDBF2/GPR1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_ZDBF2_GPR1.sep <- cSplit(fi5000comCGdata2_line2st_ZDBF2_GPR1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_ZDBF2_GPR1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(206249859, 206271820), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_ZDBF2_GPR1.sep.svg", width=25, height=5, units="cm", dpi=96)
#NAP1L5
fi5000comCGdata2_line2st_NAP1L5 <- fi5000comCGdata2_line2st[grep("NAP1L5", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_NAP1L5.sep <- cSplit(fi5000comCGdata2_line2st_NAP1L5, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_NAP1L5.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(88697033, 88698086), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_NAP1L5.sep.svg", width=25, height=5, units="cm", dpi=96)
#VTRNA2_1
fi5000comCGdata2_line2st_VTRNA2_1 <- fi5000comCGdata2_line2st[grep("VTRNA2-1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_VTRNA2_1.sep <- cSplit(fi5000comCGdata2_line2st_VTRNA2_1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_VTRNA2_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(136079113, 136080956), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_VTRNA2_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#FAM50B
fi5000comCGdata2_line2st_FAM50B <- fi5000comCGdata2_line2st[grep("FAM50B", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_FAM50B.sep <- cSplit(fi5000comCGdata2_line2st_FAM50B, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_FAM50B.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3848848, 3850125), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_FAM50B.sep.svg", width=25, height=5, units="cm", dpi=96)
#PLAGL1
fi5000comCGdata2_line2st_PLAGL1 <- fi5000comCGdata2_line2st[grep("PLAGL1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_PLAGL1.sep <- cSplit(fi5000comCGdata2_line2st_PLAGL1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_PLAGL1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(144006941, 144008751), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_PLAGL1.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2R
fi5000comCGdata2_line2st_IGF2R <- fi5000comCGdata2_line2st[grep("IGF2R", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_IGF2R.sep <- cSplit(fi5000comCGdata2_line2st_IGF2R, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_IGF2R.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(160005526, 160006529), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_IGF2R.sep.svg", width=25, height=5, units="cm", dpi=96)
#VDR27
fi5000comCGdata2_line2st_VDR27 <- fi5000comCGdata2_line2st[grep("VDR27", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_VDR27.sep <- cSplit(fi5000comCGdata2_line2st_VDR27, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_VDR27.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(169654408, 169655522), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_VDR27.sep.svg", width=25, height=5, units="cm", dpi=96)
#GRB10
fi5000comCGdata2_line2st_GRB10 <- fi5000comCGdata2_line2st[grep("GRB10", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_GRB10.sep <- cSplit(fi5000comCGdata2_line2st_GRB10, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_GRB10.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(50781029, 50783615), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_GRB10.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG10
fi5000comCGdata2_line2st_PEG10 <- fi5000comCGdata2_line2st[grep("PEG10", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_PEG10.sep <- cSplit(fi5000comCGdata2_line2st_PEG10, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_PEG10.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(94656225, 94658648), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_PEG10.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEST
fi5000comCGdata2_line2st_MEST <- fi5000comCGdata2_line2st[grep("MEST", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_MEST.sep <- cSplit(fi5000comCGdata2_line2st_MEST, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_MEST.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(130490281, 130494547), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_MEST.sep.svg", width=25, height=5, units="cm", dpi=96)
#SVOPL
fi5000comCGdata2_line2st_SVOPL <- fi5000comCGdata2_line2st[grep("SVOPL", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_SVOPL.sep <- cSplit(fi5000comCGdata2_line2st_SVOPL, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_SVOPL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(138663373, 138664324), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_SVOPL.sep.svg", width=25, height=5, units="cm", dpi=96)
#HTR5A
fi5000comCGdata2_line2st_HTR5A <- fi5000comCGdata2_line2st[grep("HTR5A", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_HTR5A.sep <- cSplit(fi5000comCGdata2_line2st_HTR5A, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_HTR5A.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(155071009, 155071672), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_HTR5A.sep.svg", width=25, height=5, units="cm", dpi=96)
#ERILN2
fi5000comCGdata2_line2st_ERILN2 <- fi5000comCGdata2_line2st[grep("ERILN2", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_ERILN2.sep <- cSplit(fi5000comCGdata2_line2st_ERILN2, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_ERILN2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(37747474, 37748570), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_ERILN2.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG13
fi5000comCGdata2_line2st_PEG13 <- fi5000comCGdata2_line2st[grep("PEG13", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_PEG13.sep <- cSplit(fi5000comCGdata2_line2st_PEG13, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_PEG13.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(140098048, 140100982), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_PEG13.sep.svg", width=25, height=5, units="cm", dpi=96)
#FANCC
fi5000comCGdata2_line2st_FANCC <- fi5000comCGdata2_line2st[grep("FANCC", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_FANCC.sep <- cSplit(fi5000comCGdata2_line2st_FANCC, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_FANCC.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(95313118, 95313462), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_FANCC.sep.svg", width=25, height=5, units="cm", dpi=96)
#INPP5F
fi5000comCGdata2_line2st_INPP5F <- fi5000comCGdata2_line2st[grep("INPP5F", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_INPP5F.sep <- cSplit(fi5000comCGdata2_line2st_INPP5F, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_INPP5F.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(119818534, 119819215), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_INPP5F.sep.svg", width=25, height=5, units="cm", dpi=96)
#H19_IGF2
fi5000comCGdata2_line2st_H19_IGF2 <- fi5000comCGdata2_line2st[grep("H19/IGF2", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_H19_IGF2.sep <- cSplit(fi5000comCGdata2_line2st_H19_IGF2, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_H19_IGF2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(1997582, 2003510), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_H19_IGF2.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2_ex9
fi5000comCGdata2_line2st_IGF2_ex9 <- fi5000comCGdata2_line2st[grep("IGF2_ex9", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_IGF2_ex9.sep <- cSplit(fi5000comCGdata2_line2st_IGF2_ex9, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_IGF2_ex9.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2132761, 2133882), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_IGF2_ex9.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF2_altTSS
fi5000comCGdata2_line2st_IGF2_altTSS <- fi5000comCGdata2_line2st[grep("IGF2_altTSS", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_IGF2_altTSS.sep <- cSplit(fi5000comCGdata2_line2st_IGF2_altTSS, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_IGF2_altTSS.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2147103, 2148538), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_IGF2_altTSS.sep.svg", width=25, height=5, units="cm", dpi=96)
#KCNQ1OT1
fi5000comCGdata2_line2st_KCNQ1OT1 <- fi5000comCGdata2_line2st[grep("KCNQ1OT1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_KCNQ1OT1.sep <- cSplit(fi5000comCGdata2_line2st_KCNQ1OT1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_KCNQ1OT1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(2698718, 2701029), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_KCNQ1OT1.sep.svg", width=25, height=5, units="cm", dpi=96)
#RB1
fi5000comCGdata2_line2st_RB1 <- fi5000comCGdata2_line2st[grep("RB1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_RB1.sep <- cSplit(fi5000comCGdata2_line2st_RB1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_RB1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(48318205, 48321627), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_RB1.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG3_DLK1
fi5000comCGdata2_line2st_MEG3_DLK1 <- fi5000comCGdata2_line2st[grep("MEG3/DLK1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_MEG3_DLK1.sep <- cSplit(fi5000comCGdata2_line2st_MEG3_DLK1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_MEG3_DLK1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100809090, 100811721), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_MEG3_DLK1.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG3
fi5000comCGdata2_line2st_MEG3 <- fi5000comCGdata2_line2st[grep("MEG3", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_MEG3.sep <- cSplit(fi5000comCGdata2_line2st_MEG3, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_MEG3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100824187, 100827641), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_MEG3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MEG8
fi5000comCGdata2_line2st_MEG8 <- fi5000comCGdata2_line2st[grep("MEG8", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_MEG8.sep <- cSplit(fi5000comCGdata2_line2st_MEG8, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_MEG8.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(100904404, 100905082), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_MEG8.sep.svg", width=25, height=5, units="cm", dpi=96)
#MKRN3
fi5000comCGdata2_line2st_MKRN3 <- fi5000comCGdata2_line2st[grep("MKRN3", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_MKRN3.sep <- cSplit(fi5000comCGdata2_line2st_MKRN3, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_MKRN3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23561939, 23567348), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_MKRN3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MAGEL2
fi5000comCGdata2_line2st_MAGEL2 <- fi5000comCGdata2_line2st[grep("MAGEL2", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_MAGEL2.sep <- cSplit(fi5000comCGdata2_line2st_MAGEL2, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_MAGEL2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23647278, 23648882), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_MAGEL2.sep.svg", width=25, height=5, units="cm", dpi=96)
#NDN
fi5000comCGdata2_line2st_NDN <- fi5000comCGdata2_line2st[grep("NDN", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_NDN.sep <- cSplit(fi5000comCGdata2_line2st_NDN, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_NDN.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(23686304, 23687612), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_NDN.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_alt
fi5000comCGdata2_line2st_SNRPN_alt <- fi5000comCGdata2_line2st[grep("SNRPN_alt", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_SNRPN_alt.sep <- cSplit(fi5000comCGdata2_line2st_SNRPN_alt, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_SNRPN_alt.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24823417, 24824334), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_SNRPN_alt.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_int1
fi5000comCGdata2_line2st_SNRPN_int1 <- fi5000comCGdata2_line2st[grep("SNRPN_int1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_SNRPN_int1.sep <- cSplit(fi5000comCGdata2_line2st_SNRPN_int1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_SNRPN_int1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24847861, 24848682), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_SNRPN_int1.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNRPN_int2
fi5000comCGdata2_line2st_SNRPN_int2 <- fi5000comCGdata2_line2st[grep("SNRPN_int2", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_SNRPN_int2.sep <- cSplit(fi5000comCGdata2_line2st_SNRPN_int2, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_SNRPN_int2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24877880, 24878758), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_SNRPN_int2.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNURF
fi5000comCGdata2_line2st_SNURF <- fi5000comCGdata2_line2st[grep("SNURF", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_SNURF.sep <- cSplit(fi5000comCGdata2_line2st_SNURF, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_SNURF.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(24954857, 24956829), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_SNURF.sep.svg", width=25, height=5, units="cm", dpi=96)
#IGF1R
fi5000comCGdata2_line2st_IGF1R <- fi5000comCGdata2_line2st[grep("IGF1R", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_IGF1R.sep <- cSplit(fi5000comCGdata2_line2st_IGF1R, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_IGF1R.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(98865267, 98866421), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_IGF1R.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF597_1
fi5000comCGdata2_line2st_ZNF597_1 <- fi5000comCGdata2_line2st[grep("ZNF597_1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_ZNF597_1.sep <- cSplit(fi5000comCGdata2_line2st_ZNF597_1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_ZNF597_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3431801, 3432388), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_ZNF597_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF597_tss
fi5000comCGdata2_line2st_ZNF597_tss <- fi5000comCGdata2_line2st[grep("ZNF597_tss", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_ZNF597_tss.sep <- cSplit(fi5000comCGdata2_line2st_ZNF597_tss, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_ZNF597_tss.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(3442828, 3444463), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_ZNF597_tss.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF331_1
fi5000comCGdata2_line2st_ZNF331_1 <- fi5000comCGdata2_line2st[grep("ZNF331_1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_ZNF331_1.sep <- cSplit(fi5000comCGdata2_line2st_ZNF331_1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_ZNF331_1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(53537256, 53538958), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_ZNF331_1.sep.svg", width=25, height=5, units="cm", dpi=96)
#ZNF331_2
fi5000comCGdata2_line2st_ZNF331_2 <- fi5000comCGdata2_line2st[grep("ZNF331_2", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_ZNF331_2.sep <- cSplit(fi5000comCGdata2_line2st_ZNF331_2, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_ZNF331_2.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(53553832, 53555171), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_ZNF331_2.sep.svg", width=25, height=5, units="cm", dpi=96)
#PEG3
fi5000comCGdata2_line2st_PEG3 <- fi5000comCGdata2_line2st[grep("PEG3", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_PEG3.sep <- cSplit(fi5000comCGdata2_line2st_PEG3, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_PEG3.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(56837125, 56841903), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_PEG3.sep.svg", width=25, height=5, units="cm", dpi=96)
#MCTS2P
fi5000comCGdata2_line2st_MCTS2P <- fi5000comCGdata2_line2st[grep("MCTS2P", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_MCTS2P.sep <- cSplit(fi5000comCGdata2_line2st_MCTS2P, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_MCTS2P.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(31546860, 31548130), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_MCTS2P.sep.svg", width=25, height=5, units="cm", dpi=96)
#NNAT
fi5000comCGdata2_line2st_NNAT <- fi5000comCGdata2_line2st[grep("NNAT", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_NNAT.sep <- cSplit(fi5000comCGdata2_line2st_NNAT, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_NNAT.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(37520202, 37522126), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_NNAT.sep.svg", width=25, height=5, units="cm", dpi=96)
#L3MBTL1
fi5000comCGdata2_line2st_L3MBTL1 <- fi5000comCGdata2_line2st[grep("L3MBTL1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_L3MBTL1.sep <- cSplit(fi5000comCGdata2_line2st_L3MBTL1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_L3MBTL1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(43513725, 43515400), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_L3MBTL1.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_NESP
fi5000comCGdata2_line2st_GNAS_NESP <- fi5000comCGdata2_line2st[grep("GNAS-NESP", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_GNAS_NESP.sep <- cSplit(fi5000comCGdata2_line2st_GNAS_NESP, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_GNAS_NESP.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58838984, 58843557), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_GNAS_NESP.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_AS1
fi5000comCGdata2_line2st_GNAS_AS1 <- fi5000comCGdata2_line2st[grep("GNAS-AS1", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_GNAS_AS1.sep <- cSplit(fi5000comCGdata2_line2st_GNAS_AS1, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_GNAS_AS1.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58850594, 58852978), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_GNAS_AS1.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNAS_XL
fi5000comCGdata2_line2st_GNAS_XL <- fi5000comCGdata2_line2st[grep("GNAS-XL", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_GNAS_XL.sep <- cSplit(fi5000comCGdata2_line2st_GNAS_XL, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_GNAS_XL.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58853850, 58856408), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_GNAS_XL.sep.svg", width=25, height=5, units="cm", dpi=96)
#GNASA_B
fi5000comCGdata2_line2st_GNASA_B <- fi5000comCGdata2_line2st[grep("GNASA/B", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_GNASA_B.sep <- cSplit(fi5000comCGdata2_line2st_GNASA_B, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_GNASA_B.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(58888210, 58890146), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_GNASA_B.sep.svg", width=25, height=5, units="cm", dpi=96)
#WRB
fi5000comCGdata2_line2st_WRB <- fi5000comCGdata2_line2st[grep("WRB", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_WRB.sep <- cSplit(fi5000comCGdata2_line2st_WRB, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_WRB.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(39385584, 39386350), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_WRB.sep.svg", width=25, height=5, units="cm", dpi=96)
#SNU13
fi5000comCGdata2_line2st_SNU13 <- fi5000comCGdata2_line2st[grep("SNU13", fi5000comCGdata2_line2st$Csites),]
fi5000comCGdata2_line2st_SNU13.sep <- cSplit(fi5000comCGdata2_line2st_SNU13, "Csites", "%")
ggplot(fi5000comCGdata2_line2st_SNU13.sep, aes(x=as.numeric(Csites_3), y=Value, group=Group))+
  geom_line(aes(color=Group))+
  geom_vline(xintercept = c(41681770, 41682869), colour = "black", linetype="solid")+
  geom_point(aes(color=Group), shape=20) + ylim(c(0,1))+
  scale_color_manual(values=c("darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey","#EC7063","#922B21","#148F77","#48C9B0","#1E8449","#52BE80"))+
  theme_classic()
ggsave("Linechart_fi5000comCGdata2_line2st_SNU13.sep.svg", width=25, height=5, units="cm", dpi=96)


WGBS_int <- c("HTR5A","IGF2_altTSS","IGF2_ex9","IGF2R","INPP5F","MKRN3","FAM50B","ZNF331_2","DIRAS3_TSS","ERILN2","ZNF597_3DMR","H19/IGF2","MAGEL2","NDN","MEG3_1")
Array_int <- c("ERILN2","MKRN3","ZNF331_2","H19/IGF2","DIRAS3_ex2","IGF2_ex9","HTR5A","MEG3_1","L3MBTL1_1","ZNF597_3DMR","MAGEL2")
library(VennDiagram)
venn.diagram(
  x = list(WGBS_int, Array_int),
  category.names = c("WGBS_int" , "Array_int"),
  filename = 'wgbs_array_venn_diagramm.png',
  output=TRUE
)

head(MERGE_myiCombat3reavg_human_ICR)
library(ggstatsplot)
#DIRAS3_ex2
MERGE_myiCombat3reavg_human_ICR_DIRAS3_ex2 <- MERGE_myiCombat3reavg_human_ICR[grep("DIRAS3_ex2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_DIRAS3_ex2$Allcontrol, MERGE_myiCombat3reavg_human_ICR_DIRAS3_ex2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_DIRAS3_ex2$PG, MERGE_myiCombat3reavg_human_ICR_DIRAS3_ex2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_DIRAS3_ex2$PG, MERGE_myiCombat3reavg_human_ICR_DIRAS3_ex2$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_DIRAS3_ex2[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_DIRAS3_ex2.svg", width=15, height=15, units="cm", dpi=96)


#ERILN2
MERGE_myiCombat3reavg_human_ICR_ERILN2 <- MERGE_myiCombat3reavg_human_ICR[grep("ERILN2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_ERILN2$Allcontrol, MERGE_myiCombat3reavg_human_ICR_ERILN2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_ERILN2$PG, MERGE_myiCombat3reavg_human_ICR_ERILN2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_ERILN2$PG, MERGE_myiCombat3reavg_human_ICR_ERILN2$C50, paired = T, alternative = "two.sided")$p.value


ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_ERILN2[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_ERILN2.svg", width=15, height=15, units="cm", dpi=96)

#H19_IGF2
MERGE_myiCombat3reavg_human_ICR_H19_IGF2 <- MERGE_myiCombat3reavg_human_ICR[grep("H19/IGF2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_H19_IGF2$Allcontrol, MERGE_myiCombat3reavg_human_ICR_H19_IGF2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_H19_IGF2$PG, MERGE_myiCombat3reavg_human_ICR_H19_IGF2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_H19_IGF2$PG, MERGE_myiCombat3reavg_human_ICR_H19_IGF2$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_H19_IGF2[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_H19_IGF2.svg", width=15, height=15, units="cm", dpi=96)

#HTR5A
MERGE_myiCombat3reavg_human_ICR_HTR5A <- MERGE_myiCombat3reavg_human_ICR[grep("HTR5A", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_HTR5A$Allcontrol, MERGE_myiCombat3reavg_human_ICR_HTR5A$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_HTR5A$PG, MERGE_myiCombat3reavg_human_ICR_HTR5A$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_HTR5A$PG, MERGE_myiCombat3reavg_human_ICR_HTR5A$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_HTR5A[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_HTR5A.svg", width=15, height=15, units="cm", dpi=96)


#IGF2_ex9
MERGE_myiCombat3reavg_human_ICR_IGF2_ex9 <- MERGE_myiCombat3reavg_human_ICR[grep("IGF2_ex9", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_IGF2_ex9$Allcontrol, MERGE_myiCombat3reavg_human_ICR_IGF2_ex9$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_IGF2_ex9$PG, MERGE_myiCombat3reavg_human_ICR_IGF2_ex9$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_IGF2_ex9$PG, MERGE_myiCombat3reavg_human_ICR_IGF2_ex9$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_HTR5A[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_HTR5A.svg", width=15, height=15, units="cm", dpi=96)


#L3MBTL1
MERGE_myiCombat3reavg_human_ICR_L3MBTL1 <- MERGE_myiCombat3reavg_human_ICR[grep("L3MBTL1", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_L3MBTL1$Allcontrol, MERGE_myiCombat3reavg_human_ICR_L3MBTL1$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_L3MBTL1$PG, MERGE_myiCombat3reavg_human_ICR_L3MBTL1$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_L3MBTL1$PG, MERGE_myiCombat3reavg_human_ICR_L3MBTL1$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_L3MBTL1[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_L3MBTL1.svg", width=15, height=15, units="cm", dpi=96)

#MAGEL2
MERGE_myiCombat3reavg_human_ICR_MAGEL2 <- MERGE_myiCombat3reavg_human_ICR[grep("MAGEL2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MAGEL2$Allcontrol, MERGE_myiCombat3reavg_human_ICR_MAGEL2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MAGEL2$PG, MERGE_myiCombat3reavg_human_ICR_MAGEL2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MAGEL2$PG, MERGE_myiCombat3reavg_human_ICR_MAGEL2$C50, paired = T, alternative = "two.sided")$p.value
ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_MAGEL2[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_MAGEL2.svg", width=15, height=15, units="cm", dpi=96)

#MEG3
MERGE_myiCombat3reavg_human_ICR_MEG3 <- MERGE_myiCombat3reavg_human_ICR[grep("MEG3", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MEG3$Allcontrol, MERGE_myiCombat3reavg_human_ICR_MEG3$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MEG3$PG, MERGE_myiCombat3reavg_human_ICR_MEG3$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MEG3$PG, MERGE_myiCombat3reavg_human_ICR_MEG3$C50, paired = T, alternative = "two.sided")$p.value
ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_MEG3[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_MEG3.svg", width=15, height=15, units="cm", dpi=96)

#MKRN3
MERGE_myiCombat3reavg_human_ICR_MKRN3 <- MERGE_myiCombat3reavg_human_ICR[grep("MKRN3", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MKRN3$Allcontrol, MERGE_myiCombat3reavg_human_ICR_MKRN3$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MKRN3$PG, MERGE_myiCombat3reavg_human_ICR_MKRN3$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_MKRN3$PG, MERGE_myiCombat3reavg_human_ICR_MKRN3$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_MKRN3[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_MKRN3.svg", width=15, height=15, units="cm", dpi=96)


#NNAT
MERGE_myiCombat3reavg_human_ICR_NNAT <- MERGE_myiCombat3reavg_human_ICR[grep("NNAT", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_NNAT$Allcontrol, MERGE_myiCombat3reavg_human_ICR_NNAT$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_NNAT$PG, MERGE_myiCombat3reavg_human_ICR_NNAT$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_NNAT$PG, MERGE_myiCombat3reavg_human_ICR_NNAT$C50, paired = T, alternative = "two.sided")$p.value
ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_NNAT[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_NNAT.svg", width=15, height=15, units="cm", dpi=96)

#ZNF331_2
MERGE_myiCombat3reavg_human_ICR_ZNF331_2 <- MERGE_myiCombat3reavg_human_ICR[grep("ZNF331_2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_ZNF331_2$Allcontrol, MERGE_myiCombat3reavg_human_ICR_ZNF331_2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_ZNF331_2$PG, MERGE_myiCombat3reavg_human_ICR_ZNF331_2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_ZNF331_2$PG, MERGE_myiCombat3reavg_human_ICR_ZNF331_2$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_ZNF331_2[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_ZNF331_2.svg", width=15, height=15, units="cm", dpi=96)

#NNAT
MERGE_myiCombat3reavg_human_ICR_NNAT <- MERGE_myiCombat3reavg_human_ICR[grep("NNAT", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(MERGE_myiCombat3reavg_human_ICR_NNAT$Allcontrol, MERGE_myiCombat3reavg_human_ICR_NNAT$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_NNAT$PG, MERGE_myiCombat3reavg_human_ICR_NNAT$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(MERGE_myiCombat3reavg_human_ICR_NNAT$PG, MERGE_myiCombat3reavg_human_ICR_NNAT$C50, paired = T, alternative = "two.sided")$p.value
ggbetweenstats(data.frame(stack(as.matrix(MERGE_myiCombat3reavg_human_ICR_NNAT[,c(5,17,20,22)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PG_myiCombat3reavg_human_ICR_NNAT.svg", width=15, height=15, units="cm", dpi=96)


library(ggstatsplot)
#DIRAS3_ex2
PR_myiCombat3reavg_human_ICR_DIRAS3_ex2 <- MERGE_myiCombat3reavg_human_ICR[grep("DIRAS3_ex2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_DIRAS3_ex2$Allcontrol, PR_myiCombat3reavg_human_ICR_DIRAS3_ex2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_DIRAS3_ex2$PG, PR_myiCombat3reavg_human_ICR_DIRAS3_ex2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_DIRAS3_ex2$PG, PR_myiCombat3reavg_human_ICR_DIRAS3_ex2$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_DIRAS3_ex2[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_DIRAS3_ex2.svg", width=15, height=15, units="cm", dpi=96)


#ERILN2
PR_myiCombat3reavg_human_ICR_ERILN2 <- MERGE_myiCombat3reavg_human_ICR[grep("ERILN2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_ERILN2$Allcontrol, PR_myiCombat3reavg_human_ICR_ERILN2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_ERILN2$PG, PR_myiCombat3reavg_human_ICR_ERILN2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_ERILN2$PG, PR_myiCombat3reavg_human_ICR_ERILN2$C50, paired = T, alternative = "two.sided")$p.value


ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_ERILN2[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_ERILN2.svg", width=15, height=15, units="cm", dpi=96)

#H19_IGF2
PR_myiCombat3reavg_human_ICR_H19_IGF2 <- MERGE_myiCombat3reavg_human_ICR[grep("H19/IGF2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_H19_IGF2$Allcontrol, PR_myiCombat3reavg_human_ICR_H19_IGF2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_H19_IGF2$PG, PR_myiCombat3reavg_human_ICR_H19_IGF2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_H19_IGF2$PG, PR_myiCombat3reavg_human_ICR_H19_IGF2$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_H19_IGF2[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_H19_IGF2.svg", width=15, height=15, units="cm", dpi=96)

#HTR5A
PR_myiCombat3reavg_human_ICR_HTR5A <- MERGE_myiCombat3reavg_human_ICR[grep("HTR5A", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_HTR5A$Allcontrol, PR_myiCombat3reavg_human_ICR_HTR5A$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_HTR5A$PG, PR_myiCombat3reavg_human_ICR_HTR5A$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_HTR5A$PG, PR_myiCombat3reavg_human_ICR_HTR5A$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_HTR5A[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_HTR5A.svg", width=15, height=15, units="cm", dpi=96)


#IGF2_ex9
PR_myiCombat3reavg_human_ICR_IGF2_ex9 <- MERGE_myiCombat3reavg_human_ICR[grep("IGF2_ex9", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_IGF2_ex9$Allcontrol, PR_myiCombat3reavg_human_ICR_IGF2_ex9$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_IGF2_ex9$PG, PR_myiCombat3reavg_human_ICR_IGF2_ex9$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_IGF2_ex9$PG, PR_myiCombat3reavg_human_ICR_IGF2_ex9$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_HTR5A[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_HTR5A.svg", width=15, height=15, units="cm", dpi=96)

#MAGEL2
PR_myiCombat3reavg_human_ICR_MAGEL2 <- MERGE_myiCombat3reavg_human_ICR[grep("MAGEL2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_MAGEL2$Allcontrol, PR_myiCombat3reavg_human_ICR_MAGEL2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_MAGEL2$PG, PR_myiCombat3reavg_human_ICR_MAGEL2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_MAGEL2$PG, PR_myiCombat3reavg_human_ICR_MAGEL2$C50, paired = T, alternative = "two.sided")$p.value
ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_MAGEL2[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_MAGEL2.svg", width=15, height=15, units="cm", dpi=96)

#MKRN3
PR_myiCombat3reavg_human_ICR_MKRN3 <- MERGE_myiCombat3reavg_human_ICR[grep("MKRN3", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_MKRN3$Allcontrol, PR_myiCombat3reavg_human_ICR_MKRN3$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_MKRN3$PG, PR_myiCombat3reavg_human_ICR_MKRN3$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_MKRN3$PG, PR_myiCombat3reavg_human_ICR_MKRN3$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_MKRN3[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_MKRN3.svg", width=15, height=15, units="cm", dpi=96)


#SNRPN_alt
PR_myiCombat3reavg_human_ICR_SNRPN_alt <- MERGE_myiCombat3reavg_human_ICR[grep("SNRPN_alt", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_SNRPN_alt$Allcontrol, PR_myiCombat3reavg_human_ICR_SNRPN_alt$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_SNRPN_alt$PG, PR_myiCombat3reavg_human_ICR_SNRPN_alt$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_SNRPN_alt$PG, PR_myiCombat3reavg_human_ICR_SNRPN_alt$C50, paired = T, alternative = "two.sided")$p.value
ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_SNRPN_alt[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_SNRPN_alt.svg", width=15, height=15, units="cm", dpi=96)

#ZNF331_2
PR_myiCombat3reavg_human_ICR_ZNF331_2 <- MERGE_myiCombat3reavg_human_ICR[grep("ZNF331_2", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_ZNF331_2$Allcontrol, PR_myiCombat3reavg_human_ICR_ZNF331_2$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_ZNF331_2$PG, PR_myiCombat3reavg_human_ICR_ZNF331_2$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_ZNF331_2$PG, PR_myiCombat3reavg_human_ICR_ZNF331_2$C50, paired = T, alternative = "two.sided")$p.value

ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_ZNF331_2[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_ZNF331_2.svg", width=15, height=15, units="cm", dpi=96)

#SNRPN_alt
PR_myiCombat3reavg_human_ICR_SNRPN_alt <- MERGE_myiCombat3reavg_human_ICR[grep("SNRPN_alt", MERGE_myiCombat3reavg_human_ICR$DMR),]
wilcox.test(PR_myiCombat3reavg_human_ICR_SNRPN_alt$Allcontrol, PR_myiCombat3reavg_human_ICR_SNRPN_alt$PG, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_SNRPN_alt$PG, PR_myiCombat3reavg_human_ICR_SNRPN_alt$C13, paired = T, alternative = "two.sided")$p.value
wilcox.test(PR_myiCombat3reavg_human_ICR_SNRPN_alt$PG, PR_myiCombat3reavg_human_ICR_SNRPN_alt$C50, paired = T, alternative = "two.sided")$p.value
ggbetweenstats(data.frame(stack(as.matrix(PR_myiCombat3reavg_human_ICR_SNRPN_alt[,c(5,18,19,21)]))), 
               x=col, 
               y=value,
               type = "nonparametric",
               p.adjust.method="BH",
               pairwise.display = "all",
               pairwise.comparisons = TRUE)
ggsave("statistics_PR__myiCombat3reavg_human_ICR_SNRPN_alt.svg", width=15, height=15, units="cm", dpi=96)


#RNA-Seq ICFs
MM_PK_Imp_gene_Exp_FPKM <- read.table("/media/ankitv/Archivio2/ankit/Array21/ICF_RNASeq/MM_PK_Imp_gene_Exp_FPKM.txt", header = T, stringsAsFactors = F)
MM_PK_Imp_gene_Exp_FPKM["LogavgUN"] <- log2((MM_PK_Imp_gene_Exp_FPKM$Un2 + MM_PK_Imp_gene_Exp_FPKM$UN)/2)
MM_PK_Imp_gene_Exp_FPKM["LogpG"] <- log2(MM_PK_Imp_gene_Exp_FPKM$pG)
MM_PK_Imp_gene_Exp_FPKM["Logc13"] <- log2(MM_PK_Imp_gene_Exp_FPKM$c13)
MM_PK_Imp_gene_Exp_FPKM["Logc50"] <- log2(MM_PK_Imp_gene_Exp_FPKM$c50)
MM_PK_Imp_gene_Exp_FPKM["LogpR"] <- log2(MM_PK_Imp_gene_Exp_FPKM$pR)
MM_PK_Imp_gene_Exp_FPKM["Logc7"] <- log2(MM_PK_Imp_gene_Exp_FPKM$c7)
head(MM_PK_Imp_gene_Exp_FPKM)


MM_PK_Imp_gene_Exp_FPKM_MEG3 <- MM_PK_Imp_gene_Exp_FPKM[grep("MEG3", MM_PK_Imp_gene_Exp_FPKM$gene_name),]
MM_PK_Imp_gene_Exp_FPKM_MEG3Log <- data.frame(stack(t(MM_PK_Imp_gene_Exp_FPKM_MEG3[,10:15])))

ggplot(MM_PK_Imp_gene_Exp_FPKM_MEG3Log)+
  geom_bar(aes(x=row, y=value, col=row, fill=row),stat = "identity")+
  ggtitle("MEG3")+
  theme_bw()+ scale_color_manual(values=c("darkgrey","red","green","#009A17","darkred","darkgreen"))+ 
  scale_fill_manual(values=c("white","white","white","white","white","white"))
ggsave("MM_PK_Imp_gene_Exp_FPKM_MEG3Log.svg", width=15, height=15, units="cm", dpi=96)
barplot( height=MM_PK_Imp_gene_Exp_FPKM_MEG3Log$value, names=MM_PK_Imp_gene_Exp_FPKM_MEG3Log$row , density=c(15,15,15,15,15,15,15) , angle=c(45,45,45,45,45,45,45) , col=c("darkgrey","#EC7063","#48C9B0","#52BE80","#922B21","#148F77"),xlab = "Samples", ylab = "Log2(FPKM)",main = "MEG3",width=c(5,5,5,5,5,5))
#MM_PK_Imp_gene_Exp_FPKM_MEG3Log.png
MM_PK_Imp_gene_Exp_FPKM_H19 <- MM_PK_Imp_gene_Exp_FPKM[grep("H19", MM_PK_Imp_gene_Exp_FPKM$gene_name),]
MM_PK_Imp_gene_Exp_FPKM_H19Log <- data.frame(stack(t(MM_PK_Imp_gene_Exp_FPKM_H19[,10:15])))

ggplot(MM_PK_Imp_gene_Exp_FPKM_H19Log)+
  geom_bar(aes(x=row, y=value, col=row, fill=row),stat = "identity")+
  ggtitle("H19")+
  theme_bw()+ scale_color_manual(values=c("darkgrey","red","green","#009A17","darkred","darkgreen"))+ 
  scale_fill_manual(values=c("white","white","white","white","white","white"))

ggsave("MM_PK_Imp_gene_Exp_FPKM_H19Log.svg", width=15, height=15, units="cm", dpi=96)
barplot( height=MM_PK_Imp_gene_Exp_FPKM_H19Log$value, names=MM_PK_Imp_gene_Exp_FPKM_H19Log$row , density=c(12,12,12,12,12,12,12) , angle=c(45,45,45,45,45,45,45) , col=c("darkgrey","#EC7063","#48C9B0","#52BE80","#922B21","#148F77"),xlab = "Samples", ylab = "Log2(FPKM)",main = "H19")
#MM_PK_Imp_gene_Exp_FPKM_H19Log.png

MM_PK_Imp_gene_Exp_FPKM_MKRN3 <- MM_PK_Imp_gene_Exp_FPKM[grep("MKRN3", MM_PK_Imp_gene_Exp_FPKM$gene_name),]
MM_PK_Imp_gene_Exp_FPKM_MKRN3Log <- data.frame(stack(t(MM_PK_Imp_gene_Exp_FPKM_MKRN3[,c(10:15)])))

ggplot(MM_PK_Imp_gene_Exp_FPKM_MKRN3Log)+
  geom_bar(aes(x=row, y=value, col=row, fill=row),stat = "identity")+
  ggtitle("MKRN3")+
  theme_bw()+ scale_color_manual(values=c("darkgrey","darkred","darkgreen","red","green","#009A17"))+ 
  scale_fill_manual(values=c("white","white","white","white","white","white"))

ggsave("MM_PK_Imp_gene_Exp_FPKM_MKRN3Log.svg", width=15, height=15, units="cm", dpi=96)
barplot( height=MM_PK_Imp_gene_Exp_FPKM_MKRN3Log$value, names=MM_PK_Imp_gene_Exp_FPKM_MKRN3Log$row , density=c(12,12,12,12,12,12,12) , angle=c(45,45,45,45,45,45,45) , col=c("darkgrey","#EC7063","#48C9B0","#52BE80","#922B21","#148F77"),xlab = "Samples", ylab = "Log2(FPKM)",main = "MKRN3")
#MM_PK_Imp_gene_Exp_FPKM_MKRN3Log.png

