### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(reshape2)
library(gtools)
library(cowplot)
library(RColorBrewer)


seu_cleaned<-readRDS(here("scratch/HLiCA/healthy_RNA_merged_harmonized_gamma.rds"))

DotPlot(seu_cleaned, features=c("XIST","RPS4Y1"), group.by = "orig.ident")


coexp<-FetchData(seu_cleaned, vars = c("XIST","RPS4Y1") )
meta<-seu_cleaned@meta.data
identical(rownames(coexp), rownames(meta))

meta_exp<-cbind(meta, coexp)

ggplot(meta_exp, aes(SAMPLE, XIST))+geom_boxplot(outlier.size=0.25)+facet_wrap(~donor_sex, scales = "free_y")+coord_flip()+theme_bw()
ggplot(meta_exp, aes(SAMPLE, RPS4Y1))+geom_boxplot(outlier.size=0.25)+facet_wrap(~donor_sex, scales = "free_y")+coord_flip()+theme_bw()

# correlation plot
meta_exp$colour<-"Neither Expressed"
meta_exp$colour[which(meta_exp$XIST>0 )]<-"XIST only"
meta_exp$colour[which(meta_exp$RPS4Y1 >0)]<-"RPS4Y1 only"
meta_exp$colour[which(meta_exp$XIST>0 & meta_exp$RPS4Y1 >0)]<-"Both Expressed"

plot_count<-meta_exp %>%
  count(SAMPLE, colour) %>%
  arrange(SAMPLE, colour)

ggplot(plot_count, aes(x = SAMPLE, y = n, fill = colour)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Cell count", fill = "Colour category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = c(
    "red",  # light blue
    "#1f78b4",  # strong blue
    "forestgreen",  # light red/pink
    "grey"   # strong red
  ))


ggplot(meta_exp, aes(reorder(SAMPLE, XIST), XIST, fill=donor_sex))+
  geom_boxplot(outlier.size=0.25)+   theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))

ggplot(meta_exp, aes(reorder(SAMPLE, RPS4Y1), RPS4Y1, fill=donor_sex))+
  geom_boxplot(outlier.size=0.25)+   theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))

ggplot(meta_exp, aes(reorder(SAMPLE, RPS4Y1), 1, fill=donor_sex))+
  geom_tile(color="black")+   theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))


XIST_box<-ggplot(meta_exp, aes(reorder(SAMPLE, XIST), XIST, fill=donor_sex))+
  geom_boxplot(outlier.size=0.25)+   theme_bw()+
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))+
  ylab("XIST Expression")+xlab("Sample ID")

boxes<-ggplot(meta_exp, aes(reorder(SAMPLE, XIST), 1, fill=donor_sex))+
  geom_tile(color="black")+   theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))+
  xlab("")+ylab("")+  theme(axis.text.y = element_blank())

ggsave(plot_grid(XIST_box,boxes,ncol=1, align="v", rel_heights = c(2,1)),
       file="scratch/HLiCA/XIST_expression.png", width=20, h=10)



RPS4Y1_box<-ggplot(meta_exp, aes(reorder(SAMPLE, RPS4Y1), RPS4Y1, fill=donor_sex))+
  geom_boxplot(outlier.size=0.25)+   theme_bw()+
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))+
  ylab("RPS4Y1 Expression (Y chromosome)")+xlab("Sample ID")

boxes<-ggplot(meta_exp, aes(reorder(SAMPLE, RPS4Y1), 1, fill=donor_sex))+
  geom_tile(color="black")+   theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))+
  xlab("")+ylab("")+  theme(axis.text.y = element_blank())

ggsave(plot_grid(RPS4Y1_box,boxes,ncol=1, align="v", rel_heights = c(2,1)),
       file="scratch/HLiCA/RPS4Y1_expression.png", width=20, h=10)

# percent both per sample
samples<-unique(meta_exp$SAMPLE)

percent_both_samples<-lapply(1:length(samples), function(x){
  sample_exp<-meta_exp[which(meta_exp$SAMPLE==samples[x]),]
  data.frame(sample=samples[x], both_percent=(nrow(sample_exp[which(sample_exp$XIST!=0 & sample_exp$RPS4Y1!=0),])/nrow(sample_exp))*100)
})

percent_both_samples<-do.call(rbind, percent_both_samples)



#### Y chr score
# downloaded list of all y chr genes
y_genes <- c(
  "SRY", "CSF2RA", "CD99", "SHOX", "IL3RA", "CRLF2", "VAMP7", "SLC25A6", "TSPY1", "DAZ1",
  "USP9Y", "DDX3Y", "IL9R", "UTY", "KDM5D", "ASMT", "PCDH11Y", "RBMY1A1", "ZFY", "NLGN4Y",
  "DAZ2", "RPS4Y1", "EIF1AY", "DAZ4", "ZBED1", "CDY1", "DAZ3", "AKAP17A", "P2RY8", "AMELY",
  "SPRY3", "PPP2R3B", "TBL1Y", "GTPBP6", "PLCXD1", "CDY1B", "HSFY1", "BPY2", "CDY2A", "RPS4Y2",
  "TMSB4Y", "ASMTL", "DHRSX", "CDY2B", "VCY", "HSFY2", "RBMY1F", "TGIF2LY", "RBMY1B", "VCY1B",
  "RBMY1J", "RBMY1D", "TSPY10", "BPY2B", "TSPY9", "RBMY1E", "BPY2C", "TSPY2", "TSPY4", "TSPY3", "TSPY8"
)

ygenes_variable<- list(c(rownames(seu_cleaned)[which(rownames(seu_cleaned)%in%y_genes)]))

seu_cleaned    <- AddModuleScore(
  object = seu_cleaned,
  features = ygenes_variable,
  ctrl = 5,
  name = 'Y_chr'
)

meta_score<-seu_cleaned@meta.data


score_box<-ggplot(meta_score, aes(reorder(SAMPLE, Y_chr1), Y_chr1, fill=donor_sex))+
  geom_boxplot(outlier.size=0.25)+   theme_bw()+
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))+
  ylab("Y_chr score")+xlab("Sample ID")

boxes<-ggplot(meta_score, aes(reorder(SAMPLE, Y_chr1), 1, fill=donor_sex))+
  geom_tile(color="black")+   theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c("#b2182b","cornflowerblue"))+
  xlab("")+ylab("")+  theme(axis.text.y = element_blank())

ggsave(plot_grid(score_box,boxes,ncol=1, align="v", rel_heights = c(2,1)),
       file="scratch/HLiCA/Y_chr_score.png", width=20, h=10)



# C59 sspecific differences
# C59<-subset(d10x, subset = donor == "C59")
# C59_blend<-FeaturePlot(C59, features = c("XIST", "RPS4Y1"), blend = TRUE)
# ggsave(C59_blend, file="figures/jpeg/C59_XIST_RPS4Y1.jpeg", width = 10, height=3)
#
# FeaturePlot(C59, features = c("XIST", "RPS4Y1"), blend = TRUE)
# FeaturePlot(C59, features = c("XIST", "RPS4Y1"), blend = TRUE, split.by = "assay")
# VlnPlot(C59, features = c("XIST", "RPS4Y1"), group.by = "assay", split.by = "assay")

