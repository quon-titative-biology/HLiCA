### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(colorspace)
library(cowplot)
library(SeuratDisk)
library(harmony)


source("scripts/00_pretty_plots.R")
source("scripts/00_fanciest_UMAP.R")
source("scripts/colour_palette.R")






# load peds
load(here("~/links/scratch/HLiCA/Adult_ped_integrated.rds"))
DimPlot(d10x.combined_healthy)



#######################
## mesenchyme
#######################
# subset just to peds becasue some samples are in HLiCA
peds_healthy<-subset(d10x.combined_healthy, subset = age_condition =="Ped Healthy")
peds_healthy_mes<-subset(peds_healthy, subset = CellType_refined %in% c("HSC"))


# load HLiCA mesenchyme only
load(here("~/links/scratch/HLiCA/seu_mesenchyme.RData"))


#batch key is sample/library so make sure it is in both objects
peds_healthy_mes$batch<-peds_healthy_mes$sample
seu_cleaned$batch<-seu_cleaned$library_uuid

peds_healthy_mes$dataset<-"Pediatric"
seu_cleaned$dataset<-"HLiCA"


# merge
combined <- merge(seu_cleaned, y = peds_healthy_mes, add.cell.ids = c("cell_barcode", "cell"), project = "HarmonyIntegration")

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 50)

save(combined, file="~/links/scratch/HLiCA/mesenchyme_peds_hlica_combined.RData")

load(here("/media/redgar/Seagate Portable Drive/HLiCA/peds_integration/mesenchyme_peds_hlica_combined.RData"))

# # Harmony by dataset
# combined <- RunHarmony(combined, group.by.vars = "dataset", reduction = "pca", assay.use = "RNA", reduction.save = "harmony")

# # Harmony by library
combined <- RunHarmony(combined, group.by.vars = "batch", reduction = "pca", assay.use = "RNA", reduction.save = "harmony")

# plots
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.2)

save(combined, file="/media/redgar/Seagate Portable Drive/HLiCA/peds_integration/mesenchyme_peds_hlica_integrated.RData")



######################
load(here("/media/redgar/Seagate Portable Drive/HLiCA/peds_integration/mesenchyme_peds_hlica_integrated.RData"))

DimPlot(combined, reduction = "umap", group.by = "dataset")
DimPlot(combined, reduction = "umap", label = TRUE) + NoLegend()


DimPlot(combined, reduction = "umap", group.by = "Gamma.Annotation")+colscale_rnd

table(combined$seurat_clusters, combined$dataset)
#83 NRXN1 high cells in ped map (2,210 mesenchyme cells in peds)

FeaturePlot(combined, features = "NRXN1", split.by = "dataset")


## fancy UMAP
UMAP <- as.data.frame(Embeddings(combined, reduction = "umap"))
meta<-combined@meta.data

UMAP$cell<-rownames(UMAP)
meta$cell<-rownames(meta)
plt_umap<-merge(meta, UMAP, by="cell")

plt_umap$Gamma.Annotation[which(plt_umap$dataset=="Pediatric")]<-"Pediatric mesenchyme"

len_x_bar<-((range(plt_umap$umap_1))[2]-(range(plt_umap$umap_1))[1])/10
len_y_bar<-((range(plt_umap$umap_2))[2]-(range(plt_umap$umap_2))[1])/10
arr <- list(x = min(plt_umap$umap_1)-2, y = min(plt_umap$umap_2)-2, x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(fill=Gamma.Annotation),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_rnd+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)


fanciest_UMAP<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(color=Gamma.Annotation),size=0.05)+
  geom_point(data=plt_umap[which(plt_umap$dataset=="Pediatric"),], size = 0.06, colour= "black", stroke = 0.75)+
  geom_point(aes(color=Gamma.Annotation),data=plt_umap[which(plt_umap$dataset=="Pediatric"),], size=0.05)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_rnd+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")
fanciest_UMAP



fanciest_UMAP_cluster<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(color=seurat_clusters),size=0.05)+
  geom_point(data=plt_umap[which(plt_umap$seurat_clusters=="6"),], size = 0.06, colour= "black", stroke = 0.75)+
  geom_point(aes(color=seurat_clusters),data=plt_umap[which(plt_umap$seurat_clusters=="6"),], size=0.05)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")
fanciest_UMAP_cluster

celltype_label<-plt_umap %>%  group_by(seurat_clusters) %>%
  summarize(
    mean_umap1 = mean(umap_1, na.rm = TRUE),
    mean_umap2 = mean(umap_2, na.rm = TRUE)  )

fanciest_UMAP_cluster <- fanciest_UMAP_cluster+
  geom_text(aes(mean_umap1, mean_umap2, label=seurat_clusters), 
            data=celltype_label,size=2.5, color="black",fontface = "bold")
fanciest_UMAP_cluster

fanciest_UMAP_cluster 

save_plts(plot_grid(plot_grid(fanciest_UMAP_cluster,fanciest_UMAP, ncol=2, align="vh"), nice_legend, rel_widths = c(4,1)), 
          "peds_hlica_integration_umap", w=8, h=3)


#########
# peds clustering alone
#########
# load peds
load(here("/media/redgar/Seagate Portable Drive/HLiCA/peds_integration/Adult_ped_integrated.rds"))

# subset just to peds becasue some samples are in HLiCA
peds_healthy<-subset(d10x.combined_healthy, subset = age_condition =="Ped Healthy")
rm(d10x.combined_healthy)
gc()
peds_healthy_mes<-subset(peds_healthy, subset = CellType_refined %in% c("HSC"))
rm(peds_healthy)
gc()

peds_healthy_mes <- NormalizeData(peds_healthy_mes)
peds_healthy_mes <- FindVariableFeatures(peds_healthy_mes, selection.method = "vst", nfeatures = 3000)
peds_healthy_mes <- ScaleData(peds_healthy_mes)
peds_healthy_mes <- RunPCA(peds_healthy_mes, npcs = 50)

peds_healthy_mes <- RunUMAP(peds_healthy_mes, reduction = "pca", dims = 1:20)
peds_healthy_mes <- FindNeighbors(peds_healthy_mes, reduction = "pca", dims = 1:20)
peds_healthy_mes <- FindClusters(peds_healthy_mes, resolution = 0.1)

## fancy UMAP
UMAP <- as.data.frame(Embeddings(peds_healthy_mes, reduction = "umap"))
meta<-peds_healthy_mes@meta.data

UMAP$cell<-rownames(UMAP)
meta$cell<-rownames(meta)
plt_umap<-merge(meta, UMAP, by="cell")


gene_exp<-FetchData(peds_healthy_mes, vars=c("NRXN1"))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt_umap<-merge(plt_umap, gene_exp, by="cell")

plt_umap<-plt_umap[order(plt_umap$value),]

len_x_bar<-((range(plt_umap$umap_1))[2]-(range(plt_umap$umap_1))[1])/10
len_y_bar<-((range(plt_umap$umap_2))[2]-(range(plt_umap$umap_2))[1])/10
arr <- list(x = min(plt_umap$umap_1)-2, y = min(plt_umap$umap_2)-2, x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(fill=value),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_gradientn(colours = c("grey92", "#bcd4ff", "#0448bd"))+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)


fanciest_UMAP<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(color=value),size=0.05)+
  geom_point(data=plt_umap[which(plt_umap$seurat_clusters=="5"),], size = 0.06, colour= "black", stroke = 0.75)+
  geom_point(aes(color=value),data=plt_umap[which(plt_umap$seurat_clusters=="5"),], size=0.05)+
  scale_color_gradientn(colours = c("grey92", "#bcd4ff", "#0448bd")) +  # Adjust colors
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")
fanciest_UMAP

fanciest_UMAP_cluster<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(color=seurat_clusters),size=0.05)+
  geom_point(data=plt_umap[which(plt_umap$seurat_clusters=="5"),], size = 0.06, colour= "black", stroke = 0.75)+
  geom_point(aes(color=seurat_clusters),data=plt_umap[which(plt_umap$seurat_clusters=="5"),], size=0.05)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")
fanciest_UMAP_cluster

celltype_label<-plt_umap %>%  group_by(seurat_clusters) %>%
  summarize(
    mean_umap1 = mean(umap_1, na.rm = TRUE),
    mean_umap2 = mean(umap_2, na.rm = TRUE)  )

fanciest_UMAP_cluster <- fanciest_UMAP_cluster+
  geom_text(aes(mean_umap1, mean_umap2, label=seurat_clusters), 
            data=celltype_label,size=1.5, color="black",fontface = "bold")
fanciest_UMAP_cluster

save_plts(plot_grid(plot_grid(fanciest_UMAP_cluster,fanciest_UMAP, ncol=2, align="vh"), nice_legend, rel_widths = c(10,1)), 
          "peds_NRXN1_UMAP", w=5, h=2)



table(peds_healthy_mes$seurat_clusters)
table(peds_healthy_mes$seurat_clusters, peds_healthy_mes$sample)
# 85 NRNX1, 69 from C115


cell_cluster_count<-peds_healthy_mes@meta.data %>%  group_by(seurat_clusters,sample %>% eval) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

cell_cluster_count<-as.data.frame(cell_cluster_count)
cell_cluster_count_NRXN<-cell_cluster_count[which(cell_cluster_count$seurat_clusters=="5"),]
colnames(cell_cluster_count_NRXN)[2]<-"Donor"

ped_NRXN1_bar<-ggplot() + 
  geom_bar(aes(fill=Donor, y=n, x=seurat_clusters),cell_cluster_count_NRXN, position="stack", stat="identity", color="black")+
  theme_classic()+th+ylab("Cell Count")+xlab("Cluster")+coord_flip()+theme(legend.position = "bottom")
ped_NRXN1_bar
save_plts(ped_NRXN1_bar, "peds_NRXN1_bar", w=6, h=2)










################################
## lymphocytes
################################
# subset just to peds becasue some samples are in HLiCA
peds_healthy<-subset(d10x.combined_healthy, subset = age_condition =="Ped Healthy")
peds_healthy_lymphocytes<-subset(peds_healthy, subset = CellType_refined %in% c("CD3+ T-cells","gd T-cells",
                                                                                "Mature B-cells","NK-like cells",
                                                                                "Plasma cells","Cycling Plasma"))
# load HLiCA lympohocytes only
load(here("~/links/scratch/HLiCA/seu_cleaned_lymphocytes.RData"))

#batch key is sample/library so make sure it is in both objects
peds_healthy_lymphocytes$batch<-peds_healthy_lymphocytes$sample
seu_cleaned$batch<-seu_cleaned$library_uuid

peds_healthy_lymphocytes$dataset<-"Pediatric"
seu_cleaned$dataset<-"HLiCA"



load(here("/media/redgar/Seagate Portable Drive/HLiCA/peds_integration/lymphocytes_peds_hlica_combined.RData"))

# Harmony by dataset
# combined <- RunHarmony(combined, group.by.vars = "dataset", reduction = "pca", assay.use = "RNA", reduction.save = "harmony")

# # Harmony by library
combined <- RunHarmony(combined, group.by.vars = "batch", reduction = "pca", assay.use = "RNA", reduction.save = "harmony")

# plots
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.4)

save(combined, file="/media/redgar/Seagate Portable Drive/HLiCA/peds_integration/lymphocytes_peds_hlica_integrated.RData")



#############################################
load(here("/media/redgar/Seagate Portable Drive/HLiCA/peds_integration/lymphocytes_peds_hlica_integrated.RData"))

DimPlot(combined, reduction = "umap", group.by = "dataset")
DimPlot(combined, reduction = "umap", label = TRUE) + NoLegend()


DimPlot(combined, reduction = "umap", group.by = "Gamma.Annotation")+colscale_rnd

table(combined$seurat_clusters, combined$dataset)
#69 T regs in ped map

FeaturePlot(combined, features = "IL2RA", split.by = "dataset")


## fancy UMAP
UMAP <- as.data.frame(Embeddings(combined, reduction = "umap"))
meta<-combined@meta.data

UMAP$cell<-rownames(UMAP)
meta$cell<-rownames(meta)
plt_umap<-merge(meta, UMAP, by="cell")

plt_umap$Gamma.Annotation[which(plt_umap$dataset=="Pediatric")]<-"Pediatric lymphocytes"

len_x_bar<-((range(plt_umap$umap_1))[2]-(range(plt_umap$umap_1))[1])/10
len_y_bar<-((range(plt_umap$umap_2))[2]-(range(plt_umap$umap_2))[1])/10
arr <- list(x = min(plt_umap$umap_1)-2, y = min(plt_umap$umap_2)-2, x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(fill=Gamma.Annotation),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  fillscale_rnd+theme_bw()+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)


fanciest_UMAP<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(color=Gamma.Annotation),size=0.05)+
  geom_point(data=plt_umap[which(plt_umap$dataset=="Pediatric"),], size = 0.06, colour= "black", stroke = 0.75)+
  geom_point(aes(color=Gamma.Annotation),data=plt_umap[which(plt_umap$dataset=="Pediatric"),], size=0.05)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  colscale_rnd+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")
fanciest_UMAP



fanciest_UMAP_cluster<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(color=seurat_clusters),size=0.05)+
  geom_point(data=plt_umap[which(plt_umap$seurat_clusters=="12"),], size = 0.06, colour= "black", stroke = 0.75)+
  geom_point(aes(color=seurat_clusters),data=plt_umap[which(plt_umap$seurat_clusters=="12"),], size=0.05)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")
fanciest_UMAP_cluster

celltype_label<-plt_umap %>%  group_by(seurat_clusters) %>%
  summarize(
    mean_umap1 = mean(umap_1, na.rm = TRUE),
    mean_umap2 = mean(umap_2, na.rm = TRUE)  )

fanciest_UMAP_cluster <- fanciest_UMAP_cluster+
  geom_text(aes(mean_umap1, mean_umap2, label=seurat_clusters), 
            data=celltype_label,size=2.5, color="black",fontface = "bold")
fanciest_UMAP_cluster

save_plts(plot_grid(plot_grid(fanciest_UMAP_cluster,fanciest_UMAP, ncol=2, align="vh"), nice_legend, rel_widths = c(4,1)), 
          "peds_hlica_integration_umap_lymphocytes", w=10, h=5)


#############################################
# peds clustering alone

# load peds
load(here("/media/redgar/Seagate Portable Drive/HLiCA/peds_integration/Adult_ped_integrated.rds"))

# subset just to peds becasue some samples are in HLiCA
peds_healthy<-subset(d10x.combined_healthy, subset = age_condition =="Ped Healthy")
rm(d10x.combined_healthy)
gc()
peds_healthy_lymphocytes<-subset(peds_healthy, subset = CellType_refined %in% c("CD3+ T-cells","gd T-cells",
                                                                                "Mature B-cells","NK-like cells",
                                                                                "Plasma cells","Cycling Plasma"))
rm(peds_healthy)
gc()

peds_healthy_lymphocytes <- NormalizeData(peds_healthy_lymphocytes)
peds_healthy_lymphocytes <- FindVariableFeatures(peds_healthy_lymphocytes, selection.method = "vst", nfeatures = 3000)
peds_healthy_lymphocytes <- ScaleData(peds_healthy_lymphocytes)
peds_healthy_lymphocytes <- RunPCA(peds_healthy_lymphocytes, npcs = 50)

peds_healthy_lymphocytes <- RunUMAP(peds_healthy_lymphocytes, reduction = "pca", dims = 1:20)
peds_healthy_lymphocytes <- FindNeighbors(peds_healthy_lymphocytes, reduction = "pca", dims = 1:20)
peds_healthy_lymphocytes <- FindClusters(peds_healthy_lymphocytes, resolution = 0.1)


DimPlot(peds_healthy_lymphocytes, reduction = "umap", label = TRUE) + NoLegend()
FeaturePlot(peds_healthy_lymphocytes, features = "IL2RA")
DotPlot(peds_healthy_lymphocytes, features = "IL2RA")


## fancy UMAP
UMAP <- as.data.frame(Embeddings(peds_healthy_lymphocytes, reduction = "umap"))
meta<-peds_healthy_lymphocytes@meta.data

UMAP$cell<-rownames(UMAP)
meta$cell<-rownames(meta)
plt_umap<-merge(meta, UMAP, by="cell")

gene_exp<-FetchData(peds_healthy_lymphocytes, vars=c("IL2RA"))
gene_exp$cell<-rownames(gene_exp)
gene_exp<-melt(gene_exp)

plt_umap<-merge(plt_umap, gene_exp, by="cell")

plt_umap<-plt_umap[order(plt_umap$value),]

len_x_bar<-((range(plt_umap$umap_1))[2]-(range(plt_umap$umap_1))[1])/10
len_y_bar<-((range(plt_umap$umap_2))[2]-(range(plt_umap$umap_2))[1])/10
arr <- list(x = min(plt_umap$umap_1)-2, y = min(plt_umap$umap_2)-2, x_len = len_x_bar, y_len = len_y_bar)

forlegned_plot<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(aes(fill=value),size=2, shape=21)+xlab("UMAP 1")+ylab("UMAP 2")+
  scale_fill_gradientn(colours = c("grey92", "#bcd4ff", "#0448bd"))+
  theme(legend.text = element_text(size=5),
        legend.title = element_text(size=6))
nice_legend<-get_leg(forlegned_plot)


fanciest_UMAP<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(data=plt_umap, size = 0.06, colour= "black", stroke = 0.75)+
  geom_point(aes(color=value),data=plt_umap, size=0.05)+
  scale_color_gradientn(colours = c("grey92", "#bcd4ff", "#0448bd")) +  # Adjust colors
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")
fanciest_UMAP

fanciest_UMAP_cluster<-ggplot(plt_umap, aes(umap_1,umap_2))+
  geom_point(data=plt_umap, size = 0.06, colour= "black", stroke = 0.75)+
  geom_point(aes(color=seurat_clusters),data=plt_umap, size=0.05)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), size=0.25,color="black",
           arrow = arrow(type = "closed", length = unit(2, 'pt'))) +
  theme_void()+theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"),
                     axis.title.x = element_text(size=5,hjust = 0.05),
                     axis.title.y = element_text(size=5,hjust = 0.05,angle = 90),
                     legend.position = "none")
fanciest_UMAP_cluster

celltype_label<-plt_umap %>%  group_by(seurat_clusters) %>%
  summarize(
    mean_umap1 = mean(umap_1, na.rm = TRUE),
    mean_umap2 = mean(umap_2, na.rm = TRUE)  )

fanciest_UMAP_cluster <- fanciest_UMAP_cluster+
  geom_text(aes(mean_umap1, mean_umap2, label=seurat_clusters), 
            data=celltype_label,size=1.5, color="black",fontface = "bold")
fanciest_UMAP_cluster

save_plts(plot_grid(plot_grid(fanciest_UMAP_cluster,fanciest_UMAP, ncol=2, align="vh"), nice_legend, rel_widths = c(10,1)), 
          "peds_treg_UMAP", w=5, h=2)
