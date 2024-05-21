# Load packages and functions for analysis. Henderson Labs package SeuratPipe is used but code elements can be pulled out and isolated for this analysis if needed.
library(Seurat)
seu <- readRDS("project/Liver_Cell_Atlas/Data/First_Download/all_harmonized_seurat.rds") # Load data

# Load in Atlas markers
# Atlas_marker_list is a copy of the Supplementarty table 2 that was downloaded from the google shared drive on the 5th Sept.
# Atlas marker list collapsed is the same list as the atlas marker list but with all the extra markers groups have added interdigitated into the first three columns.
# Updated collapsed marker list was generated after talking with Colleagues and changed a few of the markers around that I had added in for the Henderson group.
atlas_markers <- read.csv("~/project/Liver_Cell_Atlas/Analysis_01/00_Atlas_Markers/marker_lists/Atlas_Marker_List_Collapsed_Updated.csv")
atlas_markers <- atlas_markers[,1:3]

# Turn the dataframe into an appropriate list for the module score analysis
lst <- lapply(split(atlas_markers[,c(1,3)]['Gene'], atlas_markers$General.Type), transform, Gene = factor(Gene)) 
lst <- lapply(lst, as.vector) # Turn list of dataframes into list of vectors 
lst <- unlist(lst, recursive = FALSE) # Unlist the last vector list
names(lst) <- paste0("General_",c("B_cell", "Cell_Cycle", "Cholangiocyte", "Endothelial", "Hepatocyte",  "Lymphocyte",  "Macrophage_Myeloid",
                                  "PDC", "RBC", "Stellate_Mesenchymal")) # Rename as appropriate

lst_general <- lst
lst_general <- lapply(lst_general, as.character) # Change to character to remove levels
lst_general <- lapply(lst_general, unique) # Make sure each list has only one copy of the marker

# Repeat above pattern on the markers in the specific list as opposed to general list
lst <- lapply(split(atlas_markers[,c(1,3)]['Gene'], atlas_markers$Specific.Type), transform, Gene = factor(Gene)) 
lst <- lapply(lst, as.vector)
lst <- unlist(lst, recursive = FALSE)
names(lst) <- paste0("Specific_",c("abT_incl.CD4.and.CD8", "ActivatedMac", "AP1.Stellate",  "ApoLipoChol", "Arterial", "B.cell", "C.Hepato", "CD4.T", "CD8.T",
                                   "Cellcycle", "Cholangio", "cNK", "cvEndo", "cvLSEC", "Endothelial", "gdT", "Hepatocyte", "KeratinsChol", "Kupffer", "LAM.like", "lrNK", 
                                   "Macrophage", "Mast", "Mesenchyme", "MHCII_cDC", "Monocyte.derived", "Mononuclear.Phagocyte", "Mucus.Secreting", "Myeloid", "Naive.B.cell", 
                                   "Neutrophil", "P.Hepato", "PDC", "Plama.B.cell", "ppLSEC", "Quiescent.Stellate", "RBC", "ResidentMac", "Stellate.Fiber", "SynapseMac","T.cell", "VSMC"))
lst_specific <- lst
lst_specific <- lapply(lst_specific, as.character)
lst_specific <- lapply(lst_specific, unique)

# Give Henderson markers a defining name to allow to easy distinquishing of figures
# Henderson Human Liver lineage Markers
human_liver_lineage_markers <- list(
  lin_endo = c("PECAM1", "CDH5", "ERG", "ENG", "CLEC4G", "ICAM2", "KDR", "STAB1", "CD34"),
  lin_mesenchyme = c("PDGFRB", "PDGFRA", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "DES", "DCN"),
  lin_hepatocyte = c("ALB", "TTR", "TF", "HP", "HNF4A", "CYP2A6"),
  lin_cholangiocyte = c("EPCAM", "KRT19", "CD24"),
  lin_mp = c("CD68", "ITGAM", "ITGAX", "HLA-DRA", "CSF1R", "CD14"),
  lin_tcell = c("CD3D", "CD3E", "CD3G", "CD8A", "CD4"),
  lin_bcell = c("CD79A", "CD79B", "CD19", "MS4A1"),
  lin_ilc = c("KLRF1", "KLRC1", "GZMA", "GZMB", "NKG7"),
  lin_pdc = c("LILRA4", "CLEC4C", "GZMB"),
  lin_mast = c("KIT", "TPSAB1", "TPSB2"),
  lin_plasma = c("CD79A", "JCHAIN", "IGHA2"),
  lin_rbc = c("HBA1", "HBA2", "HBB"),
  lin_cycling = c("MKI67", "CCNA2", "CCNB2", "STMN1")
)

names(human_liver_lineage_markers) <- gsub("lin_", "Hendo_", names(human_liver_lineage_markers))
human_liver_lineage_markers <- lapply(human_liver_lineage_markers, unique)



# Check which genes are not present in the dataset and remove from each list
unlist(lst_specific)[!unlist(lst_specific) %in% rownames(seu)]
#Specific_Arterial3 Specific_Cholangio23 Specific_Hepatocyte7   Specific_Kupffer12 Specific_Macrophage3     Specific_ppLSEC1 
#"PLPP"             "GNB2L1"               "G6PC"              "SEPP1"               "CD54"              "ADRIF" 

tmp <- lst_specific
tmp$Specific_Arterial <- tmp$Specific_Arterial[-3]
tmp$Specific_Cholangio <- tmp$Specific_Cholangio[-23]
tmp$Specific_Hepatocyte <- tmp$Specific_Hepatocyte[-7]
tmp$Specific_Kupffer <- tmp$Specific_Kupffer[-12]
tmp$Specific_Macrophage <- tmp$Specific_Macrophage[-3]
tmp$Specific_ppLSEC <- tmp$Specific_ppLSEC[-1]
unlist(tmp)[!unlist(tmp) %in% rownames(seu)]
#named character(0)
lst_specific <- tmp



unlist(lst_general)[!unlist(lst_general) %in% rownames(seu)]
#General_Cholangiocyte27         General_Endothelial3        General_Endothelial35         General_Hepatocyte22 General_Macrophage_Myeloid23 General_Macrophage_Myeloid38 
#"GNB2L1"                       "PLPP"                      "ADRIF"                       "G6PC"                      "SEPP1"                       "CD54" 

tmp <- lst_general
tmp$General_Cholangiocyte <- tmp$General_Cholangiocyte[-27]
tmp$General_Endothelial <- tmp$General_Endothelial[-c(3,35)]
tmp$General_Hepatocyte <- tmp$General_Hepatocyte[-22]
tmp$General_Macrophage_Myeloid <- tmp$General_Macrophage_Myeloid[-c(23,38)]
unlist(tmp)[!unlist(tmp) %in% rownames(seu)]
lst_general <- tmp

unlist(human_liver_lineage_markers)[!unlist(human_liver_lineage_markers) %in% rownames(seu)]
# 0


# List of Module score lists
modules <- list(Atlas_general = lst_general,
                Atlas_specific = lst_specific,
                Hendo = human_liver_lineage_markers)


saveRDS(modules, "project/Liver_Cell_Atlas/Analysis_01/00_Atlas_Markers/modules.rds")
