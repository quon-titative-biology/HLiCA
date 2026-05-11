

myColors_celltype <- c(
  "#994676",  # Kupffer Cells
  "#ce1256",  # Classical Monocytes
  "#d9667f",  # Type 2 cDCs
  "#ce1256",  # Non-classical Monocytes
  "#b80783",  # TREM2+ Macrophages
  "#b01629",  # MAMLD1+ Trans Monocytes
  "#d9667f",  # DCmac
  "#914e84",  # Neutrophils
  "#b01629",  # Activated Monocytes
  "#d9667f",  # migDCs
  "#d9667f",  # Type 1 cDCs
  "#ce1256",  # BAG3+ Monocytes
  
  "#FFA500",  # Small ApoLipo 
  "#FF8C00",  # Small Keratin 
  "#FFB347",  # LAMC2+ Small 
  "#FF7F50",  # CXCL8+ Small Keratin
  "#FFAA33",   # Large Mucus Secreting 
  
  "#FFA500",  # ApoLipo 
  "#FF8C00",  # Keratin 
  "#FFB347",  # LAMC2+ Cholangiocyte 
  "#FF7F50",  # CXCL8+  Keratin
  "#FFAA33",  # Mucus Secreting 
  
  "#8B4513",  # Periportal LSEC 
  "#A0522D",  # Portal Vein 
  "#CD853F",  # Central Venous LSEC 
  "#B5651D",  # Lymphatic 
  "#DEB887",  # Hepatic Artery 
  "#8B3E2F",  # LSEC 
  "#A52A2A",   # Vascular Endothelial Cell 
  
  "#F0E68C",  # Hepatic Stellate Cell 
  "#FFD700",  # Vascular Smooth Muscle Cell 
  "#FFFACD",  # Portal Fibroblast
  "#F5DEB3",   # CUX2+ Hepatic Stellate Cell
  "#f09d02",  "#f09d02",
  
  "#663399",  "#800080",  '#9370DB',  '#9851b0',
  
  "#006400", "#339500", "#4CAF00", "#66C700",  "#b0db93", "#5f9c4e", "#4c9638",
  "#6495ED", "#1d91c0", "#74a9cf", "#9ecae1", "#0570b0", "#08519c", 
  
  "#eb150c",  # pDC
  "#000066",  # Cycling
  "grey",   # Erythrocyte
  "grey50","grey70","grey90"
)
color_possibilities_celltype<-c(
  "Kupffer Cells", "Classical Monocytes", "Type 2 cDCs",
  "Non-classical Monocytes", "TREM2+ Macrophages", "MAMLD1+ Trans Monocytes", "DCmac",
  "Neutrophils", "Activated Monocytes", "migDCs", "Type 1 cDCs",
  "BAG3+ Monocytes", 
  
  "Small ApoLipo", "Small Keratin", "LAMC2+ Small","CXCL8+ Small Keratin", "Large Mucus Secreting", 
  "ApoLipo", "Keratin", "LAMC2+ Cholangiocyte","CXCL8+ Keratin", "Mucus Secreting", 
  
  "Periportal LSEC", "Portal Vein", "Central Venous LSEC", "Lymphatic", "Hepatic Artery", "Vascular Endothelial Cell", "LSEC", "Hepatic Stellate Cell", "Vascular Smooth Muscle Cell",
  "Portal Fibroblast", "CUX2+ Hepatic Stellate Cell", "NRXN1+ Hepatic Stellate Cell","NRXN1+ Stromal Cell",
  
  "Plasma B Cell", "B Cell", "IgA B cells", "IgG B cells",
  "Bright NK Cell", "Dim NK Cell", "Helper T Cell", "CD8 T Cell", "MAIT T Cell", 
  "Exhausted CD8 T Cell",  "Regulatory T Cell",
  
  "Ribosomal+ Hepatocyte","UGT+ Hepatocyte" ,"Mito+ Hepatocyte" , "Pericentral Hepatocyte", "Periportal Hepatocyte","SERPINE1+ Hepatocyte",
  
  "pDC","Cycling", "Erythrocyte",
  "Doublet seen in beta annotation","Doublet seen in alpha annotation","Clustered by technical factors once sorted into alpha lineage"
)


names(myColors_celltype) <- color_possibilities_celltype
fillscale_cellType <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)
colscale_cellType <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype, drop = T, limits=force)







myColors_celltype_alpha <- c(
  "#ce1256",  # Myeloid Cell
  "#FFA500",  # Cholangiocyte
  "#CD853F",  # Endothelia
  "#FFD700",  # Mesenchyme
  "#339500", #Lymphocyte
  "#6495ED", #Hepatocyte
  "pink",  # Erythrocyte
  "#eb150c",  # pDC
  "#000066",  # Cycling
  "grey90",  # Potential Doublet
  "grey", #Doublet
  "grey", #removed
  "grey" #TBC
  
)
color_possibilities_celltype_alpha<-c(
  "Myeloid Cell", "Cholangiocyte","Endothelia",
  "Mesenchyme","Lymphocyte","Hepatocyte","Erythrocyte",
  "pDC","Cycling","Potential Doublet","Doublet","removed","TBC")


names(myColors_celltype_alpha) <- color_possibilities_celltype_alpha
fillscale_alpha <- scale_fill_manual(name="Cell Type",
                                        values = myColors_celltype_alpha, drop = T, limits=force)
colscale_alpha <- scale_color_manual(name="Cell Type",
                                        values = myColors_celltype_alpha, drop = T, limits=force)



myColors_rnd <-  c(
    "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#654522", 
    "#8DB600", "#008856", "#E68FAC", "#0067A5", "#F99379", "#604E97",
    # Group separator
    "#F6A600", "#0067A5", "#DCD300", "#8DB600","#eb150c",
    "#F6A600", "#0067A5", "#DCD300", "#8DB600","#eb150c",
    # Group separator
    "#0067A5", "#E25822", "#B3446C", "#F3C300", "#875692", "#F38400", "#8DB600",
    # Group separator
    "#A1CAF1", "#BE0032", "#8DB600", "#DCD300", "#67038f", "#67038f", 
    # Group separator
    "#7ecfba", "#0067A5", "#F99379", "#6A5ACD", 
    "#F6A600", "#B3446C", "#8510b0", "#3b7823", "#A1CAF1", "#BE0032", 
    "#DCD300", "#882D17", "#8510b0", "#F38400", "#3b7823", "#0067A5", 
    "#B3446C",  "#eb150c",  # pDC
    "#000066",  # Cycling
    # Group separator
    "grey", "grey50", "grey70", "grey90")

names(myColors_rnd) <- color_possibilities_celltype
fillscale_rnd <- scale_fill_manual(name="Cell Type",
                                        values = myColors_rnd, drop = T, limits=force)
colscale_rnd <- scale_color_manual(name="Cell Type",
                                        values = myColors_rnd, drop = T, limits=force)




###################
### for sankey2
###################
myColors_sankey <- c(
  "#ce1256",  # Myeloid Cell
  "#FFA500",  # Cholangiocyte
  "#CD853F",  # Endothelia
  "#FFD700",  # Mesenchyme
  "#339500",  # Lymphocyte
  "#6495ED",  # Hepatocyte
  "pink",     # Erythrocyte
  "#000066",  # Cycling
  "grey90",   # Potential Doublet
  "grey",     # Doublet
  "grey",     # removed
  "grey",     # TBC

  "#994676",  # Kupffer Cells
  "#ce1256",  # Classical Monocytes
  "#d9667f",  # Type 2 cDCs
  "#ce1256",  # Non-classical Monocytes
  "#b80783",  # TREM2+ Macrophages
  "#b01629",  # MAMLD1+ Trans Monocytes
  "#d9667f",  # DCmac
  "#914e84",  # Neutrophils
  "#b01629",  # Activated Monocytes
  "#d9667f",  # migDCs
  "#d9667f",  # Type 1 cDCs
  "#ce1256",  # BAG3+ Monocytes

  "#FFA500",  # Small ApoLipo
  "#FF8C00",  # Small Keratin
  "#FFB347",  # LAMC2+ Small
  "#FF7F50",  # CXCL8+ Small Keratin
  "#FFAA33",  # Large Mucus Secreting
  
  "#FFA500",  # Small ApoLipo
  "#FF8C00",  # Small Keratin
  "#FFB347",  # LAMC2+ Small
  "#FF7F50",  # CXCL8+ Small Keratin
  "#FFAA33",  # Large Mucus Secreting

  "#8B4513",  # Periportal LSEC
  "#A0522D",  # Portal Vein
  "#CD853F",  # Central Venous LSEC
  "#B5651D",  # Lymphatic
  "#DEB887",  # Hepatic Artery
  "#CD853F",  # LSEC
  "#CD853F",  # Vascular Endothelial Cell

  "#F0E68C",  # Hepatic Stellate Cell
  "#FFD700",  # Vascular Smooth Muscle Cell
  "#FFFACD",  # Portal Fibroblast
  "#F5DEB3",  # CUX2+ Hepatic Stellate Cell
  "#f09d02",  # NRXN1+ Hepatic Stellate Cell
  "#f09d02",  # NRXN1+ Hepatic Stellate Cell

  "#663399", "#800080",  '#9370DB',  '#9851b0',

  "#006400", "#339500", "#4CAF00", "#66C700",  "#b0db93", "#5f9c4e", "#4c9638",
  "#6495ED", "#1d91c0", "#74a9cf", "#9ecae1", "#0570b0", "#08519c",

  "#eb150c",  # pDC
  "#000066",  # Cycling
  "grey",     # Erythrocyte
  "grey50", "grey70", "grey90",

  # betas
  "#000066",  # Cycling Cells
  "#eb150c", # "pDCs"
  "#ce1256",  # Non Classical Monocytes
  "#873b57",  # "LAM-Like"
  "#d9667f",# "Macrophages (cDCs?)"
  "#63314b",# "Activated Macrophages"
  "#d9667f",# "Type 1 cDCs"
  "#4CAF00",# "T Cells?"
  "#FF8C00",# "Keratin / Small / Ribosomal / Cycling"
  "#FFA500",# "ApoLipo / Small"
  "#FFB347",# "LAMC2+ / Small"
  "#FF8C00",# "Small (Keratin)"
  "#b8ad93",# "Stressed?"
  "#FFAA33", # "Mucus Secreting / Large / Keratin"
  "#FF8C00",# "Keratin / CXCL8+"
  "#CD853F",# "Central Vein"
  '#800080',# "Naive B Cell"
  "#006400", # "cNK / Bright NK Cell"
  "#339500",# "Liver Resident / Dim NK Cell"
  "#66C700",# "CD4+ / Helper T Cell"
  "#5f9c4e", # "CD8+ / Cytotoxic T Cell"
  "#b0db93",# "CD4+ / MAIT T Cell"
  "#4c9638",# "CD4+ / Regulatory T Cell"
  "#8ea389",# "Stressed Cells"
  "#8ea389",# "Poor Quality"
  
  "#6495ED", # "Ribosomal+ / APOA2+"
  "#1d91c0",# "ORM2+ / APOE+"
  "#74858c",# "Potentially Technical"
  "#74a9cf",# "Mito+"
  "#6495ED",# "Ribosomal+ / FTL+"
  "#9ecae1",# "Central"
  "#0570b0",# "Portal"
  "#08519c" # "SERPINE1+"
)

color_possibilities_sankey <- c(
  "Myeloid Cell", "Cholangiocyte", "Endothelia",
  "Mesenchyme", "Lymphocyte", "Hepatocyte", "Erythrocyte",
  "Cycling", "Potential Doublet", "Doublet", "removed", "TBC",

  "Kupffer Cells", "Classical Monocytes", "Type 2 cDCs",
  "Non-classical Monocytes", "TREM2+ Macrophages", "MAMLD1+ Trans Monocytes", "DCmac",
  "Neutrophils", "Activated Monocytes", "migDCs", "Type 1 cDCs",
  "BAG3+ Monocytes", 
  
  "Small ApoLipo", "Small Keratin", "LAMC2+ Small","CXCL8+ Small Keratin", "Large Mucus Secreting", 
  "ApoLipo", "Keratin", "LAMC2+ Cholangiocyte","CXCL8+ Keratin", "Mucus Secreting", 
  
  
  "Periportal LSEC", "Portal Vein", "Central Venous LSEC", "Lymphatic", "Hepatic Artery", "Vascular Endothelial Cell", "LSEC", "Hepatic Stellate Cell", "Vascular Smooth Muscle Cell",
  "Portal Fibroblast", "CUX2+ Hepatic Stellate Cell", "NRXN1+ Hepatic Stellate Cell","NRXN1+ Stromal Cell",
  
  "Plasma B Cell", "B Cell", "IgA B cells", "IgG B cells",
  "Bright NK Cell", "Dim NK Cell", "Helper T Cell", "CD8 T Cell", "MAIT T Cell", 
  "Exhausted CD8 T Cell",  "Regulatory T Cell",
  
  "Ribosomal+ Hepatocyte","UGT+ Hepatocyte" ,"Mito+ Hepatocyte" , "Pericentral Hepatocyte", "Periportal Hepatocyte","SERPINE1+ Hepatocyte",
  
  "pDC","Cycling", "Erythrocyte",
  "Doublet seen in beta annotation","Doublet seen in alpha annotation","Clustered by technical factors once sorted into alpha lineage",
  
  #betas
  "Cycling Cells", "pDCs", "Non Classical Monocytes", "LAM-Like",
  "Macrophages (cDCs?)", "Activated Macrophages", "Type 1 cDCs", "T Cells?",
  "Keratin / Small / Ribosomal / Cycling", "ApoLipo / Small", "LAMC2+ / Small", "Small (Keratin)",
  "Stressed?", "Mucus Secreting / Large / Keratin", "Keratin / CXCL8+", "Central Vein",
  "Naive B Cell", "cNK / Bright NK Cell", "Liver Resident / Dim NK Cell", "CD4+ / Helper T Cell",
  "CD8+ / Cytotoxic T Cell", "CD4+ / MAIT T Cell", "CD4+ / Regulatory T Cell", "Stressed Cells",
  "Poor Quality", "Ribosomal+ / APOA2+", "ORM2+ / APOE+", "Potentially Technical",
  "Mito+", "Ribosomal+ / FTL+", "Central", "Portal",
  "SERPINE1+"
  
)


names(myColors_sankey) <- color_possibilities_sankey
fillscale_sankey <- scale_fill_manual(name="Cell Type",
                                      values = myColors_sankey, drop = T, limits=force)
colscale_sankey <- scale_color_manual(name="Cell Type",
                                      values = myColors_sankey, drop = T, limits=force)


##########
## GSEA
##########
myColors_gsea <- c(
  "#6495ED", 
  "#BE0032"
)
color_possibilities_gsea<-c(
  "Down-regulated", "Up-regulated")


names(myColors_gsea) <- color_possibilities_gsea
fillscale_gsea <- scale_fill_manual(name="Cell Type",
                                    values = myColors_gsea, drop = T, limits=force)
colscale_gsea <- scale_color_manual(name="Cell Type",
                                    values = myColors_gsea, drop = T, limits=force)




##########
## disease
##########
myColors_disease <- c(
  "gray70",  "#D73027", "#CDE077","#CDE077",
  "#005E73","#005E73", "#003C4C","#003C4C",
  "#9470B3", "#6B4E9B", 
  "#E69F00",   "#E69F00", 
  "#B86B00",   "#B86B00", 
  "#C94FB3", "#C94FB3",
  "#FFD700", "#228B22","#556B2F"
)
color_possibilities_disease<-c(
  "healthy" ,
  "Cholestasis",
  "Primary.sclerosing.cholangitis" ,  "Primary sclerosing cholangitis" ,
  "Intrahepatic.cholangiocarcinoma",  "Intrahepatic cholangiocarcinoma",
  "Extrahepatic.cholangiocarcinoma" ,  "Extrahepatic cholangiocarcinoma" ,
  "Fibrosis" ,
  "Cirrhosis" ,
  "Chronic.HBV" ,  "Chronic HBV" ,
  "nonA.E.hepatitis.acute.liver.failure",  "nonA-E hepatitis acute liver failure",
  "acetaminophen.induced.acute.liver.failure" ,  "acetaminophen-induced acute liver failure" ,
  "NAFLD" ,
  "HCC",
  "Hepatoblastoma")


names(myColors_disease) <- color_possibilities_disease
fillscale_disease <- scale_fill_manual(name="Disease",
                                    values = myColors_disease, drop = T, limits=force)
colscale_disease <- scale_color_manual(name="Disease",
                                    values = myColors_disease, drop = T, limits=force)


