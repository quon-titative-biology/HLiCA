library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

getwd()

sampleName<-"FRS1"
sampleFolder<-paste0(sampleName,"/")

##add some subfolders
dir.create(paste0(sampleFolder,"results"))
dir.create(paste0(sampleFolder,"Robjects"))

########################################
##### Functions
########################################
source('functions.R')


################################################################################
########## LOAD DATA
################################################################################
rawData <- Read10X(paste0(sampleFolder,"rawData/filtered_feature_bc_matrix/"))
dim(rawData)

colnames(rawData)<-gsub('-1',paste0('-',sampleName),colnames(rawData))

################################################################################
########## QC: CELLS
################################################################################

########################################
########## Calculate metrics
########################################

##### Create object #####
sce<-SingleCellExperiment(list(counts=rawData))
dim(sce)

##### Get mitochondrial genes #####
is.mito <- grepl("^MT-", rownames(sce))
sum(is.mito)
rownames(sce)[is.mito]

##### Calculate QC metrics #####
sce <- perCellQCMetrics(sce, subsets=list(Mt=is.mito))

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr"=colnames(rawData),"orig.ident"='tmp',
                     "nGene"=sce$detected,"nUMI"=sce$sum,"percent.mito"=sce$subsets_Mt_percent, 
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1

metaData$orig.ident<-stringr::str_split(rownames(metaData),"-") %>% purrr::map_chr(~.[[2]])
table(metaData$orig.ident)

########################################
########## Get outliers
########################################
nmad_low_feature<-3
nmad_high_feature<-3

nmad_low_UMI<-3
nmad_high_UMI<-3

nmad_high_mito<-5
nmad_low_mito<-3

##### Aim: remove cells with low library sizes, low numbers of expressed features and with high mitochondrial proportions
##same as nGene in Seurat pipeline
feature.drop.low <- isOutlier(sce$detected, nmads=nmad_low_feature, type="lower", log=TRUE)
sum(feature.drop.low)

feature.drop.high <- isOutlier(sce$detected, nmads=nmad_high_feature, type="higher", log=TRUE)
sum(feature.drop.high)

feature.drop<-as.logical(feature.drop.low + feature.drop.high)
sum(feature.drop)

##same as UMI in Seurat pipeline
libsize.drop.low <- isOutlier(sce$sum, nmads=nmad_low_UMI, type="lower", log=TRUE)
sum(libsize.drop.low)

libsize.drop.high <- isOutlier(sce$sum, nmads=nmad_high_UMI, type="higher", log=TRUE)
sum(libsize.drop.high)

libsize.drop<-as.logical(libsize.drop.low+libsize.drop.high)
sum(libsize.drop)

##% mitochondrial genes
mito.drop.low <- isOutlier(sce$subsets_Mt_percent, nmads=nmad_low_mito, type="lower")
sum(mito.drop.low)

mito.drop.high <- isOutlier(sce$subsets_Mt_percent, nmads=nmad_high_mito, type="higher")
sum(mito.drop.high)

mito.drop<-as.logical(mito.drop.low+mito.drop.high)
sum(mito.drop)


##### add to metaData matrix #####
metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop


########################################
########## Create violinPlots
########################################

### Before filtering
toPlot<-metaData
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/1a_beforeFiltering.png"))
#Split per sample
drawVlnPlot_split(toPlot, fileName = paste0(sampleFolder,"results/1b_beforeFiltering_splitted.png"))

### After filtering
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/2a_afterFiltering.png"))
# Split per sample
drawVlnPlot_split(toPlot, fileName = paste0(sampleFolder,"results/2b_afterFiltering_splitted.png"))





########################################
########## Remove outliers
########################################

goodCells <- rownames(sce)[!(libsize.drop | feature.drop | mito.drop)]
length(goodCells)

### Number of cells removed
nrow(metaData)-length(goodCells)


################################################################################
########## FINALIZE QC
################################################################################

rawDataFiltered<-rawData[,goodCells]
dim(rawDataFiltered)
saveRDS(goodCells, file=paste0(sampleFolder,"Robjects/cells_afterFiltering.rds"))

### Remove some variables
rm(sce)
rm(rawData)


################################################################################
########## CREATE SEURAT OBJECT
################################################################################

seuratObj<-CreateSeuratObject(rawDataFiltered, min.cells = 3, min.features = 200)
dim(seuratObj)

########## Add metadata ##########
seuratObj@meta.data$orig.ident<-as.character(seuratObj@meta.data$orig.ident)
seuratObj@meta.data$orig.ident<-stringr::str_split(rownames(seuratObj@meta.data),"-") %>% purrr::map_chr(~.[[2]])

########## Visualize vlnPlot ##########
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

########################################
##### Normalize
########################################
seuratObj<-NormalizeData(seuratObj)

########################################
##### HVG
########################################
seuratObj<-FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(seuratObj), 10)


########################################
##### Scaling
########################################
seuratObj<-ScaleData(seuratObj)

########################################
##### PCA
########################################
seuratObj<-RunPCA(seuratObj, features = VariableFeatures(seuratObj))
ElbowPlot(seuratObj, ndims = 40)

########################################
##### Clustering
########################################
dimsToTry<-c(20,25,30,35)
resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30)
  umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, label.size = 4)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"results/3_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), 
         width = 20, height = 6)
  
}

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_",sampleName,".rds"))


##### Test gene
FeaturePlot(object = seuratObj, features = c("MKI67"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')


################################################################################
########## DOUBLET DETECTION
################################################################################

########################################
##### Doublet Finder
########################################

##### Function doubletFinder
runDoubletFinder <- function(object, minPCT = 1, maxPCT = 10, pN = 0.25, 
                             reuse.pANN = F, PCs=1:25){
  if(!require('DoubletFinder')){
    devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
  }
  require("dplyr")
  require("stringr")
  nDoubletsMin <- (ncol(object@assays$RNA@data) * (minPCT/100) ) %>% round()
  nDoubletsMax <- (ncol(object@assays$RNA@data) * (maxPCT/100) ) %>% round()
  
  if( reuse.pANN ) {
    expectedColumn <- paste0("pANN_", pN, "_")
    
    if( grepl(pattern = expectedColumn, colnames(metaData) %>% any() )) {
      
      columnToUse <- grep( pattern = expectedColumn, x = colnames(metaData), value = T )
      pKValue <- str_split_fixed(columnToUse, "_", n = 4)[,3] %>% as.numeric()
      object <- doubletFinder(object, pN, pKValue, nExp = nDoubletsMin, reuse.pANN = columnToUse)
      object <- doubletFinder(object, pN, pKValue, nExp = nDoubletsMax, reuse.pANN = columnToUse)
      
    }
    else{
      return(
        message("No pANN found for given pN in given object, please retry with reusepANN = F")
      )
    }
  }
  else{
    findpK <- paramSweep_v3(object) %>%
      summarizeSweep() %>%
      find.pK()
    maxScore <- findpK %>% pull('BCmetric') %>% which.max()
    pKValue <- findpK[maxScore, 'pK'] %>% as.character() %>% as.numeric()
    
    object <- doubletFinder_v3(seu = object, pN = pN, pK = pKValue, nExp = nDoubletsMin, reuse.pANN = F, PCs=PCs)
    object <- doubletFinder_v3(seu = object, pN = pN, pK = pKValue, nExp = nDoubletsMax, reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_", nDoubletsMin),PCs=PCs)
    
  }
  
  
  object@meta.data$DFPrediction <- object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_", nDoubletsMax)]
  object@meta.data$DFPrediction[ object@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_", nDoubletsMin)] == 'Doublet'] <- "High Confidence"
  object@meta.data$DFPrediction <- gsub('Doublet', "Low Confidence", object@meta.data$DFPrediction)
  
  return(object)
}

### Run doubletFinder
seuratObjNew <- runDoubletFinder(seuratObj)

### Create plot
seuratObjNew@meta.data$DFPrediction<-factor(seuratObjNew@meta.data$DFPrediction,
                                            levels=c('Singlet','Low Confidence','High Confidence'))
p <- DimPlot(seuratObjNew, group.by = "DFPrediction", cols = c('#C9C9C9', 'yellow2','red'), order=T)
ggsave(p, file=paste0(sampleFolder,'results/4a_doubletFinder.png'),width = 8, height = 6)

table(seuratObjNew@meta.data$DFPrediction)


########################################
##### scDblFinder
########################################
library("scDblFinder")

### create object
sce <- SingleCellExperiment(list(counts = rawData[,colnames(seuratObj)]))

### test for several doublet rates
dbr <- tibble::lst(0.10)

for (rate in dbr){
  ## Run test (doublet rate 1-10%)
  scePer <- scDblFinder(sce, samples = seuratObj@meta.data$orig.ident, dbr=rate)
  
  ## Print output
  print(paste0("Doublet rate ",rate,":"))
  print(table(scePer@colData@listData$scDblFinder.class))
  
  ## Add to seuratObj
  seuratObj@meta.data$scDblFinder<-scePer@colData[rownames(seuratObj@meta.data),'scDblFinder.class']

  ## Create plot  
  p <- DimPlot(seuratObj, group.by = 'scDblFinder', cols = c('#C9C9C9', 'red'), order=T)
  ggsave(p, file=paste0(sampleFolder,'results/4b_scdblfinder.png'),width = 8, height = 6)
}




