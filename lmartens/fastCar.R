# devtools::install_git("https://git.web.rug.nl/P278949/FastCAR")

library('Seurat')
library('Matrix')
library('qlcMatrix')
library('FastCAR')
library('ggplot2')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

getwd()

sampleName<-"."
sampleFolder<-paste0(sampleName,"/")

##add some subfolders
dir.create(paste0(sampleFolder,"results_fastCAR"))


########################################
##### Functions
########################################

##### Function do calculations for ambient plot
calculateForPlot<-function(fullMatrix, stopLimit, byStep){
  print('Getting ambient profile')
  ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = fullMatrix, 
                                             start = 0, 
                                             stop = stopLimit, 
                                             by = byStep, 
                                             contaminationChanceCutoff = 0.05)
  return(ambProfile)
}


##### Function make ambient plot
makePlot<-function(ambProfile, theCutOff, sampleName){
  ### Make plot
  print('Making plot')
  png(file=paste0(sampleFolder,'results_fastCAR/plotAmbientProfile_',sampleName,'.png'), width = 900, height = 800)
  par(mfrow = c(3, 1))
  plot(as.numeric(rownames(ambProfile)), ambProfile[,1], 
       main = "Total number of empty droplets at cutoffs", 
       xlab = "empty droplet UMI cutoff", ylab = "Number of empty droplets")
  abline(v=theCutOff, col="red")
  
  plot(as.numeric(rownames(ambProfile)), ambProfile[,2], 
       main = "Number of genes in ambient RNA", xlab = "empty droplet UMI cutoff", 
       ylab = "Genes in empty droplets")
  abline(v=theCutOff, col="red")
  
  plot(as.numeric(rownames(ambProfile)), ambProfile[,3], 
       main = "number of genes to correct", xlab = "empty droplet UMI cutoff", 
       ylab = "Genes identified as contamination")
  abline(v=theCutOff, col="red")
  dev.off()
}

##### Function get ambient profile
getAmbientProfile<-function(fullMatrix, cellMatrix, emptyDropletCutoff, contaminationChanceCutoff){
  ### Determine ambient RNA
  ambientProfile = determine.background.to.remove(fullMatrix, cellMatrix, emptyDropletCutoff, contaminationChanceCutoff)
  topGenes<-tail(sort(ambientProfile), 50)
  
  print(paste0(sum(ambientProfile>0),' genes using for correction'))
  
  ##Create barplot of top genes
  toPlot<-as.data.frame(topGenes)
  toPlot$gene<-rownames(toPlot)
  p<-ggplot(data=toPlot, aes(x=reorder(gene, -topGenes), y=topGenes)) +
    geom_bar(stat="identity",color='gray20',fill='steelblue') +
    coord_flip() +
    theme_classic() +
    theme(axis.text.y =  element_text(size=6))
  ggsave(p, file=paste0(sampleFolder,"results_fastCAR/barplot_topAmbientGenes_",sampleName,".png"))
  
  ##Write all genes used for correction
  corrGenes<-ambientProfile[ambientProfile>0]
  corrGenes<-corrGenes[order(corrGenes, decreasing = T)]
  write.table(as.matrix(corrGenes), 
              file=paste0(sampleFolder,'results_fastCAR/correctionGenes_',sampleName,'.txt'),
              sep='\t', col.names = F)
  
  ###Return ambientProfile
  return(ambientProfile)
}



################################################################################
########## RUN FASTCAR
################################################################################
sampleName<-'CS049'

########## Read data ##########
cellMatrix     = read.cell.matrix(paste0("rawData/",sampleName,"/filtered_feature_bc_matrix"))
dim(cellMatrix)
# 32285  2999

fullMatrix     = read.full.matrix(paste0("rawData/",sampleName,"/raw_feature_bc_matrix"))
dim(fullMatrix)
# 32285 6794880

### Filter raw folder
min(colSums(cellMatrix))
# 500

the_colSums<-colSums(fullMatrix)
neededCells<-names(the_colSums[the_colSums > 10])
fullMatrix<-fullMatrix[,neededCells]
dim(fullMatrix)
# 32285 78295

########## Make plot ##########
ambProfile<-calculateForPlot(fullMatrix, 700, 10)

theCutOff<-120
makePlot(ambProfile, theCutOff, sampleName)

### How many empty droplets
tmp<-names(the_colSums[the_colSums < theCutOff])
length(intersect(tmp, neededCells))
# 7786

########## Get ambient genes ##########
ambientProfile = getAmbientProfile(fullMatrix, cellMatrix, theCutOff, 0.05)
#76 genes

########## Remove background ##########
cellMatrix_new     = remove.background(cellMatrix, ambientProfile)
dim(cellMatrix_new)
# 32285  2999
saveRDS(cellMatrix_new, file=paste0(sampleFolder,'results_fastCAR/outputFastCAR/cellMatrixNew_',sampleName,'.rds'))

##### Remove some variables
rm(cellMatrix)
rm(cellMatrix_new)
rm(fullMatrix)






