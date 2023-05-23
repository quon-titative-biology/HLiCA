##### Function drawVlnPlot
drawVlnPlot<-function(toPlot, fileName, width=12, height=8){
  toPlot<-toPlot[order(toPlot[,'nGene.drop']),]
  p_nGene <- ggplot(toPlot, aes(staticNr, nGene)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic() +
    theme(legend.position = "none")
  
  toPlot<-toPlot[order(toPlot[,'nUMI.drop']),]
  p_nUMI <- ggplot(toPlot, aes(staticNr, nUMI)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic() +
    theme(legend.position = "none")
  
  toPlot<-toPlot[order(toPlot[,'mito.drop']),]
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic() +
    theme(legend.position = "none")
  
  grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
  if(fileName != ""){
    ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file=fileName, dpi=200,
           width=width, height=height)
  }
}

##### Function drawVlnPlot_split
drawVlnPlot_split<-function(toPlot, fileName, width=15, height=8){
  toPlot<-toPlot[order(toPlot[,'nGene.drop']),]
  p_nGene <- ggplot(toPlot, aes(orig.ident, nGene, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop), alpha=0.5) +
    geom_violin() + 
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=14), axis.text.y = element_text(size = 14))
  
  toPlot<-toPlot[order(toPlot[,'nUMI.drop']),]
  p_nUMI <- ggplot(toPlot, aes(orig.ident, nUMI, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop), alpha=0.5) +
    geom_violin() + 
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=14), axis.text.y = element_text(size = 14))
  
  toPlot<-toPlot[order(toPlot[,'mito.drop']),]
  p_mito <- ggplot(toPlot, aes(orig.ident, percent.mito, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes(col=mito.drop), alpha=0.5) +
    geom_violin() + 
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14))
  
  thePlot<-grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
  if(fileName != ""){
    ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file=fileName, dpi=200, 
           width=width, height=height)
  }
  else{
    return(thePlot)
  }
}




