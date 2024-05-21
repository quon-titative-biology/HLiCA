
## Flatter nested Lists for Data Presentation Titles
flatten <- function(lst) {
  do.call(c, lapply(lst, function(x) if(is.list(x)) flatten(x) else list(x)))
}



## Gene set analysis
geneset_analysis <- function(
  gene_markers, out_dir, species, topn_genes = NULL, cluster_pval_adj = 0.05,
  quick_analysis = FALSE, ontologies = c("CC", "MF", "BP"), gs_pval_adj = 0.05, gs_qval = 0.1) {
  
  # Which species do we perform geneset analysis
  if (species == "human") {
    orgdb <- "org.Hs.eg.db"
    keggdb <- "hsa"
    pathwaydb <- "human"
  } else if (species == "mouse") {
    orgdb <- "org.Mm.eg.db"
    keggdb <- "mmu"
    pathwaydb <- "mouse"
  } else {
    stop("Unrecognised 'species' argument.")
  }
  
  # Search higher levels of Ontology only in quick_analysis mode
  if (quick_analysis) {
    levs <- seq(4, 6)
  } else {
    levs <- seq(1, 6)
  }
  
  if (is.null(ontologies)) {
    ontologies <- c("CC", "MF", "BP")
  }
  
  if (is.character(gene_markers)) {
    # Load marker genes from each cluster, e.g. output of Seurat::FindAllMarkers
    markers <- read.csv(gene_markers, head = TRUE) |>
      dplyr::filter(p_val_adj < cluster_pval_adj)
    # In case we want to retain topn marker genes per cluster
    if (!is.null(topn_genes)) {
      markers <- markers |> dplyr::group_by(cluster) |>
        dplyr::slice_max(n = topn_genes, order_by = avg_log2FC)
    }
  } else {
    markers <- gene_markers
  }
  
  # Create a list for each geneset
  geneset <- list()
  for (i in unique(markers$cluster)) {
    geneset[[paste0("Cluster", i)]] <- markers[markers$cluster == i, ] |> tibble::remove_rownames()
  }
  
  # Add entrez ID column
  for (m in names(geneset)) {
    tmp <- clusterProfiler::bitr(geneset[[m]]$gene, fromType = "SYMBOL", toType = "ENTREZID",
                                 OrgDb = orgdb, drop = FALSE) |>
      dplyr::mutate(order = 1:n()) |>
      dplyr::group_by(SYMBOL) |>
      dplyr::slice_head(n = 1) |> dplyr::ungroup() |>
      dplyr::arrange(order)
    geneset[[m]]$entrez <- tmp$ENTREZID
    geneset[[m]] <- na.omit(geneset[[m]]) %>% tibble::remove_rownames()
  }
  
  ### Gene Ontology Analyses --------------------------------------------------------------------
  # MF = Molecular Function, BP = Biological Process, CC = Cellular Compartment
  
  ##
  # 1. Check for overlaps between each gene and the gene set associated
  #    with each GO term at a given semantic level.
  ##
  if (!dir.exists(paste0(out_dir, "/go_membership"))) {
    dir.create(paste0(out_dir, "/go_membership"), recursive = TRUE)
  }
  if (!quick_analysis) {
    lapply(X = names(geneset), FUN = function(set) {
      level = 3
      while (level < 7) {
        for (d in ontologies) {
          go_group <- clusterProfiler::groupGO(
            gene = geneset[[set]]$gene, keyType = "SYMBOL", OrgDb = orgdb,
            ont = d, level = level, readable = FALSE)
          go_group <- dplyr::filter(go_group@result, Count >= 1)
          go_group <- go_group[order(-go_group$Count), ]
          write.csv(go_group, paste0(out_dir, "/go_membership/",
                                     set, "_GO_", d, "_group_level", level, ".csv"))
        }
        level = level + 1
      }
    })
  }
  
  ##
  # 2. GO term enrichment analysis: identifies statistical enrichment
  #    of genes belonging to GO terms.
  ##
  out_dir_enr <- paste0(out_dir, "/go_enrichment/")
  if (!dir.exists(paste0(out_dir_enr, "/outs/"))) {
    dir.create(paste0(out_dir_enr, "/outs/"), recursive = TRUE)
  }
  for (set in names(geneset)) {
    genes_fc <- geneset[[set]]$avg_log2FC
    names(genes_fc) <- geneset[[set]]$gene
    
    for (d in ontologies) {
      go_enr <- clusterProfiler::enrichGO(
        gene = geneset[[set]]$gene, OrgDb = orgdb, keyType = "SYMBOL", ont = d,
        pAdjustMethod = "BH", pvalueCutoff = gs_pval_adj, qvalueCutoff = gs_qval)
      
      for (l in levs) {
        go_enr_level <- clusterProfiler::gofilter(go_enr, level = l)
        if (length(go_enr_level@result$geneID) > 0) {
          write.csv(go_enr_level@result,
                    paste0(out_dir_enr, "/outs/", set, "_GO_", d, "_enrichment_level", l, ".csv"))
          
          # filter for terms which are statistical significant and have more than 2 genes
          tmp <- dplyr::filter(go_enr_level, Count > 2)
          if (length(tmp@result$geneID) > 1 && length(tmp@result$p.adjust) > 1 &&
              tmp@result$p.adjust[[1]] < gs_pval_adj) {
            png(paste0(out_dir_enr, "/", set, "_GO_", d, "_level", l, "_ConceptNetwork.png"),
                width = 10, height = 7, res = 300, units = "in")
            p1 <- enrichplot::cnetplot(tmp, foldChange = genes_fc, colorEdge = TRUE) +
              ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
              ggplot2::ggtitle(paste0("GO Enrichment ", d, ":semantic level ", l))
            print(p1)
            dev.off()
            
            png(paste0(out_dir_enr, "/", set, "_GO_", d, "_level", l, "_DotPlot.png"),
                width = 13, height = 13, res = 300, units = "in")
            p2 <- enrichplot::dotplot(tmp, showCategory = 20) +
              ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
              ggplot2::theme_classic() +
              ggplot2::ggtitle(paste0("GO Enrichment ", d, ":semantic level ", l))
            print(p2)
            dev.off()
            
            png(paste0(out_dir_enr, "/", set, "_GO_", d, "_level", l, "_Heatmap.png"),
                width = 15, height = 8, res = 300, units = "in")
            p3 <- enrichplot::heatplot(tmp, foldChange = genes_fc) +
              ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
              ggplot2::ggtitle(paste0("GO Enrichment ", d, ":semantic level ", l))
            print(p3)
            dev.off()
          }
        }
      }
    }
  }
  
  
  ### Cluster Comparison ----------------------------------------------------------------------------
  if (!dir.exists(paste0(out_dir, "/cluster_comparison"))) {
    dir.create(paste0(out_dir, "/cluster_comparison"), recursive = TRUE)
  }
  # Keep only a list of genes (entrez IDs) for each cluster
  genes <- list()
  for (set in names(geneset)) {
    genes[[set]] <- geneset[[set]]$entrez
  }
  
  # Compute enriched genesets
  for (d in ontologies) {
    ego <- clusterProfiler::compareCluster(
      geneCluster = genes, fun = "enrichGO", OrgDb = orgdb, ont = d,
      pvalueCutoff = gs_pval_adj, qvalueCutoff = gs_qval)
    write.csv(ego@compareClusterResult, paste0(out_dir, "/cluster_comparison/enrichGO_", d, ".csv" ))
    
    png(paste0(out_dir, "/cluster_comparison/enrichGO_", d, ".png"),
        width = 20, height = 25, res = 300, units = "in")
    p <- enrichplot::dotplot(ego) +
      ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
      ggplot2::theme_classic()
    print(p)
    dev.off()
  }
  ekegg <- NULL
  try(ekegg <- clusterProfiler::compareCluster(
    geneCluster = genes, fun = "enrichKEGG", organism = keggdb,
    pvalueCutoff = gs_pval_adj, qvalueCutoff = gs_qval))
  if(!is.null(ekegg)){
    write.csv(ekegg@compareClusterResult, paste0(out_dir, "/cluster_comparison/enrichKEGG.csv" ))
    png(paste0(out_dir, "/cluster_comparison/enrichKEGG.png"),
        width = 20, height = 25, res = 300, units = "in")
    p <- enrichplot::dotplot(ekegg) +
      ggplot2::scale_colour_distiller(palette = "RdYlBu", direction = 1) +
      ggplot2::theme_classic()
    print(p)
    dev.off()
  }
}



##' wraping long string to multiple lines
##'
##'
##' @title str_wrap
##' @param string input string
##' @param width the maximum number of characters before wrapping to a new line
##' @return update strings with new line character inserted
##' @export
##' @author Guangchuang Yu and Erqiang Hu
str_wrap <- function(string, width = getOption("width")) {
  result <- vapply(string,
                   FUN = function(st) {
                     words <- list()
                     i <- 1
                     while(nchar(st) > width) {
                       if (length(grep(" ", st)) == 0) break
                       y <- gregexpr(' ', st)[[1]]
                       n <- nchar(st)
                       y <- c(y,n)
                       idx <- which(y < width)
                       # When the length of first word > width
                       if (length(idx) == 0) idx <- 1
                       # Split the string into two pieces
                       # The length of first piece is small than width
                       words[[i]] <- substring(st, 1, y[idx[length(idx)]] - 1)
                       st <- substring(st, y[idx[length(idx)]] + 1, n)
                       i <- i + 1
                     }
                     words[[i]] <- st
                     paste0(unlist(words), collapse="\n")
                   },
                   FUN.VALUE = character(1)
  )
  names(result) <- NULL
  result
}
#' default_labeller
#'
#' default labeling function that uses the
#' internal string wrapping function `yulab.utils::str_wrap`
#'
#' @importFrom yulab.utils str_wrap
default_labeller <- function(n) {
  function(str){
    str <- gsub("_", " ", str)
    str_wrap(str, n)
  }
}



plot_enrichment <- function(df, terms = NULL, pval_thresh = 0.05, max_terms = 25,
                            text_size = 1.0, dot_size = 5.0, label_format = 30) {
  
  # Split long terms to multiple lines
  label_func <- default_labeller(label_format)
  
  # Sanity checks
  stopifnot(is.numeric(pval_thresh))
  # Filter out pathways/terms
  df <- df %>% dplyr::filter(p.adjust <= pval_thresh)
  if (NROW(df) == 0) {
    cat("No siginificant terms at the specified 'pval_thresh' threshold")
    return(ggplot())
  }
  
  df <- df %>% dplyr::slice_min(n = max_terms, order_by = p.adjust) %>%
    dplyr::mutate(logp = -log10(p.adjust + 1e-100)) %>%
    dplyr::mutate(start = 0) %>%
    dplyr::arrange(p.adjust)
  
  # Order according to Significance
  df$Description <- factor(df$Description <- df$Description,
                           levels = df$Description[order(df$p.adjust, decreasing = TRUE)])
  
  if (!is.null(terms)) {
    df <- df |> dplyr::filter(Description %in% terms)
  }
  
  p <- ggplot(df, aes_string(x = "Description", y = "logp")) +
    geom_point(size = dot_size) +
    geom_hline(yintercept = -log10(pval_thresh), linetype="longdash") +
    geom_segment(aes_string(xend = "Description", yend = "start")) +
    ylab("-log10 adjusted p-value") +
    scale_x_discrete(labels = label_func) +
    coord_flip() +
    theme(
      axis.text.y = element_text(size = rel(text_size), hjust = 1, color = 'black'),
      axis.text.x = element_text(size = rel(1.2), vjust = 0.5, color = 'black'),
      axis.title.y = element_blank(),
      legend.position = 'none',
      panel.background = element_blank()
    )
  return(p)
}



dot_plot_gene <-  function(seu, features, group.by = NULL, labels = NULL,
                           xlab = "Signature", ylab = "Cluster",
                           legend.position = "right", col_pal = NULL, ...) {
  dot_params <- rlang::list2(...)
  params_seu <- dot_params[which(names(dot_params) %in%
                                   methods::formalArgs(Seurat::DotPlot))]
  params_gg <- dot_params[which(names(dot_params) %in%
                                  c(methods::formalArgs(ggplot2::scale_colour_distiller),
                                    "limits"))]
  params_seu <- params_seu[which(!(names(params_seu) %in% c("features", "group.by") ) )]
  params_gg <- params_gg[which(!(names(params_gg) %in% c("name", "type", "palette") ) )]
  
  
  # If we want to have different labels plotted
  if (!is.null(labels)) {
    labels <- labels[idx]
    if (length(labels) != length(features)) { labels <- features }
  } else {
    labels <- features
  }
  if (is.null(col_pal)) { col_pal = "RdYlBu" }
  
  p <- eval(rlang::expr(Seurat::DotPlot(seu, features = features,
                                        group.by = group.by, !!!params_seu))) +
    ggplot2::coord_flip() +
    eval(rlang::expr(ggplot2::scale_colour_distiller(name = NULL, type = "div",
                                                     palette = col_pal, !!!params_gg))) +
    ggplot2::scale_y_discrete(name = ylab) +
    ggplot2::scale_x_discrete(name = xlab, labels = labels) +
    ggplot2::theme(legend.position = legend.position)
  return(p)
}





percent_contribution_hist <- function(metadata = seu@meta.data, value = "orig.ident", clusters = "seurat_clusters", col_pal = opts$discrete_col_pal){
  cl <- metadata |> dplyr::group_by(!!rlang::sym(clusters), !!rlang::sym(value)) |>
    dplyr::summarise(n = dplyr::n()) |> dplyr::mutate(freq = n / sum(n))
  return(ggplot2::ggplot(data = cl,
                       ggplot2::aes(x = !!rlang::sym(clusters), y = freq, fill =!!rlang::sym(value))) +
         ggplot2::geom_bar(stat = "identity", color = "black") +
         ggplot2::theme_classic() +
         ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)) +
         ggplot2::xlab(NULL) + ggplot2::ylab("Cluster Contribution Percentage") + ggplot2::ggtitle(NULL) +
         ggplot2::scale_fill_manual(values = col_pal))
}

relative_contribution_hist <- function(metadata = seu@meta.data, value = "orig.ident", clusters = "seurat_clusters", col_pal = opts$discrete_col_pal){
  cl <- metadata |> dplyr::group_by(!!rlang::sym(value), !!rlang::sym(clusters)) |>
    dplyr::summarise(n = dplyr::n()) |> dplyr::mutate(freq = n / sum(n))
  return(ggplot2::ggplot(data = cl,
                       ggplot2::aes(x = !!rlang::sym(clusters), y = freq, fill =!!rlang::sym(value))) +
         ggplot2::geom_bar(stat = "identity", color = "black", position = ggplot2::position_dodge()) +
         ggplot2::theme_classic() +
         ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)) +
         ggplot2::xlab(NULL) + ggplot2::ylab("Relative Contribution Percentage") + ggplot2::ggtitle(NULL) +
         ggplot2::scale_fill_manual(values = col_pal))
}



presentation_create <- function(template = "~/project/Functions/Data_Presentation_template_v3.Rmd",
                                title = "Data Presenation",
                                obj_path = NULL,
                                results_path = NULL,
                                final_cluster_param = "seurat_clusters",
                                output_file = NULL){
  if (is.null(obj_path)){
    obj_path <- paste0(getwd(), "/", io$out_dir,"/",grep(paste0("*rds"), list.files(io$out_dir), value = TRUE))
  }
  if (is.null(results_path)){
    results_path = paste0(getwd(), "/",io$out_dir)
  }
  if (is.null(output_file)){
    output_file <- paste0(getwd(),"/Data_Presenation.html")
  }
  rmarkdown::render(template, params = list(title = title, obj_path = obj_path, results_path = results_path, final_cluster_param = final_cluster_param), output_file = output_file)
}


# Code written by Zoe Clarke from Gary Bader's group to threshold cells based on Malat1 expression. A gene she has found to be assocated with poor cell quality.
# Was found by correlating gene expression to dropQC results.

# This script filters cells based on MALAT1 expression
# Code is written assuming the use of a Seurat object, but
# should be able to be applied to any single-cell object


# Choose rough x value that surrounds the minimum between 0 and where slope starts to rise again; this is max_counts
# Counts is vector of MALAT1 counts, min and max are range of count values to look between
# Raise or lower bw  (e.g. to 0.1 or 0.001) if minimum doesn't look totally accurate
# Change lwd to change thickness of red line
# Change breaks to change bins of histogram
define_malat1_threshold <- function(counts, min_counts = 0, bw = 0.01, lwd = 2, breaks = 100, max_counts) {
  # Visualise MALAT1 histogram
  hist(counts, breaks = breaks, freq = FALSE)
  # Calculate the density values
  density_data <- density(counts, bw = bw)
  # Isolate x value at minimum
  range <- density_data$x >= min_counts & density_data$x <= max_counts
  threshold <- density_data$x[range][which(density_data$y[range] == min(density_data$y[range]))]
  threshold <- threshold[!is.na(threshold)]
  if (length(threshold) > 1){
    threshold <- threshold[length(threshold)]
  }
  # Visualise on histogram
  abline(v = threshold, col = "red", lwd = lwd)
  return(threshold)
}


## Example of how to use the function.
## Look at histogram of normalized MALAT1 expression
# hist(sobj@assays$RNA@data["MALAT1"], freq = FALSE, breaks = 100)

## If there is a peak (even a small peak) at zero, followed by a dip, then a more
## normal distribution, run the following code.


## Run this function on MALAT1 reads, making sure to include max_counts eg:
# threshold <- define_malat1_threshold(sobj@assays$RNA@data["MALAT1",])

## Use this value to subset your cells
# cells_subset <- WhichCells(sobj, expression = MALAT1 > threshold, slot = "data")
# sobj_subset <- subset(sobj, cells = cells_subset)

## OR assign TRUE/FALSE values to determine which cells would pass the filter
# malat1_threshold <- sobj@assays$RNA@data['MALAT1',] > threshold
# sobj$malat1_threshold <- malat1_threshold
# DimPlot(sobj, group.by = "malat1_threshold")
