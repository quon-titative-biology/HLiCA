
# Module Score analysis function taken from SeuratPipe and adapted to produce violin and dot plots.
module_score_analysis <- function(
    seu, modules_group, plot_dir = NULL, reduction = "umap", max.cutoff = "q98",
    min.cutoff = NA, legend.position = "top", col_pal = NULL, dims_plot = c(1, 2),
    seed = 1, pt.size = 1.4, fig.res = 200, alpha = c(0.1, 0.9),
    pt.size.factor = 1.1, spatial_col_pal = "inferno", crop = FALSE, ...) {
  # If no modules are given, return Seurat object
  if (is.null(modules_group)) {
    return(seu)
  }

  # Extract default assay
  assay <- Seurat::DefaultAssay(object = seu)
  
  # Iterate over the group of modules
  for (mg in names(modules_group)) {
    # Extract list of modules within the group
    modules <- modules_group[[mg]]
    # For each module compute module scores and plot
    for (m in names(modules)) {
      # Extract genes that are present in data
      features <- modules[[m]][modules[[m]] %in% rownames(seu)]
      if (length(features) == 0) { next }
      seu <- SeuratPipe::compute_module_score(seu = seu, features = features, name = m,
                                  seed = seed, ...)
      
      if (!is.null(plot_dir)) {
        plot_dim <- .plot_dims(feat_len = length(features))
        png(paste0(plot_dir, "01_markers_", m, "_feature.png"), width = plot_dim$width,
            height = plot_dim$height, res = 200, units = "in")
        print(SeuratPipe::feature_plot(seu = seu, reduction = reduction,
                          features = features,
                          max.cutoff = max.cutoff, min.cutoff = min.cutoff,
                          ncol = plot_dim$ncols,
                          legend.position = legend.position, col_pal = col_pal,
                          dims_plot = dims_plot, pt.size = pt.size,
                          combine = TRUE, ...) & Seurat::NoAxes())
        dev.off()
        
        plot_dim <- .plot_dims_2(feat_len = length(features))
        png(paste0(plot_dir, "01_markers_", m, "_violin.png"), width = plot_dim$width,
            height = plot_dim$height, res = 200, units = "in")
        print(VlnPlot(seu, features = features, ncol = plot_dim$ncols, pt.size = 0, cols = c(SeuratPipe:::discrete_col_pal, "slategrey", "maroon")) & 
               ggplot2::geom_jitter(height = 0, size = 0.5, show.legend = FALSE, alpha=0.1))
        dev.off()
        
        png(paste0(plot_dir, "01_markers_", m, "_dot.png"), width = 20,
            height = length(features), res = 200, units = "in")
        print(dot_plot_gene(seu, features = features, scale.min = 0, scale.max = 100))
        dev.off()
        
        

      }
    }
    if (!is.null(plot_dir)) {
      # Module scores in one plot
      plot_dim <- .plot_dims(feat_len = length(modules))
      png(paste0(plot_dir, "02_score_", mg, "_feature.png"), width = plot_dim$width,
          height = plot_dim$height, res = 200, units = "in")
      print(SeuratPipe::feature_plot(seu = seu, reduction = reduction,
                        features = names(modules),
                        max.cutoff = max.cutoff, min.cutoff = min.cutoff, ncol = plot_dim$ncols,
                        legend.position = legend.position, col_pal = col_pal,
                        dims_plot = dims_plot, pt.size = pt.size,
                        combine = TRUE, ...) & Seurat::NoAxes())
      dev.off()
      plot_dim <- .plot_dims_2(feat_len = length(modules))
      png(paste0(plot_dir, "02_score_", mg, "_violin.png"), width = plot_dim$width,
          height = plot_dim$height, res = 200, units = "in")
      print(VlnPlot(seu, features = names(modules), ncol = plot_dim$ncols, pt.size = 0, cols = c(SeuratPipe:::discrete_col_pal, "slategrey", "maroon")) & 
             ggplot2::geom_jitter(height = 0, size = 0.5, show.legend = FALSE, alpha=0.1))
      dev.off()
      
      png(paste0(plot_dir, "02_score_", mg, "_dot.png"), width = 20,
          height = length(names(modules)), res = 200, units = "in")
      print(SeuratPipe::dot_plot(seu, features = names(modules), scale.min = 0, scale.max = 100))
      dev.off()
      
    }
  }
  return(seu)
}



# Print Module scores order = TRUE
module_score_print <- function(
    seu, modules_group, plot_dir = NULL, reduction = "umap", max.cutoff = "q98",
    min.cutoff = NA, legend.position = "top", col_pal = NULL, dims_plot = c(1, 2),
    seed = 1, pt.size = 1.4, fig.res = 200, alpha = c(0.1, 0.9),
    pt.size.factor = 1.1, order = TRUE, ...) {
  # If no modules are given, return Seurat object
  if (is.null(modules_group)) {
    return(seu)
  }
  
  # Extract default assay
  assay <- Seurat::DefaultAssay(object = seu)
  
  # Iterate over the group of modules
  for (mg in names(modules_group)) {
    # Extract list of modules within the group
    modules <- modules_group[[mg]]
    # For each module compute module scores and plot
    for (m in names(modules)) {
      # Extract genes that are present in data
      features <- modules[[m]][modules[[m]] %in% rownames(seu)]
      if (length(features) == 0) { next }
      
      if (!is.null(plot_dir)) {
        plot_dim <- .plot_dims(feat_len = length(features))
        png(paste0(plot_dir, "01_markers_", m, "_feature.png"), width = plot_dim$width,
            height = plot_dim$height, res = 200, units = "in")
        print(SeuratPipe::feature_plot(seu = seu, reduction = reduction,
                                       features = features,
                                       max.cutoff = max.cutoff, min.cutoff = min.cutoff,
                                       ncol = plot_dim$ncols,
                                       legend.position = legend.position, col_pal = col_pal,
                                       dims_plot = dims_plot, pt.size = pt.size,
                                       combine = TRUE, order = order, ...) & Seurat::NoAxes())
        dev.off()
      }
    }
    if (!is.null(plot_dir)) {
      # Module scores in one plot
      plot_dim <- .plot_dims(feat_len = length(modules))
      png(paste0(plot_dir, "02_score_", mg, "_feature.png"), width = plot_dim$width,
          height = plot_dim$height, res = 200, units = "in")
      print(SeuratPipe::feature_plot(seu = seu, reduction = reduction,
                                     features = names(modules),
                                     max.cutoff = max.cutoff, min.cutoff = min.cutoff, ncol = plot_dim$ncols,
                                     legend.position = legend.position, col_pal = col_pal,
                                     dims_plot = dims_plot, pt.size = pt.size,
                                     combine = TRUE, order = order, ...) & Seurat::NoAxes())
      dev.off()
    }
  }
  return(seu)
}


# Plot_dim function taken from Seurat Pipe (needed for module score)
.plot_dims <- function(feat_len) {
  if (!is.numeric(feat_len)) {
    message("Invalid argument:  Setting 'feat_len' for plotting to 2")
    feat_len <- 2
  } else if (feat_len == 0){
    feat_len <- 2
  }
  if (feat_len == 1) {
    width = 9.5
    height = 7
    ncols = feat_len
  } else if (feat_len <= 3 & feat_len > 1) {
    width = 7.5 * feat_len
    height = 6
    ncols = feat_len
  } else {
    width = 30
    height = 6 * ceiling(feat_len / 4)
    ncols = 4
  }
  return(list(width = width, height = height, ncols = ncols))
}


# Adapated plot dim function for the violin plots
.plot_dims_2 <- function(feat_len) {
  if (!is.numeric(feat_len)) {
    message("Invalid argument:  Setting 'feat_len' for plotting to 2")
    feat_len <- 2
  } else if (feat_len == 0){
    feat_len <- 2
  }
  if (feat_len == 1) {
    width = 9.5
    height = 7
    ncols = feat_len
  } else if (feat_len <= 3 & feat_len > 1) {
    width = 45 
    height = 12
    ncols = 3
  } else {
    width = 60
    height = 6 * ceiling(feat_len / 4)
    ncols = 4
  }
  return(list(width = width, height = height, ncols = ncols))
}


# Dot plot function taken from seuratpipe and adpated to print genes not just metadata
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