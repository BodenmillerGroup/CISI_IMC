# Load libraries
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggsci)
library(metan)
library(circlize)
library(reshape2)
library(ComplexHeatmap)
library(stringr)
library(gridGraphics)
library(zellkonverter)
library(SingleCellExperiment)
library(cytomapper)
library(jaccard)


## General helper fnc.

# Define fnc to read in all modules and save into dataframe
read_U <- function(file, type){
  # Depending which test run different helper columns are added
  if(type=="size"){
    k_name <- gsub("k_", "", str_split(file, "/")[[1]][8])
    rep_name <- str_split(file, "/")[[1]][10]
    size_name <- str_split(file, "/")[[1]][9]
  } else if(type=="par"){
    k_name <- str_split(file, "/")[[1]][9]
    rep_name <- str_split(file, "/")[[1]][10]
  } else if(type=="training") {
    k_name <- gsub("k_", "", str_split(file, "/")[[1]][8])
    rep_name <- str_split(file, "/")[[1]][5]
  }
  
  temp.df <- read.csv(file)
  if(ncol(temp.df)!=0){
    temp.df <- temp.df %>%
      dplyr::rename(protein=X) %>%
      dplyr::rename_with(~gsub("X", "", .x, fixed=TRUE)) %>%
      melt(id="protein") %>%
      mutate(rep=rep_name)
    
    if(type=="par" | type=="training"){
      temp.df <- temp.df %>%
        mutate(k=k_name)
    } else if (type=="size"){
      temp.df <- temp.df %>%
        mutate(k=k_name,
               size=size_name)
    }
    
    temp.df
  } else {
    data.frame()
  }
}


# Function that reads in results .csv file and adjust its just that all results
# file for the comparison of lung and tonsil can be put into one dataframe
read_results <- function(file, type, voi="k"){
  if(grepl("training", file)){
    dataset_name <- str_split(file, "/")[[1]][5]
    training_name <- dataset_name
    k_name <- gsub("*._", "", str_split(file, "/")[[1]][8])
    
    if(type=="res"){
      datasize_name <- gsub("*._", "", str_split(file, "/")[[1]][7])
    }
  } else {
    dataset_name <- gsub(".*_", "", str_split(file, "/")[[1]][6])
    training_name <- stringi::stri_replace_all_regex(str_split(file, "/")[[1]][5],
                                                     pattern=c(dataset_name, "_vs_"),
                                                     replacement=c("", ""), vectorize=FALSE)
    k_name <- gsub("*._", "", str_split(file, "/")[[1]][7])
    if(type=="res"){
      datasize_name <- NA
    }
  }
  
  if(type=="res"){
    res <- read_tsv(file, show_col_types=FALSE) %>%
      mutate(dataset=dataset_name,
             training=training_name,
             simulation=ifelse(grepl("composite", file), "composite", "noisy"),
             !!as.symbol(voi):=k_name,
             datasize=datasize_name)
    res[which(res=="none")] <- NA
    
  } else if(type=="x") {
    res <- readH5AD(file)
    metadata(res) <- list(dataset=dataset_name, training=training_name,  
                          ground_truth=ifelse(grepl("X_test", file), "ground_truth", "simulated"))
    metadata(res)[[voi]] <- k_name
  }
  
  res
}


# Function that binarizes U by calculating which proteins explain up to thresh
# much of the weighted sum of squared entries in one module
compute_module_membership <- function(m.col, thresh=0.9){
  cs <- cumsum(-sort(-m.col^2))
  cs <- cs / cs[length(cs)]
  cs <- ifelse(cs<=thresh, 1, 0)
  cs[which(cs==0)[1]] <- 1
  cs <- cs[order(match(names(cs), names(m.col)))]
  cs
}


# Function to compare correlation matrices of proteins in U (all plots using
# this fnc. are not in the final output .html, since the mantel test is used
# instead)
compare_cor <- function(A, B){
  sum(abs(A - B))
}


# Compute confusion matrix using pre-trained RF classifier for a list of SCE
compute_confusionMatrix <- function(sce, rffit){
  cur_mat <- t(assay(sce, "exprs"))
  cur_mat <- cbind(as.data.frame(colData(sce)) %>%
                     dplyr::select(cell_id, subbatch) %>%
                     mutate(dummy=1,
                            subbatch=paste0("subbatch", subbatch)) %>%
                     dcast(cell_id~subbatch, value.var="dummy", fill=0) %>%
                     dplyr::select(-cell_id))
                     
  # Predict cell phenotypes in test data
  cur_pred <- predict(rffit, 
                      newdata=cur_mat)
  # cur_pred <- as.character(predict.train(rffit.lung,
  #                                        newdata=cur_mat,
  #                                        type="raw"))
  cm <- confusionMatrix(data=cur_pred,
                        reference=factor(sce$cell_labels),
                        mode="everything")
  cm
}


## Plotting fncs.

# Plots U for different iter_var as heatmaps with multiple heatmaps per iter_var
# next to each other according to repetition variable
plot_U <- function(df, iter_var, repetition){
  # Empty list to save correlation matrices
  df.cor <- list()
  
  for (i in unique(df[[iter_var]])) {
    cat('#####', iter_var, ' = ', i, '\n')
    temp_list <-  list()
    
    # Create empty HeatmapList and add all repetitions as individual heatmaps
    ht_list = NULL
    names.repetition <- c()
    for(r in unique(df[[repetition]])) {
      res.temp <- df %>%
        dplyr::filter(!!as.symbol(repetition)==r & !!as.symbol(iter_var)==i) %>%
        dplyr::select(-c(repetition, iter_var)) %>%
        pivot_wider(names_from = variable, values_from = value) %>%
        column_to_rownames(var="protein") %>%
        as.matrix()
      
      if (nrow(res.temp)!=0){
        # Calculate correlation matrix
        temp_list[[length(temp_list)+1]] <- cor(t(res.temp))
        
        res.temp <- apply(res.temp, 2, compute_module_membership)
        
        ht_list = ht_list + Heatmap(res.temp, column_title = paste0(str_to_title(repetition), 
                                                                    " ", r, 
                                                                    "\n(dictionary size: ", 
                                                                    ncol(res.temp), ")"),
                                    col=structure(c("grey", pal_npg("nrc")("1")[1]), 
                                                  names=c("0", "1")), 
                                    show_heatmap_legend=FALSE, 
                                    show_row_dend=FALSE, row_names_gp=gpar(fontsize = 8),
                                    column_names_gp=gpar(fontsize=8),
                                    column_title_gp=gpar(fontsize=10),
                                    rect_gp=gpar(col="white", lwd=1))
        names.repetition <- c(names.repetition, r)
      }
    }
    
    # Draw all heatmaps together
    draw(ht_list, heatmap_legend_list=list(Legend(labels=c("not active", "active"), 
                                                  title="Activity",
                                                  legend_gp=gpar(fill=c("grey", 
                                                                        pal_npg("nrc")("1")[1])),
                                                  at=c("0", "1"))),
         column_title=paste0(iter_var, " = ", i))
    
    names(temp_list) <- names.repetition
    df.cor[[length(df.cor)+1]] <- temp_list
    cat('\n\n')
  }
  
  names(df.cor) <- unique(df[[iter_var]])
  df.cor  
}


# Plots U for different iter_var as heatmaps where two U's are concatenated and 
# coloured by which dataset the module belongs to
# Shared modules are only shown once and coloured accordingly
plot_U_membership <- function(df, iter_var, repetition){
  
  for (i in unique(df[[iter_var]])) {
    cat('#####', iter_var, ' = ', i, '\n')
    temp_list <-  list()
    
    for(r in unique(df[[repetition]])) {
      res.temp <- df %>%
        dplyr::filter(!!as.symbol(repetition)==r & !!as.symbol(iter_var)==i) %>%
        dplyr::select(-c(repetition, iter_var)) %>%
        pivot_wider(names_from = variable, values_from = value) %>%
        column_to_rownames(var="protein") %>%
        as.matrix()
      
      temp_list[length(temp_list)+1] <- list(apply(res.temp, 2, compute_module_membership))
    }
    names(temp_list) <- unique(df[[repetition]])
    
    modules.both <- lapply(1:ncol(temp_list[[1]]), function(i){
      unlist(lapply(1:ncol(temp_list[[2]]), function(j){
        if(all(temp_list[[1]][, i]==temp_list[[2]][, j])){
          data.frame(x=c(i), y=c(j))
        }
      }))
    }) %>% bind_rows()
    
    modules.matrix <- temp_list[[1]][, modules.both$x] %>%
      bind_cols(2 * temp_list[[1]][, -(modules.both$x)]) %>%
      bind_cols(3 * temp_list[[2]][, -(modules.both$y)]) %>%
      as.matrix()
    colnames(modules.matrix) <- NULL
    rownames(modules.matrix) <- rownames(temp_list[[1]])
    
    
    modules.heatmap <- Heatmap(modules.matrix, 
                               col=structure(c("grey", pal_npg("nrc")("3")[1:3]), 
                                             names=c("0", "1", "2", "3")), 
                               show_heatmap_legend=TRUE, 
                               show_row_dend=FALSE, row_names_gp=gpar(fontsize = 8),
                               column_names_gp=gpar(fontsize=8),
                               column_title_gp=gpar(fontsize=10),
                               rect_gp=gpar(col="white", lwd=1),
                               clustering_distance_columns = function(x, y) 1 - jaccard(ifelse(x>0, 1, 0), 
                                                                                        ifelse(y>0, 1, 0)),
                               clustering_distance_rows=function(x, y) 1 - jaccard(ifelse(x>0, 1, 0),
                                                                                   ifelse(y>0, 1, 0)),
                               heatmap_legend_param=list(at=c("0", "1", "2", "3"),
                                                         labels=c("not active", "active in both",
                                                                  paste0("active in ", 
                                                                         unique(df[[repetition]]))),
                                                         title="Proteins",
                                                         title_gp=gpar(fontsize=10, 
                                                                       fontface="bold"),
                                                         labels_gp=gpar(fontsize=8)))
    print(modules.heatmap)
    cat('\n\n')
  }
}


# Compute pairwise compare_cor between U matrices and plot as boxplots
plot_dist_boxplot <- function(df, x_lab, to_int=TRUE){
  # Compute all pairwise distances between correlation matrices of repetions for
  # each k
  df.dist <- data.frame((sapply((lapply(df, function(l){
    temp_comp <- sapply(l, function(x) sapply(l, function(y) compare_cor(x, y)))
    temp_comp[upper.tri(temp_comp)]
  })), c))) %>%
    dplyr::rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    melt() %>%
    dplyr::rename(Distance=value)
  
  if(to_int){
    df.dist <- df.dist %>%
      mutate(variable=factor(df.dist$variable, levels=sort(strtoi(unique(df.dist$variable)))))
  } 
  
  df.boxplot <- ggplot(df.dist, aes(x=variable, y=Distance, fill=variable)) +
    geom_boxplot() +
    scale_fill_npg() +
    theme_cowplot() +
    theme(legend.position = "none") +
    labs(x=x_lab) 
  
  df.boxplot
}


# Compute mantel test per var_name and use metans pairs_mantel to plot
plot_metan_mantel <- function(df, var_name){
  for (n in names(df)) {
    cat('##### ', var_name, ' = ', n, '\n')
    mantel_plot <- as.lpcor(df[[n]][[1]],
                            df[[n]][[2]],
                            df[[n]][[3]],
                            df[[n]][[4]],
                            df[[n]][[5]],
                            df[[n]][[6]],
                            df[[n]][[7]],
                            df[[n]][[8]],
                            df[[n]][[9]],
                            df[[n]][[10]]) %>%
      pairs_mantel(diag = TRUE,
                   pan.spacing = 0,
                   shape.point = 21,
                   col.point = "black",
                   fill.point = "red",
                   size.point = 1.5,
                   alpha.point = 0.6,
                   alpha = 0.2)
    print(mantel_plot)
    cat('\n\n')
  }
}


# Compute pairwise mantel test between U matrices and plot results as boxplots
plot_mantel_boxplot <- function(df, var_name, to_int=TRUE){
  df.mantel <- data.frame((sapply((lapply(df, function(l){
    temp_comp <- sapply(l, function(x) sapply(l, function(y) mantel_test(x, y)$mantel_r))
    temp_comp[upper.tri(temp_comp)]
  })), c))) %>%
    dplyr::rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    melt() %>%
    dplyr::rename(mantel=value)
  
  if(to_int){
    df.mantel <- df.mantel %>%
      mutate(variable=factor(df.mantel$variable, levels=sort(strtoi(unique(df.mantel$variable)))))
  }
  
  df.mantelBoxplot <- ggplot(df.mantel, aes(x=variable, y=mantel, fill=variable)) +
    geom_boxplot() +
    scale_fill_npg() +
    theme_cowplot() +
    guides(fill=FALSE) +
    labs(x=var_name, y=expression(cor["Mantel"])) 
  
  df.mantelBoxplot  
}


# Plot results of CISI vs ground truth for protein of interest poi and show for
# cell segmentation
plot_cells <- function(sce.list, masks.list, poi, layer="none"){
  
  colour.cells <- list(el=c("#FFFFFF", pal_npg("nrc")("8")[8]))
  names(colour.cells) <- poi
  
  img.list <- list()
  idx <- c()
  for(el in sce.list){
    # Select if normal counts or transformed counts should be plotted
    if(layer=="none"){
      p <- plotCells(mask=masks.list, object=el,
                     cell_id="ObjectNumber", img_id="sample_id", colour_by=poi,
                     return_plot=TRUE,  image_title=list(cex=1),
                     colour=colour.cells, display="single",
                     scale_bar=list(cex=1, lwidth=5),
                     legend = list(colour_by.title.cex=0.7, colour_by.labels.cex=0.7))
    } else {
      p <- plotCells(mask=masks.list, object=el,
                     cell_id="ObjectNumber", img_id="sample_id", colour_by=poi,
                     return_plot=TRUE,  image_title=list(cex=1),
                     colour=colour.cells, display="single",
                     scale_bar=list(cex=1, lwidth=5),
                     legend = list(colour_by.title.cex=0.7, colour_by.labels.cex=0.7),
                     exprs_values=layer)
    }
    p.gg <- lapply(p$plot, function(x){ggdraw(x, clip="on") })
    img.list <- append(img.list, p.gg)
    idx <- c(idx, seq_along(p.gg))
  }
  
  names(img.list) <- make.unique(names(img.list))
  idx <- order(idx)
  img.list[idx]
}


# Function that plots exprs (asinh counts) of protein_x vs protein_y for both
# ground truth and decomposed data coloured by celltype
plot_exprs <- function(sce.list, celltype_col, protein_x, protein_y, layer="exprs"){
  # Remove special characters in protein names to make plotting of such proteins
  # possible
  sce.list <- lapply(sce.list, function(sce){
    rownames(sce) <- gsub("/|-", "_", rownames(sce))
    sce
    })
  protein_x <- gsub("/|-", "_", protein_x)
  protein_y <- gsub("/|-", "_", protein_y)

  if ("celltype" %in% names(colData(sce.list[[2]]))){
    col_cells <- rep("grey", length(unique(colData(sce.list[[1]])$celltype)))
    names(col_cells) <- unique(colData(sce.list[[1]])$celltype)
    if(celltype_col!=""){
      col_cells[celltype_col] <- "red"
    }
    
    alpha_cells <- rep(0.2, length(unique(colData(sce.list[[1]])$celltype)))
    names(alpha_cells) <- unique(colData(sce.list[[1]])$celltype)
    alpha_cells[celltype_col] <- 1
    
    plot.sim <- ggcells(sce.list[[1]], 
                        aes(x=!!as.symbol(protein_x), y=!!as.symbol(protein_y), 
                            colour=celltype, alpha=celltype), 
                        exprs_values=layer) +
      geom_point(size=0.3) +
      theme_cowplot(10) +
      scale_color_manual(values=col_cells) +
      guides(color=FALSE) +
      scale_alpha_manual(guide='none', values=alpha_cells)
    
    plot.true <- ggcells(sce.list[[2]], 
                         aes(x=!!as.symbol(protein_x), y=!!as.symbol(protein_y), 
                             colour=celltype, alpha=celltype), 
                         exprs_values=layer) +
      geom_point(size=0.4) +
      theme_cowplot(10) +
      scale_color_manual(values=col_cells) +
      scale_alpha_manual(guide='none', values=alpha_cells)
  } else {
    plot.sim <- ggcells(sce.list[[1]], 
                        aes(x=!!as.symbol(protein_x), y=!!as.symbol(protein_y)), 
                        exprs_values=layer) +
      geom_point(size=0.3) +
      theme_cowplot(10) +
      scale_color_manual(values=col_cells) +
      guides(color=FALSE) +
      scale_alpha_manual(guide='none', values=alpha_cells)
    
    plot.true <- ggcells(sce.list[[2]], 
                         aes(x=!!as.symbol(protein_x), y=!!as.symbol(protein_y)), 
                         exprs_values=layer) +
      geom_point(size=0.4) +
      theme_cowplot(10) +
      scale_color_manual(values=col_cells) +
      scale_alpha_manual(guide='none', values=alpha_cells)
  }
  
  plot_grid(plot.sim, plot.true, 
            labels=names(sce.list), 
            label_size=15, hjust=c(-2, -1.5), vjust=1)
}