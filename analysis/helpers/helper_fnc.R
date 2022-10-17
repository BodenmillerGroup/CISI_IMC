################################################################################
### Helper fncs.
# This file contains helper fncs. used in the analysis of CISI run on IMC data
# (Master Thesis)
################################################################################



## Load libraries

library(caret)
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
library(ggpattern)
library(scater)
library(viridis)



## Plot defaults

# Set fontsizes used throughout this script
title.fontsize <- 10
axis_title.fontsize <- 8



## General helper fnc.

# Define fnc to read in a single dictionary U from "file" into dataframe in long
# format (columns: protein, module, membership=how a specific protein is part
# of a specific module)
read_single_U <- function(file){
  # Read in file
  temp.df <- read.csv(file)
  # Catch the error of CISI failing, e.g. U being empty by returning emtpy dataframe
  if (ncol(temp.df)!=0){
    # Convert data frame to long format for easier handling later on
    temp.df <- temp.df %>%
      dplyr::rename(protein=X) %>%
      dplyr::rename_with(~gsub("X", "", .x, fixed=TRUE)) %>%
      melt(id="protein") %>%
      dplyr::rename(module=variable, membership=value)
    
    temp.df
  } else {
    # Return empty data frame if U is empty
    data.frame()
  }
}


# Define fnc to read in all modules and save into dataframe
# type: "size", "par" or "training"
read_U <- function(file, type){
  # Depending what kind of test was run (type) different helper columns are added
  # Possible columns added are: k used, which repetition iteration data is form
  # and if datasize CISI was run on
  if (type=="size"){
    k_name <- gsub("k_", "", str_split(file, "/")[[1]][8])
    rep_name <- str_split(file, "/")[[1]][10]
    size_name <- str_split(file, "/")[[1]][9]
  } else if (type=="par"){
    k_name <- str_split(file, "/")[[1]][9]
    rep_name <- str_split(file, "/")[[1]][10]
  } else if (type=="training") {
    k_name <- gsub("k_", "", str_split(file, "/")[[1]][8])
    rep_name <- str_split(file, "/")[[1]][5]
  }
  
  # Read in file
  temp.df <- read_single_U(file)
  
  # Catch the error of CISI failing, e.g. U being empty by not including these files
  if (ncol(temp.df)!=0){
    # Convert data frame to long format to combine all files into one big data frame 
    temp.df <- temp.df %>%
      mutate(rep=rep_name,
             k=k_name)
    
    # Add column size if type=="size"
    if (type=="size"){
      temp.df <- temp.df %>%
        mutate(size=size_name)
    }
  } 
  
  temp.df
}


# Define fnc to read in a single experiment matrix A/Phi from "file" into dataframe 
# in long format (columns: channel, one column per protein)
read_single_A <- function(file){
  # Read in A/Phi 
  a <- read.table(file, row.names=1, sep="\t") %>%
    dplyr::rename_all(~ (stringr::str_replace_all( ., "\\.", "-" ))) %>%
    rownames_to_column(var="channel")
  
  a
}


# Define fnc to read in all simulation designs (A/Phi) and save into dataframe
# voi: name of added column (folder name in which A/Phi is saved, should give information
#      on what the parameter or experiment of interest for that CISI ran was)
# use_voi: TRUE or FALSE if an extra column voi should be added
read_A <- function(file, voi="k", use_voi=TRUE){
  # Depending if the file is somwhere in the training folder or not (as is the
  # case where we looked at cross-dataset results) the extra columns are taking 
  # from different parts of the file path
  # Added columns: dataset CISI was tested on, trained on and what the voi was
  if (grepl("training", file)){
    dataset_name <- str_split(file, "/")[[1]][5]
    training_name <- dataset_name
    
    if (use_voi){
      k_name <- gsub("*._", "", str_split(file, "/")[[1]][8])
    }
  } else {
    dataset_name <- gsub(".*_", "", str_split(file, "/")[[1]][6])
    training_name <- stringi::stri_replace_all_regex(str_split(file, "/")[[1]][5],
                                                     pattern=c(dataset_name, "_vs_"),
                                                     replacement=c("", ""), vectorize=FALSE)
    if (use_voi){
      k_name <- gsub("*._", "", str_split(file, "/")[[1]][7])
    }
  }
  
  # Read in A/Phi and add a colun for the channel and add above additional columns
  a <- read_single_A(file) %>%
    mutate(dataset=dataset_name,
           training=training_name)
  
  if (use_voi){
    a <- a %>%
      mutate(!!as.symbol(voi):=k_name)
  }
  
  a
}


# Define fnc to read in a single correlation results from "file" into dataframe 
# and add column specifying if the results come from noisy or no_noise simulations
read_single_correlation_results <- function(file){
  # Read in correlation results
  res <- read_tsv(file, show_col_types=FALSE) %>%
    mutate(simulation=ifelse((grepl("composite", file) | grepl("no_noise", file)), 
                             "no_noise", "noisy"))
  
  res
}


# Define fnc to read in a single anndata object from "file" and add to the metadata
# if the object contains ground_truth or simulated expressions
read_single_anndata <- function(file){
  # Read in anndata object
  res <- readH5AD(file)
  metadata(res) <- list(ground_truth=ifelse(grepl("X_test", file), 
                                            "ground_truth", 
                                            "simulated"))
  
  res
}


# Function that reads in results .csv file and adjust its just that all results
# file for the comparison of lung and tonsil can be put into one dataframe
# type: res (results.txt files with correlations) or x (X_*.h5ad files with anndata
# object)
# voi: name of added column (folder name in which A/Phi is saved, should give information
#      on what the parameter or experiment of interest for that CISI ran was)
# use_voi: TRUE or FALSE if an extra column voi should be added
read_results <- function(file, type, voi="k", use_voi=TRUE){
  # Depending if the file is somwhere in the training folder or not (as is the
  # case where we looked at cross-dataset results) the extra columns are taking 
  # from different parts of the file path
  # Added columns: dataset CISI was tested on, trained on, what the voi was and
  # dataset size (e.g. full set of proteins or only subset of protein was used)
  if (grepl("training", file)){
    dataset_name <- str_split(file, "/")[[1]][5]
    training_name <- dataset_name
    
    if (use_voi){
      k_name <- gsub(".*_", "", str_split(file, "/")[[1]][8])
    }
    
    if (type=="res"){
      datasize_name <- str_split(file, "/")[[1]][7]
    }
  } else {
    dataset_name <- gsub(".*_", "", str_split(file, "/")[[1]][6])
    training_name <- stringi::stri_replace_all_regex(str_split(file, "/")[[1]][5],
                                                     pattern=c(dataset_name, "_vs_"),
                                                     replacement=c("", ""), vectorize=FALSE)
    if (use_voi){
      k_name <- gsub(".*_", "", str_split(file, "/")[[1]][7])
    }
    
    if (type=="res"){
      datasize_name <- NA
    }
  }
  
  if (type=="res"){
    # If type="res" we read in file containing correlation results
    # Column gets added to show if results comes from noisy or no noise simulations
    res <- read_single_correlation_results(file) %>%
      mutate(dataset=dataset_name,
             training=training_name,
             datasize=datasize_name)
    
    if (use_voi){
      res <- res %>%
        mutate(!!as.symbol(voi):=k_name)
    }
    
    res[which(res=="none")] <- NA
    
  } else if (type=="x") {
    # If type="x", then anndata object gets read into SCE containing either
    # ground truth or simulated/decomposed expressions (normalized and subseted
    # according to what was specified/used in the CISI run)
    res <- read_single_anndata(file)
    metadata(res) <- append(metadata(res),
                            list(dataset=dataset_name, training=training_name))
    
    if (use_voi){
      metadata(res)[[voi]] <- k_name
    }
  }
  
  res
}


# Function that binarizes U by calculating which proteins explain up to thresh
# much of the weighted sum of squared entries in one module
# (From original paper)
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
# !Deprecated!
compare_cor <- function(A, B){
  sum(abs(A - B))
}


# Train a random forest classifier on arcsinh transformed counts
# using the gating labels of the training data of a pretrained RF classifier
# rffit.original
train_rf <- function(sce, rffit.original){
  # Extract training data
  train_mat <- t(assay(sce, "exprs"))
  # Add subbatch information to training data and filter for cells used in the
  # pretrained classifier for training
  train_mat <- cbind(train_mat,
                   as.data.frame(colData(sce)) %>%
                     dplyr::select(cell_id, subbatch) %>%
                     mutate(dummy=1,
                            subbatch=paste0("subbatch", subbatch)) %>%
                     dcast(cell_id~subbatch, value.var="dummy", fill=0) %>%
                     dplyr::select(-cell_id)) %>%
    dplyr::filter(rownames(.) %in% rownames(rffit.original$trainingData))
  # Add gating labels from pretrained classifier
  train_mat <- train_mat[colnames(train_mat) %in% colnames(rffit.original$trainingData)] %>%
    merge(rffit.original$trainingData %>%
            dplyr::select(".outcome"), by=0, all.x=TRUE)
  # Convert Row.names column to rownames of matrix
  rownames(train_mat) <- train_mat$Row.names
  train_mat <- train_mat %>%
    dplyr::select(-Row.names)
  
  # Specify training parameters for 5-fold cross validation
  fitControl <- trainControl(method="cv",
                             number=5)
  
  # Train a random forest model for predicting cell labels
  # This call also performs parameter tuning
  rffit <- train(x=train_mat %>%
                   dplyr::select(-".outcome"), 
                 y=factor(train_mat$.outcome),
                 method="rf", ntree=1000,
                 tuneLength=5,
                 trControl=fitControl)
  
  rffit
}


# Compute confusion matrix using pre-trained RF classifier for a SCE
compute_confusionMatrix <- function(sce, rffit){
  # Extract test data to predict celltypes on
  cur_mat <- t(assay(sce, "exprs"))
  # Add subbatch information to training data and filter for cells NOT used in the
  # pretrained classifier for training or having celltype "uncertain"
  cur_mat <- cbind(cur_mat,
                   as.data.frame(colData(sce)) %>%
                     dplyr::select(cell_id, subbatch) %>%
                     mutate(dummy=1,
                            subbatch=paste0("subbatch", subbatch)) %>%
                     dcast(cell_id~subbatch, value.var="dummy", fill=0) %>%
                     dplyr::select(-cell_id)) %>%
    dplyr::filter(!(rownames(.) %in% rownames(rffit$trainingData)) &
                  !(rownames(.) %in% rownames(colData(sce)[colData(sce)$celltype=="uncertain",])))
  cur_mat <- cur_mat[colnames(cur_mat) %in% colnames(rffit$trainingData)]
                     
  # Predict cell phenotypes in test data
  cur_pred <- predict(rffit, 
                      newdata=cur_mat)

  # Extract ground truth labels
  gt_labels <- colData(sce) %>%
    as.data.frame() %>%
    dplyr::filter(!(rownames(.) %in% rownames(rffit$trainingData)) &
                 (celltype!="uncertain"))
  gt_labels <- gt_labels[match(rownames(cur_mat), rownames(gt_labels)), ] %>%
    dplyr::pull(celltype)
  
  # Compute confusion matrix
  cm <- confusionMatrix(data=cur_pred,
                        reference=gt_labels,
                        mode="everything")
  cm
}


# Compute celltype probabilities using pre-trained RF classifier for a SCE
compute_celltypeProb <- function(sce, rffit){
  # Extract test data to predict celltypes on
  cur_mat <- t(assay(sce, "exprs"))
  # Add subbatch information to training data and filter for cells NOT used in the
  # pretrained classifier for training or having celltype "uncertain"
  cur_mat <- cbind(cur_mat,
                   as.data.frame(colData(sce)) %>%
                     dplyr::select(cell_id, subbatch) %>%
                     mutate(dummy=1,
                            subbatch=paste0("subbatch", subbatch)) %>%
                     dcast(cell_id~subbatch, value.var="dummy", fill=0) %>%
                     dplyr::select(-cell_id)) %>%
    dplyr::filter(!(rownames(.) %in% rownames(rffit$trainingData)) &
                    !(rownames(.) %in% rownames(colData(sce)[colData(sce)$celltype=="uncertain",])))
  cur_mat <- cur_mat[colnames(cur_mat) %in% colnames(rffit$trainingData)]
  
  # Predict cell phenotypes probabilities in test data
  cur_pred <- predict(rffit, 
                      newdata=cur_mat, 
                      type="prob")
  
  # Extract ground truth labels
  gt_labels <- colData(sce) %>%
    as.data.frame() %>%
    dplyr::filter(!(rownames(.) %in% rownames(rffit$trainingData)) &
                    (celltype!="uncertain"))
  gt_labels <- gt_labels[match(rownames(cur_mat), rownames(gt_labels)), ] %>%
    dplyr::pull(celltype)
  
  # Add ground truth labels
  cur_pred$truth <- factor(gt_labels)
  
  cur_pred
}



## Plotting fncs.

# Plot single U heatmap with specified "title" from res dataframe as is read in
# read_single_U fnc
plot_single_U <- function(res, title){
  # Pivot dataframe to long format
  res <- res %>%
    pivot_wider(names_from=module, values_from=membership) %>%
    column_to_rownames(var="protein") %>%
    as.matrix()
  
  # Compute module membership for all proteins
  res.membership <- apply(res, 2, compute_module_membership)
  
  # Create heatmap of module membership
  res.heatmap <- Heatmap(res.membership, 
                         column_title=title,
                         col=structure(c("grey", pal_npg("nrc")("1")[1]), 
                                       names=c("0", "1")), 
                         show_heatmap_legend=FALSE, 
                         show_row_dend=FALSE, 
                         row_names_gp=gpar(fontsize=axis_title.fontsize),
                         column_names_gp=gpar(fontsize=axis_title.fontsize),
                         column_title_gp=gpar(fontsize=title.fontsize),
                         rect_gp=gpar(col="white", lwd=1))
  
  res.heatmap
}


# Plots U for different iter_var as heatmaps with multiple heatmaps per iter_var
# next to each other according to repetition variable
plot_U <- function(df, iter_var, repetition){
  # Empty list to save correlation matrices
  df.cor <- list()
  
  for (i in unique(df[[iter_var]])) {
    cat('#####', iter_var, ' = ', i, '\n')
    temp_list <-  list()
    
    # Create empty HeatmapList and add all repetitions as individual heatmaps
    ht_list <-NULL
    names.repetition <- c()
    for (r in unique(df[[repetition]])) {
      # Filter for current iter_var and repetition and convert to matrix used
      # by ComplexHeatmaps
      res.temp <- df %>%
        dplyr::filter(!!as.symbol(repetition)==r & !!as.symbol(iter_var)==i) %>%
        dplyr::select(-c(repetition, iter_var))
      
      if (nrow(res.temp)!=0){
        # Calculate correlation matrix
        temp_list[[length(temp_list)+1]] <- cor(t(res.temp %>%
                                                    pivot_wider(names_from=module, values_from=membership) %>%
                                                    column_to_rownames(var="protein") %>%
                                                    as.matrix()))
        
        # Create heatmap of module membership
        ht_list <- ht_list + plot_single_U(res.temp, paste0(str_to_title(repetition), 
                                                            " ", 
                                                            r, 
                                                            "\n(dictionary size: ", 
                                                            length(unique(res.temp$module)), 
                                                            ")"))
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
    
    # Add correlations to final returned df
    names(temp_list) <- names.repetition
    df.cor[[length(df.cor)+1]] <- temp_list
    
    cat('\n\n')
  }
  
  names(df.cor) <- unique(df[[iter_var]])
  
  df.cor  
}


# Plots U as heatmaps where two U's are concatenated and coloured by which 
# dataset the module belongs to
# Shared modules are only shown once and coloured accordingly
# Shown as multiple heatmaps per iter_var next to each other according to 
# repetition variable
plot_U_membership <- function(df, iter_var, repetition){
  
  for (i in unique(df[[iter_var]])) {
    cat('#####', iter_var, ' = ', i, '\n')
    temp_list <-  list()
    
    for (r in unique(df[[repetition]])) {
      # Filter for current iter_var and repetition and convert to matrix used
      # by ComplexHeatmaps
      res.temp <- df %>%
        dplyr::filter(!!as.symbol(repetition)==r & !!as.symbol(iter_var)==i) %>%
        dplyr::select(-c(repetition, iter_var)) %>%
        pivot_wider(names_from=module, values_from=membership) %>%
        column_to_rownames(var="protein") %>%
        as.matrix()
      
      # Compute module memberships
      temp_list[length(temp_list)+1] <- list(apply(res.temp, 2, compute_module_membership))
    }
    names(temp_list) <- unique(df[[repetition]])
    
    # Compare modules and retain only modules present in both U's
    modules.both <- lapply(1:ncol(temp_list[[1]]), function(i){
      unlist(lapply(1:ncol(temp_list[[2]]), function(j){
        if (all(temp_list[[1]][, i]==temp_list[[2]][, j])){
          data.frame(x=c(i), y=c(j))
        }
      }))
    }) %>% bind_rows()
    
    # Add missing modules only present in one of the two U's
    modules.matrix <- temp_list[[1]][, modules.both$x] %>%
      bind_cols(2 * temp_list[[1]][, -(modules.both$x)]) %>%
      bind_cols(3 * temp_list[[2]][, -(modules.both$y)]) %>%
      as.matrix()
    colnames(modules.matrix) <- NULL
    rownames(modules.matrix) <- rownames(temp_list[[1]])
    
    # Create heatmap of combined module membership and cluster according to jaccard distance
    modules.heatmap <- Heatmap(modules.matrix, 
                               col=structure(c("grey", pal_npg("nrc")("3")[1:3]), 
                                             names=c("0", "1", "2", "3")), 
                               show_heatmap_legend=TRUE, 
                               show_row_dend=FALSE, row_names_gp=gpar(fontsize=axis_title.fontsize),
                               column_names_gp=gpar(fontsize=axis_title.fontsize),
                               column_title_gp=gpar(fontsize=title.fontsize),
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
                                                         title_gp=gpar(fontsize=title.fontsize, 
                                                                       fontface="bold"),
                                                         labels_gp=gpar(fontsize=axis_title.fontsize)))
    print(modules.heatmap)
    cat('\n\n')
  }
}


# Plot single A heatmap with specified "title" from df dataframe as is read in
# read_single_A fnc
plot_single_A <- function(a, title){
  # Transform into matrix with channels as rownames
  a.matrix <- a %>%
    column_to_rownames("channel") %>%
    as.matrix()
  
  a.heatmap <- Heatmap(a.matrix, 
                       column_title=title,
                       col=colorRamp2(c(0, 1), colors=c("grey", pal_npg("nrc")("1"))), 
                       show_heatmap_legend=FALSE, 
                       show_row_dend=FALSE, 
                       row_names_gp=gpar(fontsize=axis_title.fontsize),
                       column_names_gp=gpar(fontsize=axis_title.fontsize),
                       column_title_gp=gpar(fontsize=title.fontsize),
                       rect_gp=gpar(col="white", lwd=1))
  
  a.heatmap
}


# Plots A for different iter_var as heatmaps with multiple heatmaps per iter_var
# next to each other according to repetition variable
plot_A <- function(df, iter_var, repetition){
  for (i in unique(df[[iter_var]])) {
    cat('#####', iter_var, ' = ', i, '\n')
    
    # Create empty HeatmapList and add all repetitions as individual heatmaps
    ht_list <-NULL
    names.repetition <- c()
    for (r in unique(df[[repetition]])) {
      # Filter for current iter_var and repetition and convert to matrix used
      # by ComplexHeatmaps
      a.temp <- df %>%
        dplyr::filter(!!as.symbol(repetition)==r & !!as.symbol(iter_var)==i) %>%
        dplyr::select(-c(repetition, iter_var))
      
      if (nrow(a.temp)!=0){
        # Create heatmap and add to list of heatmaps
        ht_list <- ht_list + plot_single_A(a.temp, paste0(str_to_title(repetition), 
                                                          " ", r))

        names.repetition <- c(names.repetition, r)
      }
    }
    
    # Draw all heatmaps together
    draw(ht_list, 
         heatmap_legend_list=list(
           Legend(col_fun=colorRamp2(c(0, 1), colors=c("grey", pal_npg("nrc")("1"))),
                  title="Chosen proteins", at=c(0, 0.25, 0.5, 0.75, 1))),
         column_title=paste0(iter_var, " = ", i))

    cat('\n\n')
  }
}


# Compute pairwise compare_cor between U matrices and plot as boxplots
# x_lab: name of x axis
# to_int: TRUE or FALSE if the variable of interest (e.g. k) in df, which is used
# as the y-axis (grouping) should be converted to int to order y-axis
# !Deprecated!
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
  
  # Convert to int if necessary to order groups in boxplot
  if (to_int){
    df.dist <- df.dist %>%
      mutate(variable=factor(df.dist$variable, levels=sort(strtoi(unique(df.dist$variable)))))
  } 
  
  # Create boxplot
  df.boxplot <- ggplot(df.dist, aes(x=variable, y=Distance, fill=variable)) +
    geom_boxplot() +
    scale_fill_npg() +
    theme_cowplot(title.fontsize) +
    theme(legend.position = "none") +
    labs(x=x_lab) 
  
  df.boxplot
}


# Compute mantel test per var_name and use metans pairs_mantel to plot
# Assumes 10 repetitions (since unfortunately a list can't be given, it had to be
# hardcoded like this)
plot_metan_mantel <- function(df, var_name){
  for (n in names(df)) {
    cat('##### ', var_name, ' = ', n, '\n')
    # Call mantel test fnc
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
# var_name: name of x axis
# to_int: TRUE or FALSE if the variable of interest (e.g. k) in df, which is used
# as the y-axis (grouping) should be converted to int to order y-axis
plot_mantel_boxplot <- function(df, var_name, to_int=TRUE){
  # Call the mantel test for all possible combinations of modules in x and y and
  # only retain only upper triangle since the mantel test is unidirectional and we 
  # don't want to plot the same values twice
  df.mantel <- data.frame((sapply((lapply(df, function(l){
    temp_comp <- sapply(l, function(x) sapply(l, function(y) mantel_test(x, y)$mantel_r))
    temp_comp[upper.tri(temp_comp)]
  })), c))) %>%
    dplyr::rename_with(~ gsub("X", "", .x, fixed = TRUE)) %>%
    melt() %>%
    dplyr::rename(mantel=value)
  
  # Convert to int if necessary to order groups in boxplot
  if (to_int){
    df.mantel <- df.mantel %>%
      mutate(variable=factor(df.mantel$variable, levels=sort(strtoi(unique(df.mantel$variable)))))
  }
  
  # Plot boxplot of mantel values
  df.mantelBoxplot <- ggplot(df.mantel, aes(x=variable, y=mantel, fill=variable)) +
    geom_boxplot() +
    scale_fill_npg() +
    theme_cowplot(title.fontsize) +
    guides(fill=FALSE) +
    labs(x=var_name, y=expression(cor["Mantel"])) 
  
  df.mantelBoxplot  
}


# Plot results from CISI for data: grouped by group on x-axis, measure is used as
# y-value and fill variable is used for coloring, results are striped according to
# noisy and no-noise results
plot_cisi_results <- function(df, group, measure, fill){
  # Rename composite to no_noise in simulation column to account for old unspecific
  # naming of these results
  df <- df %>%
    mutate(simulation=ifelse(simulation=="composite", "no_noise", simulation))
  
  # Create barplot
  barplot <- ggplot(df, aes(x=!!sym(group), y=!!sym(measure),
                            fill=!!sym(fill), pattern=simulation)) +
    geom_bar_pattern(stat="identity",
                     position=position_dodge(preserve="single"),
                     width=0.6,
                     color="black", 
                     pattern_fill="black",
                     pattern_angle=45,
                     pattern_density=0.3,
                     pattern_spacing=0.01,
                     pattern_key_scale_factor=0.6) +
    scale_y_continuous(limits=c(0, 1)) +
    scale_fill_npg() +
    scale_pattern_manual(values=c(no_noise="stripe", noisy="none"), 
                         labels=c("No noise", "Noisy (SNR=5)")) +
    labs(x=str_to_title(group), y=str_to_title(measure), pattern="Simulation type") + 
    guides(pattern=guide_legend(override.aes=list(fill="white")),
           fill=guide_legend(override.aes=list(pattern="none"))) +
    theme_cowplot(title.fontsize)
  
  barplot
}


# Plot results of expression values in "layer" of CISI vs ground truth for 
# protein of interest poi and segmented according to masks.list
plot_cells <- function(sce.list, masks.list, poi, layer="none"){
  
  # Set colours for poi
  colour.cells <- list(el=c("#FFFFFF", pal_npg("nrc")("8")[8]))
  names(colour.cells) <- poi
  
  # Since all the images and legends are returned individually, we have to
  # subset and reorder them such that the simulated and ground truth are next to
  # each other and only one legend is shown
  img.list <- list()
  idx <- c()
  for (el in sce.list){
    # Select if normal counts or transformed counts should be plotted
    if (layer=="none"){
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
    # Add plot to image list
    p.gg <- lapply(p$plot, function(x){ggdraw(x, clip="on") })
    img.list <- append(img.list, p.gg)
    idx <- c(idx, seq_along(p.gg))
  }
  
  # Reorder image list
  names(img.list) <- make.unique(names(img.list))
  idx <- order(idx)
  img.list[idx]
}


# Function that plots exprs (stored in "layer") of protein_x vs protein_y for both
# ground truth and decomposed data coloured by celltype_col
plot_exprs <- function(sce.list, celltype_col, protein_x, protein_y, layer="exprs"){
  # Remove special characters in protein names to make plotting of such proteins
  # possible
  sce.list <- lapply(sce.list, function(sce){
    rownames(sce) <- gsub("/|-", "_", rownames(sce))
    sce
    })
  protein_x <- gsub("/|-", "_", protein_x)
  protein_y <- gsub("/|-", "_", protein_y)
  
  # If "celltype" column is present in ground_truth SCE and celltype_col is given
  # then every cell specified as a celltype_col type cell will be plotted in red
  if ("celltype" %in% names(colData(sce.list[[2]]))){
    col_cells <- rep("grey", length(unique(colData(sce.list[[1]])$celltype)))
    names(col_cells) <- unique(colData(sce.list[[1]])$celltype)
    if (celltype_col!=""){
      col_cells[celltype_col] <- "red"
    }
    
    # Grey cells are made somewhat transparent to allow easier spotting of red cells
    # since there is a lot of overlap
    alpha_cells <- rep(0.2, length(unique(colData(sce.list[[1]])$celltype)))
    names(alpha_cells) <- unique(colData(sce.list[[1]])$celltype)
    alpha_cells[celltype_col] <- 1
    
    # Plot simualted and ground_truth cells
    plot.sim <- ggcells(sce.list[[1]], 
                        aes(x=!!as.symbol(protein_x), y=!!as.symbol(protein_y), 
                            colour=celltype, alpha=celltype), 
                        exprs_values=layer) +
      geom_point(size=0.3) +
      theme_cowplot(title.fontsize) +
      scale_color_manual(values=col_cells) +
      guides(color=FALSE) +
      scale_alpha_manual(guide='none', values=alpha_cells)
    plot.true <- ggcells(sce.list[[2]], 
                         aes(x=!!as.symbol(protein_x), y=!!as.symbol(protein_y), 
                             colour=celltype, alpha=celltype), 
                         exprs_values=layer) +
      geom_point(size=0.4) +
      theme_cowplot(title.fontsize) +
      scale_color_manual(values=col_cells) +
      scale_alpha_manual(guide='none', values=alpha_cells)
    
  } else {
    # If no celltype column is present, then all cells are black
    plot.sim <- ggcells(sce.list[[1]], 
                        aes(x=!!as.symbol(protein_x), y=!!as.symbol(protein_y)), 
                        exprs_values=layer) +
      geom_point(size=0.3) +
      theme_cowplot(title.fontsize)  
    
    plot.true <- ggcells(sce.list[[2]], 
                         aes(x=!!as.symbol(protein_x), y=!!as.symbol(protein_y)), 
                         exprs_values=layer) +
      geom_point(size=0.4) +
      theme_cowplot(title.fontsize) 
  }
  
  # Add plots next to each other
  plot_grid(plot.sim, plot.true, 
            labels=names(sce.list), 
            label_size=15, hjust=c(-2, -1.5), vjust=1)
}


# Plot ground_truth and simulated image "im" where cells are colored according 
# to celltypes using the colours "celltype.col"
plot_celltype_images <- function(sce, masks, im, celltype.col){
  # Create empty plot list
  plot.list <- list()
  
  # Create ground truth image with coloured celltypes
  gt.image <- plotCells(masks[im],
                        object=sce, 
                        cell_id="ObjectNumber", 
                        img_id="sample_id",
                        colour_by="celltype_NA",
                        colour=list(celltype_NA=celltype.col),
                        return_plot=TRUE,
                        image_title=list(cex=1.5),
                        legend=list(colour_by.title.cex=0.001, colour_by.labels.cex=0.001),
                        scale_bar=list(length=2, cex=1, lwidth=2),
                        display="single")
  # Create simulated image with coloured celltypes
  sim.image <- plotCells(masks[im],
                         object=sce, 
                         cell_id="ObjectNumber", 
                         img_id="sample_id",
                         colour_by="celltype_pred",
                         colour=list(celltype_pred=celltype.col),
                         return_plot=TRUE,
                         image_title=list(cex=1.5),
                         legend=list(colour_by.title.cex=0.001, colour_by.labels.cex=0.001),
                         scale_bar=list(length=2, cex=1, lwidth=2),
                         display="single")
  # Append both images together
  plot.list <- append(append(plot.list, 
                             lapply(gt.image$plot, function(x){ggdraw(x, clip="on")})), 
                      lapply(sim.image$plot, function(x){ggdraw(x, clip="on")}))
  
  # Remove legend plots (formatted weirdly) and add new legend 
  plot.list <- append(plot.list[names(plot.list)!="legend"], list(
    legend=grid.grabExpr(draw(Legend(labels=c(names(celltype.col)[1:(length(celltype.col)-1)], 
                                              "Background/NA"), 
                                     title="Celltypes",
                                     legend_gp=gpar(fill=celltype.col),
                                     grid_height=unit(5, "mm"),
                                     grid_width=unit(5, "mm"))))))
  
  # Put two images and legend into a single row of plots
  images <- plot_grid(plotlist=plot.list, nrow=1, labels=c("Ground truth", "Simulated"),
                      label_size=15, hjust=c(-2, -1.5), 
                      vjust=1, scale=0.9)
  
  images
}


# Plot barplot of correlations per protein
plot_protein_cor <- function(X.cor){
  # Create barplot with proteins on y-axis and their correlations on the y-axis
  protein.plot <- ggplot(X.cor %>%
           dplyr::select(protein, correlation) %>%
           ungroup() %>%
           distinct(),
         aes(x=reorder(protein, correlation), y=correlation, fill=protein)) +
    geom_bar(stat="identity") +
    theme_cowplot(title.fontsize) +
    scale_fill_igv() +
    guides(fill=FALSE) +
    xlab("protein") +
    coord_flip()
  
  protein.plot
}


# Plot boxplots of protein expr per dataset coloured by protein correlation results
# (performance) with "lower_limit.colour" specifying where colour gradient starts
# (for better visualizations)
plot_protein_dist <- function(X.cor, lower_limit.colour=0.7){
  # Create boxplot for each protein coloured by performance per dataset
  protein.dist <- ggplot(X.cor, 
                        aes(x=protein, y=ground_truth, color=correlation)) +
    facet_wrap(~dataset, ncol=1) +
    geom_boxplot() +
    theme_cowplot(title.fontsize) +
    scale_color_viridis(direction=-1, limits=c(lower_limit.colour, 1), oob=scales::squish) +
    stat_summary(fun=mean, geom="point", shape=8) +
    labs(y="protein expression")
  
  protein.dist
}