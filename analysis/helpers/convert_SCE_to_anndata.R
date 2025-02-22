# Import libraries
library(zellkonverter)
library(imcRtools)

## Function that reads in SCE data and converts/saves them to pyhton anndata objects

# Write SingleCellExperiment saved in input.path to Anndata object in out.path
# for use in python
sce_to_anndata <- function(input.path, out.path){
  sce <- readRDS(input.path)
  writeH5AD(SingleCellExperiment(assays(sce),
                                 colData=colData(sce),
                                 rowData=rowData(sce)), file=out.path)
}


input.path <- "/mnt/bb_dqbm_volume/data/Immucan_lung/Lung_sce.rds"
out.path <- "/mnt/bb_dqbm_volume/data/Immucan_lung/Lung_sce_original.h5ad"
sce_to_anndata(input.path, out.path)

# sce <- readRDS(input.path)
# 
# steinbock.path <- "/mnt/bb_dqbm_volume/data/Tonsil_th152/steinbock"
# spe <- read_steinbock(steinbock.path)
