#!/usr/bin/env bash

# change directory
BASEDIR="Tonsil_th152"
cd "${BASEDIR}"

# setup steinbock alias
shopt -s expand_aliases
alias steinbock="docker run -e STEINBOCK_MASK_DTYPE=uint32  -v ${BASEDIR}:/data \
-u $(id -u):$(id -g) ghcr.io/bodenmillergroup/steinbock:0.14.1"

# preprocessing
steinbock preprocess imc panel --imcpanel config/TH152_panel.csv --mcd mcd/1 \
-o steinbock/panel.csv
steinbock preprocess imc images --hpf 50 --mcd mcd/1 --panel steinbock/panel.csv \
--imgout steinbock/img --infoout steinbock/images.csv

# classification using existing classifier
# steinbock classify ilastik prepare --cropsize 500 --seed 123
# rm pixel_classifier.ilp && mv IMCWorkflow.ilp pixel_classifier.ilp
# rm -r ilastik_crops && mv analysis/crops ilastik_crops
# steinbock classify ilastik fix --no-backup
# steinbock classify ilastik run
#
# # random forest-based segmentation
# steinbock segment cellprofiler prepare
# steinbock segment cellprofiler run -o masks_ilastik

# deep learning-based segmentation
steinbock segment deepcell --minmax --panel steinbock/panel.csv --img steinbock/img \
-o steinbock/masks_deepcell

# measurement
steinbock measure intensities --panel steinbock/panel.csv --img steinbock/img \
--masks steinbock/masks_deepcell -o steinbock/intensities
steinbock measure regionprops --img steinbock/img --masks steinbock/masks_deepcell \
-o steinbock/regionprops
steinbock measure neighbors --masks steinbock/masks_deepcell --type expansion \
--dmax 4 -o steinbock/neighbors

# export
# steinbock export ome
# steinbock export histocat --masks masks_deepcell
# steinbock export csv intensities regionprops -o cells.csv
# steinbock export csv intensities regionprops --no-concat -o cells_csv
# steinbock export fcs intensities regionprops -o cells.fcs
# steinbock export fcs intensities regionprops --no-concat -o cells_fcs
# steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors -o cells.h5ad
# steinbock export anndata --intensities intensities --data regionprops --neighbors neighbors --no-concat -o cells_h5ad
# steinbock export graphs --data intensities

# zip -r cells_csv.zip cells_csv
# zip -r cells_fcs.zip cells_fcs
# zip -r cells_h5ad.zip cells_h5ad
# zip -r graphs.zip graphs
# zip -r histocat.zip histocat
# zip -r ilastik_crops.zip ilastik_crops
# zip -r ilastik_img.zip ilastik_img
# zip -r ilastik_probabilities.zip ilastik_probabilities
# zip -r img.zip img
# zip -r intensities.zip intensities
# zip -r masks_deepcell.zip masks_deepcell
# zip -r masks_ilastik.zip masks_ilastik
# zip -r neighbors.zip neighbors
# zip -r ome.zip ome
# zip -r regionprops.zip regionprops
