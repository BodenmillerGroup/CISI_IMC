#!/usr/bin/env bash

# print start
echo "----------------------------------------------------"
echo $(date +"%D %T") "Init steinbock..."
echo "----------------------------------------------------"

# change directory
# BASEDIR="/mnt/bb_dqbm_volume/data/Tonsil_th152/steinbock"
BASEDIR="/mnt/bb_dqbm_volume/data/20221108_TsH_LSS_cisiabmix1_179/steinbock"
cd "${BASEDIR}"

# setup steinbock alias
shopt -s expand_aliases
alias steinbock="docker run -e STEINBOCK_MASK_DTYPE=uint32 -v ${BASEDIR}:/data \
-u $(id -u):$(id -g) ghcr.io/bodenmillergroup/steinbock:0.14.2"


# print starting steinbock
echo "----------------------------------------------------"
echo $(date +"%D %T") "Start preprocessing...."
echo "----------------------------------------------------"

# preprocessing
steinbock preprocess imc panel
steinbock preprocess imc images --hpf 50


# print steinbock segmentation starts
echo "----------------------------------------------------"
echo $(date +"%D %T") "Start segmentation..."
echo "----------------------------------------------------"

# deep learning-based segmentation
steinbock segment deepcell --minmax -o masks_deepcell


# print steinbock meausurements starts
echo "----------------------------------------------------"
echo $(date +"%D %T") "Start computing measurements..."
echo "----------------------------------------------------"

# measurement
steinbock measure intensities --masks masks_deepcell 
steinbock measure regionprops --masks masks_deepcell
steinbock measure neighbors --masks masks_deepcell --type expansion --dmax 4


# print done
echo "----------------------------------------------------"
echo $(date +"%D %T") "Steinbock pipeline is done."
echo $(date +"%D %T") "Files can be found at: $BASEDIR"
