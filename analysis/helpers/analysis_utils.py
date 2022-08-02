# Import libraries
import anndata as ad
import tifffile
from pathlib import Path
import os
import numpy as np
import pandas as pd


'''
Helper Functions
(Contains functions used in segmentation training)
'''

# Function that checks that file_path is an existing file and has a certain extension
def is_valid_file(file_path, extension):
    file = Path(file_path)
    if file.exists() and file.is_file() and file.suffix in extension:
        return True
    else:
        return False
    
    
# Function that checks if directory exists and contains at least one file
# with certain extension
def is_valid_directory(directory_path): 
    # Checking if the directory exists 
    if os.path.exists(directory_path): 
        # Checking if the directory is empty or not
        if any([is_valid_file(os.path.join(directory_path, f), '.tiff') 
                for f in directory_path.iterdir() if not f.name.startswith('.')]) == True:
            return True
        else:
            return False
    else:
        return  False
    
    
# Function that reads in tiff files and panel file into anndata object 
def anndata_from_tiff(tiffs_path, panel_path):
    # Read in TIFF images and flatten into one numpy array (channels x pixels)
    image_names = [f.stem for f in tiffs_path.iterdir() if not f.name.startswith('.')]
    images_unflattened = [tifffile.imread(os.path.join(tiffs_path, f)) 
                          for f in tiffs_path.iterdir() if not f.name.startswith('.')]
    images_list = [img.reshape(img.shape[0], (img.shape[1]*img.shape[2])) for img in images_unflattened]

    # Read in panel metadata
    images_panel = pd.read_csv(panel_path)

    # Create anndata from images
    # Add image intensities per pixel matrix
    images = ad.AnnData(np.transpose(np.hstack(images_list)))

    # Add observation and variable names
    num_pixels = [img.shape[1] for img in images_list]
    images.obs_names = [(ele + '_' + str(j)) for i, ele in enumerate(image_names) for j in range(num_pixels[i])]
    images.var_names = images_panel['name']

    # Add observation metadata
    images.obs = pd.DataFrame({
        'sample_id': [ele for i, ele in enumerate(image_names) for j in range(num_pixels[i])],
        'ObjectNumber': [j for i, ele in enumerate(image_names) for j in range(num_pixels[i])],
        'ROI': [ele.replace('20220520_TsH_th152_cisi1_00', '') for i, ele in enumerate(image_names) 
                for j in range(num_pixels[i])]
    })
    images.obs.index = images.obs['sample_id'].tolist()

    # Add panel data as variable metadata
    images.var = images_panel
    images.var.index = images.var['name'].tolist()
    
    return images