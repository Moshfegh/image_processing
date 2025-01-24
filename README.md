# image-processing-measurement-analysis-ml
preprocess and segment 384-well plates of 4-channel imaged iNeurons

notes coming soon...

____________________________________________________________________________________________________________________

### image processing

--> all python scripts run on sp3.yml requirements

--> preprocess.py - includes image QC, flatfield correction, background subtraction, and clipping & rescaling of raw tiffs
   
   * requirements are path to raw image files, plate number, and a treatments .txt file with 2 columns 'guide' & 'well'   ***NB*** well needs to be 3 characters (e.g. M05)

    example call:
    python preprocess.py /mnt/IXM/1014E/20241217-sp3-neuron-cpscreen-P1/2024-12-17/633/TimePoint_1/ treatments.txt 1

   * preprocess will use qc.py

   * output folders/files:

    ALL FOLDERS HAVE SUBFOLDERS STRATIFIED BY PLATE
    
    raw_imgs - copies over all (5 channels) raw tiffs to local folder - image names should include the following for each channel:
          tubulin: _w1
          mitochondria: _w2
          lysosomes: _w3
          nuclear: _w4
          brightfield: _w5
    unused - image stacks (all 5 channels) that didn't pass QC  **includes .png images of removed images for viewing and text file with list of FOVs (i.e. well_site)
    imgs_corrected - flatfield corrected and background subtracted tiffs
    rescaled_imgs - final processed images, i.e., rescaled based on min-max percentiles (0.1, 99.9 respectively) based on plate-wise (and channel-wise) intensity histograms, saved as tiffs; these are input for cellprofiler
    processed - (images subfolder) FOV stacks saved as .npy; these are image input for machine learning, CNN, transformers, image analysis, etc.

### cell segmentation

--> make_mask.py - uses Cellpose 3.0 with cyto3 model + mitochondria channel to create and save segmented masks

      example call:
      python make_mask.py 4 treatments.txt

   * output folders/files:

     masks - cellular masks for each FOV saved as .npy
     masks_tif - masks filtered for only BFP+ cells, saved as .tif; these are also an input for cellprofiler
     bfp_masks - nuclear masks, used for calculating cell counts per well and checking lenti incorporation

--> process_mask.py - uses output masks from make_mask to isolate single cell image stacks and saves .npy in processed/cells/ folder

   * used by make_mask.py but can also be run separately:

      example call:
      python process_mask.py 4 treatments.txt


