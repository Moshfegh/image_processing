# image_processing
preprocess and segment 384-well plates of 4-channel imaged iNeurons

notes coming soon...





--> output folders/files:

    ALL FOLDERS HAVE SUBFOLDERS STRATIFIED BY PLATE
    raw_imgs - copies over all (5 channels) raw tiffs to local folder
    unused - image stacks that didn't pass QC
      **includes
    imgs_corrected - flatfield corrected images that passed QC
    rescaled_imgs - final processed images, i.e., rescaled based on min-max percentiles (0.1, 99.9 respectively) based on plate-wise (and channel-wise) intensity histograms
