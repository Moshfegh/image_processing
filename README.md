# image-processing-measurement-analysis-ml
preprocess and segment 384-well plates of 4-channel imaged iNeurons

notes coming soon...

____________________________________________________________________________________________________________________

### image processing

--> all .py scripts run on sp3.yml requirements

--> preprocess.py - includes image QC, flatfield correction, background subtraction, and clipping & rescaling of raw tiffs
   
   requirements are path to raw image files, plate number, and a treatments file with 2 columns 'guide' & 'well'   ***NB*** well needs to be 3 characters (e.g. M05)

    example call:
    python preprocess.py /mnt/IXM/1014E/20241217-sp3-neuron-cpscreen-P1/2024-12-17/633/TimePoint_1/ treatments.txt 1

   preprocess will use qc.py

--> output folders/files:

    ALL FOLDERS HAVE SUBFOLDERS STRATIFIED BY PLATE
    
    raw_imgs - copies over all (5 channels) raw tiffs to local folder
    unused - image stacks that didn't pass QC
      **includes
    imgs_corrected - flatfield corrected images that passed QC
    rescaled_imgs - final processed images, i.e., rescaled based on min-max percentiles (0.1, 99.9 respectively) based on plate-wise (and channel-wise) intensity histograms
