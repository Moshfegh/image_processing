# python process_mask.py plate_no treatments_file

import cellpose
from cellpose import utils, io, models, core
from cellpose import plot
import numpy
import matplotlib.pyplot as plt
import skimage.io
import shutil
import os
import cv2
import glob
from scipy import ndimage
import time
import pandas
import sys


treatments = pandas.read_csv(sys.argv[2], sep='\t')
p = str(sys.argv[1])

# os.system('mkdir processed/plate' + p)
os.system('mkdir processed/plate' + p + '/cells')
for g in treatments.guide.unique():
    os.system('mkdir processed/plate' + p + '/cells/' + g)


path = 'rescaled_imgs/plate'
w = {'tubulin': 0, 'mito': 1, 'lysosome': 2, 'dapi': 3, 'brightfield': 4}

dic1 = dict()
dic2 = dict()
dic3 = dict()
dic4 = dict()

def make_dic(dic, path):

    for x in glob.glob(path + '*.tif'):
        if x.split('_')[-3] + '_' + x.split('_')[-2] in dic.keys():
            dic[x.split('_')[-3] + '_' + x.split('_')[-2]].append(x)
        else: dic[x.split('_')[-3] + '_' + x.split('_')[-2]] = [x]
    print(len(dic))
    
    for x, y in dic.items():
        dic[x] = sorted(y)


# make_dic(dic1, path+'1/')
# make_dic(dic2, path+'2/')
# make_dic(dic3, path+'3/')
# make_dic(dic4, path+'4/')


dicdict = {'1':dic1, '2':dic2, '3':dic3, '4':dic4}

make_dic(dicdict[p], path + p + '/')

# processing images w/ masks & saving single cell stacks:
def process_mask(dic, well):
    if dic == dic1:
        pl = '1/'
    elif dic == dic2:
        pl = '2/'
    elif dic == dic3:
        pl = '3/'
    elif dic == dic4:
        pl = '4/'
    
    tubu = cv2.imread(dic[well][0], cv2.IMREAD_GRAYSCALE)
    mito = cv2.imread(dic[well][1], cv2.IMREAD_GRAYSCALE)
    lyso = cv2.imread(dic[well][2], cv2.IMREAD_GRAYSCALE)
    dapi = cv2.imread(dic[well][w['dapi']], cv2.IMREAD_GRAYSCALE)

    # stack = cv2.merge((tubu, mito, lyso, dapi))
    # numpy.save('processed/plate' + p + 'images/' + well, stack)
    
    
    mask = numpy.load('masks/plate' + pl + 'mask_' + well + '.npy')
    r = 60
    d = mask.shape[0]
    
    mask[mask==numpy.min(mask)] = 0
    
    if numpy.max(mask) > 0:
        
        for x in list(set(mask[mask != 0])):
            m = numpy.where(mask != x, 0, mask)
            if numpy.max(m)>0:
                m[m > 1] = 1
                im = mito * m
                it = tubu * m
                il = lyso * m
                ida = dapi * m
                centroid = ndimage.center_of_mass(m)
                if (r < centroid[0] < d-r) & (r < centroid[1] < d-r):
                    
                    tcell = it[round(centroid[0])-r:round(centroid[0])+r, round(centroid[1])-r:round(centroid[1])+r]#.astype(numpy.uint8)
            
                    mcell = im[round(centroid[0])-r:round(centroid[0])+r, round(centroid[1])-r:round(centroid[1])+r]#.astype(numpy.uint8)
                
                    lcell = il[round(centroid[0])-r:round(centroid[0])+r, round(centroid[1])-r:round(centroid[1])+r]#.astype(numpy.uint8)

                    dcell = ida[round(centroid[0])-r:round(centroid[0])+r, round(centroid[1])-r:round(centroid[1])+r]#.astype(numpy.uint8)

                    if numpy.max(mcell) > 0:
                        stcell = cv2.merge((tcell, mcell, lcell, dcell))
                        numpy.save('processed/plate' + pl + 'cells/' + treatments[treatments.well == well[:3]].guide.item() + '/' + well + '_' + str(x), stcell)


t0 = time.time()
for x in dicdict[p].keys():
    process_mask(dicdict[p], x)

print('process_mask time: ' + str((time.time()-t0)/60) + ' minutes')
