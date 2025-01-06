# python make_mask.py plate_no treatments_file
# ml environment


import os
import sys
os.environ["CUDA_VISIBLE_DEVICES"] = str(sys.argv[1])

import cellpose
from cellpose import utils, io, models, core
from cellpose import plot
import numpy
import matplotlib.pyplot as plt
import skimage.io
import shutil
import cv2
import glob
from scipy import ndimage
import time
import warnings
warnings.filterwarnings('ignore')
import pandas


path = 'raw_imgs/plate'
p = str(sys.argv[1])

# masks folders created in preprocess:
os.system('mkdir masks/plate' + p)
os.system('mkdir masks_tif/plate' + p)

treatments = pandas.read_csv(sys.argv[2], sep='\t')

# os.system('mkdir processed/plate' + p)
os.system('mkdir processed/plate' + p + '/cells')
for g in treatments.guide.unique():
    os.system('mkdir processed/plate' + p + '/cells/' + g)

dic1 = dict()
dic2 = dict()
dic3 = dict()
dic4 = dict()

def make_dic(dic, path):

	for x in glob.glob(path + '*.tif'):
		if x.split('_')[2] + '_' + x.split('_')[3] in dic.keys():
			dic[x.split('_')[2] + '_' + x.split('_')[3]].append(x)
		else: dic[x.split('_')[2] + '_' + x.split('_')[3]] = [x]
	print(len(dic))
	
	for x, y in dic.items():
		dic[x] = sorted(y)


dicdict = {'1':dic1, '2':dic2, '3':dic3, '4':dic4}
print(dicdict)

# print(dicdict[p])

make_dic(dicdict[p], path + p + '/')

w = {'tubulin': 0, 'mito': 1, 'lysosome': 2, 'dapi': 3, 'brightfield': 4}
print(w)


# generating & saving masks:
# mito + cyto3 model from cellpose3 for ineurons - 20X
# NONE diameter
### with gpu:
model = models.Cellpose(gpu=True, model_type='cyto3')
# model = models.CellposeModel(gpu=True,  model_type='/home/ymoshfegh/1014e/20230926/20240125/CP_GCaMP')


def make_mask(dc, channel, well, diam):
	if dc == dic1:
		pl = '1/'
	elif dc == dic2:
		pl = '2/'
	elif dc == dic3:
		pl = '3/'
	elif dc == dic4:
		pl = '4/'
		
	dapi = cv2.imread(dc[well][3], cv2.IMREAD_GRAYSCALE)
	# mito = cv2.imread(dc[well][1], cv2.IMREAD_GRAYSCALE)
	img = cv2.imread(dc[well][w[channel]], cv2.IMREAD_GRAYSCALE)
	stack = cv2.merge((img, dapi))
	

	mask, flows, styles, diams = model.eval(stack, diameter=diam, channels=[0,0])
	# mask, flows, styles, diams = model.eval(img, diameter=diam, channels=[1,3])
	# mask, flows, styles = model.eval(mito, diameter=125, channels=[1,3])
	# skimage.io.imsave(path + '/mask.png', mask)
	# matplotlib.image.imsave('masks/mask_' + well + '.png', mask)
	numpy.save('masks/plate' + pl + 'mask_' + well, mask)
	skimage.io.imsave('masks_tif/plate' + pl + 'mask_' + well + '.tif', mask)
	#print(well)


for x in dicdict[p].keys():
	make_mask(dicdict[p], 'mito', x, None)





# process_mask:
import process_mask


# def process_mask(dic, well):
    # if dic == dic1:
        # pl = '1/'
    # elif dic == dic2:
        # pl = '2/'
    # elif dic == dic3:
        # pl = '3/'
    # elif dic == dic4:
        # pl = '4/'
    
    # tubu = cv2.imread(dic[well][0], cv2.IMREAD_GRAYSCALE)
    # mito = cv2.imread(dic[well][1], cv2.IMREAD_GRAYSCALE)
    # lyso = cv2.imread(dic[well][2], cv2.IMREAD_GRAYSCALE)
    # dapi = cv2.imread(dic[well][w['dapi']], cv2.IMREAD_GRAYSCALE)

    # # stack = cv2.merge((tubu, mito, lyso, dapi))
    # # numpy.save('processed/plate' + pl + 'images/' + well, stack)
    
    
    # mask = numpy.load('masks/plate' + pl + 'mask_' + well + '.npy')
    # r = 75
    # d = mask.shape[0]
    
    # mask[mask==numpy.min(mask)] = 0
    
    # if numpy.max(mask) > 0:
        
        # for x in list(set(mask[mask != 0])):
            # m = numpy.where(mask != x, 0, mask)
            # if numpy.max(m)>0:
                # m[m > 1] = 1
                # im = mito * m
                # it = tubu * m
                # il = lyso * m
                # ida = dapi * m
                # centroid = ndimage.center_of_mass(m)
                # if (r < centroid[0] < d-r) & (r < centroid[1] < d-r):
                    
                    # tcell = it[round(centroid[0])-r:round(centroid[0])+r, round(centroid[1])-r:round(centroid[1])+r]#.astype(numpy.uint8)
            
                    # mcell = im[round(centroid[0])-r:round(centroid[0])+r, round(centroid[1])-r:round(centroid[1])+r]#.astype(numpy.uint8)
                
                    # lcell = il[round(centroid[0])-r:round(centroid[0])+r, round(centroid[1])-r:round(centroid[1])+r]#.astype(numpy.uint8)

                    # dcell = ida[round(centroid[0])-r:round(centroid[0])+r, round(centroid[1])-r:round(centroid[1])+r]#.astype(numpy.uint8)

                    # if numpy.max(mcell) > 0:
                        # stcell = cv2.merge((tcell, mcell, lcell, dcell))
                        # numpy.save('processed/plate' + pl + 'cells/' + treatments[treatments['well'] == well[:3]].guide.item() + '/' + well + '_' + str(x), stcell)


t0 = time.time()
#for x in dicdict[p].keys():
#    process_mask(dicdict[p], x)

print('process_mask time:')
print(time.time()-t0)
