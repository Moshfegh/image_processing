# python make_mask.py plate_no treatments_file path_to_imgs
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


# path = 'rescaled_imgs/plate'
path = sys.argv[3]
p = str(sys.argv[1])

if os.path.isdir('calcein_masks'):
    os.system('mkdir calcein_masks/plate' + p)
    # os.system('mkdir masks_tif/plate' + p)
    os.system('mkdir calcein_bfp/plate' + p)
else:
    os.system('mkdir calcein_masks')
    # os.system('mkdir masks_tif')
    os.system('mkdir calcein_masks/plate' + p)
    # os.system('mkdir masks_tif/plate' + p)
    os.system('mkdir calcein_bfp')
    os.system('mkdir calcein_bfp/plate' + p)

treatments = pandas.read_csv(sys.argv[2], sep='\t')

# os.system('mkdir processed/plate' + p + '/cells')
# for g in treatments.guide.unique():
    # os.system('mkdir processed/plate' + p + '/cells/' + g)

dic1 = dict()
dic2 = dict()
dic3 = dict()
dic4 = dict()

def make_dic(dic, path):

	for x in glob.glob(path + '*.tif'):
		if x.split('_')[-3] + '_' + x.split('_')[-2] in dic.keys():
			dic[x.split('_')[-3] + '_' + x.split('_')[-2]].append(x)
		else: dic[x.split('_')[-3] + '_' + x.split('_')[-2]] = [x]
	# print(len(dic))
	
	for x, y in dic.items():
		dic[x] = sorted(y)

	remove = pandas.read_csv('unused/plate' +  p + '/remove.txt', header=None)
	for i in list(remove[0]):
		if i in dic.keys():
			del dic[i]
	rmv = []
	for x in dic.keys():
		if x[:3] not in list(treatments.well):
			rmv.append(x)
	for i in rmv:
		del dic[i]


dicdict = {'1':dic1, '2':dic2, '3':dic3, '4':dic4}
# print(dicdict)


make_dic(dicdict[p], path)# + p + '/')

# w = {'tubulin': 0, 'mito': 1, 'lysosome': 2, 'dapi': 3, 'brightfield': 4}
w = {'calcein': 0, 'bfp': 1}


# generating & saving masks:
# mito + cyto3 model from cellpose3 for ineurons - 20X
# NONE diameter
### with gpu:
# model = models.Cellpose(gpu=True, model_type='cyto3')
model = models.CellposeModel(gpu=True, model_type='/home/ymoshfegh/1014e/20230926/CP_GCaMP')


def make_mask(dc, channel, well, diam):
	if dc == dic1:
		pl = '1/'
	elif dc == dic2:
		pl = '2/'
	elif dc == dic3:
		pl = '3/'
	elif dc == dic4:
		pl = '4/'
		
	dapi = cv2.imread(dc[well][w['bfp']], cv2.IMREAD_GRAYSCALE)
	# mito = cv2.imread(dc[well][1], cv2.IMREAD_GRAYSCALE)
	img = cv2.imread(dc[well][w[channel]], cv2.IMREAD_GRAYSCALE)
	stack = cv2.merge((img, dapi))

# for calcein images w/ fine-tuned model:
	# mask, flows, styles = model.eval(stack, diameter=diam, channels=[0,0])
	mask, flows, styles = model.eval(img, diameter=diam, channels=[1,3])
	bfpmask, bfpflows, bfpstyles = model.eval(dapi, diameter=None, channels=[1,3])

	# mask, flows, styles, diams = model.eval(stack, diameter=diam, channels=[0,0])
	# bfpmask, bfpflows, bfpstyles, bfpdiams = model.eval(dapi, diameter=None, channels=[1,3])
	# numpy.save('masks/plate' + pl + 'mask_' + well, mask)
	# numpy.save('bfp_masks/plate' + pl + 'mask_' + well, bfpmask)
	numpy.save('calcein_bfp/plate' + pl + 'mask_' + well, bfpmask)
# saving masks_tif (for cellprofiler) with only BFP+ cells:
# keeping .npy masks as everything - can use BFP- cells as controls in future
	bfpmask[bfpmask > 0] = 1
	intersection = mask * bfpmask

	filt = numpy.unique(intersection[intersection > 0])
	for i in numpy.unique(mask):
		if i not in filt:
			mask[mask==i] = 0

	skimage.io.imsave('calcein_masks/plate' + pl + 'mask_' + well + '.tif', mask)
	#print(well)



time0 = time.time()

for x in dicdict[p].keys():
	make_mask(dicdict[p], 'calcein', x, 35)

print('make_mask time: ' + str((time.time()-time0)/60) + ' minutes')


# process_mask:
import calcein_process_mask

# time1 = time.time()
calcein_process_mask

# print('process_mask time: ' + str((time.time()-time1)/60) + ' minutes')
