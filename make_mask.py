# python make_mask.py plate_no treatments_file path_to_calcein_imgs
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


if os.path.isdir('masks'):
	os.system('mkdir masks/plate' + p)
	os.system('mkdir masks_tif/plate' + p)
	os.system('mkdir bfp_masks/plate' + p)
else:
	os.system('mkdir masks')
	os.system('mkdir masks_tif')
	os.system('mkdir masks/plate' + p)
	os.system('mkdir masks_tif/plate' + p)
	os.system('mkdir bfp_masks')
	os.system('mkdir bfp_masks/plate' + p)

treatments = pandas.read_csv(sys.argv[2], sep='\t')

# os.system('mkdir processed/plate' + p + '/cells')
# for g in treatments.guide.unique():
	# os.system('mkdir processed/plate' + p + '/cells/' + g)

dic1 = dict()
dic2 = dict()
dic3 = dict()
dic4 = dict()

def make_dic(dic, path):

	for x in [x for x in glob.glob(path + '*.tif') if 'thumb' not in x]:
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
	# print(len(dic))
	# print(list(dic.keys())[:5])


dicdict = {'1':dic1, '2':dic2, '3':dic3, '4':dic4}
# print(dicdict)

# print(dicdict[p])

make_dic(dicdict[p], path)# + p + '/')

# w = {'tubulin': 0, 'mito': 1, 'lysosome': 2, 'dapi': 3, 'brightfield': 4}
w = {'calcein': 0, 'bfp': 1}


# generating & saving masks:
# calcein + Cellpose4 built-in SAM model
# 57 diameter
### with gpu:
bfpmodel = models.CellposeModel(gpu=True)
# model = models.CellposeModel(gpu=True, model_type='/home/ymoshfegh/1014e/20230926/CP_GCaMP')
model = models.CellposeModel(gpu=True)
w_mtx = [numpy.array([[1,0,2.13594],[0,1,-34.15589]], dtype=numpy.float32)]

def make_mask(dc, channel, well, diam):
	if dc == dic1:
		pl = '1/'
	elif dc == dic2:
		pl = '2/'
	elif dc == dic3:
		pl = '3/'
	elif dc == dic4:
		pl = '4/'
		
	bfp = cv2.imread(dc[well][w['bfp']], cv2.IMREAD_GRAYSCALE)
	dapi = numpy.load('processed/plate' + p + '/images/' + treatments[treatments['well']==well[:3]]['guide'].item() + '/' + well + '.npy')[:,:,3].astype(numpy.float32)
	img = cv2.imread(dc[well][w[channel]], cv2.IMREAD_GRAYSCALE)
	stack = cv2.merge((img, bfp))

# for calcein images w/ Cellpose4 model:
	# mask, flows, styles, diams = model.eval(stack, diameter=diam, channels=[0,0])
	mask, flows, styles = model.eval(stack, diameter=diam)#, channels=[1,3])
	bfpmask, bfpflows, bfpstyles = bfpmodel.eval(bfp, diameter=None)#, channels=[1,3])
	
	
# aligning calcein masks to processed images:
	size = dapi.shape
	warp_mode = cv2.MOTION_TRANSLATION
	warp_matrix = numpy.eye(2,3, dtype=numpy.float32)
	criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, 5000, 1e-10)
	try:
		(cc, warp_matrix) = cv2.findTransformECC (dapi, bfp.astype(numpy.float32), warp_matrix, warp_mode, criteria)
		w_mtx[0] = warp_matrix
	except:
		warp_matrix = w_mtx[0]#numpy.array([[1,0,2.13594],[0,1,-34.15589]], dtype=numpy.float32)
	mask_aligned = cv2.warpAffine(mask, warp_matrix, (size[1],size[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)
	bfpmask_aligned = cv2.warpAffine(bfpmask, warp_matrix, (size[1],size[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)


	numpy.save('masks/plate' + pl + 'mask_' + well, mask_aligned)
	numpy.save('bfp_masks/plate' + pl + 'mask_' + well, bfpmask_aligned)
# saving masks_tif (for cellprofiler) with only BFP+ cells:
# keeping .npy masks as everything - can use BFP- cells as controls in future
	bfpmask_aligned[bfpmask_aligned > 0] = 1
	intersection = mask_aligned * bfpmask_aligned

	filt = numpy.unique(intersection[intersection > 0])
	for i in numpy.unique(mask_aligned):
		if i not in filt:
			mask_aligned[mask_aligned==i] = 0

	skimage.io.imsave('masks_tif/plate' + pl + 'mask_' + well + '.tif', mask_aligned)



time0 = time.time()

for x in dicdict[p].keys():
	make_mask(dicdict[p], 'calcein', x, 57)


print('make_mask time: ' + str(round((time.time()-time0)/60, 4)) + ' minutes')


# process_mask:
import process_mask

# time1 = time.time()
process_mask

# print('process_mask time: ' + str(round((time.time()-time1)/60, 4)) + ' minutes')



# post-segmentation QC:
import post_qc
post_qc

