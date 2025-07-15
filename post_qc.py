# 2nd QC filtering after cell segmentation
# 3 strategies - filters out low intensity cells (in all channels), filter for unaligned cropped cells, & removes clusters of cells
# python post_qc.py plate_no treatments_file


import pandas
import numpy
import matplotlib.pyplot as plt
import glob
import os
import skimage.io
from cellpose import models
import sys
import time
import json
from skimage.segmentation import clear_border
import warnings


warnings.filterwarnings('ignore')
os.environ["CUDA_VISIBLE_DEVICES"] = str(sys.argv[1])

time0 = time.time()

treatments = pandas.read_csv(sys.argv[2], sep='\t')
# w = {'tubulin': 0, 'mito': 1, 'lysosome': 2, 'dapi': 3, 'brightfield': 4}

# path = sys.argv[3]
p = str(sys.argv[1])


rmv_dic = dict()

for g in treatments.guide.unique():
    rmv_dic[g] = [0, []]

for m in glob.glob('masks_tif/plate' + p + '/*.tif'):
    name = m.split('/')[-1][:11]
    well = name.split('_')[1]
    g = treatments[treatments['well'] == well]['guide'].item()
    mask = plt.imread(m)
    images = [numpy.load(x) for x in sorted(glob.glob('processed/plate' + p + '/cells/' + g + '/' + name[5:] + '*.npy'))]
    labels = [x.split('/')[-1].split('.')[0] for x in sorted(glob.glob('processed/plate' + p + '/cells/' + g + '/' + name[5:] + '*.npy'))]
    idxs = [int(x.split('_')[-1]) for x in labels]

    tub = []
    mit = []
    lys = []
    dap = []
    cellmasks = []
    clear_borders = []

    for cell in images:
        tub.append(numpy.max(cell[:,:,0]))
        mit.append(numpy.max(cell[:,:,1]))
        lys.append(numpy.max(cell[:,:,2]))
        dap.append(numpy.max(cell[:,:,3]))
        cellmask = cell[:,:,1].copy()
        cellmask[cellmask > 1] = 1
        cm = clear_border(cellmask)
        cellmasks.append(cellmask)
        clear_borders.append(cm)

    remove = []
    for idx in numpy.arange(len(images)):
        rmv_dic[g][0] += 1
# filter out cells with low intensity:
        # if (((tub[idx] < 26) & (mit[idx] < 26) & (lys[idx] < 26)) | (dap[idx] < 26)):
        if ((numpy.max(images[idx][:,:,:3]) < 26) | (dap[idx] < 26)):
            remove.append(idxs[idx])
            rmv_dic[g][1].append(idxs[idx])
# filter out cells that weren't aligned well and got cropped in process_mask.py:
        elif (cellmasks[idx] != clear_borders[idx]).any():
            remove.append(idxs[idx])
            rmv_dic[g][1].append(idxs[idx])

    for i in numpy.unique(mask):
        if i not in idxs:
            mask[mask == i] = 0
        elif i in remove:
            mask[mask == i] = 0
            os.system('mv processed/plate' + p + '/cells/' + g + '/' + name[5:] + '_' + str(i) + '.npy' + ' unused/plate' + p + '/')
            os.system('mv processed/plate' + p + '/cells/' + g + '/' + name[5:] + '_' + str(i) + '.tif' + ' unused/plate' + p + '/')
    skimage.io.imsave('masks_tif/plate' + p + '/' + name + '.tif', mask)


with open('unused/plate' + p + '/remove_intensity.json', 'w') as file:
    json.dump(rmv_dic, file)

with open('unused/plate' + p + '/intensity_stats.txt', 'w') as f:
    f.write('\t'.join(['guide', 'total_cells', 'filtered_cells', 'pct_removed']) + '\n')
    for g in rmv_dic.keys():
        if rmv_dic[g][0] == 0:
            pct = '0'
        else: pct = str(round(len(rmv_dic[g][1])/rmv_dic[g][0], 4))
        f.write('\t'.join([g, str(rmv_dic[g][0]), str(len(rmv_dic[g][1])), pct]) + '\n')


# Cellpose-SAM model from cellpose4 on nuclear channel to remove cell crops with more than one nucleus:

def make_dic(dc, plate):
    for g in treatments.guide.unique():
        dc[g] = ([numpy.load(x)[:,:,3] for x in sorted(glob.glob('processed/plate' + plate + '/cells/' + g + '/*.npy'))], [x.split('/')[-1].split('.')[0] for x in sorted(glob.glob('processed/plate' + plate + '/cells/' + g + '/*.npy'))])
    print(len(dic))


dic = dict()


make_dic(dic, p)


#cellpose4 built-in model:
model = models.CellposeModel(gpu=True)


def get_mask(dc, name, guide, diam):

    img = dc[guide][0][dc[guide][1].index(name)]

    mask, flows, styles = model.eval(img, diameter=diam, cellprob_threshold=-2, flow_threshold=.95)

    return len(numpy.unique(mask))


masks = [plt.imread(x) for x in glob.glob('masks_tif/plate' + p + '/*.tif')]
labels = [x.split('.')[0][-6:] for x in glob.glob('masks_tif/plate' + p + '/*.tif')]


bfpmasks = dict()
remove_0 = []
remove_clump = []

# using diameter=45 for nuclei

for g in treatments.guide.unique():
    bfpmasks[g] = []
    for i in dic[g][1]:
        name = i[:6]
        m = get_mask(dic, i, g, 45)
        bfpmasks[g].append(m)
        if m == 1:
            remove_0.append(i)
            idx = int(i.split('_')[-1])
            masks[labels.index(name)][masks[labels.index(name)] == idx] = 0

        elif m > 2:
            remove_clump.append(i)

            idx = int(i.split('_')[-1])
            masks[labels.index(name)][masks[labels.index(name)] == idx] = 0


with open('unused/plate' + p + '/nuclei.json', 'w') as file:
    json.dump(bfpmasks, file)

with open('unused/plate' + p + '/remove_clumps.txt', 'w') as f:
    f.write('\n'.join(x for x in remove_clump))

with open('unused/plate' + p + '/remove_no_nuclei.txt', 'w') as f:
    f.write('\n'.join(x for x in remove_0))

with open('unused/plate' + p + '/clumping_stats.txt', 'w') as f:
    f.write('\t'.join(['guide', 'total_cells', 'filtered_cells', 'pct_clumps']) + '\n')
    for g in bfpmasks.keys():
        if len(bfpmasks[g]) > 0:
       # f.write('\t'.join([g, str(len(bfpmasks[g])), str(len(bfpmasks[g]) - bfpmasks[g].count(2)), str(round(1 - bfpmasks[g].count(2)/len(bfpmasks[g]), 4))]) + '\n')
            f.write('\t'.join([g, str(len(bfpmasks[g]) - bfpmasks[g].count(1)), str(len(bfpmasks[g]) - bfpmasks[g].count(1) - bfpmasks[g].count(2)), str(round(1 - bfpmasks[g].count(2)/(len(bfpmasks[g]) - bfpmasks[g].count(1)), 4))]) + '\n')


for x in bfpmasks.keys():
    plt.hist(bfpmasks[x], bins=10)
plt.savefig('unused/plate' + p + '/clumping_histogram_p' + p + '.png')


for g in treatments.guide.unique():
    for i in range(len(dic[g][1])):
        name = dic[g][1][i][:6]
        if bfpmasks[g][i] != 2:
            idx = int(dic[g][1][i].split('_')[-1])
            os.system('mv processed/plate' + p + '/cells/' + g + '/' + name + '_' + str(idx) + '.npy' + ' unused/plate' + p + '/')
            os.system('mv processed/plate' + p + '/cells/' + g + '/' + name + '_' + str(idx) + '.tif' + ' unused/plate' + p + '/')

            idx = int(dic[g][1][i].split('_')[-1])
            masks[labels.index(name)][masks[labels.index(name)] == idx] = 0


for cm in range(len(masks)):
    skimage.io.imsave('masks_tif/plate' + p + '/mask_' + labels[cm] + '.tif', masks[cm])


print('post-segmentation QC time: ' + str((time.time()-time0)/60) + ' minutes')
