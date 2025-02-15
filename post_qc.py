# 2nd QC filtering after cell segmentation
# 2 strategies - filters out low intensity cells (in all channels) & removes clusters of cells
# python post_qc.py plate_no treatments_file path_to_imgs


import warnings
warnings.filterwarnings('ignore')
import pandas
import numpy
import matplotlib.pyplot as plt
import glob
import os
import skimage.io
import cellpose
from cellpose import utils, io, models, core
import sys
import time
import json
os.environ["CUDA_VISIBLE_DEVICES"] = str(sys.argv[1])

time0 = time.time()

treatments = pandas.read_csv(sys.argv[2], sep='\t')
# w = {'tubulin': 0, 'mito': 1, 'lysosome': 2, 'dapi': 3, 'brightfield': 4}

path = sys.argv[3]
p = str(sys.argv[1])


rmv_dic = dict()

for g in treatments.guide.unique():
    rmv_dic[g] = [0, []]

for m in glob.glob('masks_tif/plate' + p + '/*.tif'):
    name = m.split('/')[-1][:11]
    well = name.split('_')[1]
    g = treatments[treatments['well']==well]['guide'].item()
    mask = plt.imread(m)
    images = [numpy.load(x) for x in sorted(glob.glob(path + 'processed/plate' + p + '/cells/' + g + '/' + name[5:] + '*.npy'))]
    labels = [x.split('/')[-1].split('.')[0] for x in sorted(glob.glob(path + 'processed/plate' + p + '/cells/' + g + '/' + name[5:] + '*.npy'))]
    idxs = [int(x.split('_')[-1]) for x in labels]

    tub = []
    mit = []
    lys = []
    dap = []

    for cell in images:
        tub.append(numpy.max(cell[:,:,0]))
        mit.append(numpy.max(cell[:,:,1]))
        lys.append(numpy.max(cell[:,:,2]))
        dap.append(numpy.max(cell[:,:,3]))

    remove = []
    for idx in numpy.arange(len(cells)): 
        rmv_dic[g][0] += 1
        if (tub[idx] < 26) | (mit[idx] < 26) | (lys[idx] < 26) | (dap[idx] < 26):
            remove.append(idxs[idx])
            rmv_dic[g][1].append(idxs[idx])



    for i in numpy.unique(mask):
        if i not in idxs:
            mask[mask==i] = 0
        elif i in remove:
            mask[mask==i] = 0
            os.system('mv processed/plate' + p + '/cells/' + g + '/' + name[5:] + '_' + str(i) + '.npy' + ' unused/plate' + p + '/')
            os.system('mv processed/plate' + p + '/cells/' + g + '/' + name[5:] + '_' + str(i) + '.tif' + ' unused/plate' + p + '/')
    skimage.io.imsave('masks_tif/plate' + p + '/' + name + '.tif', mask)


with open('unused/plate' + p + '/remove_intensity.json', 'w') as file:
    json.dump(rmv_dic, file)

with open('unused/plate' + p + '/intensity_stats.txt', 'w') as f:
    f.write('\t'.join(['guide', 'total_cells', 'filtered_cells', 'pct_low_intensity']) + '\n')
    for g in rmv_dic.keys():
        f.write('\t'.join([g, str(rmv_dic[g][0]), str(len(rmv_dic[g][1])), str(round(len(rmv_dic[g][1])/rmv_dic[g][0], 4))]) + '\n')



# cyto3 model from cellpose on nuclear channel to remove cell crops with more than one nucleus:

def make_dic(dc, plate):
    for g in treatments.guide.unique():
        dc[g] = ([numpy.load(x)[:,:,3] for x in sorted(glob.glob('processed/plate' + plate + '/cells/' + g + '/*.npy'))], [x.split('/')[-1].split('.')[0] for x in sorted(glob.glob('processed/plate' + plate + '/cells/' + g + '/*.npy'))])
    print(len(dic))

dic = dict()

make_dic(dic, p)

def get_mask(dc, name, guide, diam):

    #cellpose3 built-in model:
    model = models.Cellpose(gpu=True, model_type='cyto3')
    
    img = dc[guide][0][dc[guide][1].index(name)]

    try:
        mask, flows, styles, diams = model.eval(img, diameter=diam, channels=[1,3])
    except:
        return -1
    return len(numpy.unique(mask))



masks = [plt.imread(x) for x in glob.glob('masks_tif/plate' + p + '/*.tif')]
labels = [x.split('.')[0][-6:] for x in glob.glob('masks_tif/plate' + p + '/*.tif')]


bfpmasks = dict()
remove = []

for g in treatments.guide.unique():
    # print(g)
    # t0 = time.time()
    bfpmasks[g] = []
    for i in dic[g][1]:
        name = i[:6]
        m = get_mask(dic, i, g, None)
        bfpmasks[g].append(m)
        if m != 2:
            remove.append(i)
            
            idx = int(i.split('_')[-1])
            masks[labels.index(name)][masks[labels.index(name)] == idx] = 0


with open('unused/plate' + p + '/nuclei.json', 'w') as file:
    json.dump(bfpmasks, file)

with open('unused/plate' + p + '/remove_clumps.txt', 'w') as f:
    f.write('\n'.join(x for x in remove))

with open('unused/plate' + p + '/clumping_stats.txt', 'w') as f:
    f.write('\t'.join(['guide', 'total_cells', 'filtered_cells', 'pct_clumps']) + '\n')
    for g in bfpmasks.keys():
        f.write('\t'.join([g, str(len(bfpmasks[g])), str(bfpmasks[g].count(2)), str(round(1 - bfpmasks[g].count(2)/len(bfpmasks[g]), 4))]) + '\n')

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
