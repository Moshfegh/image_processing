# python preprocess.py path_to_images treatment_file plate_no
# ml environment

from matplotlib import pyplot as plt
import pandas
import glob
import numpy
from basicpy import BaSiC
import os
import scipy.stats as stats
import cv2
import skimage
import pyclesperanto_prototype as cle
import sys
import time
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

os.environ["CUDA_VISIBLE_DEVICES"] = str(sys.argv[3])

time0 = time.time()

treatments = pandas.read_csv(sys.argv[2], sep='\t')
w = {'tubulin': 0, 'mito': 1, 'lysosome': 2, 'dapi': 3, 'brightfield': 4}

path = Path(sys.argv[1].replace('\\', ''))
p = str(sys.argv[3])

if os.path.isdir('raw_imgs'):
    os.system('mkdir raw_imgs/plate' + p)
    os.system('mkdir imgs_corrected/plate' + p)
    os.system('mkdir processed/plate' + p)
    os.system('mkdir processed/plate' + p + '/images')
    os.system('mkdir rescaled_imgs/plate' + p)
else:
    os.system('mkdir raw_imgs')
    os.system('mkdir raw_imgs/plate' + p)
    os.system('mkdir imgs_corrected')
    os.system('mkdir imgs_corrected/plate' + p)
    os.system('mkdir processed')
    os.system('mkdir processed/plate' + p)
    os.system('mkdir processed/plate' + p + '/images')
    os.system('mkdir rescaled_imgs')
    os.system('mkdir rescaled_imgs/plate' + p)

for g in treatments.guide.unique():
    os.system('mkdir processed/plate' + p + '/images/' + g)


r = [_file.absolute().as_posix() for _file in path.iterdir()]
raw_imgs = [f for f in r if ('tif' in f) & ('thumb' not in f)]
print('# images: ' + str(len(raw_imgs)))
for f in raw_imgs:
    name = f.split('/')[-1]
    tmp = name.split('_')[-3]
    # os.system('ln -s ' + f + ' raw_imgs/')
    if tmp in list(treatments.well):
        os.system('cp ' + f + ' raw_imgs/plate' + p + '/' + name)


# remove bad images:
import qc
qc


tubu = [plt.imread(x) for x in sorted(glob.glob('raw_imgs/plate' + p + '/*.tif')) if '_w1' in x]
tubu_labels = [x.split('/')[-1] for x in sorted(glob.glob('raw_imgs/plate' + p + '/*.tif')) if '_w1' in x]
print(len(tubu))


dapi = [plt.imread(x) for x in sorted(glob.glob('raw_imgs/plate' + p + '/*.tif')) if '_w4' in x]
dapi_labels = [x.split('/')[-1] for x in sorted(glob.glob('raw_imgs/plate' + p + '/*.tif')) if '_w4' in x]
print(len(dapi))


mito = [plt.imread(x) for x in sorted(glob.glob('raw_imgs/plate' + p + '/*.tif')) if '_w2' in x]
mito_labels = [x.split('/')[-1] for x in sorted(glob.glob('raw_imgs/plate' + p + '/*.tif')) if '_w2' in x]
print(len(mito))


lyso = [plt.imread(x) for x in sorted(glob.glob('raw_imgs/plate' + p + '/*.tif')) if '_w3' in x]
lyso_labels = [x.split('/')[-1] for x in sorted(glob.glob('raw_imgs/plate' + p + '/*.tif')) if '_w3' in x]
print(len(lyso))


tubu = numpy.array(tubu)
mito = numpy.array(mito)
lyso = numpy.array(lyso)
dapi = numpy.array(dapi)


def correct(channel, names):
    basic = BaSiC(get_darkfield=True, smoothness_flatfield=1)
    basic.fit(channel)
    images_transformed = basic.transform(channel)
    return (images_transformed, basic)


tub_flatfield, b = correct(tubu, tubu_labels)
# plot the fit results:
fig, axes = plt.subplots(1, 3, figsize=(9, 3))
im = axes[0].imshow(b.flatfield)
fig.colorbar(im, ax=axes[0])
axes[0].set_title("Flatfield")
im = axes[1].imshow(b.darkfield)
fig.colorbar(im, ax=axes[1])
axes[1].set_title("Darkfield")
axes[2].plot(b.baseline)
axes[2].set_xlabel("Frame")
axes[2].set_ylabel("Baseline")
fig.tight_layout()
fig.savefig('imgs_corrected/plate' + p + '/tubulin_flatfield_p' + p + '.png')

mit_flatfield, m = correct(mito, mito_labels)
# plot the fit results:
fig, axes = plt.subplots(1, 3, figsize=(9, 3))
im = axes[0].imshow(m.flatfield)
fig.colorbar(im, ax=axes[0])
axes[0].set_title("Flatfield")
im = axes[1].imshow(m.darkfield)
fig.colorbar(im, ax=axes[1])
axes[1].set_title("Darkfield")
axes[2].plot(m.baseline)
axes[2].set_xlabel("Frame")
axes[2].set_ylabel("Baseline")
fig.tight_layout()
fig.savefig('imgs_corrected/plate' + p + '/mito_flatfield_p' + p + '.png')

lys_flatfield, l = correct(lyso, lyso_labels)
# plot the fit results:
fig, axes = plt.subplots(1, 3, figsize=(9, 3))
im = axes[0].imshow(l.flatfield)
fig.colorbar(im, ax=axes[0])
axes[0].set_title("Flatfield")
im = axes[1].imshow(l.darkfield)
fig.colorbar(im, ax=axes[1])
axes[1].set_title("Darkfield")
axes[2].plot(l.baseline)
axes[2].set_xlabel("Frame")
axes[2].set_ylabel("Baseline")
fig.tight_layout()
fig.savefig('imgs_corrected/plate' + p + '/lysosome_flatfield_p' + p + '.png')

dap_flatfield, d = correct(dapi, dapi_labels)
# plot the fit results:
fig, axes = plt.subplots(1, 3, figsize=(9, 3))
im = axes[0].imshow(d.flatfield)
fig.colorbar(im, ax=axes[0])
axes[0].set_title("Flatfield")
im = axes[1].imshow(d.darkfield)
fig.colorbar(im, ax=axes[1])
axes[1].set_title("Darkfield")
axes[2].plot(d.baseline)
axes[2].set_xlabel("Frame")
axes[2].set_ylabel("Baseline")
fig.tight_layout()
fig.savefig('imgs_corrected/plate' + p + '/nuclei_flatfield_p' + p + '.png')


# background subtraction

def bg_subtract(img, rad=100):
    result = cle.top_hat_sphere(img, radius_x=rad, radius_y=rad)
    result = numpy.asarray(result)
    return result


tub_corrected = []
mit_corrected = []
lys_corrected = []
dap_corrected = []

for i in tub_flatfield:
    tub_corrected.append(bg_subtract(numpy.array(i)))

for i in mit_flatfield:
    mit_corrected.append(bg_subtract(numpy.array(i)))

for i in lys_flatfield:
    lys_corrected.append(bg_subtract(numpy.array(i)))

for i in dap_flatfield:
    dap_corrected.append(bg_subtract(numpy.array(i)))


# for channel, label in [(tub_corrected, tubu_labels), (mit_corrected, mito_labels), (lys_corrected, lyso_labels), (dap_corrected, dapi_labels)]:
    # for i in range(len(channel)):
        # name = label[i]
        # img = (channel[i]/255).astype(numpy.uint8)
        # # img = channel[i]
        # skimage.io.imsave('imgs_corrected/plate' + p + '/' + name, img)


tub_corrected = numpy.array(tub_corrected)
mit_corrected = numpy.array(mit_corrected)
lys_corrected = numpy.array(lys_corrected)
dap_corrected = numpy.array(dap_corrected)


# rescaling based on 16-bit corrected images:

pmin0, pmax0 = stats.scoreatpercentile([tub_corrected[x] for x in range(len(tub_corrected))], (0.1, 99.9))
print(pmin0, pmax0)
pmin1, pmax1 = stats.scoreatpercentile([mit_corrected[x] for x in range(len(mit_corrected))], (0.1, 99.9))
print(pmin1, pmax1)
pmin2, pmax2 = stats.scoreatpercentile([lys_corrected[x] for x in range(len(lys_corrected))], (0.1, 99.9))
print(pmin2, pmax2)
pmin3, pmax3 = stats.scoreatpercentile([dap_corrected[x] for x in range(len(dap_corrected))], (0.1, 99.9))
print(pmin3, pmax3)


images = []
labels = []

for idx in range(len(tub_corrected)):
    t0 = tub_corrected[idx]
    l0 = tubu_labels[idx]
    t1 = mit_corrected[idx]
    l1 = mito_labels[idx]
    t2 = lys_corrected[idx]
    l2 = lyso_labels[idx]
    t3 = dap_corrected[idx]
    l3 = dapi_labels[idx]

    imgs = cv2.merge((t0, t1, t2, t3))
    l = [l0, l1, l2, l3]
    images.append(imgs)
    labels.append(l)


for img in range(len(images)):
    tmp0 = (skimage.exposure.rescale_intensity(images[img][:,:,0], in_range=(pmin0, pmax0))*255).astype(numpy.uint8)
    skimage.io.imsave('rescaled_imgs/plate' + p + '/' + labels[img][0], tmp0)

    tmp1 = (skimage.exposure.rescale_intensity(images[img][:,:,1], in_range=(pmin1, pmax1))*255).astype(numpy.uint8)
    skimage.io.imsave('rescaled_imgs/plate' + p + '/' + labels[img][1], tmp1)

    tmp2 = (skimage.exposure.rescale_intensity(images[img][:,:,2], in_range=(pmin2, pmax2))*255).astype(numpy.uint8)
    skimage.io.imsave('rescaled_imgs/plate' + p + '/' + labels[img][2], tmp2)

    tmp3 = (skimage.exposure.rescale_intensity(images[img][:,:,3], in_range=(pmin3, pmax3))*255).astype(numpy.uint8)
    skimage.io.imsave('rescaled_imgs/plate' + p + '/' + labels[img][3], tmp3)

    stack = cv2.merge((tmp0, tmp1, tmp2, tmp3))

    name = labels[img][0].split('_')[1] + '_' + labels[img][0].split('_')[2]

    numpy.save('processed/plate' + p + '/images/' + treatments[treatments.well == name[:3]].guide.item() + '/' + name, stack)

print('preprocess done')
print('time: ' + str(round((time.time()-time0)/60, 4)) + ' minutes')
