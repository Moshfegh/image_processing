import sys
import numpy
import matplotlib.pyplot as plt
import skimage.io as io
import os
import glob
import warnings
warnings.filterwarnings('ignore')

p = str(sys.argv[3])

if os.path.isdir('unused'):
    os.system('mkdir unused/plate' + p)
else:
    os.system('mkdir unused')
    os.system('mkdir unused/plate' + p)


# function to display images:
def view_img(dc, well):
    t = io.imread(dc[well][w['tubulin']])
    m = io.imread(dc[well][w['mito']])
    l = io.imread(dc[well][w['lysosome']])

    # fig = plt.figure(figsize=(36, 20))
    fig = plt.figure()

    fig, ax = plt.subplots(nrows=1, ncols=3)

    ax[0].imshow(t, cmap='gray')
    ax[0].set_title('tubulin')
    ax[0].axis('off')

    ax[1].imshow(m, cmap='gray')
    ax[1].set_title('mitochondria')
    ax[1].axis('off')

    ax[2].imshow(l, cmap='gray')
    ax[2].set_title('lysosome')
    ax[2].axis('off')

    plt.suptitle(well, y=.75, fontsize=12)
    fig.tight_layout()
    # plt.show()
    plt.savefig('unused/plate' + p + '/' + well + '.png')


w = {'tubulin': 0, 'mito': 1, 'lysosome': 2, 'dapi': 3, 'brightfield': 4}


def make_dic(dic, path):
    for x in glob.glob(path + '*.tif'):# + glob.glob('unused_imgs/*.tif'):
        if x.split('_')[2] + '_' + x.split('_')[3] in dic.keys():
            dic[x.split('_')[2] + '_' + x.split('_')[3]].append(x)
        else: dic[x.split('_')[2] + '_' + x.split('_')[3]] = [x]
    print(len(dic))

    for x, y in dic.items():
        dic[x] = sorted(y)


dic = dict()
make_dic(dic, 'raw_imgs/plate' + p + '/')


def get_img(dc, well, channel):
    img = io.imread(dc[well][w[channel]])
    return img


tubu = []
mito = []
lyso = []
for d in dic.keys():
    tubu.append(numpy.std(get_img(dic, d, 'tubulin')))
    mito.append(numpy.std(get_img(dic, d, 'mito')))
    lyso.append(numpy.std(get_img(dic, d, 'lysosome')))

keys = []
for d in dic.keys():
    keys.append(d)
len(keys)

outliers = []
for x in tubu:
    if (x >= numpy.mean(tubu) + 2*numpy.std(tubu)) | (x <= numpy.mean(tubu) - 2*numpy.std(tubu)):
        outliers.append(keys[tubu.index(x)])
for x in mito:
    if (x >= numpy.mean(mito) + 2*numpy.std(mito)) | (x <= numpy.mean(mito) - 2*numpy.std(mito)):
        outliers.append(keys[mito.index(x)])
for x in lyso:
    if (x >= numpy.mean(lyso) + 2*numpy.std(lyso)) | (x <= numpy.mean(lyso) - 2*numpy.std(lyso)):
        outliers.append(keys[lyso.index(x)])

outliers = list(set(outliers))
len(outliers)

with open('unused/plate' + p + '/remove.txt', 'w') as f:
    f.write('\n'.join(x for x in outliers))

for x in outliers:
    view_img(dic, x)


for x in glob.glob('raw_imgs/plate' + p + '/*.tif'):
    for r in outliers:
        if r in x:
            os.system('mv ' + x + ' unused/plate' + p + '/')
