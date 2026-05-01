#   torch-gpu env

# pca plot - aggreggated by well and grouped by gene, plus distance heatmap

# python pca.py <cellprofiler table> <metadata_col>


import sys
import pandas
import matplotlib
import matplotlib.pyplot as plt
import numpy
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler
from scipy import stats
from scipy.spatial.distance import cdist
import seaborn
import warnings
warnings.filterwarnings('ignore')


cp = pandas.read_csv(sys.argv[1], sep='\t')
column = sys.argv[2]



def pca_plot(table, feature, mark='o', c='tab20b'):

    cols = table.select_dtypes(include='number').columns
    
    table = table.groupby(['well',feature]).mean().reset_index()

    norm = StandardScaler().fit_transform(table[cols].values)
    
    pca = PCA(n_components=2)
    pc = pca.fit_transform(norm)
    
    # print('Explained variation per principal component: {}'.format(pca.explained_variance_ratio_))

    pca_df = pandas.DataFrame(data = pc, columns = ['pc1', 'pc2'])# + ['na']*98)
    pca_df['well'] = list(table['well'])
    # pca_df['guide'] = table['guide']
    pca_df[feature] = list(table[feature])

    plt.figure(figsize=(8,5))
    # seaborn.jointplot(
    seaborn.scatterplot(
    x="pc1", y="pc2",
    hue=feature, hue_order=sorted(list(table[feature].unique())),
    palette=seaborn.color_palette(c, table[feature].nunique()),
    data=pca_df,
    legend="full",
    # marker="$\circ$", ec="face", # <-- for unfilled circles
    marker = mark,#'+',
    alpha=1)
    plt.title('PCA plot of CP features\n' + 'Explained variation per principal component: ' + str(pca.explained_variance_ratio_))#, fontsize=20)
    plt.legend(bbox_to_anchor=(1, 1), ncol=max([1, int(numpy.floor(table[feature].nunique()/25))]))
    return pca_df
    plt.tight_layout()
    plt.show()



def cluster_dist_heatmap(table, feature='guide', cm='Blues_r', f=(15,10)):

    tmp = pandas.DataFrame(columns=list(table[feature].unique()), index=list(table[feature].unique()))


    for x in tmp.columns:
        for y in tmp.columns:
            a = table[table[feature]==x][[table.columns[0], table.columns[1]]]
            b = table[table[feature]==y][[table.columns[0], table.columns[1]]]


            if x == y:
                between = 0.0
        
            else:
            # ctrl = numpy.ravel(cdist(a,a, 'euclidean')).mean()
                between = numpy.ravel(cdist(a,b, 'euclidean')).mean()

            tmp.loc[x, y] = between
    tmp = tmp.apply(pandas.to_numeric, errors='coerce')
    ######
    tmp.sort_index(axis=0, ascending=True, inplace=True)
    tmp = tmp[list(tmp.index)]
############

    plt.figure(figsize=f)
    
    seaborn.heatmap(data=tmp, cmap=cm, annot_kws={'fontsize':30}, square=True)
    # hue_order=sorted(list(table[feature].unique()))
    # plt.legend()
    plt.xticks(rotation=90)
    plt.title('euclidean distances between treatments', fontsize=f[0])
    plt.tight_layout()
    plt.show()


df = pca_plot(cp, column)

cluster_dist_heatmap(df, feature=column)
