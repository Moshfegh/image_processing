#   torch-gpu env

# "final" volcano plot of hits

# python hits.py <cellprofiler table> <metadata_col> <control>


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
control = sys.argv[3]


def pcadf(table, feature):

    cols = table.select_dtypes(include='number').columns
    
    table = table.groupby(['well',feature]).mean().reset_index()

    norm = StandardScaler().fit_transform(table[cols].values)
    
    pca = PCA(n_components=2)
    pc = pca.fit_transform(norm)

    pca_df = pandas.DataFrame(data = pc, columns = ['pc1', 'pc2'])
    pca_df['well'] = table['well']    
    pca_df[feature] = table[feature]

    return pca_df


def cluster_dist(g, table, ntc, feature='guide'):

    distances = []
    a = table[table[feature]==ntc][[table.columns[0], table.columns[1]]]
    ctrl = numpy.ravel(cdist(a,a, 'euclidean'))

    b = table[table[feature]==g][[table.columns[0], table.columns[1]]]
    guide = numpy.ravel(cdist(a, b, 'euclidean'))


    if b.shape[0] > 0:
        u = stats.mannwhitneyu(ctrl, guide, alternative='two-sided')[1]
        t = stats.ttest_ind(ctrl, guide, alternative='two-sided', equal_var=False, trim=0)[1]
        ks = stats.ks_2samp(ctrl, guide, alternative='two-sided')[1]

        pvalues = [u, t, ks]
        fdr = stats.false_discovery_control(pvalues)
                
        p = fdr.max()
        # if p < 1e-100:
            # p = 1e-20

    else: p=1

    return numpy.abs(guide.mean()-ctrl.mean()), p


def make_dist(table, ntc, feature='guide'):
    guides = []
    distances = []
    pvalues = []

    for guide in table[feature].unique():
        a, b = cluster_dist(guide, table[table.guide.isin([ntc, guide])], ntc)
        guides.append(guide)
        distances.append(a)
        pvalues.append(b)

    dist = pandas.DataFrame({'well':guides, 'distance':distances, 'pvalue':[-numpy.log10(x) for x in pvalues]})
    dist = dist[dist.pvalue > 0]
    return dist


def dist_plot(table, ntc, measure, p, d):
    t = measure

    df = table.copy()
    
    plt.figure(figsize=(7, 7), dpi=100)
    plt.scatter(x=df.distance,y=df.pvalue,s=200,label="not significant", color='grey', alpha=.5)
    
    up = df[(df['distance']>d)&(df['pvalue']>p)]
    
    plt.scatter(x=up['distance'],y=up['pvalue'],s=200,label="significant distance",color="lightblue")
    for i,r in up.iterrows():
        plt.text(x=r['distance'],y=r['pvalue'],s=r.well)
    
    plt.xlabel(f"distance to {ntc}")
    plt.ylabel("-log10(padj)")
    # plt.axvline(.8,color="grey",linestyle="--")
    plt.axvline(d,color="grey",linestyle="--")
    plt.axhline(p,color="grey",linestyle="--")
    plt.title(t, fontweight='bold')
    plt.show()



df = pcadf(cp, column)

dist_plot(make_dist(df, control, column), control, f'pca distances to {control}', 2, 2)