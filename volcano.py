# torch-gpu env

# volcano plot of individual perturbations

# python volcano.py <cellprofiler table> <control> <gene> <metadata_col>

import sys
import pandas
import matplotlib.pyplot as plt
import numpy
from scipy.stats import ttest_ind, mannwhitneyu, ks_2samp
from scipy import stats
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler
import warnings
warnings.filterwarnings('ignore')

table = pandas.read_csv(sys.argv[1], sep='\t')

metadata = table[table.select_dtypes(exclude='number').columns]
table = table[table.select_dtypes(include='number').columns]

scaler = MinMaxScaler()

scaled = scaler.fit_transform(table)

table = pandas.DataFrame(scaled, columns = table.columns)
table = pandas.concat([table, metadata], axis=1)

guide = sys.argv[3]
ntc = sys.argv[2]
column = sys.argv[4]

def volcano(g, blank, cp, meta, alpha=1e-15):

    features = []
    logFC = []
    FC = []
    stat = []
    FDR = []
    sig_features = []
    # g = 'DNM1L_1'
    d = cp.copy()
    for col in d.select_dtypes(include='number').columns:
        # print(col)
        tmp1 = d[d[meta]==g][col]
        tmp2 = d[d[meta]==blank][col]
        # s1 = numpy.sum(numpy.square(tmp1))/len(tmp1)
    
        s1 = numpy.mean(tmp1)
        s2 = numpy.mean(tmp2)
    
        # if len(tmp1)>10:
        stat1, pt = ttest_ind(tmp1, tmp2, alternative='two-sided', equal_var=False, trim=0)
        stat2, pu = mannwhitneyu(tmp1, tmp2, alternative='two-sided')
        stat3, pk = ks_2samp(tmp1, tmp2, alternative='two-sided')

        pvalues = [pt, pu, pk]
        pvalues = [x for x in pvalues if str(x) != 'nan']
        fdr = stats.false_discovery_control(pvalues)
        
        if len(fdr) > 0:
            p = fdr.max()
            # if p < 1e-100:
            #     p = 1e-100
            # 
            # print('stat=%.3f, p=%.3f' % (stat, p))
        
            features.append(col)
            logFC.append(numpy.log2(s1/s2))
            FC.append(s1/s2)
            stat.append(max(pvalues))
            FDR.append(p)
            
            if p < alpha:
                # print('Probably the same distribution')
            # else:
                # print('Probably different distributions ************************************************')
                sig_features.append(col)

    print(f"#significant features: {len(sig_features)}")


    logp = [-numpy.log10(x) for x in FDR]

    df = pandas.DataFrame({'feature':features, 'pvalue':stat, 'padj':FDR, 'log10padj':logp, 'FC':FC, 'logFC':logFC})
    

    plt.scatter(x=df.logFC,y=df.log10padj,s=20,label="not significant", color='grey', alpha=.5)
    
    # # highlight down- or up- regulated genes
    # down = df[(df['FC']<=.75)&(df['stat']<=0.01)]
    # up = df[(df['FC']>=1.25)&(df['stat']<=0.01)]

    down = df[(df['logFC']<-.5)&(df['padj']<alpha)]
    up = df[(df['logFC']>.5)&(df['padj']<alpha)]
    
    plt.scatter(x=down['logFC'],y=down['log10padj'],s=20,label=f"lower in {g}",color="skyblue")
    plt.scatter(x=up['logFC'],y=up['log10padj'],s=20,label=f"higher in {g}",color="royalblue")
    # for i,r in up.iterrows():
    #     plt.text(x=r['FC'],y=r['pvalue'],s=r.feature)
    # for i,r in down.iterrows():
    #     plt.text(x=r['FC'],y=r['pvalue'],s=r.feature)
    
    plt.xlabel(f"log2FC - {g}/{blank}")
    plt.ylabel("-log10(padj)")
    plt.axvline(-.5,color="grey",linestyle="--")
    plt.axvline(.5,color="grey",linestyle="--")
    plt.axhline(-numpy.log10(alpha),color="grey",linestyle="--")
    plt.legend(bbox_to_anchor=(1, 1))
    plt.title(g, fontweight='bold')
    # plt.tight_layout()
    plt.show()


volcano(guide, ntc, table, column)
