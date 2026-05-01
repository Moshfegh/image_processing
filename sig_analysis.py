# torch-gpu env
# python sig_analysis.py <cellprofiler table>

import pandas
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn
import sys
import json
import time
import warnings
warnings.filterwarnings('ignore')


table = pandas.read_csv(sys.argv[1], sep='\t', index_col=0)



def pval_test(ctl, treatment):
    d = dict(zip(treatment.columns, [[] for x in range(len(treatment.columns))]))
    measurements = []
    for x in treatment.columns[:-4]:
        dtmp = ctl[x]
        tmp = treatment[x]
        if (dtmp.shape[0] > 10) & (tmp.shape[0] > 10):
            u = stats.mannwhitneyu(dtmp, tmp, alternative='two-sided')
            if u[1] < 10e-15:
                d[x].append(u[1])
            t = stats.ttest_ind(dtmp, tmp, alternative='two-sided', equal_var=False, trim=0)
            if t[1] < 10e-15:
                d[x].append(t[1])
            ks = stats.ks_2samp(dtmp, tmp, alternative='two-sided')
            if ks[1] < 10e-15:
                d[x].append(ks[1])

    for k, v in d.items():
        if len(v) == 3:
            measurements.append(k)

    return measurements



dic = dict()
for nt in table.guide.unique():
    dic[nt] = dict()
    for g in table.guide.unique():
        dic[nt][g] = pval_test(table[table.guide==nt], table[table.guide==g])



def make_data_dic(d):
    data = dict()
    for x in d.keys():
        data[x] = []
        for k, v in d[x].items():
            data[x].append(len(v))
    return data


def make_heatmap(d, cm=None):
    df = pandas.DataFrame([x for x in d.values()], columns=list(d.keys()), index=list(d.keys()))

    sort_dic = dict()
    for x, y in enumerate(sorted(df.index)):
        sort_dic[y] = x
    
    sorter = [sort_dic[x] for x in df.index]
    df['sort'] = sorter
    
    df = df.sort_values(["sort"])[sorted(df.index)]

    f = (12,9)
    
    plt.figure(figsize=f)
    seaborn.heatmap(data=df, cmap=cm)
    plt.show()


make_heatmap(make_data_dic(dic))
