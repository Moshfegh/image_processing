# bar plot of total recovered cells/treatment

# python recovered_cells.py <cellprofiler table>

import sys
import pandas
import matplotlib.pyplot as plt

table = pandas.read_csv(sys.argv[1], sep='\t')

counts = pandas.DataFrame(table.guide.value_counts())
counts.plot(kind='bar', figsize=(40,10), grid=False, color='k', fontsize=20)
plt.title('total recovered cells/treatment', fontsize=50, fontweight='bold')
plt.tight_layout()
plt.show()
