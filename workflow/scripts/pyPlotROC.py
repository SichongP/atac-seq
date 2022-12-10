import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
import pybedtools as pbt
        
atac = snakemake.input['atac_peaks']
h3k4me1 = snakemake.input['h3k4me1_peaks']
h3k4me3 = snakemake.input['h3k4me3_peaks']
h3k27me3 = snakemake.input['h3k27me3_peaks']
h3k27ac = snakemake.input['h3k27ac_peaks']
chrom_sizes = snakemake.input['chrom_sizes']

out_png = snakemake.output['png']
out_csv = snakemake.output['csv']

atac_peaks = pd.read_csv(atac, sep = '\t', names = ['chr', 'start', 'end', 'strand', 'score'])

h3k4me1_peaks = pbt.BedTool(h3k4me1).slop(b = 200, g = chrom_sizes)
h3k4me3_peaks = pbt.BedTool(h3k4me3).slop(b = 200, g = chrom_sizes)
h3k27ac_peaks = pbt.BedTool(h3k27ac).slop(b = 200, g = chrom_sizes)
h3k27me3_peaks = pbt.BedTool(h3k27me3).slop(b = 200, g = chrom_sizes)

active = (h3k4me1_peaks.merge(i=h3k4me3_peaks) + h3k27ac_peaks - h3k27me3_peaks).sort().merge()
inactive = h3k27me3_peaks.sort().merge()
RP_count = active.count()
RN_count = inactive.count()
ratio_count_df = pd.DataFrame(columns=['cutoff', 'ratio', 'value'])
for cutoff in list(range(0,int(atac_peaks['score'].quantile(.99)),1))+list(range(int(atac_peaks['score'].quantile(.9)), int(max(atac_peaks['score'])),int(max(atac_peaks['score'])/50))):
    filtered = atac_peaks[atac_peaks['score'] > cutoff]
    active_overlap = active.intersect(pbt.bedtool.BedTool.from_dataframe(filtered), u = True).count()
    inactive_overlap = inactive.intersect(pbt.bedtool.BedTool.from_dataframe(filtered), u = True).count()
    ratio_count_df = pd.concat([ratio_count_df, pd.DataFrame(data = {'cutoff': [cutoff], 'ratio': 'active', 'value': [active_overlap]})])
    ratio_count_df = pd.concat([ratio_count_df, pd.DataFrame(data = {'cutoff': [cutoff], 'ratio': 'inactive', 'value': [inactive_overlap]})])
ratio_count_df['cutoff'] = ratio_count_df['cutoff'].astype(int)
ratio_count_df['value'] = ratio_count_df['value'].astype(int)
ratio_count_df_pivot = ratio_count_df.pivot_table(index = 'cutoff', columns = 'ratio', values = 'value')
ratio_count_df_pivot.columns = ['TP', 'FP']
ratio_count_df_pivot['precision'] = ratio_count_df_pivot['TP']/(ratio_count_df_pivot['TP'] + ratio_count_df_pivot['FP'])
ratio_count_df_pivot['TPR'] = ratio_count_df_pivot['TP']/RP_count
ratio_count_df_pivot['FPR'] = ratio_count_df_pivot['FP']/RN_count
ratio_count_df_pivot.to_csv(out_csv)


fig, ax = plt.subplots(figsize = (6, 4))
sns.lineplot(data = ratio_count_df_pivot, x = 'FPR', y = 'TPR', ax = ax)
lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),
    np.max([ax.get_xlim(), ax.get_ylim()]),
]

ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
ax.set_aspect('equal')
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_xlabel('FPR', fontsize = 20)
ax.set_ylabel('TPR', fontsize = 20)
ax.set_title('ROC Curve', fontsize = 30)
fig.savefig(out_png, dpi = 300)