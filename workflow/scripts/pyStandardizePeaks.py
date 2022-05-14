import numpy as np
import pandas as pd
import glob, os

peak_file = snakemake.input['peaks']
chrom_sizes = snakemake.input['chrom_sizes']
outdir = snakemake.params['outdir']
caller = snakemake.wildcards['caller']

def parse_single_sample(file, chrom_sizes_df, out_dir = '.'):
    bio_sample, condition = os.path.basename(file).split('.')[0].split('_')[:2]
    if caller == 'macs3':
        peaks = pd.read_csv(file, sep = '\t', index_col = False, comment = '#', names = ['chr', 'start', 'abs_summit', 'name', '-log10(pvalue)'])
        standard_peak = peaks[['chr', '-log10(pvalue)', 'abs_summit']]
    else:
        peaks = pd.read_csv(file, sep = '\t', index_col = False, comment = '#', names = ['chr', 'start', 'end', 'name', 'score', 'strand', 'auc', '-log10(pvalue)', '-log10(qvalue)', 'summit_pos'])
        standard_peak = peaks[['chr', '-log10(pvalue)']]
        standard_peak['abs_summit'] = peaks['start'] + peaks['summit_pos']
    standard_peak = standard_peak.assign(start = standard_peak['abs_summit'] - 250)
    standard_peak = standard_peak.assign(end = standard_peak['abs_summit'] + 250)
    standard_peak = standard_peak.merge(chrom_sizes_df)
    standard_peak = standard_peak[(standard_peak['start'] >= 0) & (standard_peak['end'] <= standard_peak['chrom_size'])]
    i=0
    while i < standard_peak.shape[0]:
        interval = pd.arrays.IntervalArray.from_arrays(left = standard_peak.iloc[i+1:]['start'], right = standard_peak.iloc[i+1:]['end'], closed = 'both')
        standard_peak = standard_peak[~np.insert(interval.overlaps(pd.Interval(left = standard_peak.iloc[i]['start'], right = standard_peak.iloc[i]['end'], closed = 'both')), 0, [False]*(i+1))]
        i+=1
    standard_peak.loc[:,'normalized_score'] = standard_peak['-log10(pvalue)'] / (standard_peak['-log10(pvalue)'].sum() / 1e6)
    standard_peak.to_csv(os.path.join(outdir, f"{bio_sample}_{condition}.unique_501bp_peaks.txt"), index = False)

chrom_sizes_df = pd.read_csv(chrom_sizes, names = ['chr', 'chrom_size'], sep = '\t')
parse_single_sample(peak_file, chrom_sizes_df, outdir)
exit(0)