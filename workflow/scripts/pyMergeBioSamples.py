import numpy as np
import pandas as pd
import glob, os

peak_files = snakemake.input['peak_files']
outfile = snakemake.output[0]

def merge_bio_samples(files, outfile = '.'):
    peaks = pd.DataFrame()
    for file in files:
        peaks = pd.concat([peaks, pd.read_csv(file, index_col = False)])
    peaks = peaks.sort_values('normalized_score', ascending = False)
    while i < peaks.shape[0]:
        interval = pd.arrays.IntervalArray.from_arrays(left = peaks.iloc[i+1:]['start'], right = peaks.iloc[i+1:]['end'], closed = 'both')
        peaks = peaks[~np.insert(interval.overlaps(pd.Interval(left = peaks.iloc[i]['start'], right = peaks.iloc[i]['end'], closed = 'both')), 0, [False]*(i+1))]
        i+=1
    norm_score = peaks['normalized_score'] / (peaks['normalized_score'].sum() / 1e6)
    if 'chrom_size' in peaks.columns:
        peaks = peaks.drop(['chrom_size', '-log10(pvalue)'], axis = 1)
    peaks.loc[:,'normalized_score'] = norm_score
    peaks.to_csv(outfile, index = False)

merge_bio_samples(peak_files, outfile)
exit(0)