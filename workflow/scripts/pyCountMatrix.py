import os
import numpy as np
import pandas as pd
from multiprocessing import Pool

ncores = snakemake.resources['cpus']
bed_files = snakemake.input['bed_se']
union_peak_file = snakemake.input['peaks']
out_file_raw = snakemake.output['matrix_raw']
out_file_norm = snakemake.output['matrix_norm']


def count_transpositions_by_chr(reads, peaks, chrom, name):
    sub_reads = reads.loc[reads['chr'] == chrom, :]
    sub_peaks = union_peaks.loc[union_peaks['chr'] == chrom].copy()
    pos = np.concatenate([sub_reads[sub_reads['name'].str.endswith('1')]['start'].values, sub_reads[sub_reads['name'].str.endswith('2')]['end'].values])
    sub_peaks = sub_peaks.sort_values('start', ascending = True)
    bins = pd.Series(pd.cut(pos, sub_peaks.index.get_level_values(1))).value_counts()
    bins.index = pd.MultiIndex.from_arrays([[chrom]*bins.shape[0], bins.index])
    return bins
    
def count_transpositions(bed_file, peaks, ncores):
    chrs = union_peaks['chr'].unique()
    name = os.path.basename(bed_file).replace('.bed', '')
    print(name)
    reads = pd.read_csv(bed_file, sep = '\t', names = ['chr', 'start', 'end', 'name', 'score', 'strand'], dtype = {'chr': 'category', 'start': np.int32, 'end': np.int32, 'score': np.int32, 'strand': 'category'})
    with Pool(ncores) as p:
        results = p.starmap(count_transpositions_by_chr, [(reads, peaks, chrom, name) for chrom in chrs])
    counts = pd.Series(dtype = np.int32)
    for i in results:
        counts = pd.concat([counts, pd.Series(i)])
    peaks.loc[:, name] = counts

union_peaks = pd.read_csv(union_peak_file)
union_peaks = union_peaks.sort_values(['chr', 'start', 'end'])
peak_intervals = pd.arrays.IntervalArray.from_arrays(left = union_peaks['start'], right = union_peaks['end'], closed = 'both')
union_peaks.index = pd.MultiIndex.from_arrays([union_peaks['chr'], pd.IntervalIndex(peak_intervals)])

for bed_file in bed_files:
    count_transpositions(bed_file, union_peaks, ncores)

union_peaks.to_csv(out_file_raw)
matrix = union_peaks.iloc[:, 5:]
matrix_norm = np.log((matrix / (matrix.sum(axis = 0) + 1)) + 1)
matrix_norm.to_csv(out_file_norm)