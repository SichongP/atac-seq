import os
import numpy as np
import pandas as pd
from multiprocessing import Pool

ncores = snakemake.resources['cpus']
bed_files = snakemake.input['bed_se']
union_peak_file = snakemake.input['peaks']
out_file = snakemake.output['matrix']

union_peaks = pd.read_csv(union_peak_file)
#union_peaks = union_peaks.sort_values(['chr', 'start', 'end'])

def find_peak(read, peaks):
    peaks_interval = pd.arrays.IntervalArray.from_arrays(left = peaks['start'], right = peaks['end'], closed = 'both')
    cutting_site = read['start'] if read['read'] == 1 else read['end']
    return np.all([peaks_interval.contains(cutting_site), peaks['chr'] == read['chr']], axis = 0)

def count_file(bed_file, peaks):
    name = os.path.basename(bed_file).replace('.bed', '')
    full_reads = pd.read_csv(bed_file, sep = '\t', names = ['chr', 'start', 'end', 'name', 'mapq', 'strand'], engine = 'c')
    full_reads['read'] = np.where(full_reads['name'].str.endswith('1'), 1, 2)
    mat = full_reads.apply(find_peak, peaks = union_peaks, axis = 1)
    return (name, np.array([i for i in mat]).sum(axis = 0))
    
with Pool(ncores) as p:
    results = p.starmap(count_file, [(bed_file, union_peaks) for bed_file in bed_files])
    for name, result in results:
        union_peaks.loc[:, name] = result
        
union_peaks['peak'] = union_peaks['chr'] + '_' + union_peaks['start'].astype(str) + '_' + union_peaks['end'].astype(str)
matrix = union_peaks.iloc[:, 5:].set_index('peak')
matrix = np.log((matrix / (matrix.sum(axis = 1) + 1)) + 1)
matrix.to_csv(out_file)