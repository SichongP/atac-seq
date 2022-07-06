rule bamCoverage:
    input:
        bam = f"{OUTDIR}/shifted_bam/{{sample, [\w\d_]+}}.bam",
        size_factor = f"{OUTDIR}/deseq2/size_factors.txt"
    output: f"{OUTDIR}/bigwig/{{sample}}.bw"
    conda: "../envs/deeptools.yaml"
    params: partition = 'med2'
    threads: 5
    resources: cpus=5, time_min=1200, mem_mb=20000, cpus_bmm=5, mem_mb_bmm=20000
    shell:
     """
     scale=$(grep {wildcards.sample} {input.size_factor} | cut -d ',' -f 2)
     bamCoverage -b {input.bam} -o {output} -of bigwig -bs 50 -p {resources.cpus} --effectiveGenomeSize {config[effective_size]} --scaleFactor $scale  --ignoreDuplicates --minMappingQuality 20
     """
     
rule multiBigWig:
    input: expand(f"{OUTDIR}/bigwig/{{sample}}.bw", sample = pep.sample_table['sample_name'])
    output:
        npz = f"{OUTDIR}/deeptools/all.npz"
    conda: "../envs/deeptools.yaml"
    params: partition = 'med2'
    threads: 15
    resources: cpus=15, time_min=6600, mem_mb=30000, cpus_bmm=15, mem_mb_bmm=30000, partition = 'med2'
    shell:
     """
     multiBigwigSummary bins -b {input} -o {output.npz} --chromosomesToSkip chrX chrM MSY --smartLabels -bs 50 -p {resources.cpus} 
     """
     
rule plotPCA:
    input: f"{OUTDIR}/deeptools/all.npz"
    output:
        png = f"{OUTDIR}/figures/PCA.png",
        tab = f"{OUTDIR}/figures/PCA.tab"
    conda: "../envs/deeptools.yaml"
    params: partition = 'med2'
    threads: 1
    resources: cpus=1, time_min=600, mem_mb=4000, cpus_bmm=1, mem_mb_bmm=4000, partition = 'med2'
    shell:
     """
     plotPCA  --transpose -in {input} -o {output.png} -T "PCA, Normalized Genome Coverage" --outFileNameData {output.tab} --ntop 0 
     """
     
rule plotCorrelation:
    input: f"{OUTDIR}/deeptools/all.npz"
    output:
        png = f"{OUTDIR}/figures/corr.png",
        tab = f"{OUTDIR}/figures/corr.tab"
    conda: "../envs/deeptools.yaml"
    params: partition = 'med2'
    threads: 1
    resources: cpus=1, time_min=600, mem_mb=4000, cpus_bmm=1, mem_mb_bmm=4000, partition = 'med2'
    shell:
     """
     plotCorrelation -c pearson -p heatmap -o {output.png} --skipZeros -T "Pearson Correlation of Genome Coverage" --removeOutliers
     """
     
rule computeMatrix:
    input:
        bw = expand(f"{OUTDIR}/bigwig/{{sample}}.bw", sample = pep.sample_table['sample_name']),
        gtf = {config['gff']}
    output: mtx = f"{OUTDIR}/deeptools/mtx.gz"
    conda: "../envs/deeptools.yaml"
    params: partition = 'med2'
    threads: 15
    resources: cpus=15, time_min=6600, mem_mb=30000, cpus_bmm=15, mem_mb_bmm=30000, partition = 'med2'
    shell:
     """
     computeMatrix scale-regions -S {input.bw} -R {input.gtf} -o {output.mtx} -b 3000 -a 3000 -m 5000 --skipZeros --smartLabels -p {resources.cpus} --metagene
     """
     
rule plotHeatmap:
    input: f"{OUTDIR}/deeptools/mtx.gz"
    output:
        png = f"{OUTDIR}/figures/heatmap.png",
        tab = f"{OUTDIR}/deeptools/heatmap.tab",
        bed = f"{OUTDIR}/deeptools/mtx.sorting.bed"
    conda: "../envs/deeptools.yaml"
    params: partition = 'med2'
    threads: 1
    resources: cpus=1, time_min=600, mem_mb=3000, cpus_med2=1, mem_mb_med2=3000, partition = 'med2'
    shell:
     """
     plotHeatmap --dpi 300 -m {input} -o {output.png} --outFileNameMatrix {output.tab} --outFileSortedRegions {output.bed} --colorMap 'Set2' --endLabel "TTS" -T "Heatmap, ATAC-seq coverage around genes"
     """     