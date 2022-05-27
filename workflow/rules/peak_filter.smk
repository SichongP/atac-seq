rule plot_roc:
    input:
        atac_peaks = f"{OUTDIR}/peaks/final/{{caller}}/{{sample}}.unique_501bp_peaks.bed",
        h3k4me1_peaks = "resources/histone_peaks/H3K4me1_{sample}_Combined_Peaks.bed",
        h3k4me3_peaks = "resources/histone_peaks/H3K4me3_{sample}_Combined_Peaks.bed",
        h3k27me3_peaks = "resources/histone_peaks/H3K27me3_{sample}_Combined_Peaks.bed",
        h3k27ac_peaks = "resources/histone_peaks/H3K27ac_{sample}_Combined_Peaks.bed",
        chrom_sizes = config['chrom_sizes']
    output:
        csv = f"{OUTDIR}/peak_filter/{{caller}}/{{sample}}.csv",
        png = f"{OUTDIR}/peak_filter/{{caller}}/{{sample}}.png"
    conda: "../envs/pybedtools.yaml"
    resources: cpus=1, cpus_med2=1, time_min=120, mem_mb=4000, mem_mb_bmm=4000, partition="med2"
    log: "logs/plot_roc/{sample}_{caller}.log"
    script: "../scripts/pyPlotROC.py"