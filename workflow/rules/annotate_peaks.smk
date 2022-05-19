rule make_annotations:
    input:
        gtf = {config['gff']}
    output:
        annotation = f"{OUTDIR}/annotation/annotations.txt",
        final = f"{OUTDIR}/annotation/annotations.final.txt",
        stat = f"{OUTDIR}/annotation/annotations.stat.txt"
    conda: "../envs/homer.yaml"
    resources: cpus=1, time_min=60, mem_mb=2000, partition="low2"
    log: "logs/annotations/make_annotations.log"
    shell:
     """
     parseGTF.pl {input.gtf} ann > {output.annotation}
     assignGenomeAnnotation {output.annotation} {output.annotation} -prioritize {output.final} > {output.stat}
     """

rule annotate_peaks:
    input:
        peaks = f"{OUTDIR}/peaks/final/{{caller}}/{{name}}.unique_501bp_peaks.bed",
        genome = {config['ref']},
        ann = f"{OUTDIR}/annotation/annotations.final.txt"
    output:
        annotated = f"{OUTDIR}/annotation/{{caller}}/{{name}}.bed"
    conda: "../envs/homer.yaml"
    resources: cpus=1, time_min=60, mem_mb=2000, partition="low2"
    log: "logs/annotations/annotate_peaks_{name}_{caller}.log"
    shell:
     """
     annotatePeaks.pl {input.peaks} {input.genome} -ann {input.ann} > {output.annotated}
     """