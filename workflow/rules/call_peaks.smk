rule make_bed_se:
    # We generate single-ended bed files from paired-end bam files
    # Because both ends of a fragment represent a Tn5 insertion site each
    # And we are specifically looking for Tn5 insertion sites
    input: f"{OUTDIR}/shifted_bam/{{sample}}.bam"
    output: f"{OUTDIR}/bed_se/{{sample}}.bed"
    conda: "../envs/bedtools.yaml"
    resources: cpus=1, cpus_med2=1, time_min=120, mem_mb=4000, mem_mb_bmm=4000, partition="med2"
    log: "logs/make_bed_se/{sample}.log"
    shell:
     """
     bedtools bamtobed -i {input} 2>{log} 1>{output} 
     """