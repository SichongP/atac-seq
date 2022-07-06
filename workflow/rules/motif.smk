rule find_motif:
    input: 
        bed = f"{OUTDIR}/peaks/final/macs3/common.unique_501bp_peaks.bed",
        ref = "resources/ref/ref.fa"
    output: f"{OUTDIR}/motifs/{{name}}/homerResults.html"
    conda: "../envs/homer.yaml"
    params: outdir = lambda w, output: os.path.split(output[0])[0]
    resources: cpus=1, time_min=120, mem_mb=2000, partition="med2"
    shell:
     """
     findMotifsGenome.pl {input.bed} {input.ref} {params.outdir} -size 250
     """