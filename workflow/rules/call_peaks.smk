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
     
rule macs3_merged:
# Here we shift reads toward 5' direction by 75bp and then extend to a fixed 150bp fragment towards 3' direction
# This centers our read (which is now considered as fragment by peak caller) around the Tn5 cutting sites and smooths the signal with a 75bp window on both sides
    input: f"{OUTDIR}/bed_se/{{sample}}.bed"
    output: f"{OUTDIR}/peaks/macs3/{{sample}}/{{sample}}_peaks.narrowPeak"
    conda: "../envs/macs3.yaml"
    params: outdir = lambda w, output: os.path.split(output[0])[0]
    log: "logs/macs3/{sample}.log"
    resources: cpus = 1, time_min=60, mem_mb=4000, cpus_bmm=1, mem_mb_bmm=4000, partition = 'med2'
    shell:
     """
     macs3 callpeak -p 0.01 -t {input} -f BED -g {config[effective_size]} --outdir {params.outdir} -n {wildcards.sample} -B --shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all
     """

rule sort_by_name:
    input: f"{OUTDIR}/shifted_bam/{{sample}}.bam"
    output: f"{OUTDIR}/shifted_bam/name_sorted/{{sample}}.bam"
    params: "-n"
    log: "logs/sort_by_name/{sample}.log"
    threads: 2
    wrapper:
        "v1.4.0/bio/sambamba/sort"

rule genrich:
    input: lambda w: expand(f"{OUTDIR}/shifted_bam/name_sorted/{{sample}}.bam", sample = pep.sample_table.loc[(pep.sample_table['condition'] == w.condition) & (pep.sample_table['sex'] == w.sex)]['sample_name'])
    output:
        peak=f"{OUTDIR}/peaks/genrich/{{condition}}_{{sex}}/{{condition}}_{{sex}}.narrowPeak",
        flog=f"{OUTDIR}/peaks/genrich/{{condition}}_{{sex}}/{{condition}}_{{sex}}.f.bedgraph",
        klog=f"{OUTDIR}/peaks/genrich/{{condition}}_{{sex}}/{{condition}}_{{sex}}.k.bedgraph"
    conda: "../envs/genrich.yaml"
    params: outdir = lambda w, output: os.path.split(output[0])[0]
    resources: cpus = 1, time_min=60, mem_mb=4000, cpus_bmm=1, mem_mb_bmm=4000, partition = 'med2'
    shell:
     """
     chr_exl=$(samtools view -H {input} | grep -E "chrUn|chrM" | cut -d ':' -f 2 | sed -z 's/\\n/,/g')
     genrich -t \"{input}\" -o {output.peak} -f {output.flog} -k {output.klog} -e $chr_exl -j -D -d 150 -p 0.01 -g 50
     """