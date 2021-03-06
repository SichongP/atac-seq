rule get_ref_genome:
    input: {config['ref']}
    output: "resources/ref/ref.fa"
    shell:
     """
     cp {input} {output}
     """

rule bwa_index:
    input:
        "resources/ref/ref.fa",
    output:
        "resources/ref/ref.0123",
        "resources/ref/ref.amb",
        "resources/ref/ref.ann",
        "resources/ref/ref.bwt.2bit.64",
        "resources/ref/ref.pac",
    log:
        "logs/bwa-mem2_index/ref.log",
    params: prefix="resources/ref/ref"
    resources: cpus=1, cpus_bmm=1, time_min=360, mem_mb=5000, mem_mb_bmm=5000, partition="bmh"
    wrapper:
        "v1.3.2/bio/bwa-mem2/index"
        
rule bwa_mem:
# Map reads using bwa-mem2, mark duplicates by samblaster and sort and index by sambamba
    input:
        reads = [f"{OUTDIR}/trimmed_reads/{{sample}}.1.fastq", f"{OUTDIR}/trimmed_reads/{{sample}}.2.fastq"],
        idx=multiext("resources/ref/ref", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    output:
        bam=f"{OUTDIR}/mapped/{{sample}}.bam",
        index=f"{OUTDIR}/mapped/{{sample}}.bam.bai",
    log: "logs/bwa_mem2_sambamba/{sample}.log"
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort_extra="-q"
    threads: 8
    resources: cpus=8, cpus_bmm=8, time_min=360, mem_mb=10000, mem_mb_bmm=10000, partition="bmm"
    wrapper:
        "v1.3.2/bio/bwa-mem2/mem-samblaster"
        
rule filter:
    input:
        bam = f"{OUTDIR}/mapped/{{sample}}.bam",
        index = f"{OUTDIR}/mapped/{{sample}}.bam.bai"
    output:
        bam = f"{OUTDIR}/filtered_bam/{{sample}}.bam",
        index = f"{OUTDIR}/filtered_bam/{{sample}}.bam.bai"
    log: "logs/filter/{sample}.log"
    resources: cpus=4, cpus_bmm=4, time_min=360, mem_mb=10000, mem_mb_bmm=10000, partition="bmm"
    conda: "../envs/samtools.yaml"
    shell:
     """
     samtools idxstats {input.bam} | cut -f 1 | grep -v "chrUn" | grep -v "chrM" | xargs samtools view -bh -@ {resources.cpus} -F 3844 -q 30 {input.bam} > {output.bam}
     samtools index -@ {resources.cpus} {output.bam}
     """
        
rule shift_reads:
    input:
        bam = f"{OUTDIR}/filtered_bam/{{sample}}.bam",
        index = f"{OUTDIR}/filtered_bam/{{sample}}.bam.bai"
    output:
        bam = f"{OUTDIR}/shifted_bam/{{sample, [\w\d_]+}}.bam",
        index = f"{OUTDIR}/shifted_bam/{{sample, [\w\d_]+}}.bam.bai"
    resources: cpus=4, cpus_bmm=4, time_min=240, mem_mb=8000, mem_mb_bmm=8000, partition="med2"
    conda: "../envs/deeptools.yaml"
    log: "logs/deeptools/{sample}.shift.log"
    shell:
     """
     alignmentSieve -b {input.bam} -o {output.bam}.tmp -p {resources.cpus} --ATACshift
     samtools sort -o {output.bam} -@ {resources.cpus} {output.bam}.tmp
     rm {output.bam}.tmp
     samtools index -@ {resources.cpus} {output.bam}
     """
