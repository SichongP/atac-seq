rule bwa_mem:
# Map reads using bwa-mem2, mark duplicates by samblaster and sort and index by sambamba
    input:
        reads = ["results/trimmed_reads/{sample}.1.fastq", "results/trimmed_reads/{sample}.2.fastq"],
        idx={config['ref_idx']}
    output:
        bam="results/mapped/{sample}.bam",
        index="results/mapped/{sample}.bam.bai",
    log:
        "logs/bwa_mem2_sambamba/{sample}.log"
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort_extra="-q",
        partition="bmm"
    threads: 8
    resources: cpus=8, cpus_bmm=8, time_min=360, mem_mb=10000, mem_mb_bmm=10000
    wrapper:
        "v1.3.2/bio/bwa-mem2/mem-samblaster"
        
rule filter:
    input:
        bam="results/mapped/{sample}.bam",
        index="results/mapped/{sample}.bam.bai"
    output:
        bam="results/filtered_bam/{sample}.bam",
        index="results/filtered_bam/{sample}.bam.bai"
    log:
        "logs/filter/{sample}.log"
    resources: cpus=4, cpus_bmm=4, time_min=360, mem_mb=10000, mem_mb_bmm=10000
    params: partition="bmm"
    conda: "../envs/samtools.yaml"
    shell:
     """
     samtools idxstats {input.bam} | cut -f 1 | grep -v "chrUn" | grep -v "chrM" | xargs samtools view -bh -@ {resources.cpus} -F 3844 -q 30 {input.bam} > {output.bam}
     samtools index -@ {resources.cpus} {output.bam}
     """
        