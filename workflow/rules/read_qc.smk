"""
A collection of rules for quality control.
This includes QC metrics of raw, trimmed reads, alignment, and called peaks.
"""

rule fastqc:
    input:
        "resources/raw_reads/{sample}_{read}.fastq.gz"
    output:
        html="results/qc/fastqc/raw/{sample}_{read}_fastqc.html",
        zip="results/qc/fastqc/raw/{sample}_{read}_fastqc.zip"
    conda: "../envs/fastqc.yaml"
    params: out_dir="results/qc/fastqc/raw/", partition="low2"
    log:
        "logs/fastqc/{sample}_{read}.log"
    threads: 1
    resources: mem_mb = 500, time_min = 60, cpus = 1
    wrapper:
        "v1.3.2/bio/fastqc"
        
rule multiqc_raw:
    input:
        expand("results/qc/fastqc/raw/{sample}_{read}_fastqc.zip", sample=pep.sample_table['sample_name'], read = ['R1', 'R2'])
    output:
        "results/qc/fastqc/raw/multiqc_report.html"
    params:
        in_dir="results/qc/fastqc/raw/",
        out_dir="results/qc/fastqc/raw/",
        partition="low2"
    resources: mem_mb = 500, time_min = 60, cpus = 1
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.3.2/bio/multiqc"
     
rule fastp:
    input: sample=expand("results/qc/fastqc/raw/{{sample}}_{read}_fastqc.zip", read = ['R1', 'R2'])
    output:
        trimmed=["results/trimmed_reads/{sample}.1.fastq", "results/trimmed_reads/{sample}.2.fastq"],
        # Unpaired reads separately
        unpaired1="results/trimmed_reads/{sample}.u1.fastq",
        unpaired2="results/trimmed_reads/{sample}.u2.fastq",
        failed="results/trimmed_reads/{sample}.failed.fastq",
        html="results/qc/trimmed/{sample}.fastp.html",
        json="results/qc/trimmed/{sample}.fastp.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 2
    params: partition="bmm"
    resources: time_min=1000, mem_mb=8000, mem_mb_bmm=8000, cpus_bmm=2, cpus=2
    wrapper:
        "v1.3.2/bio/fastp"
        
rule multiqc_trimmed:
    input:
        expand("results/qc/trimmed/{sample}.fastp.json", sample=pep.sample_table['sample_name'])
    output:
        "results/qc/trimmed/multiqc_report.html"
    params:
        in_dir="results/qc/trimmed/",
        out_dir="results/qc/trimmed/",
        partition="bmm"
    log:
        "logs/multiqc_trimmed.log"
    wrapper:
        "v1.3.2/bio/multiqc"