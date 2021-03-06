"""
A collection of rules for quality control.
This includes QC metrics of raw, trimmed reads, alignment, and called peaks.
"""

rule fastqc:
    input: lambda w: pep.sample_table.loc[w.sample, w.read]
    output:
        html = f"{OUTDIR}/qc/fastqc/raw/{{sample}}_{{read}}_fastqc.html",
        zip = f"{OUTDIR}/qc/fastqc/raw/{{sample}}_{{read}}_fastqc.zip"
    params: "--quiet"
    log: "logs/fastqc/{sample}_{read}.log"
    threads: 1
    resources: mem_mb = 500, time_min = 60, cpus = 1, partition="low2"
    wrapper:
        "v1.3.2/bio/fastqc"
        
rule multiqc_raw:
    input: expand(f"{OUTDIR}/qc/fastqc/raw/{{sample}}_{{read}}_fastqc.zip", sample=pep.sample_table['sample_name'], read = ['R1', 'R2'])
    output: f"{OUTDIR}/qc/fastqc/raw/multiqc_report.html"
    params:
        in_dir = lambda w, input: os.path.splitext(input[0])[0],
        out_dir = lambda w, output: os.path.splitext(output[0])[0]
    resources: mem_mb = 500, time_min = 60, cpus = 1, partition="low2"
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.3.2/bio/multiqc"
     
rule fastp:
    input: sample=lambda w: pep.sample_table.loc[w.sample, ['R1', 'R2']]
    output:
        trimmed=[f"{OUTDIR}/trimmed_reads/{{sample}}.1.fastq", f"{OUTDIR}/trimmed_reads/{{sample}}.2.fastq"],
        # Unpaired reads separately
        unpaired1=f"{OUTDIR}/trimmed_reads/{{sample}}.u1.fastq",
        unpaired2=f"{OUTDIR}/trimmed_reads/{{sample}}.u2.fastq",
        failed=f"{OUTDIR}/trimmed_reads/{{sample}}.failed.fastq",
        html=f"{OUTDIR}/qc/trimmed/{{sample}}.fastp.html",
        json=f"{OUTDIR}/qc/trimmed/{{sample}}.fastp.json"
    log: "logs/fastp/{sample}.log"
    threads: 2
    resources: time_min=1000, mem_mb=8000, mem_mb_med=8000, cpus_bmm=2, cpus=2, partition="bmm"
    wrapper:
        "v1.3.2/bio/fastp"
        
rule multiqc_trimmed:
    input:
        expand(f"{OUTDIR}/qc/trimmed/{{sample}}.fastp.json", sample=pep.sample_table['sample_name'])
    output: f"{OUTDIR}/qc/trimmed/multiqc_report.html"
    params:
        in_dir = lambda w, input: os.path.splitext(input[0])[0],
        out_dir = lambda w, output: os.path.splitext(output[0])[0]
    log: "logs/multiqc_trimmed.log"
    resources: mem_mb = 500, time_min = 60, cpus = 1, partition="low2"
    wrapper:
        "v1.3.2/bio/multiqc"