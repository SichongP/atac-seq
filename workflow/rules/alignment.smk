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
    wrapper:
        "v1.3.2/bio/bwa-mem2/index"
        
rule bwa_mem:
# Map reads using bwa-mem2, mark duplicates by samblaster and sort and index by sambamba
    input:
        reads=lambda wildcards: [pep.sample_table.loc[wildcards.sample, 'read1'], pep.sample_table.loc[wildcards.sample, 'read2']],
        idx=multiext("resources/ref/ref", ".fa", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    output:
        bam="results/mapped/{sample}.bam",
        index="results/mapped/{sample}.bam.bai",
    log:
        "logs/bwa_mem2_sambamba/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort_extra="-q",  # Extra args for sambamba.
    threads: 8
    wrapper:
        "v1.3.2/bio/bwa-mem2/mem-samblaster"