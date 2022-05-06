rule frag_size_pe:
    input: f"{OUTDIR}/mapped/{{sample}}.bam"
    output:
        fig = f"{OUTDIR}/figures/frag_size/{{sample}}.png",
        stat = f"{OUTDIR}/metrics/frag_size/{{sample}}.tsv",
        raw = f"{OUTDIR}/metrics/frag_size/counts_{{sample}}.txt"
    conda: "../envs/deeptools.yaml"
    resources: cpus=8, cpus_bmm=8, mem_mb=4000, mem_mb_bmm=4000, time_min=60, partition='bmm'
    log: "logs/frag_size/{sample}.log"
    shell:
     """
     bamPEFragmentSize -hist {output.fig} -p {resources.cpus} -T "Fragment Size Distribution" -bs 10000 --maxFragmentLength 2000 -n 100000 --table {output.stat} --outRawFragmentLengths {output.raw} -b {input}
     """

rule sambamba_flagstat:
    input: f"{OUTDIR}/mapped/{{sample}}.bam"
    output: f"{OUTDIR}/mapped/stats/{{sample}}.flagstats.txt"
    log: "logs/sambamba-flagstat/{sample}.log"
    params: extra = ""
    resources: cpus=1, cpus_bmm=1, mem_mb=3000, mem_mb_bmm=3000, time_min=300, partition='bmm'
    threads: 1
    wrapper:
        "v1.3.2/bio/sambamba/flagstat"

rule tss_from_gff3:
    input: {config['gff']}
    output: "resources/annotation/tss.bed"
    resources: cpus=1, cpus_bmm=1, mem_mb=3000, mem_mb_bmm=3000, time_min=90, partition='low2'
    log: "logs/tss_from_gff3/tss_from_gff3.txt"
    shell:
     """
     cat {input} | awk '{{if ($3=="mRNA") print $0}}' | awk '{{if ($7=="+") {{print $1"\t"$4-1"\t"$4"\t"$7}} else {{print $1"\t"$5-1"\t"$5"\t"$7}}}}' 1> {output} 2> {log}
     """

rule nuclear_mit_ratio:
    input:
        bam = f"{OUTDIR}/mapped/{{sample}}.bam",
        filtered = f"{OUTDIR}/filtered_bam/{{sample}}.bam"
    output: f"{OUTDIR}/mapped/stats/{{sample}}.nucl.mit.ratio"
    conda: "../envs/sambamba.yaml"
    resources: cpus=1, mem_mb = 3000, cpus_bmm=1, mem_mb_bmm=3000, time_min=180, partition='bmm'
    log: "logs/nuclear_mit_ratio/{sample}.log"
    shell:
     """
     chrM_count=$(sambamba view -c -t {resources.cpus} -F "proper_pair and not duplicate and not (supplementary or secondary_alignment)" {input.bam} chrM)
     nucl_count=$(sambamba view -c -t {resources.cpus} {input.filtered})
     echo $chrM_count"\t"$nucl_count > {output}
     """
     
rule tss_enrichment:
    input:
        bam = f"{OUTDIR}/filtered_bam/{{sample}}.bam",
        bed = "resources/annotation/tss.bed"
    output: f"{OUTDIR}/tss_enricment/{{sample}}.TssEnrichment"
    conda: "../envs/tss_enrichment.yaml"
    resources: cpus=1, mem_mb = 5000, cpus_bmm=1, mem_mb_bmm=5000, time_min = 360, partition='bmm'
    log: "logs/nuclear_mit_ratio/{sample}.log"
    shell:
     """
     python workflow/scripts/pyTssEnrichment.py -a {input.bam} -b {input.bed} -o {output} -p ends -z -v -s 4 -u 2000 -d 2000
     """
