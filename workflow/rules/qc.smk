# Illumina QC
# Chunyu Zhao

import subprocess
from collections import OrderedDict
import pandas

# Quality-control reads
TARGET_CLEAN = expand(
    str(QC_FP/'02_bbduk'/'{sample}_{rp}.fastq.gz'),
    sample = Samples.keys(), rp = Pairs)


rule all_qc:
    """Runs trimmomatic and fastqc on all input files."""
    input:
        TARGET_CLEAN


ruleorder: trimmomatic_paired > trimmomatic_unpaired
ruleorder: entropy_filter_paired > entropy_filter_unpaired
ruleorder: align_to_host_paired > align_to_host_unpaired
ruleorder: filter_host_reads > filter_host_reads_unpaired


rule trimmomatic_unpaired:
    input:
        lambda wildcards: Samples[wildcards.sample]['1']
    output:
        str(QC_FP/'01_trimmomatic'/'{sample}_1.fastq.gz')
    log:
        str(QC_FP/'log'/'trimmomatic'/'{sample}.out')
    params:
        adapter = CONDA_PATH / f"share/trimmomatic-0.39-2/adapters/{ADAPTER_FASTA}",
        jar = str(CONDA_PATH/'share/trimmomatic-0.39-2/trimmomatic.jar')
    threads: 8
    shell:
        """
        java -Xmx1028m -jar {params.jar} \
        SE -threads {threads} -phred33 \
        {input} {output} \
        ILLUMINACLIP:{params.adapter}:2:30:10:8:true \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 \
        > >(tee -a {log}) 2> >(tee -a {log} >&2)
        """


rule trimmomatic_paired:
    input:
        r1 = lambda wildcards: Samples[wildcards.sample]['1'],
        r2 = lambda wildcards: Samples[wildcards.sample]['2']
    output:
        pair_r1 = str(QC_FP/'01_trimmomatic'/'{sample}_1.fastq.gz'),
        pair_r2 = str(QC_FP/'01_trimmomatic'/'{sample}_2.fastq.gz'),
        unpair_r1 = temp(str(QC_FP/'01_trimmomatic'/'unpaired'/'{sample}_1_unpaired.fastq.gz')),
        unpair_r2 = temp(str(QC_FP/'01_trimmomatic'/'unpaired'/'{sample}_2_unpaired.fastq.gz'))
    log:
        str(QC_FP/'log'/'trimmomatic'/'{sample}.out'),
    params:
        adapter = CONDA_PATH / f"share/trimmomatic-0.39-2/adapters/{ADAPTER_FASTA}",
        jar = str(CONDA_PATH/'share/trimmomatic-0.39-2/trimmomatic.jar')
    threads: 8
    shell:
        """
        java -Xmx1028m -jar {params.jar} \
        PE -threads {threads} -phred33 \
        {input.r1} {input.r2} \
        {output.pair_r1} {output.unpair_r1} \
        {output.pair_r2} {output.unpair_r2} \
        ILLUMINACLIP:{params.adapter}:2:30:10:8:true \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 \
        > >(tee -a {log}) 2> >(tee -a {log} >&2)
        """


rule entropy_filter_unpaired:
    input:
        r1 = str(QC_FP/'01_trimmomatic'/'{sample}_1.fastq.gz'),
    output:
        r1 = str(QC_FP/'02_bbduk'/'{sample}_1.fastq.gz'),
    log:
        str(QC_FP/'log'/'bbduk'/'{sample}.log'),
    threads: 8
    shell:
        """
        bbduk.sh -Xmx1028m in={input.r1} out={output.r1} \
        k=23 hdist=1 entropy=0.5 entropywindow=50 entropyk=5 \
        minlen=50 threads={threads} 2>&1 | tee {log}
        """


rule entropy_filter_paired:
    input:
        r1 = str(QC_FP/'01_trimmomatic'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'01_trimmomatic'/'{sample}_2.fastq.gz'),
    output:
        r1 = str(QC_FP/'02_bbduk'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'02_bbduk'/'{sample}_2.fastq.gz'),
    log:
        str(QC_FP/'log'/'bbduk'/'{sample}.log'),
    threads: 8
    shell:
        """
        bbduk.sh -Xmx1028m in={input.r1} in2={input.r2} \
        out={output.r1} out2={output.r2} \
        k=23 hdist=1 entropy=0.5 entropywindow=50 entropyk=5 \
        minlen=50 threads={threads} 2>&1 | tee {log}
        """


rule align_to_host_unpaired:
    input:
        r1 = str(QC_FP/'02_bbduk'/'{sample}_1.fastq.gz'),
        index = str(Cfg['qc']['host_fp']/'{host}.rev.2.bt2')
    output:
        str(QC_FP/'03_decontam'/'ids'/'{host}'/'{sample}')
    conda:
        "../envs/decontam.yml"
    threads: 8
    params:
        index_fp = str(Cfg['qc']['host_fp'])
    shell:
        """
        ## turn off bash strict mode
        set +o pipefail

        bowtie2 --no-unal --local --very-sensitive-local \
        -X 500.0 -x {params.index_fp}/{wildcards.host} \
        --threads {threads} -q -U {input.r1} | \
        samtools view --threads {threads} -S - | awk '{{print $1}}' | sort | uniq > {output}
        """


rule align_to_host_paired:
    input:
        r1 = str(QC_FP/'02_bbduk'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'02_bbduk'/'{sample}_2.fastq.gz'),
        index = str(Cfg['qc']['host_fp']/'{host}.rev.2.bt2')
    output:
        str(QC_FP/'03_decontam'/'ids'/'{host}'/'{sample}')
    conda:
        "../envs/decontam.yml"
    threads: 8
    params:
        index_fp = str(Cfg['qc']['host_fp'])
    shell:
        """
        ## turn off bash strict mode
        set +o pipefail

        bowtie2 --no-unal --local --very-sensitive-local \
        -X 500.0 -x {params.index_fp}/{wildcards.host} \
        --threads {threads} -q -1 {input.r1} -2 {input.r2} |
        samtools view --threads {threads} -S - | awk '{{print $1}}' | sort | uniq > {output}
        """


rule gather_host_ids:
    input:
        expand(str(QC_FP/'03_decontam'/'ids'/'{host}'/'{{sample}}'),
            host=HostGenomes.keys())
    output:
        str(QC_FP/'03_decontam'/'hostreads'/'{sample}')
    shell:
        """
        cat {input} | sort | uniq > {output}
        """


rule filter_host_reads_unpaired:
    input:
        ids = str(QC_FP/'03_decontam'/'hostreads'/'{sample}'),
        r1 = str(QC_FP/'01_trimmomatic'/'{sample}_1.fastq.gz'),
    output:
        r1 = str(QC_FP/'03_decontam'/'{sample}_1.fastq.gz'),
    conda:
        "../envs/decontam.yml"
    shell:
        """
        seqkit grep -v -f {input.ids} {input.r1} -o {output.r1}
        """


rule filter_host_reads:
    input:
        ids = str(QC_FP/'03_decontam'/'hostreads'/'{sample}'),
        r1 = str(QC_FP/'01_trimmomatic'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'01_trimmomatic'/'{sample}_2.fastq.gz'),
    output:
        r1 = str(QC_FP/'03_decontam'/'{sample}_1.fastq.gz'),
        r2 = str(QC_FP/'03_decontam'/'{sample}_2.fastq.gz')
    conda:
        "../envs/decontam.yml"
    shell:
        """
        seqkit grep -v -f {input.ids} {input.r1} -o {output.r1}
        seqkit grep -v -f {input.ids} {input.r2} -o {output.r2}
        """
