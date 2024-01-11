import os 

samples = ['PMC-01-PDOX1', 'PMC-01-PDOX2']

base_dir = config['base_dir']
tmp_dir = config['tmp_dir']

rule all:
    input:
        expand(os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam'), sample=samples),
        expand(os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam.bai'), sample=samples),
        expand(os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.fastq.gz'), sample=samples),

rule extract_reads:
    input:
        bam = os.path.join(base_dir, 'data/{sample}.bam'),
        intervals = os.path.join(base_dir, 'metadata/regions.intervals'),
    output:
        bam = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam'),
    singularity:
        '/juno/work/shah/users/chois7/singularity/sif/gatk4.sif',
        #'docker://broadinstitute/gatk',
    params:
        tmp_dir = tmp_dir,
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
        "PrintReads -I {input.bam} -L {input.intervals} -O {output.bam} "
        "--disable-tool-default-read-filters true"

rule samtools_index:
    input:
        bam = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam'),
    output:
        bai = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam.bai'),
    shell:
        'samtools index {input.bam}'

rule bedtofastq:
    input:
        bam = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam'),
    output:
        fastq = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.fastq.gz'),
    params:
        fastq = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.fastq'),
    shell:
        'bedtools bamtofastq -i {input.bam} -fq {params.fastq} && '
        'gzip {params.fastq}'

