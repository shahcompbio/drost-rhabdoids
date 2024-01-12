import os 

samples = ['PMC-01-PDOX1', 'PMC-01-PDOX2']

base_dir = config['base_dir']
tmp_dir = config['tmp_dir']

rule all:
    input:
        os.path.join(base_dir, 'metadata/reference/GRCh38_mm10m.fa.fai'),
        expand(os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam'), sample=samples),
        expand(os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam.bai'), sample=samples),
        expand(os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.fastq.gz'), sample=samples),
        expand(os.path.join(base_dir, 'data/alignment/bam/default/{sample}.svs.bam'), sample=samples),

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

rule reheader_mouse_reference:
    input:
        mouse_fa = os.path.join(base_dir, 'metadata/reference/mm10_build38_mouse.fasta'),
        mouse_fai = os.path.join(base_dir, 'metadata/reference/mm10_build38_mouse.fasta.fai'),
    output:
        mouse_fa = os.path.join(base_dir, 'metadata/reference/mm10m.fa'),
        mouse_fai = os.path.join(base_dir, 'metadata/reference/mm10m.fa.fai'),
    shell:
        "cat {input.mouse_fa} | sed 's/^>/>m/' > {output.mouse_fa} && "
        "cat {input.mouse_fai} | sed 's/^>/>m/' > {output.mouse_fai}"

rule create_joint_reference:
    input:
        human_fa = os.path.join(base_dir, 'metadata/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa'),
        mouse_fa = os.path.join(base_dir, 'metadata/reference/mm10m.fasta'),
    output:
        fa = os.path.join(base_dir, 'metadata/reference/GRCh38_mm10m.fa'),
        fai = os.path.join(base_dir, 'metadata/reference/GRCh38_mm10m.fa.fai'),
    shell:
        "cat {input.human_fa} {input.mouse_fa} > {output.fa} && "
        "samtools faidx {output.fa}"

rule minimap_index:
    input:
        fa = os.path.join(base_dir, 'metadata/reference/GRCh38_mm10m.fa'),
    output:
        mmi = os.path.join(base_dir, 'metadata/reference/GRCh38_mm10m.fa.mmi'),
    singularity:
        "/juno/work/shah/users/chois7/singularity/sif/minimap2__v2.15dfsg-1-deb_cv1.sif",
        #"docker://biocontainers/minimap2:v2.15dfsg-1-deb_cv1",
    threads: 10,
    shell:
        'minimap2 -ax map-ont -t {threads} -I8g -d {output.mmi} {input.fa}'

rule minimap_align:
    input:
        fastq = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.fastq.gz'),
        mmi = os.path.join(base_dir, 'metadata/reference/GRCh38_mm10m.fa.mmi'),
    output:
        sam = os.path.join(base_dir, 'data/alignment/bam/default/{sample}.svs.sam'),
    singularity:
        "/juno/work/shah/users/chois7/singularity/sif/minimap2__v2.15dfsg-1-deb_cv1.sif",
        #"docker://biocontainers/minimap2:v2.15dfsg-1-deb_cv1",
    threads: 20,
    shell:
        'minimap2 -ax map-ont -t {threads} -I8g ' # increase index chunk to 8Gbp
        '{input.mmi} {input.fastq} > {output.sam}'

rule sam_to_bam:
    input:
        sam = os.path.join(base_dir, 'data/alignment/bam/default/{sample}.svs.sam'),
    output:
        bam = os.path.join(base_dir, 'data/alignment/bam/default/{sample}.svs.bam'),
    threads: 6,
    shell:
        'samtools view -b -h -O BAM -@ {threads} -o {output.bam} {input.sam}'
    
#rule nextflow_minimap:
#    input:
#        fastq = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.fastq.gz'),
#    output:
#        bam = 
