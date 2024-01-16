import os 

samples = ['PMC-01-PDOX1', 'PMC-01-PDOX2']

base_dir = config['base_dir']
tmp_dir = config['tmp_dir']
genome_version = config['genome_version']
genome_fas = config['genome_fa']
genome_mmis = config['genome_mmi']
genome_fa = genome_fas[genome_version]
genome_mmi = genome_mmis[genome_version]
result_version = config['result_version']
minimap_I_option = ''
if genome_version == 'GRCh38_mm10m':
    minimap_I_option = '-I8g'


rule all:
    input:
        expand(os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.bam'), sample=samples),
        expand(os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.sorted.bam.bai'), sample=samples),

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
        fa = genome_fa,
    output:
        mmi = genome_mmi,
    singularity:
        "/juno/work/shah/users/chois7/singularity/sif/minimap2__v2.15dfsg-1-deb_cv1.sif",
        #"docker://biocontainers/minimap2:v2.15dfsg-1-deb_cv1",
    threads: 10,
    params: I_option = minimap_I_option,
    shell:
        'minimap2 -ax map-ont -t {threads} {params.I_option} -d {output.mmi} {input.fa}'

rule minimap_align:
    input:
        fastq = os.path.join(base_dir, 'data/alignment/fastq/{sample}_R1_001.fastq.gz'),
        mmi = genome_mmi,
    output:
        sam = os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.sam'),
    singularity:
        "/juno/work/shah/users/chois7/singularity/sif/minimap2__v2.15dfsg-1-deb_cv1.sif",
        #"docker://biocontainers/minimap2:v2.15dfsg-1-deb_cv1",
    threads: 20,
    params: I_option = minimap_I_option,
    shell:
        'minimap2 -ax map-ont -t {threads} {params.I_option} ' # increase index chunk to 8Gbp
        '{input.mmi} {input.fastq} > {output.sam}'

rule sam_to_bam:
    input:
        sam = os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.sam'),
    output:
        bam = os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.bam'),
    threads: 6,
    shell:
        'samtools view -b -h -O BAM -@ {threads} -o {output.bam} {input.sam}'

rule samtools_sort:
    input:
        bam = os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.bam'),
    output:
        bam = os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.sorted.bam'),
        bai = os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.sorted.bam.bai'),
    params:
        prefix = os.path.join(base_dir, f'data/alignment/bam/{result_version}/{{sample}}.svs'),
    threads: 6,
    shell:
        'samtools sort -@ {threads} -o {output.bam} -T {params.prefix} {input.bam} && '
        'samtools index {output.bam}'

