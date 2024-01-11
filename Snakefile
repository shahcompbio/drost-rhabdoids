import os 

samples = ['PMC-01-PDOX1', 'PMC-01-PDOX2']

base_dir = config['base_dir']

rule all:
    input:
        expand(os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam'), sample=samples),

rule extract_reads:
    input:
        bam = os.path.join(base_dir, 'data/{sample}.bam'),
        intervals = ''
    output:
        bam = os.path.join(base_dir, 'data/alignment/extracted/{sample}.svs.bam'),
    shell:
        'PrintReads.sh {input.bam} {input.intervals} {output.bam}'
