#!/bin/bash

[ $# -ne 3 ] && { echo -e "\nUsage: $(basename $0) <in.bam> <in.intervals> <out.bam>\n" 1>&2 && exit 1; }
in_bam=$1
in_intervals=$2
out_bam=$3
[ ! -f $in_bam ] && { echo "[ERROR] $in_bam does not exist" 1>&2 && exit 1; }
[ ! -f $in_intervals ] && { echo "[ERROR] $in_intervals does not exist" 1>&2 && exit 1; }

image="docker://broadinstitute/gatk"

tmp_dir="/juno/work/shah/users/chois7/tmp"
[ ! -d $tmp_dir ] && { echo "[ERROR] $tmp_dir does not exist" 1>&2 && exit 1; }

cmd="singularity run -B /juno -B /home -B /ifs $image"
cmd="$cmd gatk --java-options '-Xmx4g -Djava.io.tmpdir=${tmp_dir}'"
cmd="$cmd PrintReads"
cmd="$cmd -I ${in_bam}"
cmd="$cmd -L ${in_intervals}"
cmd="$cmd -O ${out_bam}"
cmd="$cmd --disable-tool-default-read-filters true;"

echo $cmd
eval $cmd

src_bai=${out_bam%%\.bam}.bai
dst_bai=${out_bam}.bai
[[ ! -f $src_bai ]] && { echo "[ERROR] $src_bai does not exist" 1>&2 && exit 1; }
[[ -f $src_bai ]] && [[ ! -f $dst_bai ]] && { cmd="mv $src_bai $dst_bai"; echo $cmd; eval $cmd; }
[[ ! -f $src_bai ]] && [[ ! -f $dst_bai ]] && { cmd="samtools $in_bam"; echo $cmd; eval $cmd; }

