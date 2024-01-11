---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.5
  kernelspec:
    display_name: Python 3.7.12 64-bit (conda)
    language: python
    name: python3712jvsc74a57bd0dc0fe05456373cce17991fe9e2e9264df4f4d1e972d77814a9700f21c9e7a8e2
---

# Settings

```python
%load_ext autoreload
%autoreload 2
```

```python
import os
```

```python
import pandas as pd
import numpy as np
```

```python
from collections import defaultdict
from tqdm import tqdm
```

```python
import matplotlib.pyplot as plt
```

```python
import re
```

```python
from venny4py.venny4py import venny4py
```

```python
from matplotlib_venn import venn2, venn3
```

```python
import pysam
```

```python
import logging
logger = logging.getLogger()
```

```python
logger.setLevel(logging.INFO)
```

# Data

```python
vcf_path = '/juno/work/shah/users/chois7/tickets/drostpmc/results/survivor/union.lenient.vcf'
vcf_cols = ['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'FILTER', 'info', 'format']
sample_cols = ['PMC-01-ORG1', 'PMC-01-PDOX1', 'PMC-01-PDOX2', 'PMC-01-T1']
vcf_cols += sample_cols
```

```python
vcf = pd.read_table(vcf_path, names=vcf_cols, comment='#')
```

# Proc

```python
vcf[sample_cols].iloc[10].values
```

```python
df = pd.DataFrame()
pattern = '[^:]+:([^:]+):[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+'
col_map = {'NaN':0, 'NA':1}
for sample in sample_cols:
    match = vcf[sample].str.extract(pattern).squeeze()
    match = match.replace(col_map)
    df[sample] = match
```

```python
vcf['ID'].str.replace('_[12]', '', regex=True).value_counts()
```

```python
_tmp = vcf['alt'].str.extract('([^\[^\]]+):([^\[^\]]+)')
_tmp.columns = ['chrom2', 'pos2']
vcf = vcf.join(_tmp)
```

```python
chroms = ['chr'+str(x) for x in range(1, 22+1)] + ['chrX', 'chrY']
```

```python
vcf = vcf[
    vcf['chrom'].isin(chroms) &
    vcf['chrom2'].isin(chroms)
]
```

```python
samples = ['PMC-01-ORG1', 'PMC-01-PDOX1', 'PMC-01-PDOX2', 'PMC-01-T1',]
```

```python
vcf.loc[2].values
```

## get svs set

```python
svs = defaultdict(set)
for _, row in tqdm(vcf.iterrows(), total=vcf.shape[0]):
    chrom1 = row['chrom']
    chrom2 = row['chrom2']
    pos1 = int(row['pos'])
    pos2 = int(row['pos2'])
    sv_coords = f'{chrom1}:{pos1}-{chrom2}:{pos2}'
    fmt = row['format'].split(':')
    for sample in samples:
        gts = row[sample].split(':')
        gts = dict(zip(fmt, gts))
        if gts['CO'] != 'NAN':
            svs[sample].add(sv_coords)
```

```python
bams = {
    'PMC-01-N1': pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-N1.bam', 'rb'),
    'PMC-01-T1': pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-T1.bam', 'rb'),
    'PMC-01-ORG1': pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-ORG1.bam', 'rb'),
    'PMC-01-PDOX1': pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-PDOX1.bam', 'rb'),
    'PMC-01-PDOX2': pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-PDOX2.bam', 'rb'),
}
```

```python
list(svs['PMC-01-PDOX2'])[:10]
```

## extract comprehensive bam stats


### modules

```python
def extract_stats_from_reads(reads):
    """
    Create a pandas DataFrame from a list of reads, filtering for common tags.

    This function takes a list of read objects and extracts specific tags from each read.
    It then creates a pandas DataFrame where each row corresponds to a read and columns
    correspond to the common tags found in these reads.

    Parameters:
    reads (list): A list of pysam AlignedSegment object

    Returns:
    pandas.DataFrame: A DataFrame where each row represents a read and each column represents
                      one of the common tags. The DataFrame is indexed by the 'qname' of each read.
      
    """
    common_tags = set(['NM', 'ms', 'AS', 'nn', 'cm', 's1', 'de', 'r1'])
    cigar_order = ['match', 'ins', 'del', 'skip', 'softclip', 'hardclip', 
                   'pad', 'qual', 'diff', 'back', 'nmtag']
    cigar_result_ix = 1 # 0: number of bases for each cigar, 1: number of blocks
    qnames = [read.qname for read in reads]
    stats = pd.DataFrame(
        [{k:v for (k,v) in read.tags if k in common_tags} for read in reads],
        index=qnames,
    )
    stats['mapq'] = [read.mapq for read in reads]
    cigars = pd.DataFrame( 
        [dict(zip(cigar_order, read.get_cigar_stats()[cigar_result_ix])) for read in reads],
        index=qnames,
    )
    stats = stats.join(cigars)
    if 'nmtag' in stats.columns: 
        stats = stats.drop('nmtag', axis=1)
    
    return stats
```

```python
def get_svstats(bams, samples, margin=1000):
    nbam = bams['PMC-01-N1']
    likely_tp = svs['PMC-01-T1'] | svs['PMC-01-ORG1']
    likely_fp = ((svs['PMC-01-PDOX1'] - svs['PMC-01-PDOX2']) - likely_tp) | ((svs['PMC-01-PDOX2'] - svs['PMC-01-PDOX1']) - likely_tp)

    cols = ['NM', 'ms', 'AS', 'nn', 'cm', 's1', 'de', 'mapq', 'match', 'ins',
            'del', 'skip', 'softclip', 'hardclip', 'pad', 'qual', 'diff', 'back']
    cols = [s+'_t' for s in cols] + [s+'_n' for s in cols] + ['TP']
    svstats = pd.DataFrame(columns=cols)
    for sample in samples:
        ssvs = svs[sample] # sample svs
        tbam = bams[sample]

        for region in tqdm(list(ssvs)):
            flt_tag = 0
            flt_tag += (region in valid_tp)
            flt_tag -= (region in likely_fp)

            reg1, reg2 = region.split('-')
            chrom1, pos1 = reg1.split(':')
            chrom2, pos2 = reg2.split(':')
            pos1, pos2 = int(pos1), int(pos2)

            for chrom, pos in (chrom1, pos1), (chrom2, pos2):
                nreads = [read for read in nbam.fetch(chrom, pos-margin, pos+margin)]
                if len(nreads) == 0: continue
                nstats = extract_stats_from_reads(nreads)
                nstatsdf = nstats.mean().to_frame().T
                treads = [read for read in tbam.fetch(chrom, pos-margin, pos+margin)]
                if len(treads) == 0: continue
                tstats = extract_stats_from_reads(treads)
                tstatsdf = tstats.mean().to_frame().T
                field = tstatsdf.join(nstatsdf, lsuffix="_t", rsuffix="_n").squeeze().tolist()
                field += [flt_tag]
                svstats.loc[svstats.shape[0]] = field
                # raise ValueError
            # break
```

```python
def get_svsupporting_sv_stats(bams, samples, margin=100):
    nbam = bams['PMC-01-N1']
    likely_tp = svs['PMC-01-T1'] | svs['PMC-01-ORG1']
    likely_fp = ((svs['PMC-01-PDOX1'] - svs['PMC-01-PDOX2']) - likely_tp) | ((svs['PMC-01-PDOX2'] - svs['PMC-01-PDOX1']) - likely_tp)

    cols = ['NM', 'ms', 'AS', 'nn', 'cm', 's1', 'de', 'mapq', 'match', 'ins',
            'del', 'skip', 'softclip', 'hardclip', 'pad', 'qual', 'diff', 'back']
    cols = [s+'_t' for s in cols] + [s+'_n' for s in cols] + ['TP']
    svstats = pd.DataFrame(columns=cols)
    for sample in samples:
        ssvs = svs[sample] # sample svs
        tbam = bams[sample]

        for region in tqdm(list(ssvs)):
            flt_tag = 0
            flt_tag += (region in valid_tp)
            flt_tag -= (region in likely_fp)

            reg1, reg2 = region.split('-')
            chrom1, pos1 = reg1.split(':')
            chrom2, pos2 = reg2.split(':')
            pos1, pos2 = int(pos1), int(pos2)
            nreads = [read for read in nbam.fetch(chrom1, pos1-margin, pos1+margin)]
            if len(nreads) == 0: continue
            treads = [read for read in tbam.fetch(chrom1, pos1-margin, pos1+margin)]
            treads = [r for r in get_sv_supporting_reads(treads, chrom1, pos1, chrom2, pos2)]
            if len(treads) == 0: continue
            
            nstats = extract_stats_from_reads(nreads)
            nstatsdf = nstats.mean().to_frame().T
            
            tstats = extract_stats_from_reads(treads)
            tstatsdf = tstats.mean().to_frame().T
            field = tstatsdf.join(nstatsdf, lsuffix="_t", rsuffix="_n").squeeze().tolist()
            field += [flt_tag]
            svstats.loc[svstats.shape[0]] = field
    return svstats
```

### valid TP

```python
valid_tp = {
    "chr16:72806711-chr16:72685741",
    "chr16:6512983-chr16:6645209",
    "chr16:6489642-chr16:6420750",
    "chr16:6366747-chr16:6328918",
    "chr13:86179166-chr13:85396707",
    "chr10:131885267-chr10:131884020",
    "chr10:66681858-chr10:66597638",
    "chr7:40600415-chr7:40643054",
    "chr4:90909961-chr4:90596261",
    "chr4:80131039-chr4:80201094",
    "chr3:79339417-chr3:79170916",
    "chr3:79119950-chr3:79007918",
    "chr2:96071256-chr2:96072063",
}
```

probably I should add stats for only those reads supporting the variant...?


## extract cigars and see SV evidences

```python
def get_secondary_end_pos(chrom, start, strand, cigarstr):
    ref_consumers = 'MDN=X'
    end_pos = int(start)
    for pattern in re.finditer('(\d+)([MIDNSHPX=])', cigarstr):
        assert pattern, 'CIGAR string pattern mismatch'
        clen, cstr = pattern.groups()
        clen = int(clen)
        if cstr in ref_consumers:
            end_pos += clen
    return end_pos
```

#### unit-test: get_secondary_end_pos


- 6ee83486-5f36-43e6-afe7-509f421ddc73 chr1,152379398,+,6560S280M3I9S,60,6;
- f7242fda-e52e-4acb-889b-1da13da63061 chr16,6643764,-,7S4005M114D4187S,60,407;
- f7242fda-e52e-4acb-889b-1da13da63061 chr16,6643764,+,42S4151M10D4006S,60,58;
- 8d9521b1-2ca3-403f-ba51-7bbf9f83892d chr16,6507185,+,37S5798M3I309S,60,86;
- 81ac77f6-a95d-4626-be2c-c023891953d1 chr16,6512037,+,42S946M2000S,60,8;
- 162cdb2e-f7b8-4039-b803-d6fe16811c7a chr16,6498279,-,9S14700M4D7019S,60,297;
- eebf0d9e-e5d6-4770-8e5e-5e233bdad613 chr16,6504535,-,10S8418M30D320S,60,195;
- 673cfc41-5867-442f-93ac-18fd0474a155 chr16,6507185,-,7S5802M3I329S,60,112;

```python
get_secondary_end_pos('chr16', 6498279, '-', '9S14700M4D7019S') == 6512983
```

```python
get_secondary_end_pos('chr16', 6512037, '+', '42S946M2000S') == 6512983
```

```python
get_secondary_end_pos('chr16', 6507185, '+', '37S5798M3I309S') == 6512983
```

## get sv supporting reads

```python
class FilterConfigs:
    def __init__(self, tumor_wobble, tumor_neighbor, normal_wobble):
        self.tumor_wobble = tumor_wobble
        self.tumor_neighbor = tumor_neighbor
        self.normal_wobble = normal_wobble
        
flt = FilterConfigs(
    tumor_wobble = 10,
    tumor_neighbor = 50,
    normal_wobble = 100,
)
```

```python
def get_tag_read(read, tags, read_nickname, chrom2, pos2, wobble, yielded, debug):
    """Iterate SA coordinates and return read if coordinate within wobble for breakpoint 2
    """
    logger.setLevel(logging.INFO)
    if debug: logger.setLevel(logging.DEBUG)
    for sa_tag in tags['SA'].rstrip(';').split(';'):
        tag_chrom2, tag_start2, tag_strand2, tag_cigarstr, tag_mapq, _ = sa_tag.split(',')
        tag_start2 = int(tag_start2)
        tag_end2 = get_secondary_end_pos(tag_chrom2, tag_start2, tag_strand2, tag_cigarstr)
        if tag_chrom2 != chrom2:
            logging.debug(f'[x-3-2: {read_nickname}] tag_chrom2:{tag_chrom2} != chrom2:{chrom2}')
        else:
            if abs(tag_start2 - pos2) < wobble: # ori of opposite side of clip is '-'
                if read.qname not in yielded:
                    logging.debug(f'[o-3-2(-): {read_nickname}] |(tag_start2-pos2):{tag_start2}-{pos2}| < {wobble} => yield {read_nickname}')
                    yielded.add(read.qname)
                    return read
            elif abs(tag_end2 - pos2) < wobble: # ori of opposite side of clip is '+'
                if read.qname not in yielded:
                    logging.debug(f'[o-3-2(+): {read_nickname}] |(tag_start2-pos2):{tag_end2}-{pos2}| < {wobble} => yield {read_nickname}')
                    yielded.add(read.qname)
                    return read
    return None
```

```python
def get_sv_supporting_reads(reads, chrom1, pos1, chrom2, pos2, wobble=10, min_var_size=100,
                            debug=False):   
    logger.setLevel(logging.INFO)
    if debug: logger.setLevel(logging.DEBUG)
    cigar_order = 'MIDNSHP=X'
    query_consumers = {cigar_order.index(x) for x in ['M', 'I', 'S', '=', 'X']}
    ref_consumers = {cigar_order.index(x) for x in ['M', 'D', 'N', '=', 'X']}
    yielded = set()
    for read in reads:
        read_nickname = read.qname[:8] # shortened version of qname
        cigars = read.cigar
        n_cigars = len(cigars)
        brk_pos = read.reference_start
        
        for cix, (cnum, clen) in enumerate(cigars):
            cstr = cigar_order[cnum]
            # logging.debug(f'cix:{cix} {cstr} length {clen}')
            # if cnum in query_consumer:

            if abs(brk_pos - pos1) < wobble and clen > min_var_size:
                logging.debug(f'[init] cstr:{cstr} clen:{clen}>{min_var_size}; |{brk_pos} - {pos1}| = {abs(brk_pos - pos1)} < {wobble}')
                
                if abs(brk_pos+clen - pos2) < wobble:
                    if read.qname not in yielded:
                        logging.debug(f'[o-3-1: {read_nickname}] |brk_pos+clen {brk_pos+clen} - pos2 {pos2}| = {abs(brk_pos+clen - pos2)} < {wobble} => yield {read_nickname}')
                        yielded.add(read.qname)
                        yield read
                else:
                    logging.debug(f'[x-3-1: {read_nickname}] |brk_pos+clen {brk_pos+clen} - pos2 {pos2}| = {abs(brk_pos+clen - pos2)} >= {wobble}')
                    
                    tags = {k:v for (k,v) in read.tags}
                    if 'SA' in tags:
                        tag_read = get_tag_read(read, tags, read_nickname, chrom2, pos2, wobble, yielded, debug)
                        if tag_read: yield read
                    else: 
                        logging.debug(f'[x-3-3: {read_nickname}] "SA" not in tags: {list(tags.keys())}')
            elif brk_pos - pos1 >= wobble:
                logging.debug(f'[x-2-1: {read_nickname}] {brk_pos} - {pos1} >= {wobble}')
                break
            else:
                # logging.debug(f'[x-2-2] |{brk_pos} - {pos1}| = {brk_pos - pos1} >= {wobble}')
                pass

            if cnum in ref_consumers:
                brk_pos += clen
                # logging.debug(f'[x-1] {cstr} in ref_consumers; added to {brk_pos}')
```

```python
def get_neighbor_supporting_reads(reads, chrom1, pos1, chrom2, pos2, wobble=10, wobble_min_size=200,
                            debug=False):   
    logger.setLevel(logging.INFO)
    if debug: logger.setLevel(logging.DEBUG)
    cigar_order = 'MIDNSHP=X'
    query_consumers = {cigar_order.index(x) for x in ['M', 'I', 'S', '=', 'X']}
    ref_consumers = {cigar_order.index(x) for x in ['M', 'D', 'N', '=', 'X']}
    var_cigars = {'I', 'D', 'S', 'N'}
    yielded = set()
    for read in reads:
        read_nickname = read.qname[:8] # shortened version of qname
        cigars = read.cigar
        n_cigars = len(cigars)
        brk_pos = read.reference_start
        tags = {k:v for (k,v) in read.tags}
        
        for cix, (cnum, clen) in enumerate(cigars):
            cstr = cigar_order[cnum]

            if abs(brk_pos - pos1) < wobble and cstr in var_cigars and clen > wobble_min_size:
                logging.debug(f'[neighbor {read_nickname}] cstr:{cstr} clen:{clen} |{brk_pos} - {pos1}| = {abs(brk_pos - pos1)} < {wobble}')
                
                if abs(brk_pos+clen - pos2) >= wobble:
                    if read.qname not in yielded:
                        logging.debug(f'[wo-3-1: {read_nickname}] clen:{clen} & |brk_pos+clen {brk_pos+clen} - pos2 {pos2}| = {abs(brk_pos+clen - pos2)} >= {wobble} => yield {read_nickname}')
                        yielded.add(read.qname)
                        yield read
                else:
                    logging.debug(f'[wo pos2 within wobble: {read_nickname}] |brk_pos+clen {brk_pos+clen} - pos2 {pos2}| = {abs(brk_pos+clen - pos2)} < {wobble}')
                    pass # dist1 < wobble and dist2 < wobble -> NOT neighbor
                # logging.debug(f'[neighbor diffchrom {chrom1} != {chrom2}')
                
                if 'SA' in tags:
                    tag_read = get_tag_read(read, tags, read_nickname, chrom2, pos2, wobble, yielded, debug)
                    
                    if not tag_read: 
                        if read.qname not in yielded:
                            logging.debug(f'[wo-3-2: No non-neighbor tag_read returned: {tag_read} => yield {read_nickname}')
                            yielded.add(read.qname)
                            yield read
                else: 
                    logging.debug(f'[wx-3-3: {read_nickname}] "SA" not in tags: {list(tags.keys())}')
                    pass

                if brk_pos - pos1 >= wobble:
                    logging.debug(f'[wx-2-1: {read_nickname}] {brk_pos} - {pos1} >= {wobble}')
                    break
                else:
                    # logging.debug(f'[wx-2-2] |{brk_pos} - {pos1}| = {brk_pos - pos1} >= {wobble}')
                    pass

            if cnum in ref_consumers:
                brk_pos += clen
                # logging.debug(f'[wx-1] {cstr} in ref_consumers; added to {brk_pos}')
```

#### neighbor-test: PMC-01-PDOX1 chr19:54964908 chr19:54965999 (+pack)

```python
test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-PDOX1.bam'
test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
chrom1, pos1, chrom2, pos2 = 'chr19', 54964908, 'chr19', 54965999
margin = 100
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_twreads1 = get_neighbor_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, wobble=10, wobble_min_size=200, debug=True) # t wiggle 1
sv_twreads2 = get_neighbor_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, wobble=10, wobble_min_size=200, debug=True) # t wiggle 2
sv_twreads1 = [read for read in sv_twreads1]
sv_twreads2 = [read for read in sv_twreads2]
```

```python
[r.qname for r in sv_twreads1]
```

```python
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=True)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False)
sv_treads1 = [read for read in sv_treads1]
sv_treads2 = [read for read in sv_treads2]
```

```python
[r.qname for r in sv_treads1]
```

#### neighbor-test: PMC-01-PDOX1 chr4:58878793-chr6:31427491 (+pack)

```python
test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-PDOX1.bam'
test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
chrom1, pos1, chrom2, pos2 = 'chr4', 58878793, 'chr6', 31427491
margin = 100
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=False)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False)
sv_twreads1 = get_neighbor_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=False) # t wiggle 1
sv_twreads2 = get_neighbor_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False) # t wiggle 2
sv_treads1 = [read for read in sv_treads1]
sv_treads2 = [read for read in sv_treads2]
sv_twreads1 = [read for read in sv_twreads1]
sv_twreads2 = [read for read in sv_twreads2]
```

```python
sv_twread_qnames1 = {
    '0cfd0362-11aa-41b3-8be6-a21f02b877ff',
    '6dd5f766-8c89-472c-831e-277324ce5d0c',
    'facad054-a0a5-4aeb-aeed-bcd0e12a51db'
}
```

```python
set([read.qname for read in sv_twreads1]) == sv_twread_qnames1
```

#### bigger wobble test: PMC-01-PDOX1 and N1 chr20:29045066 chr21:5216246

```python
test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-N1.bam'
test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
chrom1, pos1, chrom2, pos2 = 'chr20', 29045066, 'chr21', 5216246
margin = 100
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, wobble=100, debug=True)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, wobble=100, debug=False)
sv_treads1 = [read for read in sv_treads1]
# sv_treads2 = [read for read in sv_treads2]
```

```python
sv_twread_qnames1 = {
    '0cfd0362-11aa-41b3-8be6-a21f02b877ff',
    '6dd5f766-8c89-472c-831e-277324ce5d0c',
    'facad054-a0a5-4aeb-aeed-bcd0e12a51db'
}
```

```python
set([read.qname for read in sv_twreads1]) == sv_twread_qnames1
```

#### unit-test (DEL): PMC-01-PDOX1 chr16:6512983-chr16:6645209

```python
test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/test_data/chr16_6512983__chr16_6645209.bam'
test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
region, chrom1, pos1, chrom2, pos2 = ('chr16:6512983-chr16:6645209', 'chr16', 6512983, 'chr16', 6645209)
margin = 100
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=False)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False)
sv_treads1 = [read for read in sv_treads1]
sv_treads2 = [read for read in sv_treads2]
```

```python
sv_tread_qnames1 = {
    '162cdb2e-f7b8-4039-b803-d6fe16811c7a',
    '673cfc41-5867-442f-93ac-18fd0474a155',
    '81ac77f6-a95d-4626-be2c-c023891953d1',
    '8d9521b1-2ca3-403f-ba51-7bbf9f83892d',
    'eebf0d9e-e5d6-4770-8e5e-5e233bdad613'
}
```

```python
sv_tread_qnames2 = {
    '162cdb2e-f7b8-4039-b803-d6fe16811c7a',
    '673cfc41-5867-442f-93ac-18fd0474a155',
    '81ac77f6-a95d-4626-be2c-c023891953d1',
    '8d9521b1-2ca3-403f-ba51-7bbf9f83892d',
    'eebf0d9e-e5d6-4770-8e5e-5e233bdad613'
}
```

```python
set([read.qname for read in sv_treads1]) == sv_tread_qnames1
```

```python
set([read.qname for read in sv_treads2]) == sv_tread_qnames2
```

#### unit-test (DUP): PMC-01-PDOX1 chr3:79339417 chr3:79170916

```python
test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/test_data/chr3_79339417__chr3_79170916.bam'
test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
chrom1, pos1, chrom2, pos2 = 'chr3', 79339417, 'chr3', 79170916
margin = 100
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=False)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False)
sv_treads1 = [read for read in sv_treads1]
sv_treads2 = [read for read in sv_treads2]
```

```python
sv_tread_qnames1 = {
    '0ebf68ee-ceff-49aa-bc68-fcb93dd2735c',
    '159ad899-74a9-41f5-9fce-05f0521600be',
    '3330e509-c5fa-4fc3-ba7c-e74954e8c86d',
    '3d1be97c-fe2d-49a8-b149-77fa733e5141',
    '592deca1-d59c-4bf8-961b-419ce9d834f1',
    '6002887e-89ff-440a-8ca3-d5595b1161ee',
    '6816b1ce-c9a5-460f-988c-46318276a55a',
    '751eb2f4-5be5-4858-b15f-c16cec5b5b69',
    '7aa30c22-d265-4407-8787-9fabaa5c7452',
    'fd60ad61-1f7e-4abf-9342-443b90d52dda'
}
```

```python
set([read.qname for read in sv_treads1]) == sv_tread_qnames1
```

#### unit-test ("D" cigar DEL): PMC-01-PDOX1 chr2:96071256-chr2:96072063

```python
test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/test_data/chr2_96071256__chr2_96072063.bam'
test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
chrom1, pos1, chrom2, pos2 = 'chr2', 96071256, 'chr2', 96072063
margin = 100
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=False)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False)
sv_treads1 = [read for read in sv_treads1]
sv_treads2 = [read for read in sv_treads2]
```

```python
sv_tread_qnames1 = {
    '6f1bb5aa-74ac-462c-9c52-402caf7139dd',
    '84fbc76b-b4b3-4f06-8534-b32bf9f478b2'
}
```

```python
set([read.qname for read in sv_treads1]) == sv_tread_qnames1
```

#### unit-test (TRA): PMC-01-T1 chr2:103362006-chr22:47379382 (SV from T1)

```python
test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/test_data/chr2_103362006__chr22_47379382.bam'
test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
chrom1, pos1, chrom2, pos2 = 'chr2', 103362006, 'chr22', 47379382
margin = 100
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=False)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False)
sv_treads1 = [read for read in sv_treads1]
sv_treads2 = [read for read in sv_treads2]
```

```python
sv_tread_qnames2 = {
    '1323da6b-6eaa-4934-8bdb-cc02c165581c',
    '14f8489d-c129-427b-b526-9210688bfa4b',
    '263049d3-f63e-4eb5-9c5b-1f46a21bd622',
    '421a1b8d-f53e-46a7-a603-723e8843e81b',
    '4a33e454-64ca-4b40-9629-622241d6b3f7',
    '65319648-3e3d-408a-8e6a-58dc9d34a0c7',
    '75a627e7-c290-4ff1-80df-e2b9343fa9b5'
}
```

```python
set([read.qname for read in sv_treads2]) == sv_tread_qnames2
```

#### unit-test (TRA): PMC-01-PDOX1 chr2:201284717-chr19:46833494 (SV from ORG1)

```python
test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/test_data/chr2_201284717__chr19_46833494.bam'
# test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-N1.bam'
test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
region, chrom1, pos1, chrom2, pos2 = ('chr2:201284717-chr19:46833494', 'chr2', 201284717, 'chr19', 46833494)
margin = 100
treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=False)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False)
sv_treads1 = [read for read in sv_treads1]
sv_treads2 = [read for read in sv_treads2]
```

```python
sv_tread_qnames1 = {
    '29eb9055-bdca-45e6-82d1-3937e280ed47',
    '74c4aadc-9ab6-430e-b9f1-5d2beae90e64',
    'df75cf49-a20b-47e4-a934-d5cc4ac53fb9'
}
```

```python
sv_tread_qnames2 = {
    '29eb9055-bdca-45e6-82d1-3937e280ed47',
    '74c4aadc-9ab6-430e-b9f1-5d2beae90e64',
    'df75cf49-a20b-47e4-a934-d5cc4ac53fb9'
}
```

```python
set([read.qname for read in sv_treads1]) == sv_tread_qnames1
```

```python
set([read.qname for read in sv_treads2]) == sv_tread_qnames2
```

### test: PMC-01-PDOX1 chr16:6512983-chr16:6645209

```python
nbam = bams['PMC-01-N1']
likely_tp = svs['PMC-01-T1'] | svs['PMC-01-ORG1']
likely_fp = ((svs['PMC-01-PDOX1'] - svs['PMC-01-PDOX2']) - likely_tp) | ((svs['PMC-01-PDOX2'] - svs['PMC-01-PDOX1']) - likely_tp)

cols = ['NM', 'ms', 'AS', 'nn', 'cm', 's1', 'de', 'mapq', 'match', 'ins',
        'del', 'skip', 'softclip', 'hardclip', 'pad', 'qual', 'diff', 'back']
cols = [s+'_t' for s in cols] + [s+'_n' for s in cols] + ['TP']
svstats = pd.DataFrame(columns=cols)
for sample in samples:
    ssvs = svs[sample] # sample svs
    tbam = bams[sample]

    for region in list(ssvs):
        # if region != 'chr16:6512983-chr16:6645209': continue # unittest1
        if region != 'chr3:79339417-chr3:79170916': continue
        
        flt_tag = 0
        flt_tag += (region in valid_tp)
        flt_tag -= (region in likely_fp)

        reg1, reg2 = region.split('-')
        chrom1, pos1 = reg1.split(':')
        chrom2, pos2 = reg2.split(':')
        pos1, pos2 = int(pos1), int(pos2)

        nreads1 = [read for read in nbam.fetch(chrom1, pos1-1, pos1)]
        nreads2 = [read for read in nbam.fetch(chrom2, pos2-1, pos2)]
        treads1 = [read for read in tbam.fetch(chrom1, pos1-1, pos1)]
        treads2 = [read for read in tbam.fetch(chrom2, pos2-1, pos2)]
        
        sv_nreads1 = get_sv_supporting_reads(nreads1, chrom1, pos1, chrom2, pos2)
        sv_nreads2 = get_sv_supporting_reads(nreads2, chrom2, pos2, chrom1, pos1)
        sv_treads1 = get_sv_supporting_reads(treads1, chrom1, pos1, chrom2, pos2)
        sv_treads2 = get_sv_supporting_reads(treads2, chrom2, pos2, chrom1, pos1)
        raise ValueError
```

```python
tbam = bams['PMC-01-PDOX2']
```

```python
region, chrom1, pos1, chrom2, pos2 = 'chr3:79339417-chr3:79170916', 'chr3', 79339417, 'chr3', 79170916
```

```python
margin = 100
```

```python
nreads1 = [read for read in nbam.fetch(chrom1, pos1-margin, pos1+margin)]
nreads2 = [read for read in nbam.fetch(chrom2, pos2-margin, pos2+margin)]
treads1 = [read for read in tbam.fetch(chrom1, pos1-margin, pos1+margin)]
treads2 = [read for read in tbam.fetch(chrom2, pos2-margin, pos2+margin)]

sv_nreads1 = get_sv_supporting_reads(nreads1,  chrom1, pos1, chrom2, pos2, debug=True)
sv_nreads2 = get_sv_supporting_reads(nreads2,  chrom2, pos2, chrom1, pos1, debug=True)
sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=True)
sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=True)
```

```python
for read in treads1:
    tags = {k:v for (k,v) in read.tags}
    if 'SA' in tags:
        print(read.qname, '-' if read.is_reverse else '+', tags['SA'].split(','))
```

```python
for read in treads2:
    tags = {k:v for (k,v) in read.tags}
    if 'SA' in tags:
        print(read.qname, tags['SA'])
```

```python
sv_nreads1 = [read for read in sv_nreads1]
```

```python
sv_treads1 = [read for read in sv_treads1]
```

```python
%debug
```

```python
sv_treads2 = [read for read in sv_treads2]
```

```python
sv_treads1
```

```python
treads1
```

```python
sv_treads2
```

```python
read = nreads1[0]
```

```python
sv_reads = list(sv_nreads1)
```

```python
sv_read = sv_reads[0]
```

```python
sv_read
```

```python
chrom1, pos1
```

## extract SV-supporting read stats

```python
# test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/test_data/chr2_201284717__chr19_46833494.bam'
# # test_bam_path = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-N1.bam'
# test_bam = pysam.AlignmentFile(test_bam_path, 'rb')
# region, chrom1, pos1, chrom2, pos2 = ('chr2:201284717-chr19:46833494', 'chr2', 201284717, 'chr19', 46833494)
# margin = 100
# treads1 = [read for read in test_bam.fetch(chrom1, pos1-margin, pos1+margin)]
# treads2 = [read for read in test_bam.fetch(chrom2, pos2-margin, pos2+margin)]
# sv_treads1 = get_sv_supporting_reads(treads1,  chrom1, pos1, chrom2, pos2, debug=False)
# sv_treads2 = get_sv_supporting_reads(treads2,  chrom2, pos2, chrom1, pos1, debug=False)
# sv_treads1 = [read for read in sv_treads1]
# sv_treads2 = [read for read in sv_treads2]
```

```python
# margin = 10
# fsvs = defaultdict(set)
# nbam = bams['PMC-01-N1']
# for sample in samples:
#     nsvs = 0
#     small_sv, low_ndepth, low_tdepth, nmapq_low, tmapq_low = 0, 0, 0, 0, 0

#     # if sample != 'PMC-01-PDOX1': continue
#     psvs = svs[sample]
#     tbam = bams[sample]
#     for region in tqdm(list(psvs)):
#         reg1, reg2 = region.split('-')
#         chrom1, pos1 = reg1.split(':')
#         chrom2, pos2 = reg2.split(':')
#         pos1, pos2 = int(pos1), int(pos2)
            
#         nreads = nbam.fetch(chrom1, pos1-margin, pos1+margin)
#         nmapqs = [read.mapq for read in nreads]
#         treads = tbam.fetch(chrom1, pos1-margin, pos1+margin)
#         treads = [r for r in get_sv_supporting_reads(treads, chrom1, pos1, chrom2, pos2)]
#         tmapqs = [read.mapq for read in treads]
        
#         if len(nmapqs) < 10: # read depth filter
#             low_ndepth += 1
#             continue
        
#         if len(tmapqs) < 10:
#             low_tdepth += 1
#             continue
            
#         if nmapqs.count(0) / len(nmapqs) > 0.5: 
#             nmapq_low += 1
#             continue
        
#         if tmapqs.count(0) / len(nmapqs) > 0.5: 
#             tmapq_low += 1
#             continue
        
#         nsvs += 1
        
#         fsvs[sample].add(region)
#     print(sample, small_sv, low_ndepth, low_tdepth, nmapq_low, tmapq_low, nsvs)
```

## umap for sv stats

```python
import umap
import seaborn as sns
from sklearn.preprocessing import StandardScaler
```

```python
_backup = svstats.copy()
```

```python
f_cols = ['TP'] + svstats.columns[:-1].tolist()
flt_cols = []
for col in f_cols:
    n_vals = svstats[col].unique().shape[0]
    if n_vals == 1:
        continue
    flt_cols.append(col)
f_cols = flt_cols
```

```python
data = svstats.loc[:, f_cols].values.astype(float)
scaler = StandardScaler()
reducer = umap.UMAP(random_state=42)
scaled = scaler.fit_transform(data)
embedding = reducer.fit_transform(scaled)
data = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])
data = data.join(svstats[f_cols].astype(float))
```

```python
log10cols = ['NM_t', 'NM_n']
for col in log10cols:
    data[col] = np.log10(data[col])
```

```python
data['match_t'].clip(inplace=True, upper=210)
data['match_n'].clip(inplace=True, upper=210)
data['ins_t'].clip(inplace=True, upper=110)
data['ins_n'].clip(inplace=True, upper=110)
data['del_t'].clip(inplace=True, upper=150)
data['del_n'].clip(inplace=True, upper=150)
data['softclip_t'].clip(inplace=True, lower=.9)
data['softclip_n'].clip(inplace=True, lower=.9)
data['hardclip_t'].clip(inplace=True, upper=1)
data['hardclip_n'].clip(inplace=True, upper=1)
```

```python
nrow=4; ncol=7
fig, axes = plt.subplots(nrow, ncol)
fig.set_figheight(nrow*2)
fig.set_figwidth(ncol*2)
for i in range(nrow*ncol):
    x = i % ncol
    y = i // ncol
    ax = axes[y][x]
    if i < len(f_cols): 
        feature = f_cols[i]
    else:
        ax.axis('off')
        continue
    color = 'tab:red' if feature.endswith('_t') else 'tab:blue'
    if feature == 'TP': color='tab:green'
    data[feature].hist(ax=ax, color=color)
    ax.set_title(feature)
plt.tight_layout()
```

```python
nrow=5; ncol=6
fig, axes = plt.subplots(nrow, ncol)
fig.set_figheight(nrow*4)
fig.set_figwidth(ncol*4)

suptitle = 'SV features'
for i in range(nrow*ncol):
    x = i % ncol
    y = i // ncol
    # print(f'x, y = {x}, {y}')
    ax = axes[y][x]
    if i < len(f_cols): 
        feature = f_cols[i]
    else:
        ax.axis('off')
        continue
    print(f'plotting umap for {feature}')
    palette = 'viridis' if feature != 'TP' else 'tab20_r'
    sns.scatterplot(ax=ax, data=data, x='UMAP1', y='UMAP2', 
                    alpha=0.7, s=10, palette=palette, hue=feature)
    # break
fig.suptitle(suptitle)
plt.tight_layout()
```

# Genotype SVs across sample


## set data


### make (ordered) union sv set

```python
chroms = ['chr'+str(c) for c in range(1, 22+1)] + ['chrX', 'chrY']
union_svs = set()
for sample in samples:
    saved = set()
    sample_svs = set()
    _svs = svs[sample]
    for sv in _svs:
        reg1, reg2 = sv.split('-')
        chrom1, pos1 = reg1.split(':')
        chrom2, pos2 = reg2.split(':')
        pos1, pos2 = int(pos1), int(pos2)
        if chrom1 == chrom2:
            if pos1 > pos2:
                chrom1, pos1, chrom2, pos2 = chrom2, pos2, chrom1, pos1
            if abs(pos2 - pos1) < 1000: continue # likely FP
            # continue
        else:
            if chroms.index(chrom1) > chroms.index(chrom2):
                chrom1, pos1, chrom2, pos2 = chrom2, pos2, chrom1, pos1
        sv = f'{chrom1}:{pos1}-{chrom2}:{pos2}'
        svid = (chrom1, pos1, chrom2, pos2)
        if svid in saved: continue
        saved.add(svid)
        sample_svs.add(sv)

    union_svs.update(sample_svs)
    print(sample, len(sample_svs))
union_svs = list(union_svs)
```

### make sv evidence table (~5m)

```python
nbam = pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-N1.bam', 'rb')
```

```python
sv_cols = ['svid', 'chrom1', 'pos1', 'chrom2', 'pos2', 'N1_depth', 'N1_support']
for sample in samples:
    sample = sample.replace('PMC-01-', '')
    sv_cols += [f'{sample}_depth', f'{sample}_support']
svdf = pd.DataFrame(columns=sv_cols)

margin = 10 # for fetching reads
for ix, sv in tqdm(enumerate(union_svs), total=len(union_svs)): # 'chr6:29925879-chr6:29925917'
    reg1, reg2 = sv.split('-')
    chrom1, pos1 = reg1.split(':')
    chrom2, pos2 = reg2.split(':')
    pos1, pos2 = int(pos1), int(pos2)
    
    nreads = [read for read in nbam.fetch(chrom1, pos1-margin, pos1+margin)]
    sv_nreads = [r for r in get_sv_supporting_reads(nreads, chrom1, pos1, chrom2, pos2, wobble=50)]
    field = [sv, chrom1, pos1, chrom2, pos2,  len(nreads), len(sv_nreads)]
    for sample in samples:
        tbam = bams[sample]
        treads = [read for read in tbam.fetch(chrom1, pos1-margin, pos1+margin)]
        sv_treads = [r for r in get_sv_supporting_reads(treads, chrom1, pos1, chrom2, pos2, wobble=30)]
        
        field += [len(treads), len(sv_treads)]
    svdf.loc[svdf.shape[0]] = field
    # if ix == 10: break
```

```python
svdf.shape
```

```python
out_path = '/juno/work/shah/users/chois7/tickets/drostpmc/tables/sv_support.tsv'
if os.path.exists(out_path):
    svdf = pd.read_table(out_path)
else:
    svdf.to_csv(out_path, sep='\t', index=False)
```

### label somatic or germline

```python
new_cols = ['svid', 'chrom1', 'pos1', 'chrom2', 'pos2', 
            'T1_depth', 'T1_support', 'ORG1_depth', 'ORG1_support', 
            'PDOX1_depth', 'PDOX1_support', 'PDOX2_depth', 'PDOX2_support', 
            'N1_depth', 'N1_support',]
svdf = svdf[new_cols]
```

### make somatic sv set

```python
somatic_sv_min_alt_cutoff = 1
germline_sv_max_alt_cutoff = 0
normal_short_sample = 'N1'

svset = {s: set() for s in samples}
for sample in samples:
    short_sample = sample.replace('PMC-01-', '')
    _svdf = svdf[
        (svdf[f'{short_sample}_support'] >= somatic_sv_min_alt_cutoff) &
        (svdf[f'{normal_short_sample}_support'] <= germline_sv_max_alt_cutoff)
    ]
    svset[sample] = set(_svdf['svid'])
```

```python
union = svset['PMC-01-T1'] | svset['PMC-01-ORG1'] | svset['PMC-01-PDOX1'] | svset['PMC-01-PDOX2']
```

```python
union - (svset['PMC-01-T1'] | svset['PMC-01-PDOX1'] | svset['PMC-01-PDOX2'])
```

```python
svset['PMC-01-T1'] & svset['PMC-01-ORG1'] & svset['PMC-01-PDOX1'] & svset['PMC-01-PDOX2']
```

## plot venn diagram for sv set

```python
import matplotlib.font_manager as font_manager
arial_path = "/home/chois7/.local/share/fonts/Arial.ttf"
# arial_italic_path = "/path/to/Arial-Italic.ttf"
arial = font_manager.FontProperties(fname=arial_path)
# arial_italic = font_manager.FontProperties(fname=arial_italic_path)
plt.rcParams['font.family'] = arial.get_name()
```

```python
fig, ax = plt.subplots(figsize=(6,6))

group1 = 'PMC-01-T1'
group2 = 'PMC-01-ORG1'
group3 = 'PMC-01-PDOX1'
group4 = 'PMC-01-PDOX2'
venny4py(sets={
    group1 + ' (tumor)': svset[group1], 
    group2 + ' (organoid)': svset[group2], 
    group3 + ' (PDX1)': svset[group3], 
    group4 + ' (PDX2)': svset[group4],
}, size=5, asax=ax,)
```

## create IGV targets based on SV comparison

```python
import os
```

```python
def make_igv_region_and_cmd(igvset, tag, igv_dir):
    # tag = 'Shared_in_PDOX'
    # igv_dir = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/svset_comparison'
    bam_files = 'PMC-01-T1.bam PMC-01-ORG1.bam PMC-01-PDOX1.bam PMC-01-PDOX2.bam PMC-01-N1.bam'
    margin = 500
    saved_coords = set()
    with open(f'{igv_dir}/region_{tag}.txt', 'w') as out:
        for coord in list(igvset):
            reg1, reg2 = coord.split('-')
            chrom1, pos1 = reg1.split(':')
            chrom2, pos2 = reg2.split(':')
            pos1, pos2 = int(pos1), int(pos2)
            line = f'{chrom1}:{pos1-margin}-{pos1+margin} {chrom2}:{pos2-margin}-{pos2+margin} {tag}'
            if (chrom1, pos1) in saved_coords or (chrom2, pos2) in saved_coords:
                continue
            else:
                saved_coords.add((chrom1, pos1))
                saved_coords.add((chrom2, pos2))
                out.write(line + '\n')
    with open(f'{igv_dir}/cmd_{tag}.sh', 'w') as cmdfile:
        cmd = "singularity run -B /juno -B /home docker://shahcompbio/igver igver.py "
        cmd += f'--bam {bam_files} -r region_{tag}.txt --config view_breakpoint.batch -o png_{tag} -g hg38'
        cmdfile.write(cmd + '\n')
```

```python
uniqsets = {}
igv_dir = '/juno/work/shah/users/chois7/tickets/drostpmc/igv/svset_comparison'
for ix in range(1, 15+1):
    targets = ['PMC-01-T1', 'PMC-01-ORG1', 'PMC-01-PDOX1', 'PMC-01-PDOX2']
    som = svset['PMC-01-T1'] | svset['PMC-01-ORG1'] | svset['PMC-01-PDOX1'] | svset['PMC-01-PDOX2'] # somatic svs
    bit = f'{ix:04b}'
    i0, i1, i2, i3 = bit
    i0, i1, i2, i3 = int(i0), int(i1), int(i2), int(i3)
    if i0 == 0:
        targets.remove('PMC-01-T1')
        som -= svset['PMC-01-T1']
    else:
        som &= svset['PMC-01-T1']
    if i1 == 0:
        targets.remove('PMC-01-ORG1')
        som -= svset['PMC-01-ORG1']
    else:
        som &= svset['PMC-01-ORG1']
    if i2 == 0:
        targets.remove('PMC-01-PDOX1')
        som -= svset['PMC-01-PDOX1']
    else:
        som &= svset['PMC-01-PDOX1']
    if i3 == 0:
        targets.remove('PMC-01-PDOX2')
        som -= svset['PMC-01-PDOX2']
    else:
        som &= svset['PMC-01-PDOX2']
    target_str = ' & '.join(targets)
    target_tag = '__'.join(targets)
    uniqsets[target_str] = som
    make_igv_region_and_cmd(som, target_tag, igv_dir)
    os.system(f'mkdir -p {igv_dir}/png_{target_tag}')
    print(target_str, len(som))
```

```python
svdf[svdf.svid.str.count('66039803')>0]
```

### plot tree

```python
import seaborn as sns
```

```python
ordered_samples = ['PMC-01-T1', 'PMC-01-ORG1', 'PMC-01-PDOX1', 'PMC-01-PDOX2']
```

```python
union = svset['PMC-01-T1'] | svset['PMC-01-ORG1'] | svset['PMC-01-PDOX1'] | svset['PMC-01-PDOX2'] # somatic svs
```

```python
phylodf = pd.DataFrame(columns=ordered_samples, index=union)
phylodf.loc[:, :] = 0
union_set = sv
for sample in ordered_samples:
    _svs = list(svset[sample])
    phylodf.loc[_svs, sample] = 1
    phylodf[sample] = phylodf[sample].astype(int)
```

```python
phylodf
```

```python
# fig, ax = plt.subplots(figsize=(4, 20))
g = sns.clustermap(phylodf, row_cluster=True, col_cluster=True, figsize=(5,7), cbar_pos=None, cmap="Reds")
ax = g.ax_heatmap
ax.set_yticks([])
```

# Venn diagrams

```python
samples
```

```python
group1 = 'PMC-01-T1'
group2 = 'PMC-01-ORG1'
venn2(subsets=[svset[group1], svset[group2]], set_labels=[group1, group2])
```

```python
group1 = 'PMC-01-PDOX1'
group2 = 'PMC-01-PDOX2'
venn2(subsets=[svset[group1], svset[group2]], set_labels=[group1, group2])
```

```python
group1 = 'PMC-01-T1'
group2 = 'PMC-01-ORG1'
venn2(subsets=[fsvs[group1], fsvs[group2]], set_labels=[group1, group2])
```

```python
fig, ax = plt.subplots(figsize=(4,4))
group1 = 'PMC-01-T1'
group2 = 'PMC-01-ORG1'
venn2(subsets=[svs[group1], svs[group2]], set_labels=[group1, group2], ax=ax)
```

```python
group1 = 'PMC-01-PDOX1'
group2 = 'PMC-01-PDOX2'
venn2(subsets=[svs[group1], svs[group2]], set_labels=[group1, group2])
```

```python
group1 = 'PMC-01-PDOX1'
group2 = 'PMC-01-PDOX2'
venn2(subsets=[svs[group1], svs[group2]], set_labels=[group1, group2])
```

```python
group1 = 'PMC-01-PDOX1'
group2 = 'PMC-01-PDOX2'
venn2(subsets=[fsvs[group1], fsvs[group2]], set_labels=[group1, group2])
```

```python
?venny4py
```

```python
svs[group1]
```

```python
fig, ax = plt.subplots(figsize=(6,6))

group1 = 'PMC-01-T1'
group2 = 'PMC-01-ORG1'
group3 = 'PMC-01-PDOX1'
group4 = 'PMC-01-PDOX2'
venny4py(sets={
    group1 + ' (tumor)': svs[group1], 
    group2 + ' (organoid)': svs[group2], 
    group3 + ' (PDX1)': svs[group3], 
    group4 + ' (PDX2)': svs[group4],
}, size=5, asax=ax,)
```

```python
fig, ax = plt.subplots(figsize=(6,6))

group1 = 'PMC-01-T1'
group2 = 'PMC-01-ORG1'
group3 = 'PMC-01-PDOX1'
group4 = 'PMC-01-PDOX2'
venny4py(sets={
    group1 + ' (tumor)': fsvs[group1], 
    group2 + ' (organoid)': fsvs[group2], 
    group3 + ' (PDX1)': fsvs[group3], 
    group4 + ' (PDX2)': fsvs[group4],
}, size=5, asax=ax,)
```

```python
match.replace({'NaN':0, 'NA':1}).value_counts()
```

```python
val.split(':')
```

```python
pattern = '[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:([^:]+):[^:]+:[^:]+:[^:]+:[^:]+:([^:]+)'
re.search(pattern, val).groups()
```

# Extract coordinates for IGV

```python
results_dir = '/juno/work/shah/users/chois7/tickets/savana/pipeline/results'
```

```python
samples
```

```python
samples_with_normals = samples + ['PMC-01-N1']
```

```python
sample_svs = {}
```

```python
shared_in_pdox = svs['PMC-01-PDOX1'] & svs['PMC-01-PDOX2']
union_in_human = svs['PMC-01-T1'] & svs['PMC-01-ORG1']
```

```python
igvset = union_in_human #
tag = 'Union_in_Human'

igv_dir = '/juno/work/shah/users/chois7/tickets/drostpmc'
bam_files = 'PMC-01-T1.bam PMC-01-ORG1.bam PMC-01-PDOX1.bam PMC-01-PDOX2.bam PMC-01-N1.bam'
margin = 1000
saved_coords = set()
with open(f'{igv_dir}/igv/region_{tag.lower()}.txt', 'w') as out:
    for coord in list(igvset):
        reg1, reg2 = coord.split('-')
        chrom1, pos1 = reg1.split(':')
        chrom2, pos2 = reg2.split(':')
        pos1, pos2 = int(pos1), int(pos2)
        line = f'{chrom1}:{pos1-margin}-{pos1+margin} {chrom2}:{pos2-margin}-{pos2+margin} {tag}'
        if (chrom1, pos1) in saved_coords or (chrom2, pos2) in saved_coords:
            continue
        else:
            saved_coords.add((chrom1, pos1))
            saved_coords.add((chrom2, pos2))
            out.write(line + '\n')
with open(f'{igv_dir}/igv/cmd_{tag.lower()}.sh', 'w') as cmdfile:
    cmd = "singularity run -B /juno -B /home docker://shahcompbio/igver igver.py "
    cmd += f'--bam {bam_files} -r region_{tag.lower()}.txt --config view_breakpoint.batch -o png -g hg38'
    cmdfile.write(cmd + '\n')
```

```python
igvset = shared_in_pdox #
tag = 'Shared_in_PDOX'

igv_dir = '/juno/work/shah/users/chois7/tickets/drostpmc'
bam_files = 'PMC-01-T1.bam PMC-01-ORG1.bam PMC-01-PDOX1.bam PMC-01-PDOX2.bam PMC-01-N1.bam'
margin = 1000
saved_coords = set()
with open(f'{igv_dir}/igv/region_{tag.lower()}.txt', 'w') as out:
    for coord in list(igvset):
        reg1, reg2 = coord.split('-')
        chrom1, pos1 = reg1.split(':')
        chrom2, pos2 = reg2.split(':')
        pos1, pos2 = int(pos1), int(pos2)
        line = f'{chrom1}:{pos1-margin}-{pos1+margin} {chrom2}:{pos2-margin}-{pos2+margin} {tag}'
        if (chrom1, pos1) in saved_coords or (chrom2, pos2) in saved_coords:
            continue
        else:
            saved_coords.add((chrom1, pos1))
            saved_coords.add((chrom2, pos2))
            out.write(line + '\n')
with open(f'{igv_dir}/igv/cmd_{tag.lower()}.sh', 'w') as cmdfile:
    cmd = "singularity run -B /juno -B /home docker://shahcompbio/igver igver.py "
    cmd += f'--bam {bam_files} -r region_{tag.lower()}.txt --config view_breakpoint.batch -o png -g hg38'
    cmdfile.write(cmd + '\n')
```

```python
igvset = likely_tp #
tag = 'T1_or_ORG1'

igv_dir = '/juno/work/shah/users/chois7/tickets/drostpmc'
bam_files = 'PMC-01-T1.bam PMC-01-ORG1.bam PMC-01-PDOX1.bam PMC-01-PDOX2.bam PMC-01-N1.bam'
margin = 1000
saved_coords = set()
with open(f'{igv_dir}/igv/region_{tag.lower()}.txt', 'w') as out:
    for coord in list(igvset):
        reg1, reg2 = coord.split('-')
        chrom1, pos1 = reg1.split(':')
        chrom2, pos2 = reg2.split(':')
        pos1, pos2 = int(pos1), int(pos2)
        line = f'{chrom1}:{pos1-margin}-{pos1+margin} {chrom2}:{pos2-margin}-{pos2+margin} {tag}'
        if (chrom1, pos1) in saved_coords or (chrom2, pos2) in saved_coords:
            continue
        else:
            saved_coords.add((chrom1, pos1))
            saved_coords.add((chrom2, pos2))
            out.write(line + '\n')
with open(f'{igv_dir}/igv/cmd_{tag.lower()}.sh', 'w') as cmdfile:
    cmd = "singularity run -B /juno -B /home docker://shahcompbio/igver igver.py "
    cmd += f'--bam {bam_files} -r region_{tag.lower()}.txt --config view_breakpoint.batch -o png -g hg38'
    cmdfile.write(cmd + '\n')
```

# Feature extraction of mouse reads

```python
import pysam
```

```python
n1_bam = pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-N1.bam', 'rb')
t1_bam = pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-T1.bam', 'rb')
org1_bam = pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-ORG1.bam', 'rb')
pdox1_bam = pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-PDOX1.bam', 'rb')
pdox2_bam = pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-PDOX2.bam', 'rb')
```

```python
def make_tag_df(bam, chrom, pos):
    margin = 500
    df_cols = ['NM', 'ms', 'AS', 'nn', 'tp', 
               'cm', 's1', 's2', 'de', 'rl', 
               'mapq', 'alen']
    tag_df = pd.DataFrame(columns=df_cols)
    qnames = []
    for i, read in enumerate(bam.fetch(chrom, pos-margin, pos+margin, until_eof=True)):
        nm_tag = read.tags[0]
        assert nm_tag[0] == 'NM', read.tags
        field = [0] * len(df_cols)
        for tag_name, tag_stat in read.tags:
            if tag_name in set(df_cols):
                ix = df_cols.index(tag_name)
            else:
                continue
            field[ix] = tag_stat
        field[ix+1] = read.mapq
        field[ix+2] = read.alen
        qnames.append(read.qname)
        tag_df.loc[tag_df.shape[0]] = field
    tag_df.index = qnames
    return tag_df
```

```python
# chrom, pos = 'chr7', 142653049 
chrom, pos = 'chr6', 32639731
# chrom, pos = 'chr1', 223558935 # TP

n1_df = make_tag_df(n1_bam, chrom, pos)
t1_df = make_tag_df(t1_bam, chrom, pos)
org1_df = make_tag_df(t1_bam, chrom, pos)
pdox1_df = make_tag_df(pdox1_bam, chrom, pos)
pdox2_df = make_tag_df(pdox2_bam, chrom, pos)
```

```python
n1_df
```

```python
pdox1_df
```

```python
bins = np.arange(0, 0.3, 0.01)
alpha = 0.3
fig, ax = plt.subplots(figsize=(4,3))
n1_df['de'].hist(ax=ax, bins=bins, alpha=alpha, label='normal')
t1_df['de'].hist(ax=ax, bins=bins, alpha=alpha, label='tumor')
org1_df['de'].hist(ax=ax, bins=bins, alpha=alpha * 2, label='organoid')
pdox1_df['de'].hist(ax=ax, bins=bins, alpha=alpha * 2, label='PDOX1', color='black')
pdox2_df['de'].hist(ax=ax, bins=bins, alpha=alpha * 2, label='PDOX2', color='brown')
ax.set_xlabel('sequence divergence'); ax.set_ylabel('density'); 
ax.set_title(f'Alignments in {chrom}:{pos-margin}-{pos+margin}');
ax.legend()
```

```python
n1_df.sort_values(['de'])
```

```python
t1_df.sort_values(['de'])
```

```python
pdox1_df.sort_values(['de'])
```

```python
pdox2_df.sort_values(['de'])
```

```python
pd.Series(nm_fracs).hist()
```

```python
pd.Series(nm_fracs).hist()
```

# Call SNVs nearby SV positions


## create SV position intervals

```python
svs.keys()
```

```python
chroms = ['chr'+str(c) for c in range(1, 22+1)] + ['chrX', 'chrY']
margin = 200
saved_breakpoints = set()
for sample in samples:
    svset = svs[sample]
    for regions in list(svset):
        reg1, reg2 = regions.split('-')
        chrom1, pos1 = reg1.split(':')
        chrom2, pos2 = reg2.split(':')
        pos1, pos2 = int(pos1), int(pos2)
        preg1 = f'{chrom1}:{pos1-margin}-{pos1+margin}'
        preg2 = f'{chrom2}:{pos2-margin}-{pos2+margin}'
        if preg1 not in saved_breakpoints:
            saved_breakpoints.add(preg1)
        if preg2 not in saved_breakpoints:
            saved_breakpoints.add(preg2)
```

```python
chroms.index(list(saved_breakpoints)[0].split(':')[0])
```

```python
saved_breakpoints = sorted(
    list(saved_breakpoints), key=lambda x: 
    (chroms.index(x.split(':')[0]),
    int(x.split(':')[1].split('-')[0]))
)
```

```python
with open('/juno/work/shah/users/chois7/tickets/drostpmc/snv/perisvs.intervals', 'w') as out:
    for preg in list(saved_breakpoints):
        out.write(preg + '\n')
```

# Pileup scan snv

```python
n1_bam = pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-N1.bam', 'rb')
pdox1_bam = pysam.AlignmentFile('/juno/work/shah/users/chois7/tickets/drostpmc/igv/PMC-01-PDOX1.bam', 'rb')

chrom = 'chr6'
pos = 32639731
margin = 500
for i, pileup in enumerate(pdox1_bam.pileup(chrom, pos-margin, pos+margin)):
    if pileup.pos == pos: break
    pass
    # nm_tag = read.tags[0]
    # assert nm_tag[0] == 'NM', read.tags
    # field = [0] * len(df_cols)
    # for tag_name, tag_stat in read.tags:
    #     if tag_name in set(df_cols):
    #         ix = df_cols.index(tag_name)
    #     else:
    #         continue
    #     field[ix] = tag_stat
    # tag_df.loc[tag_df.shape[0]] = field
```

```python
_ = pileup.get_num_aligned
```

```python
_.
```

```python
[s for s in dir(pileup) if not s.startswith('__')]
```

```python
# pileup.get_query_sequences(mark_matches=False) # don't!
```

# Phylo


## test

```python
import seaborn as sns
```

```python
# fig, ax = plt.subplots(figsize=(4, 20))
sns.clustermap(df, row_cluster=True, col_cluster=True, figsize=(4, 10))
```

```python
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

# Create a synthetic distance matrix (replace this with your actual data)
data = np.random.rand(4, 10)
data = df.copy()

# Create a dendrogram for hierarchical clustering
# dendrogram = sch.dendrogram(sch.linkage(data, method='ward'))

# Reorder the data matrix based on the hierarchical clustering
order = dendrogram['leaves']
data_ordered = data.iloc[order, :]#[:, order]

# Create the heatmap
plt.figure(figsize=(8, 6))
plt.imshow(data_ordered, cmap='viridis', aspect='auto', origin='lower')
plt.colorbar(label='Distance')
plt.title('Phylogenetic Tree Heatmap')
plt.xlabel('Species')
plt.ylabel('Species')

# Show the plot
plt.show()
```
