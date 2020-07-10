'''
Alevin returns a gene-by-cell count matrix, with the genes as ensembl IDs.
These are unwieldy for a hands-on analysis, so we want to swap them out for
more standard gene symbols. The gold-standard way to convert between IDs is
the R package biomaRt, but in my experience it doesn't have a mapping for
all ids in our matrix, so instead we're going back to the source of the gene
list--the gtf and fasta files.

This script should only need to be run once for any
new reference transcriptome version.
'''

import os
import re
import gzip
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', help='location of index files')
parser.add_argument('-f', '--fasta', help='fa.gz file of ref transcriptome')
parser.add_argument('-g', '--gtf', help='gtf.gz file of ref transcriptome')
args = parser.parse_args()

os.chdir(args.directory)

gtf_transcripts = []
with gzip.open(args.gtf, 'r') as f:
    for line in f:
        line = line.decode('utf-8')
        if line[0] == "#":
            continue
        line = line.strip()
        lineparts = line.split('\t')
        if lineparts[2] != "transcript":
            continue
        attributes = {}
        for keyvaluepair in lineparts[8].split(';'):
            keyplusvalue = keyvaluepair.strip().split()
            if len(keyplusvalue) == 2:
                key = keyplusvalue[0]
                value = re.sub('"', '', keyplusvalue[1])
                attributes[key] = value
        gtf_transcripts.append(attributes)

gtf_transcripts = pd.DataFrame(gtf_transcripts)
gtf_transcripts = gtf_transcripts[['gene_id', 'gene_name']]

fasta_transcripts = []
with gzip.open(args.fasta, 'r') as f:
    for line in f:
        line = line.decode('utf-8')
        if line[0] != ">":
            continue
        attributes = {}
        line = re.sub(">", "", line)
        lineparts = line.split("|")
        attributes['gene_name'] = lineparts[5]
        attributes['gene_id'] = lineparts[1]
        fasta_transcripts.append(attributes)

fasta_transcripts = pd.DataFrame(fasta_transcripts)

all_transcripts = (
    pd.merge(gtf_transcripts, fasta_transcripts, how='outer')
      .drop_duplicates())
all_transcripts.to_csv("gene2symbol.tsv", sep="\t", index=False, header=False)
