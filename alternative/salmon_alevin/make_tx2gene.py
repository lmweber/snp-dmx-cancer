'''
One of the requirements of every alevin run is an indexed reference
transcriptome (eg gencode), plus a mapping between transcript IDs
and gene IDs. This script makes that mapping file. It hypothetically
only needs a GTF file, but if there are fasta sequences not found in
the GTF file, alevin will throw an error and quit. So to work around
this, we'll process both GTF and fasta files and return a union of
transcript/gene mappings from both.

This script should only need to be run once for any new
reference transcriptome version.
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
gtf_transcripts = gtf_transcripts[['transcript_id', 'gene_id']]

print(gtf_transcripts.head())

fasta_transcripts = []
with gzip.open(args.fasta, 'r') as f:
    for line in f:
        line = line.decode('utf-8')
        if line[0] != ">":
            continue
        attributes = {}
        line = re.sub(">", "", line)
        lineparts = line.split("|")
        attributes['transcript_id'] = lineparts[0]
        attributes['gene_id'] = lineparts[1]
        fasta_transcripts.append(attributes)

fasta_transcripts = pd.DataFrame(fasta_transcripts)

all_transcripts = pd.merge(gtf_transcripts, fasta_transcripts, how='outer')
all_transcripts.to_csv("tx2gene.tsv", sep="\t", index=False, header=False)
