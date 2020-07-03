# Bash script to run salmon index


# runtime: 5 min

#qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=200G run_salmon_index.sh


# (1) download files (if required)

wget -nv -P ../data/index ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.pc_transcripts.fa.gz
wget -nv -P ../data/index ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz


# (2) run Python scripts (from Ariel Hippen Anderson)

# note: assumes valid python3 installation in path, with some required packages e.g. pandas

# make gene by transcript mapping file
python3 make_tx2gene.py -d ../data/index -f gencode.v34.pc_transcripts.fa.gz -g gencode.v34.annotation.gtf.gz

# run make_gene2symbol.py
python3 make_gene2symbol.py -d ../data/index -f gencode.v34.pc_transcripts.fa.gz -g gencode.v34.annotation.gtf.gz


# (3) run salmon index

# load required modules on JHPCE server
module load gcc/5.5.0
module load cmake/3.15.4

salmon index -i ../data/index -k 23 --gencode -p 10 -t ../data/index/gencode.v34.pc_transcripts.fa.gz


