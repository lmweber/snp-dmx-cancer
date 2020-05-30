###########################
# Snakefile to run workflow
###########################

# directories

dir_scripts = "scripts"
dir_data = "../data"
dir_outputs = "../outputs"
dir_timestamps = "../timestamps"


# variables

sample_ids = ["16030X2", "16030X3", "16030X4"]

dirs_fastq = {"16030X2": dir_data + "/16030R/Fastq/16030X2_HJVMLDMXX", 
              "16030X3": dir_data + "/16030R/Fastq/16030X3_HJTWLDMXX", 
              "16030X4": dir_data + "/16030R/Fastq/16030X4_HJTWLDMXX"}

dir_ref = dir_data + "/GRCh38/refdata-cellranger-GRCh38-3.0.0"


# --------------
# run on cluster
# --------------

# command to run pipeline on cluster

# snakemake --cluster "qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=200G" -j 3 --local-cores 20


# -----
# rules
# -----

# default rule (run all)

rule all:
  input:
    expand(dir_outputs + "/{sample}/outs/possorted_genome_bam.bam", sample = sample_ids)


# run cellranger count

rule run_cellranger:
  input:
    script_cellranger = dir_scripts + "/run_cellranger.sh"
  output:
    dir_outputs + "/{sample}/outs/possorted_genome_bam.bam"
  params:
    dir_fastq = lambda wildcards: dirs_fastq[wildcards.sample]
  shell:
    "bash {input.script_cellranger} {wildcards.sample} {params.dir_fastq} {dir_ref} {dir_outputs}"



