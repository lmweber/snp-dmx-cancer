###########################
# Snakefile to run workflow
###########################

# directories

dir_scripts = "scripts"
dir_data = "../data"
dir_out = "../outputs"
dir_timestamps = "../timestamps"


# variables

sample_ids = ["16030X2", "16030X3", "16030X4"]

dirs_fastq = {"16030X2": dir_data + "/16030R/Fastq/16030X2_HJVMLDMXX", 
              "16030X3": dir_data + "/16030R/Fastq/16030X3_HJTWLDMXX", 
              "16030X4": dir_data + "/16030R/Fastq/16030X4_HJTWLDMXX"}

dir_ref = dir_data + "/GRCh38/refdata-cellranger-GRCh38-3.0.0"


# number of cores and memory usage for parameters inside scripts

n_cores_cellranger = 10
n_cores_cellsnp = 20
mem_cellranger = 50

# number of cores and memory usage for Linux cluster settings (total memory = memory * cores)

local_cellranger = 10
mem_free_cellranger = "10G"
h_vmem_cellranger = "20G"
h_fsize_cellranger = "200G"


# note: shell commands in Snakemake rules use "qsub" syntax for Linux cluster job scheduling



# default rule (run all)

rule all:
  input:
    expand(dir_timestamps + "/timestamp_cellranger_{sample}.txt", sample = sample_ids)


# run cellranger count

rule run_cellranger:
  input:
    script_cellranger = dir_scripts + "/run_cellranger.sh"
  output:
    dir_timestamps + "/timestamp_cellranger_{sample}.txt"
  params:
    sample = "{wildcards.sample}", 
    dir_fastq = lambda wildcards: DIRS_FASTQ[wildcards.sample], 
    dir_ref = dir_ref, 
    cores = n_cores_cellranger, 
    mem = mem_cellranger, 
    path_timestamp = dir_timestamps + "/timestamp_cellranger_{wildcards.sample}.txt"
  shell:
    """
    CWD=$(pwd)
    mkdir -p dir_out && cd dir_out && mkdir -p cellranger && cd cellranger
    qsub -V -cwd -pe local {local_cellranger} -l mem_free={mem_free_cellranger},h_vmem={h_vmem_cellranger},h_fsize={h_fsize_cellranger} {input.script_cellranger} {params.sample} {params.dir_fastq} {params.dir_ref} {params.cores} {params.mem}
    cd CWD & mkdir -p dir_timestamps
    date > {params.path_timestamp}
    """


