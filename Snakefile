###########################
# Snakefile to run pipeline
###########################

# variables

dir_scripts = "scripts"
dir_data = "../data"
dir_outputs = "../outputs"
dir_runtimes = "../runtimes"
dir_timestamps = "../timestamps"


sample_ids_HGSOC = ["16030X2_HJVMLDMXX", "16030X3_HJTWLDMXX", "16030X4_HJTWLDMXX"]

fastq_HGSOC = {"16030X2_HJVMLDMXX": [dir_data + "/HGSOC/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R1_001.fastq.gz", 
                                     dir_data + "/HGSOC/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R2_001.fastq.gz"], 
               "16030X3_HJTWLDMXX": [dir_data + "/HGSOC/16030R/Fastq/16030X3_HJTWLDMXX/16030X3_HJTWLDMXX_S1_L001_R1_001.fastq.gz", 
                                     dir_data + "/HGSOC/16030R/Fastq/16030X3_HJTWLDMXX/16030X3_HJTWLDMXX_S1_L001_R2_001.fastq.gz"], 
               "16030X4_HJTWLDMXX": [dir_data + "/HGSOC/16030R/Fastq/16030X4_HJTWLDMXX/16030X4_HJTWLDMXX_S2_L001_R1_001.fastq.gz", 
                                     dir_data + "/HGSOC/16030R/Fastq/16030X4_HJTWLDMXX/16030X4_HJTWLDMXX_S2_L001_R2_001.fastq.gz"]}


# ------------
# run pipeline
# ------------

# command to run pipeline on cluster

# snakemake --cluster "qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=200G" -j 6 --local-cores 20


# --------------
# pipeline rules
# --------------

# default rule

rule all:
  input:
    expand(dir_timestamps + "/HGSOC/alevin/timestamp_alevin_{sample}.txt", sample = sample_ids_HGSOC)


# run salmon alevin

rule run_alevin:
  input:
    script_alevin = dir_scripts + "/run_salmon_alevin.sh", 
    salmon_index_timestamp = dir_timestamps + "/salmon_index/timestamp_salmon_index.txt"
  output:
    dir_timestamps + "/HGSOC/alevin/timestamp_alevin_{sample}.txt"
  params:
    fastq_1 = lambda wildcards: fastq_HGSOC[wildcards.sample][0], 
    fastq_2 = lambda wildcards: fastq_HGSOC[wildcards.sample][1]
  shell:
    "bash {input.script_alevin} {wildcards.sample} {dir_runtimes} {dir_timestamps} 10 {params.fastq_1} {params.fastq_2} {dir_data}/salmon_index {dir_outputs}/HGSOC/{wildcards.sample}"


# run salmon index

rule run_salmon_index:
  input:
    script_salmon_index = dir_scripts + "/run_salmon_index.sh"
  output:
    dir_timestamps + "/salmon_index/timestamp_salmon_index.txt"
  shell:
    "bash {input.script_salmon_index} NA {dir_runtimes} {dir_timestamps} 10 {dir_data}/salmon_index"

