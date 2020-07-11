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

short_sample_ids_HGSOC = {"16030X2_HJVMLDMXX": "X2", 
                          "16030X3_HJTWLDMXX": "X3", 
                          "16030X4_HJTWLDMXX": "X4"}


# ------------
# run pipeline
# ------------

# command to run pipeline on cluster

# snakemake --cluster "qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G" -j 6 --local-cores 20


# --------------
# pipeline rules
# --------------

# default rule

rule all:
  input:
    dir_timestamps + "/HGSOC/merge_and_index_BAM/timestamp_merge_and_index_BAM.txt"
    


# merge and index BAM files

rule merge_and_index_BAM:
  input:
    expand(dir_timestamps + "/HGSOC/convert_BAM/timestamp_convert_BAM_{sample}.txt", sample = sample_ids_HGSOC), 
    script_merge_and_index_BAM = dir_scripts + "/merge_and_index_BAM.sh"
  output:
    dir_timestamps + "/HGSOC/merge_and_index_BAM/timestamp_merge_and_index_BAM.txt"
  params:
    sample_id_1 = lambda wildcards: sample_ids_HGSOC[0], 
    sample_id_2 = lambda wildcards: sample_ids_HGSOC[1], 
    sample_id_3 = lambda wildcards: sample_ids_HGSOC[2]
  shell:
    "bash {input.script_merge_and_index_BAM} {dir_runtimes} {dir_timestamps} NA {dir_outputs}/HGSOC {params.sample_id_1} {params.sample_id_2} {params.sample_id_2}"


# convert SAM to BAM files

rule convert_BAM:
  input:
    dir_timestamps + "/HGSOC/parse_SAM_barcodes/timestamp_parse_SAM_barcodes_{sample}.txt", 
    script_convert_BAM = dir_scripts + "/convert_BAM.sh"
  output:
    dir_timestamps + "/HGSOC/convert_BAM/timestamp_convert_BAM_{sample}.txt"
  shell:
    "bash {input.script_convert_BAM} {dir_runtimes} {dir_timestamps} NA {wildcards.sample} {dir_outputs}/HGSOC/{wildcards.sample}/alevin_mappings"


# parse SAM files to add sample IDs to cell barcodes

rule parse_SAM_barcodes:
  input:
    dir_timestamps + "/HGSOC/alevin/timestamp_alevin_{sample}.txt", 
    script_parse_SAM_barcodes = dir_scripts + "/parse_SAM_cell_barcodes.sh"
  output:
    dir_timestamps + "/HGSOC/parse_SAM_barcodes/timestamp_parse_SAM_barcodes_{sample}.txt"
  params:
    short_sample_id = lambda wildcards: short_sample_ids_HGSOC[wildcards.sample]
  shell:
    "bash {input.script_parse_SAM_barcodes} {dir_runtimes} {dir_timestamps} NA {wildcards.sample} {params.short_sample_id} {dir_outputs}/HGSOC/{wildcards.sample}/alevin_mappings"


# run salmon alevin

rule run_alevin:
  input:
    dir_timestamps + "/salmon_index/timestamp_salmon_index.txt", 
    script_alevin = dir_scripts + "/run_salmon_alevin.sh"
  output:
    dir_timestamps + "/HGSOC/alevin/timestamp_alevin_{sample}.txt"
  params:
    fastq_1 = lambda wildcards: fastq_HGSOC[wildcards.sample][0], 
    fastq_2 = lambda wildcards: fastq_HGSOC[wildcards.sample][1]
  shell:
    "bash {input.script_alevin} {dir_runtimes} {dir_timestamps} 10 {wildcards.sample} {params.fastq_1} {params.fastq_2} {dir_data}/salmon_index {dir_outputs}/HGSOC/{wildcards.sample}"


# run salmon index

rule run_salmon_index:
  input:
    script_salmon_index = dir_scripts + "/run_salmon_index.sh"
  output:
    dir_timestamps + "/salmon_index/timestamp_salmon_index.txt"
  shell:
    "bash {input.script_salmon_index} {dir_runtimes} {dir_timestamps} 10 NA {dir_data}/salmon_index"

