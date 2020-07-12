###########################
# Snakefile to run pipeline
###########################

# ---------
# variables
# ---------

dir_scripts = "scripts"
dir_data = "../data"
dir_data_HGSOC = dir_data + "/HGSOC"
dir_outputs = "../outputs"
dir_outputs_HGSOC = dir_outputs + "/HGSOC"
dir_runtimes = "../runtimes"
dir_runtimes_HGSOC = dir_runtimes + "/HGSOC"
dir_timestamps = "../timestamps"
dir_timestamps_HGSOC = dir_timestamps + "/HGSOC"


sample_ids_HGSOC = ["16030X2_HJVMLDMXX", "16030X3_HJTWLDMXX", "16030X4_HJTWLDMXX"]

sample_ids_HGSOC_short = ["X2", "X3", "X4"]

sample_ids_HGSOC_short_named = {"16030X2_HJVMLDMXX": "X2", 
                                "16030X3_HJTWLDMXX": "X3", 
                                "16030X4_HJTWLDMXX": "X4"}

fastq_HGSOC = {"16030X2_HJVMLDMXX": [dir_data_HGSOC + "/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R1_001.fastq.gz", 
                                     dir_data_HGSOC + "/16030R/Fastq/16030X2_HJVMLDMXX/16030X2_HJVMLDMXX_S1_L002_R2_001.fastq.gz"], 
               "16030X3_HJTWLDMXX": [dir_data_HGSOC + "/16030R/Fastq/16030X3_HJTWLDMXX/16030X3_HJTWLDMXX_S1_L001_R1_001.fastq.gz", 
                                     dir_data_HGSOC + "/16030R/Fastq/16030X3_HJTWLDMXX/16030X3_HJTWLDMXX_S1_L001_R2_001.fastq.gz"], 
               "16030X4_HJTWLDMXX": [dir_data_HGSOC + "/16030R/Fastq/16030X4_HJTWLDMXX/16030X4_HJTWLDMXX_S2_L001_R1_001.fastq.gz", 
                                     dir_data_HGSOC + "/16030R/Fastq/16030X4_HJTWLDMXX/16030X4_HJTWLDMXX_S2_L001_R2_001.fastq.gz"]}


# ------------
# run pipeline
# ------------

# command to run pipeline on cluster

# snakemake --cluster "qsub -V -cwd -pe local 10 -l mem_free=10G,h_vmem=20G,h_fsize=300G" -j 6 --local-cores 20 --latency-wait 10


# --------------
# pipeline rules
# --------------

# default rule

rule all:
  input:
    dir_timestamps_HGSOC + "/parse_and_merge_barcodes/timestamp_parse_and_merge_barcodes.txt"


# parse and merge cell barcode files

rule parse_and_merge_barcodes:
  input:
    dir_timestamps_HGSOC + "/merge_and_index_BAM/timestamp_merge_and_index_BAM.txt", 
    script_parse_and_merge_barcodes = dir_scripts + "/parse_and_merge_barcodes.sh"
  output:
    dir_timestamps_HGSOC + "/parse_and_merge_barcodes/timestamp_parse_and_merge_barcodes.txt"
  params:
    sample_id_1 = lambda wildcards: sample_ids_HGSOC[0], 
    sample_id_2 = lambda wildcards: sample_ids_HGSOC[1], 
    sample_id_3 = lambda wildcards: sample_ids_HGSOC[2], 
    sample_id_1_short = lambda wildcards: sample_ids_HGSOC_short[0], 
    sample_id_2_short = lambda wildcards: sample_ids_HGSOC_short[1], 
    sample_id_3_short = lambda wildcards: sample_ids_HGSOC_short[2]
  shell:
    "bash {input.script_parse_and_merge_barcodes} {dir_runtimes_HGSOC} {dir_timestamps_HGSOC} NA "
    "{dir_outputs_HGSOC} "
    "{params.sample_id_1} {params.sample_id_2} {params.sample_id_3} "
    "{params.sample_id_1_short} {params.sample_id_2_short} {params.sample_id_3_short}"


# merge and index BAM files

rule merge_and_index_BAM:
  input:
    expand(dir_timestamps_HGSOC + "/convert_BAM/timestamp_convert_BAM_{sample}.txt", sample = sample_ids_HGSOC), 
    script_merge_and_index_BAM = dir_scripts + "/merge_and_index_BAM.sh"
  output:
    dir_timestamps_HGSOC + "/merge_and_index_BAM/timestamp_merge_and_index_BAM.txt"
  params:
    sample_id_1 = lambda wildcards: sample_ids_HGSOC[0], 
    sample_id_2 = lambda wildcards: sample_ids_HGSOC[1], 
    sample_id_3 = lambda wildcards: sample_ids_HGSOC[2]
  shell:
    "bash {input.script_merge_and_index_BAM} {dir_runtimes_HGSOC} {dir_timestamps_HGSOC} NA "
    "{dir_outputs_HGSOC} {params.sample_id_1} {params.sample_id_2} {params.sample_id_3}"


# convert SAM to BAM files

rule convert_BAM:
  input:
    dir_timestamps_HGSOC + "/parse_SAM_barcodes/timestamp_parse_SAM_barcodes_{sample}.txt", 
    script_convert_BAM = dir_scripts + "/convert_BAM.sh"
  output:
    dir_timestamps_HGSOC + "/convert_BAM/timestamp_convert_BAM_{sample}.txt"
  shell:
    "bash {input.script_convert_BAM} {dir_runtimes_HGSOC} {dir_timestamps_HGSOC} NA "
    "{wildcards.sample} {dir_outputs_HGSOC}"


# parse SAM files to add sample IDs to cell barcodes

rule parse_SAM_barcodes:
  input:
    dir_timestamps_HGSOC + "/alevin/timestamp_alevin_{sample}.txt", 
    script_parse_SAM_barcodes = dir_scripts + "/parse_SAM_cell_barcodes.sh"
  output:
    dir_timestamps_HGSOC + "/parse_SAM_barcodes/timestamp_parse_SAM_barcodes_{sample}.txt"
  params:
    short_sample_id = lambda wildcards: sample_ids_HGSOC_short_named[wildcards.sample]
  shell:
    "bash {input.script_parse_SAM_barcodes} {dir_runtimes_HGSOC} {dir_timestamps_HGSOC} NA "
    "{wildcards.sample} {params.short_sample_id} {dir_outputs_HGSOC}"


# run salmon alevin

rule run_alevin:
  input:
    dir_timestamps + "/salmon_index/timestamp_salmon_index.txt", 
    script_alevin = dir_scripts + "/run_salmon_alevin.sh"
  output:
    dir_timestamps_HGSOC + "/alevin/timestamp_alevin_{sample}.txt"
  params:
    fastq_1 = lambda wildcards: fastq_HGSOC[wildcards.sample][0], 
    fastq_2 = lambda wildcards: fastq_HGSOC[wildcards.sample][1]
  shell:
    "bash {input.script_alevin} {dir_runtimes_HGSOC} {dir_timestamps_HGSOC} 10 "
    "{wildcards.sample} {params.fastq_1} {params.fastq_2} "
    "{dir_data}/salmon_index {dir_outputs_HGSOC}"


# run salmon index

rule run_salmon_index:
  input:
    script_salmon_index = dir_scripts + "/run_salmon_index.sh"
  output:
    dir_timestamps + "/salmon_index/timestamp_salmon_index.txt"
  shell:
    "bash {input.script_salmon_index} {dir_runtimes} {dir_timestamps} 10 NA {dir_data}/salmon_index"

