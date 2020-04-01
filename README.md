# HGSOC-Vireo

Repository containing scripts for Vireo analyses to demultiplex combined scRNA-seq samples in our HGSOC dataset.


## Scripts

Order of scripts:

- Preprocessing:

    - `run_cellranger_16030X2.sh`, `run_cellranger_16030X3.sh`, `run_cellranger_16030X4.sh`: run Cell Ranger to generate BAM file for each scRNA-seq sample
    - `run_STAR_index.sh`, `run_STARsolo_16030X2.sh`, `run_STARsolo_16030X3.sh`, `run_STARsolo_16030X4.sh`: alternative scripts to run STARsolo instead of Cell Ranger due to problems when running Cell Ranger on JHPCE server (now redundant since Cell Ranger scripts fixed, but may be useful in the future since STARsolo is much faster)
    - `combine_BAM_files_CellRanger.sh`: combine BAM files from Cell Ranger into one multiplexed BAM file
    - `combine_BAM_files_STARsolo.sh`: alternatively, combine BAM files from STARsolo (now redundant)

- Vireo:

    - `run_cellSNP.sh`: run CellSNP to genotype cells


