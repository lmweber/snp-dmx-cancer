# snp-dmx-cancer

This repository contains a `Snakemake` workflow for demultiplexing single-cell RNA sequencing samples using genetic variation based demultiplexing tools (`Vireo` and `cellSNP`). We demonstrate that these tools can be applied to cancer data (HGSOC and lung adenocarcinoma).


## Contents

#### Main pipeline

- [Snakefile](Snakefile): `Snakefile` for main `Snakemake` workflow
- [scripts/](scripts/): directory containing scripts for main `Snakemake` workflow


#### Additional scripts

- [filter_vcf/](filter_vcf/): to do
- [evaluations/](evaluations/): scripts to evaluate performance on our datasets
- [download_EGA/](download_EGA/): to do
- [salmon_alevin/](salmon_alevin/): alternative scripts for running `salmon alevin` instead of `Cell Ranger` (using `salmon alevin` in the pipeline currently does not work, since `cellSNP` expects a genomic BAM instead of transcriptomic BAM; however we have kept these scripts here in case they are useful in the future to update the pipeline to use `salmon alevin`)
- [convert_VCF/](convert_VCF/): alternative scripts to convert the VCF file for `cellSNP` from genomic coordinates to transcriptomic coordinates (these scripts may be useful in the future to update the pipeline to use `salmon alevin`)


## Links

- Vireo: [documentation](https://vireosnp.readthedocs.io/en/latest/index.html)
- Vireo paper: [Huang et al. (2019), Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1865-2)
- cellSNP: [documentation](https://github.com/single-cell-genetics/cellSNP)
- Cell Ranger: [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome)

