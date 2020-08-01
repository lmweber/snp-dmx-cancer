# snp-dmx-cancer

This repository contains a `Snakemake` workflow for demultiplexing single-cell RNA sequencing samples using genetic variation based demultiplexing tools (`Vireo` and `cellSNP`). We demonstrate that these tools can be applied to cancer data (HGSOC and lung adenocarcinoma).


## Contents

### Pipeline

Files for the main workflow are saved in the directory [pipeline/](pipeline/).

- [Snakefile](pipeline/Snakefile): `Snakefile` to run the workflow using `Snakemake`
- [scripts/](pipeline/scripts/): directory containing scripts for the individual steps in the `Snakemake` workflow


### Additional scripts

Scripts for additional steps outside the main workflow.

- [evaluations/](evaluations/): R scripts for performance evaluations on our datasets
- [filter_vcf/](filter_vcf/): to do
- [download_EGA/](download_EGA/): to do


### Alternatives

Scripts for alternative options or tools to use within the workflow. This mainly includes scripts to use `salmon alevin` instead of `Cell Ranger`. We used `Cell Ranger` in the final pipeline due to incompatibility of the `salmon alevin` outputs with `cellSNP` and `Vireo` (`salmon alevin` generates a transcriptomic BAM, while `cellSNP` and `Vireo` expect a genomic BAM and VCF), so these scripts are not used within the main pipeline. However, we have kept them here in case they are useful in the future or for other work.

These scripts are saved in the directory [alternative](alternative/).

- [salmon_alevin/](alternative/salmon_alevin/): alternative scripts to run `salmon alevin` instead of `Cell Ranger`
- [convert_VCF/](alternative/convert_VCF/): alternative scripts to convert VCF file from genomic coordinates to transcriptomic coordinates


## Links

- Vireo: [documentation](https://vireosnp.readthedocs.io/en/latest/index.html)
- Vireo paper: [Huang et al. (2019), Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1865-2)
- cellSNP: [documentation](https://github.com/single-cell-genetics/cellSNP)
- Cell Ranger: [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome)

