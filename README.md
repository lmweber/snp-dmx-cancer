# snp-dmx-cancer

This repository contains reproducible code for the benchmark evaluations and `Snakemake` workflow in our paper evaluating the performance of genetic variation-based demultiplexing tools (in particular `Vireo`) for pooled single-cell RNA sequencing samples in cancer (high-grade serous ovarian cancer (HGSOC) and lung adenocarcinoma).

The benchmark evaluations cover a range of simulated doublet proportions (up to 30% for "super-loading" experimental designs) and choice of genotype reference list of SNPs (e.g. 1000 Genomes Project, or genotyping from matched bulk RNA-seq samples).

The `Snakemake` workflow implements a complete workflow for one dataset and doublets scenario (HGSOC dataset, 20% doublets), using the best-performing set of tools (bulk RNA-seq samples genotype from `bcftools`; demultiplexing using `cellSNP`/`Vireo`). The workflow is modular, and can be adapted to substitute alternative tools.

For more details, see our paper (preprint to be posted on bioRxiv soon).


## Contents

### Workflow

Scripts for the `Snakemake` workflow are saved in [workflow/](workflow/):

- [workflow/Snakefile](workflow/Snakefile): `Snakefile` defining the `Snakemake` workflow
- [workflow/scripts/](workflow/scripts/): scripts for the individual steps in the workflow


### Benchmark evaluations

Scripts for the benchmark evaluations are saved in [benchmarking/](benchmarking/):

- [benchmarking/Snakefile](benchmarking/Snakefile): additional `Snakefile` for the initial steps in the workflow for the benchmark evaluations
- [benchmarking/parse_BAM_doublets/](benchmarking/parse_BAM_doublets/): scripts for doublets scenarios
- [benchmarking/run_scenarios/](benchmarking/run_scenarios/): scripts to run tools for each benchmark scenario
- [evaluations/](evaluations/): scripts to calculate performance evaluations (precision and recall) and generate plots (precision-recall, runtimes, estimated cost savings)


### Additional scripts

Scripts for additional steps outside the main workflow and benchmark evaluations:

- [alternative](alternative/): scripts for alternative tools that were not used in the final workflow, which may be useful in the future (e.g. `salmon alevin` instead of `Cell Ranger`)
- [download_EGA/](download_EGA/): script to download data files for lung adenocarcinoma dataset from European Genome-phenome Archive (EGA) (requires access to the controlled access data repository)
- [filter_vcf/](filter_vcf/): script to filter 1000 Genomes Project genotype VCF file to retain only SNPs in 3' untranslated region (UTR), for much faster runtime
- [genotype/](genotype/): scripts to run different options of tools to generate custom genotype VCF file, including from matched bulk RNA-seq samples (using either `bcftools` or `cellSNP`), or directly from single-cell RNA-seq samples (using `cellSNP`)


## Links

- Vireo paper: [Huang et al. (2019), Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1865-2)
- Vireo: [documentation](https://vireosnp.readthedocs.io/en/latest/index.html)
- cellSNP: [documentation](https://github.com/single-cell-genetics/cellSNP)

