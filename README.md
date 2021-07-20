# snp-dmx-cancer

This repository contains reproducible code for the `Snakemake` workflow, benchmark evaluations, and supplementary analyses in our paper evaluating genetic variation-based demultiplexing tools (in particular `Vireo`) for pooled single-cell RNA sequencing samples in cancer (high-grade serous ovarian cancer (HGSOC) and lung adenocarcinoma).


## Paper

For details on the analyses, see our paper: [Weber et al. (2021), "Genetic demultiplexing of pooled single-cell RNA-sequencing samples in cancer facilitates effective experimental design", bioRxiv](https://www.biorxiv.org/content/10.1101/2020.11.06.371963v3)


## Contents

### Workflow

The `Snakemake` workflow implements a complete workflow for one dataset and doublets simulation scenario (HGSOC dataset, 20% doublets), using the best-performing set of tools (bulk RNA-seq samples genotype from `bcftools`, demultiplexing using `cellSNP`/`Vireo`) from the benchmark. The workflow is modular, and can be adapted to substitute alternative tools.

Scripts for the `Snakemake` workflow are saved in [workflow/](workflow/):

- [workflow/Snakefile](workflow/Snakefile): `Snakefile` defining the `Snakemake` workflow
- [workflow/scripts/](workflow/scripts/): scripts for the individual steps in the workflow


### Benchmark evaluations

Scripts for the benchmark evaluations are saved in [benchmarking/](benchmarking/):

- [benchmarking/Snakefile](benchmarking/Snakefile): `Snakefile` for the initial steps in the workflow for the benchmark evaluations
- [benchmarking/parse_BAM_doublets/](benchmarking/parse_BAM_doublets/): scripts for doublets scenarios
- [benchmarking/run_scenarios/](benchmarking/run_scenarios/): scripts to run tools for each benchmark scenario
- [evaluations/](evaluations/): scripts to calculate performance evaluations (precision and recall) and generate plots (precision-recall, runtimes, estimated cost savings)


### Supplementary analyses

Scripts for the supplementary analyses are saved in the following directories:

- [supplementary_debris/](supplementary_debris/): scripts for additional scenarios including 10%, 20%, or 40% ambient RNA from simulated cell debris or lysed cells
- [supplementary_snparray/](supplementary_snparray/): scripts for additional scenarios using subset of SNPs from SNP array as genotype reference
- [supplementary_healthy/](supplementary_healthy/): scripts for additional scenarios for healthy (non-cancer) cell line dataset
- [supplementary_downstream/](supplementary_downstream/): scripts to run downstream doublet detection tool


### Additional scripts

Scripts for additional steps outside the main workflow and benchmark evaluations:

- [alternative/](alternative/): scripts for alternative tools that were not used in the final workflow, which may be useful in the future (e.g. `salmon alevin` instead of `Cell Ranger`)
- [download_EGA/](download_EGA/): script to download data files for lung adenocarcinoma dataset (Kim et al. 2020) from European Genome-phenome Archive (EGA) (requires access to the controlled access data repository)
- [download_souporcell/](download_souporcell/): scripts to download data files for healthy (non-cancer) iPSC cell line dataset from souporcell paper (Heaton et al. 2020) from European Nucleotide Archive (ENA)
- [filter_vcf/](filter_vcf/): script to filter 1000 Genomes Project genotype VCF file to retain only SNPs in 3' untranslated regions (UTRs), for faster runtime
- [genotype/](genotype/): scripts to run different options of tools to generate custom genotype VCF file, including from matched bulk RNA-seq samples (using either `bcftools` or `cellSNP`), or directly from single-cell RNA-seq samples (using `cellSNP`)


## Links

- `Vireo` paper: [Huang et al. (2019), Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1865-2)
- [Vireo documentation](https://vireosnp.readthedocs.io/en/latest/index.html)
- [cellsnp-lite documentation](https://github.com/single-cell-genetics/cellsnp-lite)

