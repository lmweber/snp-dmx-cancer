# Run Vireo

# see https://vireosnp.readthedocs.io/en/latest/manual.html


# notes:
# run Vireo using option 1 (without any sample genotyping)
# using cell genotyping output from run_cellSNP.sh (also in mode 1)


#qsub -V -cwd -l mem_free=100G,h_vmem=200G,h_fsize=200G run_vireo.sh

vireo -c out_cellSNP -N 3 -o out_vireo


