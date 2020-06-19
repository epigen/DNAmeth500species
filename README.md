# DNAmeth500species
Comparative analysis of DNA methylation across more than 500 animal species

1. The init file src/00.0_init.R is user specific and needs to be slightly adapted (one path section)

2. The src/00.0_init.R script to initialize all further scripts, it will standardize setting, make annotations availble, and set paths to data and scripts. First numbers after . indicates the order, in which scripts should be run. 
In general, scripts can be divided into two categories - ones, that have to be ran in HPC mode on each species in parallel, and summary/vizualization scripts, that generate summary tables and figures. In first case the corresponding sbatch script with reccomended parametrs to run the job is provided. (those need to be updated, depending on the comutational setup)

## Reproducting the project step-by-step
### 00: pre-analysis

Preparing to run RefFreeDMA, if you want to repeat the steps from the raw data
```bash
Rscript 00.1_prepareRefFreeDMA.R 
```

the quality control, based on the unconverted samples (for several steps):
```bash
sbatch sbatch/00.21.sbatch
#vizualize
Rscript 00.22__plot_QC.R
```
Generation of the phylogenetic tree for iTOl, it's optimization for future analysis and exporting the phylogenetic order of species for future use:

```bash
Rscript 00.3_create_ITOL.R
Rscript 00.31_species_order.R
```
