# DNAmeth500species
Comparative analysis of DNA methylation across more than 500 animal species

In general, scripts can be divided into two categories - ones, that have to be ran in HPC mode on each species in parallel, and summary/vizualization scripts, that generate summary tables and figures. In first case the corresponding sbatch script with reccomended parametrs to run the job is provided. (those need to be updated, depending on the comutational setup)

## Reproducing the project step-by-step
### Adapting to your setup
00.0_init.R is the script, to initialize all further scripts, it will standardize setting, make annotations availble, and set paths to data and scripts. Implicit paths should be changed in the user - specific manner in the corresponding section.
Additionaly, two enviromnental variables are required:  
CODEBASE
LOGDIR
### 00: preparing to run the pipeline on raw data

Preparing to run RefFreeDMA, if you want to repeat the steps from the raw data
```bash
Rscript 00.1_prepareRefFreeDMA.R 
```
### 01: calculating different metrics
Following scripts will generate some additional data metrics, that we can use in the future.
Creating the big summary annotation table and other initial stats:
```bash
Rscript 01.1_collectStats.R
RScript 01.2_taxonomicAnnotation.R. #generates the stats annot table
RScript 01.3_quality_stratification.R
```

Calcualating the frequency of CGs for future analysis:
```
Rscript 01.4_dedRef_CpG_count.R  ##takes a while
```
Calculating the frequency of kmers, using the MEME tool. 
!Importantly, the path in 01.5.sbatch has to be updated additionally
```bash
sbatch sbatch/01.5.1.sbatch ##prefilteing
sbatch sbatch/01.5.sbatch
Rscript 01.5.2_kmer_count_summary.R  ##summary for further use
```

Generation of the phylogenetic tree for iTOl, it's optimization for future analysis and exporting the phylogenetic order of species for future use:
```bash
Rscript 01.6_create_ITOL.R
Rscript 01.7_species_order.R
```

Calculating the CpG island criteria:
```bash
sbatch sbatch/01.8.sbatch 
```

### 02: Main statisctics
Following scripts perform statistical analysis and generate summary plots,giving us an overview of the data
02.1_stats_uc.R  
02.2_stats_general.R  
02.3_stats_detail.R  
02.4_stats_location.R  
02.5_stats_methylation.R  

Saving the supplementary tables and the stats for the data
```bash
Rscript 02.6_save_annot.R
```
### 03: Recinstructing phylogeny
### 04: Explaining mean DNA methylation
### 05: Prediction of DNA methylation from sequence
With-in species DNA methylation prediction
```bash
sbatch sbatch/05.1.sbatch: repeated on k 1:10, k = 3, k = 4, k = 7
Rscript 05.2_summarize_within_sp_prediction.R
Rscript 05.3_feature_weights_analysis.R
```
Additional check with more stringent thresholds:
```bash
sbatch sbatch/05.1.1.sbatch
Rscript 05.2.1_summary_alt.R
```
Cross-species prediction:
```bash
sbatch sbatch/05.4.sbatch
Rscript 05.5_summary_cross_species_prediction.R
```

### 06: Investigating the invertiness
### 07: TFBS analysis across species and tissues
### 08: Tissue-specific DNA methylation

