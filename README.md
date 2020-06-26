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


### 03: Explaining mean DNA methylation
Step-by-step exploring features, that explain mean DNA methylation level per species
```bash
Rscript 03.1_prediction_features.R ##creating a table of features
Rscript 03.2_kmer_predictability_by_group.R ##prediction values for each group, takes a couple of hours
Rscript 03.3_kmer3_stability.R ##stepwise selection on 3mers, takes a while
```
### 04: Reconstructing phylogeny
Using kmer frequencies to create a hierarchical clustering of species, based on their sequence composition
```bash
Rscript 04.1_phylogeny_reconstruction.R
```
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
Identifying the species, that we classify as "inverted":
```bash
Rscript 06.1_inverted.R
```
Than, we confirm the observed effect with additional experiments on subsets of data in a tissue-specific and via performing subsampling:
```bash
sbatch sbatch/06.2.1.sbatch <candidate species>
sbatch sbatch/06.2.2.sbatch <candidate species>
sbatch sbatch/06.2.3.sbatch <candidate species>
Rscript 06.2_verify_inv_summary.R
```
Analysing the frequencies of kmers specifically in high/low methylated sequences:
!!first two scripts require additional specification of the data_dir path
```bash
python 06.3.1_to_fasta_for_kmer_count.py
sbatch 06.3.2.sbatch #runs on all the species, where classes contain inverted
RScript 06.3_freq_kmers.R
```
And summarizing all the knowledge, oriented on the WHH species
```
Rscript 06.4_explore_inverted_final_WHH.R
```

### 07: TFBS analysis across species and tissues
As a preparation, the JASPAR matrix of all the TF profiles has to be downloaded: http://jaspar.genereg.net/downloads/
We used the non-redundant vertebrate TFBS set.
First, for each species we calculate frequencies of each TFBS in deduced reference fragments.
```bash
sbatch sbatch/07.1.sbatch
sbatch sbatch/07.2.sbatch
Rscript 07.3_TFBS_summary.R
```

Than we run differential methylation analysis and motif enrichment in tissue-specific manner:
(make sure you have ame installed)
```bash
sbatch sbatch/07.4.1.sbatch
sbatch sbatch/07.4.2.sbatch
Rscript 07.4.ame_summary.R
```

Than on the enriched we are running specific analysis:
For integrating the methylation preference specificity step, we used the data from (reference and how we share it? link?)
```bash
Rscript 07.5_ame_selex.R
Rscript 07.6_GO.R
Rscript 07.7_GO.R

### 08: Tissue-specific DNA methylation
<Johanna?>
