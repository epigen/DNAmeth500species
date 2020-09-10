# DNAmeth500species
Comparative analysis of DNA methylation across more than 500 animal species

In general, scripts can be divided into two categories: Those that have to be rn in HPC mode on each species in parallel, andthose that are used for summarizing and vizualizing results by generating summary tables and figure panels. For the HPC scrips sbatch scripts with reccomended parametrs to run the job are provided but need to be adjusted if a sheduler other than SLURM is to be used.

## Reproducing the project step-by-step
### Adapting to your setup
00.0_init.R is the script, to initialize all further scripts, it will standardize setting, make annotations availble, and set paths to data and scripts. Few absolute paths need to be adjusted and two enviromnental variables are required:  
CODEBASE (location of this github repository)
LOGDIR (location for pipelen log files)


### 00: preparing to run the pipeline on raw data
Preparing to run RefFreeDMA, if you want to start from the raw data. This will prepare the run folders and shell scripts containing the commands to run RefFeeDMA.
If run in interactive mode stats can be checked and parameters can be adjusted.
```bash
Rscript 00.1_prepareRefFreeDMA.R 
```
### 01: creating a comprehensive sample annotation table
Following scripts will generate some additional data metrics, that are used in most subsequent analysis.
Here we create the big summary annotation table and other initial stats:
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
!Importantly, the path in 01.5.sbatch has to be adjusted
```bash
sbatch sbatch/01.5.1.sbatch ##prefilteing
sbatch sbatch/01.5.sbatch
Rscript 01.5.2_kmer_count_summary.R  ##summary for further use
```

Generation of the phylogenetic tree for iTOL, it's optimization for future analysis and exporting the taxonomimc order of species for use in subsequent analysis:
```bash
Rscript 01.6_create_ITOL.R
Rscript 01.7_species_order.R
```

Calculating the CpG island properties:
```bash
sbatch sbatch/01.8.sbatch 
```

### 02: Main statisctics
Following scripts create summary statistics and generate summary plots, giving us an overview of the data  
02.1_stats_uc.R  
02.2_stats_general.R  
02.3_stats_detail.R  
02.4_stats_location.R  
02.5_stats_methylation.R  

Saving the supplementary tables and the stats for the paper
```bash
Rscript 02.6_save_annot.R
```


### 03: Explaining mean DNA methylation
Step-by-step exploration of features, that explain mean DNA methylation level per species
```bash
Rscript 03.1_prediction_features.R ##creating a table of features
Rscript 03.2_kmer_predictability_by_group.R ##prediction values for each group, takes a couple of hours
Rscript 03.3_kmer3_stability.R ##stepwise selection on 3mers, takes a while
```
### 04: Reconstructing phylogeny
Using kmer frequencies to create a hierarchical clustering of species, based on the sequence composition of the consensus genomes
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
Additional check with more stringent thresholds
```bash
sbatch sbatch/05.1.1.sbatch
Rscript 05.2.1_summary_alt.R
```
Cross-species prediction
```bash
sbatch sbatch/05.4.sbatch
Rscript 05.5_summary_cross_species_prediction.R
```

### 06: Investigating the inverted prediction phenotype
Identifying the species, that we classify as "inverted":
```bash
Rscript 06.1_inverted.R
```
Than, we confirm the observed effect with additional experiments on subsets of data (by tissue, individual, and random subsampling)
```bash
sbatch sbatch/06.2.1.sbatch <candidate species>
sbatch sbatch/06.2.2.sbatch <candidate species>
sbatch sbatch/06.2.3.sbatch <candidate species>
Rscript 06.2_verify_inv_summary.R
```
Analysing the frequencies of kmers specifically in highly/lowly methylated sequences:
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
Calculate and report frequencies of each TFBS in consensus reference of each species.
```bash
sbatch sbatch/07.1.sbatch
sbatch sbatch/07.2.sbatch
Rscript 07.3_TFBS_summary.R
```

Identification of differentially methylated fragemnts between heart and liver tissues and motif enrichment analysis:
(make sure you have ame installed)
```bash
sbatch sbatch/07.4.1.sbatch
sbatch sbatch/07.4.2.sbatch
Rscript 07.4.ame_summary.R
```

Follow-up on enriched TFBS:
For integrating the methylation preference specificity step, we used the data from (reference and how we share it? link?)
```bash
Rscript 07.5_ame_selex.R
Rscript 07.6_GO.R
Rscript 07.7_GO.R
```
### 08: Investigate sources of within-species DNA methylation variability
Multidimensional scaling analysis for selected species
```bash
Rscript 08.1_tissue_MDS.R
```
Comparison of DNA methylation variability due to tissue and individual
```bash
Rscript 08.2_tissue_correlation.R
Rscript 08.3_tissue_distance.R
```
08.4_tissue_diffMeth.R
