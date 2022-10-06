# DNAmeth500species
Comparative analysis of DNA methylation across more than 500 animal species

In general, scripts can be divided into two categories: Those that have to be run in HPC mode on each species in parallel, and those that are used for summarizing and vizualizing results by generating summary tables and figure panels. For the HPC scrips sbatch scripts with reccomended parametrs to run the job. Those scripts have to be adjusted to the HPC and scheduler settings.

## Reproducing the project step-by-step

### Adapting to your setup
00.0_init.R is the script, to initialize all further scripts, it will standardize setting, make annotations availble, and set paths to data and scripts. Few absolute paths need to be adjusted and two enviromnental variables are required:  
```
CODEBASE (location of this github repository) 
LOGDIR (location for pipeline log and tmp files) 
```
#### Preparation
Make sure, you have all the nessesary packages in R:
1. For data processing, set up RefFreeDMA and bisulfiteBlast,a s described in corresponding repositories:
https://github.com/jklughammer/RefFreeDMA 
https://github.com/jklughammer/bisulfiteBlast

2. for the analysis we use two versions of R (3.6.1 for the most part and 4.1 for the SVM modelling). We provide the conda environmets with listed packages as yaml files in the envs folders. The missing packages that were not available via conda and must be installed within R are provided as csv tables.

### 00: Preparing to run the pipeline on raw data

For the first processing step the RefFreeDMA pipeline[https://github.com/jklughammer/RefFreeDMA] is used. The 00.1 script will prepare the run folders and shell scripts containing the commands to run RefFeeDMA. For setting up RefFreeDMA please follow the instructions in the corresponding repository.

```bash
Rscript 00.1_prepareRefFreeDMA.R 
```
### 01: Creating a comprehensive sample annotation table

Following scripts will generate some additional data metrices, that are used in most subsequent analysis.
Here we create the big summary annotation table and other initial stats:
```bash
Rscript 01.1_collectStats.R
RScript 01.2_taxonomicAnnotation.R. #generates the stats annot table
RScript 01.3_quality_stratification.R
```

Calcualating the frequency of CGs for subsequent analysis:
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

### 02: Main statistics
Following scripts create summary statistics and generate summary plots, giving us an overview of the data  
```bash
Rscript 02.1_stats_uc.R  
Rscript 02.2_stats_general.R  
Rscript 02.3_stats_detail.R  
Rscript 02.4_stats_location.R  
Rscript 02.5_stats_methylation.R  
```

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
### 05: Predicting DNA methylation from sequence
With-in species DNA methylation prediction
```bash
sbatch sbatch/05.1.sbatch ##repeated on k 1:10, k = 3, k = 4, k = 7
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
Confirm the observed effect with additional experiments on subsets of data (by tissue, individual, and random sub-sampling)
```bash
sbatch sbatch/06.2.1.sbatch <candidate species>
sbatch sbatch/06.2.2.sbatch <candidate species>
sbatch sbatch/06.2.3.sbatch <candidate species>
Rscript 06.2_verify_inv_summary.R
```
Analysing the frequencies of kmers/repeats specifically in highly/lowly methylated sequences:
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
We used the non-redundant vertebrate TFBS set (additionally download the JASPAR annotation table)
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
For integrating the methylation preference for TF binding, we used the data from (reference and how we share it? link?)
```bash
Rscript 07.5_ame_selex.R
Rscript 07.6_GO.R
Rscript 07.7_GO.R
```
### 08: Investigating sources of within-species DNA methylation variability
Multidimensional scaling analysis for selected species
```bash
Rscript 08.1_tissue_MDS.R
```
Comparison of DNA methylation variability due to tissue and individual
```bash
Rscript 08.2_tissue_correlation.R
Rscript 08.3_tissue_distance.R
```
Create correlation plots between liver and heart DNA methylattion levels for selected species
```bash
08.4_tissue_diffMeth.R
```

## Validation on on reference genomes


### 01 CrossMapping RefFreeDMA analysis to reference genomes

determining which genomes to process, cleaning them up and concatting 
```bash
Rscript 01.1_find_mapping_match.R ## finding pairs of genomes to match
Rscript 01.1.2_axolotl.R ## rearrange and remove small chromosomes of the genome of the axolotl
sbatch 01.2_concat_genome.sbatch ## running genome concat
```
Generating the config files and running the cross-mapping:
```bash
Rscript 01.3_save_configs.R
```
afterwards run the job-subbmitting bash script same way as for original RefFreeDMA

Summarise the results of crossMapping:
```bash
Rscript 01.4_crossMapping_stats.R
sbatch sbatch/01.5_cross-mapping_analysis.sbatch
Rscript 01.5.2_crossMappingAnalysis_summary.R
Rscript 01.5.3_crossMapping_plot.R
```

### 02 Insilico RRBS simulation

In interactive jupyer notebook:
Run the preparation notebook for genomes:
02.1_insilico_digest_prepare_genomes.R

```bash
sbatch 02.2_insilicoDigest.sbatch
Rscript 02.3_insilicoDigest_summary.R
```
### 03 WGBS analysis

WGBS data download is described at 03.1_download_WGBS.ipynb

```bash
Rscript 03.2.process_WGBS_data.R
Rscript 03.2.2_compare_mean_meth.R

sbatch 03.3_WGBS_fragment_CpGmeth.sbatch ##preparing the data
```

##### SVM-modeling in WGBS

```bash
sbatch 03.4.sbatch
sbatch 03.5.sbatch
sbatch 03.6.sbatch

Rscript 03.4.1_WGBS_prediction_summary.R
Rscript 03.5_summary.R
Rscript 03.6_summary.R
```
