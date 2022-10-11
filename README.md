# DNAmeth500species
## Comparative analysis of genome-scale, base-resolution DNA methylation profiles across 580 animal species.

Methylation of cytosines is the prototypic epigenetic modification of the DNA. It has been implicated in various regulatory mechanisms throughout the animal kingdom and particularly in vertebrates. We mapped DNA methylation in 580 animal species (535 vertebrates, 45 invertebrates), resulting in 2443 genome-scale DNA methylation profiles of multiple organs. Reference-genome independent analysis of this comprehensive dataset quantified the association of DNA methylation with the underlying genomic DNA sequence throughout vertebrate evolution. We observed a broadly conserved link with two major transitions – once in the first vertebrates and again with the emergence of reptiles. Cross-species com-parisons focusing on individual organs supported a deeply conserved association of DNA methylation with tissue type, and cross-mapping analysis of DNA methylation at gene promoters revealed evolu-tionary changes for orthologous genes with conserved DNA methylation patterns. In summary, this study establishes a large resource of vertebrate and invertebrate DNA methylomes, it showcases the power of reference-free epigenome analysis in species for which no reference genomes are available, and it contributes an epigenetic perspective to the study of vertebrate evolution.

More supplementary data and interactive data exploration  are available on the [**supplementary website**](https://medical-epigenomics.org/papers/DNAmeth500species/).  

This repository provides a set of scripts used to analyze the data. Below please find step-by-step instructions on how to reproduce the analysis in the paper. There are two categories of scripts: those that have to be run in HPC mode on each species in parallel, and those used for summarizing and visualizing results by generating summary tables and figure panels. The provided HPC scripts are written for a SLURM scheduling system and would have to be adjusted for your own computational environment and scheduler.

## Linking figures and scripts:
|Panel|Script|
| --- | --- |
|1A	|Overview, no script | 
|1B	|02.3_stats_detail.R|
|1C	|02.6_create_ITOL.R|
|1D	|02.3_stats_detail.R|
|1E	|02.5_stats_methylation.R|
|1F	|03.1_predict_average_meth_levels.R & 03.2_kmer_stability.R|
|1G	|04.1_phylogeny_reconstruction.R|
|1H	|03.3_phyloglm.R|
|1I	|02.7_speciesMeta.R|
|1J	|02.7_speciesMeta.R|
|1K	|02.7_speciesMeta.R|
|2A	|Overview, no script|
|2B	|05.2_summarize_within_sp_prediction.R|
|2C	|05.2_summarize_within_sp_prediction.R|
|2D	|05.2_summarize_within_sp_prediction.R|
|2E	|05.3_feature_weight_analysis.R|
|2F	|05.3_feature_weight_analysis.R|
|3A	|05.5_summary_cross_species.R|
|3B	|05.4_predict_crossSp.R|
|3C	|05.5_summary_cross_species.R|
|3D	|06.1_inverted.R|
|3E	|06.4_explore_inverted_final_WHH.R|
|3F	|06.4_explore_inverted_final_WHH.R|
|4A	|08.2_tissue_correlation.R|
|4B	|Overview, no script|
|4C	|07.5_ame_selex.R|
|4D	|07.8_followup_analysis.ipynb & cytoscape| 
|4E |07.8_followup_analysis.ipynb & cytoscape|
|5A|validation/01.5.3_crossMapping_plot.R| 
|5B|validation/01.5.3_crossMapping_plot.R| 
|5C|validation/01.5.3_crossMapping_plot.R| 
|5D|validation/01.5.3_crossMapping_plot.R| 
|5E|validation/01.5.3_crossMapping_plot.R| 

## Reproducing the project step-by-step:

### Setup.
For computational analysis, we use R version 3.6.1 with all the corresponding packages. The SVM module requires an independent block with R version 4.1.
We provide a full list of packages installed in each of the subsets with their versions in the envs/ folder as csv files. For your convinience, we also provide the mirrors of conda environments as yaml files, generated with ```conda env export```. Please note that not all packages were installed within conda, so additionally you would still have to install missing R packages from the csv tables if needed for the corresponding part of the analysis. 
 
The project folder contains the following folders:
-src 
-meta
-raw_data (can be downloaded from SRA SRP357738)
-results_pipeline
-results_analysis

The *results_pipeline* folder should contain the per-species output of the RefFreeDMA. Most relevant files (deduced reference genome and methylation profiles) are shared and could be downloaded from the GEO repository under accession number GSE195869 or from the [supplementary website](https://medical-epigenomics.org/papers/DNAmeth500species/).

The *results_analysis* folder contains the output of most of the scripts below. We also have provided the most relevant summary tables for download on the [supplementary website](https://medical-epigenomics.org/papers/DNAmeth500species/).

00.0_init.R is the script, to initialize all further scripts, it will standardize settings, make annotations available, and set paths to data and scripts. A few absolute paths need to be adjusted and two environmental variables are required (they have to be visible by both Rscript and bash): 

```
CODEBASE (location of this github repository) 
LOGDIR (location for pipeline log and tmp files) 
```

### 00: Preparing to run the pipeline on raw data.

For data processing you would have to use [RefFreeDMA](https://github.com/jklughammer/RefFreeDMA) and [bisulfiteBlast](https://github.com/jklughammer/bisulfiteBlast).
The 00.1 script will prepare the run folders and shell scripts containing the commands to run RefFeeDMA. For setting up RefFreeDMA please follow the instructions in the [corresponding repository](https://github.com/jklughammer/bisulfiteBlast).

```bash
Rscript 00.1_prepareRefFreeDMA.R 
```
### 01: Creating a comprehensive sample annotation table.

Following scripts will generate some additional data metrices, that are used in most subsequent analysis.
Here we create the big summary annotation table and other initial stats:
```bash
Rscript 01.1_collectStats.R
Rscript 01.2_taxonomicAnnotation.R. #generates the stats annot table
Rscript 01.3_quality_stratification.R
```

Calcualating the frequency of CGs for subsequent analysis:
```
Rscript 01.4_dedRef_CpG_count.R  ##takes a while
```
Calculating the frequency of kmers, using the MEME tool:
!!Importantly, the path in 01.5.sbatch has to be adjusted

```bash
sbatch sbatch/01.5.1.sbatch ##prefilteing
sbatch sbatch/01.5.2.sbatch
Rscript 01.5.3_kmer_count_summary.R  ##summary for further use
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
Coverage statistics and sequence features:
```bash
Rscript 01.9_coverage_analysis.R
Rscript 01.10_seq_features.R
```
### 02: Main statistics.

The following scripts create summary statistics and generate summary plots, giving us an overview of the data:  
```bash
Rscript 02.1_stats_uc.R  
Rscript 02.2_stats_general.R  
Rscript 02.3_stats_detail.R  
Rscript 02.4_stats_location.R  
Rscript 02.5_stats_methylation.R 
Rscript 02.7_speicesMeta.R
```

Saving the supplementary tables and the stats for the paper:
```bash
Rscript 02.6_save_annot.R
```


### 03: Explaining mean DNA methylation.
Step-by-step exploration of features, that explain mean DNA methylation level per species:

```bash
Rscript 03.1_predict_average_meth_levels.R ##creating a table of features & prediction values for each group (takes a couple of hours on a large dataset)
Rscript 03.2_kmer3_stability.R ##stepwise selection on 3mers
Rscript 03.3_phyloglm.R ##phylogenetic models
```
### 04: Reconstructing phylogeny.
Using kmer frequencies to create a hierarchical clustering of species, based on the sequence composition of the consensus genomes:
```bash
Rscript 04.1_phylogeny_reconstruction.R
```
### 05: Predicting DNA methylation from sequence.
With-in species DNA methylation prediction:
```bash
sbatch sbatch/05.1.sbatch ##repeated on k 1:10, k = 3, k = 4, k = 7
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

### 06: Investigating the inverted prediction phenotype.
Identifying the species, that we classify as "inverted":
```bash
Rscript 06.1_inverted.R
```
Confirm the observed effect with additional experiments on subsets of data (by tissue, individual, and random sub-sampling):
```bash
sbatch sbatch/06.2.1.sbatch <candidate species>
sbatch sbatch/06.2.2.sbatch <candidate species>
sbatch sbatch/06.2.3.sbatch <candidate species>
Rscript 06.2_verify_inv_summary.R
```
Analysing the frequencies of kmers/repeats specifically in highly/lowly methylated sequences:
!!first two python scripts require additional specification of the data_dir path

```bash
python 06.3.1_to_fasta_for_kmer_count.py
sbatch 06.3.2.sbatch #runs on all the species, where classes contain inverted
RScript 06.3_freq_kmers.R
```
And summarizing all the knowledge, oriented on the WHH species
```
Rscript 06.4_explore_inverted_final_WHH.R
```

### 07: TFBS analysis across species and tissues. 

As a preparation, the JASPAR matrix of all the TF profiles has to be downloaded: http://jaspar.genereg.net/downloads/
We used the non-redundant vertebrate TFBS set (additionally download the JASPAR annotation table)
Calculate and report frequencies of each TFBS in consensus reference of each species:

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
For integrating the methylation preference for TF binding, we used the data from [Yin et. al, 2017](https://pubmed.ncbi.nlm.nih.gov/28473536/).

```bash
Rscript 07.5_ame_selex.R
Rscript 07.6_GO_enrichment.R
Rscript 07.7_GRN.R
```
Interactively: 07.8_followup_analysis.ipynb


### 08: Investigating sources of within-species DNA methylation variability.

Multidimensional scaling analysis for selected species:
```bash
Rscript 08.1_tissue_MDS.R
```
Comparison of DNA methylation variability due to tissue and individual:
```bash
Rscript 08.2_tissue_correlation.R
Rscript 08.3_tissue_distance.R
```
Create correlation plots between liver and heart DNA methylattion levels for selected species:
```bash
08.4_tissue_diffMeth.R
```

## Validation on publicly available reference genomes.
### 01: CrossMapping RefFreeDMA analysis to reference genomes.

determining which genomes to process, cleaning them up and concatting:

```bash
Rscript 01.1_find_mapping_match.R ## finding pairs of genomes to match
Rscript 01.1.2_axolotl.R ## rearrange and remove small chromosomes of the genome of the axolotl
sbatch 01.2_concat_genome.sbatch ## running genome concat
```
The download procedure is described downstream in 02.1_insilico_digest_prepare_genomes_and_tracks.ipynb.
Generating the config files and running the cross-mapping:
```bash
Rscript 01.3_save_configs.R
```
afterwards run the job-subbmitting bash script same way as for original RefFreeDMA.

Summarise the results of crossMapping:
```bash
Rscript 01.4_crossMapping_stats.R
sbatch sbatch/01.5_cross-mapping_analysis.sbatch
Rscript 01.5.2_crossMappingAnalysis_summary.R
Rscript 01.5.3_crossMapping_plot.R
```

### 02: Insilico RRBS simulation.

In interactive jupyer notebook run the preparation notebook for genomes:
02.1_insilico_digest_prepare_genomes.R

```bash
sbatch 02.2_insilicoDigest.sbatch
Rscript 02.3_insilicoDigest_summary.R
```
### 03: WGBS analysis.

WGBS data download is described at 03.1_download_WGBS.ipynb:
```bash
Rscript 03.2.process_WGBS_data.R
Rscript 03.2.2_compare_mean_meth.R

sbatch 03.3_WGBS_fragment_CpGmeth.sbatch ##preparing the data
```

##### SVM-modeling in WGBS.

```bash
sbatch 03.4.sbatch
sbatch 03.5.sbatch
sbatch 03.6.sbatch

Rscript 03.4.1_WGBS_prediction_summary.R
Rscript 03.5_summary.R
Rscript 03.6_summary.R
```
