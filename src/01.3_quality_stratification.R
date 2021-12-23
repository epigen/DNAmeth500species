source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(pheatmap)

wd=file.path(analysis_dir,"01_basicStats/01.3_quality_stratification")
dir.create(wd)
setwd(wd)


quality_measures=c("mapping_efficiency","coveredCpGs","conversionRate","k1_unmeth" ,"k3_meth","others", "fragments_uncovered_perc","cont_rat","max_lowest_common" )

for (quality_measure in quality_measures){
  print(quality_measure)
  #converted samples
  if(quality_measure %in% c("mapping_efficiency","coveredCpGs","conversionRate","k3_meth","max_lowest_common" )){
    stats_annot[conversion_type=="converted",paste0(quality_measure,"_qual_tier"):=as.numeric(as.character(cut(get(quality_measure),breaks=quantile(get(quality_measure),c(0,0.05,0.5,0.95,1),na.rm=TRUE),labels=c(4,3,2,1),include.lowest=TRUE))),]
    }else {
    stats_annot[conversion_type=="converted",paste0(quality_measure,"_qual_tier"):=as.numeric(as.character(cut(get(quality_measure),breaks=quantile(get(quality_measure),c(0,0.05,0.5,0.95,1),na.rm=TRUE),labels=c(1,2,3,4),include.lowest=TRUE))),] 
  }
  #unconverted samples
  if(quality_measure %in% c("mapping_efficiency","coveredCpGs","k1_unmeth","k3_meth","max_lowest_common" )){
    stats_annot[conversion_type=="unconverted",paste0(quality_measure,"_qual_tier"):=as.numeric(as.character(cut(get(quality_measure),breaks=quantile(get(quality_measure),c(0,0.05,0.5,0.95,1),na.rm=TRUE)+c(0,0.000000001,0,-0.000000001,0),labels=c(4,3,2,1),include.lowest=TRUE))),]
  }else {
    stats_annot[conversion_type=="unconverted",paste0(quality_measure,"_qual_tier"):=as.numeric(as.character(cut(get(quality_measure),breaks=quantile(get(quality_measure),c(0,0.05,0.5,0.95,1),na.rm=TRUE)+c(0,0.000000001,0,-0.000000001,0),labels=c(1,2,3,4),include.lowest=TRUE))),] 
  }
}

stats_annot[,max_qual_tier:=max(unlist(mget(paste0(quality_measures,"_qual_tier"))),na.rm=TRUE),by=1:nrow(stats_annot)]
stats_annot[,mean_qual_tier:=mean(unlist(mget(paste0(quality_measures,"_qual_tier"))),na.rm=TRUE),by=1:nrow(stats_annot)]


cor_mat_conv=cor(stats_annot[conversion_type=="converted",mget(paste0(quality_measures,"_qual_tier")),],use="complete.obs")
cor_mat_unconv=cor(stats_annot[conversion_type=="unconverted",mget(paste0(quality_measures,"_qual_tier")),],use="complete.obs")

pdf("qual_tier_correlation.pdf",height=6,width=6.5)
pheatmap(cor_mat_conv,main="Correlation matrix converted")
pheatmap(cor_mat_unconv,main="Correlation matrix unconverted")
dev.off()

#check potentially discordant species
pdf("lowest_common_count_distrib.pdf", height=3,width=5.5)
ggplot(stats_annot,aes(x=blast_count1,fill=max_lowest_common>=26))+geom_density(alpha=0.5)+xlim(c(0,50))+geom_vline(xintercept=8,lty=2)
dev.off()

pdf("lowest_common_check.pdf", height=20,width=15)
ggplot(stats_annot[scientific_name%in%scientific_name[max_lowest_common<26&blast_count1>8]],aes(x=Sample_Name,y=max_lowest_common,fill=blast_count1))+geom_bar(stat="identity")+geom_hline(yintercept=10,lty=2)+rotate_labels()+facet_wrap(~abbreviation_sp+scientific_name,scale="free_x")+ggtitle("blast_count > 8")

ggplot(stats_annot[scientific_name%in%scientific_name[max_lowest_common<26&blast_count1<=8]],aes(x=Sample_Name,y=max_lowest_common,fill=blast_count1))+geom_bar(stat="identity")+geom_hline(yintercept=10,lty=2)+rotate_labels()+facet_wrap(~abbreviation_sp+scientific_name,scale="free_x")+ggtitle("blast_count <= 8")
dev.off()


#ABU_2_L  viral infection
#AFB_1_L_uc possible mixup (bird --> fish)
#*AN all samples affected (snake --> sheep, Python third)
#*BF all samples affected (Frog --> fish, Xenopus third)
#BL_1 all samples affected (aeromonas infection)
#BL_2_H lizard --> sheep, Anolis third
#*BP_1_H  snake --> sheep, Python third
#*BTF all samples affected (frog --> fish), origin from pathology vienna --> probably ok
#CAT_2  all samples affected (cat --> chicken), possibly mixup
#CD_2_S Crocodile --> mouse, possibly mixup
#CHD_1_H  contamination with e coli
#CTL_1_L  Lizard --> kangaroo, possibly mixup
#*CW_1_A_uc probably due to missing reference
#*ECC_1 sea star --> mussle, possibly due to contamination with mussle (food?)
#FPF_1_G_uc fish --> bird, possibly mixup
#*LSK_1_G ray --> carp, rey third
#NCL_1_H & NCL_1_H_uc lizard --> dog, possibly mixup
#*NCL_2_LU lizard --> turtle, possibly due to missing reference
#*PA nearly all samples affected (lizard --> sheep), probably due to missing reference
#PLF all samples affected (pseudomonas infection)
#PO_2_S klebsiella infection
#*PS all except 2 affected (lizard --> human/fish), probably due to missing reference
#RED_1_L  deer --> duck, probably mixup
#RI_2_H & RI_2_L  lizard -->human/turtle, possibly sample mixup
#RS_1_H snake --> mouse, probably sample mixup/food?
#*RSS all somewhat affected (sea squirt --> sea star, sea squirt 6th), possibly due to missing reference
#*SGS all somewhat affected (sea squirt --> clam, sea squirt later), possibly due to missing reference
#SME_1_G eel --> shark, probably mixup
#VIS_2_H & VIS_1_L  viscacha --> chicken/pig, probably mixup
#WA_3_L wallaby --> salmon, probably sample mixup
#WBD_2_F vibrio infection
#WG_1_H e coli infection
#*WIS multiple affected rey --> zebrafish, rey third
#YS_1_K aeromonas infection
#*ACT all affected probably due to missing reference. NOTE: Repl. 1 and 2 seem to be technical replicates (not biologicla)
#*ASC all affecte, probably due to missing reference
#*BD_2_L lizard --> sheep, lizard third
#*BF all affected,frog thisrd, missing reference
#*BLS all affected, missing reference
#*BRS all affected, missing reference
#*BST_1_TF_uc, sea star third
#*CCS all affected
#*CNR all affected
#*COL_2_L lizard third
#*CSD several affected, missing reference
#*CT all affected, missing reference
#*ESB_1_H_uc, mayby bacterial infection, also missing reference
#*FBS_1, missing reference
#*GB, missing reference
#*GE_1_L, sheep, Gekko later
#*GF, missing reference
#*GKS, missing reference
#*HT, all affected, toad --> carp, xenopus later, possibly missing reference
#*ISS, all affected, missing reference
#*LIS_1_H, squid third
#LOW_1 possibly something else (no good hits)
#*LOW_2 owl --> chicken (has eagle hits later)
#*MF nearly all affected, missing reference
#*MOL_1_H lizard --> sheep, lizard later
#NCL_1_H lizard --> dog, possibly sample mixup
#NCL_1_H_uc lizard --> dog, possibly sample mixup
#*NSS all affected, possibly missing reference
#*OG_1_H gekko -->sheep, gekko later
#*OG_1_L_uc gekko -->sheep, gekko later
#*OSC_2_TF cucumaria third
#*PA lizard -->sheep, lizard later
#*PCA_2_H chameleon --> sheep, cameleon third
#PR_2_L possibly something different
#*PS most affected, maybe missing reference
#*PSS probably missing reference
#*RAY probably missing reference
#RBF_2_G possibly something different
#*RIN_1_H_uc lizard --> crocodile, lizard third
#*RLF pssibly missing reference or sample mix: lantern fish --> medusa (low target complexity)
#*RSS possibly missing reference
#*SBS all affected, missing reference
#*SCC missing reference
#*SEC missing reference
#*SEP missing reference
#*SES missing reference
#*SGO_1_L lizard --> sheep, lizard later
#*SGS individuals might all be slightly different species
#*STL_1_H lizard --> sheep, lizard later
#*STT_2_P missing reference
#*TRF missing reference
#*VC_1_L missing reference
#*VC_2_L missing reference
#*WSC_2_DT missing reference
#*WSC_3_DT missing reference
#*XEN missing reference, lizard --> sheep, lizard later
#*YFS_1_M missing reference 


#check potentially problematic samples
susp_df <- stats_annot[scientific_name%in%scientific_name[max_lowest_common<26&blast_count1>8],list(Sample_Name,scientific_name, blast_species2, blast_species1,  ncbi_name, ncbi_order, ncbi_class,   ncbi_group, ncbi_name_blastS1, ncbi_name_blastS2, lowest_common_blastS1, lowest_common_blastS2, max_lowest_common, cont_rat, blast_count1, blast_count2)]

susp_df[max_lowest_common<26&blast_count1>8, mismatch:="true",] ## marking problematic sample

## save the suspisious samples into a table and 
write.csv(susp_df, "problematic_samples.tsv")
print("Stop here and manually explore the table!")
#stop()

## by hand checking and using the annitation
susp_df_curated <- fread("problematic_samples_annotated.csv")
## uploading the new version
#problematic_samples_old=c("ABU_2_L$","AFB_1_L_uc","BL_1","BL_2_H$","CAT_2","CD_2_S$","CHD_1_H$","CTL_1_L$","FPF_1_G_uc","NCL_1_H$","NCL_1_H_uc","PLF","PO_2_S$","RED_1_L$","RI_2_H$", "RI_2_L$","RS_1_H$","SME_1_G$","VIS_2_H$","VIS_1_L$","WA_3_L$","WBD_2_F$","WG_1_H$","YS_1_K$","LOW_1","NCL_1_H$","NCL_1_H_uc","PR_2_L$","RBF_2_G$")

problematic_samples <- susp_df_curated[summary==1]$Sample_Name
##added problematic, based on prev runs:
problematic_samples <- c(problematic_samples, 
                         "LOW_1_H", "LOW_1_L", "RBF_2_G",  "SME_1_G")

#LOW_1* has an extremly low blust count and was mapped as problematic last time
#SME_1_G quite high mapping rate and was also mapped to a distant with a higher rate last time
#RFB_2_G  low blust count and was mapped as problematic last time

#check
sort(stats_annot[grep(paste0(problematic_samples,collapse="|^"),Sample_Name)]$Sample_Name)
sort(problematic_samples)

#add column
stats_annot[,species_check:=ifelse(Sample_Name%in%grep(paste0(problematic_samples,collapse="|^"),Sample_Name,value=TRUE),"fail","pass"),]

my_wt(stats_annot,"stats_annot.tsv")
save(stats_annot,file="stats_annot.RData")
