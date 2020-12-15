source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))
library(tableHTML)

wd=file.path(analysis_dir,"01_basicStats")
dir.create(wd,recursive=TRUE)
setwd(wd)

stats_annot=fread("all_stats.tsv")

#find samples for unconverted sequencing
stats_optim=stats_annot[!grepl("_uc",Sample_Name)][order(`fragments_uncovered_perc`)]
stats_optim[,rank:=c(1:length(`fragments_uncovered_perc`)),by="species"]
stats_optim=stats_optim[order(rank,`pre-BC_CT`)]
#stat_missing_uc[Sample_Name==sel]

stats_optim[,pool_group:=as.integer(NA),]
#now make pooling suggestions
group_size=20
group=1
row=1
max_perc_cov=30
max_ct_diff=3
#exclude samples from Experiments before AK15
exclude_exp=c("AK02","AK09","AK10","AK11","AK12","AK13","AK14")
exclude_samples=c("ATK_1_M","MAS_1_GO","OG_1_H","ACO_1_G","MT_1_B", "CRF_1_BL", "RES_1_H", "WSD_1_G", "MBI_2_L","CTL_2_LU","LYN_2_H")

stats_optim[`Experiment ID`%in%exclude_exp,pool_group:=-2,]
stats_optim[`Sample_Name`%in%exclude_samples,pool_group:=-2,]

#when samples have already been selected as UC remove all samples from that species from the list before creating new groups
done_sp=unique(stats_optim[!is.na(unconverted)]$abbreviation_sp)
stats_optim[abbreviation_sp%in%done_sp,pool_group:=-3,]


while (nrow(stats_optim[is.na(pool_group)])>0){
  
  if (!is.na(stats_optim[row]$pool_group)){
    row=row+1
    print(paste0(nrow(stats_optim[is.na(pool_group)])," - ", group," - ",row))
    next
  }
  
  while (nrow(stats_optim[pool_group==group])<group_size){
    taken_adapters=stats_optim[pool_group==group]$Adapter
   # ct_mean=mean(stats_optim[pool_group==group]$`pre-BC_CT`)
   ct_min=min(c(stats_optim[pool_group==group]$`pre-BC_CT`,Inf))
   ct_max=max(c(stats_optim[pool_group==group]$`pre-BC_CT`,-Inf))
   if (ct_min==Inf){ct_min=stats_optim[row]$`pre-BC_CT`}
   if (ct_max==-Inf){ct_max=stats_optim[row]$`pre-BC_CT`}
   
    if (!stats_optim[row]$Adapter%in%taken_adapters&is.na(stats_optim[row]$pool_group)&(stats_optim[row]$fragments_uncovered_perc<=max_perc_cov|stats_optim[row]$rank==min(stats_optim[abbreviation_sp==stats_optim[row]$abbreviation_sp&is.na(pool_group)]$rank))&(abs(stats_optim[row]$`pre-BC_CT`-ct_min)<=max_ct_diff)&(abs(stats_optim[row]$`pre-BC_CT`-ct_max)<=max_ct_diff)){
      stats_optim[row,pool_group:=group]
      stats_optim[species==stats_optim[row]$species&is.na(pool_group),pool_group:=-1,]
    }
    if (row < nrow(stats_optim)-1){
      row=row+1
      print(paste0(nrow(stats_optim[is.na(pool_group)])," - ", group," - ",row))
    }else{break}
  }
  group=group+1
  row=1
  
  exclude_sp=unique(stats_optim[`Experiment ID`%in%exclude_exp,]$abbreviation_sp)
  if (nrow(stats_optim[!abbreviation_sp%in%exclude_sp&is.na(pool_group)])==0){
    stats_optim[abbreviation_sp%in%exclude_sp,pool_group:=-2]
  }
  print(paste0(nrow(stats_optim[is.na(pool_group)])," - ", group," - ",row))
}

table(stats_optim$pool_group) #For Tasmanian devil the first sample of TD_3_Tm was selected for unconverted, but later the second sample was kept (automated selection of the second sample if there are multiple) --> fix by excluding second sample --> rerun analysis
stats_optim[pool_group==2,list(rank, pool_group,fragments_uncovered_perc,`pre-BC_CT`,Sample_Name,Adapter),]


stats_optim_30=copy(stats_optim)

stats_optim_30_sel=stats_optim[pool_group>0,list(rank, pool_group,fragments_uncovered_perc,`pre-BC_CT`,Sample_Name,species,Adapter,`Experiment ID`,group_count=.N),by=pool_group]

write.table(stats_optim_30_sel,"stats_optim_30_sel_new4.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#only needed for sharing
#write.table(stats_optim_30_sel,"/data/groups/lab_bock/public_html/jklughammer/compEpi_auto/stats_optim_30_sel.tsv",sep="\t",quote=FALSE,row.names=FALSE)
#cat(paste0("<a href='stats_optim_30_sel.tsv'>download table</a>\n",tableHTML(stats_optim_30_sel,collapse="separate",theme="rshiny-blue",spacing = "4px")),file="/data/groups/lab_bock/public_html/jklughammer/compEpi_auto/unconv_sel.html")


#find species/samples for which uc is still missing

stats_optim_30_sel_red=stats_optim_30_sel[pool_group<7]

#use after selecting new pool groups
stat_missing_uc=stats_annot[!grepl("_uc",Sample_Name)][!abbreviation_sp%in%abbreviation_sp[!is.na(unconverted)]&!abbreviation_sp%in%stats_optim_30_sel_red$species]

#use in final round, when in principle all species should be covered 
#should be empty
stat_missing_uc=stats_annot[!grepl("_uc",Sample_Name)][!abbreviation_sp%in%abbreviation_sp[!is.na(unconverted)]]
#check if _ucs have been sequenced
stats_annot[is.na(lane_uc)&!is.na(unconverted)]


stat_missing_uc[,rank_others:=rank(others,ties.method="min"),by=abbreviation_sp]
stat_missing_uc[,rank_uncov:=rank(fragments_uncovered_perc,ties.method="min"),by=abbreviation_sp]
stat_missing_uc[,rank_comb:=pmax(rank_others,rank_uncov),]
stat_missing_uc[,sel:=Sample_Name[which.min(rank_comb)],by=abbreviation_sp]


stat_missing_uc_sel=stat_missing_uc[Sample_Name==sel]
stat_missing_uc_sel[,pool_group:=0]

stats_optim_30_sel_red_mean_ct=stats_optim_30_sel_red[,list(mean_ct=median(`pre-BC_CT`),diff_max_min=max(`pre-BC_CT`)-min(`pre-BC_CT`)),by=pool_group]

for (pg in stats_optim_30_sel_red_mean_ct[order(diff_max_min)]$pool_group){
  pool=stats_optim_30_sel_red[pool_group==pg]
  gc_diff=18-unique(pool$group_count)
  mean_ct=mean(pool$`pre-BC_CT`)
  stat_missing_uc_sel[,cd_diff:=abs(`pre-BC_CT`-mean_ct),]
  stat_missing_uc_sel=stat_missing_uc_sel[order(cd_diff)]
  sel_samples=stat_missing_uc_sel[pool_group==0][1:gc_diff]$Sample_Name
  stat_missing_uc_sel[Sample_Name%in%sel_samples, pool_group:=pg]
  
}

stat_missing_uc_sel_re=stat_missing_uc_sel[,c("pool_group","fragments_uncovered_perc","pre-BC_CT","Sample_Name","abbreviation_sp","Experiment ID","Box","Well_y","Well_x","DNA (ng/µl)","Volume (µl)"),with=FALSE]
setnames(stat_missing_uc_sel_re,"abbreviation_sp","species")

stats_optim_30_sel_combined=rbindlist(list(stats_optim_30_sel_red,stat_missing_uc_sel_re),use.names=TRUE,fill=TRUE)
stats_optim_30_sel_combined=stats_optim_30_sel_combined[order(pool_group)]
stats_optim_30_sel_combined[,group_count_new:=.N,by="pool_group"]

write.table(stats_optim_30_sel_combined,"stats_optim_30_sel_new4_combined.tsv",sep="\t",quote=FALSE,row.names=FALSE )
