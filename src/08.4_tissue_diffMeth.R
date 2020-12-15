source(file.path(Sys.getenv("CODEBASE"),"DNAmeth500species/src/00.0_init.R"))

wd=file.path(analysis_dir,"04_tissue")
dir.create(wd)
setwd(wd)


all_species=c("BU","WIF","NBI","ABU","JAK","DAS","NBI","EH","YS")

for (species in all_species){
print(species)
diff_file=paste0(data_dir,"/results_pipeline/",species,"/toSelf_filtered_0.08mm_final_concat/diffMeth_cpg/",species,"_diff_meth.tsv")
if(!file.exists(diff_file)){
  print("Not available")
  next
}  
  
diffMeth=fread(diff_file)

heart_file=paste0(data_dir,"/results_pipeline/",species,"/motifAnalysis/ded_top500_cov2_Heart.fa")
liver_file=paste0(data_dir,"/results_pipeline/",species,"/motifAnalysis/ded_top500_cov2_Liver.fa")
if(!file.exists(heart_file)|!file.exists(liver_file)){
  print("Not available")
  next
}

heart=readLines(heart_file)
liver=readLines(liver_file)

heart_sel=sub(">","",grep(">",heart,value=TRUE))
liver_sel=sub(">","",grep(">",liver,value=TRUE))


cov_trsh_l=c(2)
cov_trsh_u=c(Inf)
i=1

diffMeth_filt=diffMeth[meth.cov_mean_g1>cov_trsh_l[i]&meth.cov_mean_g2>cov_trsh_l[i]]
top_diff_data=diffMeth_filt[dedRef_ID%in%c(heart_sel,liver_sel)]


ggp_frag=ggplot(diffMeth_filt,aes(x=meth.meth_mean_g1,y=meth.meth_mean_g2))+ 
  stat_binhex(bins=30,col="white") +
  geom_point(data=top_diff_data, aes(x=meth.meth_mean_g1,y=meth.meth_mean_g2),position=position_jitter(width=3,height=3),alpha=0.4,col="green",size=1,shape=21) + 
  ggtitle(,label=paste0("r = ",round(cor(diffMeth_filt[,meth.meth_mean_g1],diffMeth_filt[,meth.meth_mean_g2]),3),"\n cov>= ",cov_trsh_l[i], "\n cov <=",cov_trsh_u[i] ,"\nN = ",dim(diffMeth_filt)[1])[1])+scale_fill_gradient(limits=c(0,nrow(diffMeth_filt)/200),high="blue",low="white",na.value="blue")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+xlab(paste0("% methylation Heart"))+ylab(paste0("% methylation Liver"))+xlim(-3,103)+ylim(-3,103)+coord_fixed()


pdf(paste0(species,"_diffMeth.pdf"),height=4,width=4.5)
print(ggp_frag)
dev.off()

}
