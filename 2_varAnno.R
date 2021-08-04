###var annotation score from LINGSIGHT, PRIMATEAI, EXPECTO

library(as.color)
library(plyr)
library(gplots)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
##############################unique risk region,bed
#setwd("./2_varFunAnno/1_var")
if(!dir.exists("./2_varFunAnno")){
    dir.create("./2_varFunAnno")
}

if(dir.exists("./2_varFunAnno/1_var")){
    setwd("./2_varFunAnno/1_var")
}else{
    dir.create("./2_varFunAnno/1_var")
    setwd("./2_varFunAnno/1_var")
}

##1 retrieve the variants under risk region from gwas or 1000G
cmd1<-c("bedtools intersect -b ../../1_summary_region_gene/0_riskRegion.bed.txt -a /Users/zhenwang/Documents/Database/23_AD_summaryVriants_2018/ALZ_NatGenet_2018_ref_bed.txt -wa > 1_AD_highRiskRegion.gwas.var.txt")
cmd2<-c("bedtools intersect -b ../../1_summary_region_gene/0_riskRegion.bed.txt -a /Users/zhenwang/Documents/Database/4_1k_genome_project/1_vcf/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -wa > 1_AD_highRiskRegion.1000G.var.txt")

system(cmd1)
system(cmd2)

##########################################################################################################################################################################
##2 merge gwas and 1000G
G1000<-read.table("1_AD_highRiskRegion.1000G.var.txt",head=F,sep="\t")
gwas<-read.table("1_AD_highRiskRegion.gwas.var.txt",head=F,sep="\t")

G1000$var_id<-paste(G1000$V1,G1000$V2,sep="_")
gwas$var_id<-paste(gwas$V1,gwas$V2,sep="_")

G1000_gwas<-G1000[G1000$var_id %in% gwas$var_id,]
G1000_gwas$varSource<-c("GWAS_1000G")

G1000_uniq<-G1000[!(G1000$var_id %in% gwas$var_id),]
G1000_uniq$varSource<-c("1000G")
gwas_uniq<-gwas[!(gwas$var_id %in% G1000$var_id),]
gwas_uniq$varSource<-c("GWAS")

G1000_gwas_bed<-G1000_gwas[,c("V1","V2","V2","varSource","V3")]
G1000_uniq.bed<-G1000_uniq[,c("V1","V2","V2","varSource","V3")]
gwas_uniq_bed<-gwas_uniq[,c("V1","V2","V3","varSource","V4")]

colnames(G1000_gwas_bed)<-c("var_chr","var_start","var_end","var_source","var_id")
colnames(G1000_uniq.bed)<-c("var_chr","var_start","var_end","var_source","var_id")
colnames(gwas_uniq_bed)<-c("var_chr","var_start","var_end","var_source","var_id")

var.bed<-rbind(gwas_uniq_bed,G1000_gwas_bed,G1000_uniq.bed)
var.uniq.bed<-var.bed[!duplicated(var.bed[,c(1,2,3)]),]

var.uniq.chr.bed<-var.uniq.bed
var.uniq.chr.bed$var_chr<-paste("chr",var.uniq.chr.bed$var_chr,sep="")

write.table(var.uniq.bed,file="2_AD_highRiskRegion.gwas_1000G_merge.bed.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(var.uniq.chr.bed,file="2_AD_highRiskRegion.gwas_1000G_merge.chr.bed.txt",quote=F,row.names=F,col.names=F,sep="\t")
##########################################################################################################################################################################
###3, annot the Linsight to variants
#setwd("../2_varLinsight")
if(dir.exists("../2_varLinsight")){
    setwd("../2_varLinsight")
}else{
    dir.create("../2_varLinsight")
    setwd("../2_varLinsight")
}
cmd3<-c("bedtools intersect -b /Users/zhenwang/Documents/Project5_AD/1_dataColleciton/15_linsight/LINSIGHT.bed -a ../1_var/2_AD_highRiskRegion.gwas_1000G_merge.chr.bed.txt -wa -wb > 1_Linsight_AD_highRiskRegion.gwas_1000G.txt")
system(cmd3)

var.chr.bed<-read.table("../1_var/2_AD_highRiskRegion.gwas_1000G_merge.chr.bed.txt", head=F,sep="\t")
colnames(var.chr.bed)<-c("var_chr","var_pos","var_end","var_source","varID")

var.linsight<-read.table("1_Linsight_AD_highRiskRegion.gwas_1000G.txt",sep="\t",head=F)
var.linsight.new<-subset(var.linsight, V2 != V7)
var.linsight.score<-var.linsight.new[,c("V1","V2","V9")]
colnames(var.linsight.score)<-c("var_chr","var_pos","linsight")

var.bed.linsight<-merge(var.chr.bed,var.linsight.score,by=c("var_chr","var_pos"),all.x=T)
var.bed.linsight$var_chr<-gsub("chr","",var.bed.linsight$var_chr)

write.table(var.bed.linsight,file="2_AD_var_all_linsight.bed.txt",sep="\t",quote=F,col.names=T,row.names=F)
write.table(var.bed.linsight,file="2_AD_var_all_linsight.noHead.bed.txt",sep="\t",quote=F,col.names=F,row.names=F)

##########################################################################################################################################################################
###4, annot the PrimateAI to variants
##create the file
if(dir.exists("../3_varPrimateAI")){
    setwd("../3_varPrimateAI")
}else{
    dir.create("../3_varPrimateAI")
    setwd("../3_varPrimateAI")
}
getwd()

ucscGene<-read.table("/Users/zhenwang/Documents/Database/16_PrimateAI/ucsc_geneID_2_geneSymbol.txt",head=T,sep="\t")
PrimateAIgene<-read.table("/Users/zhenwang/Documents/Database/16_PrimateAI/PrimateAI_scores_v0.2.tsv.geneName.txt",head=T,sep="\t")
PrimateAIgene.symbol<-merge(PrimateAIgene,ucscGene,by.x="UCSC_gene",by.y="kgID")
PrimateAIscore<-read.table(gzfile("/Users/zhenwang/Documents/Database/16_PrimateAI/PrimateAI_scores_v0.2.tsv.gz"),header=T)

high<-read.table("../../../3_summary_gene/3_1_ALS.risk.gene.high.score.list.txt",sep="\t",head=F)
high.ucsc<-merge(high,PrimateAIgene.symbol,by.x="V1",by.y="geneSymbol")
colnames(high.ucsc)<-c("geneSymbol","UCSC_gene")
PrimateAIscore.high<-merge(high.ucsc,PrimateAIscore,by.x="UCSC_gene")
write.table(PrimateAIscore.high,file="1_PrimateAIscore.high.txt",row.names=F,col.names=T,sep="\t",quote=F)

high.without.primateAI<-high[!(high$V1 %in% high.ucsc$geneSymbol),]
write.table(high.without.primateAI,file="1_PrimateAI_nonIncluding_highGenes.txt",row.names=F,col.names=F,sep="\t",quote=F)

##########################################################################################################################################################################
###5, annot the ExPecto to variants
##create the file
if(dir.exists("../4_varExPecto")){
    setwd("../4_varExPecto")
}else{
    dir.create("../4_varExPecto")
    setwd("../4_varExPecto")
}
getwd()

##start from the risk gene and all the promoter effectsize of each cell line
#5-1 merge all the effect size of cells into one file

effect<-read.table("/Users/zhenwang/Documents/Database/14_ExPecto/tissueName/effects_coor_AD.txt",head=T)
cell.file<-read.table("/Users/zhenwang/Documents/Database/14_ExPecto/tissueName/tissueName.txt",sep="\t",head=F)

for (i in 1:dim(cell.file)[1]){
cell<-cell.file[i,1]
print(i)
print(cell)
cell.name<-gsub("effects_pergene_mat_","",cell)
cell.name<-gsub(".txt","",cell.name)
effect.cell<-read.table(paste("/Users/zhenwang/Documents/Database/14_ExPecto/tissueName/1_AD_cellLine/AD_",cell,sep=""),head=F)
effect.cell<-data.frame(as.vector(t(effect.cell)))
colnames(effect.cell)<-cell.name
effect<-cbind(effect,effect.cell)
}
write.table(effect,file="1_AD_ExPecto_allCell_mat.txt",quote=F,col.names=T,row.names=F,sep="\t")

####5-2 match only the variants from gwas or 1000G to reduce the size of the file
effect.size<-read.table("1_AD_ExPecto_allCell_mat.txt",sep="\t",head=T)
effect.size$var_id<-paste(effect.size$chr,effect.size$pos,effect.size$ref,effect.size$alt,sep="_")
effect.varid<-data.frame(var_id=unique(effect.size$var_id))

G1000<-read.table("../1_var/1_AD_highRiskRegion.1000G.var.txt",head=F,sep="\t")
gwas<-read.table("../1_var/1_AD_highRiskRegion.gwas.var.txt",head=F,sep="\t")

G1000$var_id<-paste(G1000$V1,G1000$V2,G1000$V4,G1000$V5,sep="_")
gwas$var_id<-paste(gwas$V1,gwas$V2,gwas$V5,gwas$V6,sep="_")

G1000_gwas<-data.frame(var_id=effect.varid[(effect.varid$var_id %in% gwas$var_id) & (effect.varid$var_id %in% G1000$var_id),])
G1000_gwas<-merge(G1000_gwas,gwas[,c("var_id","V4")],by="var_id")
G1000_gwas<-merge(G1000_gwas,G1000[,c("var_id","V3")],by="var_id")
G1000_gwas$varSource<-c("GWAS_1000G")
colnames(G1000_gwas)<-c("varID","snpID_gwas","snpID_1000G","varSource")

G1000_uniq<-data.frame(var_id=effect.varid[(!effect.varid$var_id %in% gwas$var_id) & (effect.varid$var_id %in% G1000$var_id),])
G1000_uniq$snpID_gwas<-c("-")
G1000_uniq<-merge(G1000_uniq,G1000[,c("var_id","V3")],by="var_id")
G1000_uniq$varSource<-c("1000G")
colnames(G1000_uniq)<-c("varID","snpID_gwas","snpID_1000G","varSource")

gwas_uniq<-data.frame(var_id=effect.varid[(effect.varid$var_id %in% gwas$var_id) & (!effect.varid$var_id %in% G1000$var_id),])
gwas_uniq<-merge(gwas_uniq,gwas[,c("var_id","V4")],by="var_id")
gwas_uniq$snpID_1000G<-c("-")
gwas_uniq$varSource<-c("GWAS")
colnames(gwas_uniq)<-c("varID","snpID_gwas","snpID_1000G","varSource")

var.all<-rbind(G1000_gwas,gwas_uniq,G1000_uniq)
effect.gwas.1000G<-merge(var.all,effect.size,by.x="varID",by.y="var_id")
write.table(effect.gwas.1000G,file="2_AD_ExPecto_allCell_mat.gwas_1000G.txt",quote=F,col.names=T,row.names=F,sep="\t")

#####5-3 get the highiest score of each variants
nerous.system<-c("astrocyte","bipolar.spindle.neuron","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Germinal_Matrix","Brain_Hippocampus","Brain_Hippocampus_Middle","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c1","Brain_Substantia_nigra","diencephalon","Fetal_Brain_Female","frontal.cortex","Nerve_Tibial","neural.cell","neural.progenitor.cell","occipital.lobe","parietal.lobe","Pituitary","spinal.cord","temporal.lobe")
ex<-read.table("2_AD_ExPecto_allCell_mat.gwas_1000G.txt",head=T,sep="\t")
###1 return the abs max value and cell line
ex$max_all<-apply(ex[,11:228],1,function(x)max(abs(x)))
ex$max_all_cell<-apply(ex[,11:228],1,function(x)names(ex[,11:228])[which.max(abs(x))])

ex$max_neuro<-apply(ex[,nerous.system],1,function(x)max(abs(x)))
ex$max_neuro_cell<-apply(ex[,nerous.system],1,function(x)names(ex[,nerous.system])[which.max(abs(x))])

ex.max<-ex[,c("varID","snpID_gwas","snpID_1000G","varSource","chr","pos","ref","alt","gene","max_all","max_all_cell","max_neuro","max_neuro_cell")]
ex.max.gwas<-ex.max[!(ex.max$varSource=="1000G"),]

write.table(ex.max,file="3_AD_ExPecto_maxScore_cell.gwas_1000G.txt",sep="\t",quote=F,col.names=T,row.names=F)
write.table(ex.max.gwas,file="3_AD_ExPecto_maxScore_cell.gwas.txt",sep="\t",quote=F,col.names=T,row.names=F)

setwd("../../")

