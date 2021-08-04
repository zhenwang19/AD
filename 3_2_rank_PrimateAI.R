
library(as.color)
library(plyr)
library(gplots)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
##############################unique risk region,bed

if(dir.exists("./3_varRank/2_PrimateAI_rank")){
    setwd("./3_varRank/2_PrimateAI_rank")
}else{
    dir.create("./3_varRank/2_PrimateAI_rank")
    setwd("./3_varRank/2_PrimateAI_rank")
}


####1) convert to bed format for gene and their all primateAI score
primateAI.all<-read.table("../../2_varFunAnno/3_varPrimateAI/1_PrimateAIscore.high.txt",head=T,sep="\t")
primateAI.all.bed<-primateAI.all
primateAI.all.bed$var_end<-primateAI.all.bed$pos
primateAI.all.bed$chr<-gsub("chr","",primateAI.all.bed$chr)
primateAI.all.bed<-primateAI.all.bed[,c("chr","pos","var_end","geneSymbol","UCSC_gene","ref","alt","refAA","altAA","strand_1pos_0neg","trinucleotide_context","ExAC_coverage","primateDL_score")]
write.table(primateAI.all.bed,file="1_PrimateAIscore.high_all.bed.txt",sep="\t",quote=F,col.names=F,row.names=F)

cmd1<-c("bedtools intersect -b ../../1_summary_region_gene/0_riskRegion.bed.txt -a 1_PrimateAIscore.high_all.bed.txt -wa -wb > 2_PrimateAIscore.high_under_riskRegion.txt ")
system(cmd1)

primateAI.var.risk<-read.table("2_PrimateAIscore.high_under_riskRegion.txt",sep="\t",head=F)
colnames(primateAI.var.risk)<-c(colnames(primateAI.all.bed),"risk_chr","risk_start","risk_end","riskID","risk_length")
primateAI.var.risk$geneRiskId<-paste(primateAI.var.risk$geneSymbol,primateAI.var.risk$riskID,sep="_")

risk.region.gene<-read.table("../../1_summary_region_gene/0_high.risk.gene.loci.tss_full.txt",head=T,sep="\t")
primateAI.var.risk.match<-subset(primateAI.var.risk,geneRiskId %in% risk.region.gene$geneRiskId)
primateAI.var.risk.match$varID<-paste(primateAI.var.risk.match$chr,primateAI.var.risk.match$pos,primateAI.var.risk.match$ref,primateAI.var.risk.match$alt,sep="_")
#write.table(primateAI.var.risk.match,sep="\t",col.names=T,row.names=F,file="3_PrimateAIscore.high_under_riskRegion_match.txt",quote=F)

###2, filter the variants that out of the gwas or 1000G
threshold<-0.7
var.gwas<-read.table("../../2_varFunAnno/1_var/1_AD_highRiskRegion.gwas.var.txt")
var.1000G<-read.table("../../2_varFunAnno/1_var/1_AD_highRiskRegion.1000G.var.txt")

var.gwas$varID<-paste(var.gwas$V1,var.gwas$V2,var.gwas$V5,var.gwas$V6,sep="_")
var.1000G$varID<-paste(var.1000G$V1,var.1000G$V2,var.1000G$V4,var.1000G$V5,sep="_")

primateAI.var.gwas<-primateAI.var.risk.match[primateAI.var.risk.match$varID %in% var.gwas$varID,]
primateAI.var.1000G<-primateAI.var.risk.match[primateAI.var.risk.match$varID %in% var.1000G$varID,]

primateAI.var.gwas.filter<-subset(primateAI.var.gwas,primateDL_score >= threshold)
primateAI.var.1000G.filter<-subset(primateAI.var.1000G,primateDL_score >= threshold)

write.table(primateAI.var.gwas,sep="\t",col.names=T,row.names=F,file="3_PrimateAIscore.high_under_riskRegion_match_gwas.txt",quote=F)
write.table(primateAI.var.1000G,sep="\t",col.names=T,row.names=F,file="3_PrimateAIscore.high_under_riskRegion_match_1000G.txt",quote=F)
write.table(primateAI.var.gwas.filter,sep="\t",col.names=T,row.names=F,file="4_PrimateAIscore.high_under_riskRegion_match_gwas_filter.txt",quote=F)
write.table(primateAI.var.1000G.filter,sep="\t",col.names=T,row.names=F,file="4_PrimateAIscore.high_under_riskRegion_match_1000G_filter.txt",quote=F)

###output the casual, gene and casual
primateAI.var.casual.gwas<-primateAI.var.gwas.filter
primateAI.var.casual.gwas$varIdPrimateAI<-paste(primateAI.var.casual.gwas$varID,primateAI.var.casual.gwas$primateDL_score,sep="_")
primateAI.var.casual.gwas.out<-primateAI.var.casual.gwas %>%  group_by(geneSymbol) %>% summarise(casualVarsGwas = paste(varIdPrimateAI, collapse="; "))

primateAI.var.casual.1000G<-primateAI.var.1000G.filter
primateAI.var.casual.1000G$varIdPrimateAI<-paste(primateAI.var.casual.1000G$varID,primateAI.var.casual.1000G$primateDL_score,sep="_")
primateAI.var.casual.1000G.out<-primateAI.var.casual.1000G %>%  group_by(geneSymbol) %>% summarise(casualVars1000G = paste(varIdPrimateAI, collapse="; "))

#num, enh gene num overall
num.gene.list<-c("Term","gene_num","var_num")
num.gene.list<-rbind(num.gene.list,c("gene-var-gwas",length(unique(primateAI.var.gwas$geneSymbol)),length(unique(primateAI.var.gwas$varID))))
num.gene.list<-rbind(num.gene.list,c("gene-var-1000G",length(unique(primateAI.var.1000G$geneSymbol)),length(unique(primateAI.var.1000G$varID))))
num.gene.list<-rbind(num.gene.list,c("gene-var-gwas-casual",length(unique(primateAI.var.gwas.filter$geneSymbol)),length(unique(primateAI.var.gwas.filter$varID))))
num.gene.list<-rbind(num.gene.list,c("gene-var-1000G-casual",length(unique(primateAI.var.1000G.filter$geneSymbol)),length(unique(primateAI.var.1000G.filter$varID))))
write.table(num.gene.list,sep="\t",col.names=F,row.names=F,file="5_Summary_gene_PrimateAI_var_number_overall.txt",quote=F)

#####num, var number for each gene

var.gene.gwas.uniq<-unique(primateAI.var.gwas[,c("geneSymbol","varID")])
num.var.gene.gwas<-table(var.gene.gwas.uniq$geneSymbol)
var.gene.1000G.uniq<-unique(primateAI.var.1000G[,c("geneSymbol","varID")])
num.var.gene.1000G<-table(var.gene.1000G.uniq$geneSymbol)

var.gene.gwas.filter.uniq<-unique(primateAI.var.gwas.filter[,c("geneSymbol","varID")])
num.var.gene.gwas.filter<-table(var.gene.gwas.filter.uniq$geneSymbol)
var.gene.1000G.filter.uniq<-unique(primateAI.var.1000G.filter[,c("geneSymbol","varID")])
num.var.gene.1000G.filter<-table(var.gene.1000G.filter.uniq$geneSymbol)

num.enh.var<-merge(num.var.gene.gwas,num.var.gene.gwas.filter,by="Var1",all=T)
num.enh.var<-merge(num.enh.var,num.var.gene.1000G,by="Var1",all=T)
num.enh.var<-merge(num.enh.var,num.var.gene.1000G.filter,by="Var1",all=T)
colnames(num.enh.var)<-c("geneSymbol","gwas_var_num","gwas_var_filter_num","1000G_var_num","1000G_var_filter_num")

num.enh.var<-merge(num.enh.var,primateAI.var.casual.gwas.out,by="geneSymbol",all=T)
num.enh.var<-merge(num.enh.var,primateAI.var.casual.1000G.out,by="geneSymbol",all=T)
write.table(num.enh.var,sep="\t",col.names=T,row.names=F,file="7_Summary_gene_enh_var_number.txt",quote=F)

setwd("../../")

