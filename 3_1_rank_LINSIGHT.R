
library(as.color)
library(plyr)
library(gplots)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
##############################unique risk region,bed

if(!dir.exists("./3_varRank")){
    dir.create("./3_varRank")
}


if(dir.exists("./3_varRank/1_Linsight_rank")){
    setwd("./3_varRank/1_Linsight_rank")
}else{
    dir.create("./3_varRank/1_Linsight_rank")
    setwd("./3_varRank/1_Linsight_rank")
}


####1) retrive enhancer of each risk gene across the whole genome
high.gene<-read.table("../../../3_summary_gene/3_1_ALS.risk.gene.high.score.list.txt",sep="\t",head=F)
colnames(high.gene)<-c("gene")
enh2gene.all<-read.table("/Users/zhenwang/Documents/Project5_AD/1_dataColleciton/17_enhTarget/enhId2gene.txt",sep="\t",head=F)
colnames(enh2gene.all)<-c("enhID","gene")
gene.enhID<-merge(high.gene,enh2gene.all,by="gene")

enh<-read.table("/Users/zhenwang/Documents/Project5_AD/1_dataColleciton/17_enhTarget/enh.cell.target.merge.txt",sep="\t",head=F)
colnames(enh)<-c("enhID","enh_chr","enh_start","enh_end","enh_cell","enh_source","enh_target")

gene.enh.infor<-merge(gene.enhID,enh,by="enhID")
gene.enh.bed<-gene.enh.infor[,c("enh_chr","enh_start","enh_end","enhID","gene","enh_cell","enh_source")]
write.table(gene.enh.bed,sep="\t",col.names=F,row.names=F,file="1_alzGene_enhancer_infor_wholeGenome.txt",quote=F)

###bedtools, enhancer under risk region,
cmd1<-c("bedtools intersect -b ../../1_summary_region_gene/0_riskRegion.bed.txt -a 1_alzGene_enhancer_infor_wholeGenome.txt -wa -wb > 2_alzGeneEnh_under_riskRegion.txt ")
system(cmd1)

###match gene-enhancer-riskRegion, and variants with linsight score
gene.enh.risk<-read.table("2_alzGeneEnh_under_riskRegion.txt",sep="\t",head=F)
colnames(gene.enh.risk)<-c("enh_chr","enh_start","enh_end","enhID","gene","enh_cell","enh_source","risk_chr","risk_start","risk_end","risk_length","risk_id")
gene.enh.risk$geneRiskId<-paste(gene.enh.risk$gene,gene.enh.risk$risk_chr,gene.enh.risk$risk_start,gene.enh.risk$risk_end,sep="_")

risk.region.gene<-read.table("../../1_summary_region_gene/0_high.risk.gene.loci.tss_full.txt",head=T,sep="\t")
gene.enh.risk.match<-subset(gene.enh.risk,geneRiskId %in% risk.region.gene$geneRiskId)
length(unique(gene.enh.risk.match$gene))
#gene.enh.risk.unmatch<-subset(gene.enh.risk,!(geneRiskId %in% risk.region.gene$geneRiskId))
write.table(gene.enh.risk.match,sep="\t",col.names=F,row.names=F,file="3_alzGeneEnh_under_riskRegion_match.txt",quote=F)

cmd2<-c("bedtools intersect -b ../../2_varFunAnno/2_varLinsight/2_AD_var_all_linsight.noHead.bed.txt -a 3_alzGeneEnh_under_riskRegion_match.txt -wa -wb > 4_alzGeneEnh_under_riskRegion_var_linsight.txt ")
system(cmd2)
##########

##########################################################################################################################################################################
#Plot the score distribution, data filter
##########################################################################################################################################################################
###Summary the AD gene, enhancer num, variants num distribution,gene-enhancer, variants,
##gene enhancer connection under risk region
threshold<-0.9
gene.enh.risk.match<-read.table("3_alzGeneEnh_under_riskRegion_match.txt",head=F,sep="\t")
colnames(gene.enh.risk.match)<-c("enh_chr","enh_start","enh_end","enhID","gene","enh_cell","enh_source","risk_chr","risk_start","risk_end","risk_length","risk_id","geneRiskId")
gene.enh.risk.var.linsight<-read.table("4_alzGeneEnh_under_riskRegion_var_linsight.txt",sep="\t",head=F)
colnames(gene.enh.risk.var.linsight)<-c("enh_chr","enh_start","enh_end","enhID","gene","enh_cell","enh_source","risk_chr","risk_start","risk_end","risk_length","risk_id","geneRiskID","var_chr","var_pos","var_pos2","var_source","varID","linsight")

gene.enh.risk.var.linsight.gwas<-subset(gene.enh.risk.var.linsight, var_source !="1000G")
gene.enh.risk.var.linsight.1000G<-subset(gene.enh.risk.var.linsight, var_source !="GWAS")

gene.enh.risk.var.linsight.gwas.filter<-subset(gene.enh.risk.var.linsight.gwas,linsight>= threshold )
gene.enh.risk.var.linsight.1000G.filter<-subset(gene.enh.risk.var.linsight.1000G,linsight>=threshold)

write.table(gene.enh.risk.var.linsight.gwas.filter,sep="\t",col.names=T,row.names=F,file="5_gene.enh.risk.var.linsight.gwas.filter.txt",quote=F)
write.table(gene.enh.risk.var.linsight.1000G.filter,sep="\t",col.names=T,row.names=F,file="5_gene.enh.risk.var.linsight.1000G.filter.txt",quote=F)

###output the casual, gene and casual

gene.var.casual.gwas<-gene.enh.risk.var.linsight.gwas.filter
gene.var.casual.gwas$varIdLinsight<-paste(gene.var.casual.gwas$var_chr,gene.var.casual.gwas$var_pos,gene.var.casual.gwas$varID,gene.var.casual.gwas$linsight,sep="_")
gene.var.casual.gwas.out<-gene.var.casual.gwas %>%  group_by(gene) %>% summarise(casualVarsGwas = paste(varIdLinsight, collapse="; "))

gene.var.casual.1000G<-gene.enh.risk.var.linsight.1000G.filter
gene.var.casual.1000G$varIdLinsight<-paste(gene.var.casual.1000G$var_chr,gene.var.casual.1000G$var_pos,gene.var.casual.1000G$varID,gene.var.casual.1000G$linsight,sep="_")
gene.var.casual.1000G.out<-gene.var.casual.1000G %>%  group_by(gene) %>% summarise(casualVars1000G = paste(varIdLinsight, collapse="; "))

######################################
#num, enh gene num overall
num.gene.list<-c("Term","gene_num","enh_or_var_num")
num.gene.list<-rbind(num.gene.list,c("gene-enh",length(unique(gene.enh.risk.match$gene)),length(unique(gene.enh.risk.match$enhID))))

num.gene.list<-rbind(num.gene.list,c("gene-var-gwas",length(unique(gene.enh.risk.var.linsight.gwas$gene)),length(unique(gene.enh.risk.var.linsight.gwas$varID))))
num.gene.list<-rbind(num.gene.list,c("gene-var-1000G",length(unique(gene.enh.risk.var.linsight.1000G$gene)),length(unique(gene.enh.risk.var.linsight.1000G$varID))))

num.gene.list<-rbind(num.gene.list,c("gene-var-gwas-casual",length(unique(gene.enh.risk.var.linsight.gwas.filter$gene)),length(unique(gene.enh.risk.var.linsight.gwas.filter$varID))))
num.gene.list<-rbind(num.gene.list,c("gene-var-1000G-casual",length(unique(gene.enh.risk.var.linsight.1000G.filter$gene)),length(unique(gene.enh.risk.var.linsight.1000G.filter$varID))))
write.table(num.gene.list,sep="\t",col.names=F,row.names=F,file="6_Summary_gene_enh_var_number_overall.txt",quote=F)

#num, enh number for eahc gene
enh.gene.uniq<-unique(gene.enh.risk.match[,c("gene","enhID")])
num.enh.gene<-table(enh.gene.uniq$gene)

var.gene.gwas.uniq<-unique(gene.enh.risk.var.linsight.gwas[,c("gene","varID")])
num.var.gene.gwas<-table(var.gene.gwas.uniq$gene)
var.gene.1000G.uniq<-unique(gene.enh.risk.var.linsight.1000G[,c("gene","varID")])
num.var.gene.1000G<-table(var.gene.1000G.uniq$gene)

var.gene.gwas.filter.uniq<-unique(gene.enh.risk.var.linsight.gwas.filter[,c("gene","varID")])
num.var.gene.gwas.filter<-table(var.gene.gwas.filter.uniq$gene)
var.gene.1000G.filter.uniq<-unique(gene.enh.risk.var.linsight.1000G.filter[,c("gene","varID")])
num.var.gene.1000G.filter<-table(var.gene.1000G.filter.uniq$gene)

num.enh.var<-merge(num.enh.gene,num.var.gene.gwas,by="Var1",all=T)
num.enh.var<-merge(num.enh.var,num.var.gene.gwas.filter,by="Var1",all=T)
num.enh.var<-merge(num.enh.var,num.var.gene.1000G,by="Var1",all=T)
num.enh.var<-merge(num.enh.var,num.var.gene.1000G.filter,by="Var1",all=T)
colnames(num.enh.var)<-c("gene","enh_num","gwas_var_num","gwas_var_filter_num","1000G_var_num","1000G_var_filter_num")

num.enh.var<-merge(num.enh.var,gene.var.casual.gwas.out,by="gene",all=T)
num.enh.var<-merge(num.enh.var,gene.var.casual.1000G.out,by="gene",all=T)
write.table(num.enh.var,sep="\t",col.names=T,row.names=F,file="7_Summary_gene_enh_var_number.txt",quote=F)

setwd("../../")
getwd()







