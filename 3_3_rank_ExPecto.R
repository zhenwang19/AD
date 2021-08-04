##Risk variants ranked by ExPecto Score
library(as.color)
library(plyr)
library(gplots)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
##############################unique risk region,bed

if(dir.exists("./3_varRank/3_ExPecto_rank")){
    setwd("./3_varRank/3_ExPecto_rank")
}else{
    dir.create("./3_varRank/3_ExPecto_rank")
    setwd("./3_varRank/3_ExPecto_rank")
}
###########################################################################################################################


ex.max.gwas<-read.csv("../../2_varFunAnno/4_varExPecto/3_AD_ExPecto_maxScore_cell.gwas.txt",sep="\t",head=T)
ex.max<-read.table("../../2_varFunAnno/4_varExPecto/3_AD_ExPecto_maxScore_cell.gwas_1000G.txt",sep="\t",head=T)

###########################################################################################################################
####PLOT threshould
threshold.gwas<-log10(1.5)
threshold.1000G<-0.3

vline_gwas<-geom_hline(aes(yintercept =threshold.gwas),linetype="dashed",color="red", size=0.2)
vline_1000G<-geom_hline(aes(yintercept =threshold.1000G),linetype="dashed",color="red", size=0.2)

##top ten effect size
top.all.1000G.10<-NULL
top.neuro.1000G.10<-NULL
top.all.gwas.10<-NULL
top.neuro.gwas.10<-NULL

top.all.1000G.10<-sort(ex.max$max_all,decreasing=TRUE)[10]
top.neuro.1000G.10<-sort(ex.max$max_neuro,decreasing=TRUE)[10]
ex.all.1000G.text<-geom_text_repel(data =subset(ex.max, max_all>= top.all.1000G.10),aes(x=gene,y=max_all,label = paste(gene,max_all_cell,sep="_")),colour = "black",size =2,box.padding = unit(0.25, "lines"),point.padding = unit(0.3, "lines"))
ex.neuro.1000G.text<-geom_text_repel(data =subset(ex.max, max_neuro >= top.neuro.1000G.10),aes(x=gene,y=max_neuro,label = paste(gene,max_neuro_cell,varID,sep="_")),colour = "black",size = 2,box.padding = unit(0.25, "lines"),point.padding = unit(0.3, "lines"))

top.all.gwas.10<-sort(ex.max.gwas$max_all,decreasing=TRUE)[10]
top.neuro.gwas.10<-sort(ex.max.gwas$max_neuro,decreasing=TRUE)[10]
ex.all.gwas.text<-geom_text_repel(data =subset(ex.max.gwas, max_all>= top.all.gwas.10),aes(x=gene,y=max_all,label = paste(gene,max_all_cell,sep="_")),colour = "black",size =2,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))
ex.neuro.gwas.text<-geom_text_repel(data =subset(ex.max.gwas, max_neuro >= top.neuro.gwas.10),aes(x=gene,y=max_neuro,label = paste(gene,max_neuro_cell,sep="_")),colour = "black",size = 2,box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines"))

p.all.1000G<-ggplot(data=ex.max) +geom_jitter(aes(x=gene,y=max_all,col=gene),alpha=0.4)+theme_classic()+labs(title="",x="Genes(93/all_cell_line/1000G)", y = "Effect size")+vline_1000G+theme(axis.text.x = element_text(angle = 45, vjust = 0.5,size= 4),legend.position="none")+ex.all.1000G.text
p.neuro.1000G<-ggplot(data=ex.max) +geom_jitter(aes(x=gene,y=max_neuro,col=gene),alpha=0.4)+theme_classic()+labs(title="",x="Genes(93/neuron_cell_line/1000G)", y = "Effect size")+vline_1000G+theme(axis.text.x = element_text(angle = 45, vjust = 0.5,size= 4),legend.position="none")+ex.neuro.1000G.text
p.all.gwas<-ggplot(data=ex.max.gwas) +geom_jitter(aes(x=gene,y=max_all,col=gene),alpha=0.4)+theme_classic()+labs(title="",x="Genes(88/all_cell_line/GWAS)", y = "Effect size")+vline_gwas+theme(axis.text.x = element_text(angle = 45, vjust = 0.5,size= 4),legend.position="none")+ex.all.gwas.text
p.neuro.gwas<-ggplot(data=ex.max.gwas) +geom_jitter(aes(x=gene,y=max_neuro,col=gene),alpha=0.4)+theme_classic()+labs(title="",x="Genes(88/neuron_cell_line/GWAS)", y = "Effect size")+vline_gwas+theme(axis.text.x = element_text(angle = 45, vjust = 0.5,size= 4),legend.position="none")+ex.neuro.gwas.text

pdf("ExPecto_effect_size_AD_all_neruo_distribution.pdf", width = 15, height = 7.5)
grid.arrange(p.all.1000G,p.all.gwas,p.neuro.1000G,p.neuro.gwas,nrow =2)
dev.off()

###2, filter the variants that out of the gwas or 1000G
ex.var.gwas.all.filter<-subset(ex.max.gwas,max_all >= threshold.gwas)
ex.var.gwas.neuro.filter<-subset(ex.max.gwas,max_neuro >= threshold.gwas)
ex.var.1000G.all.filter<-subset(ex.max,max_all >= threshold.1000G)
ex.var.1000G.neuro.filter<-subset(ex.max,max_neuro >= threshold.1000G)

write.table(ex.var.gwas.all.filter,sep="\t",col.names=T,row.names=F,file="1_ExPecto.high_under_riskRegion_gwas_all_filter.txt",quote=F)
write.table(ex.var.gwas.neuro.filter,sep="\t",col.names=T,row.names=F,file="1_ExPecto.high_under_riskRegion_gwas_neuro_filter.txt",quote=F)
write.table(ex.var.1000G.neuro.filter,sep="\t",col.names=T,row.names=F,file="1_ExPecto.high_under_riskRegion_1000G_neuro_filter.txt",quote=F)
write.table(ex.var.1000G.all.filter,sep="\t",col.names=T,row.names=F,file="1_ExPecto.high_under_riskRegion_1000G_all_filter.txt",quote=F)

########################################################################################summary
###output the casual, gene and casual
ex.var.casual.gwas.all<-ex.var.gwas.all.filter
ex.var.casual.gwas.all$varIdExPecto<-paste(ex.var.casual.gwas.all$varID,ex.var.casual.gwas.all$max_all,sep="_")
ex.var.casual.gwas.all.out<-ex.var.casual.gwas.all %>%  group_by(gene) %>% summarise(casualVarsGwasAll = paste(varIdExPecto, collapse="; "))

ex.var.casual.gwas.neuro<-ex.var.gwas.neuro.filter
ex.var.casual.gwas.neuro$varIdExPecto<-paste(ex.var.casual.gwas.neuro$varID,ex.var.casual.gwas.neuro$max_neuro,sep="_")
ex.var.casual.gwas.neuro.out<-ex.var.casual.gwas.neuro %>%  group_by(gene) %>% summarise(casualVarsGwasNeuro = paste(varIdExPecto, collapse="; "))

ex.var.casual.1000G.all<-ex.var.1000G.all.filter
ex.var.casual.1000G.all$varIdExPecto<-paste(ex.var.casual.1000G.all$varID,ex.var.casual.1000G.all$max_all,sep="_")
ex.var.casual.1000G.all.out<-ex.var.casual.1000G.all %>%  group_by(gene) %>% summarise(casualVars1000GAll = paste(varIdExPecto, collapse="; "))

ex.var.casual.1000G.neuro<-ex.var.1000G.neuro.filter
ex.var.casual.1000G.neuro$varIdExPecto<-paste(ex.var.casual.1000G.neuro$varID,ex.var.casual.1000G.neuro$max_all,sep="_")
ex.var.casual.1000G.neuro.out<-ex.var.casual.1000G.neuro %>%  group_by(gene) %>% summarise(casualVars1000GNeuro = paste(varIdExPecto, collapse="; "))

#num, enh gene num overall
num.gene.list<-c("Term","gene_num","var_num")
num.gene.list<-rbind(num.gene.list,c("gene-var-gwas",length(unique(ex.max.gwas$gene)),length(unique(ex.max.gwas$varID))))
num.gene.list<-rbind(num.gene.list,c("gene-var-1000G",length(unique(ex.max$gene)),length(unique(ex.max$varID))))

num.gene.list<-rbind(num.gene.list,c("gene-var-gwas-casual-all",length(unique(ex.var.gwas.all.filter$gene)),length(unique(ex.var.gwas.all.filter$varID))))
num.gene.list<-rbind(num.gene.list,c("gene-var-gwas-casual-neuro",length(unique(ex.var.gwas.neuro.filter$gene)),length(unique(ex.var.gwas.neuro.filter$varID))))

num.gene.list<-rbind(num.gene.list,c("gene-var-1000G-casual-all",length(unique(ex.var.1000G.all.filter$gene)),length(unique(ex.var.1000G.all.filter$varID))))
num.gene.list<-rbind(num.gene.list,c("gene-var-1000G-casual-neuro",length(unique(ex.var.1000G.neuro.filter$gene)),length(unique(ex.var.1000G.neuro.filter$varID))))
write.table(num.gene.list,sep="\t",col.names=F,row.names=F,file="2_Summary_gene_ExPecto_var_number_overall.txt",quote=F)

#####num, var number for each gene

var.gene.gwas.uniq<-unique(ex.max.gwas[,c("gene","varID")])
num.var.gene.gwas<-table(var.gene.gwas.uniq$gene)
var.gene.1000G.uniq<-unique(ex.max[,c("gene","varID")])
num.var.gene.1000G<-table(var.gene.1000G.uniq$gene)

##MAX ALL
var.gene.gwas.filter.uniq<-unique(ex.var.gwas.all.filter[,c("gene","varID")])
num.var.gene.gwas.filter<-table(var.gene.gwas.filter.uniq$gene)
var.gene.1000G.filter.uniq<-unique(ex.var.1000G.all.filter[,c("gene","varID")])
num.var.gene.1000G.filter<-table(var.gene.1000G.filter.uniq$gene)

##MAX NEURON
var.gene.gwas.neuro.filter.uniq<-unique(ex.var.gwas.neuro.filter[,c("gene","varID")])
num.var.gene.gwas.neuro.filter<-table(var.gene.gwas.neuro.filter.uniq$gene)
var.gene.1000G.neuro.filter.uniq<-unique(ex.var.1000G.neuro.filter[,c("gene","varID")])
num.var.gene.1000G.neuro.filter<-table(var.gene.1000G.neuro.filter.uniq$gene)

####
num.enh.var<-merge(num.var.gene.gwas,num.var.gene.gwas.filter,by="Var1",all=T)
num.enh.var<-merge(num.enh.var,num.var.gene.gwas.neuro.filter,by="Var1",all=T)

num.enh.var<-merge(num.enh.var,num.var.gene.1000G,by="Var1",all=T)
num.enh.var<-merge(num.enh.var,num.var.gene.1000G.filter,by="Var1",all=T)
num.enh.var<-merge(num.enh.var,num.var.gene.1000G.neuro.filter,by="Var1",all=T)

colnames(num.enh.var)<-c("gene","gwas_var_num","gwas_var_all_filter_num","gwas_var_neuro_filter_num","1000G_var_num","1000G_var_all_filter_num","1000G_var_neuro_filter_num")

num.enh.var<-merge(num.enh.var,ex.var.casual.gwas.all.out,by="gene",all=T)
num.enh.var<-merge(num.enh.var,ex.var.casual.gwas.neuro.out,by="gene",all=T)
num.enh.var<-merge(num.enh.var,ex.var.casual.1000G.all.out,by="gene",all=T)
num.enh.var<-merge(num.enh.var,ex.var.casual.1000G.neuro.out,by="gene",all=T)

write.table(num.enh.var,sep="\t",col.names=T,row.names=F,file="3_Summary_gene_ExPecto_var_number.txt",quote=F)
setwd("../../")

