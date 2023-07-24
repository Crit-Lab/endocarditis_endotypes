#########CLUSTERING#########

# Carga de datos ----

library(DESeq2)
library(pheatmap)
library(tidyverse)
library(genefilter)
library(ensembldb)
library(rjson)
library(readr)
library(RColorBrewer)
library(AnnotationHub)
library(tximport)
library(clusterProfiler)
library(umap)
library(factoextra)
library(org.Hs.eg.db)

#Colors:
library(critcolors) #To install, use: devtools::install_github("cecilomar6/critcolors")
cluster1<-"#486D87"
cluster2<-"#FFB838"
clusters_colors<-c(cluster1, cluster2)
sig_colors<-critcolors(100)
exp_colors<-colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

#Clinical data
samples<-read_csv("anonimized_data.csv")

#Load transcriptomes

ah<-AnnotationHub()
edb<-ah[["AH73986"]]

txi<-readRDS("transcripts_d1_github.rds") #Local file
#Removing case without RNAseq (all patiens must have its corresponding fastq file)
samples<-samples[samples$nhc %in% colnames(txi$abundance),]

#Sorting clinical data
samples<-samples[match(colnames(txi$abundance), samples$nhc),]

#Count matrix normalized 
mat<-txi$abundance

#log2 
mat <- log2(mat +1)

#Select highly variable features

mat.filter <- varFilter(mat, var.cutoff = 0.95, filterByQuantile = TRUE)
#With varFilter(), remove features that exhibit little variation between samples.
#Such non-specific filtering can be advantageous for downstream data analysis.
#Prevents occurrence of "Curse od dimensionality" --> The common theme of these problems is that when dimensionality increases,
#the volume of space increases so fast that available data becomes scarce.
#To get a reliable result, the amount of data needed often grows exponentially with dimensionality.
#In addition, the organization and search of data is often based on the detection of areas where objects form groups with similar properties;
#However, in high-dimensional data, all objects appear to be sparse and different in many respects, which prevents common data organization strategies from being efficient.

dim(mat)#34528
dim(mat.filter)#1727

genenames<-AnnotationDbi::select(edb, keys=rownames(mat.filter), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(mat.filter)<-genenames$SYMBOL 

#Distance matrices
distance <- dist(t(mat.filter), method = "euclidean")


# Creación de Clusters ----

clusters <- hclust(distance, method = "ward.D") 

#Hierarchical tree of the samples
#K=2
hc<-hcut(distance, k=2, hc.method="ward.D", hc_func = "hclust")
fviz_dend(hc, 
          k=2,
          show_labels = FALSE,
          rect=TRUE,
          rect_fill=TRUE,
          k_colors = rev(clusters_colors),
          main=NULL,
          cex=32,
          ggtheme=theme_classic(base_size = 24))+
  theme(aspect.ratio=1.618)
clustercut <- as.factor(cutree(clusters, 2)) #Uses a hclust(clusters) to create a tree with a given number of clusters
table(clustercut,useNA="always")
#clustercut
# 1    2 <NA> 
# 18   14    0 

samples$cluster<-clustercut[match(samples$nhc, names(clustercut))]

#UMAP ----

set.seed(789)
my.umap <- umap::umap(t(mat.filter))
dat <- as.data.frame(my.umap$layout)
dat$clust <- samples$cluster[match(rownames(dat), samples$nhc)]

ggplot(dat, aes(x=V1,y=V2,col = clust))+
  geom_point(size = 2.5)+
  scale_color_manual(values = clusters_colors)+
  theme_gray(base_size = 24)+
  theme(aspect.ratio = 1, legend.position = "none")+
  labs(x="UMAP1", y="UMAP2")


#Heatmap of genes used for clustering
sc.mat.filter <- t(apply(X = mat.filter,MARGIN = 1,FUN = function(y){
  (y - mean(y))/sd(y)}))
samples$colors<-clusters_colors[as.numeric(samples$cluster)]

df <- data.frame(samples$cluster)
colnames(df)[1] <- "Cluster"
rownames(df)<-colnames(mat)

pheatmap(sc.mat.filter,
         #fontsize = 4,
         cellwidth = 4,
         cellheight = 0.2,
         clustering_method="ward.D",
         cluster_cols=TRUE,
         cluster_rows = TRUE,
         show_rownames = F,
         show_colnames = F,
         annotation_colors=list(Cluster=c("1"=clusters_colors[1], "2"=clusters_colors[2])),
         #border_color = "#979a9a",
         annotation_col = df)

#Clinical data at admission
median_iqr<-function(x) {
  percentiles<-round(quantile(x, probs=c(0.25, 0.75), na.rm=TRUE),2)
  paste(round(median(x, na.rm=TRUE),2), " (", percentiles[1]," - ", percentiles[2],")", sep="")
}

samples$nl_ratio<-samples$neutrofilos_1/samples$linfocitos_1

quant<-c("edad", "imc", "charlson", "delay",
         "leucos_1", 
         "neutrofilos_1", "linfocitos_1", "nl_ratio", "monocitos_1",  "hb_1", 
         "plaq_1", 
         "cr_1", 
         "ast_1", 
         "alt_1", "ldh_1", 
         "pcr_1", 
         "tp_1", 
         "ttpa_1", 
         "qx_euroscoreii %")

quant_results<-data.frame(variable=character(), C1=character(), C2=character(), C3=character(), p_value=numeric())

for(i in quant) {
  values<-aggregate(as.numeric(samples[[i]])~samples$cluster, FUN=median_iqr)[,2]
  p_val<-round(wilcox.test(as.numeric(samples[[i]])~samples$cluster)$p.value,3)
  df.line<-data.frame(variable=i, CTP1=values[1], CTP2=values[2], p_value=p_val)
  quant_results<-rbind(quant_results, df.line)
  
}
#View(quant_results)

qual<-c("sexo", "clinica_fecha_global", "ei_tipo", "ei_Ao", "ei_Mi", "ei_tric", "ei_pulm", "ei_marcap", "ei_multiple",
        "ei_etiologia", "ei_qx", "fumador", "alcohol", "enfpulmonar", "coronaria", 
        "icc", "dm", "hta", "dislipemia", "evp", "acv", "tejconectivo",
        "inmunodeprimido", "neo", "neo_tipo", "insufrenal", "hepatopatia", "advp", 
        "charlson", "ei_origen", "ei_origen_tipo", "hc_germen", "cultivo_germen", 
        "comp_icc", "comp_nrl", "comp_abd_renal", "comp_abd_hep", "comp_insufrenal",
        "qx", "clinica_inicio_fiebremeg", "clinica_inicio_cardiaca", "clinica_inicio_nrl", 
        "clinica_inicio_pulmonar", "clinica_inicio_renal", "clinica_inicio_reumatica")

qual_results<-data.frame(variable=character(), 
                         values=character(), 
                         C1=numeric(), 
                         C2=numeric(),
                         p.value=character())

for(i in qual) {
  values<-table(samples[[i]], samples$cluster)
  p_val<-round(chisq.test(values)$p.value,3)
  df.chunk<-data.frame(variable=c(i, rep("", nrow(values)-1)),
                       values=rownames(values),
                       C1=as.numeric(values[,"1"]), 
                       C2=as.numeric(values[,"2"]),
                       p.value=c(p_val, rep("", nrow(values)-1)))
  qual_results<-rbind(qual_results, df.chunk)
  
}

#View(qual_results)




#Surgery
table(samples$qx, samples$cluster) #0: No indication; 1: Not done; 2: Scheduled; 3: Urgent; 4: Emergent

#Differential expression
dds_cluster<-DESeqDataSetFromTximport(txi, colData = samples, design = ~cluster) 

#Prefiltering
dim(dds_cluster)#34528
keep<-rowSums(counts(dds_cluster)>0)>=10
dds_cluster<-dds_cluster[keep,]
dim(dds_cluster)#16772

#DE-analysis
dds_cluster<-DESeq(dds_cluster)
genenames<-AnnotationDbi::select(edb, keys=rownames(dds_cluster), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(dds_cluster)<-genenames$SYMBOL 

#Dispersion-plot
plotDispEsts(dds_cluster)

#Results
res<-results(dds_cluster, alpha=0.05)
summary(res)
resOrdered_cluster <- res[order(res$log2FoldChange),]
sig.genes<-res[!is.na(res$padj) & res$padj<0.01,]
res$gene <- rownames(res)

genes<-as.data.frame(res)
genes$genename<-rownames(genes)


#GSEA between clusters
dds_gsea<-DESeqDataSetFromTximport(txi, colData = samples, design = ~cluster) 
keep<-rowSums(counts(dds_gsea)>0)>=10
dds_gsea<-dds_gsea[keep,]
dds_gsea<-DESeq(dds_gsea)
df<-results(dds_gsea, alpha = 0.05) 

original_gene_list<-df$stat
names(original_gene_list)<-rownames(df)

#omit NAs
gene_list<-na.omit(original_gene_list)

#sort gene list in decreasing order (required for clusterProfiler)
gene_list<-sort(gene_list, decreasing = T)

## Gene Set Enrichment ----

set.seed(789)
gse_dx<-gseGO(geneList = gene_list,
           ont = "BP",
           keyType = "ENSEMBL",
           pvalueCutoff = 0.01,
           minGSSize = 10,
           maxGSSize = 800,
           verbose = T,
           OrgDb = org.Hs.eg.db,
           pAdjustMethod = "BH")

GO<-as.data.frame(gse_dx)

library(DOSE)
library(scales)
library(munsell)

dim(gse_dx)

dotplot(gse_dx,
        showCategory= 76,
        split=".sign", font.size=9, label_format=100)+
  facet_grid(.~.sign)+
  theme(aspect.ratio = 2.5)

library(enrichplot)
treeplot(enrichplot::pairwise_termsim(gse_dx),
         nCluster=6,
         showCategory=50,
         font.size=5,
         label_format = 100,
         color="enrichmentScore")+
  scale_color_gradientn(name = "fold change",
                        colours = exp_colors,
                        limits= c(-1, 1), 
                        breaks=c(-1 , 0, 1) ) +
  theme(aspect.ratio = 1)

#Adding labels to volcano plot

go_labels_tcell<-gse_dx@result$ID[grepl("T cell", gse_dx@result$Description)]
go_labels_stat<-gse_dx@result$ID[grepl("STAT", gse_dx@result$Description)]

go_labels_genes_tcell<-unique(unlist(strsplit(gse_dx@result$core_enrichment[gse_dx@result$ID %in% go_labels_tcell], "/")))
go_labels_genes_stat<-unique(unlist(strsplit(gse_dx@result$core_enrichment[gse_dx@result$ID %in% go_labels_stat], "/")))
go_labels_genes_tcell<-bitr(go_labels_genes_tcell, fromType="ENSEMBL", toType="SYMBOL", OrgDb = "org.Hs.eg.db")
go_labels_genes_stat<-bitr(go_labels_genes_stat, fromType="ENSEMBL", toType="SYMBOL", OrgDb = "org.Hs.eg.db")
go_labels_genes_tcell<-unique(go_labels_genes_tcell$SYMBOL)
go_labels_genes_stat<-unique(go_labels_genes_stat$SYMBOL)

go_labels_genes_tcell<-go_labels_genes_tcell[go_labels_genes_tcell %in% rownames(sig.genes)]
go_labels_genes_stat<-go_labels_genes_stat[go_labels_genes_stat %in% rownames(sig.genes)]

genes$Significant <- ifelse(genes$padj <= 0.01 , "FDR < 0.01", "Not Sig")
genes$Significant[genes$Significant=="FDR < 0.01" & genes$genename %in% go_labels_genes_tcell]<-"FDR < 0.01, T cell pathway"
genes$Significant[genes$Significant=="FDR < 0.01" & genes$genename %in% go_labels_genes_stat]<-"FDR < 0.01, STAT pathway"


library(ggrepel)

ggplot(genes, aes(x = log2FoldChange, y = -log10(padj))) +
  #geom_point(aes(color = Significant), shape=1, alpha=(0.99-0.2*(genes$padj>0.01))) +
  geom_point(aes(color = Significant), 
             shape=1+15*((rownames(genes) %in% c(go_labels_genes_stat, go_labels_genes_tcell))), 
             alpha=(0.1+
                                                         0.24*(genes$padj>0.01)+
                                                         0.75*(rownames(genes) %in% c(go_labels_genes_stat, go_labels_genes_tcell)))) +
  geom_hline(yintercept = -log10(0.01), linetype=2)+
  theme_light(base_size = 24) + 
  theme(legend.position = "none")+
  geom_text_repel(
    data = genes[(genes$genename %in% go_labels_genes_stat) & genes$padj<0.01,],
    aes(label = genename),
    size = 4,
    point.padding = unit(2, "lines"),
    max.overlaps=50)+
  geom_text_repel(
    data = genes[(genes$genename %in% go_labels_genes_tcell) & genes$padj<0.0001,],
    aes(label = genename),
    size = 4,
    point.padding = unit(2, "lines"),
    max.overlaps=50)+
  xlim(-4, 4)+
  theme(aspect.ratio = 1/1.618)+
  scale_color_manual(values=c("#07044BFF", "#486D87", "#FFB838", "lightgrey"))

#Deconvolution
library(MetaIntegrator)
library(data.table)
library(dplyr)
library(ggpubr)

immunoStatesMatrix2<-immunoStatesMatrix[,!colnames(immunoStatesMatrix) %in% c("MAST_cell", "macrophage_m0", "macrophage_m1", "macrophage_m2")]
bulk<-counts(dds_cluster, normalized=T) #Counts matrix (with annotate's genes)
sum(rownames(bulk) %in% rownames(immunoStatesMatrix2)) #297 genes from our bulk is in immunoStatesMatrix

#Filtering bulk to leave only genes presents in immunoStatesMatrix
bulk<-bulk[rownames(bulk) %in% rownames(immunoStatesMatrix2),] 

#unique genes
bulk<-rowsum(bulk, rownames(bulk)) 
dim(bulk)#288 genes

#run immunoStates
outDT <-as.data.table(MetaIntegrator:::iSdeconvolution(immunoStatesMatrix2, bulk), keep.rownames = T)
#add new variables (compiling the different variables of each cell type)
outDT[,natural_killer_cell:=CD56bright_natural_killer_cell+CD56dim_natural_killer_cell]
outDT[,monocyte:=CD14_positive_monocyte+CD16_positive_monocyte]
outDT[,B_cell:=naive_B_cell+memory_B_cell+plasma_cell]
outDT[,T_cell:=CD8_positive_alpha_beta_T_cell+CD4_positive_alpha_beta_T_cell+gamma_delta_T_cell]
outDT[,granulocyte:=neutrophil+eosinophil+basophil]
outDT[,lymp:=B_cell+T_cell]

keep<-apply(outDT[,2:17],2, function(x) sum(x>0)>8) #Only cell populations present in more than 8 samples
cell_types<-colnames(outDT[,2:17])[keep]

outDT_scaled<-as.data.frame(outDT)
outDT_scaled<-outDT_scaled[,cell_types]
outDT_scaled<-as.data.frame(t(apply(outDT_scaled,1, FUN=function(x) x/sum(x))))
rowSums(outDT_scaled)

cell_types<-c("neutrophil", "eosinophil", "basophil",
              "CD14_positive_monocyte", "plasmacytoid_dendritic_cell", "CD56dim_natural_killer_cell", 
              "naive_B_cell", "memory_B_cell", "plasma_cell",
              "CD4_positive_alpha_beta_T_cell", "CD8_positive_alpha_beta_T_cell")

pcalc<-function(x) {
  wt<-wilcox.test(as.numeric(x)~as.factor(samples$cluster))
  return(round(wt$p.value, 3))
}

p_val<-apply(outDT_scaled[cell_types], 2, pcalc)
p_val<-paste("p=", p_val, sep="")
p_val[p_val=="p=0"]<-"p<0.001"
p_val<-as.factor(p_val)

cyto_plot<-function(celltype, titulo, subtitulo) {
  cytoquine<-as.numeric(outDT_scaled[[celltype]])
  genot<-as.factor(samples$cluster)
  ggplot(NULL, aes(y=cytoquine, x=genot, col=genot, fill=genot))+
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
    geom_jitter(width=0.15, height=0, stroke=0, size=3.5)+
    theme_grey(base_size = 12)+
    theme(aspect.ratio=1.618, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
          strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
    labs(x=NULL, y=NULL, title=titulo)+
    scale_x_discrete(labels=c("EE1", "EE2"))+
    scale_y_continuous(labels=scales::percent)+
    scale_color_manual(values=clusters_colors)+
    scale_fill_manual(values=clusters_colors)+
    facet_wrap(subtitulo)
}

titulos<-cell_types
titulos<-c("neutrophil", "eosinophil", "basophil",
           "CD14_positive_monocyte", "plasmacytoid_dendritic_cell", "CD56dim_natural_killer_cell", 
           "naive_B_cell", "memory_B_cell", "plasma_cell",
           "CD4_positive_alpha_beta_T_cell", "CD8_positive_alpha_beta_T_cell")

lista<-list()
for(i in 1:length(titulos)) {
  lista[[i]]<-cyto_plot(titulos[[i]], titulos[[i]], p_val[[i]])
  #ggsave(filename=paste(titulos[[i]],".pdf", sep=""), plot=lista[[i]])
}

library(cowplot)
plot_grid(plotlist=lista, align="hv", nrow=4)



#Effect of surgery: changes in gene expression and cell populations

samples2<-read.csv("anonimized_data2.csv")
table(samples2$cluster)/2
summary(samples2$delay)
aggregate(samples2$delay~samples2$cluster, FUN=summary)
wilcox.test(samples2$delay~samples2$cluster)

txi<-readRDS("transcripts_qx_github.rds")

#Sorting clinical data
samples2<-samples2[match(colnames(txi$abundance), samples2$samples_name),]
samples2$condition<-factor(samples2$condition, levels=c("pre_qx", "post_qx"))

samples_cluster1 <- samples2 %>% filter(cluster == 1)
samples_cluster1$condition
rownames(samples_cluster1) <- samples_cluster1$samples_name

txi_cluster1 <- txi
txi_cluster1$abundance <- txi_cluster1$abundance[,colnames(txi_cluster1$abundance) %in% samples_cluster1$samples_name]
txi_cluster1$counts <- txi_cluster1$counts[,colnames(txi_cluster1$counts) %in% samples_cluster1$samples_name]
txi_cluster1$length <- txi_cluster1$length[,colnames(txi_cluster1$length) %in% samples_cluster1$samples_name]

all(rownames(samples_cluster1) == colnames(txi_cluster1$abundance))

ddstxi_cluster1 <- DESeqDataSetFromTximport(txi_cluster1, colData = samples_cluster1, design = ~condition)
dim(ddstxi_cluster1)#34528

keep<-rowSums(counts(ddstxi_cluster1)>0)>=7

ddstxi_cluster1 <- ddstxi_cluster1[keep,]
dim(ddstxi_cluster1)

dds_cluster1 <- DESeq(ddstxi_cluster1)
res_clust1 <- results(dds_cluster1)
summary(res_clust1, alpha=0.05)

res_clust1 <- as.data.frame(res_clust1)
res_clust1$padj[is.na(res_clust1$padj)]<-1

original_gene_list <- res_clust1$stat
names(original_gene_list) <- rownames(res_clust1)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

set.seed(987)
gse_cluster1 <- gseGO(geneList = gene_list,
                      ont = "BP",
                      keyType = "ENSEMBL",
                      minGSSize = 10,
                      maxGSSize = 800,
                      pvalueCutoff = 1,
                      verbose = T,
                      OrgDb = org.Hs.eg.db,
                      pAdjustMethod = "BH")

treeplot(enrichplot::pairwise_termsim(gse_cluster1),
         nCluster=5,
         showCategory=gse_cluster1@result$Description[1:50],
         font.size=5,
         label_format = 100,
         color="enrichmentScore")+
  scale_color_gradientn(name = "fold change",
                        colours = exp_colors,
                        limits= c(-1, 1), 
                        breaks=c(-1 , 0, 1) )



samples_cluster2 <- samples2 %>% filter(cluster == 2)
samples_cluster2$condition
rownames(samples_cluster2) <- samples_cluster2$samples_name

txi_cluster2 <- txi
txi_cluster2$abundance <- txi_cluster2$abundance[,colnames(txi_cluster2$abundance) %in% samples_cluster2$samples_name]
txi_cluster2$counts <- txi_cluster2$counts[,colnames(txi_cluster2$counts) %in% samples_cluster2$samples_name]
txi_cluster2$length <- txi_cluster2$length[,colnames(txi_cluster2$length) %in% samples_cluster2$samples_name]

all(rownames(samples_cluster2) == colnames(txi_cluster2$abundance))

ddstxi_cluster2 <- DESeqDataSetFromTximport(txi_cluster2, colData = samples_cluster2, design = ~condition)
dim(ddstxi_cluster2)#34528

keep<-rowSums(counts(ddstxi_cluster2)>0)>=5

ddstxi_cluster2 <- ddstxi_cluster2[keep,]
dim(ddstxi_cluster2)

dds_cluster2 <- DESeq(ddstxi_cluster2)
res_clust2 <- results(dds_cluster2)
summary(res_clust2, alpha=0.05)

res_clust2 <- as.data.frame(res_clust2)
res_clust2$padj[is.na(res_clust2$padj)]<-1
sum(rownames(res_clust2[res_clust2$padj<0.05,]) %in% rownames(res_clust1[res_clust1$padj<0.05,]))



original_gene_list <- res_clust2$stat
names(original_gene_list) <- rownames(res_clust2)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

set.seed(987)
gse_cluster2 <- gseGO(geneList = gene_list,
                      ont = "BP",
                      keyType = "ENSEMBL",
                      minGSSize = 10,
                      maxGSSize = 800,
                      pvalueCutoff = 1,
                      verbose = T,
                      OrgDb = org.Hs.eg.db,
                      pAdjustMethod = "BH")

treeplot(enrichplot::pairwise_termsim(gse_cluster2),
         nCluster=5,
         showCategory=gse_cluster2@result$Description[1:50],
         font.size=5,
         label_format = 100,
         color="enrichmentScore")+
  scale_color_gradientn(name = "fold change",
                        colours = exp_colors,
                        limits= c(-1, 1), 
                        breaks=c(-1 , 0, 1) )


go_cluster1<-gse_cluster1@result[gse_cluster1@result$p.adjust<0.05,]
go_cluster2<-gse_cluster2@result[gse_cluster2@result$p.adjust<0.05,]

go_of_interest_cluster1<-gse_cluster1@result[gse_cluster1@result$ID %in% c(go_cluster1$ID, go_cluster2$ID),]
go_of_interest_cluster2<-gse_cluster2@result[gse_cluster2@result$ID %in% c(go_cluster1$ID, go_cluster2$ID),]

go_of_interest_cluster1$ETP<-"Cluster 1"
go_of_interest_cluster2$ETP<-"Cluster 2"

go_of_interest<-rbind(go_of_interest_cluster1, go_of_interest_cluster2)
go_of_interest<-go_of_interest[go_of_interest$ID %in% names(table(go_of_interest$ID)>1)[table(go_of_interest$ID)>1],]

divergent_go<-tapply(go_of_interest$enrichmentScore, go_of_interest$ID, prod)
divergent_go<-names(divergent_go[divergent_go<0])

#Simplify GOs
library(GOSemSim)
d<-godata("org.Hs.eg.db", ont="BP", computeIC=FALSE)
simtable<-mgoSim(divergent_go, divergent_go, semData=d, combine=NULL )
redundant_gos<-apply(simtable, 2, FUN=function(x) dimnames(simtable)[[1]][x>0.7])

go_go<-function(x) {
    if(length(x)==1) return(x)
    simtable_subset<-simtable[x,x]
    chosen<-which(simtable_subset==min(simtable_subset))
    return(na.omit(dimnames(simtable_subset)[[1]][chosen][[1]]))
  }

divergent_go_simple<-unique(unlist(lapply(redundant_gos, FUN=go_go)))


ggplot(go_of_interest[go_of_interest$ID %in% divergent_go_simple,], aes(x = enrichmentScore, y = Description, color = p.adjust, size = setSize))+
  geom_point()+
  geom_vline(xintercept=0, linetype=2)+
  theme_dose(font.size=12)+
  facet_wrap(~ETP)+
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  scale_color_gradientn(colours = sig_colors,
                        trans="log",
                        limits= c(1e-8, 0.05), 
                        breaks=c(1e-8,1e-5, 1e-3, 0.05) ) +
  labs(size="Number of genes", color="Adjusted p value")


#Cell deconvolution by cluster and surgery

samples2$cluster<-as.factor(samples2$cluster)
dds_cluster<-DESeqDataSetFromTximport(txi, colData = samples2, design = ~cluster) 

genenames<-AnnotationDbi::select(edb, keys=rownames(dds_cluster), keytype = "GENEID", columns=c("SYMBOL", "GENEID") )
rownames(dds_cluster)<-genenames$SYMBOL 

dds_cluster<-estimateSizeFactors(dds_cluster)
bulk<-counts(dds_cluster, normalized=T) #Counts matrix (with annotate's genes)
sum(rownames(bulk) %in% rownames(immunoStatesMatrix2))
bulk<-bulk[rownames(bulk) %in% rownames(immunoStatesMatrix2),] 

#unique genes
bulk<-rowsum(bulk, rownames(bulk)) 
dim(bulk)

outDT <-as.data.table(MetaIntegrator:::iSdeconvolution(immunoStatesMatrix2, bulk), keep.rownames = T)
#add new variables (compiling the different variables of each cell type)
outDT[,natural_killer_cell:=CD56bright_natural_killer_cell+CD56dim_natural_killer_cell]
outDT[,monocyte:=CD14_positive_monocyte+CD16_positive_monocyte]
outDT[,B_cell:=naive_B_cell+memory_B_cell+plasma_cell]
outDT[,T_cell:=CD8_positive_alpha_beta_T_cell+CD4_positive_alpha_beta_T_cell+gamma_delta_T_cell]
outDT[,granulocyte:=neutrophil+eosinophil+basophil]
outDT[,lymp:=B_cell+T_cell]
cell_types<-colnames(outDT)[2:17]
keep<-apply(outDT[,2:17],2, function(x) sum(x>0)>4)
cell_types<-cell_types[keep]

outDT_scaled<-as.data.frame(outDT)
outDT_scaled<-outDT_scaled[,cell_types]
outDT_scaled<-as.data.frame(t(apply(outDT_scaled,1, FUN=function(x) x/sum(x))))
rowSums(outDT_scaled)

samples2_cells<-cbind(samples2, outDT_scaled[match(samples2$samples_name, outDT$rn),])

cyto_plot_qx<-function(celltype, titulo, subtitulo) {
  cytoquine<-as.numeric(outDT_scaled[[celltype]])
  genot<-as.factor(samples2_cells$cluster)
  ggplot(NULL, aes(y=cytoquine, x=samples2_cells$condition, col=genot, fill=genot))+
    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3)+
    geom_jitter(stroke=0, size=3.5, position=position_jitterdodge())+
    theme_grey(base_size = 12)+
    theme(aspect.ratio=1, axis.text.x = element_text(angle=45, hjust = 1), legend.position="none", 
          strip.background = element_rect(colour=NA, fill=NA), strip.text=element_text(colour="black"))+
    labs(x=NULL, y=NULL, title=titulo)+
    scale_x_discrete(limits=c("pre_qx", "post_qx"), labels=c("Pre-surgery", "After surgery"))+
    scale_y_continuous(labels=scales::percent)+
    scale_color_manual(values=clusters_colors)+
    scale_fill_manual(values=clusters_colors)+
    facet_wrap(subtitulo)
}

titulos<-c("neutrophil", "eosinophil", "basophil",
           "CD14_positive_monocyte", "plasmacytoid_dendritic_cell", "CD56dim_natural_killer_cell", 
           "naive_B_cell", "memory_B_cell", "plasma_cell",
           "CD4_positive_alpha_beta_T_cell", "CD8_positive_alpha_beta_T_cell")


#ANCOVA
ancova<-function(celltype) {
  reference<-samples2_cells[samples2_cells$condition=="pre_qx", celltype]-mean(samples2_cells[samples2_cells$condition=="pre_qx",celltype])
  delta<-tapply(samples2_cells[[celltype]], samples2_cells$nhc, FUN=function(x) x[1]-x[2])
  delta<-delta[match(samples2_cells[samples2_cells$condition=="pre_qx", "nhc"], names(delta))]
  model_ancova<-summary(lm(delta~reference+samples2_cells$cluster[1:13]))
  model_ancova[["coefficients"]][3,4]
}

p_val<-round(unlist(lapply(titulos, ancova)),3)
p_val<-paste("p=", p_val, sep="")
p_val[p_val=="p=0"]<-"p<0.001"
p_val<-as.factor(p_val)



lista<-list()
for(i in 1:length(titulos)) {
  lista[[i]]<-cyto_plot_qx(titulos[[i]], titulos[[i]], p_val[[i]])
  #ggsave(filename=paste(titulos[[i]],".pdf", sep=""), plot=lista[[i]])
}

plot_grid(plotlist=lista, align="hv", nrow=4)


#Outcome

library(survival)
library(ggfortify)

clust_var <- as.factor(samples$cluster)
qx<-as.factor(samples$qx>1)
table(samples$exitus, clust_var, useNA = "always")
samples$exitus<-factor(samples$exitus, levels=c("-", "0", "1"))
table(samples$exitus)

surv.icu<-Surv(samples$follow_up, samples$exitus)

autoplot(survfit(surv.icu~clust_var)[,2], conf.int = FALSE, surv.size=2, ylim=c(0,1))+
  scale_color_manual(values=clusters_colors)+
  theme_gray(base_size = 24)+
  theme(legend.position = "none")+
  scale_x_continuous(breaks=c(0,10,20,30, 40, 50, 60))+
  labs(x="Time (days)", y="Probability of hospital discharge alive")

autoplot(survfit(surv.icu~clust_var)[,3], conf.int = FALSE, surv.size=2,  ylim=c(0,1))+
  scale_color_manual(values=clusters_colors)+
  theme_gray(base_size = 24)+
  theme(legend.position = "none")+
  scale_x_continuous(breaks=c(0,10,20,30, 40, 50, 60))+
  labs(x="Time (days)", y="Probability of death")

summary(survfit(surv.icu~clust_var), times=c(0,10,20,30, 40, 50, 60))
model<-coxph(surv.icu~clust_var+edad+sexo, data=samples, id=nhc) 
summary(model)

exitus.table<-table(samples$cluster, samples$exitus)
exitus.table
prop.table(exitus.table, margin=1)
fisher.test(exitus.table)



