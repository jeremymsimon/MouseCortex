#####################################
# Introduction
#####################################

# This is an R markdown script that outlines the key steps involved in identifying clusters of developing mouse cortex cell types from Drop-seq dataset published in
# Loo et al., "Single-cell transcriptomic catalog of mouse cortical development", bioRxiv, 2017
# Much of this was derived from Shekhar et al., "Comprehensive classification of retinal bipolar cells using single-cell transcriptomics", Cell, 2016, using the version pushed to Github on October 19th, 2016

library(genefilter)
library(sva)
library(igraph)
library(ggplot2)
library(Matrix)
library(gmodels)
library(RANN)
library(reshape)
library(cluster)
library(SLICER)
library(gplots)


# Please ensure that the accessory file class.R is available in your working directory

# First change to the directory containing the above file, and set it as your working directory. Next, load required packages and functions

source("class.R")

#####################################
# Data preprocessing
#####################################

# Load the data file containing the raw gene expression matrix

e14_dge = read.table("E14_combined_matrix.txt",header=T,sep="\t",row.names=1)

print(dim(e14_dge))

[1] 21313 11069

e14_dge[1:5, 1:2]
              e14.WT10_AAAAAGCAAGAA e14.WT10_AAAAATCTCTCC
0610005C13Rik                     0                     0
0610007N19Rik                     0                     0
0610007P14Rik                     0                     1
0610009B14Rik                     0                     0
0610009B22Rik                     0                     0

# The raw matrix thus consists of 21,313 genes and 11,069 cells. 

# Next, remove cells that contain more than 10% mitochondrially derived transcripts,

mt.genes = grep("mt-", rownames(e14_dge), value = TRUE)
cells.use = colnames(e14_dge)[colSums(e14_dge[mt.genes, ])/colSums(e14_dge) < 0.1]
e14_dge = e14_dge[, cells.use]
dim(e14_dge)
[1] 21313 11024

# Generate histogram in Supplementary Figure 1a-b
pdf("raw_colSums.pdf")
hist(colSums(e14_dge))
dev.off()

# Initialize single cell data as an S4 class object. 
# Only cells where > 500 genes are detected are considered. 
# Among the selected cells, only genes that are present in > 10 cells and those having > 60 transcripts summed across all the selected cells are considered,

dsq.e14ctx=scDrop(count.data=e14_dge)
dsq.e14ctx=initialize(dsq.e14ctx, min.genes = 500, min.cells = 10, min.counts=60)

[1] "Initializing S4 object"
[1] "Median normalizing counts and log-transforming"
[1] "z-scoring each gene"

print(dim(dsq.e14ctx@data))

[1] 12982 10931

# Thus the filtered matrix contains 12982 genes and 10931 cells. 

# Let’s examine the number of cells from each of the samples,

table(dsq.e14ctx@meta$sample)

 e14.WT10  e14.WT11 e14.WT8.1 e14.WT8.2 e14.WT9.1 e14.WT9.2 
     2376      2938      1389       866      1434      1928   

# Note that for publication, we renamed these as e14.WT1 through e14.WT6, respectively

# We can examine the genes and transcripts per cell (Supplementary Figure 1c-d)
pdf("e14_violin.pdf")
violinplot(dsq.e14ctx,c("num.genes","num.trans"))
dev.off()


#####################################
# Batch correction and PCA
#####################################

# Next, we perform batch correction using ComBat (Johnson et al., Biostatistics, 2007) and the batch-corrected expression matrix is then reduced using PCA, and the top 100 PCs are stored.

unique(as.character(dsq.e14ctx@meta$sample))
[1] "e14.WT10"  "e14.WT11"  "e14.WT8.1" "e14.WT8.2" "e14.WT9.1" "e14.WT9.2"
 
batchname = as.character(dsq.e14ctx@meta$sample)
batchid = rep(1,length(batchname))
batchid[batchname=="e14.WT11"] = 2
batchid[batchname=="e14.WT8.1"] = 3
batchid[batchname=="e14.WT8.2"] = 4
batchid[batchname=="e14.WT9.1"] = 5
batchid[batchname=="e14.WT9.2"] = 6

dsq.e14ctx=doBatchCorrection(dsq.e14ctx, batch.cov=batchid)

dsq.e14ctx=doPCA(dsq.e14ctx,pcs.store=150)

# Both steps above are computation and/or memory intensive, and we recommend that they be run on a cluster with at least 30 GB memory. 
# The batch corrected expression matrix is stored in dsq.e14ctx@batch.data. 
# The scores, loadings, and eigenvalues from the PCA are stored in dsq.e14ctx@pca.scores, dsq.e14ctx@pca.load, and dsq.e14ctx@pca.eigenvalues


# Visualize PC-scores as a scatter plot

pdf("e14_PCAplots_batchCorrected.pdf")
data.plot = dsq.e14ctx@pca.scores
data.plot$group = dsq.e14ctx@group
ggplot(data.plot, aes(x = PC1, y = PC2)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC1, y=PC2)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC1, y=PC3)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC1, y=PC4)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC1, y=PC5)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC2, y=PC3)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC2, y=PC4)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC2, y=PC5)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC3, y=PC4)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC3, y=PC5)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
ggplot(data.plot, aes(x=PC4, y=PC5)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
dev.off()


# Write out data table for permutations
write.table(dsq.e14ctx@scale.data,"dsq.e14ctx.scale.WT.batchcorrected.data.txt",quote=F,sep="\t",col.names=NA)


#####################################
# Cell-wise permutation to determine number of significant eigenvalues and therefore significant PCs to look at
# This takes a very, very long time, therefore we recommend saving the below as a standalone script and run as 50 parallelized jobs of 10.
#####################################

library(gmodels)
nperm=10
randEV=rep(NA,nperm)
data=read.table("dsq.e14ctx.scale.WT.batchcorrected.data.txt",header=T,sep="\t",row.names=1)
data=as.matrix(data)

for (i in 1:nperm) {
	perm.mat=matrix(NA,nrow=length(rownames(data)),ncol=length(colnames(data)))

	for (j in (1:length(rownames(data)))) {
		s = sample(1:length(data[j,]))
		perm.mat[j,]=rbind(data[j,s])	
	}

	rownames(perm.mat) = rownames(data)
	colnames(perm.mat) = c(1:length(colnames(data)))

	data.use=perm.mat
	pc.genes = rownames(perm.mat)
         
	#Remove genes with zero variation

	pc.genes.var = apply(data.use[pc.genes,],1,function(x) var(x))
	genes.use = pc.genes[pc.genes.var>0]
	pc.data = data.use[genes.use,]
            
	pca.obj = fast.prcomp(t(pc.data),center=FALSE, scale=FALSE)
	ev=pca.obj$sdev^2	
	randEV[i] = max(ev)
	print(max(ev))
}
print(randEV)

#############


# Determine how many "significant" PCs there are based on maximum eigenvalue from above permutation
length(dsq.e14ctx@pca.eigenvalues$ev[dsq.e14ctx@pca.eigenvalues$ev>4.38])
[1] 89


# Visualize PC-scores for the top X "significant" PCs as a scatter plot

pdf("e14_PCAplot_batchCorrected_postReduction.pdf")
data.plot = dsq.e14ctx@pca.scores[,1:89]
data.plot$group = dsq.e14ctx@group
ggplot(data.plot, aes(x = PC1, y = PC2)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
dev.off()

# Write out PCA scores to load into t-SNE (rows are cell IDs, columns are top X significant PCs)
write.table(dsq.e14ctx@pca.scores[,1:89],"DropSeqAnalysis_e14_batchCorrected_PCAdimReduced_top89PCs.txt",quote=F,sep="\t",col.names=NA)


#####################################
# Visualize in tSNE
#####################################

# Now that we have reduced the number of dimensions of our dataset down to 89 (PCs) x 10,931 (cells), we can visualize them using t-SNE
# We use the scikit-learn implementation in Python, but this could be done in R as well if users prefer
# See "E14_tSNE.py" for code to run

# We can specify different values for perplexity and learning_rate to get different results, as illustrated here: https://distill.pub/2016/misread-tsne/
# We therefore iterated the perplexity parameter, testing the following values: 5, 10, 20, 25, 30, 50, 60, 100, and 150.
# To be clear, this has no bearing on the downstream analyses, only for visualization purposes
# We set perplexity=50 for this dataset

# It is also good practice to iterate through the learning_rate values, we tested the following: 100, 250, 500, 600, 700, 750, 800, 900, and 1000.
# Again, this has no bearing on the downstream analyses, only for visualization purposes
# We set learning_rate=750 for this dataset


# Once that is complete, we directly read in the embedded values and store in the slot dsq.e14ctx@tsne.y

tsne.y = read.table("DropSeqAnalysis_e14_batchCorrected_PCAdimReduced_top89PCs_perp50_learn750_tSNE.txt",sep="\t",header=T,row.names=1)
colnames(tsne.y) = c("tSNE1","tSNE2")
dsq.e14ctx@tsne.y = tsne.y

# Next, we visualize single-cell expression based on tSNE coordinates, as in Supplementary Figure 2.

pdf("DropSeqAnalysis_e14_batchCorrected_PCAdimReduced_top89PCs_tSNE_perp50_learn750_tSNE_unclustered.pdf")
data.plot = dsq.e14ctx@tsne.y
data.plot$group = dsq.e14ctx@group
ggplot(data.plot, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
dev.off()


#####################################
Louvain-Jaccard Clustering
#####################################

# Next, we cluster the cells based on their PC scores using the Louvain-Jaccard method. The basic command recommended by Shekhar et al. is: 
# dsq.e14ctx = doGraph_clustering(dsq.e14ctx, pcs.use = 1:89, num.nn = 30, do.jaccard = TRUE, method = "Louvain")

# Except determining the proper number of nearest neighbors to use is non-trivial, and does make a large difference in the final cluster assignments and cluster robustness
# We therefore want to optimize this parameter to find the best #NN to use that produces the optimal baseline cluster assignments

# To do this, we iterate through #NNs from 10-100, plotting mean silhouette width as a measure of cluster robustness for each iteration. 
# Silhouette widths here are computed using 1-cor() on values for all genes in PC space. Typically spearman correlations work best, however this may be changed in favor of "pearson" or "kendall"

filePrefix = "DropSeqAnalysis_e14_batchCorrected_PCAdimReduced_top89PCs_tSNE_perp50_learn750"
corMethod = "spearman"
pcs.use = 1:89
print("Optimizing number of nearest neighbors (NN) from 10 through 100")

sils=rep(NA,100)
numclust=rep(NA,100)
for (i in 10:100) {		# Range of NNs to use is from 10 through 100. First 9 values will be set to NA in final sils vector
	dsq.e14ctx = doGraph_clustering(dsq.e14ctx, pcs.use = pcs.use, num.nn = i, do.jaccard = TRUE, method = "Louvain")
	numclust[i] = length(table(dsq.e14ctx@group))
	data.use=dsq.e14ctx@pca.scores[,pcs.use]
	data.use=t(data.use)
	sub<-as.numeric(dsq.e14ctx@group)
	dm<-as.dist(1-cor(data.use,method=corMethod))
	si<-silhouette(sub,dm,data.use)
	sils[i] = mean(si[,3])
}
dsq.e14ctx@sils = sils
dsq.e14ctx@numclust = numclust

min.sils = round(min(dsq.e14ctx@sils,na.rm=T),2) - 0.01
max.sils = round(max(dsq.e14ctx@sils,na.rm=T),2) + 0.01
min.clust = round(min(dsq.e14ctx@numclust,na.rm=T),2) - 5
max.clust = round(max(dsq.e14ctx@numclust,na.rm=T),2) + 5
opt.NN = which.max(dsq.e14ctx@sils)

pdfName = paste0(filePrefix,"_sils_numclust_iter_10-100_",corMethod,".pdf")
pdf(file=pdfName)
print(ggplot(as.data.frame(dsq.e14ctx@sils), aes(x=seq(1:100),y=dsq.e14ctx@sils),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.5) + scale_x_continuous(breaks=c(seq(0,100,10)),limits=c(0,100),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.sils,max.sils,0.01)),limits=c(min.sils,max.sils),minor_breaks=NULL))
print(ggplot(as.data.frame(dsq.e14ctx@numclust), aes(x=seq(1:100),y=dsq.e14ctx@numclust),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,100,10)),limits=c(0,100),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.clust,max.clust,5)),limits=c(min.clust,max.clust),minor_breaks=NULL))
dev.off()

# This process determines that the optimal NN value to use for this dataset is 34.


# Next, we could generate cluster assignments by now running "doGraph_clustering" using that optimized NN value, however these baseline cluster assignments are still imperfect
# Most notably, there are actually many statistical outliers per cluster.
# Instead, we created a new function, "perform_refined_clustering", that runs the clustering iteratively
# During each iteration, we compute Silhouette widths to assess cluster robustness and look for outlier cells
# This iteration refines these cluster assignments, plus allows for new clusters to form comprised of outlier cells
# New clusters can form if at least 10% of cells (and >10 cells) within a cluster were considered outliers AND that they themselves formed a new robust grouping. 
# Otherwise no new cluster will be formed and they instead will get a chance to join their closest neighbor in the subsequent iteration.

# This function can also run the above NN optimization first, if desired, by specifying "optimize.NN=TRUE"
# Since we already determined the optimal NN, we will specify "optimize.NN=FALSE" and provide "opt.NN=34" to simply refine the existing clusters
# Plots containing Silhouette widths for each cell and cluster will be output during the iteration process for the user to monitor changes to cluster assignments

dsq.e14ctx = perform_refined_clustering(dsq.e14ctx, pcs.use=1:89, corMethod="spearman",filePrefix="DropSeqAnalysis_e14_batchCorrected_PCAdimReduced_top89PCs_tSNE_perp50_learn750", min.newcluster.size=10, min.proportion.outliers=0.10, optimize.NN=FALSE,opt.NN=34)

# Visualize final clusters on tSNE plot (Figure 1)

pdf("DropSeqAnalysis_e14_batchCorrected_PCAdimReduced_top89PCs_tSNE_perp50_learn750_34NNclustered_tSNE.pdf")
plot.tsne(dsq.e14ctx)
dev.off()


#####################################
# Find and plot marker genes of each cell cluster
#####################################

# We are now ready to identify marker genes of each cluster to aid in determining their biological identity. 
# First we generate matrices corresponding to the raw counts Count.mat and median normalized counts TPM.mat

filePrefix = "DropSeqAnalysis_e14_batchCorrected_PCAdimReduced_top89PCs_tSNE_perp50_learn750_cellStates_cluster"
TPM.mat = exp(dsq.e14ctx@data) - 1
percExp = rep(NA,length(rownames(dsq.e14ctx@data)))
for (j in 1:length(levels(dsq.e14ctx@group))) {
	markers = paste("markers", j, sep = "")
	assign(markers, markers.binom(dsq.e14ctx, clust.1 = j, effect.size = log(2), TPM.mat = TPM.mat, Count.mat = Count.mat))
	print(head(get(markers),10))
	write.table(get(markers),paste0(filePrefix,"_markerGenes_clust",j,".txt"),quote=F,sep="\t",col.names=NA)
}


# Make dotplot to illustrate cell type specificity of marker gene expression

e14genes=c("Dcx","Vim","Nes","Slc17a6","Bcl11b","Fezf2","Tbr1","Reln","Trp73","Lhx1","Lhx5","Crym","Ndrg1","Mc4r","Nfe2l3","Nxph3","Sla","Sybu","Tfap2d","Rwdd3","Ppp1r14c","Cux2","Satb2","Tiam2","Pdzrn3","Eomes","Nhlh1","Unc5d","Btbd17","Mki67","Top2a","Sgol2","Spag5","Hes1","Hes5","Ednrb","Tk1","Tcf19","Pkmyt1","Arhgef39","Rspo1","Dkk3","Dlx1","Dlx2","Gad1","Gad2","Pnoc","Lhx6","Nxph1","Isl1","Syt6","Gpr88","Rxrg","Tcf7l2","Shox2","Syt13","Ak7","Folr1","Foxj1","Kcne2","Otx2","Trem2","Mpeg1","Igfbp7","Epas1","Pecam1","Col3a1")
pdf("DropSeqAnalysis_e14_WTonly_batchCorrected_WTonly_PCAdimReduced_top89PCs_tSNE_perp50_learn750_34NNclustered_dotplot_named.pdf",width=30,height=10)
a=c(17,5,13,3,7,2,4,11,15,8,14,21,10,6,1,12,16,9,19,22,20,18)
groups=c("LayerI","LayerV-VI","LayerV-VI","LayerV-VI","LayerV-VI","LayerV-VI","SVZ1(migrating)","SVZ2(migrating)","SVZ3(proliferating)","RG1","RG2","RG3 (cortical hem)","RG4","Ganglionic eminences","Int1","Int2","Striatal inh1","Striatal inh2","Thalamic","Choroid plexus","Microglia","Endothelial")
dot.plot(dsq.e14ctx, features.use = e14genes,group.use=a,group.names=groups,max.val.perc = 0.9, max.val.exp = 5, max.size = 10)
dev.off()




#####################################
# Mean collapsing of clusters, for web-based visualization tools
#####################################

dsq.e14ctx.scaled.collapsed = matrix(NA,nrow=length(rownames(dsq.e14ctx@scale.data)),ncol=length(levels(dsq.e14ctx@group)))
rownames(dsq.e14ctx.scaled.collapsed) = rownames(dsq.e14ctx@scale.data)
colnames(dsq.e14ctx.scaled.collapsed) = paste0("Cluster",seq(1:22))

for (i in 1:length(levels(dsq.e14ctx@group))) {
	samples.use = unique(as.character(names(dsq.e14ctx@group[which(dsq.e14ctx@group==i)])))
	submat = dsq.e14ctx@scale.data[,which(colnames(dsq.e14ctx@scale.data) %in% samples.use)]
	submat.collapse = apply(submat,1,mean)
	dsq.e14ctx.scaled.collapsed[,i] = submat.collapse
}

write.table(dsq.e14ctx.scaled.collapsed,"dsq.e14ctx.scaled.collapsed.txt",quote=F,sep="\t",col.names=NA)




#####################################
# Identification of cell sub-clusters
#####################################

# To identify cellular sub-clusters, we focus on genes whose expression is heterogeneous within a cluster
# To do this, we first identify genes expressed in at least 25% but less than 75% of cells
# Then, to refine this further, we use a feature-selection tool as part of the SLICER package, which--though still unsupervised--identifies genes with somewhat coherent expression patterns among sub-clusters of cells
# If there were at least 20 genes identified this way, we then perform hierarchical clustering on the expression of those genes for all cells within a cluster to identify putative sub-clusters
# These putative sub-clusters may be driven by only 1-2 genes or very small subsets of cells, so we inspect each cluster and putative sub-clusters and cull those that do not represent robust groups
# We also include a function below for median-centering the expression values, which we perform for each cluster

medianCtr = function(x){
	namesAll = dimnames(x)
	medians = apply(x,1,median,na.rm=T)
	x = t(scale(t(x),center=medians,scale=F))
	dimnames(x) = namesAll
	return(x)
}

percExp = rep(NA,length(rownames(dsq.e14ctx@data)))
TPM.mat = exp(dsq.e14ctx@data) - 1

filePrefix = "CellSubclusters/E14_cellSubcluster"
dir.create("CellSubclusters", showWarnings = FALSE)

for (j in 1:length(levels(dsq.e14ctx@group))) {
	samples.use = unique(as.character(names(dsq.e14ctx@group[which(dsq.e14ctx@group==j)])))
	for (i in 1:length(rownames(dsq.e14ctx@data))) {
		percExp[i] = sum(dsq.e14ctx@data[i,samples.use]>0)/length(dsq.e14ctx@data[i,samples.use])
	}
	sub = TPM.mat[which(percExp>.25 & percExp <.75),samples.use]
	sub.ctr = medianCtr(sub)
	write.table(sub.ctr,paste0(filePrefix,j,"_TPMdata_medCtr.txt"),quote=F,sep="\t",col.names=NA)

	sub.ctr.t=t(sub.ctr)
	genes=select_genes(sub.ctr.t)

	if(length(genes)>=20){
		sub.slicer=t(sub.ctr.t[,genes])
		write.table(sub.slicer,paste0(filePrefix,j,"_TPMdata_medCtr_subclusters.txt"),quote=F,sep="\t",col.names=NA)
		png(paste0(filePrefix,j,"_TPMdata_medCtr_subclusters.png"))
		heatmap.2(sub.slicer,hclust=function(x) hclust(x,method="average"),distfun=function(x) as.dist(1-cor(t(x))),labCol=NA,col=colorRampPalette(c("blue", "white", "red"))(n = 50),breaks=c(seq(-3,3,0.12)),symbreaks=T,key=F,trace="none",density.info="none")
		dev.off()
	}
}
