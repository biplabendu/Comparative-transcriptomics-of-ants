install.packages(c("fields", "impute", "dynamicTreeCut", "qvalue") )
install.packages(c("flashClust", "Hmisc"))
library(flashClust)
library(cutreeStatic)

source("http://bioconductor.org/biocLite.R")
biocLite("impute")

#PART 1:  Loading and Cleaning the Data
library(WGCNA)
library(dynamicTreeCut)
options(stringsAsFactors = FALSE);

# Read in the data set
dat1=read.csv("dataset.csv", header=T)

# Take a quick look at what is in the data set
dim(dat1)
names(dat1)
datExpr0 = as.data.frame(t(dat1[, -c(1)]));
names(datExpr0) = dat1$cluster;
rownames(datExpr0) = names(dat1)[-c(1)];

# Check genes with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster the samples to inspect for outlier arrays.
sampleTree = flashClust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut (can trim outliers)
abline(h = 200000, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 300000, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Load trait data
traitData = read.csv("traits.csv", header=T);
dim(traitData)
names(traitData)

# Remove columns that hold information we do not need
allTraits = traitData;
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits
WWDSDRSamples = rownames(datExpr);
traitRows = match(WWDSDRSamples, allTraits$Samples);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();

# Re-cluster the samples
sampleTree2 = flashClust(dist(datExpr), method = "average")

# Convert traits to a color representation; white means low, red means high, grey means missing
traitColors = numbers2colors(datTraits, signed = FALSE);


# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

# Choose the soft-thresholding power for analysis of network topology
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 1;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex.axis = 3,col="green",cex.lab =20,cex.main =2)

# Now calculate the adjacencies using the soft thresholding power (i.e., softPower = beta)
softPower = 8;
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency matrix into topological overlap matrix; then convert to dissimilarity matrix (1-TOM).
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function; flastClust is much faster clustering routine than hclust
geneTree = flashClust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);

# Module identification using dynamic tree cut
minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, cutHeight = 0.97, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)


# Convert numeric labels into colors and plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

# Choose a cutoff of 0.2
MEDissThres = 0.2

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Correlate eigengenes with external traits and look for the most significant associations
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, method="pearson");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitCor=cor(MEs, datTraits, method="pearson");
moduleTraitPvalue = glm(MEs ~ datTraits)

# Graphical representation of modules-phenotypes correlations
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE,cex.text = 0.5, zlim = c(-1,1), main = paste("Module-trait relationships"))
table(moduleColors)

