
# LOAD DEs ---------------------------

rm(list = ls())
df_all <- readRDS("df_wgcna.rds")

# LOAD LIBRARIES ---------------------------

lib <-  c("openxlsx", "heatmap3", "WGCNA", "dplyr", "readxl")

invisible(lapply(lib, library, character.only = TRUE))


# PREP DATA ---------------------------

allDat <- df_all %>%
  mutate(Protein.IDs = rownames(df_all)) %>% 
  dplyr::select(Protein.IDs, everything())


rownames(allDat) <- allDat$Protein.IDs
Group <- read_excel("groups_wgcna.xlsx")
Group <- Group$Group


# DEFINE FUNCTIONS ---------------------------

plotErrorBarsLines <- function (v, barSizes, lines, labels = NULL, col = "blue", 
                                ylim = c(min(lines), max(lines)), ...) 
{
  barSizes[is.na(barSizes)] <- 0
  topBars <- v + 0.5 * barSizes
  bottomBars <- v - 0.5 * barSizes
  N <- length(v)
  if (is.null(labels)) 
    labels <- 1:N
  ylims <- c(min(bottomBars, ylim[1], min(lines)), max(topBars, 
                                                       ylim[2], max(lines)))
  par(pch = 19, xaxt = "n")
  plot(as.numeric(labels), v, ylim = ylims, col = col, type = "b", 
       lwd = 3, ...)
  par(xaxt = "s")
  
  for (i in 1:N) {
    lines(c(i, i), c(topBars[i], bottomBars[i]))
  }
  for (i in 1:ncol(lines)) {
    lines(as.numeric(labels), lines[, i], lwd = 0.5, lty = "dotted", 
          col = "gray")
  }
}

plotClusterProfileWGCNA <- function(cluster.data, moduleColors, group, MEs=NULL, 
                                    ylab="Abundance", 
                                    file="ClusterPatterns.png", ...) {
  
  gp = group
  noClusters <- nlevels(as.factor(moduleColors))
  
  r.temp <- aggregate(t(cluster.data), by=list(gp=gp), FUN=mean)
  ag.sample <- r.temp[,-1]
  rownames(ag.sample) <- r.temp[,1]
  
  ag.sample <- ag.sample[c(1,2,5,3,4),]
  
  ag.genes <- aggregate(t(ag.sample), by=list(Cluster=moduleColors), FUN=mean)
  ag.sd <- aggregate(t(ag.sample), by=list(Cluster=moduleColors), FUN=sd)
  ag.matrix <- as.matrix(ag.genes[,-1])
  
  if(!is.null(MEs) ) {		
    r.temp <- aggregate(MEs, by=list(gp=gp), FUN=mean)
    ag.matrix <- t(r.temp[,-1])
    colnames(ag.matrix) <- r.temp[,1]
  }
  ag.counts <- summary(as.factor(moduleColors))
  ag.bars <- as.matrix(ag.sd[,-1])
  
  fScale = max(8,noClusters)/8
  
  
  png(file, 2000, 3000*fScale, res=300)
  par(fg=gray(0.3), mar= c(8, 6, 2, 1) + 0.1, col.main="black", col.sub="black", col.lab="black", col.axis="black")
  layout(matrix(1:(ceiling(noClusters/2)*2), ncol=2, byrow=TRUE))
  NSig <- noClusters
  cols = c('gray20', 'salmon4', 'lightgreen', 'gray65', 'magenta3', 'pink', 'indianred1')
  for(i in 1:NSig) {
    gname <-  paste(levels(as.factor(moduleColors))[i], "(", ag.counts[i], "proteins )")
    lines <- ag.sample[, moduleColors==levels(as.factor(moduleColors))[i], drop=FALSE]
    plotErrorBarsLines(ag.matrix[i,], 2*ag.bars[i,], lines, 
                       labels=1:ncol(ag.matrix), 
                       col=cols[i],  main=gname, # bgcol="gray", split=split,
                       ylab=ylab, xlab="", xaxt="none",
                       ylim=c(min(ag.matrix), max(ag.matrix)), ...)
    text(1:length(unique(gp)), par("usr")[3]-0.5, 
         srt = 60, adj = 1, xpd = TRUE,
         labels = c("Control","IVIg","Day 7 - Passive ITP", "Active ITP", "Day 3 - Passive ITP"), cex = 1.1, col = "black")
    # axis(1,at=1:ncol(ag.matrix), las=2, labels=colnames(ag.matrix), col="black", ...)
    abline(h=mean(as.matrix(cluster.data)), lty="dotted")
  }
  
  dev.off()
  
}


# DEF PARAMETERS ---------------------------
RCutoff = .85
MCutheight = .1
PowerUpper = 30
minModuleSize = 10

IgnoreCols = 1 	


enableWGCNAThreads()

# Only abundance data - samples as rows and proteins as columns
datExpr <- as.data.frame(t(allDat[,-c(IgnoreCols)]) )

rownames(datExpr) <- colnames(allDat)[-c(IgnoreCols)]
colnames(datExpr) <- rownames(allDat)


# BUILD WGCNA NETWORK ---------------------------

# Choose a set of soft-thresholding powers
powers <- c(c(2:10), seq(from = 12, to=PowerUpper, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = RCutoff, verbose = 5)

# Plot the results
#png("ScaleFreeTopology.png", 3000,2000,res=300)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
h = 0.85
abline(h=h,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")	 

#dev.off()

# Get the soft power
softPower <- sft$powerEstimate

# Build the adjacency table - use "signed" for proteomics data
adjacency <- adjacency(datExpr, power = softPower, type="signed")


# Turn adjacency into topological overlap distance
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM


# Clustering using TOM-based dissimilarity

proTree <- hclust(as.dist(dissTOM), method = "average")


# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = proTree, distM = dissTOM, 
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

print("Dynamic tree cut results:")
print(table(dynamicMods))


# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# c('gray20','skyblue2', 'salmon4', 'lightgreen', 'gray65', 'magenta3', 'pink', 'indianred1', 'turquoise3', 'lightgoldenrod1')

# change to more pretty colors
dycol <- dynamicColors
dycol <- ifelse(dycol == "black", "gray20", dycol)
dycol <- ifelse(dycol == "blue", "skyblue2", dycol)
dycol <- ifelse(dycol == "brown", "salmon4", dycol)
dycol <- ifelse(dycol == "green", "lightgreen", dycol)
dycol <- ifelse(dycol == "grey", "gray65", dycol)
dycol <- ifelse(dycol == "magenta", "magenta3", dycol)
dycol <- ifelse(dycol == "red", "indianred1", dycol)
dycol <- ifelse(dycol == "turquoise", "turquoise3", dycol)
dycol <- ifelse(dycol == "yellow", "lightgoldenrod1", dycol)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(proTree, dycol, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Protein dendrogram and module colors")

# Merge clusters
mergedClust <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MCutheight, verbose = 3)

mergedColors <- mergedClust$colors

mergedMEs <- mergedClust$newMEs


#############################################
# Calculate and plot module eigenproteins
# - Dendrogram
# - Heatmap
# - Boxplot
# - KME
#############################################


# Rename to moduleColors
moduleColors <- mergedColors


print("Modules after merging:")
print(table(moduleColors))

# Change to more pretty colors
dycol2 <- mergedColors
dycol2 <- ifelse(dycol2 == "black", "gray20", dycol2)
dycol2 <- ifelse(dycol2 == "brown", "salmon4", dycol2)
dycol2 <- ifelse(dycol2 == "green", "lightgreen", dycol2)
dycol2 <- ifelse(dycol2 == "grey", "gray65", dycol2)
dycol2 <- ifelse(dycol2 == "magenta", "magenta3", dycol2)
dycol2 <- ifelse(dycol2 == "red", "indianred1", dycol2)

# Plot dendrogram

#png("DendroColorMergedClust.png", 2000, 2000, res=300)
plotDendroAndColors(proTree, cbind(dycol, dycol2),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Get the module eigenproteins
MEs <- mergedMEs

rownames(MEs) <- rownames(datExpr)


# reorder MEs by color names of modules
MEs <- MEs[,order(colnames(MEs))]

# Plot module profiles with eigenproteins overlaid
WGCNAClusterID <- moduleColors

# plotClusterProfileWGCNA(t(datExpr), WGCNAClusterID, Group,  MEs= MEs,
#                         ylab="Average log ratio", file="WGCNAClusterPattenME.png",
#                         cex.main=1.3, cex.lab=1.2, cex.axis=1.1)
# 


# plotClusterProfileWGCNA(t(datExpr), WGCNAClusterID, Group,
#                         ylab="Average log ratio", file="WGCNAClusterPattenAve.pdf",
#                         cex.main=1.3, cex.lab=1.2, cex.axis=1.1)


# dendrogram and heatmap for eigenproteins
#png("Dendrogram eigenproteins.png", 2000,2000,res=300)					
plotEigengeneNetworks(MEs, "Eigenprotein Network", 
                      marHeatmap = c(3,4,2,2), 
                      marDendro = c(3,4,2,5),
                      plotDendrograms = TRUE, 
                      xLabelsAngle = 90,
                      heatmapColors=blueWhiteRed(50),
                      excludeGrey = TRUE)	#The grey module comprises the set of genes which have not been clustered in any module
#dev.off()
#
#removeGreyME()

MEs2 <- MEs
rownames(MEs2) <- c("Day 7 - Control 1",
                    "Day 3 - Control 2",
                    "Day 7 - Control 2",
                    "Day 3 - Control 3",
                    "Day 7 - Control 3",
                    "Day 3 - Passive ITP 4",
                    "Day 3 - Passive ITP 5",
                    "Day 3 - Passive ITP 6",
                    "Day 7 - Passive ITP 4",
                    "Day 7 - Passive ITP 5",
                    "Day 7 - Passive ITP 6",
                    "Active ITP 2",
                    "Active ITP 3",
                    "Active ITP 5",
                    "IVIg 10",
                    "IVIg 6",
                    "IVIg 7",
                    "IVIg 8",
                    "IVIg 9")

Group_n <- Group

Group_n <- ifelse(Group_n == "PTI_passive_3", "Day 3 - Passive ITP", Group_n)
Group_n <- ifelse(Group_n == "PTI_passive_7", "Day 7 - Passive ITP", Group_n)
Group_n <- ifelse(Group_n == "PTI_active", "Active ITP", Group_n)


# png("Heatmap_eigenproteins.png", 1000,700,res=100)
heatmap3(t(MEs2), margins=c(10,10), distfun = function(x) dist(x, method="euclidean"), 
         ColSideColors=ggsci::pal_locuszoom("default", alpha = 1)(nlevels(as.factor(Group_n)))[as.factor(Group_n)],
         method = "average")
legend("topright", fill=ggsci::pal_locuszoom("default", alpha = 1)(nlevels(as.factor(Group_n)))[1:nlevels(as.factor(Group_n))],
       legend=levels(as.factor(Group_n)), cex=.5, xpd=TRUE, inset=-.0003 )
# dev.off()


kmes <- signedKME(datExpr, MEs)


# separate results by modules, order by kME, hub proteins on top

dat.res <- data.frame(allDat, moduleColors , kmes)

list.cluster.dat <- lapply(levels(as.factor(moduleColors)), 
                           function(x) {dtemp = dat.res[dat.res$moduleColors == x,];
                           dtemp[order(dtemp[,paste0('kME',x)==colnames(dtemp)], decreasing=TRUE),
                                 -setdiff(grep("^kME", colnames(dtemp)), which(paste0('kME',x)==colnames(dtemp)))]} )

names(list.cluster.dat) <- levels(as.factor(moduleColors))



# Boxplot for eigenproteins

ag.temp <- aggregate(MEs, by=list(Group=Group), FUN=mean)
ag.eigengenes <- t(ag.temp[,-1])
colnames(ag.eigengenes) <- ag.temp[,1]

fScale <- max(8,nlevels(as.factor(moduleColors)))/8

#png("Boxplot_eigenproteins.png", 2000, 3000*fScale, res=300)

par(mar= c(7, 4, 2, 1) + 0.1)
layout(matrix(1:(ceiling(nlevels(as.factor(moduleColors))/2)*2), ncol=2, byrow=TRUE))
cols = c('gray20', 'salmon4', 'lightgreen', 'gray65', 'magenta3', 'pink', 'indianred1')
for(ii in 1:ncol(MEs)){
  boxplot(MEs[,ii] ~ Group, col=cols[ii], ylab = "log ratio",
          main=paste(colnames(MEs)[ii], table(moduleColors)[ii]), cex.main=1.1, cex.lab=1, xlab="", xaxt="none")
  text(1:length(unique(Group)), par("usr")[3]-0.05, 
       srt = 60, adj = 1, xpd = TRUE,
       labels = c("Control","IVIg","Active ITP", "Day 3 - Passive ITP", "Day 7 - Passive ITP"), cex = 0.9, font = 2)}

#dev.off()

# Boxplot for top 6 hub proteins

for(ii in 1:length(list.cluster.dat)) {
  
  # png(paste0("Boxplot hub proteins - ", names(list.cluster.dat)[ii], ".png"), 2000, 2500, res=300)
  par(oma= c(5, 2, 2, 1) + 0.1)
  layout(matrix(c(1:6), ncol=2))
  for(jj in c(1:6)){
    boxplot(as.numeric(t(list.cluster.dat[[ii]][jj,2:20])) ~ Group, 
            main=paste(list.cluster.dat[[ii]][jj,1],"\nkME=", round(list.cluster.dat[[ii]][jj,ncol(list.cluster.dat[[ii]])],2)), 
            col=rainbow(nlevels(as.factor(Group))), ylab="Log ratio", cex.main=1.5, las=2,cex.lab=1.2, alpha = 0.3)
  }
  dev.off()
}


# OUTPUT ---------------------------

# wb <- createWorkbook()
# 
# addWorksheet(wb, "AllData")
# 
# writeData(wb, "AllData", dat.res)
# 
# # write modules only tabs
# 
# for(ii in 1:length(list.cluster.dat)) {
#   addWorksheet(wb, names(list.cluster.dat)[ii])
#   writeData(wb, names(list.cluster.dat)[ii], list.cluster.dat[[ii]])
# 
# }
# 
# saveWorkbook(wb, "ResultsWGCNA.xlsx", overwrite=TRUE)
# 
# # output the eigenprotein
# write.csv(MEs,"Module_eigenproteins.csv")



# TRAITS ---------------------------
Groups <- read_excel("groups_traits.xlsx")

datTraits <- Groups[, -1]
rownames(datTraits) <- Groups$Sample

# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)


MEs0 <- orderMEs(MEs)
modTraitCor <- cor(MEs0, datTraits, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nSamples)


textMatrix <- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", 
                    sep = "")
dim(textMatrix) <- dim(modTraitCor)

#png("corTraits2.png", 2000, 2000, res=300)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
coll <- paste("ME", sep = "", c('salmon4', 'lightgreen', 'pink', 'magenta3', 'gray20', 'indianred1', 'gray65'))

labeledHeatmap(Matrix = modTraitCor, xLabels = c("IVIg", "Day 7 - Passive ITP","Active ITP", "Day 3 - Passive ITP"), yLabels = coll, 
               ySymbols = names(MEs0), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = TRUE, cex.text = 0.5, main = paste("Module-trait relationships"), cex.lab = 0.7)

#dev.off()
