library(ape)
library(dplyr)
library(tidyr)
library(phytools)
library(phylosignal)

#Vascular plant tree
phyall <- read.tree(file = "Data/PhylogeneticResources/Vascular_Plants_rooted.dated.tre")

#Germination warming risk data set filtered for individual species
MSBP.tree <- MSBP.meta %>%
  filter(SEbreadth<200) %>%
  mutate(FWT = Future.Temp - Tmax) %>%
  distinct(Taxon_ID, .keep_all = T)

MSBP.tree.table <- MSBP.tree

rownames(MSBP.tree) <- MSBP.tree$Taxon_ID
MSBP.Species <- as.vector(MSBP.tree$Taxon_ID)
#Prune full tree for matches
PrunedMSBP <- drop.tip(phyall, phyall$tip.label[na.omit(-match(MSBP.Species, phyall$tip.label))])
MSBP.tree <- rename(MSBP.tree, "PrunedMSBP$tip.label" = Taxon_ID)
MSBP.tree <- right_join(MSBP.tree, as.data.frame(PrunedMSBP$tip.label), by = "PrunedMSBP$tip.label")
MSBP.tree <- rename(MSBP.tree, Taxon_ID="PrunedMSBP$tip.label")
MSBP.Species <- as.vector(MSBP.tree$Taxon_ID)
MSBP.tree[, "phylo"] <- MSBP.Species
tree <- compute.brlen(PrunedMSBP)
cor <- vcv(tree, cor = T)

#List of taxa

MSBP4d <- phylo4d(tree, MSBP.tree$Topt)

teststat <- phyloSignal(MSBP4d, rep =999)

#Tests for pylogenetic signal: Blomberg's K, and Pagel's lamda
phylosig(tree, MSBP.tree$Tmax, method = "K", test = T, se = MSBP.tree$SEmax)
phylosig(tree, MSBP.tree$Tmax, method = "lambda", test = T, se = MSBP.tree$SEmax)

phylosig(tree, MSBP.tree$Tmax, method = "K", test = T, se = MSBP.tree$SEmax)
phylosig(tree, MSBP.tree$Tmax, method = "lambda", test = T, se = MSBP.tree$SEmax)

phylosig(tree, MSBP.tree$Future.Temp, method = "K", test = T)
phylosig(tree, MSBP.tree$Future.Temp, method = "lambda", test = T)

phylosig(tree, (MSBP.tree$FWT), method = "K", test = T, se = MSBP.tree$SEmax)

phylosig(tree, (MSBP.tree$FWT), method = "lambda", test = T, se = MSBP.tree$SEmax)

phylosig(tree, (MSBP.tree$FWT), method = "K", test = T)
phylosig(tree, (MSBP.tree$FWT), method = "lambda", test = T)

#Test to see the values that pop out with "dummy variables"
phylosig(tree, (MSBP.tree$Longitude), method = "K", test = T)
phylosig(tree, (MSBP.tree$AccessionNumber), method = "lambda", test = T)


source("https://raw.githubusercontent.com/liamrevell/phytools/master/R/plotTree.wBars.R")


xtree <-pbtree(n=100,b=1,d=0.5)
x<-fastBM(xtree)

rownames(MSBP.tree) <- MSBP.tree$Taxon_ID

tree[4] <- NULL 

plotTree.wBars(xtree, x, type = "fan") #Random brownian motion 


rownames(MSBP.tree) <- MSBP.tree$Taxon_ID
WTvector <- as.vector(MSBP.tree$WT)
names(WTvector) <- MSBP.tree$Taxon_ID

#why doesn't this work? But above does? Something to do with it being ultrmetric

plotTree.wBars(tree, WTvector,type = "fan", scale = 0.002) #Actual data



#Very slow
MSBP.cg <- phyloCorrelogram(MSBP4d)

plot(MSBP.cg)
pdf("phylo.cg.pdf", width = 5, height =5)

#looking for local areas of strong signal
#i.e. groups with high correlation of future warming tolerance
local.i <- lipaMoran(MSBP4d, prox.phylo = "nNodes", as.p4d = TRUE)
points.col <- lipaMoran(MSBP4d, prox.phylo = "nNodes")$p.value
points.col <- ifelse(points.col < 0.01, "red", "black")
dotplot.phylo4d(local.i, dot.col = points.col)
pdf("testphysignalplot.pdf", width = 15, height = 20)


mat.e <- matrix(MSBP.tree$SEmax, ncol = 1,
                dimnames = list(tipLabels(local.i), names(tdata(local.i))))

barplot(local.i, error.bar.sup = mat.e, error.bar.inf = mat.e)