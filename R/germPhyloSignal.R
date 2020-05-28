#+ germ-phylo-signal, eval = FALSE

library(dplyr)
require(tidyr)
require(ggtree)
require(viridis)
require(phylobase)
library("phylosignal")
require(ape)
require(phytools)


#core dataset
MSBP <- read.csv(file = "./Data/MSBP_Clean") #core dataset

#Prepare data set
MSBP.tree <- MSBP %>%
             mutate(vi = SEmax*SEmax) %>% #Standard error for future warming risk
             filter(!is.na(Tmax) & #Filter where Tmax and Future.Hot.Quart is not NA
                    !is.na(Future.Hot.Quart) &
                    SEmax < 200) %>% #and SEmax is lower than 200
             mutate(WR = Future.Hot.Quart - SEmax) %>% #Create FWT variable
             distinct(Taxon_ID, .keep_all = T)



#Vascular plant tree
phyall <- read.tree(file = "./Data/Phylo/ALLMB.tre")

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

MSBP4d <- phylo4d(tree, tip.data = MSBP.tree$WR)

WRSignal <- phyloSignal(MSBP4d, rep =999)

WRSignal
