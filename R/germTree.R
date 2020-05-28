#+ germ-tree, eval = FALSE

require(ggplot2)
require(dplyr)
require(tidyr)
require(ape)
require(ggtree)
require(viridis)

MSBP <- read.csv(file = "./Data/MSBP_Clean") #core dataset

#Prepare data set, if SE's are too high, the analyses don't run (but get weighted out anyway)
MSBP.tree.WT <- MSBP %>%
  mutate(vi = SEmax*SEmax) %>% #Standard error for warming risk
  filter(!is.na(Tmax) &
         !is.na(Future.Hot.Quart)) %>% #Filter NAs
  distinct(Taxon_ID, .keep_all = T) #Remove duplicate species

#Seed plant phylogeny - Smith and Brown 2019
phyall <- read.tree(file = "./Data/Phylo/ALLMB.tre")

#List of taxa
rownames(MSBP.tree.WT) <- MSBP.tree.WT$Taxon_ID

MSBP.Species.WT <- as.vector(MSBP.tree.WT$Taxon_ID)

#Prune full tree for matches
PrunedMSBP <- drop.tip(phyall, phyall$tip.label[na.omit(-match(MSBP.Species.WT, phyall$tip.label))])
MSBP.tree.WT <- rename(MSBP.tree.WT, "PrunedMSBP$tip.label" = Taxon_ID)
MSBP.tree.WT <- right_join(MSBP.tree.WT, as.data.frame(PrunedMSBP$tip.label), by = "PrunedMSBP$tip.label")
MSBP.tree.WT <- rename(MSBP.tree.WT, Taxon_ID="PrunedMSBP$tip.label")
MSBP.Species.WT <- as.vector(MSBP.tree.WT$Taxon_ID)
MSBP.tree.WT[, "phylo"] <- MSBP.Species.WT

tree.df <- data.frame(Taxon_ID = PrunedMSBP$tip.label)

tree.df <- data.frame(left_join(tree.df, MSBP.tree.WT, by = "Taxon_ID"))

pruned.tree <- ggtree(PrunedMSBP)


p.tree <- pruned.tree %<+% filter(tree.df) +
  geom_segment2(aes(subset=isTip, 
                    yend = y, 
                    xend = x+70, 
                    color = (Future.Hot.Quart) - Tmax, 
                    alpha=(1/log(SEmax))), 
                size = 2.5)+
  scale_colour_viridis(option = "plasma")+
  theme(legend.position=c(0,0.7))+
  guides(alpha=FALSE) #removes alpha key

pnames <- geom_tiplab2(size = 2, offset = 30)

pclades <- p.tree + 
  theme(plot.margin = unit(c(2,5,2,2), "cm"))+

coord_cartesian(xlim=c(0,1000))

pclades

ggsave(filename = "./Outputs/MSBPtreelabels.png", plot = pclades, height = 15, width =10)

