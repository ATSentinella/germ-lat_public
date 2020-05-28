#+ get-phylo-correlation, eval = FALSE

getPhyloCor <- function(data) {
  
  require(ape)
  require(dplyr)
  require(tidyr)
  
  MSBP.tree <- data
  phyall <- read.tree(file = "./Data/Phylo/ALLMB.tre") #full phylo tree for plants
  
  rownames(MSBP.tree) <- MSBP.tree$Taxon_ID
  
  MSBP.Species <- as.vector(MSBP.tree$Taxon_ID)
  
  #Prune full tree for matches
  PrunedMSBP <- drop.tip(phyall, phyall$tip.label[na.omit(-match(MSBP.Species, phyall$tip.label))])
  
  MSBP.tree <- rename(MSBP.tree, "PrunedMSBP$tip.label" = Taxon_ID) %>%
               right_join(as.data.frame(PrunedMSBP$tip.label), by = "PrunedMSBP$tip.label") %>%
               rename(Taxon_ID="PrunedMSBP$tip.label")
  
  MSBP.Species <- as.vector(MSBP.tree$Taxon_ID)
  
  nSpp <- nrow(MSBP.Species)
  
  MSBP.tree[, "phylo"] <- MSBP.Species
  
  tree <- compute.brlen(PrunedMSBP)
  
  Vphy <- vcv(tree, cor = T) #correlation matrix

  return(list(MSBP.tree = MSBP.tree, cor = Vphy))
  
}
