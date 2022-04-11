#' This is code to run analayses and create outputs for the effect of
#' Longevity and growth form on current climate mismatch
#+ woodiness-longevity-climate-mismatch-current, eval = FALSE


getLifeWoodCM <- function() {
  
  require(dplyr)
  require(tidyr)
  require(metafor)
  require(ggplot2)
  source("./R/getPhyloCor.R")


longevity <- read.csv(file = "./Data/Longevity/TRY.longevity.csv") #longevity dataset

MSBP <- read.csv(file = "./Data/MSBP_Clean") #core dataset

#Prepare data set, if SE's are too high, the analyses don't run (but get weighted out anyway)
MSBP <- MSBP %>%
        mutate(vi = SEopt.upp*SEopt.upp) %>% #Standard error for climate mismatch
        filter(!is.na(Topt.upp)          &    
               !is.na(Current.Hot.Quart) &
               !is.na(SeedAge) &
               !is.na(Altitude)) #Remove NAs

#Take longevity data and seperate into "annual" and "perennial"
longevityfilt <- longevity %>%
  mutate(
    ann_per = case_when(
      OrigValueStr == "annual" |
        OrigValueStr == "Annual" |
        OrigValueStr == "annuals" |
        OrigValueStr == "summer annuals" |
        OrigValueStr == "always annual" |
        OrigValueStr == "winter annuals" |
        OrigValueStr == "annual-winter annual" |
        OrigValueStr == "winter annual" |
        (OriglName == "Life history" &
           OrigValueStr == "1") |
        (OriglName == "Plant phenology: Annual" &
           OrigValueStr == "yes")     ~   "annual", #annuals
      OrigValueStr == "perennial" |
        OrigValueStr == "Perennial" |
        OrigValueStr == "perennials" |
        OrigValueStr == "always pluriennial-pollakanthic" |
        (OriglName == "Plant phenology: Biennial" &
           OrigValueStr == "yes") |
        OrigValueStr == "perennial < 20 years" |
        OrigValueStr == "woody" |
        OrigValueStr == "perennial/woody" |
        OrigValueStr == "perennial > 20 years" |
        OrigValueStr == "poly-annuals > 50 years (long-lived perennials)" |
        OrigValueStr == "always biennial, always pluriennial-hapaxanthic" |
        OrigValueStr == "always biennial, always pluriennial-pollakanthic" |
        OrigValueStr == "tree" |
        OrigValueStr == "shrub" |
        OrigValueStr == "always pluriennial-hapaxanthic, always pluriennial-pollakanthic" |
        OrigValueStr == "always pluriennial-hapaxanthic" |
        OrigValueStr == "biennial" |
        OrigValueStr == "annual/biennial" |
        OrigValueStr == "poly-annuals < 5 years (short-lived perennials)" |
        OrigValueStr == "Biennial" |
        OrigValueStr == "biennial/perennial" |
        OrigValueStr == "always biennial" |
        OrigValueStr == "biennial-perennial" |
        OrigValueStr == "sometimes biennial, always pluriennial-hapaxanthic, sometimes pluriennial-pollakanthic" |
        OrigValueStr == "sometimes biennial, sometimes pluriennial-hapaxanthic, always pluriennial-pollakanthic" |
        OrigValueStr == "biennial/perennial/woody" |
        OrigValueStr == "sometimes biennial, always pluriennial-pollakanthic" |
        OrigValueStr == "poly-annuals 5-50 years (medium-lived perennials)" |
        (OriglName == "Plant phenology: Perennial" &
           OrigValueStr == "yes") |
        (OriglName == "Plant phenology: Annual" &
           OrigValueStr == "no")    ~   "perennial"  #perennials
    )
  )

longevityfilt <- data.frame(lapply(longevityfilt, function(x) {gsub(" ", "_", x)})) %>%
                 distinct(SpeciesName, .keep_all = T) %>%
                 select(SpeciesName, ann_per)

#Add Longevity to MSBP core dataset
MSBPlong <- left_join(MSBP, longevityfilt, by = c("Taxon_ID" ="SpeciesName"))

##Get woodiness information
woody <- read.csv(file = "./Data/GlobalWoodinessDatabase.csv")

woody <- data.frame(lapply(woody, function(x) {gsub(" ", "_", x)}))%>%
          distinct(gs, .keep_all = T) %>%
         select(gs, woodiness)

MSBPlong <- left_join(MSBPlong, woody, by = c("Taxon_ID" ="gs"))

#Assign any taxa with "w" (woody) to be perennial
MSBPlong <- MSBPlong %>%
            mutate(ann_per=replace(ann_per, woodiness=="W", "perennial"))

#Filter dataset for analysis (without phylogeny)
MSBP.meta.wood <- mutate(MSBPlong, vi = SEopt.upp*SEopt.upp) %>%
  filter(woodiness == "W"|woodiness == "H")

#Filter dataset for analysis (with phylogeny)
MSBP.meta.wood.phy <- mutate(MSBPlong, vi = SEopt.upp*SEopt.upp) %>%
                      filter(woodiness == "W"|woodiness == "H")  %>%
                      distinct(Taxon_ID, .keep_all = T) # needed for phlogentic analysis

#Get Phylogentic Matrix
MSBP.tree.cor.wood <- getPhyloCor(MSBP.meta.wood.phy)


## Current Climate Mismatch ~ Woodiness only (with phylogeny)
wood.mv.CM.phy <- rma.mv(yi = Current.Hot.Quart - Topt.upp, 
                         V = vi, 
                         mod = ~woodiness + log10(SeedAge) + altitude,
                         random = list(~1 | Taxon_ID, ~1 | grid.ll, ~1 | phylo),
                         R = list(phylo = MSBP.tree.cor.wood$cor),
                         method = "REML", 
                         data = MSBP.tree.cor.wood$MSBP.tree, verbose=F, digits=5, 
                         control=list(optimizer = "optim", optmethod = "Nelder-Mead", maxit= 10000))

## Current Climate Mismatch ~ Woodiness only (without phylogeny)
wood.mv.CM <- rma.mv(yi = Current.Hot.Quart - Topt.upp, 
                     V = vi, 
                     mod = ~woodiness + log10(SeedAge) + altitude,
                     random = list(~1 | Taxon_ID, ~1 | grid.ll),
                     method = "REML", 
                     data = MSBP.meta.wood)


## Current Climate Mismatch ~ Woodiness and lat (with phylogeny)
wood.mv.CM.phy.lat <- rma.mv(yi = Current.Hot.Quart - Topt.upp, 
                         V = vi, 
                         mod = ~AbsLat + woodiness + woodiness:AbsLat + log10(SeedAge) + altitude,
                         random = list(~1 | Taxon_ID, ~1 | grid.ll, ~1 | phylo),
                         R = list(phylo = MSBP.tree.cor.wood$cor),
                         method = "REML", 
                         data = MSBP.tree.cor.wood$MSBP.tree, verbose=F, digits=5, 
                         control=list(optimizer = "optim", optmethod = "Nelder-Mead", maxit= 10000))

## Current Climate Mismatch ~ Woodiness  and lat (without phylogeny)
wood.mv.CM.lat <- rma.mv(yi = Current.Hot.Quart - Topt.upp, 
                     V = vi, 
                     mod = ~AbsLat + woodiness + woodiness:AbsLat + log10(SeedAge) + altitude,
                     random = list(~1 | Taxon_ID, ~1 | grid.ll),
                     method = "REML", 
                     data = MSBP.meta.wood)


#Boxplot comparing woody and herbaceous species
woodbox <- ggplot(data = MSBP.meta.wood, aes(x = woodiness, y = Current.Hot.Quart - Topt.upp)) +
            geom_boxplot(fill = "springgreen1") + 
            theme_classic()+
            ylab(expression("Current Environment Temperature - \n Upper Optimal Germination Temperature ",(degree*C)))+
            xlab("") +
            scale_x_discrete(breaks=c("W", "H"),  labels=c("Woody", "Herbaceous"))+
            ylim(-35, 10) +
            theme(legend.position="none")

#Filter dataset for analysis (with phylogeny)
MSBP.meta.long <- mutate(MSBPlong, vi = SEopt.upp*SEopt.upp) %>%
                  filter(ann_per == "annual"|ann_per == "perennial")

#Filter dataset for analysis (without phylogeny)
MSBP.meta.long.phy <- mutate(MSBPlong, vi = SEopt.upp*SEopt.upp) %>%
                      filter(ann_per == "annual"|ann_per == "perennial")%>%
                      distinct(Taxon_ID, .keep_all = T) # needed for phlogentic analysis

#Get Phylogentic Matrix
MSBP.tree.cor <- getPhyloCor(MSBP.meta.long.phy)


## Current Climate Mismatch  ~ Longevity only (with phylogeny)
long.mv.CM.phy <- rma.mv(yi = Current.Hot.Quart - Topt.upp, 
                     V = vi, 
                     mod = ~ann_per + log10(SeedAge) + altitude,
                     random = list(~1 | Taxon_ID, ~1 | grid.ll, ~1 | phylo),
                     R = list(phylo = MSBP.tree.cor$cor),
                     method = "REML", 
                     data = MSBP.tree.cor$MSBP.tree, verbose=F, digits=5, 
                     control=list(optimizer = "optim", optmethod = "Nelder-Mead", maxit= 10000))

## Current Climate Mismatch ~ Longevity only (without phylogeny)
long.mv.CM <- rma.mv(yi = Current.Hot.Quart - Topt.upp, 
                     V = vi, 
                     mod = ~ann_per + log10(SeedAge) + altitude,
                     random = list(~1 | Taxon_ID, ~1 | grid.ll),
                     method = "REML", 
                     data = MSBP.meta.long)

## Current Climate Mismatch  ~ Longevity and AbsLat (with phylogeny)
long.mv.CM.phy.lat <- rma.mv(yi = Current.Hot.Quart - Topt.upp, 
                         V = vi, 
                         mod = ~AbsLat + ann_per + ann_per:AbsLat +log10(SeedAge) + altitude,
                         random = list(~1 | Taxon_ID, ~1 | grid.ll, ~1 | phylo),
                         R = list(phylo = MSBP.tree.cor$cor),
                         method = "REML", 
                         data = MSBP.tree.cor$MSBP.tree, verbose=F, digits=5, 
                         control=list(optimizer = "optim", optmethod = "Nelder-Mead", maxit= 10000))

## Current Climate Mismatch ~ Longevity and AbsLat (without phylogeny)
long.mv.CM.lat <- rma.mv(yi = Current.Hot.Quart - Topt.upp, 
                     V = vi, 
                     mod = ~AbsLat + ann_per + ann_per:AbsLat  + log10(SeedAge) + altitude,
                     random = list(~1 | Taxon_ID, ~1 | grid.ll),
                     method = "REML", 
                     data = MSBP.meta.long)


#Boxplot comparing woody and herbaceous species
longbox <- ggplot(data = MSBP.meta.long, aes(x = ann_per, y = Current.Hot.Quart - Topt.upp)) +
            geom_boxplot(fill = "springgreen1") + 
            theme_classic()+
            ylab(expression("Current Environment Temperature - \n Upper Optimal Germination Temperature ",(degree*C)))+
            xlab("") +
            ylim(-35, 10) +
            scale_x_discrete(breaks=c("annual", "perennial"),  labels=c("Annual", "Perennial"))+
            theme(legend.position="none")


output <- list(Wood.Model = wood.mv.CM, 
                 Wood.Model.Phy = wood.mv.CM.phy, 
                 Wood.Plot = woodbox,
                 Long.Model = long.mv.CM,
                 Long.Model.Phy = long.mv.CM.phy,
                 Long.Plot = longbox,
                 wood.mv.CM.phy.lat = wood.mv.CM.phy.lat,
                 wood.mv.CM.lat = wood.mv.CM.lat,
                 long.mv.CM.phy.lat = long.mv.CM.phy.lat,
                 long.mv.CM.lat = long.mv.CM.lat)

return(output)

}

CM.wood.life <- getLifeWoodCM()

