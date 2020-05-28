#' This is code to run analayses and create outputs for Current Climate Vulnerability
#' For Seed Germination over Latitude Project
#' 
#' @param 
#' @return 
#' @author Alex Sentinella ATSentinella@gmail.com


getLifeWoodWR <- function() {
  
  require(dplyr)
  require(tidyr)
  require(metafor)
  require(ggbeeswarm)
  source("./R/getPhyloCor.R")


longevity <- read.csv(file = "./Data/Longevity/TRY.longevity.csv")

MSBP <- read.csv(file = "./Data/MSBP_Clean") #core dataset

#Prepare data set, if SE's are too high, the analyses don't run (but get weighted out anyway)
MSBP <- MSBP %>%
  mutate(vi = SEmax*SEmax) %>% #Standard error for warming risk
  filter(Tmax > 0 & SEmax<200) %>% #Filter where Tmax is above 0 and SEmax is lower than 200
  distinct(Taxon_ID, .keep_all = T) #Is this needed?

longevityfilt <- longevity %>% 
  mutate(ann_per = case_when(OrigValueStr == "annual" | 
                               OrigValueStr == "Annual" | 
                               OrigValueStr == "annuals" |
                               OrigValueStr == "summer annuals" |
                               OrigValueStr == "always annual" |
                               OrigValueStr == "winter annuals" |
                               OrigValueStr == "annual-winter annual" |
                               OrigValueStr == "winter annual" |
                               (OriglName == "Life history" & OrigValueStr == "1" ) |
                               (OriglName == "Plant phenology: Annual" & OrigValueStr == "yes" )     ~   "annual",  ######annuals
                             OrigValueStr == "perennial" | 
                               OrigValueStr == "Perennial" | 
                               OrigValueStr == "perennials" |                                   
                               OrigValueStr == "always pluriennial-pollakanthic" | 
                               (OriglName == "Plant phenology: Biennial" & OrigValueStr == "yes" ) | 
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
                               (OriglName == "Plant phenology: Perennial" & OrigValueStr == "yes" )| 
                               (OriglName == "Plant phenology: Annual" & OrigValueStr == "no" )    ~   "perennial"  ###perennials
  ))

longevityfilt <- data.frame(lapply(longevityfilt, function(x) {gsub(" ", "_", x)}))

MSBPlong <- left_join(MSBP, longevityfilt, by = c("Taxon_ID" ="SpeciesName"))

woody <- read.csv(file = "./Data/GlobalWoodinessDatabase.csv")

woody <- data.frame(lapply(woody, function(x) {gsub(" ", "_", x)}))

MSBPlong <- left_join(MSBPlong, woody, by = c("Taxon_ID" ="gs"))

count(MSBPlong, ann_per)

MSBPlong <- MSBPlong %>%
  mutate(ann_per=replace(ann_per, woodiness=="W", "perennial"))

count(MSBPlong, woodiness)



MSBP.meta.wood.phy <- mutate(MSBPlong, vi = SEmax*SEmax) %>%
   filter(woodiness == "W"|woodiness == "H")%>%
                  distinct(Taxon_ID, .keep_all = T)

MSBP.meta.wood <- mutate(MSBPlong, vi = SEmax*SEmax) %>%
                  filter(woodiness == "W"|woodiness == "H")



MSBP.tree.cor.wood <- getPhyloCor(MSBP.meta.wood.phy)

wood.mv.WT.phy <- rma.mv(yi = Future.Temp - Tmax, 
                         V = vi, 
                         mod = ~woodiness,
                         random = list(~1 | Taxon_ID, ~1 | grid.ll, ~1 | phylo),
                         R = list(phylo = MSBP.tree.cor.wood$cor),
                         method = "REML", 
                         data = (MSBP.tree.cor.wood$MSBP.tree), verbose=F, digits=5, 
                         control=list(optimizer = "optim", optmethod = "Nelder-Mead", maxit= 10000))
summary(wood.mv.WT.phy)


wood.mv.WT <- rma.mv(yi = Future.Temp - Tmax, 
                     V = vi, 
                     mod = ~woodiness,
                     random = list(~1 | Taxon_ID, ~1 | grid.ll),
                     method = "REML", 
                     data = (MSBP.meta.wood))
summary(wood.mv.WT)


woodbox <- ggplot(data = MSBP.meta.wood, aes(x = woodiness, y = Future.Temp - Tmax)) +
  geom_boxplot(aes(fill = woodiness)) + 
  geom_quasirandom()+
  scale_fill_brewer(palette="Accent")+
  theme_classic()+
  ylab("Predicted 2070 Environment Temperature - \n Maximum Germination Temerature (°C)")+
  xlab("") +
  scale_x_discrete(breaks=c("W", "H"),  labels=c("W", "H"))+
  theme(legend.position="none")


MSBP.meta.long <- mutate(MSBPlong, vi = SEmax*SEmax) %>%
  filter(ann_per == "annual"|ann_per == "perennial")

MSBP.meta.long.phy <- mutate(MSBPlong, vi = SEmax*SEmax) %>%
  filter(ann_per == "annual"|ann_per == "perennial")%>%
  distinct(Taxon_ID, .keep_all = T)

MSBP.tree.cor <- getPhyloCor(MSBP.meta.long.phy)



long.mv.WT.phy <- rma.mv(yi = Future.Temp - Tmax, 
                     V = vi, 
                     mod = ~ann_per,
                     random = list(~1 | Taxon_ID, ~1 | grid.ll, ~1 | phylo),
                     R = list(phylo = MSBP.tree.cor$cor),
                     method = "REML", 
                     data = (MSBP.tree.cor$MSBP.tree), verbose=F, digits=5, 
                     control=list(optimizer = "optim", optmethod = "Nelder-Mead", maxit= 10000))
summary(long.mv.WT.phy)


long.mv.WT <- rma.mv(yi = Future.Temp - Tmax, 
                     V = vi, 
                     mod = ~ann_per,
                     random = list(~1 | Taxon_ID, ~1 | grid.ll),
                     method = "REML", 
                     data = MSBP.meta.long)
summary(long.mv.WT)



longbox <- ggplot(data = MSBP.meta.long, aes(x = ann_per, y = Future.Temp - Tmax)) +
  geom_boxplot(aes(fill = ann_per)) + 
  geom_quasirandom()+
  scale_fill_brewer(palette="Dark2")+
  theme_classic()+
  ylab("Predicted 2070 Environment Temperature - \n Maximum Germination Temerature (°C)")+
  xlab("") +
  scale_x_discrete(breaks=c("annual", "perennial"),  labels=c("Annual", "Perennial"))+
  theme(legend.position="none")

WoodLong <- list(Wood.Model = wood.mv.WT, 
                 Wood.Model.Phy = wood.mv.WT.phy, 
                 Wood.Plot = woodbox,
                 Long.Model = long.mv.WT,
                 Long.Model.Phy = long.mv.WT.phy,
                 Long.Plot = longbox)

return(WoodLong)

}



