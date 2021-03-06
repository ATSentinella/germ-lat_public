#' This is code to run analayses and create outputs for Current Climate Mismatch
#+ climate-mismatch-current, eval = FALSE

getCM <- function(fromfile = T) {

  if (fromfile == T) {
    output <- readRDS("./Outputs/CM_current")
  } else {
  
  
require(dplyr)
require(tidyr)
require(metafor)
require(ggplot2)
source("./R/getPhyloCor.R")

MSBP <- read.csv(file = "./Data/MSBP_Clean") #core dataset

#Prepare data set
MSBP.CM <- MSBP %>%
  mutate(vi = SEopt.upp * SEopt.upp) %>% #Standard error for Topt.upp
  filter(!is.na(Topt.upp) &
         !is.na(Current.Hot.Quart) &
         !is.na(SeedAge) &
         !is.na(Altitude)) #remove NAs

MSBP.CM.phy <- distinct(MSBP.CM, Taxon_ID, .keep_all = T) #Needed for phylogeny analysis

#Get phylogeny tree[[1]] and correlation matrix[[2]]
MSBP.CM.phy <- getPhyloCor(MSBP.CM.phy)

## Current Climate Mismatch ~ Latitude and Hemispshere (without phylogeny)
# additional fixed effects of log10(SeedAge) and altitude
# using random effect of site (grid.ll) and species (Taxon_ID)
CM.h <- rma.mv(
  yi = Current.Hot.Quart - Topt.upp,
  V = vi,
  mod = ~ AbsLat + (AbsLat:NorthTF) + log10(SeedAge) + altitude,
  random = list(~ 1 | Taxon_ID, 
                ~ 1 | grid.ll),
  method = "REML",
  data = MSBP.CM
)
CM.h #no effect of hemisphere, so repeat without

## Current Climate Mismatch ~ Latitude and Hemispshere
# additional fixed effects of log10(SeedAge) and altitude
# using random effect of site (grid.ll) and species (Taxon_ID)
# and phylogeny with corresponding correlation matrix
CM.h.phy <- rma.mv(
  yi = Current.Hot.Quart - Topt.upp,
  V = vi,
  mod = ~ AbsLat + (AbsLat:NorthTF) + log10(SeedAge) + altitude,
  random = list( ~ 1 | phylo, 
                 ~ 1 | Taxon_ID, 
                 ~ 1 | grid.ll),
  R = list(phylo = MSBP.CM.phy$cor),
  method = "REML",
  data = (MSBP.CM.phy[[1]]),
  verbose = FALSE,
  digits = 5,
  control = list(optimizer = "optim", optmethod = "Nelder-Mead")
)
CM.h.phy #no effect of hemisphere, so repeat without

## Current Climate Mismatch ~ Latitude
# additional fixed effects of log10(SeedAge) and altitude
# using random effect of site (grid.ll) and species (Taxon_ID)
# and phylogeny with corresponding correlation matrix
CM.phy <- rma.mv(
  yi = Current.Hot.Quart - Topt.upp,
  V = vi,
  mod = ~ AbsLat + log10(SeedAge) + altitude,
  random = list( ~ 1 | phylo, 
                 ~ 1 | Taxon_ID, 
                 ~ 1 | grid.ll),
  R = list(phylo = MSBP.CM.phy$cor),
  method = "REML",
  data = MSBP.CM.phy[[1]]
)
CM.phy #no interactions with abslat 

## Current Climate Mismatch ~ Latitude and Hemispshere
# additional fixed effects of log10(SeedAge) and altitude
# using random effect of site (grid.ll) and species (Taxon_ID)
CM <- rma.mv(
  yi = Current.Hot.Quart - Topt.upp,
  V = vi,
  mod = ~ AbsLat + log10(SeedAge) + altitude,
  random = list( ~ 1 | Taxon_ID, 
                 ~ 1 | grid.ll),
  method = "REML",
  data = MSBP.CM
)


CM.low <- rma.mv(
  yi = Current.Hot.Quart - Topt.low,
  V = vi,
  mod = ~ AbsLat + log10(SeedAge) + altitude,
  random = list( ~ 1 | Taxon_ID, 
                 ~ 1 | grid.ll),
  method = "REML",
  data = MSBP.CM
)

CM #no interactions with abslat 

CM.plot <- ggplot(data = MSBP.CM , aes(x = AbsLat)) +
  geom_rect(aes(
    xmin = 0, xmax = 23.5,
    ymax = Inf, ymin = -Inf
  ), fill = "grey90") +
    geom_hline(aes(yintercept = 0), size = 0.5) +
  geom_point(aes(y = Current.Hot.Quart - Topt.upp, 
                 alpha = sqrt(1 / vi)), colour = "orange") +
  geom_rug(aes(y = Current.Hot.Quart - Topt.upp),
    alpha = 1 / 2,
    position = "jitter",
    sides = "b") +
  theme_classic() +
  ylab(expression(
    "Current Environmental Temperature - \n Upper Optimal Germination Temperature  " 
    (degree * C))) +
  xlab(expression("Absolute Latitude "(degree))) +
  theme(legend.position = "none",
        plot.margin = margin(10, 10, 10, 30)) +
    scale_alpha_continuous(range = c(0.4,1))+
  scale_x_continuous(limits = c(0, 60), expand = c(0, 0))


 CM.plot 

 
#Values to get confidence intervals
mods <- cbind(AbsLat = c(1:57), 
              SeedAge= rep(log10(1479), 57),
              altitude= rep(393, 57)) #abslat - 1 to 57

### calculate predicted values from model
preds1 <- predict(CM, newmods=mods, addx=T)

ci.plot1 <- as.data.frame(
  cbind(c(seq(1, 57)), 
        preds1[[1]],
        preds1[[3]], 
        preds1[[4]])) %>%
  filter(V1 > 0)

#Plot original figure with confidence intervals from above
CM.plot <-  CM.plot +
  geom_line(data = ci.plot1,
            aes(x = V1, y = (V2)),
            colour = "orange",
            linetype = "solid") +
  geom_line(data = ci.plot1,
            aes(x = V1, y = (V3)),
            linetype = "dashed") +
  geom_line(data = ci.plot1,
            aes(x = V1, y = (V4)),
            linetype = "dashed") +
  xlab("Latitude (°)")

CM.plot

ggsave(
  CM.plot,
  filename = "./Outputs/CVSE.tiff",
  units = "in",
  width = 4,
  height = 5,
  dpi = 1000,
  compression = 'lzw'
)


output <- list(
  Model = CM,
  Model.phy = CM.phy,
  Plot = CM.plot,
  Model.h = CM.h,
  Model.h.phy = CM.h.phy,
  Model.low = CM.low
)

saveRDS(output, "./Outputs/CM_current")
  }
  
  return(output)
  
}

