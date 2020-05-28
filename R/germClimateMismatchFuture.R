#' This is code to run analayses and create outputs for Future Climate Mismatch
#' Results in supplement
#' For Seed Germination over Latitude Project
#+ climate-mismatch-future, eval = FALSE

getCMF <- function() {
  
  require(dplyr)
  require(tidyr)
  require(metafor)
  require(ggplot2)
  source("./R/getPhyloCor.R")
  
  MSBP <- read.csv(file = "./Data/MSBP_Clean") #core dataset
  
  #Prepare data set, if SE's are too high, the analyses don't run (but get weighted out of the model anyway)
  MSBP.CM <- MSBP %>%
    mutate(vi = SEopt.upp * SEopt.upp) %>% #Standard error for opt
    filter(!is.na(Topt.upp) &
           !is.na(Future.Hot.Quart))
  
  MSBP.CM.phy <- distinct(MSBP.CM, Taxon_ID, .keep_all = T) #Needed for phylogeny analysis
  
  #Get phylogeny tree[[1]] and correlation matrix[[2]]
  MSBP.CM.phy <- getPhyloCor(MSBP.CM.phy)
  
  ## Current Climate Mismatch ~ Latitude and Hemispshere (without phylogeny)
  #using random effect of site (grid.ll) and species (Taxon_ID), hemisphere
  CM.h <- rma.mv(
    yi = Future.Hot.Quart - Topt.upp,
    V = vi,
    mod = ~ AbsLat + (AbsLat:NorthTF) + log10(SeedAge) + altitude,
    random = list(~ 1 | Taxon_ID, 
                  ~ 1 | grid.ll),
    method = "REML",
    data = MSBP.CM
  )
  
  ## Current Climate Mismatch ~ Latitude and Hemispshere (with phylogeny)
  #using random effect of site (grid.ll) and species (Taxon_ID), hemisphere
  CM.h.phy <- rma.mv(
    yi = Future.Hot.Quart - Topt.upp,
    V = vi,
    mod = ~ AbsLat + (AbsLat:NorthTF)+ log10(SeedAge) + altitude,
    random = list( ~ 1 | phylo, 
                   ~ 1 | Taxon_ID, 
                   ~ 1 | grid.ll),
    R = list(phylo = MSBP.CM.phy$cor),
    method = "REML",
    data = (MSBP.CM.phy[[1]]),
    verbose = TRUE,
    digits = 5,
    control = list(optimizer = "optim", optmethod = "Nelder-Mead")
  )
  
  ## Current Climate Mismatch ~ Latitude only (without phylogeny)
  #TSM using random effect of site (grid.ll) and species (Taxon_ID), no hemisphere
  CM <- rma.mv(
    yi = Future.Hot.Quart - Topt.upp,
    V = vi,
    mod = ~ AbsLat + log10(SeedAge) + altitude,
    random = list( ~ 1 | Taxon_ID, 
                   ~ 1 | grid.ll),
    method = "REML",
    data = MSBP.CM
  )
  
  CM.low <- rma.mv(
    yi = Future.Hot.Quart - Topt.low,
    V = vi,
    mod = ~ AbsLat+ log10(SeedAge) + altitude,
    random = list( ~ 1 | Taxon_ID, 
                   ~ 1 | grid.ll),
    method = "REML",
    data = MSBP.CM
  )
  
  ## Current Climate Mismatch ~ Latitude only (with phylogeny)
  #TSM using random effect of site (grid.ll) and species (Taxon_ID), no hemisphere
  CM.phy <- rma.mv(
    yi = Future.Hot.Quart - Topt.upp,
    V = vi,
    mod = ~ AbsLat + log10(SeedAge) + altitude,
    random = list( ~ 1 | phylo, 
                   ~ 1 | Taxon_ID, 
                   ~ 1 | grid.ll),
    R = list(phylo = MSBP.CM.phy$cor),
    method = "REML",
    data = MSBP.CM.phy[[1]]
  )
  
  CM.plot <- ggplot(data = MSBP.CM , aes(x = AbsLat)) +
    geom_rect(aes(
      xmin = 0, xmax = 23.5,
      ymax = Inf, ymin = -Inf
    ), fill = "grey90") +
    geom_point(aes(y = Future.Hot.Quart - Topt.upp, 
                   alpha = log(1 / vi)), colour = "orange") +
    geom_rug(aes(y = Future.Hot.Quart - Topt.upp),
             alpha = 1 / 2,
             position = "jitter",
             sides = "b") +
    geom_hline(aes(yintercept = 0), size = 0.5) +
    theme_classic() +
    ylab(expression(
      "Future Environmental Temperature - \n Upper Optimal Germination Temperature  " 
      (degree * C))) +
    xlab(expression("Absolute Latitude "(degree))) +
    theme(legend.position = "none",
          plot.margin = margin(10, 10, 10, 30)) +
    ylim (-35, 15) +
    scale_x_continuous(limits = c(0, 60), expand = c(0, 0))
  
  #Values to get confidence intervals
  mods <- cbind(AbsLat = c(1:57), 
                SeedAge = rep(log10(median(MSBP.CM$SeedAge, na.rm =T)), 57),
                altitude= rep(median(MSBP.CM$altitude, na.rm =T), 57)) #abslat - 1 to 57
  
  ### calculate predicted values from model
  preds1 <- predict(CM, newmods=mods, addx=T)
  
  ci.plot1 <- as.data.frame(
    cbind(c(seq(1, 57)), 
          preds1[[1]], 
          preds1[[3]], 
          preds1[[4]]))
  
  #Plot original figure with confidence intervals from above
  CM.plot.f.ci <-  CM.plot +
    geom_line(data = filter(ci.plot1, V1 > -55),
              aes(x = V1, y = V2),
              colour = "orange",
              linetype = "solid") +
    geom_line(data = filter(ci.plot1, V1 > -55),
              aes(x = V1, y = V3),
              linetype = "dashed") +
    geom_line(data = filter(ci.plot1, V1 > -55),
              aes(x = V1, y = V4),
              linetype = "dashed") +
    xlab("Latitude (°)")
  
  ggsave(
  CM.plot.f.ci,
  filename = "./Outputs/CVSE_future.tiff",
  units = "in",
  width = 7,
  height = 5,
  dpi = 1000,
  compression = 'lzw')

  
  output <- list(
    Model = CM,
    Model.Phy = CM.phy,
    Plot = CM.plot.f.ci,
    Model.h = CM.h,
    Model.h.Phy = CM.h.phy,
    Model.low = CM.low
  )
  
  return(output)
}

CMF <- getCMF()
