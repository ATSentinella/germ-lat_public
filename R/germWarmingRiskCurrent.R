#' Current Warming Risk analayses and outputs
#' For Seed Germination over Latitude Project
#+ warming-risk-current, eval = FALSE


getWRC <- function() {
  
  require(dplyr)
  require(tidyr)
  require(metafor)
  require(ggplot2)
  source("./R/getPhyloCor.R")
  
  MSBP <- read.csv(file = "./Data/MSBP_Clean") #core dataset
  
  ### Prepare data set, if SE's are too high, the analyses don't run
  # (but get weighted out of the model anyway)
  MSBP.WR.data <- MSBP %>%
    mutate(vi = SEmax * SEmax) %>% #Standard error for warming risk
    filter(!is.na(Tmax) & #Remove NAs
           !is.na(Current.Hot.Quart))  #Remove NAs 
  
  MSBP.WR.data.phy <- 
    MSBP.WR.data %>%
    distinct(Taxon_ID, .keep_all = T) #Needed for phylogeny analysis
  
  # Get phylogeny tree[[1]] and correlation matrix[[2]]
  MSBP.tree.cor <- getPhyloCor(MSBP.WR.data.phy)
  
  ## Future Warming Risk ~ Latitude and Hemispshere (with phylogeny)
  #WR, using random effect of phylo, site (grid.ll) and species (Taxon_ID), and altitude
  WR.h.phy <- rma.mv(
    yi = Current.Hot.Quart - Tmax,
    V = vi,
    mod = ~ AbsLat + (AbsLat:NorthTF) + 
      log10(SeedAge) + altitude,
    random = list( ~ 1 | phylo, ~ 1 | Taxon_ID, 
                   ~ 1 | grid.ll),
    R = list(phylo = MSBP.tree.cor$cor),
    method = "REML",
    data = MSBP.tree.cor$MSBP.tree,
    verbose = TRUE,
    digits = 5,
    control = list(optimizer = "optim", optmethod = "Nelder-Mead")
  )
  
  summary(WR.h.phy)
  
  
  WR.phy <- rma.mv(
    yi = Current.Hot.Quart - Tmax,
    V = vi,
    mod = ~ AbsLat+ log10(SeedAge) + altitude,
    random = list( ~ 1 | phylo, ~ 1 | Taxon_ID, 
                   ~ 1 | grid.ll),
    R = list(phylo = MSBP.tree.cor$cor),
    method = "REML",
    data = MSBP.tree.cor$MSBP.tree,
    verbose = TRUE,
    digits = 5,
    control = list(optimizer = "optim", optmethod = "Nelder-Mead")
  )
  summary(WR.phy)
  
  
  ## Future Warming Risk ~ Latitude and Hemispshere (with phylogeny)
  #WR, using random effect of phylo, site (grid.ll) and species (Taxon_ID), and altitude  
  WR.h <- rma.mv(
    yi = Current.Hot.Quart - Tmax,
    V = vi,
    mod = ~ AbsLat + (AbsLat:NorthTF) + 
      log10(SeedAge) + altitude,
    random = list(~ 1 | Taxon_ID, ~ 1 | grid.ll),
    method = "REML",
    data = MSBP.WR.data
  )
  
  summary(WR.h)
  
  
  WR <- rma.mv(
    yi = Current.Hot.Quart - Tmax,
    V = vi,
    mod = ~ AbsLat + log10(SeedAge) + altitude,
    random = list(~ 1 | Taxon_ID, ~ 1 | grid.ll),
    method = "REML",
    data = MSBP.WR.data
  )
  
  summary(WR)
  
  #Plot of future WR with lat and hemisphere effect
  
  WR.plot.c <- ggplot(data = MSBP.WR.data , aes(x = AbsLat)) +
    geom_rect(aes(
      xmin = 0,
      xmax = 23.5,
      ymax = Inf,
      ymin = -Inf
    ), fill = "grey90") +
    geom_rug(aes(y = Current.Hot.Quart - Tmax), alpha = 1/2, position = "jitter", sides="b") +
    geom_point(aes(y = Current.Hot.Quart - Tmax, 
                   alpha = log(1 / (vi))), colour = "Red") +
    geom_hline(aes(yintercept = 0), size = 0.5) +
    theme_classic() +
    ylab(expression(
      "2070 Environmental Temperature - \n Maximum Germination Temperature   "(degree * C)
    )) +
    xlab(expression("Absolute Latitude "(degree))) +
    theme(legend.position = "none",
          plot.margin = margin(10,10,10,30))+
    ylim (-45,15)+
    scale_x_continuous(limits = c(0,60), expand = c(0,0))
  
  WR.plot.c
  
  
  #Values to get confidence intervals
  mods <- cbind(AbsLat = c(1:57), 
                SeedAge = rep(log10(median(MSBP.WR.data$SeedAge, na.rm =T)), 57), 
                altitude= rep(median(MSBP.WR.data$altitude, na.rm =T), 57))
  
  ### calculate predicted values from model
  preds1 <- predict(WR, newmods=mods, addx=T)
  
  preds1
  
  ci.plot1 <-as.data.frame(cbind(c(seq(1, 57)), 
                                 preds1[[1]], 
                                 preds1[[3]], 
                                 preds1[[4]]))
  
  #Plot original figure with confidence intervals from above
  WR.plot.c.ci <-
    WR.plot.c +
    geom_line(data = filter(ci.plot1, V1 > -55),
              aes(x = V1, y = V2),
              colour = "red",
              linetype = "solid") +
    geom_line(data = filter(ci.plot1, V1 > -55),
              aes(x = V1, y = V3),
              linetype = "dashed") +
    geom_line(data = filter(ci.plot1, V1 > -55),
              aes(x = V1, y = V4),
              linetype = "dashed") +
    xlab("Latitude (°)")
  
  WR.plot.c.ci
  
  ggsave(
    WR.plot.c.ci,
    filename = "./Outputs/WRSE_current.tiff",
    units = "in",
    width = 7,
    height = 5,
    dpi = 1000,
    compression = 'lzw')
  
  #Return results
  output <- list(
    Model = WR,
    Model.phy = WR.phy,
    Plot = WR.plot.c.ci,
    Model.h = WR.h,
    Model.h.phy = WR.h.phy
  )
  
  return(output)
  
}


WRC <- getWRC()
