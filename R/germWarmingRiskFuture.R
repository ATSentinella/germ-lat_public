#' Future Warming Risk analayses and outputs
#' For Seed Germination over Latitude Project
#+ warming-risk-future, eval = FALSE


getWR <- function(fromfile = T) {
  
  if (fromfile == T) {
    output <- readRDS("./Outputs/WR_future")
  } else {
    
  
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
    filter(!is.na(Tmax) & 
           !is.na(Future.Hot.Quart) &
           !is.na(SeedAge) &
           !is.na(Altitude)) #remove NAs
  
  MSBP.WR.data.phy <- distinct(MSBP.WR.data, Taxon_ID, .keep_all = T) #Needed for phylogeny analysis
  
  # Get phylogeny tree[[1]] and correlation matrix[[2]]
  MSBP.tree.cor <- getPhyloCor(MSBP.WR.data.phy)

  ## Future Warming Risk ~ Latitude and Hemispshere (with phylogeny)
  #WR, using random effect of phylo, site (grid.ll) and species (Taxon_ID), and altitude
  WR.h.phy <- rma.mv(
    yi = Future.Hot.Quart - Tmax,
    V = vi,
    mod = ~ AbsLat + (AbsLat:NorthTF) + log10(SeedAge) + altitude,
    random = list( ~ 1 | phylo, ~ 1 | Taxon_ID, 
                   ~ 1 | grid.ll),
    R = list(phylo = MSBP.tree.cor$cor),
    method = "REML",
    data = MSBP.tree.cor$MSBP.tree,
    verbose = F,
    digits = 5,
    control = list(optimizer = "optim", optmethod = "Nelder-Mead")
  )
  
  summary(WR.h.phy) #Hemisphere effect (p = 0.02), no interaction effects 
  
  ## Future Warming Risk ~ Latitude and Hemispshere (with phylogeny)
  #WR, using random effect of phylo, site (grid.ll) and species (Taxon_ID), and altitude  
  WR.h <- rma.mv(
    yi = Future.Hot.Quart - Tmax,
    V = vi,
    mod = ~ AbsLat + (AbsLat:NorthTF) + log10(SeedAge) + altitude,
    random = list(~ 1 | Taxon_ID, ~ 1 | grid.ll),
    method = "REML",
    data = MSBP.WR.data
  )
  
  summary(WR.h) #Hemisphere effect (p = 0.02)


#Plot of future WR with lat and hemisphere effect, no interaction effects 
  
WR.plot.f <- ggplot(data = MSBP.WR.data , aes(x = Latitude)) +
    geom_rect(aes(
      xmin = -23.5,
      xmax = 23.5,
      ymax = Inf,
      ymin = -Inf ), fill = "grey90") +
    geom_rug(aes(y = Future.Hot.Quart - Tmax), alpha = 1/2, position = "jitter", sides="b") +
    geom_point(aes(y = Future.Hot.Quart - Tmax, 
                   alpha = sqrt(1/SEmax)), colour = "Red") +
    geom_hline(aes(yintercept = 0), size = 0.5) +
    theme_classic() +
    ylab(expression(
      "Future Environmental Temperature - \n Maximum Germination Temperature   "(degree * C)
    )) +
    xlab(expression("Absolute Latitude "(degree))) +
  scale_alpha_continuous(range = c(0.4,1))+
    theme(legend.position = "none", plot.margin = margin(10,10,10,30))+
     scale_x_continuous(limits = c(-60,64), expand = c(0,0))
  
  WR.plot.f

  #Values to get confidence intervals
  mods <- cbind(AbsLat = c(0:70, 0:70), 
                SeedAge= rep(log10(1479), 142),
                altitude= rep(393, 142),
                NorthTF= c(rep(T, 71), rep(F, 71)))
  
  ### calculate predicted values from model
  preds1 <- predict(WR.h, newmods=mods, addx=T)
  
  ci.plot1 <- as.data.frame(
    cbind(c(0:70, 0:(-70)), 
          preds1[[1]],
          preds1[[3]], 
          preds1[[4]]))
  
  #Plot original figure with confidence intervals from above
  WR.plot.f <-  WR.plot.f +
    geom_line(data = filter(ci.plot1, V1 > -53 & V1 < 1),
              aes(x = V1, y = (V2)),
              colour = "red",
              linetype = "solid") + #South
    geom_line(data = filter(ci.plot1, V1 > -53 & V1 < 1),
              aes(x = V1, y = (V3)),
              linetype = "dashed") + #South
    geom_line(data = filter(ci.plot1, V1 > -53 & V1 < 1),
              aes(x = V1, y = (V4)),
              linetype = "dashed") + #South
    geom_line(data = filter(ci.plot1, V1 > -1 & V1 < 64),
              aes(x = V1, y = (V2)),
              colour = "red",
              linetype = "solid") + #North
    geom_line(data = filter(ci.plot1, V1 > -1 & V1 < 64),
              aes(x = V1, y = (V3)),
              linetype = "dashed") + #North
    geom_line(data = filter(ci.plot1, V1 > -1 & V1 < 64),
              aes(x = V1, y = (V4)),
              linetype = "dashed") + #North
    xlab("Latitude (°)")
  
  WR.plot.f

ggsave(
  WR.plot.f,
  filename = "./Outputs/WRSE.tiff",
  units = "in",
  width = 4,
  height = 5,
  dpi = 1000,
  compression = 'lzw')

  
#Return results
output <- list(
    Model.h = WR.h,
    Model.h.phy = WR.h.phy,
    Plot = WR.plot.f
  )
  
saveRDS(output, "./Outputs/WR_future")
  }
  
  return(output)
  
}



