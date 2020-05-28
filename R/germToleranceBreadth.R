#' Germination Temperature Breadth Over Latitude analyses and outputs
#' For Seed Germination over Latitude Project
#'
#' fromfile = TRUE - Use previous analysis file
#' fromfile = FALSE - Rerun analysis
#' 
#+ get-tolerance-breath, eval = FALSE


getToleranceBreadth <- function(fromfile = T) {
  if (fromfile == T) {
    #switch to not rerun analysis
    output <- readRDS("./Outputs/TB")
  } else {
    #Load packages
    require(dplyr)
    require(tidyr)
    require(metafor)
    require(ggplot2)
    source("./R/getPhyloCor.R")
    
    #Load data
    MSBP <- read.csv(file = "./Data/MSBP_Clean")
    
    #Prepare data set
    MSBP.meta.TB <- MSBP %>%
      mutate(vi = SEbreadth * SEbreadth) %>% #Variances
      filter(!is.na(Tbreadth) &
               !is.na(SeedAge) &
               !is.na(Altitude)) #Remove NA values
    
    #Remove duplicate species, for phylogeny analysis
    MSBP.meta.TB.phy <-
      distinct(MSBP.meta.TB, Taxon_ID, .keep_all = T)
    
    #Get phylogeny: tree - .[[1]] and correlation matrix - .[[2]]
    MSBP.tree.cor <- getPhyloCor(MSBP.meta.TB.phy)
    
    ### Models below:
    # breadth - Tbreadth ~ AbsLat
    # .h - with fixed effect for hemipshere
    # .phy - with random effect for phylogeny, using variance covariance matrix
    #
    # All have species and site as random factors
    # AND log10(SeedAge) and altitude as fixed effects
    
    breadth.h.phy <-  rma.mv(
      yi = TB,
      V = vi,
      mod = ~ AbsLat + AbsLat:NorthTF + log10(SeedAge) + altitude,
      random = list(~ 1 | phylo,
                    ~ 1 | Taxon_ID,
                    ~ 1 | grid.ll),
      R = list(phylo = MSBP.tree.cor$cor),
      method = "REML",
      data = MSBP.tree.cor$MSBP.tree,
      verbose = F,
      digits = 5,
      control = list(optimizer = "optim", optmethod = "Nelder-Mead")
    )
    breadth.h.phy #no hemisphere effect, therefore remove term
    
    breadth.h <- rma.mv(
      yi = TB,
      V = vi,
      mod = ~ AbsLat + AbsLat:NorthTF + log10(SeedAge) + altitude,
      random = list( ~ 1 | Taxon_ID,
                     ~ 1 | grid.ll),
      method = "REML",
      data = MSBP.meta.TB,
      control = list(optimizer = "optim", optmethod = "Nelder-Mead")
    )
    breadth.h #no hemisphere effect, therefore remove term
    
    
    breadth.phy <-  rma.mv(
      yi = TB,
      V = vi,
      mod = ~ AbsLat + log10(SeedAge) + altitude,
      random = list(~ 1 | phylo,
                    ~ 1 | Taxon_ID,
                    ~ 1 | grid.ll),
      R = list(phylo = MSBP.tree.cor$cor),
      method = "REML",
      data = MSBP.tree.cor$MSBP.tree,
      verbose = F,
      digits = 5,
      control = list(optimizer = "optim", optmethod = "Nelder-Mead")
    )
    breadth.phy
    
    
    breadth <- rma.mv(
      yi = TB,
      V = vi,
      mod = ~ AbsLat + log10(SeedAge)  + altitude,
      random = list( ~ 1 | Taxon_ID,
                     ~ 1 | grid.ll),
      method = "REML",
      data = MSBP.meta.TB,
      control = list(optimizer = "optim", optmethod = "Nelder-Mead")
    )
    breadth
    
    
    #Plot of ranges against absolute latitude with lat effect
    TB_plot <- ggplot(data = MSBP.meta.TB, aes(x = AbsLat)) +
      geom_rect(aes(
        xmin = 0,
        xmax = 23.5,
        ymax = Inf,
        ymin = -Inf
      ), fill = "grey90") +
      geom_point(aes(y = TB, alpha = sqrt(1 / SEbreadth)), colour = "green4") +
      geom_rug(
        aes(y = TB),
        alpha = 1 / 2,
        position = "jitter",
        sides = "b"
      ) +
      theme_classic() +
      ylim(0, 50) +
      ylab(
        expression(
"Maximum Germination Temperature - \n Minimum Germination Temperature  "(degree * C)
        )
      ) +
      xlab(expression("Absolute Latitude "(degree))) +
      theme(legend.position = "none",
            plot.margin = margin(10, 10, 10, 30)) +
      scale_alpha_continuous(range = c(0.5, 1)) +
      scale_x_continuous(limits = c(0, 60), expand = c(0, 0))
    
    TB_plot
    
    #Values to get confidence intervals of model fit
    mods <- cbind(
      AbsLat = c(1:57),
      #abslat - 1 to 57
      SeedAge = rep(log10(1437), 57),
      #median of SeedAge, 57 times
      altitude = rep(435, 57)
    ) #median of altitude, 57 times
    
    ### calculate predicted values from model
    preds1 <- predict.rma(breadth, newmods = mods, addx = T)
    
    ci.plot1 <- as.data.frame(cbind(c(seq(1, 57)), #1 to 57
                                    preds1[[1]],  #Model fit
                                    preds1[[3]],  #Confidence interval of model
                                    preds1[[4]])) #Confidence interval of model
    
    #Plot original figure with confidence intervals from above
    TB_plot <-  TB_plot +
      geom_line(data = filter(ci.plot1, V1 > -54 & V1 > 0),
                aes(x = V1, y = V3),
                linetype = "dashed") +
      geom_line(data = filter(ci.plot1, V1 > -54 & V1 > 0),
                aes(x = V1, y = V4),
                linetype = "dashed") +
      xlab("Latitude (°)")
    
    TB_plot
    
    ggsave(
      file = "./Outputs/LinesAbsLatRangeSE.tiff",
      units = "in",
      width = 5,
      height = 4,
      dpi = 1000,
      compression = 'lzw'
    )
    
    
    output <- list(
      Model = breadth,
      Model.phy = breadth.phy,
      Plot = TB_plot,
      Model.h = breadth.h,
      Model.h.phy = breadth.h.phy
    )
    
    saveRDS(output, "./Outputs/TB")
  }
  
  return(output)
  
}
