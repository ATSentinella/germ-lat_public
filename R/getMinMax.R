#'Maximum and Minimum Germiantion temepratures Over Latitude analyses and outputs
#' For Seed Germination over Latitude Project
#' 
#' @author Alex Sentinella ATSentinella@@gmail.com
#' @return TmaxTmin (Model, Model.Phy, Plot, pseudoR2, Model.h, Model.h.Phy)


getMinMax <- function() {
  
  require(dplyr)
  require(tidyr)
  require(metafor)
  require(ggplot2)
  source("./R/getPhyloCor.R")
  
  
  library(summarytools)
  view(dfSummary(iris))

  
  MSBP <- read.csv(file = "./Data/MSBP_Clean") 
  
  #Prepare data set, if SE's are too high, the analyses don't run (but get weighted out of the model anyway)
  MSBP.meta.min <- MSBP %>%
    mutate(vi = SEmin*SEmin) %>% #Standard error for breadth
    filter(!is.na(Tmin) & SEmin<200) #Filter where Tmax is above 0 and SEbreadth is lower than 200
  
  
  
  
  MSBP.meta.min.phy <- distinct(MSBP.meta.min, Taxon_ID, .keep_all = T) #Needed for phylogeny analysis
  
  #Get phylogeny tree[[1]] and correlation matrix[[2]]
  MSBP.tree.cor.min <- getPhyloCor(MSBP.meta.min.phy)
  
  ## Tmin ~ Latitude and Hemispshere (with phylogeny)
  #Using random effect of site (grid.ll) and species (Taxon_ID), hemisphere
  MinLat.h.phy <- rma.mv(yi = Tmin, 
                            V = vi, 
                            mod = ~AbsLat + (AbsLat:NorthTF) + altitude, 
                            random = list(~1 | phylo, ~1 | Taxon_ID, 
                                          ~1 | grid.ll),
                            R = list(phylo = MSBP.tree.cor.min$cor),
                            method = "REML", 
                            data = MSBP.tree.cor.min$MSBP.tree, verbose=TRUE, digits=5, 
                            control=list(optimizer = "optim", optmethod = "Nelder-Mead"))
  
  summary(MinLat.h.phy)
  
  ## Tmin ~ Latitude and Hemispshere (without phylogeny)
  MinLat.h <- rma.mv(yi = Tmin, 
                        V = vi, 
                        mod = ~AbsLat + (AbsLat:NorthTF)+ altitude, 
                        random = list(~1 | Taxon_ID, ~1 | grid.ll),
                        method = "REML", 
                        data = MSBP.meta.min,
                        control=list(optimizer = "optim", optmethod = "Nelder-Mead"))
  
  
  summary(MinLat.h) #hemi not significant
  
  ## Tmin ~ Latitude only (with phylogeny)
  #using random effect of site (grid.ll) and species (Taxon_ID), no hemisphere
  MinLat.phy <- rma.mv(yi = Tmin, 
                          V = vi, 
                          mod = ~AbsLat+ altitude,
                          random = list(
                            ~1 | phylo, 
                            ~1 | Taxon_ID, 
                            ~1 | grid.ll),
                          R = list(phylo = MSBP.tree.cor.min$cor),
                          method = "REML", 
                          data = MSBP.tree.cor.min$MSBP.tree, verbose=TRUE, digits=5, 
                          control=list(optimizer = "optim", optmethod = "Nelder-Mead"))
  
  summary(MinLat.phy) 
# Model Results:
#   
#   estimate       se      zval     pval     ci.lb     ci.ub     
#   intrcpt   17.27672  4.57322   3.77780  0.00016   8.31337  26.24007  ***
#   AbsLat    -0.10795  0.02368  -4.55910  <.00001  -0.15436  -0.06154  ***
#   altitude  -0.00070  0.00055  -1.27825  0.20116  -0.00177   0.00037     
  
  
  ## Tmin ~ Latitude only (without phylogeny)
  MinLat <- rma.mv(yi = Tmin, 
                      V = vi, 
                      mod = ~AbsLat + altitude, 
                      random = list(
                        ~1 | Taxon_ID, 
                        ~1 | grid.ll),
                      method = "REML", 
                      data = MSBP.meta.min,
                      control=list(optimizer = "optim", optmethod = "Nelder-Mead"))
  
  summary(MinLat)
#  Model Results:
#    
#    estimate      se     zval    pval    ci.lb    ci.ub     
#    intrcpt    20.1889  0.6750  29.9073  <.0001  18.8658  21.5119  ***
#    AbsLat     -0.1539  0.0205  -7.5089  <.0001  -0.1940  -0.1137  ***
#    altitude   -0.0005  0.0005  -0.9338  0.3504  -0.0016   0.0005     

  #Plot of ranges against absolute latitude with lat effect, not significant R^2 0.73%
  Min.plot <- ggplot(data = MSBP.meta.min , aes(x = AbsLat)) +
    geom_rect(aes(xmin=0, xmax = 23.5, ymax =Inf, ymin = -Inf), fill = "grey90")+
    geom_point(aes(y = Tmin, alpha = log(1/vi)), colour = "blue") +   
    theme_classic() +
    ylab(expression("Germination Temperature Breadth "(degree*C)))+
    xlab(expression("Absolute Latitude "(degree))) +
    geom_segment(x = 1, y = coef(MinLat)[1]+ (coef(MinLat)[2])*1, 
                 xend = 58, yend = (coef(MinLat)[1]) + (coef(MinLat)[2])*58, 
                 size = 1, colour = "blue", linetype = "solid")+
    #coord_flip() +
    theme(legend.position="none", plot.margin = margin(10,10,10,30))+
    scale_y_continuous(limits = c(-15,50), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,60), expand = c(0,0))
  
  Min.plot
  
  ggsave(file = "./Outputs/TminAbsLat.tiff", units="in", width=5, height=4, dpi=1000, compression = 'lzw')
  
  
  #Prepare data set, if SE's are too high, the analyses don't run (but get weighted out of the model anyway)
  MSBP.meta.max <- MSBP %>%
    mutate(vi = SEmax*SEmax) %>% #Standard error for breadth
    filter(!is.na(Tmax) & SEmax<200) #Filter where Tmax is above 0 and SEbreadth is lower than 200
  
  
  
  
  MSBP.meta.max.phy <- distinct(MSBP.meta.max, Taxon_ID, .keep_all = T) #Needed for phylogeny analysis
  
  #Get phylogeny tree[[1]] and correlation matrix[[2]]
  MSBP.tree.cor.max <- getPhyloCor(MSBP.meta.max.phy)
  
  ## Tmax ~ Latitude and Hemispshere (with phylogeny)
  #Using random effect of site (grid.ll) and species (Taxon_ID), hemisphere
  MaxLat.h.phy <- rma.mv(yi = Tmax, 
                         V = vi, 
                         mod = ~AbsLat + (AbsLat:NorthTF) +altitude, 
                         random = list(~1 | phylo, ~1 | Taxon_ID, 
                                       ~1 | grid.ll),
                         R = list(phylo = MSBP.tree.cor.max$cor),
                         method = "REML", 
                         data = MSBP.tree.cor.max$MSBP.tree, verbose=TRUE, digits=5, 
                         control=list(optimizer = "optim", optmethod = "Nelder-Mead"))
  
  summary(MaxLat.h.phy) #hemi not sig
  
  ## Tmax ~ Latitude and Hemispshere (without phylogeny)
  MaxLat.h <- rma.mv(yi = Tmax, 
                     V = vi, 
                     mod = ~AbsLat + (AbsLat:NorthTF) +altitude, 
                     random = list(~1 | Taxon_ID, ~1 | grid.ll),
                     method = "REML", 
                     data = MSBP.meta.max,
                     control=list(optimizer = "optim", optmethod = "Nelder-Mead"))
  
  
  summary(MinLat.h) #hemi not significant
  
  ## Tmax ~ Latitude only (with phylogeny)
  #using random effect of site (grid.ll) and species (Taxon_ID), no hemisphere
  MaxLat.phy <- rma.mv(yi = Tmax, 
                       V = vi, 
                       mod = ~AbsLat +altitude,
                       random = list(
                         ~1 | phylo, 
                         ~1 | Taxon_ID, 
                         ~1 | grid.ll),
                       R = list(phylo = MSBP.tree.cor.max$cor),
                       method = "REML", 
                       data = MSBP.tree.cor.max$MSBP.tree, verbose=TRUE, digits=5, 
                       control=list(optimizer = "optim", optmethod = "Nelder-Mead"))
  
  summary(MaxLat.phy) 
#  Model Results:
#    
#    estimate       se      zval     pval     ci.lb     ci.ub     
#    intrcpt   30.63323  3.46408   8.84310  <.00001  23.84375  37.42270  ***
#    AbsLat    -0.08335  0.03130  -2.66266  0.00775  -0.14470  -0.02200   **
#    altitude  -0.00243  0.00067  -3.64134  0.00027  -0.00373  -0.00112  ***
    
  
  ## Tmax ~ Latitude only (without phylogeny)
  MaxLat <- rma.mv(yi = Tmax, 
                   V = vi, 
                   mod = ~AbsLat + altitude, 
                   random = list(
                     ~1 | Taxon_ID, 
                     ~1 | grid.ll),
                   method = "REML", 
                   data = MSBP.meta.max,
                   control=list(optimizer = "optim", optmethod = "Nelder-Mead"))
  
  summary(MaxLat)
#  Model Results:
#    
#    estimate      se     zval    pval    ci.lb    ci.ub     
#    intrcpt    34.2852  0.8960  38.2641  <.0001  32.5290  36.0413  ***
#    AbsLat     -0.1822  0.0275  -6.6156  <.0001  -0.2362  -0.1282  ***
#    altitude   -0.0024  0.0007  -3.6287  0.0003  -0.0037  -0.0011  ***
  

  
  
  #Plot of ranges against absolute latitude with lat effect, not significant R^2 0.73%
  Max.plot <- ggplot(data = MSBP.meta.max , aes(x = AbsLat)) +
    geom_rect(aes(xmin=0, xmax = 23.5, ymax =Inf, ymin = -Inf), fill = "grey90")+
    geom_point(aes(y = Tmax, alpha = log(1/vi)), colour = "red") +   
    theme_classic() +
    ylab(expression("Maximum Germiantion Temperature "(degree*C)))+
    xlab(expression("Absolute Latitude "(degree))) +
    geom_segment(x = 1, y = coef(MaxLat)[1]+ (coef(MaxLat)[2])*1, 
                 xend = 57, yend = (coef(MaxLat)[1]) + (coef(MaxLat)[2])*57, 
                 size = 1, colour = "red", linetype = "solid")+
    #coord_flip() +
    theme(legend.position="none", plot.margin = margin(10,10,10,30))+
    scale_y_continuous(limits = c(0,60), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,60), expand = c(0,0))
  
  Max.plot
  
  ggsave(file = "./Outputs/TmaxAbsLat.tiff", units="in", width=5, height=4, dpi=1000, compression = 'lzw')

  
  #Plot of ranges against absolute latitude with lat effect, not significant R^2 0.73%
  Max.Min.plot <- ggplot(data = MSBP.meta.max , aes(x = AbsLat)) +
    geom_rect(aes(xmin=0, xmax = 23.5, ymax =Inf, ymin = -Inf), fill = "grey90")+
    geom_point(aes(y = Tmax, alpha = log(1/vi)), colour = "red") +   
    theme_classic() +
    ylab(expression("Germiantion Temperature "(degree*C)))+
    xlab(expression("Absolute Latitude "(degree))) +
    geom_segment(x = 1, y = coef(MaxLat)[1]+ (coef(MaxLat)[2])*1, 
                 xend = 57, yend = (coef(MaxLat)[1]) + (coef(MaxLat)[2])*57, 
                 size = 1, colour = "red", linetype = "solid")+
    geom_point(data = MSBP.meta.min, aes(y = Tmin, alpha = log(1/vi)), colour = "blue") +   
    geom_segment(x = 1, y = coef(MinLat)[1]+ (coef(MinLat)[2])*1, 
                 xend = 58, yend = (coef(MinLat)[1]) + (coef(MinLat)[2])*58, 
                 size = 1, colour = "blue", linetype = "solid")+
    #coord_flip() +
    theme(legend.position="none", plot.margin = margin(10,10,10,30))+
    scale_y_continuous(limits = c(-15,60), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,60), expand = c(0,0))
  
  Max.Min.plot
  
  ggsave(file = "./Outputs/TmaxTminAbsLat.tiff", units="in", width=5, height=4, dpi=1000, compression = 'lzw')
  
  TmaxTmin <- list(Max.Model = MaxLat, 
                   Max.Model.Phy = MaxLat.phy, 
                   Max.Plot = Max.plot, 
                   Max.Model.h = MaxLat.h,
                   Max.Model.h.phy = MaxLat.h.phy,
                   Min.Model = MinLat, 
                   Min.Model.Phy = MinLat.phy, 
                   Min.Plot = Min.plot, 
                   Min.Model.h = MinLat.h,
                   Min.Model.h.phy = MinLat.h.phy,
                   Max.Min.plot = Max.Min.plot
                                      )
  
    
  
  return(TmaxTmin)
  
}