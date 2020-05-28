
#' Figure 3. Germination Temperature Breadth Over Latitude
#' Viewing a subset of of the models
#' For Seed Germination over Latitude Project
#' 
#' @author Alex Sentinella ATSentinella@@gmail.com
#' Issues
#' - 
# Viewing a subset of of the groups
# 
# Issues:
#   
# - (resolved) when subsetting, old levels are retained (including levels with no data)
# - SOLUTION: use droplevels(df) to reset levels in a factor
# - good weights on one side implies better weighting on the other 

#"Harrisia_divaricata_17.5_-71.5") #good species
#"Eriosyce_kunzei_-29.5_-71.5") #Good opt, bad min/max
#"Phyllanthus_casticum_-22.5_43.5") #Good min, bad opt/max
#"Mandragora_autumnalis_32.5_35.5") #Weighted too much by 0's
#"Encelia_californica_30.5_-115.5") #Uninformative?
#"Symphyotrichum_concolor_35.5_-79.5") #Too wide
#"Ludwigia_alternifolia_35.5_-79.5") #lots of data, trend less clear
#"Scleria_depressa_11.5_-3.5" #Little data, but informative?
# "Rinorea_greveana_-21.5_43.5" # Too narrow
#"Acridocarpus_monodii_14.5_-3.5" )#Bad opt/max
#" Lithospermum_arvense_51.5_-1.5" # Why no model?
#"Prunus_cerasia_33.5_35.5")

require(tidyverse)


MSBP_MaxMinOpt <- read.csv(file = "./Outputs/MSBP_MaxMinOpt.csv")

Sample.Plot <- MSBP_MaxMinOpt %>% 
               group_by(Grid.ID) %>%
               droplevels()

Sample.Group <- sample(levels(Sample.Plot$Grid.ID), 1)# Number of groups

Sample.Plot <- Sample.Plot  %>% 
  filter(Grid.ID %in% 
           Sample.Group) #random species
#"Harrisia_divaricata_17.5_-71.5") #Good species to test
#"Ribes_nevadense_42.5_-122.5") #Looks to have a simple curve, not matched by estimates

quadbi <- function(x){
  acoef <- Sample.Plot$a[1]
  bcoef <- Sample.Plot$b[1]
  ccoef <- Sample.Plot$c[1]
  eta =(acoef*((x)^2) + bcoef*(x) + ccoef)
  exp(eta)/(1+exp(eta))
}

QuadPlot <- ggplot(data=Sample.Plot, aes(x = Test.Temp, y = propgerm, colour = Grid.ID)) + 
            geom_point(size=3) +
            geom_vline(aes(xintercept=Tmax, colour = Grid.ID), size=1, linetype = "longdash", colour = "red") +
            geom_vline(aes(xintercept=Tmin, colour = Grid.ID), size=1, linetype = "longdash", colour = "blue") +
            geom_vline(aes(xintercept=Topt, colour = Grid.ID), size=1, colour = "goldenrod2", linetype = "longdash") +
            geom_hline(aes(yintercept = 0)) + 
            facet_wrap(~Grid.ID) +guides(color = FALSE)+
            stat_function(fun=quadbi, colour="purple", size =1) +
            ylim(0,1) +
            theme_classic()+
            ylab("Proportion of seeds germinated")+
            xlab("Temperature (°C)")+
            labs(title = paste("SEmin:", signif(Sample.Plot$SEmin), 
                               ",  SEopt:", signif(Sample.Plot$SEopt), 
                               ",  SEmax:", signif(Sample.Plot$SEmax), 
                               ",  SEbreadth:", signif(Sample.Plot$SEbreadth)))


QuadPlot

t(Sample.Plot[1,])

#ggsave(filename = "../Outputs/assymetricalestimate.jpg", width = 7, height = 4)

#rm(list=ls())
