#' Figure 2. Example germination percentage plot
#' Uses Harrisia divaricata (17.5, -71.5)
#' For Seed Germination over Latitude Project
#'
#' @author Alex Sentinella ATSentinella@@gmail.com
#' @param Germ_Filtered.csv filereted germiantion data
#' @return MSBP_MaxMinOpt.csv


library(ggplot2)
library(dplyr)
library(tidyr)

MSBP_MaxMinOpt <- read.csv(file = "./Outputs/MSBP_MaxMinOpt.csv")

Sample.Plot <- MSBP_MaxMinOpt %>%
  group_by(Grid.ID)

Sample.Plot <- droplevels(Sample.Plot)

Sample.Group <- sample(levels(Sample.Plot$Grid.ID), 1) # Number of groups, in this case 1

Sample.Plot <- Sample.Plot%>%
  filter(Grid.ID %in%
           "Harrisia_divaricata_17.5_-71.5") #You can change this other species, or even multiple species

quadbi <- function(x) {
  acoef <- Sample.Plot$a[1]
  bcoef <- Sample.Plot$b[1]
  ccoef <- Sample.Plot$c[1]
  eta = (acoef * ((x) ^ 2) + bcoef * (x) + ccoef)
  exp(eta) / (1 + exp(eta))
}

germfigure <-
  ggplot(data = Sample.Plot, aes(x = Test.Temp, y = propgerm, colour = Grid.ID)) +
  geom_point(size = 3) +
  geom_vline(
    aes(xintercept = Tmax, colour = Grid.ID),
    size = 1,
    linetype = "longdash",
    colour = "red"
  ) +
  geom_vline(
    aes(xintercept = Tmin, colour = Grid.ID),
    size = 1,
    linetype = "longdash",
    colour = "blue"
  ) +
  geom_vline(
    aes(xintercept = Topt.low, colour = Grid.ID),
    size = 1,
    colour = "goldenrod2",
    linetype = "longdash"
  ) +  
  geom_vline(
    aes(xintercept = Topt.upp, colour = Grid.ID),
    size = 1,
    colour = "goldenrod2",
    linetype = "longdash"
  ) +
  geom_hline(aes(yintercept = 0)) +
  guides(color = FALSE) +
  stat_function(fun = quadbi, colour = "purple", size =
                  1) +
  theme_classic() +
  ylab("Germination Percentage (%)") +
  xlab(expression("Temperature " (degree * C))) +
  theme(plot.title = element_blank())

#Show figure
germfigure

#Save figure
ggsave(germfigure,
       filename = "./Outputs/germinationmodelfigure.png",
       width = 7,
       height = 4)
