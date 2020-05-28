library(dplyr)
library(xtable)

MSBP <- read.csv(file = "./Data/MSBP_Clean")


data_appendix <- MSBP %>%
  rename(Species.Location = Grid.ID,
         Species = Taxon_ID,
         Location = grid.ll,
         Latitude.grid = Grid.Lat,
         Longitude.grid = Grid.Long,
         n.Records = n) %>%
  select(Species.Location,
         Tmin,
         Tmax,
         Topt,
         Tbreadth,
         Topt.upp,
         Topt.low,
         Toptbreadth,
         CM.Current, WR.Future,
         CM.Future, WR.Current,
         Altitude, SeedAge,
         n.Records)

print.xtable(xtable(data_appendix), "./Outputs/data_appendix.html", type = "html")
print.xtable(xtable(data_appendix), "./Outputs/data_appendix.tex", type = "latex")


pander::pander(data_appendix)
