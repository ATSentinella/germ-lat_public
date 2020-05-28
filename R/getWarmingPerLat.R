# Load Packages

library(raster)
library(sp)
library("rgdal")
library(ggplot2)
library(dplyr)
library(tidyr)


# Load Climate Data 
# BIO1 = Annual Mean Temperature
Current <- getData("worldclim",var="bio",res=10)
Current <- Current[[1]]
names(Current) <- "Temp"

## Future (2070) Average Temperature, rcp =85, model = AC
Future85 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=70)
Future85 <- Future85[[1]]
names(Future85) <- "Temp"


poly60n <- rbind(c(-180, 60), c(180, 60), c(180, 61), c(-180, 61)) 
poly5n <- rbind(c(-180, 5), c(180,5), c(180, 6), c(-180, 6)) 
poly33s <- rbind(c(-180, -33), c(180, -33), c(180, -34), c(-180, -34))

polys <- spPolygons(poly60n, poly5n, poly33s)

CurrentMeans <- raster::extract(Current[[1]], polys, weights=FALSE, fun=mean, na.rm=TRUE)
FutureMeans <- raster::extract(Future85[[1]], polys, weights=FALSE, fun=mean, na.rm=TRUE)

Diff.Table <- cbind(CurrentMeans, FutureMeans) %>%
  as.data.frame() %>%
  mutate(ChangeTemp = (V2-V1)/10)

write.csv(Diff.Table, "./Outputs/difftable.csv")

detach("package:raster", unload=TRUE)



