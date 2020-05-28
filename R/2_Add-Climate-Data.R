#+ climate-data, eval = FALSE

getClimateData <- function() {
  # Load Packages
  require(dplyr)
  require(tidyr)
  require(raster)
  require(sp)
  require("rgdal")
  
  MSBP <- tbl_df(read.csv(file = "./Outputs/MSBP_MaxMinOpt.csv"))
  
  # Load Current Climate Data
  Current <- getData("worldclim",
                     var = "bio",
                     res = 10,
                     path = "./")
  Current <- Current[[c(1, 8)]]
  names(Current) <- c("Current.Temp", "Current.Temp.Hot.Quart")
  
  # weighted means of each polygon, removes NAs first, takes a long time
  coords <-
    dplyr::select(MSBP, one_of(c("Longitude", "Latitude"))) %>%
    rename(x = Longitude, y = Latitude)
  points <- SpatialPoints(coords, proj4string = Current@crs)
  values <- raster::extract(Current, points)
  Currentdf <- cbind.data.frame(coordinates(points), values)
  
  ## Future (2070) Average Temperature, rcp =85, model = AC
  Future85 <-  getData('CMIP5',
      var = 'bio',
      res = 10,
      rcp = 85,
      model = 'AC',
      year = 70,
      path = "./")
  Future85 <- Future85[[c(1, 8)]]
  names(Future85) <- c("Future.Temp", "Future.Temp.Hot.Quart")
  
  coords <- dplyr::select(MSBP, one_of(c("Longitude", "Latitude"))) %>%
    rename(x = Longitude, y = Latitude)
  points <- SpatialPoints(coords, proj4string = Future85@crs)
  values <- raster::extract(Future85, points)
  Future85df <- cbind.data.frame(coordinates(points), values)
  
  ## Altitude from worldclim dataset
  Altitude.dl <- getData('worldclim', var = 'alt', res = 10)
  Altitude.dl <- Altitude.dl[[1]]
  names(Altitude.dl) <- "Altitude"
  
  coords <- dplyr::select(MSBP, one_of(c("Longitude", "Latitude"))) %>%
    rename(x = Longitude, y = Latitude)
  points <- SpatialPoints(coords, proj4string = Altitude.dl@crs)
  values <- raster::extract(Altitude.dl, points)
  Altitudedf <- cbind.data.frame(coordinates(points), values)
  
  #Making one table
  
  MSBP.Climate <- cbind.data.frame(MSBP, 
                                   Currentdf[3:4], 
                                   Future85df[3:4], 
                                   Altitudedf[3])
  
  detach("package:raster", unload = TRUE)
  
  mean.na <- function (...)  {mean(..., na.rm = T)}
  
  MSBP.Climate <- MSBP.Climate %>%
    mutate(Altitude = ifelse(values < 0, 0, values)) %>%
    select(-values) %>%
    group_by(Grid.ID) %>%
    mutate(
      altitude = mean.na(Altitude),
      SeedAge = mean.na(SeedAge),
      Current.Temp = mean.na(Current.Temp) / 10,
      Future.Temp = mean.na(Future.Temp) / 10,
      Current.Hot.Quart = mean.na(Current.Temp.Hot.Quart) /  10,
      Future.Hot.Quart = mean.na(Future.Temp.Hot.Quart) / 10
    ) %>%
    mutate(
      CM.Current = Current.Hot.Quart - Topt.upp,
      WR.Future = Future.Hot.Quart - Tmax,
      CM.Future = Future.Hot.Quart - Topt.upp,
      WR.Current = Current.Hot.Quart - Tmax,
      TB = Tmax - Tmin,
      AbsLat = abs(Grid.Lat),
      Hemisphere = ifelse(Grid.Lat >= 0, "N", "S"),
      Warming = Future.Temp - Current.Temp
    ) %>%
    as_tibble() %>%
    distinct(Grid.ID,  .keep_all = T) %>%
  select(-c("X")) %>%
  unite(grid.ll, c(Grid.Lat, Grid.Long), remove =F)%>%
  mutate(NorthTF = ifelse(Hemisphere == "N", T, F)) %>%
  filter(!is.na(Tmax)|!is.na(Tmin))

write.csv(MSBP.Climate, file = "./Data/MSBP_Clean", row.names = FALSE)
  
}

getClimateData()
