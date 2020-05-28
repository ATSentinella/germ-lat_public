#' This is code to create tables for all model tests
#' For Seed Germination over Latitude Project

getGermTable <- function(model, title){
  x <- xtable(caption = title, digits = 4, 
              coef(summary(model))
  ) 
  return(x)
}

getTables <- function(models, model.name) {
  require(tidyverse)
  require(metafor)
  require(xtable)


model.tables <- list()



if (length(models)>3) {
  model.tables[[1]] <- getGermTable(models[[5]], 
                                    paste(model.name, 
                                          "over Absolute Latitude and Hemisphere (with Phylogeny)",
                                          "n =", models[[5]]$k))
  model.tables[[2]] <- getGermTable(models[[4]], 
                                    paste(model.name, 
                                          "over Absolute Latitude and Hemisphere (without Phylogeny)",
                                          "n =", models[[4]]$k))
  model.tables[[3]] <- getGermTable(models[[2]], 
                                    paste(model.name, 
                                          "over Absolute Latitude (with Phylogeny)",
                                          "n =", models[[2]]$k))
  model.tables[[4]] <- getGermTable(models[[1]], 
                                    paste(model.name, 
                                          "over Absolute Latitude (without Phylogeny)",
                                          "n =", models[[1]]$k))
} else {
  model.tables[[1]] <- getGermTable(models[[2]], 
                                    paste(model.name, 
                                        "over Absolute Latitude and Hemisphere (with Phylogeny)",
                                        "n =", models[[2]]$k))
  model.tables[[2]] <- getGermTable(models[[1]], 
                                    paste(model.name, "over Absolute Latitude and Hemisphere (with Phylogeny)",
                                        "n =", models[[1]]$k))
}

return(model.tables)
         
}

S2_Tables <- list()

S2_Tables[[1]] <- getTables(TB, "Tolerance Breadth")
S2_Tables[[2]]  <- getTables(CM, "Current Climate Mismatch")
S2_Tables[[3]]  <- getTables(WR, "Future Warming Risk")
S2_Tables[[4]]  <- getTables(WRC, "Current Warming Risk")
S2_Tables[[5]]  <- getTables(CMF, "Future Climate Mismatch")
S2_Tables[[6]] <- getTables(Tmax.lat, "Maximum Germination Temperature")
S2_Tables[[7]] <-  getTables(Topt.upp.lat, "Upper Optimum Germination Temperature")
S2_Tables[[8]] <-  getTables(Topt.low.lat, "Lower Optimum Germination Temperature")
S2_Tables[[9]] <-  getTables(Tmin.lat, "Minimum Germination Temperature")
S2_Tables[[10]] <-  getTables(Topt.range.lat, "Optimum Germination Temperature Breadth")
S2_Tables[[11]] <- getGermTable(WR.wood.life$Wood.Model.Phy, 
                                paste("Future Warming Risk by growth form (with Phylogeny)",
                                    "n =", WR.wood.life$Wood.Model.Phy$k))
S2_Tables[[12]] <- getGermTable(WR.wood.life$Wood.Model, 
                                paste("Future Warming Risk by growth form (without Phylogeny)",
                                    "n =", WR.wood.life$Wood.Model$k))
S2_Tables[[13]] <- getGermTable(WR.wood.life$Long.Model.Phy, 
                                paste("Future Warming Risk by longevity (with Phylogeny)",
                                    "n =", WR.wood.life$Long.Model.Phy$k))
S2_Tables[[14]] <- getGermTable(WR.wood.life$Long.Model, 
                                paste("Future Warming Risk by longevity (without Phylogeny)",
                                    "n =", WR.wood.life$Long.Model$k))
saveRDS(S2_Tables, "S2_Tables")

write(S2_Tables, "./Outputs/S2_Tables.txt")
