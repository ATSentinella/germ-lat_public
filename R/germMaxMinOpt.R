#'Run tests against latitude 
#+ germ-max-min-opt, eval = FALSE


lat.tests <- function(MSBP, Temp, SETemp){
require(dplyr)
require(tidyr)
require(metafor)
require(ggplot2)
source("./R/getPhyloCor.R")
  
#Prepare data sets
MSBP.temp <- MSBP %>%
  mutate(vi = get(SETemp) * get(SETemp)) %>% #Standard error for opt
  filter(!is.na(get(Temp)) &
         !is.na(SeedAge) &
         !is.na(altitude))

#Remove duplicate species for phylogeny analysis               
MSBP.temp.phy <- distinct(MSBP.temp, Taxon_ID, .keep_all = T) 

#Get phylogeny tree[[1]] and correlation matrix[[2]]
MSBP.temp.phy <- getPhyloCor(MSBP.temp.phy)

model.h <- rma.mv(
  yi = get(Temp),
  V = vi,
  mod = ~ AbsLat + (AbsLat:NorthTF) + log10(SeedAge) + altitude,
  random = list(~ 1 | Taxon_ID, 
                ~ 1 | grid.ll),
  method = "REML",
  data = MSBP.temp
)


model.h.phy <- rma.mv(
  yi = get(Temp),
  V = vi,
  mod = ~ AbsLat + (AbsLat:NorthTF) + log10(SeedAge) + altitude,
  random = list( ~ 1 | phylo, 
                 ~ 1 | Taxon_ID, 
                 ~ 1 | grid.ll),
  R = list(phylo = MSBP.temp.phy$cor),
  method = "REML",
  data = (MSBP.temp.phy[[1]]),
  verbose = TRUE,
  digits = 5,
  control = list(optimizer = "optim", optmethod = "Nelder-Mead")
)


model <- rma.mv(
  yi = get(Temp),
  V = vi,
  mod = ~ AbsLat + log10(SeedAge) + altitude,
  random = list( ~ 1 | Taxon_ID, 
                 ~ 1 | grid.ll),
  method = "REML",
  data = MSBP.temp
)

model.phy <- rma.mv(
  yi = get(Temp),
  V = vi,
  mod = ~ AbsLat + log10(SeedAge) + altitude,
  random = list( ~ 1 | phylo, 
                 ~ 1 | Taxon_ID, 
                 ~ 1 | grid.ll),
  R = list(phylo = MSBP.temp.phy$cor),
  method = "REML",
  data = MSBP.temp.phy[[1]]
)


output <- list(
  Model = model,
  Model.Phy = model.phy,
  Data = MSBP.temp,
  Model.h = model.h,
  Model.h.Phy = model.h.phy
) 
  return(output)

}

MSBP <- read.csv(file = "./Data/MSBP_Clean") #core dataset

Tmax.lat <- lat.tests(MSBP = MSBP, Temp = "Tmax", SETemp = "SEmax")
Topt.upp.lat <- lat.tests(MSBP = MSBP, Temp = "Topt.upp", SETemp = "SEopt.upp")
Topt.low.lat <- lat.tests(MSBP = MSBP, Temp = "Topt.low", SETemp = "SEopt.low")
Tmin.lat <- lat.tests(MSBP = MSBP, Temp = "Tmin", SETemp = "SEmin")
Topt.range.lat <- lat.tests(MSBP = MSBP, Temp = "Toptbreadth", SETemp = "SEoptbreadth")



Tmax.lat$Model.h.Phy #hemi not significant
Tmax.lat$Model.h #hemi significant
Tmax.lat$Model.Phy #both abslat and alititide significant
Tmax.lat$Model #both abslat and alititide significant

Topt.upp.lat$Model.h.Phy #hemi not significant
Topt.upp.lat$Model.h #hemi not significant
Topt.upp.lat$Model.Phy #abslat significant
Topt.upp.lat$Model #abslat significant

Topt.low.lat$Model.h.Phy #hemi not significant
Topt.low.lat$Model.h #hemi not significant
Topt.low.lat$Model.Phy #abslat significant
Topt.low.lat$Model #abslat significant

Tmin.lat$Model.h.Phy #hemi not significant
Tmin.lat$Model.h #hemi not significant
Tmin.lat$Model.Phy #both abslat and alititide significant
Tmin.lat$Model #both abslat and alititide significant


Topt.range.lat$Model.h.Phy
Topt.range.lat$Model.h
Topt.range.lat$Model.Phy
Topt.range.lat$Model

R2(Tmax.lat$Model.Phy)

R2(Topt.upp.lat$Model.Phy)

R2(Topt.low.lat$Model.Phy)

R2(Tmin.lat$Model.Phy)


R2(Tmax.lat$Model)

R2(Topt.upp.lat$Model)

R2(Topt.low.lat$Model)

R2(Tmin.lat$Model)

maxminlines <- ggplot(data = filter(MSBP) , aes(x = AbsLat)) +
  theme_classic() +
  geom_point(data = Tmax.lat$data,
             aes(y = Tmax, alpha = log(1 / SEmax)),
             colour = "red") +
  geom_point(data = Tmin.lat$data,
             aes(y = Tmin, alpha = log(1 /SEmin)),
             colour = "blue") +
  geom_segment(x = 0.5, y = coef(Tmax.lat$Model)[1]+
                 (coef(Tmax.lat$Model)[2])*0.5+              
                 (coef(Tmax.lat$Model)[3])*3.07+
                 (coef(Tmax.lat$Model)[4])*462, 
               xend = 65, yend = (coef(Tmax.lat$Model)[1]) + 
                 (coef(Tmax.lat$Model)[2])*65+                
                 (coef(Tmax.lat$Model)[3])*3.07+
                 (coef(Tmax.lat$Model)[4])*462, 
               size = 1.5, colour = "red")+
 geom_segment(x = 0.5, y = coef(Tmin.lat$Model)[1]+
                (coef(Tmin.lat$Model)[2])*0.5+              
                (coef(Tmin.lat$Model)[3])*3.07+
                (coef(Tmin.lat$Model)[4])*462, 
               xend = 67, yend = (coef(Tmin.lat$Model)[1]) + 
                (coef(Tmin.lat$Model)[2])*67+ 
                (coef(Tmin.lat$Model)[3])*3.06+
                (coef(Tmin.lat$Model)[4])*462.3, 
               size = 1.5, colour = "blue")+
  ylab("Temperature (°C)")+
  xlab("Absolute Latitude (°)") +
  theme(legend.position="none")+
  scale_y_continuous(limits = c(-10,50)) +
  scale_x_continuous(limits = c(0,70), expand = c(0,0))

maxminlines



optlines <- ggplot(data = filter(MSBP) , aes(x = AbsLat)) +
  theme_classic() +
  geom_point(data = Topt.upp.lat$data,
             aes(y = Topt.upp, alpha = log(1 / SEopt.upp)),
             colour = "darkorange") +
  geom_point(data = Topt.low.lat$data,
             aes(y = Topt.low, alpha = log(1 /SEopt.low)),
             colour = "orange") +
  geom_segment(x = 0.5, y = coef(Topt.upp.lat$Model)[1]+
                 (coef(Topt.upp.lat$Model)[2])*0.5+              
                 (coef(Topt.upp.lat$Model)[3])*3.15569+
                 (coef(Topt.upp.lat$Model)[4])*427, 
               xend = 62, yend = (coef(Topt.upp.lat$Model)[1]) + 
                 (coef(Topt.upp.lat$Model)[2])*62+                
                 (coef(Topt.upp.lat$Model)[3])*3.15569+
                 (coef(Topt.upp.lat$Model)[4])*427, 
               size = 1.5, colour = "darkorange")+
  geom_segment(x = 0.5, y = coef(Topt.low.lat$Model)[1]+
                 (coef(Topt.low.lat$Model)[2])*0.5+              
                 (coef(Topt.low.lat$Model)[3])*3.15569+
                 (coef(Topt.low.lat$Model)[4])*427, 
               xend = 62, yend = (coef(Topt.low.lat$Model)[1]) + 
                 (coef(Topt.low.lat$Model)[2])*62+ 
                 (coef(Topt.low.lat$Model)[3])*3.15569+
                 (coef(Topt.low.lat$Model)[4])*427, 
               size = 1.5, colour = "orange")+
  geom_segment(aes(x = AbsLat, y = Topt.low, 
                   xend = AbsLat, yend = Topt.upp), 
               size = 0.5, colour = "orange")+
  ylab("Temperature (°C)")+
  xlab("Absolute Latitude (°)") +
  theme(legend.position="none")+
  scale_y_continuous(limits = c(0,50)) +
  scale_x_continuous(limits = c(0,70), expand = c(0,0))

optlines

ggsave(
  maxminlines,
  filename = "./Outputs/maxminlines.tiff",
  units = "in",
  width = 7,
  height = 5,
  dpi = 1000,
  compression = 'lzw')


ggsave(
  optlines,
  filename = "./Outputs/optlines.tiff",
  units = "in",
  width = 7,
  height = 5,
  dpi = 1000,
  compression = 'lzw')



