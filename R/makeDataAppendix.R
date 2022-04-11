library(dplyr)
library(xtable)

MSBP <- read.csv(file = "./Data/MSBP_Clean")


longevity <- read.csv(file = "./Data/Longevity/TRY.longevity.csv")

longevityfilt <- longevity %>% 
  mutate(ann_per = case_when(OrigValueStr == "annual" | 
                               OrigValueStr == "Annual" | 
                               OrigValueStr == "annuals" |
                               OrigValueStr == "summer annuals" |
                               OrigValueStr == "always annual" |
                               OrigValueStr == "winter annuals" |
                               OrigValueStr == "annual-winter annual" |
                               OrigValueStr == "winter annual" |
                               (OriglName == "Life history" & OrigValueStr == "1" ) |
                               (OriglName == "Plant phenology: Annual" & OrigValueStr == "yes" )     ~   "annual",  ######annuals
                             OrigValueStr == "perennial" | 
                               OrigValueStr == "Perennial" | 
                               OrigValueStr == "perennials" |                                   
                               OrigValueStr == "always pluriennial-pollakanthic" | 
                               (OriglName == "Plant phenology: Biennial" & OrigValueStr == "yes" ) | 
                               OrigValueStr == "perennial < 20 years" | 
                               OrigValueStr == "woody" | 
                               OrigValueStr == "perennial/woody" | 
                               OrigValueStr == "perennial > 20 years" | 
                               OrigValueStr == "poly-annuals > 50 years (long-lived perennials)" | 
                               OrigValueStr == "always biennial, always pluriennial-hapaxanthic" | 
                               OrigValueStr == "always biennial, always pluriennial-pollakanthic" | 
                               OrigValueStr == "tree" | 
                               OrigValueStr == "shrub" | 
                               OrigValueStr == "always pluriennial-hapaxanthic, always pluriennial-pollakanthic" | 
                               OrigValueStr == "always pluriennial-hapaxanthic" | 
                               OrigValueStr == "biennial" | 
                               OrigValueStr == "annual/biennial" | 
                               OrigValueStr == "poly-annuals < 5 years (short-lived perennials)" | 
                               OrigValueStr == "Biennial" | 
                               OrigValueStr == "biennial/perennial" | 
                               OrigValueStr == "always biennial" | 
                               OrigValueStr == "biennial-perennial" | 
                               OrigValueStr == "sometimes biennial, always pluriennial-hapaxanthic, sometimes pluriennial-pollakanthic" | 
                               OrigValueStr == "sometimes biennial, sometimes pluriennial-hapaxanthic, always pluriennial-pollakanthic" | 
                               OrigValueStr == "biennial/perennial/woody" | 
                               OrigValueStr == "sometimes biennial, always pluriennial-pollakanthic" | 
                               OrigValueStr == "poly-annuals 5-50 years (medium-lived perennials)" | 
                               (OriglName == "Plant phenology: Perennial" & OrigValueStr == "yes" )| 
                               (OriglName == "Plant phenology: Annual" & OrigValueStr == "no" )    ~   "perennial"  ###perennials
  ))

longevityfilt <- data.frame(lapply(longevityfilt, function(x) {gsub(" ", "_", x)}))

MSBPlong <- left_join(MSBP, longevityfilt, by = c("Taxon_ID" ="SpeciesName"))

woody <- read.csv(file = "./Data/GlobalWoodinessDatabase.csv")

woody <- data.frame(lapply(woody, function(x) {gsub(" ", "_", x)}))

MSBPlong <- left_join(MSBPlong, woody, by = c("Taxon_ID" ="gs"))

MSBPlong <- MSBPlong %>%
  mutate(ann_per=replace(ann_per, woodiness=="W", "perennial")) %>%
  distinct(Grid.ID, .keep_all = T)



data_appendix <- MSBPlong %>%
  rename(Species.Location = Grid.ID,
         Species = Taxon_ID,
         Location = grid.ll,
         Latitude.grid = Grid.Lat,
         Longitude.grid = Grid.Long,
         n.Records = n) %>%
  select(Species.Location,
         Tmin,
         Tmax,
         Tbreadth,
         Topt.upp,
         Topt.low,
         CM.Current, WR.Future,
         CM.Future, WR.Current,
         Altitude, SeedAge,
         woodiness, ann_per, 
         n.Records)

print.xtable(xtable(data_appendix), "./Outputs/data_appendix.tex", type = "latex")
