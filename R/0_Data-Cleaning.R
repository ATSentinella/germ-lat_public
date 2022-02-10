#+ eval = FALSE

# This script details the data consolidation and cleaning process for the MSBP
# data sets

library(tidyr)
library(ggplot2)
library(dplyr)
library(readxl)

# We start wih the data sets downloaded from the MSBP website -

# Download 1: AccessionNumber, Family, Taxon, Latitude, Longitude, Cultivated (T/F)
# Download 2: AccessionNumber,  Genus, Sp1 (species), Country, Altitude
# Download 3: Accession Number, BrahmsGermTestId, StepSummary
# Download 4: BrahmsGermTestId, Rate, InspectionInterval, Interval
# Download 5: BrahmsGermTestId, Sown, Germinated, Family
# Download 6: BrahmsGermTestId, Fresh, Empty, Infested, Abnormal
# Download 7: BrahmsGermTestId, Unusable, Mouldy, TestNotes

# 1 and 2 are accession info
#
# 3-7 are germination information
#
# See MSBP Data Standards for more information (other varibles possible, many
# blanks at this stage though)

MSBP_DL1 <-
  read_excel("./Data/171117MSBPdownload/2017-11-16_231507939-Download_1.xlsx")
MSBP_DL2 <-
  read_excel("./Data/171117MSBPdownload/2017-11-16_231507939-Download_2.xlsx")
MSBP_DL3 <-
  read_excel("./Data/171117MSBPdownload/2017-11-16_233020054-Download_3.xlsx")
MSBP_DL4 <-
  read_excel("./Data/171117MSBPdownload/2017-11-16_233902116-Download_4.xlsx")
MSBP_DL5 <-
  read_excel("./Data/171117MSBPdownload/2017-11-16_234637337-Download_5.xlsx")
MSBP_DL6 <-
  read_excel("./Data/171117MSBPdownload/2017-11-16_234948032-Download_6.xlsx")
MSBP_DL7 <-
  read_excel("./Data/171117MSBPdownload/2017-11-16_235202784-Download_7.xlsx")

# Merge Data sets and create raw data files
# - This is done by row, as this was downloaded from the same large dataset.
#   Checks are done to ensure rows match up.

Accession <- bind_cols(MSBP_DL1, MSBP_DL2)

Germ <- bind_cols(MSBP_DL3, MSBP_DL4, MSBP_DL5, MSBP_DL6, MSBP_DL7)

#Save Raw, combined data - Accession
Accession_Raw <- select(Accession,-AccessionNumber...7) %>%
  rename(AccessionNumber = AccessionNumber...1)
write.csv(Accession_Raw, file = "./Data/Accession_Raw.csv")

#Save Raw, combined data - Germination
Germ_Raw <-
  select(
    Germ,
    -BrahmsGermTestId...4,
    -BrahmsGermTestId...9,
    -BrahmsGermTestId...12,
    -BrahmsGermTestId...17
  ) %>%
  rename(BrahmsGermTestId = BrahmsGermTestId...2)
write.csv(Germ_Raw, file = "./Data/Germ_Raw.csv")


#  Start to clean germination data set
#
#
#  Start with the following variables:
#    AccessionNumber, BrahmsGermTestID, StepSummary, NumSown, NumGerm
#
#  Steps:
#    - Start with:
#    AccessionNumber, BrahmsGermTestID, StepSummary, NumSown, NumGerm
#  - Remove Duplicates
#  - Make BrahmsGermTestID, NumSown, NumGerm numbers
#  - Remove Negative Germination IDs
#  - Remove rows with Number Sown = 0
#  - Rationalised duplicate GermIDs

#Starting data set 
Germ_Cleaning <- select(Germ_Raw, one_of(
  c(
    "AccessionNumber",
    "BrahmsGermTestId",
    "StepSummary",
    "NumSown",
    "NumGerm"
  )
))


#Remove duplicates 
Germ_Cleaning <- distinct(Germ_Cleaning)

#Makes columns number, rows with characters as NA
Germ_Cleaning <- transform(
  Germ_Cleaning,
  BrahmsGermTestId = as.numeric(BrahmsGermTestId),
  NumSown = as.numeric(NumSown),
  NumGerm = as.numeric(NumGerm),
  AccessionNumber = as.numeric(AccessionNumber)
)


#Remove negative IDS and NAs 
Germ_Cleaning <- filter(Germ_Cleaning, AccessionNumber > 0) %>%
  filter(BrahmsGermTestId > 0)

#Remove NumSown = 0 
Germ_Cleaning <- filter(Germ_Cleaning, NumSown > 0)

#Number of each ID (should be 1)
countgerm <- count(Germ_Cleaning, BrahmsGermTestId)

#Added a row showing where there are double IDs
Germ_Cleaning <-
  left_join(Germ_Cleaning, countgerm, by = "BrahmsGermTestId")

#Kept 1 of duplucate IDs where NumSown and NumGerm are the same
Germ_Cleaning <-
  distinct(Germ_Cleaning,
           AccessionNumber,
           BrahmsGermTestId,
           NumSown,
           NumGerm,
           .keep_all = TRUE)

countgerm <- count(Germ_Cleaning, BrahmsGermTestId)

Germ_Cleaning <-
  left_join(Germ_Cleaning, countgerm, by = "BrahmsGermTestId")

#Removed GermIds with different NumGerms, can't be sure which is correct 
Germ_Cleaning <- filter(Germ_Cleaning, n.y == 1)

Germ_Cleaning <- select(Germ_Cleaning,-n.x,-n.y)

# Add in temperature information
#
# - Find "step""
# - Find First temp value
# - Return all temps
#
# Issues:
# - From excel method at this stage, need to work out R method (resolved)
# - This leaves non filtered values with NAs (resolved)
# - These are filtered out at this stage (resolved)
# - steps aren't all included e.g. missing step 1 or 4 (sorts itself out in code)

Germ_Temps <- Germ_Cleaning

#No. of Steps
Germ_Temps <-
  mutate(Germ_Temps, Steps = lengths(regmatches(
    StepSummary, gregexpr("*Step", StepSummary)
  )))

#Remove Rows with 10 or more steps
Germ_Temps <- filter(Germ_Temps, Steps < 10)

#Temp 1, if it doesn't find "1:T" it will return everything, 
# so that should be NA
Germ_Temps <- Germ_Temps %>%
  mutate(
    Temp1 = ifelse(
      startsWith(sub(".*?1: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?1: T(.*?)(/.*|$)", "\\1", StepSummary)
    ),
    Temp2 = ifelse(
      startsWith(sub(".*?2: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?2: T(.*?)(/.*|$)", "\\1", StepSummary)
    ),
    Temp3 = ifelse(
      startsWith(sub(".*?3: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?3: T(.*?)(/.*|$)", "\\1", StepSummary)
    ),
    Temp4 = ifelse(
      startsWith(sub(".*?4: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?4: T(.*?)(/.*|$)", "\\1", StepSummary)
    ),
    Temp5 = ifelse(
      startsWith(sub(".*?5: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?5: T(.*?)(/.*|$)", "\\1", StepSummary)
    ),
    Temp6 = ifelse(
      startsWith(sub(".*?6: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?6: T(.*?)(/.*|$)", "\\1", StepSummary)
    ),
    Temp7 = ifelse(
      startsWith(sub(".*?7: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?7: T(.*?)(/.*|$)", "\\1", StepSummary)
    ),
    Temp8 = ifelse(
      startsWith(sub(".*?8: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?8: T(.*?)(/.*|$)", "\\1", StepSummary)
    ),
    Temp9 = ifelse(
      startsWith(sub(".*?9: T(.*?)(/.*|$)", "\\1", StepSummary), "*"),
      NA,
      sub(".*?9: T(.*?)(/.*|$)", "\\1", StepSummary)
    )
  )

#Count Unique temps
GermTempCount <- select(Germ_Temps, contains("Temp"))

GermTempCount$nunique <-
  apply(GermTempCount, 1, function(x)
    length(unique(na.omit(x))))

GermTempCount <- select(GermTempCount, nunique)

Germ_Temps <- bind_cols(Germ_Temps, GermTempCount)

# Now we have the first 9 temperatures and the count of unique temperatures tested
#
# Steps:
#
# - If only 1 temperature tested, that is Test.temp
# - If only 2 temperatures tested
# - THEN If the first temeperture was 60 and above
# - First different temperature after that is Test.temp
# - ELSE THEN If the first temeperture was 0 or less
# - First different temperature after that is Test.temp
# - ELSE "Too many"
# - The rest "Too many"
#
# - Any first temperature 60 and above was considered a heatshock,
# while some below could be it was too hard to determine if they were different to a summer treatment
# - Find the first non heatshock temperature, subtract 1 from nunique
# - If the only temeperature tested is 0, that is correct
# -If 0s followed by a single temperature, we assume coldshock or null value
# - We should now have a filtered list with a single temperature tested
#
# Issues:
# - heatshock defined as >=60 (in methods)
# - 0 should be ignored or considered coldshock (in methods)
# - How to handle NA  values? (removed)
# - How to check for alternating 0 values (e.g. 0, 5, 0, 5) (Ignore)


Germ_Temps <- transform(Germ_Temps, Temp1 = as.numeric(Temp1))
Germ_Temps <- transform(Germ_Temps, Temp2 = as.numeric(Temp2))
Germ_Temps <- transform(Germ_Temps, Temp3 = as.numeric(Temp3))
Germ_Temps <- transform(Germ_Temps, Temp4 = as.numeric(Temp4))
Germ_Temps <- transform(Germ_Temps, Temp5 = as.numeric(Temp5))
Germ_Temps <- transform(Germ_Temps, Temp6 = as.numeric(Temp6))
Germ_Temps <- transform(Germ_Temps, Temp7 = as.numeric(Temp7))
Germ_Temps <- transform(Germ_Temps, Temp8 = as.numeric(Temp8))
Germ_Temps <- transform(Germ_Temps, Temp9 = as.numeric(Temp9))

Germ_Temps <-
  mutate(Germ_Temps, Heatshock = ifelse(is.na(Temp1), (ifelse(
    Temp2 >= 60, TRUE, FALSE
  )),
  (ifelse(
    Temp1 >= 60, TRUE, FALSE
  ))))


Germ_Temps <- Germ_Temps %>%
  mutate(Test.Temp = ifelse(
    nunique == 1,
    ifelse(is.na(Temp1), Temp2, Temp1),
    # If 1 Temp
    ifelse(
      Heatshock == TRUE,
      #Heatshock
      ifelse(nunique == 2,
             ifelse(
               Temp2 == Temp1,
               ifelse(Temp3 == Temp2,
                      ifelse(
                        Temp4 == Temp3,
                        ifelse(Temp5 == Temp4,
                               ifelse(
                                 Temp6 == Temp5,
                                 ifelse(Temp7 == Temp6,
                                        ifelse(Temp8 == Temp7,
                                               Temp9,
                                               Temp8),
                                        Temp7),
                                 Temp6
                               ),
                               Temp5),
                        Temp4
                      ),
                      Temp3),
               Temp2
             ),
             "Too many"),
      ifelse(Temp1 <= 0, #Coldshock or error
             ifelse(
               nunique == 2,
               ifelse(Temp2 == Temp1,
                      ifelse(
                        Temp3 == Temp2,
                        ifelse(Temp4 == Temp3,
                               ifelse(
                                 Temp5 == Temp4,
                                 ifelse(Temp6 == Temp5,
                                        ifelse(
                                          Temp7 == Temp6,
                                          ifelse(Temp8 == Temp7,
                                                 Temp9,
                                                 Temp8),
                                          Temp7
                                        ),
                                        Temp6),
                                 Temp5
                               ),
                               Temp4),
                        Temp3
                      ),
                      Temp2),
               "Too many"
             ),
             "Too many")
    )
  ))

#Remove any rows with "Too many" or NA
Single_Temps <-
  filter(Germ_Temps,!Test.Temp %in% c('Too many', NA))

Single_Temps <-
  transform(Single_Temps, Test.Temp = as.numeric(Test.Temp))

Single_Temps <- filter(Single_Temps, Test.Temp < 100)

Single_Temps <-
  select(Single_Temps,-(
    c(
      Temp1,
      Temp2,
      Temp3,
      Temp4,
      Temp5,
      Temp6,
      Temp7,
      Temp8,
      Temp9,
      Steps,
      nunique,
      StepSummary
    )
  ))

# Adding in taxon information
#
# Important variables: AccessionNumber, Family, Taxon, Latitude, Longitude,
# Cultivated (not relevant)
#
# Steps:
#
# - Merge taxon information based on Accession number - For blank taxon
# information, *Family_[Accession number]  is used - Filtered out rows with
# cultivated =True (not anymore)

x <-
  distinct(select(Accession_Raw, one_of(
    c(
      "AccessionNumber",
      "Family",
      "Genus",
      "Sp1",
      "Latitude",
      "Longitude",
      "CultivatedFlag"
    )
  )))


# arrange to have lat first
x <- arrange(x, Latitude)

x <-
  distinct(x, AccessionNumber, Family, Genus, Sp1, .keep_all = TRUE)

countx <- count(x, AccessionNumber)

uniquex <- left_join(x, countx, by  = "AccessionNumber")

uniquex <- uniquex %>%
  filter (n == 1) %>%
  transform(AccessionNumber = as.numeric(AccessionNumber))

#This has filtered out all Accession numbers with unclear taxon information
Germ_Taxon <-
  inner_join(Single_Temps, uniquex, by  = "AccessionNumber")

Germ_Taxon <- select(Germ_Taxon,-(c(CultivatedFlag, n, Heatshock)))

Germ_Taxon <- unite(Germ_Taxon, Taxon_ID, c("Genus", "Sp1"))

Germ_Taxon <- Germ_Taxon %>%
  mutate(Taxon_ID = ifelse(
    Taxon_ID == "Indet_sp.",
    paste("*", Family, AccessionNumber, sep = ""),
    Taxon_ID
  ))

Germ_Taxon <- select(Germ_Taxon,-(c(Family)))

#Moves columns to the front, use -ColumnID to move to back
Germ_Taxon <- select(Germ_Taxon, Taxon_ID, everything())

# Filtering based on Latitude
#
# Steps:
# - Convert lat/long to numeric
# - Round latitudes using floor(), and add 0.5 to centre
# - Make new ID which is Grid.ID - Species_Genus_Lat_Long

Germ_Lat <- filter(Germ_Taxon,!is.na(Latitude))

Germ_Lat <- filter(Germ_Taxon,!is.na(Longitude))

Germ_Lat <- transform(Germ_Lat, Latitude = as.numeric(Latitude))
Germ_Lat <- transform(Germ_Lat, Longitude = as.numeric(Longitude))

Germ_Lat <- mutate(Germ_Lat, Grid.Lat = floor(Latitude) + 0.5)
Germ_Lat <- mutate(Germ_Lat, Grid.Long = floor(Longitude) + 0.5)

Germ_Lat <-
  mutate(Germ_Lat, Grid.ID = paste(Taxon_ID, Grid.Lat, Grid.Long, sep = "_"))

Germ_Lat <- select(Germ_Lat, Grid.ID, everything())

Germ_Lat <- transform(Germ_Lat, Test.Temp = as.numeric(Test.Temp))

Germ_Filt <- Germ_Lat %>%
  group_by(Grid.ID) %>%
  filter(n_distinct(Test.Temp) > 3)

Germ_Filt <- Germ_Filt %>%
  group_by(Grid.ID) %>%
  filter(max(NumGerm) != 0)

Germ_Test <-
  transform(Germ_Raw, BrahmsGermTestId = as.numeric(BrahmsGermTestId))
Germ_Filt <-
  left_join(Germ_Filt, Germ_Test, by = "BrahmsGermTestId")

Germ_Filt <- Germ_Filt %>%
  transform(Full = as.numeric(Full)) %>%
  mutate(
    NumNotGermFull = Full,
    NumSownFull = NumGerm.x + Full,
    NumGerm = NumGerm.x,
    NumSownOrig = NumSown.x,
    AccessionNumber = AccessionNumber.x
  )

MSBP_Collect <-
  read_excel("./Data/2019-09-26_035309887-BRAHMSOnlineData.xlsx") %>%
  transform(
    AccessionNumber = as.integer(AccessionNumber),
    DateCollected = as.Date(DateCollected)
  ) %>%
  select(-Taxon)

MSBP_Test  <-
  read_excel("./Data/2019-09-26_044526851-BRAHMSOnlineData.xlsx") %>%
  separate(DateStarted,
           sep = " ",
           into = c("Date", "Time", "AMPM")) %>%
  select(c(BrahmsGermTestId, Date)) %>%
  filter(!is.na(Date)) %>%
  transform(
    BrahmsGermTestId = as.integer(BrahmsGermTestId),
    GermDate = as.Date(Date, tryFormats = "%m/%d/%Y")
  ) %>%
  select(c(BrahmsGermTestId, GermDate))

Germ_Filt <- Germ_Filt %>%
  left_join(MSBP_Collect, by = "AccessionNumber") %>%
  left_join(MSBP_Test, by = "BrahmsGermTestId") %>%
  mutate(SeedAge = GermDate - DateCollected)

Germ_Filt <-   select(
  Germ_Filt,
  c(
    "Grid.ID"           ,
    "Taxon_ID"   ,
    "AccessionNumber" ,
    "BrahmsGermTestId"  ,
    "Test.Temp"         ,
    "Latitude"   ,
    "Longitude"         ,
    "Grid.Lat"          ,
    "Grid.Long"  ,
    "NumNotGermFull"    ,
    "NumSownFull"       ,
    "NumGerm"    ,
    "NumSownOrig",
    "SeedAge",
    "GermDate",
    "DateCollected"
  )
)

write.csv(Germ_Filt, file = "./Outputs/Germ_Filtered.csv")

rm(list = ls())
