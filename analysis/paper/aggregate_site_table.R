# Create table of site characteristics
library(here)
library(tidyverse)


df <- read_csv(here("analysis/data/derived_data/FA_all_tidied.csv"))|>
  filter(Type == "Invert")|>
  select(-(4:60))|>
  group_by(Site, EcoRegion, Sample_Date, ForTot,Hay_per, UrbTot_per, Discharge_L.s, SRP_ug.L, NO3_ug.L,Canopy_Cover_per,Lat_Long,
           Temp, pH, Conductivity_uS.cm,DO_mg.L, Peri_Chla_mg.cm2)|>
  summarize(invert_taxa = paste(list(Taxon)))|>
  ungroup()|>
  mutate(invert_taxa = str_remove_all(invert_taxa, "\""),
         invert_taxa = str_remove(invert_taxa, "c"),
         invert_taxa = str_remove_all(invert_taxa, "\\(|\\)"))|>
  mutate(EcoRegion = str_replace(EcoRegion, "_", " "),
         Sample_Date2 = as_date(Sample_Date, format = "%d-%b-%y"))|>
  mutate(Latitude = str_split_fixed(Lat_Long, ",", 2)[, 1],
         Longitude = str_split_fixed(Lat_Long, ",", 2)[, 2],
         Longitude = str_trim(Longitude),
         .after = Site)|>
  select(-Lat_Long)|>arrange(desc(EcoRegion), Sample_Date2)|>
  mutate(Sample_Date = str_remove(Sample_Date,"-19"),
         Sample_Date = str_replace(Sample_Date, "-", " "))|>
  rename(`Invert Taxa` = "invert_taxa",
         Date = "Sample_Date",
         `Canopy Cover (%)` = "Canopy_Cover_per",
         `Discharge (L/s)` = "Discharge_L.s",
         `SRP (ug/L)` = "SRP_ug.L",
         `NO3-N (ug/L)` = "NO3_ug.L",
         `Periphyton Chl-a (mg/cm2)` = "Peri_Chla_mg.cm2",
         `Catchment Forest (%)` = "ForTot",
         `Catchment Hay/Pasture (%)` = "Hay_per",
         `Catchment Urbanized (%)` = "UrbTot_per",
         `DO (mg/L)` = "DO_mg.L",
         `Specific Conductance (uS/cm)` = "Conductivity_uS.cm",
         `Water Temperature (C)` = "Temp")|>
  mutate(across(c(Latitude, Longitude), as.numeric))|>
  mutate(across(c(Latitude, `Catchment Urbanized (%)`), str_trim))|>
  select(-`Invert Taxa`, -Sample_Date2, -`Catchment Hay/Pasture (%)`,
         -`Catchment Urbanized (%)`, -`DO (mg/L)`, -pH) #Cutting these for space in the table




df|>write_csv(here("analysis/paper/figures_and_tables/Table1_site_chars.csv"))
