library(tidyverse)
library(sf)
library(ggmap)
library(ggsci)
library(here)


key <- read_delim(key_file, delim = "/t",col_names = FALSE)|>pull()
#Google's API to use, e.g., Stamen Maps, is free but must have a registered key
# I manually set 'key_file' to the location of my key on each computer I work on
# but do not want that key to be public on github.
register_google(key = key)


sites <- read_csv(here("analysis/data/raw_data/Ozarks_SiteData_Summarized.csv"))|>
  separate_wider_delim(cols = Lat_Long, names = c("lat", "lon"), delim = ",")|>
  mutate(across(c(lat,lon), as.numeric))|>
  mutate(EcoRegion = ifelse(EcoRegion == "Boston_Mountains", "Boston Mountains", "Ozark Highlands"))

er_shp <- st_read(here("analysis/data/raw_data/EcoRegion_shp/NA_CEC_Eco_Level3.shp"))
er_shp <- st_read(here("analysis/data/raw_data/EcoRegion_shp/NA_CEC_Eco_Level3.shp"))|>
  filter(NA_L3NAME %in% c("Boston Mountains", "Ozark Highlands"))|>
  st_transform(4326)

armo <- tigris::states(cb = TRUE, resolution = "20m", class = "sf") %>%
  filter(STUSPS %in% c("AR", "MO"))

state_labs <- tibble(name = c("Arkansas", "Missouri"), lat = c(36.4,36.6), lon = c(-94.2,-94.2))

#terr_map <- get_map(location = c(lat = 36, lon = -92), maptype = "satellite", zoom =7)

ggmap(terr_map)+
  geom_sf(data = er_shp, aes(color = NA_L3NAME, fill = NA_L3NAME),
          inherit.aes = FALSE, alpha = 0.1, linewidth = 2)+
  geom_sf(data = armo, inherit.aes = FALSE, alpha = 0,
          color = "black", linewidth = 1.5)+
  geom_point(data = sites, aes(y = lat, x = lon, fill = EcoRegion),
             inherit.aes = FALSE,  size = 5, shape = 21, color = "black")+
  scale_fill_jco()+
  scale_color_jco()+
  scale_y_continuous(limits = c(35,37))+
  scale_x_continuous(limits = c(-94.5,-90))+
  geom_text(data = state_labs, aes(x = lon, y = lat, label = name), size = 8)+
    theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.78,0.15),
        legend.background = element_blank(),
        legend.text = element_text(size = 18, face = "bold"))
