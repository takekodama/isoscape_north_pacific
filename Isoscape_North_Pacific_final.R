# ==============================================================================
# Pipeline: POM & Fish d15N Isoscape Modeling & Spatial Variance Decomposition
# North Pacific Oceanographic & Isotopic Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Libraries & Global Environment Setup
# ------------------------------------------------------------------------------
# Core Data Manipulation & I/O
library(tidyverse)
library(readxl)
library(lubridate)
library(magrittr)
library(data.table)

# Spatial & Oceanographic Tools
library(sf)
library(terra)
library(akima)
library(mapdata)
library(tidync)
library(metR)
library(gstat)
library(stars)
library(rnaturalearth)

# Modeling & Bayesian Inference
library(mgcv)
library(mgcViz)
library(MuMIn)
library(brms)
library(cmdstanr)
library(posterior)

# Visualization & Layout
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggrastr)
library(ggsci)
library(MetBrewer)

# Global Options & Color Palettes
select <- dplyr::select
options(na.action = "na.fail")
colors_rev <- met.brewer("Nizami", direction = -1)

# Base Map Layer Definition (World2 projection for North Pacific)
world_map <- map_data("world2")
geom_base_map <- geom_map(
  data = world_map, map = world_map,
  aes(x = long, y = lat, map_id = region),
  color = NA, fill = "gray30"
)


# ------------------------------------------------------------------------------
# 2. Custom Utility Functions
# ------------------------------------------------------------------------------

#' Remove outliers based on Grubbs-like t-test criterion
remove_outliers <- function(x, conf.level = 0.95) {
  x <- x[!is.na(x)]
  del.val <- NULL
  while (TRUE) {
    n <- length(x)
    if (n < 3) break
    r <- range(x)
    t <- abs(r - mean(x)) / sd(x)
    q2 <- (n - 2) / ((n - 1)^2 / t^2 / n - 1)
    q2[q2 < 0] <- 0
    q <- sqrt(q2)
    p <- n * pt(q, n - 2, lower.tail = FALSE)
    if (t[1] < t[2]) {
      if (p[2] < 1 - conf.level) {
        del.val <- c(del.val, r[2])
        x <- x[x != r[2]]
        next
      }
    } else {
      if (p[1] < 1 - conf.level) {
        del.val <- c(del.val, r[1])
        x <- x[x != r[1]]
        next
      }
    }
    break
  }
  return(list(x = x, del.val = del.val))
}

#' Assign oceanic region classifications based on geographic coordinates
assign_oceanic_areas <- function(df) {
  df %>%
    mutate(
      area1 = case_when(
        long <= 139 & lat <= 34.75 ~ "Kuroshio",
        long <= 141 & lat >= 34.75 ~ "SoJ",
        long <= 142 & lat >= 44 ~ "SoJ",
        long > 140 & long <= 144 & lat <= 43 & lat > 40 ~ "SoJ",
        long > 142 & long < 146 & lat >= 44 ~ "Okhotsk",
        TRUE ~ "KuroExt"
      ),
      area2 = case_when(
        long <= 139 & lat <= 34.75 ~ "Kuroshio",
        long <= 141 & lat >= 34.75 ~ "SoJ",
        long <= 142 & lat >= 44 ~ "SoJ",
        long > 140 & long <= 144 & lat <= 43 & lat > 40 ~ "SoJ",
        long > 142 & long < 146 & lat >= 44 ~ "Okhotsk",
        long >= 150 | lat > 40 ~ "KuroExt",
        TRUE ~ "Kuroshio"
      ),
      area1 = if_else(cruise_station %in% c("HK1805_P1", "HK1805_P2"), "KuroExt", area1),
      area2 = if_else(cruise_station %in% c("HK1805_P1", "HK1805_P2"), "KuroExt", area2)
    )
}

#' Perform spatial grid interpolation for a specific area and variable
interpolate_grid <- function(df, var_name, area_col, area_val, res = 0.1) {
  sub_df <- df %>%
    filter(.data[[area_col]] == area_val, !is.na(.data[[var_name]]))
  
  if (nrow(sub_df) < 3) return(NULL)
  
  min_long <- floor(min(sub_df$long))
  max_long <- ceiling(max(sub_df$long))
  min_lat  <- floor(min(sub_df$lat))
  max_lat  <- ceiling(max(sub_df$lat))
  xo_seq   <- seq(min_long, max_long, by = res)
  yo_seq   <- seq(min_lat, max_lat, by = res)
  
  akima::interp(
    x = sub_df$long, y = sub_df$lat, z = sub_df[[var_name]],
    duplicate = "mean", xo = xo_seq, yo = yo_seq
  ) %>%
    akima::interp2xyz() %>%
    as_tibble() %>%
    na.omit() %>%
    rename(long = x, lat = y, !!var_name := z)
}

#' Combine spatial grid interpolations across all defined oceanic areas
interpolate_variable <- function(df, var_name, res = 0.1) {
  areas_1 <- c("SoJ", "Okhotsk", "KuroExt")
  res_list <- list()
  for (a in areas_1) {
    res_list[[a]] <- interpolate_grid(df, var_name, "area1", a, res = res)
  }
  res_list[["Kuroshio"]] <- interpolate_grid(df, var_name, "area2", "Kuroshio", res = res)
  
  bind_rows(res_list) %>%
    mutate(long = round(long, 1), lat = round(lat, 1)) %>%
    group_by(long, lat) %>%
    summarise(!!var_name := mean(.data[[var_name]]), .groups = "drop")
}

#' Generate inverse distance weighted (IDW) raster grid with degree-based buffer masking
get_idw_grid <- function(data, target_var, buffer_deg = 2.5, dx = 0.25, dy = 0.25) {
  pts_sf <- st_as_sf(data, coords = c("long", "lat"), crs = 4326)
  
  suppressMessages(sf_use_s2(FALSE))
  buf_sf <- st_buffer(pts_sf, dist = buffer_deg) %>% st_union() %>% st_make_valid()
  suppressMessages(sf_use_s2(TRUE))
  
  grid_bbox  <- st_bbox(c(xmin = 120, xmax = 200, ymin = 20, ymax = 50), crs = 4326)
  grid_stars <- st_as_stars(grid_bbox, dx = dx, dy = dy)
  
  idw_formula <- as.formula(paste(target_var, "~ 1"))
  idw_model   <- idw(idw_formula, pts_sf, grid_stars)
  
  r_idw    <- rast(idw_model["var1.pred"])
  r_masked <- mask(r_idw, vect(buf_sf))
  
  df_grid <- as.data.frame(r_masked, xy = TRUE)
  colnames(df_grid) <- c("long", "lat", "val")
  return(df_grid)
}

#' Create standardized spatial panel ggplot component using IDW grids
create_spatial_panel <- function(data, fill_var, title_text, limits_val = NULL, breaks_val = waiver()) {
  df_grid <- get_idw_grid(data, fill_var, buffer_deg = 2.5)
  
  p <- ggplot() +
    rasterise(geom_raster(data = df_grid, aes(x = long, y = lat, fill = val)), dpi = 300) +
    scale_fill_gradientn(
      colours = colors_rev,
      limits  = limits_val,
      breaks  = breaks_val,
      oob     = scales::squish,
      guide   = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 8,
        barheight = 0.8,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    geom_point(data = data, aes(x = long, y = lat), color = "black", size = 0.8, alpha = 0.6) +
    geom_base_map + 
    coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = FALSE) +
    theme_pubr() +
    labs(x = "", y = "", title = title_text, fill = "Effect") +
    theme(
      plot.title   = element_text(size = 10, face = "bold"),
      axis.text    = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text  = element_text(size = 8)
    )
  return(p)
}


# ------------------------------------------------------------------------------
# 3. Data Ingestion & Preprocessing
# ------------------------------------------------------------------------------

# --- 3.1 POM Isotopic Database Loading & Outlier Removal ---
file_pom_main  <- "data/SI_POM_Database_230511.xlsx"
file_pom_saury <- "data/SI_POM_saury.xlsx"
file_pom_csv   <- "data/SI_POM_Database_211210.csv"

if (file.exists(file_pom_main) && file.exists(file_pom_saury)) {
  pom_SI_main  <- read_xlsx(file_pom_main, sheet = "database")
  pom_SI_saury <- read_xlsx(file_pom_saury, sheet = "database")
  
  pom_SI_saury$station <- as.character(pom_SI_saury$station)
  pom_SI_main$sal      <- as.numeric(pom_SI_main$sal)
  
  pom_SI_combined <- bind_rows(pom_SI_main, pom_SI_saury)
}

pom_surface <- read_csv(file_pom_csv) %>%
  filter(depth <= 10) %>%
  assign_oceanic_areas()

d13C_outlier_info <- remove_outliers(pom_surface$d13C)

pom_surface_rev <- pom_surface %>%
  filter(
    !d13C %in% d13C_outlier_info$del.val,
    !is.na(d13C), !is.na(sal), !is.na(temp),
    !is.na(chla), !is.na(nitrate)
  ) %>%
  mutate(
    nitrate   = if_else(nitrate < 0.01, 0.01, nitrate),
    phosphate = if_else(phosphate < 0.01, 0.01, phosphate),
    silicate  = if_else(silicate < 0.01, 0.01, silicate),
    chla      = if_else(chla < 0.01, 0.01, chla)
  )

# Summer subset (June–September) and regional classification
pom_surface_summer_rev <- pom_surface_rev %>%
  filter(month %in% 6:9) %>%
  mutate(
    area1 = if_else(station == "SJ01", "KuroExt", area1),
    area2 = if_else(station == "SJ01", "Kuroshio", area2),
    area3 = case_when(
      station == "SJ01" ~ "Kuro+Ext",
      station %in% c("P1", "P2") ~ "KuroExt",
      area1 == area2 ~ area1,
      TRUE ~ "Kuro+Ext"
    )
  )

# --- 3.2 Spatial Interpolation of POM Variables ---
vars_to_interpolate <- c("d15N", "temp", "sal", "chla", "nitrate")
interp_results <- purrr::map(vars_to_interpolate, ~ interpolate_variable(pom_surface_summer_rev, .x))

interrep_pom_surface_summer_rev <- reduce(interp_results, left_join, by = c("long", "lat")) %>%
  mutate(nitrate = if_else(nitrate < 0.01, 0.01, nitrate))

# --- 3.3 Fish Isotopic & Oceanographic Data Merging ---
fish_SI <- read_csv("./data/fish_SI.csv") 
fish_SI$longitude <- round(fish_SI$LonDeg + fish_SI$LonMin / 60, 1)
fish_SI$longitude <- if_else(fish_SI$LonCode == "W", 360 - fish_SI$longitude, fish_SI$longitude)
fish_SI$latitude  <- round(fish_SI$LatDeg + fish_SI$LatMin / 60, 1)
fish_SI$date      <- as.Date(paste(fish_SI$Year, fish_SI$Month, if_else(is.na(fish_SI$Day), 15, fish_SI$Day), sep = "-")) %>%
  as.POSIXct.Date()

glo_phys    <- read_csv("./data/points_with_temp_sal.csv")
glo_obs_chl <- read_csv("./data/points_with_CHL_data.csv")
glo_no3     <- read_csv("./data/points_with_no3_data.csv")

fish_SI <- fish_SI %>%
  left_join(glo_phys, by = c("date", "latitude", "longitude")) %>%
  left_join(glo_obs_chl, by = c("date", "latitude", "longitude")) %>%
  left_join(glo_no3, by = c("date", "latitude", "longitude")) %>%
  left_join(interrep_pom_surface_summer_rev %>% select(latitude = lat, longitude = long, d15NPOM = d15N),
            by = c("latitude", "longitude"))

fish_SI_dropna <- fish_SI %>%
  select(
    date, lat = latitude, long = longitude, BodyLength,
    temp = thetao, sal = so, chla = CHL, nitrate = no3,
    d15N, d15NPOM, Genus
  ) %>%
  drop_na()

# --- 3.4 Nitrate Isotope Database Loading ---
nitrate_pangea   <- read_csv("data/nitrate_isotope_PANGEA.csv", skip = 68)
nitrate_wang     <- read_xlsx("data/nitrate_isotope_add.xlsx", sheet = "Wang")
nitrate_horii    <- read_xlsx("data/nitrate_isotope_add.xlsx", sheet = "Horii")
nitrate_marconi  <- read_csv("data/928002_v1_gp15_nitrate_d15n_d18o_leg1.csv")
nitrate_umezawa  <- read_xlsx("data/nitrate_isotope_Japan.xlsx", sheet = "Umezawa")
nitrate_sugimoto <- read_xlsx("data/nitrate_isotope_Japan.xlsx", sheet = "Sugimoto")

nitrate_north_pacific <- bind_rows(
  nitrate_pangea %>% select(lat, long, depth, d15N_NO3),
  nitrate_wang   %>% select(lat = Latitude, long = Longitude, depth = Depth, d15N_NO3 = d15NO3) %>% mutate(d15N_NO3 = as.numeric(d15N_NO3)),
  nitrate_horii  %>% select(lat = Latitude, long = Longitude, depth = Depth, d15N_NO3 = d15NO3) %>% mutate(d15N_NO3 = as.numeric(d15N_NO3)),
  nitrate_marconi%>% select(lat = Start_Latitude, long = Start_Longitude, depth = Sample_Depth, d15N_NO3 = NITRATE_15_14_D_DELTA_BOTTLE_lnw6zq),
  nitrate_umezawa%>% select(lat = Latitude, long = Longitude, depth = Layer, d15N_NO3 = d15N),
  nitrate_sugimoto%>% select(lat, long, depth = Depth, d15N_NO3 = d15NO3)
) %>%
  mutate(long = if_else(long < 0, 360 + long, long)) %>%
  filter(depth <= 200, long >= 100, long <= 240, lat >= 0, lat <= 70)


# ------------------------------------------------------------------------------
# 4. Stage 1: Bayesian GAM for POM d15N (brms)
# ------------------------------------------------------------------------------

pom_data <- pom_surface_summer_rev %>%
  mutate(log10_chla = log10(chla))

# --- Model Selection using WAIC ---
fit1 <- brm(d15N ~ s(long, lat) + temp + sal + nitrate + log10_chla, data = pom_data, family = gaussian(), chains = 4, cores = 4, iter = 2000, seed = 123, backend = "cmdstanr")
fit2 <- brm(d15N ~ s(long, lat) + temp + nitrate + log10_chla,       data = pom_data, family = gaussian(), chains = 4, cores = 4, iter = 2000, seed = 123, backend = "cmdstanr")
fit3 <- brm(d15N ~ s(long, lat) + temp + sal + nitrate,              data = pom_data, family = gaussian(), chains = 4, cores = 4, iter = 2000, seed = 123, backend = "cmdstanr")
fit4 <- brm(d15N ~ s(long, lat) + temp + nitrate,                    data = pom_data, family = gaussian(), chains = 4, cores = 4, iter = 2000, seed = 123, backend = "cmdstanr")

fit1 <- add_criterion(fit1, "waic")
fit2 <- add_criterion(fit2, "waic")
fit3 <- add_criterion(fit3, "waic")
fit4 <- add_criterion(fit4, "waic")

loo_compare(fit1, fit2, fit3, fit4, criterion = "waic")

# --- Final Stage 1 Bayesian GAM Model ---
fit_brms_pom <- fit1

posterior_summary(fit_brms_pom, variable = c("b_temp", "b_sal", "b_nitrate", "b_log10_chla", "sigma"))
bayes_R2(fit_brms_pom)


# ------------------------------------------------------------------------------
# 5. Figure 01: Sampling Map, Observed POM Isoscape & Bayesian GAM Diagnostics
# ------------------------------------------------------------------------------

# --- Figure 01a: Ocean Currents & Sampling Locations ---
current <- tidync("./data/cmems_mod_glo_phy_my_0.083deg-climatology_P1M-m_1781494578834.nc") %>%
  hyper_tibble() %>% 
  group_by(latitude, longitude) %>%
  summarise(uo = mean(uo), vo = mean(vo), .groups = "drop") %>%
  mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

current01 <- current %>%
  na.omit() %>%
  mutate(lat = round(latitude * 2, 0) / 2, long = round(longitude * 2, 0) / 2) %>% 
  group_by(lat, long) %>%
  summarise(uo = mean(uo), vo = mean(vo), .groups = "drop") %>%
  mutate(velocity = sqrt(vo^2 + uo^2))

current01$uo[current01$velocity < 0.1] <- 0
current01$vo[current01$velocity < 0.1] <- 0

current01 <- current01 %>%
  complete(long = seq(min(long), max(long), by = 0.5), lat = seq(min(lat), max(lat), by = 0.5),
           fill = list(uo = 0, vo = 0, velocity = 0))

figure01a <- ggplot(pom_surface_summer_rev, aes(x = long, y = lat)) +
  geom_point(aes(color = area3), shape = 16) +
  geom_point(data = fish_SI_dropna %>% select(long, lat) %>% unique(), aes(x = long, y = lat), shape = 4) +
  geom_streamline(data = current01, aes(x = long, y = lat, dx = uo, dy = vo),
                  L = 100, res = 10, arrow.length = 0.3, arrow.type = "open", arrow.angle = 30,
                  color = "black", n = 12, alpha = 0.4) +
  geom_base_map +
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = 0) +
  scale_color_atlassian() +
  theme_pubr() + 
  xlab("Longitude (°E)") + ylab("Latitude (°N)")

ggsave("output/figure01a.pdf", figure01a, width = 160, height = 100, units = "mm")

# --- Figure 01b: Surface POM d15N Isoscape Map ---
fig01_d15N <- ggplot(interrep_pom_surface_summer_rev) +
  rasterise(geom_point(aes(x = long, y = lat, color = d15N), shape = 15, size = 1, alpha = 1), dpi = 300) +
  geom_point(data = pom_surface_summer_rev, aes(x = long, y = lat), size = 0.2) +
  scale_color_gradientn(breaks = seq(-2, 10, 2), colours = colors_rev) +
  geom_base_map +
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = 0) +
  theme_pubr(legend = "bottom") +
  xlab("Longitude (°E)") + ylab("Latitude (°N)") +
  guides(color = guide_colorbar(barwidth = 10, barheight = 1))

ggsave("output/Figure01b.pdf", fig01_d15N, width = 18, height = 18, units = "cm")

# --- Figure 01c-g: Partial Effects & Spatial Spline Surface ---
ce_temp    <- conditional_effects(fit_brms_pom, "temp")$temp
ce_sal     <- conditional_effects(fit_brms_pom, "sal")$sal
ce_nitrate <- conditional_effects(fit_brms_pom, "nitrate")$nitrate
ce_chla    <- conditional_effects(fit_brms_pom, "log10_chla")$log10_chla

p_temp <- ggplot(ce_temp, aes(x = temp, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "gray80", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = pom_data, aes(x = temp, y = d15N), size = 0.5, alpha = 0.4) +
  labs(x = "Temperature (°C)", y = "Partial effect") + theme_pubr()

p_sal <- ggplot(ce_sal, aes(x = sal, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "gray80", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = pom_data, aes(x = sal, y = d15N), size = 0.5, alpha = 0.4) +
  labs(x = "Salinity", y = "Partial effect") + theme_pubr()

p_nitrate <- ggplot(ce_nitrate, aes(x = nitrate, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "gray80", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = pom_data, aes(x = nitrate, y = d15N), size = 0.5, alpha = 0.4) +
  labs(x = "Nitrate (µM)", y = "Partial effect") + theme_pubr()

p_chla <- ggplot(ce_chla, aes(x = 10^log10_chla, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "gray80", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(data = pom_data, aes(x = chla, y = d15N), size = 0.5, alpha = 0.4) +
  scale_x_continuous(trans = "log10") +
  labs(x = "Chl-a (µg/L)", y = "Partial effect") + theme_pubr()

fig01_partial <- (p_temp + p_sal) / (p_nitrate + p_chla)

# Spatial spline raster with 2-degree buffer masking
pts_sf <- st_as_sf(pom_data, coords = c("long", "lat"), crs = 4326)

suppressMessages(sf_use_s2(FALSE))
buf_sf <- st_buffer(pts_sf, dist = 2) %>% st_union() %>% st_make_valid()
suppressMessages(sf_use_s2(TRUE))

world_360 <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_shift_longitude() %>%
  st_make_valid()

grid_spatial <- expand.grid(
  long       = seq(120, 200, by = 0.25),
  lat        = seq(20, 50, by = 0.25),
  temp       = mean(pom_data$temp, na.rm = TRUE),
  sal        = mean(pom_data$sal, na.rm = TRUE),
  nitrate    = mean(pom_data$nitrate, na.rm = TRUE),
  log10_chla = mean(pom_data$log10_chla, na.rm = TRUE)
)

pred_grid <- fitted(fit_brms_pom, newdata = grid_spatial)
grid_spatial$z <- pred_grid[, "Estimate"]

r_spatial <- rast(grid_spatial %>% select(long, lat, z), type = "xyz", crs = "EPSG:4326")
r_spatial_masked <- mask(r_spatial, vect(buf_sf))
r_spatial_masked <- mask(r_spatial_masked, vect(world_360), inverse = TRUE)

grid_spatial_masked <- as.data.frame(r_spatial_masked, xy = TRUE)
colnames(grid_spatial_masked) <- c("long", "lat", "z")

fig01_spatial <- ggplot() +
  geom_raster(data = grid_spatial_masked, aes(x = long, y = lat, fill = z)) +
  scale_fill_gradientn(colours = colors_rev, na.value = NA) +
  geom_point(data = pom_data, aes(x = long, y = lat), size = 0.8, color = "black") +
  geom_base_map +
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = FALSE) +
  theme_pubr(legend = "bottom") +
  labs(fill = "Spline", x = "Longitude (°E)", y = "Latitude (°N)")

Figure01c_g <- fig01_spatial / fig01_partial
ggsave("output/Figure01c_g.pdf", Figure01c_g, width = 18, height = 24, units = "cm")


# ------------------------------------------------------------------------------
# 6. Stage 2: Hierarchical Bayesian Fish Model (cmdstanr) & Pipeline Integration
# ------------------------------------------------------------------------------

# --- 6.1 Construct Master Site Reference & Standardize Environmental Data ---
master_sites <- bind_rows(
  pom_surface_summer_rev %>% select(lat, long),
  fish_SI_dropna %>% select(lat, long)
) %>%
  mutate(
    lat_round  = round(lat, 1),
    long_round = round(long, 1),
    site_key   = paste(lat_round, long_round, sep = "_")
  ) %>%
  distinct(site_key, .keep_all = TRUE) %>%
  arrange(site_key) %>%
  mutate(site_id = row_number())

N_site_total <- nrow(master_sites)

all_env_data <- bind_rows(
  pom_surface_summer_rev %>% select(lat, long, temp, nitrate, chla, sal),
  fish_SI_dropna %>% select(lat, long, temp, nitrate, chla, sal)
) %>%
  mutate(
    lat_round  = round(lat, 1),
    long_round = round(long, 1),
    site_key   = paste(lat_round, long_round, sep = "_")
  ) %>%
  group_by(site_key) %>%
  summarise(
    temp    = mean(temp, na.rm = TRUE),
    nitrate = mean(nitrate, na.rm = TRUE),
    chla    = mean(log(chla), na.rm = TRUE),
    sal     = mean(sal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    temp_s    = as.vector(scale(temp)),
    nitrate_s = as.vector(scale(nitrate)),
    chla_s    = as.vector(scale(chla)),
    sal_s     = as.vector(scale(sal))
  )

stan_data_pom_df <- pom_surface_summer_rev %>%
  filter(nitrate > 0.01, !is.na(d15N)) %>%
  mutate(
    lat_round  = round(lat, 1),
    long_round = round(long, 1),
    site_key   = paste(lat_round, long_round, sep = "_")
  ) %>%
  left_join(master_sites %>% select(site_key, site_id), by = "site_key") %>%
  left_join(all_env_data %>% select(site_key, temp_s, chla_s, sal_s, nitrate_s), by = "site_key")

# --- 6.2 Predict Latent Baseline POM d15N via Stage 1 brms Model ---
df_fish <- fish_SI_dropna %>%
  mutate(
    species_id = as.integer(as.factor(Genus)),
    BL         = BodyLength,
    log10_chla = log10(chla),
    lat_round  = round(lat, 1),
    long_round = round(long, 1),
    site_key   = paste(lat_round, long_round, sep = "_")
  ) %>%
  select(-any_of(c("site_id", "temp_s", "nitrate_s", "chla_s", "sal_s"))) %>%
  left_join(master_sites %>% select(site_key, site_id), by = "site_key") %>%
  left_join(all_env_data %>% select(site_key, temp_s, nitrate_s, chla_s, sal_s), by = "site_key")

pom_pred_brms <- fitted(fit_brms_pom, newdata = df_fish, re_formula = NA)
df_fish <- df_fish %>% mutate(d15NPOM_model1 = pom_pred_brms[, "Estimate"])

site_env <- master_sites %>%
  select(site_id, site_key) %>%
  left_join(all_env_data, by = "site_key") %>%
  arrange(site_id)

# --- 6.3 Execute Stage 2 Stan Sampling ---
stan_data_fish <- list(
  N            = nrow(df_fish),
  N_site       = N_site_total,
  d15N         = df_fish$d15N,
  POM_obs      = df_fish$d15NPOM_model1, 
  BL           = as.vector(scale(df_fish$BL)),
  temp_site    = site_env$temp_s,
  nitrate_site = site_env$nitrate_s,
  chla_site    = site_env$chla_s,
  sal_site     = site_env$sal_s,
  site         = df_fish$site_id,
  species      = df_fish$species_id,
  N_species    = length(unique(df_fish$species_id))
)

mod_fish <- cmdstan_model("./models/model_stage3_final.stan")

fit_fish <- mod_fish$sample(
  data            = stan_data_fish,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1500,
  iter_sampling   = 1500,
  adapt_delta     = 0.95,
  seed            = 123,
  init            = 0.5
)

fit_fish$summary(c(
  "beta_pom", "beta_chla", "beta_temp",
  "beta_nitrate", "beta_sal", "beta_BL", "sigma_pom"
))

# Extract Posterior Effect Means
beta_temp    <- fit_fish$summary("beta_temp")$mean
beta_nitrate <- fit_fish$summary("beta_nitrate")$mean
beta_chla    <- fit_fish$summary("beta_chla")$mean
beta_sal     <- fit_fish$summary("beta_sal")$mean
beta_pom     <- fit_fish$summary("beta_pom")$mean


# ------------------------------------------------------------------------------
# 7. Spatial Variance Decomposition & Figure 02 Generation
# ------------------------------------------------------------------------------

# --- 7.1 Site-Level Spatial Variance Decomposition ---
draws      <- as_draws_df(fit_fish$draws())
pom_cols   <- grep("^POM_true\\[", names(draws), value = TRUE)
POM_matrix <- as.matrix(draws[, pom_cols])
n_draws    <- nrow(POM_matrix)

contrib <- data.frame(
  baseline    = numeric(n_draws),
  chla        = numeric(n_draws),
  nitrate     = numeric(n_draws),
  temperature = numeric(n_draws),
  salinity    = numeric(n_draws)
)

for (i in 1:n_draws) {
  base_val  <- draws$beta_pom[i] * POM_matrix[i, ]
  chla_val  <- draws$beta_chla[i] * stan_data_fish$chla_site
  nit_val   <- draws$beta_nitrate[i] * stan_data_fish$nitrate_site
  temp_val  <- draws$beta_temp[i] * stan_data_fish$temp_site
  sal_val   <- draws$beta_sal[i] * stan_data_fish$sal_site
  
  total_val <- base_val + chla_val + nit_val + temp_val + sal_val
  
  contrib$baseline[i]    <- var(base_val) / var(total_val)
  contrib$chla[i]        <- var(chla_val) / var(total_val)
  contrib$nitrate[i]     <- var(nit_val) / var(total_val)
  contrib$temperature[i] <- var(temp_val) / var(total_val)
  contrib$salinity[i]    <- var(sal_val) / var(total_val)
}

summary_contrib <- contrib %>%
  summarise(
    across(everything(), list(mean = mean, low = ~quantile(.x, 0.025), high = ~quantile(.x, 0.975)))
  ) %>%
  pivot_longer(everything(), names_to = c("variable", ".value"), names_sep = "_")

plot_df <- summary_contrib %>%
  mutate(
    variable = factor(
      variable,
      levels = c("baseline", "chla", "nitrate", "temperature", "salinity"),
      labels = c("Baseline δ15NPOM", "Chlorophyll-a", "Nitrate", "Temperature", "Salinity")
    )
  )

p_variance_decomp <- ggplot(plot_df, aes(x = mean, y = variable)) +
  geom_errorbar(aes(xmin = low, xmax = high), orientation = "y", height = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5) +
  theme_bw(base_size = 14) +
  labs(x = "Proportion of site-level spatial variance", y = NULL) +
  theme(panel.grid.minor = element_blank())

print(p_variance_decomp)

# --- 7.2 Figure 02 Spatial Effects Panels ---
points_sf  <- st_as_sf(df_fish, coords = c("long", "lat"), crs = 4326)
grid_bbox  <- st_bbox(c(xmin = 120, xmax = 200, ymin = 20, ymax = 50), crs = 4326)
grid_stars <- st_as_stars(grid_bbox, dx = 0.25, dy = 0.25)
idw_model  <- idw(d15N ~ 1, points_sf, grid_stars)
r_idw      <- rast(idw_model["var1.pred"])

suppressMessages(sf_use_s2(FALSE))
buffer_sf  <- st_buffer(points_sf, dist = 3) %>% st_union() %>% st_make_valid()
suppressMessages(sf_use_s2(TRUE))

r_masked   <- mask(r_idw, vect(buffer_sf))
df_grid    <- as.data.frame(r_masked, xy = TRUE)
colnames(df_grid) <- c("long", "lat", "d15N")

p_d15Nfish <- ggplot() +
  geom_raster(data = df_grid, aes(x = long, y = lat, fill = d15N)) +
  geom_point(data = df_fish, aes(x = long, y = lat), color = "black", size = 0.5) +
  scale_fill_gradientn(
    colours = colors_rev,
    name = expression(delta^15*N[fish]~"(‰)"),
    limits = c(6, 14),
    oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = 16, barheight = 0.8,
      frame.colour = "black", ticks.colour = "black"
    )
  ) +
  geom_base_map +
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = FALSE) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)", tag = "(a)") +
  theme_pubr() +
  theme(
    plot.title   = element_text(size = 10, face = "bold"),
    axis.text    = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8)
  )

df_env <- df_fish %>% 
  select(long, lat, d15NPOM = d15NPOM_model1, chla_s, temp_s, sal_s, nitrate_s) %>% 
  mutate(
    env_temp             = beta_temp * temp_s,
    env_nitrate          = beta_nitrate * nitrate_s,
    env_chla             = beta_chla * chla_s,
    env_sal              = beta_sal * sal_s,
    environmental_effect = env_temp + env_nitrate + env_chla + env_sal,
    baseline             = d15NPOM * beta_pom 
  )

p_base    <- create_spatial_panel(df_env, "baseline", "Baseline Effect")
p_env     <- create_spatial_panel(df_env, "environmental_effect", "Total Environmental Effect", limits_val = c(-2.5, 2.5), breaks_val = seq(-2, 2, 1))
p_temp    <- create_spatial_panel(df_env, "env_temp", "Temperature", limits_val = c(-2.5, 2.5), breaks_val = seq(-2, 2, 1))
p_nitrate <- create_spatial_panel(df_env, "env_nitrate", "Nitrate", limits_val = c(-2.5, 2.5), breaks_val = seq(-2, 2, 1))
p_chla    <- create_spatial_panel(df_env, "env_chla", "Chlorophyll", limits_val = c(-2.5, 2.5), breaks_val = seq(-2, 2, 1))
p_sal     <- create_spatial_panel(df_env, "env_sal", "Salinity", limits_val = c(-2.5, 2.5), breaks_val = seq(-2, 2, 1))

design_layout <- "
  AA
  AA
  BB
  BB
  CC
  CC
  DE
  FG
"

Figure_environmental_effects <- wrap_plots(
  A = p_d15Nfish, B = p_base, C = p_env, D = p_chla, E = p_temp, F = p_nitrate, G = p_sal,
  design = design_layout
)

ggsave("./output/figure02.pdf", Figure_environmental_effects, width = 21, height = 60, units = "cm")


# ------------------------------------------------------------------------------
# 8. Supporting Information Figures (Figure S01 - S04)
# ------------------------------------------------------------------------------

# --- Figure S01: Environmental Property Maps ---
S01_p_temp <- ggplot(interrep_pom_surface_summer_rev) +
  rasterise(geom_point(aes(x = long, y = lat, color = temp), shape = 15, size = 1, alpha = 0.8), dpi = 300) +
  geom_point(data = pom_surface_summer_rev, aes(x = long, y = lat), size = 0.2) +
  scale_color_stepsn(colours = colors_rev, breaks = seq(5, 29, 2), limits = c(5, 29), oob = scales::squish) +
  geom_base_map +
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = 0) +
  theme_pubr(legend = "bottom") +
  guides(color = guide_colorbar(barwidth = 10, barheight = 1)) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)")

S01_p_sal <- ggplot(interrep_pom_surface_summer_rev) +
  rasterise(geom_point(aes(x = long, y = lat, color = sal), shape = 15, size = 1, alpha = 0.8), dpi = 300) +
  geom_point(data = pom_surface_summer_rev, aes(x = long, y = lat), size = 0.2) +
  scale_color_stepsn(colours = colors_rev, breaks = seq(31.5, 34.5, 0.25), limits = c(32, 34.5), oob = scales::squish) +
  geom_base_map +
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = 0) +
  theme_pubr(legend = "bottom") +
  guides(color = guide_colorbar(barwidth = 10, barheight = 1)) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)")

S01_p_chla <- ggplot(interrep_pom_surface_summer_rev) +
  rasterise(geom_point(aes(x = long, y = lat, color = chla), shape = 15, size = 1, alpha = 0.8), dpi = 300) +
  geom_point(data = pom_surface_summer_rev, aes(x = long, y = lat), size = 0.2) +
  scale_color_stepsn(colours = colors_rev, trans = "log10", limits = c(0.1, 5), oob = scales::squish, n.breaks = 10) +
  geom_base_map +
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = 0) +
  theme_pubr(legend = "bottom") +
  guides(color = guide_colorbar(barwidth = 10, barheight = 1)) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)")

S01_p_nitrate <- ggplot(interrep_pom_surface_summer_rev) +
  rasterise(geom_point(aes(x = long, y = lat, color = nitrate), shape = 15, size = 1, alpha = 0.8), dpi = 300) +
  geom_point(data = pom_surface_summer_rev, aes(x = long, y = lat), size = 0.2) +
  scale_color_stepsn(colours = colors_rev, trans = "log10", limits = c(0.01, 30), oob = scales::squish, n.breaks = 10) +
  geom_base_map +
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = 0) +
  theme_pubr(legend = "bottom") +
  guides(color = guide_colorbar(barwidth = 10, barheight = 1)) +
  xlab("Longitude (°E)") + ylab("Latitude (°N)")

FigureS01 <- S01_p_temp / S01_p_sal / S01_p_chla / S01_p_nitrate
ggsave("output/FigureS01.pdf", FigureS01, width = 16, height = 28, units = "cm")

# --- Figure S02: Temperature-Salinity Relationships vs Isotopic/Nutrient Properties ---
FigureS02 <- ggplot(pom_surface_summer_rev, aes(x = sal, y = temp)) +
  geom_point(size = 0.5, color = "gray") +
  geom_point(data = pom_surface_summer_rev %>% filter(area3 %in% c("Kuroshio", "KuroExt", "Kuro+Ext")), size = 2, aes(color = d15N)) +
  scale_color_gradientn(colours = colors_rev, limits = c(-3, 10)) +
  theme_pubr() + theme(legend.position = "right") + guides(color = guide_colorbar(barwidth = 1, barheight = 5)) +
  ggplot(pom_surface_summer_rev, aes(x = sal, y = temp)) +
  geom_point(size = 0.5, color = "gray") +
  geom_point(data = pom_surface_summer_rev %>% filter(!area3 %in% c("Kuroshio", "KuroExt", "Kuro+Ext")), size = 2, aes(color = d15N)) +
  scale_color_gradientn(colours = colors_rev, limits = c(-3, 10)) +
  theme_pubr() + theme(legend.position = "right") + guides(color = guide_colorbar(barwidth = 1, barheight = 5)) +
  ggplot(pom_surface_summer_rev, aes(x = sal, y = temp)) +
  geom_point(size = 0.5, color = "gray") +
  geom_point(data = pom_surface_summer_rev %>% filter(area3 %in% c("Kuroshio", "KuroExt", "Kuro+Ext")), size = 2, aes(color = chla)) +
  scale_color_gradientn(colours = colors_rev, trans = "log10", limits = c(0.01, 5)) +
  theme_pubr() + theme(legend.position = "right") + guides(color = guide_colorbar(barwidth = 1, barheight = 5)) +
  ggplot(pom_surface_summer_rev, aes(x = sal, y = temp)) +
  geom_point(size = 0.5, color = "gray") +
  geom_point(data = pom_surface_summer_rev %>% filter(!area3 %in% c("Kuroshio", "KuroExt", "Kuro+Ext")), size = 2, aes(color = chla)) +
  scale_color_gradientn(colours = colors_rev, trans = "log10", limits = c(0.01, 5)) +
  theme_pubr() + theme(legend.position = "right") + guides(color = guide_colorbar(barwidth = 1, barheight = 5)) +
  ggplot(pom_surface_summer_rev, aes(x = sal, y = temp)) +
  geom_point(size = 0.5, color = "gray") +
  geom_point(data = pom_surface_summer_rev %>% filter(area3 %in% c("Kuroshio", "KuroExt", "Kuro+Ext")), size = 2, aes(color = nitrate)) +
  scale_color_gradientn(colours = colors_rev, trans = "log10", limits = c(0.01, 10), oob = scales::squish) +
  theme_pubr() + theme(legend.position = "right") + guides(color = guide_colorbar(barwidth = 1, barheight = 5)) +
  ggplot(pom_surface_summer_rev, aes(x = sal, y = temp)) +
  geom_point(size = 0.5, color = "gray") +
  geom_point(data = pom_surface_summer_rev %>% filter(!area3 %in% c("Kuroshio", "KuroExt", "Kuro+Ext")), size = 2, aes(color = nitrate)) +
  scale_color_gradientn(colours = colors_rev, trans = "log10", limits = c(0.01, 10), oob = scales::squish) +
  theme_pubr() + theme(legend.position = "right") + guides(color = guide_colorbar(barwidth = 1, barheight = 5)) +
  plot_layout(ncol = 3, byrow = FALSE)

ggsave("./output/FigureS02.pdf", plot = FigureS02, width = 300, height = 150, units = "mm")

# --- Figure S03: Frequentist GAM Analysis for Fish d15N ---
fish_d15N_gam <- gam(d15N ~ s(long, lat) + d15NPOM + BodyLength + temp + sal + log10(chla) + nitrate,
                     data = fish_SI_dropna)
dredge(fish_d15N_gam)

fish_d15N_gam <- gam(d15N ~ s(long, lat) + d15NPOM + BodyLength + temp + sal + nitrate,
                     data = fish_SI_dropna)

fish_d15N_gam_viz <- getViz(fish_d15N_gam)
summary(fish_d15N_gam)

figureS03A <- plot(sm(fish_d15N_gam_viz, 1))$ggObj$data %>%
  na.omit() %>%
  ggplot(aes(x = x, y = y)) +
  geom_tile(aes(fill = z)) +
  scale_fill_gradientn(colours = colors_rev) +
  geom_point(data = fish_SI_dropna, aes(x = long, y = lat), size = 2) + 
  geom_point(data = pom_surface_summer_rev, aes(x = long, y = lat), size = 0.5, shape = 1) + 
  geom_base_map + 
  coord_sf(xlim = c(120, 200), ylim = c(20, 50), expand = 0) +
  theme_pubr(legend = "bottom") +
  guides(color = guide_colorbar(barwidth = 5, barheight = 1)) +
  labs(fill = "spline", x = "Longitude (°E)", y = "Latitude (°N)")

figureS03B <- ggplot(data = plot(pterm(fish_d15N_gam_viz, 1))$data$res, aes(x = x, y = y)) +
  geom_point(size = 0.5) + stat_smooth(method = "lm") +
  theme_pubr(legend = "right") +
  labs(x = "POM d15N (‰)", y = "Effect on Fish d15N (‰)")+ 
  theme(
    aspect.ratio = 1,
    axis.text    = element_text(size = 8),
    axis.title   = element_text(size = 9)
  )

figureS03C <- ggplot(data = plot(pterm(fish_d15N_gam_viz, 2))$data$res, aes(x = x, y = y)) +
  geom_point(size = 0.5) + stat_smooth(method = "lm") +
  theme_pubr(legend = "right") +
  labs(x = "Body Length (mm)", y = "Effect on Fish d15N (‰)")+ 
  theme(
    aspect.ratio = 1,
    axis.text    = element_text(size = 8),
    axis.title   = element_text(size = 9)
  )

figureS03D <- ggplot(data = plot(pterm(fish_d15N_gam_viz, 3))$data$res, aes(x = x, y = y)) +
  geom_point(size = 0.5) + stat_smooth(method = "lm") +
  theme_pubr(legend = "right") +
  labs(x = "Temperature (°C)", y = "Effect on Fish d15N (‰)")+ 
  theme(
    aspect.ratio = 1,
    axis.text    = element_text(size = 8),
    axis.title   = element_text(size = 9)
  )

figureS03E <- ggplot(data = plot(pterm(fish_d15N_gam_viz, 4))$data$res, aes(x = x, y = y)) +
  geom_point(size = 0.5) + stat_smooth(method = "lm") +
  theme_pubr(legend = "right") +
  labs(x = "Salinity", y = "Effect on Fish d15N (‰)")+ 
  theme(
    aspect.ratio = 1,
    axis.text    = element_text(size = 8),
    axis.title   = element_text(size = 9)
  )

figureS03F <- ggplot(data = plot(pterm(fish_d15N_gam_viz, 5))$data$res, aes(x = x, y = y)) +
  geom_point(size = 0.5) + stat_smooth(method = "lm") +
  theme_pubr(legend = "right") +
  labs(x = "Nitrate (µM)", y = "Effect on Fish d15N (‰)")　+ 
  theme(
    aspect.ratio = 1,
    axis.text    = element_text(size = 8),
    axis.title   = element_text(size = 9)
  )

figureS03 <- wrap_plots(
  A = figureS03A,
  B = figureS03B,
  C = figureS03C,
  D = figureS03D,
  E = figureS03E,
  F = figureS03F,
  design = design_layout_s03,
  heights = c(1.1, 1, 1) # Aの高さ割合を圧縮（1.1 : 1 : 1）
)

ggsave("./output/figureS03.pdf", plot = figureS03, width = 210, height = 210, units = "mm")
