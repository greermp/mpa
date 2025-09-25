library(ggplot2)
library(dplyr)
library(scales)
library(metR)
setwd("~/Anduril/interview/plots")

theme_set(
  theme_minimal() +
    theme(axis.title = element_text(size=12),
          axis.text = element_text(size=11),
          strip.text = element_text(size=13)
          )
)

cost_model <- function(mach, alt_kft, sensor_cost){
  return (50 * mach^2 - 35 * mach + 0.03 * alt_kft^2 - 0.02 * alt_kft^2 + 11 + sensor_cost)
}

endurance_model <- function(mach, alt_kft){
  return (-18.75 * mach^2 + 8.0893 * mach + 0.01 * alt_kft^2 + 0.05 * alt_kft + 9.2105)
}

# Example parameter ranges
altitude <- seq(0, 25,    by = .5)         # kft
mach     <- seq(0.4, 0.9, by = 0.05)       # degrees
fov_levels <- c(15, 30, 60)          # FOV degrees


sensors  <- c("Sensor 1" = 50000 / 1e6, 
              "Sensor 2" = 1e6   / 1e6,
              "Sensor 3" = 10e6  / 1e6)

# Create df of design points
df <- expand.grid(
  mach = mach, 
  alt_kft = altitude, 
  sensor = names(sensors))

# Map sensor cost to sensor name
df$sensor_cost <- sensors[df$sensor]


df$cost <- cost_model(df$mach, df$alt_kft, df$sensor_cost)
df$endurance <- endurance_model(df$mach, df$alt_kft)

###############################################################################
# Plot Endurance Model
df_endurance <- df[!duplicated(df[c("mach", "alt_kft")]), ]

df_endurance$cat = "hello" # Facet wrap var, for constant plot size
ggplot(df_endurance, aes(x = mach, y = alt_kft)) +
  geom_raster(aes(fill = endurance), interpolate = TRUE) +
  scale_fill_viridis_c(option = "turbo", direction = -1) +
  geom_contour(aes(z = endurance), color = "white", breaks = seq(4,16,2)) +
  geom_text_contour(  aes(z = endurance),                # labels along the lines
                      breaks = seq(4, 16, by = 2),      # same breaks as contour
                      color='white',
                      stroke.color='black',  
                      nudge_x = -.03,
                      stroke = 0.1,                     # white outline for contrast
                      size = 4, 
                      skip = 0                          # ensures every line gets at least one label
  ) +
  labs(x = "Speed (mach)", y = "Altitude (kft)", fill = "Endurance (hrs)") +
  scale_x_continuous(breaks = pretty_breaks()) 
  # facet_wrap(~cat)

ggsave("endurange_cont.png", height = 5, width = 6, dpi = 500)
###############################################################################

###############################################################################
# Plot Cost Model
df <- df %>% mutate(sensor = paste0(sensor, " ($", sensor_cost, "M)"))
ggplot(df, aes(x = mach, y = alt_kft)) +
  geom_raster(aes(fill=cost),interpolate = TRUE) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       breaks = seq(0, 50, by = 10 )) +
  geom_contour(aes(z = cost), color = "white", breaks = seq(10,40,5)) +
  geom_text_contour(  aes(z = cost),                # labels along the lines
    breaks = seq(10, 40, by = 5),      # same breaks as contour
    color='white',
    stroke.color='black',  
    nudge_x = -.03,
    stroke = 0.1,                     # white outline for contrast
    size = 4, 
    skip = 0                          # ensures every line gets at least one label
  ) +
  facet_wrap(~sensor) +
  labs(x = "Mach", y = "Altitude (kft)", fill = "Cost ($M)") +
  scale_x_continuous(breaks = pretty_breaks())
ggsave("cost.png", height = 5, width = 12.5, dpi = 500)

###############################################################################

# Compute footprint area (assuming square FOV for simplicity)
df <- expand.grid(
  alt_kft = altitude,
  fov_deg = fov_levels
)

# Compute footprint (assuming square FOV)
df <- df %>%
  mutate(
    fov_rad = fov_deg * pi / 180,
    footprint_km2 = (2 * alt_kft * 1000 * 0.3048 * tan(fov_rad/2))^2 / 1e6   # kft -> ft -> m -> km²
  )

# Line plot
ggplot(df, aes(x = alt_kft, y = footprint_km2, color = factor(fov_deg))) +
  geom_line(size = 1.2) +
  labs(
    x = "Altitude (kft)",
    y = "Sensor Footprint (km²)",
    color = "FOV (deg)",
    # title = "Sensor Footprint vs Altitude for Different FOVs",
    subtitle = "Assumes square sensor FOV"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = comma) +
  scale_color_brewer(palette = 'Dark2') 
ggsave("sensor_fov.png", height = 5, width = 6, dpi = 500)
