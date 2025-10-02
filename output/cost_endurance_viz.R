library(ggplot2)
library(dplyr)
library(scales)
library(metR)
setwd("~/Anduril/interview/output")

theme_set(
  theme_minimal(base_size = 14) +
    theme(
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(face = "bold")
          )
)

# Return cost in millions. Note - alt param is in thousands of feet, sensor cost in millions!
cost_model <- function(mach, alt_ft, sensor_cost_mm){
  alt_kft = alt_ft/1000
  return (50 * mach^2 - 35 * mach + 0.03 * alt_kft^2 - 0.2 * alt_kft + 11 + sensor_cost_mm)
}

# Return endurance in hours. Note - alt param is in thousands of feet!
endurance_model <- function(mach, alt_ft){
  alt_kft = alt_ft/1000
  return (-18.75 * mach^2 + 8.0893 * mach + 0.01 * alt_kft^2 + 0.05 * alt_kft + 9.2105)
}

# Example parameter ranges
altitude <-   seq(0, 25000,    by = 500)         # kft
mach     <-   seq(0.4, 0.9, by = 0.05)       # degrees
fov_levels <- c(15, 30, 60)          # FOV degrees


#
sensor_cost_mm  <- c("Sensor 1" = 50000 / 1e6, 
                     "Sensor 2" = 1e6   / 1e6,
                     "Sensor 3" = 10e6  / 1e6)

# Create df of design points
df <- expand.grid(
  mach = mach, 
  alt_ft = altitude, 
  sensor = names(sensor_cost_mm))

# Map sensor cost to sensor name
df$sensor_cost_mm <- sensor_cost_mm[df$sensor]


df$cost <- cost_model(df$mach, df$alt_ft, df$sensor_cost_mm)
df$endurance <- endurance_model(df$mach, df$alt_ft)

###############################################################################
# Plot Endurance Model
df_endurance <- df[!duplicated(df[c("mach", "alt_ft")]), ]

ggplot(df_endurance, aes(x = mach, y = alt_ft)) +
  geom_raster(aes(fill = endurance), interpolate = TRUE) +
  scale_fill_viridis_c(option = "turbo", direction = -1) +
  scale_y_continuous(labels=comma) +
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
  labs(x = "Speed (mach)", y = "Altitude (ft)", fill = "Endurance (hrs)") +
  scale_x_continuous(breaks = pretty_breaks()) 

ggsave("plots/endurange_cont.png", height = 5, width = 6.5, dpi = 500)
###############################################################################

###############################################################################
# Plot Cost Model
df <- df %>% mutate(sensor = paste0(sensor, " ($", sensor_cost_mm, "M)"))
ggplot(df, aes(x = mach, y = alt_ft)) +
  geom_raster(aes(fill=cost),interpolate = TRUE) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       breaks = seq(0, 50, by = 10 )) +
  geom_contour(aes(z = cost), color = "white", breaks = seq(10,40,5)) +
  scale_y_continuous(labels=comma) +
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
  labs(x = "Mach", y = "Altitude (ft)", fill = "Cost ($M)") +
  scale_x_continuous(breaks = pretty_breaks())
ggsave("plots/cost.png", height = 5, width = 13, dpi = 500)

###############################################################################

# Compute footprint area (assuming square FOV for simplicity)
df <- expand.grid(
  alt_ft = altitude,
  fov_deg = fov_levels
)

# Compute footprint (assuming square FOV)
df <- df %>%
  mutate(
    fov_rad = fov_deg * pi / 180,
    sensor_width_m = 2 * alt_ft * 0.3048 * tan(fov_rad / 2),
    footprint_km2 = (sensor_width_m^2) / 1e6
  )


# width
ggplot(df, aes(y = alt_ft, x = sensor_width_m, color = factor(fov_deg))) +
  geom_line(size = 1.2) +
  labs(
    y = "Altitude (ft)",
    x = "Sensor Coverage Width (m)",
    color = "FOV (deg)"
    # title = "Sensor Footprint vs Altitude for Different FOVs",
    # subtitle = "Assumes square sensor FOV"
  ) +
  scale_y_continuous(labels=comma) +
  scale_x_continuous(limits = c(0,9000), labels = comma, breaks=scales::pretty_breaks()) +
  scale_color_brewer(palette = 'Dark2') 
ggsave("plots/sensor_width.png", height = 5, width = 7, dpi = 500)

####################################################################

# Constants
g <- 32.174                 # ft/s²
bank_angle <- 45           # degrees
ft_to_m <- 0.3048          # ft → m conversion

# Mach speeds
mach_speeds <- seq(0.4, 0.9, by = 0.1)

# FOV values
fovs <- c(15, 30, 60)

# Altitudes
altitudes <- seq(0, 25000, length.out = 200)  # ft

# Function: speed of sound vs altitude (ISA)
speed_of_sound_at_alt <- function(alt_ft) {
  alt_m <- alt_ft * ft_to_m
  T0 <- 288.15            # sea level temp in K
  lapse_rate <- -0.0065   # K/m
  T <- T0 + lapse_rate * alt_m
  gamma <- 1.4
  R <- 287.05
  sqrt(gamma * R * T) / ft_to_m  # returns ft/s
}

# Create grid for all combinations
grid <- expand.grid(FOV = fovs, Mach = mach_speeds, Altitude = altitudes)

# Footprint calculation in meters
grid$Footprint <- 2 * grid$Altitude * tan(grid$FOV / 2 * pi / 180) * ft_to_m

# Turn radius calculation with altitude-dependent speed of sound
grid <- grid %>%
  rowwise() %>%
  mutate(speed_of_sound = speed_of_sound_at_alt(Altitude),
         TurnDiameter = ((Mach * speed_of_sound)^2 / (g * tan(bank_angle * pi/180))) * ft_to_m * 2) %>%
  ungroup()

grid$FOV = factor(grid$FOV, levels= c('60','30','15'))

# Plot
ggplot(grid, aes(y = Altitude )) +
  # Footprint ribbon
  geom_ribbon(aes(xmin = 0, xmax = Footprint, fill=FOV, group = interaction(FOV, Mach))) +
  
  # Footprint lines
  geom_line(aes(x = Footprint, group = interaction(FOV, Mach)), size = 0.8) +
  
  # Turn radius dashed lines
  geom_line(aes(x = TurnDiameter, y = Altitude , color = as.factor(Mach), group = Mach),
            linetype = "dashed", size = 1.3) +
  
  # scale_x_continuous(breaks = scales::pretty_breaks(), limits=c(0,10000),  labels=comma) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(label=comma) +
  # Facets by FOV
  # facet_wrap(~FOV, #scales = "free_x",
  #            labeller = labeller(FOV = function(x) paste0("FOV = ", x, "°"))) +
  # 
  # Colors
  scale_fill_brewer(palette = 'Pastel1', name = "Swath Width", direction = -1) +
  # scale_fill_manual(values = fov_cols, name = "Sensor FOV") +
  scale_color_brewer(palette = 'Dark2', name = "Mach Number") +
  
  # Labels
  labs(
    x = "Distance (m)",
    y = "Altitude (ft)",
    title = "At Tactical Speeds, Turn Diameter > Sensor Footprint",
    subtitle = "Dashed lines = turn circle diameter (level turn, constant mach, 45° AOB)"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical") +
  guides(color = guide_legend(nrow = 1))+
  theme(panel.ontop = TRUE,
        panel.grid.major = element_line(color = "darkgrey", size  = 0.25),
        # panel.grid.minor = element_line(color = "darkgrey", size  = 0.02),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 30, 10, 10))  

ggsave("plots/turn_radius_sensor_width.png", height = 6, width = 7, dpi = 500)



