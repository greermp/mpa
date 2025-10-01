library(dplyr)
library(ggplot2)
library(scales)

setwd("~/Anduril/interview/output")

theme_set(
  theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(face = "bold", size=12),
      strip.text = element_text(size=14, face="bold")
    )
)

# ----- Params -----
target_m    <- 15            # meters (corvette beam)
altitudes_ft<- seq(0, 25000, by = 1000)
px_required <- 4

# Sensor definitions
sensors <- data.frame(
  FOV_deg = c(15, 30, 60),
  Npx     = c(640, 1024, 1920),
  Sensor  = c("15° / 640px", "30° / 1024px", "60° / 1920px")
)

# ----- Helpers -----
ifov_from_fovN <- function(fov_deg, n_px) {
  (pi/180) * fov_deg / n_px   # radians per pixel
}

# Nadir-only pixels across target
pixels_across_nadir <- function(fov_deg, n_px, target_m, alt_ft) {
  R_m <- alt_ft * 0.3048
  if (R_m <= 0) return(Inf)
  ifov <- ifov_from_fovN(fov_deg, n_px)
  target_m / (ifov * R_m)
}

# Edge/worst-case pixels across target
pixels_across_edge <- function(fov_deg, n_px, target_m, alt_ft) {
  R_m <- alt_ft * 0.3048
  if (R_m <= 0) return(Inf)
  ifov <- ifov_from_fovN(fov_deg, n_px)
  
  # Half-FOV angle
  theta <- (fov_deg/2) * pi/180
  slant_range <- R_m / cos(theta)   # worst-case geometry at edge
  target_m / (ifov * slant_range)
}

# ----- Build dataset -----
grid <- do.call(rbind,
                lapply(1:nrow(sensors), function(i) {
                  sens <- sensors[i,]
                  data.frame(
                    alt_ft   = altitudes_ft,
                    Sensor   = sens$Sensor,
                    FOV_deg  = sens$FOV_deg,
                    Npx      = sens$Npx,
                    px_nadir = sapply(altitudes_ft, function(alt)
                      pixels_across_nadir(sens$FOV_deg, sens$Npx, target_m, alt)),
                    px_edge  = sapply(altitudes_ft, function(alt)
                      pixels_across_edge(sens$FOV_deg, sens$Npx, target_m, alt))
                  )
                })
)

# Reshape long for plotting
grid_long <- grid %>%
  tidyr::pivot_longer(cols = c(px_nadir, px_edge),
                      names_to = "Assumption", values_to = "px_across") %>%
  mutate(Assumption = recode(Assumption,
                             px_nadir = "Boresight",
                             px_edge  = "Worst-case (edge)"))

# ----- Plot -----
p <- ggplot(grid_long %>% filter(px_across<=20), aes(x = px_across, y = alt_ft,
                                                     color = Sensor, linetype = Assumption)) +
  geom_line(linewidth = 1.1, alpha=0.9) +
  geom_vline(xintercept = px_required, linetype = "dashed", size = 1.3) +
  scale_y_continuous(labels = scales::comma) +
  # scale_x_continuous(limits=c(0,20)) +
  labs(
    title = "Pixels Across vs Altitude",
    subtitle = sprintf("L ≈ %g m, boresight vs Edge GSD", target_m),
    x = "Pixels across target (px)",
    y = "Altitude (ft)",
    linetype = "𝛼"
  ) +
  scale_color_brewer(palette = "Dark2")
p

ggsave("plots/pix_alt_nadir_vs_edge.png", p, height = 5.25, width = 7, dpi = 500)
