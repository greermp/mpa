library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(gghighlight)

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

AOI_area = 10000

df <- read_csv("mpa.csv") %>%
  mutate(area_per_cost = area_sanitized_km2 / cost,
         per_mission_complete = area_sanitized_km2/AOI_area) %>% 
  mutate(speed_cat = 
        case_when(
         speed_mach <= 0.65  ~"slow",
         speed_mach >0.65    ~ "fast"))

df_big <- read_csv("mpa_big.csv") %>%
  mutate(area_per_cost = area_sanitized_km2 / cost,
         per_mission_complete = area_sanitized_km2/AOI_area) %>% 
  mutate(speed_cat = 
           case_when(
             speed_mach <= 0.65  ~"slow",
             speed_mach >0.65    ~ "fast"))
  
slow <- df %>% filter(speed_cat == "slow")

# Cost per mission complete
ggplot(df, aes(x = per_mission_complete, y =  cost, color = factor(alt_ft))) +
  geom_point(aes(shape = factor(speed_mach)), 
             size=4,
             show.legend = FALSE
  ) +
  facet_wrap(~sensor_name) +
  scale_y_continuous(limits=c(0,50), expand = expansion(mult=0))+
  scale_x_continuous(labels = percent, limits=c(0,1), expand = expansion(mult=0.05)) +
  # scale_y_reverse() +
  scale_shape_manual(guide = guide_legend(nrow = 1), 
                     values=c(16, 17, 15, 3, 4, 1)
  ) +
  labs(x="AOI Sanitized (%)",
       y="Procurement Cost ($M)",
       linetype = "",
       shape="Altitude (ft))", 
       color="Altitude (ft)") +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'bottom') +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5))



ggsave("plots/fast_stinks_cost.png", height = 3.65, width = 12.5, dpi = 500)

ggplot(slow, aes(x = per_mission_complete, y =  cost, color = factor(alt_ft))) +
  geom_point(aes(shape = factor(speed_mach)), 
             size=4,
             show.legend = FALSE
  ) +
  facet_wrap(~sensor_name) +
  scale_y_continuous(limits=c(0,50), expand = expansion(mult=0))+
  scale_x_continuous(labels = percent, limits=c(0,1), expand = expansion(mult=0.05)) +
  scale_shape_manual(guide = guide_legend(nrow = 1), 
                     values=c(16, 17, 15, 3, 4, 1)
  ) +
  labs(x="AOI Sanitized (%)",
       y="Procurement Cost ($M)",
       linetype = "",
       shape="Altitude (ft))", 
       color="Altitude (ft)") +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'bottom') +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5))


ggsave("plots/fast_stinks_filt.png", height = 3.65, width = 12.5, dpi = 500)


### Mult sorties

many <- slow %>%
  mutate(sorties_req = 
    case_when(
      area_sanitized_km2>=AOI_area ~ 1,
      TRUE ~ ceiling(AOI_area/area_sanitized_km2)
    ),
    cost_n=cost*sorties_req,
    mission_time_mult = 1/per_mission_complete*mission_time_hrs/sorties_req)

ggplot(many, aes(x = sorties_req, y =  cost_n, color = factor(alt_ft))) +
  geom_point(aes(shape = factor(speed_mach)), 
             size=4,
             show.legend = FALSE
  ) +
  facet_wrap(~sensor_name) +
  scale_y_continuous(limits=c(0,50), expand = expansion(mult=0))+
  scale_x_continuous(breaks = seq(1,8,1)) +
  scale_shape_manual(guide = guide_legend(nrow = 1), 
                     values=c(16, 17, 15, 3, 4, 1)
  ) +
  labs(x="Sorties Required",
       y="Procurement Cost * Sorties Required",
       linetype = "",
       shape="Altitude (ft))", 
       color="Altitude (ft)") +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'bottom') +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5))

ggsave("plots/mult_sorties.png", height = 3.65, width = 12.5, dpi = 500)


ggplot(many, aes(x = sorties_req, y =  cost_n, color = factor(alt_ft))) +
  geom_point(aes(shape = factor(speed_mach)), 
             size=4,
             show.legend = FALSE
  ) +
  gghighlight((sorties_req>1 & sensor_name=='EO/IR Sensor 2'), calculate_per_facet = TRUE) +
  facet_wrap(~sensor_name) +
  
  scale_y_continuous(limits=c(0,50), expand = expansion(mult=0))+
  scale_x_continuous(breaks = seq(1,8,1)) +
  scale_shape_manual(guide = guide_legend(nrow = 1), 
                     values=c(16, 17, 15, 3, 4, 1)
  ) +
  labs(x="Sorties Required",
       y="Procurement Cost * Sorties Required",
       linetype = "",
       shape="Altitude (ft))", 
       color="Altitude (ft)") +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'bottom') +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5))

ggsave("plots/mult_sorties_highight.png", height = 3.65, width = 12.5, dpi = 500)

many_filt <- many %>% filter(sensor_name=='EO/IR Sensor 2', sorties_req>1)


# Speed
 ggplot(many_filt, aes(x = mission_time_mult, y =  cost_n, color = factor(alt_ft))) +
  geom_point(aes(shape = factor(speed_mach)), 
             size=4,
             show.legend = FALSE
  ) +  
  geom_text(data=many_filt %>% filter(speed_mach>=0.6), show.legend = FALSE,
    aes(label=paste0("n =", sorties_req), x=mission_time_mult+0.16, y=cost_n+1.2)) +
   scale_y_continuous(limits=c(0,50), expand = expansion(mult=0))+
   
  labs(x="Mission Completion Time (with n sorties)",
       y="Procurement Cost * Sorties Required",
       title="Similar mission completion time",
       linetype = "",
       shape="Altitude (ft))", 
       color="Altitude (ft)") +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'bottom') +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5))



ggsave("plots/mult_sorties_speed.png", height = 3.8, width = 4.2, dpi = 500)

ggplot(many_filt, aes(x = mission_time_mult, y =  cost_n, color = factor(alt_ft))) +
  geom_point(aes(shape = factor(speed_mach)), 
             size=4,
             show.legend = FALSE
  ) +  
  scale_y_continuous(limits=c(0,50), expand = expansion(mult=0))+
  gghighlight((speed_mach>=0.6 & alt_ft==10000) | speed_mach==0.4) +
  geom_text(data=many_filt %>% filter(speed_mach==0.4), show.legend = FALSE,
            aes(label=paste0("n =", sorties_req), x=mission_time_mult-0.16, y=cost_n+1.2)) +
  
  labs(x="Mission Completion Time (with n sorties)",
       y="Procurement Cost * Sorties Required",
       title="Mass Tradeoffs",
       linetype = "",
       shape="Altitude (ft))", 
       color="Altitude (ft)") +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'bottom') +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5))

ggsave("plots/options.png", height = 3.8, width = 4.2, dpi = 500)


ggplot(slow, aes(x = per_mission_complete, y =  cost, color = factor(alt_ft))) +
  geom_point(aes(shape = factor(speed_mach)), 
             size=4,
             show.legend = FALSE
  ) +
  gghighlight(per_mission_complete==1, calculate_per_facet = TRUE) +
  facet_wrap(~sensor_name) +
  scale_y_continuous(limits=c(0,50), expand = expansion(mult=0))+
  scale_x_continuous(labels = percent, limits=c(0,1), expand = expansion(mult=0.05)) +
  scale_shape_manual(guide = guide_legend(nrow = 1), 
                     values=c(16, 17, 15, 3, 4, 1)
  ) +
  labs(x="AOI Sanitized (%)",
       y="Procurement Cost ($M)",
       linetype = "",
       shape="Altitude (ft))", 
       color="Altitude (ft)") +
  scale_color_brewer(palette = 'Dark2') +
  theme(legend.position = 'bottom') +
  theme(panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5))

ggsave("plots/mission_complete.png", height = 3.65, width = 12.5, dpi = 500)


success_big <- df_big %>% 
  filter(speed_cat == "slow") %>% 
  filter(
    (sensor_name=="EO/IR Sensor 3" & alt_ft>5000) |
    (sensor_name=="EO/IR Sensor 2" & alt_ft>10000) | 
      (sensor_name=="EO/IR Sensor 1" & alt_ft==25000)  
  )

success_small <- slow %>% filter(per_mission_complete==1)
success_small$alt_ft = factor(success_small$alt_ft)
success_small$speed_mach = factor(success_small$speed_mach)

success_big$alt_ft = factor(success_big$alt_ft)
success_big$speed_mach = factor(success_big$speed_mach)

success_small <- success_small %>% 
  mutate(km_per_m = area_sanitized_km2/cost)

ggplot(success_small, aes(x = mission_time_hrs, y =  cost)) +
  geom_point(aes(shape = speed_mach,
                 color=sensor_name,
                 fill = alt_ft), stroke=1.3, 
             size=5,
  ) +
  scale_y_continuous(limits=c(0,50), expand = expansion(mult=0))+
  scale_shape_manual(values=c(21, 24, 22)) +
  scale_color_manual(values = c('#CCC', '#FFFFFF', 'black')) +
  scale_fill_brewer(palette = 'Dark2') +
  guides(color = guide_legend(override.aes = list(shape = 23, size = 5, stroke = 1.5, fill = "white")),
         fill = guide_legend(override.aes = list( shape = 23, size = 5, stroke = 1.5, color = NA)))








pareto_points <- success_small %>%
  arrange(mission_time_hrs) %>%
  mutate(min_cost = cummin(cost)) %>%
  filter(cost == min_cost)

ggplot(success_small, aes(x = mission_time_hrs, y = cost)) +
  geom_point(aes(shape = speed_mach,
                 color = sensor_name,
                 fill  = alt_ft),
             stroke = 1.3, size = 5) +
  geom_line(data = pareto_points, color = "cyan", linewidth = 1.2) +
  geom_point(data = pareto_points, color = "black", size = 1.5) +
  scale_y_continuous(limits = c(0,50), expand = expansion(mult = 0)) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = c('#CCC', '#FFCC00', 'black')) +
  scale_fill_brewer(palette = 'Dark2') +
  guides(color = guide_legend(override.aes = list(shape = 23, size = 5, stroke = 1.5, fill = "white")),
         fill  = guide_legend(override.aes = list(shape = 23, size = 5, stroke = 1.5, color = NA))) +
  labs(x="Mission duration (hrs)",
       y="Procurement Cost ($M)",
       # linetype = "",
       title="Cost/Mission Length Pareto Front",
       subtitle="100% AOI Sanitized",
       shape="Speed (Mach))", 
       color="Sensor",
       fill="Altitude (ft)")


ggsave("plots/pareto.png", height = 5.1, width = 7, dpi = 500)


success_small <- success_small %>%
  mutate(endurance_remaining = endurance_hr - mission_time_hrs) 

pareto_points_endurance <- success_small %>% 
  arrange(desc(endurance_remaining)) %>%
  mutate(min_cost = cummin(cost)) %>%
  filter(cost == min_cost)

ggplot(success_small, aes(x = endurance_remaining, y = cost)) +
  geom_point(aes(shape = speed_mach,
                 color = sensor_name,
                 fill  = alt_ft),
             stroke = 1.3, size = 5) +
  geom_line(data = pareto_points_endurance, color = "cyan", linewidth = 1.2) +
  geom_point(data = pareto_points_endurance, color = "black", size = 1.5) +
  scale_y_continuous(limits = c(0,50), expand = expansion(mult = 0)) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_manual(values = c('#CCC', '#FFCC00', 'black')) +
  labs(x = "Endurance Remaining (hrs)", y = "Cost ($M)") +
  guides(color = guide_legend(
    override.aes = list(shape = 23, size = 5, stroke = 1.5, fill = "white")),
    fill = guide_legend(
      override.aes = list(shape = 23, size = 5, stroke = 1.5, color = NA)))+
  labs(x="Residual Endurance (hrs)",
       y="Procurement Cost ($M)",
       # linetype = "",
       title="Cost/Residual Endurance Pareto Front",
       subtitle="100% AOI Sanitized",
       shape="Speed (Mach))", 
       color="Sensor",
       fill="Altitude (ft)")


ggsave("plots/pareto2.png", height = 5.1, width = 7, dpi = 500)

