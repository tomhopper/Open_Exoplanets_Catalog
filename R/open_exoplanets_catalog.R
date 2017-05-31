library(MASS)
library(fields)
library(rvest)
library(dplyr)
library(ggplot2)
library(magrittr)
library(scales)
library(ggrepel)
library(Cairo)
library(extrafont)

#font_import()
loadfonts()

file_oec <- "https://raw.githubusercontent.com/OpenExoplanetCatalogue/oec_tables/master/comma_separated/open_exoplanet_catalogue.txt"
#file_kep <- "https://raw.githubusercontent.com/OpenExoplanetCatalogue/oec_tables/master/comma_separated/open_exoplanet_catalogue_kepler.txt"

oec_df <- read.csv(file = file_oec, header = FALSE, skip = 31, stringsAsFactors = FALSE)
oec_names <- read.csv(file = file_oec, header = FALSE, skip = 5, nrows = 25, stringsAsFactors = FALSE)[[1]]
oec_names <- sapply(strsplit(oec_names, split = ": "), "[[", 2)
colnames(oec_df) <- make.names(oec_names)

# kep_df <- read.csv(file = file_kep, header = FALSE, skip = 31, stringsAsFactors = FALSE)
# kep_names <- read.csv(file = file_kep, header = FALSE, skip = 5, nrows = 25, stringsAsFactors = FALSE)[[1]]
# kep_names <- sapply(strsplit(kep_names, split = ": "), "[[", 2)
# colnames(kep_df) <- make.names(kep_names)

oec_df$density <- oec_df$Planetary.mass..Jupiter.masses. / (4*pi*(oec_df$Radius..Jupiter.radii.^3)/3)
oec_df$density <- oec_df$density / oec_df$density[oec_df$Primary.identifier.of.planet == "Jupiter"]

oec_df$plot_density <- rep(NA, times = nrow(oec_df))
index <- which(!is.na(oec_df$Planetary.mass..Jupiter.masses.) & !is.na(oec_df$Radius..Jupiter.radii.))
oec_df$plot_density[index] <- fields::interp.surface(
  MASS::kde2d(log10(oec_df[index,]$Planetary.mass..Jupiter.masses.), log10(oec_df[index,]$Radius..Jupiter.radii.)), 
  log10(oec_df[index,c("Planetary.mass..Jupiter.masses.", "Radius..Jupiter.radii.")]))

lo <- sort(oec_df$density, decreasing = FALSE)[1]
hi <- sort(oec_df$density, decreasing = TRUE)[1]
if(lo < 1) lo_start <- ceiling(log10(lo)) else lo_start <- floor(log10(lo))
hi_end <- floor(log10(hi))
density_values <- rescale(quantile(x = sort(oec_df$density[!is.na(oec_df$density)])))
density_colors_temp <- RColorBrewer::brewer.pal(n = hi_end - lo_start, name = "PuBuGn")
density_colors_temp[1] <- muted(density_colors_temp[1])
density_colors <- colorRampPalette(c(density_colors_temp[1], density_colors_temp[length(density_colors_temp)]))(length(density_values))

sol_colors <- RColorBrewer::brewer.pal(n = 1, name = "Oranges")

sol_planets_df <- oec_df %>% 
  filter(Primary.identifier.of.planet %in% c("Mercury","Venus","Earth","Mars","Jupiter","Saturn","Neptune","Uranus"))

planets_p <- oec_df %>% 
  select(Primary.identifier.of.planet, Planetary.mass..Jupiter.masses., Radius..Jupiter.radii., density, plot_density) %>% 
  na.omit() %>% 
  ggplot() +
  aes(x = Planetary.mass..Jupiter.masses., y = Radius..Jupiter.radii.) + #alpha = 1/plot_density,
  stat_function(aes(), fun = function(x) {log10(x^(1/3))}, color = density_colors[3], alpha = 0.25) +
  geom_point(aes(color = density, alpha = rescale(1/plot_density, to = c(0.15, 1))), size = 2.5) +
  #geom_point(aes(color = density), shape = 1, alpha = 0.7, size = 2.5) +
  geom_point(data = sol_planets_df, color = sol_colors[2], size = 2) +
  geom_text_repel(data = sol_planets_df, aes(label = Primary.identifier.of.planet), alpha = 0.75, size = 3.5) +
  scale_x_log10(labels = scales::comma, name = expression(paste("Mass,"  %*% M[Jup]))) +
  scale_y_log10(labels = scales::comma, name = expression(paste("Radius,"  %*% R[Jup]))) +
  scale_alpha(guide = FALSE) +
  scale_color_gradientn(guide = FALSE, labels = scales::comma, colours = density_colors, values = density_values, name = expression(paste("Density,"  %*% rho[jup]))) +
  coord_equal() +
  labs(title = "Size and mass of discovered exoplanets",
       subtitle = paste(nrow(oec_df %>% 
                               select(Planetary.mass..Jupiter.masses., Radius..Jupiter.radii., density, plot_density) %>% 
                               na.omit()), 
                        "planets plotted of",
                        nrow(oec_df),
                        "discovered"),
       caption = "plot by Thomas Hopper 2017, data from Open Exoplanets Catalogue") +
  annotation_logticks(size = 0.1, color = "grey") +
  annotate(geom = "text", x = 0.1, y = 5, label = "less dense than Jupiter", size = 3, colour = muted(density_colors[1])) +
  annotate(geom = "text", x = 2, y = 0.1, label = "more dense than Jupiter", size = 3, colour = density_colors[3]) +
  theme_minimal() +
  theme(plot.caption = element_text(size = 6))
planets_p

png(filename = "figs/planets.png", width = 2.2*700, height = 700, units = "px", type = "cairo-png", res = 150)
print(planets_p)
dev.off()
