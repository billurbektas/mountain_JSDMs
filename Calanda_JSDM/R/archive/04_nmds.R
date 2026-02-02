# Load required packages
library(tidyverse)  # For data manipulation
library(vegan)      # For NMDS analysis
library(ggplot2)    # For plotting (part of tidyverse)

# Perform NMDS
nmds_result = metaMDS(Y, 
                      distance = "jaccard",  # Appropriate for presence/absence data
                      k = 2,                 # Number of dimensions
                      trymax = 100,          # Maximum number of random starts
                      autotransform = FALSE) # No transformation needed for binary data

# Extract NMDS site scores as a tibble for ggplot
nmds_sites = as_tibble(scores(nmds_result)$sites) %>%
  mutate(Site = rownames(Y))

# Calculate distance from origin for each species to identify influential ones
species_scores = species_scores %>%
  as_tibble() %>%
  mutate(
    distance_from_origin = sqrt(NMDS1^2 + NMDS2^2),
    influential = distance_from_origin > quantile(distance_from_origin, 0.7)  # Show top 20% influential
  )

# Improved plot with species names only
nmds_plot_with_species = ggplot() +
  # Add a better theme
  theme_bw() +
  
  # Plot sites (points only, no labels)
  geom_point(data = nmds_sites, aes(x = NMDS1, y = NMDS2), 
             size = 1, alpha = 0.7, color = "grey50") +
  
  # Add species - no points, only labels for influential species
  ggrepel::geom_text_repel(
    data = species_scores %>% filter(influential == TRUE),
    aes(x = NMDS1, y = NMDS2, label = Species),
    color = "darkred", 
    fontface = "italic",
    size = 3.5,
    max.overlaps = 15,
    box.padding = 0.5,
    segment.color = "darkred",
    segment.alpha = 0.6
  ) +
  
  # Cleaner axes and title
  labs(
    title = paste0("NMDS of Species Presence/Absence (Stress = ", 
                   round(nmds_result$stress, 3), ")"),
    x = "NMDS1", 
    y = "NMDS2",
    caption = "Based on Jaccard distance"
  ) +
  
  # Add a theme with larger text and better spacing
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  )

# Display the improved plot
print(nmds_plot_with_species)

# Display the combined plot
pdf(file = "nmds.pdf", height = 10, width = 8)
print(nmds_plot_with_species)
dev.off()
