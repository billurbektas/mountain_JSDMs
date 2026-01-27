# ==============================================================================
# TRAIT COVERAGE ASSESSMENT
# ==============================================================================
# This script assesses what percentage of species in each community have
# trait data available, both by presence/absence and weighted by relative cover
cat("\n=== Trait Coverage Assessment ===\n")

library(tidyverse)
library(dgof)  # For Cramér-von Mises test

# Load the species trait summary
trait_summary = read_csv("output/species_trait_summary.csv")%>%
  mutate(plant_species = case_when(plant_species == "Acinos alpinus" ~ "Clinopodium alpinum",
                                   plant_species == "Euphrasia rostkoviana" ~ "Euphrasia officinalis subsp. pratensis",
                                   plant_species == "Festuca pratensis" ~ "Lolium pratense",
                                   plant_species == "Helictotrichon pubescens" ~ "Avenula pubescens",
                                   plant_species == "Hieracium lactucella" ~ "Pilosella lactucella",
                                   plant_species == "Hieracium pilosella" ~ "Pilosella officinarum",
                                   plant_species == "Polygala chamaebuxus" ~ "Polygaloides chamaebuxus",
                                   plant_species == "Polygonum viviparum" ~ "Bistorta vivipara",
                                   plant_species == "Viola riviana" ~ "Viola riviniana",
                                   .default = plant_species
                                   ))

# Get list of species with trait data (at least one trait measured)
species_with_traits = trait_summary %>%
  filter(n_samples > 0) %>%
  pull(plant_species)

cat(sprintf("\nTotal species with trait data: %d\n", length(species_with_traits)))
cat(sprintf("Total species in trait summary: %d\n", nrow(trait_summary)))

# Load abundance data (veg.abund should be in the workspace or loaded from RData)
if(!exists("veg.abund")) {
  cat("\nLoading data from starter_data_25.04.25.RData...\n")
  load("output/starter_data_25.04.25.RData")
}

# Check if veg.abund exists
if(!exists("veg.abund")) {
  stop("Error: veg.abund not found in workspace or RData file")
}

cat(sprintf("\nNumber of communities (plots): %d\n", nrow(veg.abund)))
cat(sprintf("Number of species in communities: %d\n", ncol(veg.abund)))

# Convert veg.abund to long format for easier analysis
abund_long = veg.abund %>%
  as.data.frame() %>%
  rownames_to_column("plot_id") %>%
  pivot_longer(-plot_id, names_to = "species", values_to = "rel_cover") %>%
  filter(rel_cover > 0)  # Only keep species present in each plot

# Calculate coverage metrics for each community
community_coverage = abund_long %>%
  group_by(plot_id) %>%
  summarize(
    # Total species richness
    n_species = n(),

    # Species with trait data
    n_species_with_traits = sum(species %in% species_with_traits),

    # Percentage of species with traits (presence/absence)
    pct_species_with_traits = 100 * n_species_with_traits / n_species,

    # Total relative cover
    total_cover = sum(rel_cover),

    # Relative cover of species with traits
    cover_with_traits = sum(rel_cover[species %in% species_with_traits]),

    # Percentage of cover represented by species with traits
    pct_cover_with_traits = 100 * cover_with_traits / total_cover,

    .groups = "drop"
  )

# Summary statistics
cat("\n=== Coverage Summary Statistics ===\n")
cat("\n--- Species Richness-based Coverage (Presence/Absence) ---\n")
cat(sprintf("Mean: %.1f%%\n", mean(community_coverage$pct_species_with_traits)))
cat(sprintf("Median: %.1f%%\n", median(community_coverage$pct_species_with_traits)))
cat(sprintf("Range: %.1f%% - %.1f%%\n",
            min(community_coverage$pct_species_with_traits),
            max(community_coverage$pct_species_with_traits)))
cat(sprintf("SD: %.1f%%\n", sd(community_coverage$pct_species_with_traits)))

cat("\n--- Relative Cover-based Coverage ---\n")
cat(sprintf("Mean: %.1f%%\n", mean(community_coverage$pct_cover_with_traits)))
cat(sprintf("Median: %.1f%%\n", median(community_coverage$pct_cover_with_traits)))
cat(sprintf("Range: %.1f%% - %.1f%%\n",
            min(community_coverage$pct_cover_with_traits),
            max(community_coverage$pct_cover_with_traits)))
cat(sprintf("SD: %.1f%%\n", sd(community_coverage$pct_cover_with_traits)))

# Identify communities with low coverage
low_coverage_pa = community_coverage %>%
  filter(pct_species_with_traits < 80) %>%
  arrange(pct_species_with_traits)

low_coverage_cover = community_coverage %>%
  filter(pct_cover_with_traits < 80) %>%
  arrange(pct_cover_with_traits)

cat(sprintf("\n--- Communities with <80%% species coverage: %d ---\n", nrow(low_coverage_pa)))
if(nrow(low_coverage_pa) > 0) {
  print(low_coverage_pa %>% select(plot_id, n_species, n_species_with_traits, pct_species_with_traits))
}

cat(sprintf("\n--- Communities with <80%% cover coverage: %d ---\n", nrow(low_coverage_cover)))
if(nrow(low_coverage_cover) > 0) {
  print(low_coverage_cover %>% select(plot_id, total_cover, cover_with_traits, pct_cover_with_traits))
}

# Identify species that are common but lack trait data
species_missing_traits = abund_long %>%
  filter(!species %in% species_with_traits) %>%
  group_by(species) %>%
  summarize(
    n_occurrences = n(),
    total_cover = sum(rel_cover),
    mean_cover = mean(rel_cover),
    .groups = "drop"
  ) %>%
  arrange(desc(n_occurrences))

cat(sprintf("\n--- Species lacking trait data: %d ---\n", nrow(species_missing_traits)))
cat("Top 20 most common species without trait data:\n")
print(species_missing_traits %>% head(20))

# Save results
write_csv(community_coverage, "output/community_trait_coverage.csv")
write_csv(species_missing_traits, "output/species_missing_traits.csv")

cat("\n=== Output files created ===\n")
cat("- output/community_trait_coverage.csv: Coverage metrics for each community\n")
cat("- output/species_missing_traits.csv: Species lacking trait data with their prevalence\n")

# Create visualization
library(ggplot2)
library(patchwork)

# Histogram of coverage by presence/absence
p1 = ggplot(community_coverage, aes(x = pct_species_with_traits)) +
  geom_histogram(bins = 30, fill = "#4ECDC4", color = "black", alpha = 0.7) +
  geom_vline(xintercept = median(community_coverage$pct_species_with_traits),
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "% of species with trait data",
       y = "Number of communities",
       title = "Trait coverage by species richness") +
  theme_bw() +
  theme(text = element_text(size = 12))

# Histogram of coverage by relative cover
p2 = ggplot(community_coverage, aes(x = pct_cover_with_traits)) +
  geom_histogram(bins = 30, fill = "#FF6B6B", color = "black", alpha = 0.7) +
  geom_vline(xintercept = median(community_coverage$pct_cover_with_traits),
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(x = "% of relative cover with trait data",
       y = "Number of communities",
       title = "Trait coverage by relative cover") +
  theme_bw() +
  theme(text = element_text(size = 12))

# Scatter plot comparing the two metrics
p3 = ggplot(community_coverage, aes(x = pct_species_with_traits,
                                     y = pct_cover_with_traits)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "% species with trait data",
       y = "% cover with trait data",
       title = "Comparison of coverage metrics") +
  theme_bw() +
  theme(text = element_text(size = 12))

# Species richness vs coverage
p4 = ggplot(community_coverage, aes(x = n_species, y = pct_cover_with_traits)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "#4ECDC4") +
  labs(x = "Species richness",
       y = "% cover with trait data",
       title = "Coverage vs. species richness") +
  theme_bw() +
  theme(text = element_text(size = 12))

# Combine plots
p_combined = (p1 | p2) / (p3 | p4)

pdf("plot/trait_coverage_assessment.pdf", width = 12, height = 10)
print(p_combined)
dev.off()

cat("\n- plot/trait_coverage_assessment.pdf: Visualization of trait coverage\n")

# ==============================================================================
# BIAS ASSESSMENT: Check if high-coverage communities differ in variance components
# ==============================================================================
cat("\n=== Checking for bias in variance components ===\n")

# Extract variance components from JSDM results
if(!exists("res")) {
  cat("\nLoading JSDM results...\n")
  load("output/starter_data_25.04.25.RData")
}

variance_data = res$internals$Sites %>%
  as.data.frame() %>%
  rownames_to_column("plot_id") %>%
  select(plot_id, env, codist, spa) %>%
  mutate(
    total = env + codist + spa,
    env_prop = env / total,
    codist_prop = codist / total,
    spa_prop = spa / total
  )

# Merge with coverage data
combined = community_coverage %>%
  left_join(variance_data, by = "plot_id") %>%
  mutate(
    coverage_group = ifelse(pct_cover_with_traits >= 75, "≥75% coverage", "All data")
  )

# Create a dataset that includes all data as one group and high-coverage as another
all_data = combined %>%
  mutate(coverage_group = "All data")

high_coverage = combined %>%
  filter(pct_cover_with_traits >= 75) %>%
  mutate(coverage_group = "≥75% coverage")

comparison_data = bind_rows(all_data, high_coverage)

cat(sprintf("\nTotal communities: %d\n", nrow(combined)))
cat(sprintf("Communities with ≥75%% coverage: %d (%.1f%%)\n",
            sum(combined$pct_cover_with_traits >= 75),
            100 * mean(combined$pct_cover_with_traits >= 75)))

# Statistical tests comparing all data vs ≥75% coverage
cat("\n=== Statistical tests (All data vs ≥75% coverage) ===\n")

# T-tests (for means)
test_env = t.test(all_data$env, high_coverage$env)
test_codist = t.test(all_data$codist, high_coverage$codist)
test_spa = t.test(all_data$spa, high_coverage$spa)
test_env_prop = t.test(all_data$env_prop, high_coverage$env_prop)
test_codist_prop = t.test(all_data$codist_prop, high_coverage$codist_prop)
test_spa_prop = t.test(all_data$spa_prop, high_coverage$spa_prop)

# Wilcoxon tests (for medians)
wilcox_env = wilcox.test(all_data$env, high_coverage$env)
wilcox_codist = wilcox.test(all_data$codist, high_coverage$codist)
wilcox_spa = wilcox.test(all_data$spa, high_coverage$spa)
wilcox_env_prop = wilcox.test(all_data$env_prop, high_coverage$env_prop)
wilcox_codist_prop = wilcox.test(all_data$codist_prop, high_coverage$codist_prop)
wilcox_spa_prop = wilcox.test(all_data$spa_prop, high_coverage$spa_prop)

# Cramér-von Mises tests (for entire distributions)
cvm_env = cvm.test(high_coverage$env, all_data$env)
cvm_codist = cvm.test(high_coverage$codist, all_data$codist)
cvm_spa = cvm.test(high_coverage$spa, all_data$spa)
cvm_env_prop = cvm.test(high_coverage$env_prop, all_data$env_prop)
cvm_codist_prop = cvm.test(high_coverage$codist_prop, all_data$codist_prop)
cvm_spa_prop = cvm.test(high_coverage$spa_prop, all_data$spa_prop)

cat("\n--- Raw variance (t-tests for means) ---\n")
cat(sprintf("Environment: t=%.3f, p=%.4f %s\n", test_env$statistic, test_env$p.value,
            ifelse(test_env$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution: t=%.3f, p=%.4f %s\n", test_codist$statistic, test_codist$p.value,
            ifelse(test_codist$p.value < 0.05, "***", "")))
cat(sprintf("Spatial: t=%.3f, p=%.4f %s\n", test_spa$statistic, test_spa$p.value,
            ifelse(test_spa$p.value < 0.05, "***", "")))

cat("\n--- Raw variance (Wilcoxon tests for medians) ---\n")
cat(sprintf("Environment: W=%.1f, p=%.4f %s\n", wilcox_env$statistic, wilcox_env$p.value,
            ifelse(wilcox_env$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution: W=%.1f, p=%.4f %s\n", wilcox_codist$statistic, wilcox_codist$p.value,
            ifelse(wilcox_codist$p.value < 0.05, "***", "")))
cat(sprintf("Spatial: W=%.1f, p=%.4f %s\n", wilcox_spa$statistic, wilcox_spa$p.value,
            ifelse(wilcox_spa$p.value < 0.05, "***", "")))

cat("\n--- Proportional variance (t-tests for means) ---\n")
cat(sprintf("Environment prop: t=%.3f, p=%.4f %s\n", test_env_prop$statistic, test_env_prop$p.value,
            ifelse(test_env_prop$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution prop: t=%.3f, p=%.4f %s\n", test_codist_prop$statistic, test_codist_prop$p.value,
            ifelse(test_codist_prop$p.value < 0.05, "***", "")))
cat(sprintf("Spatial prop: t=%.3f, p=%.4f %s\n", test_spa_prop$statistic, test_spa_prop$p.value,
            ifelse(test_spa_prop$p.value < 0.05, "***", "")))

cat("\n--- Proportional variance (Wilcoxon tests for medians) ---\n")
cat(sprintf("Environment prop: W=%.1f, p=%.4f %s\n", wilcox_env_prop$statistic, wilcox_env_prop$p.value,
            ifelse(wilcox_env_prop$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution prop: W=%.1f, p=%.4f %s\n", wilcox_codist_prop$statistic, wilcox_codist_prop$p.value,
            ifelse(wilcox_codist_prop$p.value < 0.05, "***", "")))
cat(sprintf("Spatial prop: W=%.1f, p=%.4f %s\n", wilcox_spa_prop$statistic, wilcox_spa_prop$p.value,
            ifelse(wilcox_spa_prop$p.value < 0.05, "***", "")))

cat("\n--- Raw variance (Cramér-von Mises tests for distributions) ---\n")
cat(sprintf("Environment: CvM=%.4f, p=%.4f %s\n", cvm_env$statistic, cvm_env$p.value,
            ifelse(cvm_env$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution: CvM=%.4f, p=%.4f %s\n", cvm_codist$statistic, cvm_codist$p.value,
            ifelse(cvm_codist$p.value < 0.05, "***", "")))
cat(sprintf("Spatial: CvM=%.4f, p=%.4f %s\n", cvm_spa$statistic, cvm_spa$p.value,
            ifelse(cvm_spa$p.value < 0.05, "***", "")))

cat("\n--- Proportional variance (Cramér-von Mises tests for distributions) ---\n")
cat(sprintf("Environment prop: CvM=%.4f, p=%.4f %s\n", cvm_env_prop$statistic, cvm_env_prop$p.value,
            ifelse(cvm_env_prop$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution prop: CvM=%.4f, p=%.4f %s\n", cvm_codist_prop$statistic, cvm_codist_prop$p.value,
            ifelse(cvm_codist_prop$p.value < 0.05, "***", "")))
cat(sprintf("Spatial prop: CvM=%.4f, p=%.4f %s\n", cvm_spa_prop$statistic, cvm_spa_prop$p.value,
            ifelse(cvm_spa_prop$p.value < 0.05, "***", "")))

# Summary statistics
cat("\n=== Summary statistics by group ===\n")
summary_stats = comparison_data %>%
  group_by(coverage_group) %>%
  summarize(
    n = n(),
    env_mean = mean(env, na.rm = TRUE), env_median = median(env, na.rm = TRUE), env_sd = sd(env, na.rm = TRUE),
    codist_mean = mean(codist, na.rm = TRUE), codist_median = median(codist, na.rm = TRUE), codist_sd = sd(codist, na.rm = TRUE),
    spa_mean = mean(spa, na.rm = TRUE), spa_median = median(spa, na.rm = TRUE), spa_sd = sd(spa, na.rm = TRUE),
    env_prop_mean = mean(env_prop, na.rm = TRUE), env_prop_median = median(env_prop, na.rm = TRUE), env_prop_sd = sd(env_prop, na.rm = TRUE),
    codist_prop_mean = mean(codist_prop, na.rm = TRUE), codist_prop_median = median(codist_prop, na.rm = TRUE), codist_prop_sd = sd(codist_prop, na.rm = TRUE),
    spa_prop_mean = mean(spa_prop, na.rm = TRUE), spa_prop_median = median(spa_prop, na.rm = TRUE), spa_prop_sd = sd(spa_prop, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# Calculate mean and median differences (with NA removal)
mean_diff_env = mean(all_data$env, na.rm = TRUE) - mean(high_coverage$env, na.rm = TRUE)
mean_diff_codist = mean(all_data$codist, na.rm = TRUE) - mean(high_coverage$codist, na.rm = TRUE)
mean_diff_spa = mean(all_data$spa, na.rm = TRUE) - mean(high_coverage$spa, na.rm = TRUE)
mean_diff_env_prop = mean(all_data$env_prop, na.rm = TRUE) - mean(high_coverage$env_prop, na.rm = TRUE)
mean_diff_codist_prop = mean(all_data$codist_prop, na.rm = TRUE) - mean(high_coverage$codist_prop, na.rm = TRUE)
mean_diff_spa_prop = mean(all_data$spa_prop, na.rm = TRUE) - mean(high_coverage$spa_prop, na.rm = TRUE)

median_diff_env = median(all_data$env, na.rm = TRUE) - median(high_coverage$env, na.rm = TRUE)
median_diff_codist = median(all_data$codist, na.rm = TRUE) - median(high_coverage$codist, na.rm = TRUE)
median_diff_spa = median(all_data$spa, na.rm = TRUE) - median(high_coverage$spa, na.rm = TRUE)
median_diff_env_prop = median(all_data$env_prop, na.rm = TRUE) - median(high_coverage$env_prop, na.rm = TRUE)
median_diff_codist_prop = median(all_data$codist_prop, na.rm = TRUE) - median(high_coverage$codist_prop, na.rm = TRUE)
median_diff_spa_prop = median(all_data$spa_prop, na.rm = TRUE) - median(high_coverage$spa_prop, na.rm = TRUE)

# Store t-test, Wilcoxon, and Cramér-von Mises results for plotting
test_results = data.frame(
  component = c("Environment", "Species associations", "Spatial",
                "Environment", "Species associations", "Spatial"),
  variance_type = c(rep("Raw variance", 3), rep("Proportional variance", 3)),
  p_value_t = c(test_env$p.value, test_codist$p.value, test_spa$p.value,
                test_env_prop$p.value, test_codist_prop$p.value, test_spa_prop$p.value),
  p_value_w = c(wilcox_env$p.value, wilcox_codist$p.value, wilcox_spa$p.value,
                wilcox_env_prop$p.value, wilcox_codist_prop$p.value, wilcox_spa_prop$p.value),
  p_value_cvm = c(cvm_env$p.value, cvm_codist$p.value, cvm_spa$p.value,
                  cvm_env_prop$p.value, cvm_codist_prop$p.value, cvm_spa_prop$p.value),
  mean_diff = c(mean_diff_env, mean_diff_codist, mean_diff_spa,
                mean_diff_env_prop, mean_diff_codist_prop, mean_diff_spa_prop),
  median_diff = c(median_diff_env, median_diff_codist, median_diff_spa,
                  median_diff_env_prop, median_diff_codist_prop, median_diff_spa_prop)
) %>%
  mutate(
    sig_label_t = case_when(
      p_value_t < 0.001 ~ "***",
      p_value_t < 0.01 ~ "**",
      p_value_t < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_label_w = case_when(
      p_value_w < 0.001 ~ "***",
      p_value_w < 0.01 ~ "**",
      p_value_w < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_label_cvm = case_when(
      p_value_cvm < 0.001 ~ "***",
      p_value_cvm < 0.01 ~ "**",
      p_value_cvm < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = sprintf("Mean Δ=%.3f (t: p=%.4f %s)\nMedian Δ=%.3f (W: p=%.4f %s)\nDistribution (CvM: p=%.4f %s)",
                    mean_diff, p_value_t, sig_label_t,
                    median_diff, p_value_w, sig_label_w,
                    p_value_cvm, sig_label_cvm)
  )

# Save combined data
write_csv(combined, "output/community_coverage_variance.csv")
write_csv(summary_stats, "output/coverage_bias_summary.csv")

cat("\n- output/community_coverage_variance.csv: Combined coverage and variance data\n")
cat("- output/coverage_bias_summary.csv: Summary statistics by coverage group\n")

# ==============================================================================
# BIAS VISUALIZATIONS
# ==============================================================================
cat("\n=== Creating bias assessment plots ===\n")

# Prepare data for plotting
plot_data_long = comparison_data %>%
  select(plot_id, coverage_group, env, codist, spa, env_prop, codist_prop, spa_prop) %>%
  pivot_longer(cols = c(env, codist, spa),
               names_to = "component", values_to = "raw_value") %>%
  pivot_longer(cols = c(env_prop, codist_prop, spa_prop),
               names_to = "component_prop", values_to = "prop_value") %>%
  filter(gsub("_prop", "", component_prop) == component) %>%
  mutate(
    component_label = factor(case_when(
      component == "env" ~ "Environment",
      component == "codist" ~ "Species associations",
      component == "spa" ~ "Spatial"
    ), levels = c("Environment", "Species associations", "Spatial"))
  )

# Add test results for annotations
test_results_raw = test_results %>%
  filter(variance_type == "Raw variance") %>%
  mutate(component_label = factor(component, levels = c("Environment", "Species associations", "Spatial")))

test_results_prop = test_results %>%
  filter(variance_type == "Proportional variance") %>%
  mutate(component_label = factor(component, levels = c("Environment", "Species associations", "Spatial")))

# Calculate summary statistics for overlaying mean and median points
summary_points = plot_data_long %>%
  group_by(coverage_group, component_label) %>%
  summarize(
    mean_raw = mean(raw_value, na.rm = TRUE),
    median_raw = median(raw_value, na.rm = TRUE),
    mean_prop = mean(prop_value, na.rm = TRUE),
    median_prop = median(prop_value, na.rm = TRUE),
    .groups = "drop"
  )

# Violin plots for raw variance with mean/median points and test annotations
p_violin_raw = ggplot(plot_data_long, aes(x = coverage_group, y = raw_value, fill = coverage_group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_point(data = summary_points, aes(x = coverage_group, y = mean_raw),
             color = "red", size = 3, shape = 18, inherit.aes = FALSE) +
  geom_point(data = summary_points, aes(x = coverage_group, y = median_raw),
             color = "blue", size = 3, shape = 16, inherit.aes = FALSE) +
  facet_wrap(~ component_label, scales = "free_y") +
  geom_text(data = test_results_raw,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.2, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("All data" = "grey70", "≥75% coverage" = "#4ECDC4")) +
  labs(x = NULL, y = "Raw variance",
       title = "Variance components: All data vs ≥75% trait coverage (raw values)",
       subtitle = "Red diamonds = mean, Blue circles = median") +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Violin plots for proportional variance with mean/median points and test annotations
p_violin_prop = ggplot(plot_data_long, aes(x = coverage_group, y = prop_value, fill = coverage_group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_point(data = summary_points, aes(x = coverage_group, y = mean_prop),
             color = "red", size = 3, shape = 18, inherit.aes = FALSE) +
  geom_point(data = summary_points, aes(x = coverage_group, y = median_prop),
             color = "blue", size = 3, shape = 16, inherit.aes = FALSE) +
  facet_wrap(~ component_label, scales = "free_y") +
  geom_text(data = test_results_prop,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.2, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("All data" = "grey70", "≥75% coverage" = "#4ECDC4")) +
  labs(x = NULL, y = "Proportional variance (dominance)",
       title = "Variance components: All data vs ≥75% trait coverage (proportional)",
       subtitle = "Red diamonds = mean, Blue circles = median") +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Combine bias plots
p_bias_combined = p_violin_raw / p_violin_prop

pdf("plot/trait_coverage_bias.pdf", width = 12, height = 8)
print(p_bias_combined)
dev.off()

cat("- plot/trait_coverage_bias.pdf: Community-level bias assessment visualization\n")

# ==============================================================================
# SPECIES-LEVEL BIAS ASSESSMENT
# ==============================================================================
cat("\n=== Species-level bias assessment ===\n")

# Get species variance components
species_variance = res$internals$Species %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  mutate(
    total = env + codist + spa,
    env_prop = env / total,
    codist_prop = codist / total,
    spa_prop = spa / total
  )

# Identify species with and without trait data
species_with_data = species_with_traits
species_without_data = setdiff(species_variance$species, species_with_traits)

cat(sprintf("\nSpecies with trait data: %d\n", length(species_with_data)))
cat(sprintf("Species without trait data: %d\n", length(species_without_data)))

# Create comparison dataset for species
all_species = species_variance %>%
  mutate(trait_group = "All species")

species_with_traits_df = species_variance %>%
  filter(species %in% species_with_data) %>%
  mutate(trait_group = "Species with traits")

species_comparison = bind_rows(all_species, species_with_traits_df)

# Statistical tests for species
cat("\n=== Statistical tests (All species vs Species with traits) ===\n")

# T-tests
test_sp_env = t.test(all_species$env, species_with_traits_df$env)
test_sp_codist = t.test(all_species$codist, species_with_traits_df$codist)
test_sp_spa = t.test(all_species$spa, species_with_traits_df$spa)
test_sp_env_prop = t.test(all_species$env_prop, species_with_traits_df$env_prop)
test_sp_codist_prop = t.test(all_species$codist_prop, species_with_traits_df$codist_prop)
test_sp_spa_prop = t.test(all_species$spa_prop, species_with_traits_df$spa_prop)

# Wilcoxon tests
wilcox_sp_env = wilcox.test(all_species$env, species_with_traits_df$env)
wilcox_sp_codist = wilcox.test(all_species$codist, species_with_traits_df$codist)
wilcox_sp_spa = wilcox.test(all_species$spa, species_with_traits_df$spa)
wilcox_sp_env_prop = wilcox.test(all_species$env_prop, species_with_traits_df$env_prop)
wilcox_sp_codist_prop = wilcox.test(all_species$codist_prop, species_with_traits_df$codist_prop)
wilcox_sp_spa_prop = wilcox.test(all_species$spa_prop, species_with_traits_df$spa_prop)

# Cramér-von Mises tests
cvm_sp_env = cvm.test(species_with_traits_df$env, all_species$env)
cvm_sp_codist = cvm.test(species_with_traits_df$codist, all_species$codist)
cvm_sp_spa = cvm.test(species_with_traits_df$spa, all_species$spa)
cvm_sp_env_prop = cvm.test(species_with_traits_df$env_prop, all_species$env_prop)
cvm_sp_codist_prop = cvm.test(species_with_traits_df$codist_prop, all_species$codist_prop)
cvm_sp_spa_prop = cvm.test(species_with_traits_df$spa_prop, all_species$spa_prop)

cat("\n--- Raw variance (t-tests for means) ---\n")
cat(sprintf("Environment: t=%.3f, p=%.4f %s\n", test_sp_env$statistic, test_sp_env$p.value,
            ifelse(test_sp_env$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution: t=%.3f, p=%.4f %s\n", test_sp_codist$statistic, test_sp_codist$p.value,
            ifelse(test_sp_codist$p.value < 0.05, "***", "")))
cat(sprintf("Spatial: t=%.3f, p=%.4f %s\n", test_sp_spa$statistic, test_sp_spa$p.value,
            ifelse(test_sp_spa$p.value < 0.05, "***", "")))

cat("\n--- Raw variance (Wilcoxon tests for medians) ---\n")
cat(sprintf("Environment: W=%.1f, p=%.4f %s\n", wilcox_sp_env$statistic, wilcox_sp_env$p.value,
            ifelse(wilcox_sp_env$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution: W=%.1f, p=%.4f %s\n", wilcox_sp_codist$statistic, wilcox_sp_codist$p.value,
            ifelse(wilcox_sp_codist$p.value < 0.05, "***", "")))
cat(sprintf("Spatial: W=%.1f, p=%.4f %s\n", wilcox_sp_spa$statistic, wilcox_sp_spa$p.value,
            ifelse(wilcox_sp_spa$p.value < 0.05, "***", "")))

cat("\n--- Proportional variance (t-tests for means) ---\n")
cat(sprintf("Environment prop: t=%.3f, p=%.4f %s\n", test_sp_env_prop$statistic, test_sp_env_prop$p.value,
            ifelse(test_sp_env_prop$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution prop: t=%.3f, p=%.4f %s\n", test_sp_codist_prop$statistic, test_sp_codist_prop$p.value,
            ifelse(test_sp_codist_prop$p.value < 0.05, "***", "")))
cat(sprintf("Spatial prop: t=%.3f, p=%.4f %s\n", test_sp_spa_prop$statistic, test_sp_spa_prop$p.value,
            ifelse(test_sp_spa_prop$p.value < 0.05, "***", "")))

cat("\n--- Proportional variance (Wilcoxon tests for medians) ---\n")
cat(sprintf("Environment prop: W=%.1f, p=%.4f %s\n", wilcox_sp_env_prop$statistic, wilcox_sp_env_prop$p.value,
            ifelse(wilcox_sp_env_prop$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution prop: W=%.1f, p=%.4f %s\n", wilcox_sp_codist_prop$statistic, wilcox_sp_codist_prop$p.value,
            ifelse(wilcox_sp_codist_prop$p.value < 0.05, "***", "")))
cat(sprintf("Spatial prop: W=%.1f, p=%.4f %s\n", wilcox_sp_spa_prop$statistic, wilcox_sp_spa_prop$p.value,
            ifelse(wilcox_sp_spa_prop$p.value < 0.05, "***", "")))

cat("\n--- Raw variance (Cramér-von Mises tests for distributions) ---\n")
cat(sprintf("Environment: CvM=%.4f, p=%.4f %s\n", cvm_sp_env$statistic, cvm_sp_env$p.value,
            ifelse(cvm_sp_env$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution: CvM=%.4f, p=%.4f %s\n", cvm_sp_codist$statistic, cvm_sp_codist$p.value,
            ifelse(cvm_sp_codist$p.value < 0.05, "***", "")))
cat(sprintf("Spatial: CvM=%.4f, p=%.4f %s\n", cvm_sp_spa$statistic, cvm_sp_spa$p.value,
            ifelse(cvm_sp_spa$p.value < 0.05, "***", "")))

cat("\n--- Proportional variance (Cramér-von Mises tests for distributions) ---\n")
cat(sprintf("Environment prop: CvM=%.4f, p=%.4f %s\n", cvm_sp_env_prop$statistic, cvm_sp_env_prop$p.value,
            ifelse(cvm_sp_env_prop$p.value < 0.05, "***", "")))
cat(sprintf("Co-distribution prop: CvM=%.4f, p=%.4f %s\n", cvm_sp_codist_prop$statistic, cvm_sp_codist_prop$p.value,
            ifelse(cvm_sp_codist_prop$p.value < 0.05, "***", "")))
cat(sprintf("Spatial prop: CvM=%.4f, p=%.4f %s\n", cvm_sp_spa_prop$statistic, cvm_sp_spa_prop$p.value,
            ifelse(cvm_sp_spa_prop$p.value < 0.05, "***", "")))

# Calculate mean and median differences for species (with NA removal)
mean_diff_sp_env = mean(all_species$env, na.rm = TRUE) - mean(species_with_traits_df$env, na.rm = TRUE)
mean_diff_sp_codist = mean(all_species$codist, na.rm = TRUE) - mean(species_with_traits_df$codist, na.rm = TRUE)
mean_diff_sp_spa = mean(all_species$spa, na.rm = TRUE) - mean(species_with_traits_df$spa, na.rm = TRUE)
mean_diff_sp_env_prop = mean(all_species$env_prop, na.rm = TRUE) - mean(species_with_traits_df$env_prop, na.rm = TRUE)
mean_diff_sp_codist_prop = mean(all_species$codist_prop, na.rm = TRUE) - mean(species_with_traits_df$codist_prop, na.rm = TRUE)
mean_diff_sp_spa_prop = mean(all_species$spa_prop, na.rm = TRUE) - mean(species_with_traits_df$spa_prop, na.rm = TRUE)

median_diff_sp_env = median(all_species$env, na.rm = TRUE) - median(species_with_traits_df$env, na.rm = TRUE)
median_diff_sp_codist = median(all_species$codist, na.rm = TRUE) - median(species_with_traits_df$codist, na.rm = TRUE)
median_diff_sp_spa = median(all_species$spa, na.rm = TRUE) - median(species_with_traits_df$spa, na.rm = TRUE)
median_diff_sp_env_prop = median(all_species$env_prop, na.rm = TRUE) - median(species_with_traits_df$env_prop, na.rm = TRUE)
median_diff_sp_codist_prop = median(all_species$codist_prop, na.rm = TRUE) - median(species_with_traits_df$codist_prop, na.rm = TRUE)
median_diff_sp_spa_prop = median(all_species$spa_prop, na.rm = TRUE) - median(species_with_traits_df$spa_prop, na.rm = TRUE)

# Store species t-test, Wilcoxon, and Cramér-von Mises results
test_results_species = data.frame(
  component = c("Environment", "Species associations", "Spatial",
                "Environment", "Species associations", "Spatial"),
  variance_type = c(rep("Raw variance", 3), rep("Proportional variance", 3)),
  p_value_t = c(test_sp_env$p.value, test_sp_codist$p.value, test_sp_spa$p.value,
                test_sp_env_prop$p.value, test_sp_codist_prop$p.value, test_sp_spa_prop$p.value),
  p_value_w = c(wilcox_sp_env$p.value, wilcox_sp_codist$p.value, wilcox_sp_spa$p.value,
                wilcox_sp_env_prop$p.value, wilcox_sp_codist_prop$p.value, wilcox_sp_spa_prop$p.value),
  p_value_cvm = c(cvm_sp_env$p.value, cvm_sp_codist$p.value, cvm_sp_spa$p.value,
                  cvm_sp_env_prop$p.value, cvm_sp_codist_prop$p.value, cvm_sp_spa_prop$p.value),
  mean_diff = c(mean_diff_sp_env, mean_diff_sp_codist, mean_diff_sp_spa,
                mean_diff_sp_env_prop, mean_diff_sp_codist_prop, mean_diff_sp_spa_prop),
  median_diff = c(median_diff_sp_env, median_diff_sp_codist, median_diff_sp_spa,
                  median_diff_sp_env_prop, median_diff_sp_codist_prop, median_diff_sp_spa_prop)
) %>%
  mutate(
    sig_label_t = case_when(
      p_value_t < 0.001 ~ "***",
      p_value_t < 0.01 ~ "**",
      p_value_t < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_label_w = case_when(
      p_value_w < 0.001 ~ "***",
      p_value_w < 0.01 ~ "**",
      p_value_w < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    sig_label_cvm = case_when(
      p_value_cvm < 0.001 ~ "***",
      p_value_cvm < 0.01 ~ "**",
      p_value_cvm < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label = sprintf("Mean Δ=%.3f (t: p=%.4f %s)\nMedian Δ=%.3f (W: p=%.4f %s)\nDistribution (CvM: p=%.4f %s)",
                    mean_diff, p_value_t, sig_label_t,
                    median_diff, p_value_w, sig_label_w,
                    p_value_cvm, sig_label_cvm)
  )

# Summary statistics for species
summary_stats_species = species_comparison %>%
  group_by(trait_group) %>%
  summarize(
    n = n(),
    env_mean = mean(env, na.rm = TRUE), env_median = median(env, na.rm = TRUE), env_sd = sd(env, na.rm = TRUE),
    codist_mean = mean(codist, na.rm = TRUE), codist_median = median(codist, na.rm = TRUE), codist_sd = sd(codist, na.rm = TRUE),
    spa_mean = mean(spa, na.rm = TRUE), spa_median = median(spa, na.rm = TRUE), spa_sd = sd(spa, na.rm = TRUE),
    env_prop_mean = mean(env_prop, na.rm = TRUE), env_prop_median = median(env_prop, na.rm = TRUE), env_prop_sd = sd(env_prop, na.rm = TRUE),
    codist_prop_mean = mean(codist_prop, na.rm = TRUE), codist_prop_median = median(codist_prop, na.rm = TRUE), codist_prop_sd = sd(codist_prop, na.rm = TRUE),
    spa_prop_mean = mean(spa_prop, na.rm = TRUE), spa_prop_median = median(spa_prop, na.rm = TRUE), spa_prop_sd = sd(spa_prop, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Summary statistics by species group ===\n")
print(summary_stats_species)

# Save species data
write_csv(species_variance %>%
            mutate(has_traits = species %in% species_with_data),
          "output/species_variance_trait_status.csv")
write_csv(summary_stats_species, "output/species_bias_summary.csv")

cat("\n- output/species_variance_trait_status.csv: Species variance with trait status\n")
cat("- output/species_bias_summary.csv: Summary statistics by species group\n")

# ==============================================================================
# SPECIES-LEVEL VISUALIZATIONS
# ==============================================================================
cat("\n=== Creating species-level bias plots ===\n")

# Prepare data for plotting
plot_data_species_long = species_comparison %>%
  select(species, trait_group, env, codist, spa, env_prop, codist_prop, spa_prop) %>%
  pivot_longer(cols = c(env, codist, spa),
               names_to = "component", values_to = "raw_value") %>%
  pivot_longer(cols = c(env_prop, codist_prop, spa_prop),
               names_to = "component_prop", values_to = "prop_value") %>%
  filter(gsub("_prop", "", component_prop) == component) %>%
  mutate(
    component_label = factor(case_when(
      component == "env" ~ "Environment",
      component == "codist" ~ "Species associations",
      component == "spa" ~ "Spatial"
    ), levels = c("Environment", "Species associations", "Spatial"))
  )

# Add test results for annotations
test_results_sp_raw = test_results_species %>%
  filter(variance_type == "Raw variance") %>%
  mutate(component_label = factor(component, levels = c("Environment", "Species associations", "Spatial")))

test_results_sp_prop = test_results_species %>%
  filter(variance_type == "Proportional variance") %>%
  mutate(component_label = factor(component, levels = c("Environment", "Species associations", "Spatial")))

# Calculate summary statistics for overlaying mean and median points
summary_points_sp = plot_data_species_long %>%
  group_by(trait_group, component_label) %>%
  summarize(
    mean_raw = mean(raw_value, na.rm = TRUE),
    median_raw = median(raw_value, na.rm = TRUE),
    mean_prop = mean(prop_value, na.rm = TRUE),
    median_prop = median(prop_value, na.rm = TRUE),
    .groups = "drop"
  )

# Violin plots for raw variance with mean/median points and test annotations
p_sp_violin_raw = ggplot(plot_data_species_long, aes(x = trait_group, y = raw_value, fill = trait_group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_point(data = summary_points_sp, aes(x = trait_group, y = mean_raw),
             color = "red", size = 3, shape = 18, inherit.aes = FALSE) +
  geom_point(data = summary_points_sp, aes(x = trait_group, y = median_raw),
             color = "blue", size = 3, shape = 16, inherit.aes = FALSE) +
  facet_wrap(~ component_label, scales = "free_y") +
  geom_text(data = test_results_sp_raw,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.2, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("All species" = "grey70", "Species with traits" = "#95E1D3")) +
  labs(x = NULL, y = "Raw variance",
       title = "Species variance components: All species vs Species with trait data (raw values)",
       subtitle = "Red diamonds = mean, Blue circles = median") +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Violin plots for proportional variance with mean/median points and test annotations
p_sp_violin_prop = ggplot(plot_data_species_long, aes(x = trait_group, y = prop_value, fill = trait_group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_point(data = summary_points_sp, aes(x = trait_group, y = mean_prop),
             color = "red", size = 3, shape = 18, inherit.aes = FALSE) +
  geom_point(data = summary_points_sp, aes(x = trait_group, y = median_prop),
             color = "blue", size = 3, shape = 16, inherit.aes = FALSE) +
  facet_wrap(~ component_label, scales = "free_y") +
  geom_text(data = test_results_sp_prop,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.2, size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("All species" = "grey70", "Species with traits" = "#95E1D3")) +
  labs(x = NULL, y = "Proportional variance (dominance)",
       title = "Species variance components: All species vs Species with trait data (proportional)",
       subtitle = "Red diamonds = mean, Blue circles = median") +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Combine species bias plots
p_sp_bias_combined = p_sp_violin_raw / p_sp_violin_prop

pdf("plot/species_trait_coverage_bias.pdf", width = 12, height = 8)
print(p_sp_bias_combined)
dev.off()

cat("- plot/species_trait_coverage_bias.pdf: Species-level bias assessment visualization\n")

cat("\n=== Analysis complete ===\n")
