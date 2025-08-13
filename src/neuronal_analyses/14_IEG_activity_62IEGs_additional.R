---
  title: "Analysis of IEG Activation in tuberal hypothalamus of MYT1L mice after social operant conditioning - additional analysis on time in zone/entries in zone"
author: "Jungeun Ji"
date: "2025-08-12"
output:
  ---
  
  # Set up environment
  ```{r message=FALSE, warning=TRUE}
# Change this to your directory
.libPaths("/ref/jdlab/software/r-envs/simona/r-packages/4.2/")
setwd("/scratch/jdlab/jungeun/")

library(Seurat)
library(patchwork)
library(corrplot)
library(tidyverse)
library(scCustomize)

set.seed(420)

source("color_utils.R")
```


# Load in data
```{r load in data}
combined_neuronal <- readRDS("data/250228_combined.neuronal_1-5_5pct.rds")
combined_neuronal
```

# Use of only core list of IEGs (Fos, Jun, Arc)
```{r all iegs list}
iegs <- c("Fos", "Fosb", "Fosl1", "Fosl2", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Egr3", "Egr4", "Nr4a1", "Nr4a2", "Nr4a3", "Arc", "Homer1", "Rheb", "Rgs2", "Plk2", "Ptgs2", "Bdnf", "Inhba", "Nptx2", "Plat", "Nrn1", "Myc", "Dusp1", "Dusp5", "Dusp6", "Pcdh8", "Gadd45b", "Trib1", "Gem", "Btg2", "Ier2", "Npas4", "Rasd1", "Crem", "Mbnl2", "Arf4", "Gadd45g", "Arih1", "Nup98", "Ppp1r15a", "Fbxo33", "Per1", "Per2", "Maff", "Zfp36", "Srf", "Mcl1", "Il6", "Atf3", "Rcan1", "Ncoa7", "Cxcl2", "Bhlhe40", "Slc2a3", "Nfkbia", "Ier3", "Sgk1", "Klf6", "Klf10", "Nfkbiz", "Flg", "Gbp2b", "Tnfaip3", "Cebpd", "Hbegf", "Ldlr", "Tsc22d1", "F3", "Ccl2", "Csrnp1", "Pmaip1", "Zfp36l2", "Plau", "Ccl5", "Saa3", "Tnf", "Irf1", "Cd83", "Map3k8", "Socs3", "Il1a", "Il12b", "Il1b", "Sod2", "Pim1", "Peli1", "Tlr2", "Noct", "Bcl3", "Ifit2", "Icam1", "Ifit1", "Tnfsf9", "Ccrl2", "Cxcl10", "Gbp2", "Il10", "Clec4e", "Acod1", "Mmp13", "Cxcl11", "Il23a", "Arhgef3", "Serpine1", "Traf1", "Vcam1", "Ackr4", "Marcksl1", "Nfkbid", "Ikbke", "Ccl12", "Ifit3", "Cebpb", "Zfp36l1", "Txnip", "Nfib", "Hes1", "Pias1", "Klf2", "Cd69", "Dusp2", "Wee1", "Thbs1", "Sik1", "Gdf15", "Ier5", "Rgs1", "Id2", "Apold1")
```

``{r filtering iegs}
# Fetch expression data for all IEGs
ieg_expression <- FetchData(combined_neuronal, vars = iegs)

num_cells_per_ieg <- colSums(ieg_expression > 0)

# Create a data frame with IEG names and the number of cells expressing each IEG
ieg_cell_count_df <- data.frame(IEG = iegs, NumCells = num_cells_per_ieg)

# Filter the IEGs based on the number of cells expressing each IEG
min_cells <- 300
max_cells <- 10000
filtered_iegs <- ieg_cell_count_df$IEG[(ieg_cell_count_df$NumCells >= min_cells) & (ieg_cell_count_df$NumCells <= max_cells)]

# Print the filtered IEG list
cat(paste(filtered_iegs, collapse = '", "'))
```


# Regressions of Behavior by Pct activated
## Rewards
```{r regression rewards by pct activated}
Idents(combined_neuronal) <- combined_neuronal$cell_type

# Fetch expression data for all IEGs
ieg_expression <- FetchData(combined_neuronal, vars = iegs)

# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  TotalRewards = combined_neuronal$Avg_Total_Rewards,
  TimeinZone= combined_neuronal$Avg_Exp_Zone_Time,
  EntriesinZone = combined_neuronal$Avg_Exp_Zone_Entries,
  ieg_expression
)

# Calculate the percentage of cells at least one IEGs within each cluster, Sex, and sample
cluster_ieg_sample_percentages <- cluster_ieg_sample_df %>%
  group_by(Cluster, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Run additional regression on time in zone and entries into zone
## Time in zone
cluster_regression_time_results <- cluster_ieg_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(TimeinZone ~ Activated_Percentage, data = .x)),
    summary = map(lm_model, summary),
    coefficients = map(lm_model, coef),
    p_value = map_dbl(summary, ~ .x$coefficients[2, 4]),
    r_squared = map_dbl(summary, ~ .x$r.squared)
  ) %>%
  dplyr::mutate(
    Slope = map_dbl(coefficients, ~ .x[2]),
    Intercept = map_dbl(coefficients, ~ .x[1])
  ) %>%
  dplyr::select(Cluster, Slope, Intercept, p_value, r_squared)

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values <- p.adjust(cluster_regression_time_results$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_regression_time_results$AdjustedPValue <- adjusted_p_values
cluster_regression_time_results <- cluster_regression_time_results %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_regression_time_results)

# # Write the results to a CSV file
write.csv(cluster_regression_time_results, file = "results/IEGactivity/250812_62IEGs_timeinzone_byactivationpercent.csv", row.names = FALSE)

## Entry in zone
cluster_regression_entry_results <- cluster_ieg_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(EntriesinZone ~ Activated_Percentage, data = .x)),
    summary = map(lm_model, summary),
    coefficients = map(lm_model, coef),
    p_value = map_dbl(summary, ~ .x$coefficients[2, 4]),
    r_squared = map_dbl(summary, ~ .x$r.squared)
  ) %>%
  dplyr::mutate(
    Slope = map_dbl(coefficients, ~ .x[2]),
    Intercept = map_dbl(coefficients, ~ .x[1])
  ) %>%
  dplyr::select(Cluster, Slope, Intercept, p_value, r_squared)

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values <- p.adjust(cluster_regression_entry_results$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_regression_entry_results$AdjustedPValue <- adjusted_p_values
cluster_regression_entry_results <- cluster_regression_entry_results %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_regression_entry_results)

# # Write the results to a CSV file
write.csv(cluster_regression_entry_results, file = "results/IEGactivity/250812_62IEGs_ntriesinzone_byactivationpercent.csv", row.names = FALSE)
```

```{r scatters rewards/timeinzone/entriesinzone by pct activated}
# Time in Zone
plot_cluster_time <- function(cluster, cluster_regression_results, cluster_ieg_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- cluster_regression_results[cluster_regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$TimeinZone,
       #main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Time in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l"
  )
  
  
  # Fit line
  fit <- lm(TimeinZone ~ Activated_Percentage, data = cluster_df)
  
  # Add the fit line
  abline(fit, col = "black", lwd = 2)
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("R² =", round(cluster_results$r_squared, 3)),
    paste("p-value =", format.pval(cluster_results$p_value, digits = 3)),
    paste("Adj. p-value =", format.pval(cluster_results$AdjustedPValue, digits = 3))
  )
  
  # Add legend
  legend("topleft",
         legend = legend_labels,
         col = "black",
         cex = 1.2,
         bty = "n"
  )
  
  # Add title with cluster name and coefficients
  title_text <- paste(
    "Cluster", cluster, "\n",
    "Intercept =", round(cluster_results$Intercept, 3), "\n",
    "Slope =", round(cluster_results$Slope, 3)
  )
  title(main = title_text, cex.main = 1.5)
}

## Get all unique cluster numbers
clusters <- unique(cluster_regression_time_results$Cluster)

# Create a directory to save the PDFs if it doesn't exist
dir.create("figures/IEGactivity/TimeinZoneActivatedPercentage_62IEGs_250812/", showWarnings = FALSE)

## Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/TimeinZoneActivatedPercentage_62IEGs_250812/cluster_", cluster, "_plot.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_time(cluster, cluster_regression_time_results, cluster_ieg_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'TimeinZoneActivatedPercentage' directory.\n")

# Entries in zone
plot_cluster_entry <- function(cluster, cluster_regression_results, cluster_ieg_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- cluster_regression_results[cluster_regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$EntriesinZone,
       #main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Entries in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l"
  )
  
  
  # Fit line
  fit <- lm(EntriesinZone ~ Activated_Percentage, data = cluster_df)
  
  # Add the fit line
  abline(fit, col = "black", lwd = 2)
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("R² =", round(cluster_results$r_squared, 3)),
    paste("p-value =", format.pval(cluster_results$p_value, digits = 3)),
    paste("Adj. p-value =", format.pval(cluster_results$AdjustedPValue, digits = 3))
  )
  
  # Add legend
  legend("topleft",
         legend = legend_labels,
         col = "black",
         cex = 1.2,
         bty = "n"
  )
  
  # Add title with cluster name and coefficients
  title_text <- paste(
    "Cluster", cluster, "\n",
    "Intercept =", round(cluster_results$Intercept, 3), "\n",
    "Slope =", round(cluster_results$Slope, 3)
  )
  title(main = title_text, cex.main = 1.5)
}

## Get all unique cluster numbers
clusters <- unique(cluster_regression_entry_results$Cluster)

# Create a directory to save the PDFs if it doesn't exist
dir.create("figures/IEGactivity/EntriesinZoneActivatedPercentage_62IEGs_250812/", showWarnings = FALSE)

# Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/EntriesinZoneActivatedPercentage_62IEGs_250812/cluster_", cluster, "_plot.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_entry(cluster, cluster_regression_entry_results, cluster_ieg_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'EntriesinZoneActivatedPercentage' directory.\n")
```


### with Sex
```{r regression rewards by pct activated and sex }
# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_Sex_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  TotalRewards = combined_neuronal$Avg_Total_Rewards,
  TimeinZone = combined_neuronal$Avg_Exp_Zone_Time,
  EntriesinZone = combined_neuronal$Avg_Exp_Zone_Entries,
  ieg_expression
)
# Calculate the percentage of cells expressing at least 1 IEG within each cluster, Sex, and sample
cluster_ieg_Sex_sample_percentages <- cluster_ieg_Sex_sample_df %>%
  group_by(Cluster, Sex, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Perform linear regression for each cluster
# Additional anlaysis for TimeinZone/EntriesinZone
## Time in Zone
# Perform linear regression for each cluster
regression_results_time <- cluster_ieg_Sex_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(TimeinZone ~ Activated_Percentage + Sex, data = .x)),
    summary = map(lm_model, summary),
    coefficients = map(lm_model, coef),
    r_squared = map_dbl(summary, ~ .x$r.squared),
    p_value_activated_percentage = map_dbl(summary, ~ {
      coef_matrix <- .x$coefficients
      if (nrow(coef_matrix) >= 2 && ncol(coef_matrix) >= 4) .x$coefficients[2, 4] else NA_real_
    }),
    p_value_sex = map_dbl(summary, ~ {
      coef_matrix <- .x$coefficients
      if (nrow(coef_matrix) >= 3 && ncol(coef_matrix) >= 4) .x$coefficients[3, 4] else NA_real_
    })
  ) %>%
  dplyr::mutate(
    Slope_Activated_Percentage = map_dbl(coefficients, ~ {
      if (length(.x) >= 2) .x[2] else NA_real_
    }),
    Slope_Sex = map_dbl(coefficients, ~ {
      if (length(.x) >= 3) .x[3] else NA_real_
    }),
    Intercept = map_dbl(coefficients, ~ .x[1])
  ) %>%
  dplyr::select(
    Cluster, Slope_Activated_Percentage, Slope_Sex, Intercept,
    p_value_activated_percentage, p_value_sex, r_squared
  )

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values_activated <- p.adjust(regression_results_time$p_value_activated_percentage, method = "BH")
adjusted_p_values_sex <- p.adjust(regression_results_time$p_value_sex, method = "BH")

# Add adjusted p-values to the results data frame
regression_results_time$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results_time$adjusted_p_values_sex <- adjusted_p_values_sex

regression_results_time <- regression_results_time %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results_time)

# Write the results to a CSV file
write.csv(regression_results, file = "results/IEGactivity/250812_timeinzone_byactivationpercent_bysex_62IEGs.csv", row.names = FALSE)

## Entries in Zone
# Perform linear regression for each cluster
regression_results_entry <- cluster_ieg_Sex_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(EntriesinZone ~ Activated_Percentage + Sex, data = .x)),
    summary = map(lm_model, summary),
    coefficients = map(lm_model, coef),
    r_squared = map_dbl(summary, ~ .x$r.squared),
    p_value_activated_percentage = map_dbl(summary, ~ {
      coef_matrix <- .x$coefficients
      if (nrow(coef_matrix) >= 2 && ncol(coef_matrix) >= 4) .x$coefficients[2, 4] else NA_real_
    }),
    p_value_sex = map_dbl(summary, ~ {
      coef_matrix <- .x$coefficients
      if (nrow(coef_matrix) >= 3 && ncol(coef_matrix) >= 4) .x$coefficients[3, 4] else NA_real_
    })
  ) %>%
  dplyr::mutate(
    Slope_Activated_Percentage = map_dbl(coefficients, ~ {
      if (length(.x) >= 2) .x[2] else NA_real_
    }),
    Slope_Sex = map_dbl(coefficients, ~ {
      if (length(.x) >= 3) .x[3] else NA_real_
    }),
    Intercept = map_dbl(coefficients, ~ .x[1])
  ) %>%
  dplyr::select(
    Cluster, Slope_Activated_Percentage, Slope_Sex, Intercept,
    p_value_activated_percentage, p_value_sex, r_squared
  )

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values_activated <- p.adjust(regression_results_entry$p_value_activated_percentage, method = "BH")
adjusted_p_values_sex <- p.adjust(regression_results_entry$p_value_sex, method = "BH")

# Add adjusted p-values to the results data frame
regression_results_entry$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results_entry$adjusted_p_values_sex <- adjusted_p_values_sex

regression_results_entry <- regression_results_entry %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results_entry)

# Write the results to a CSV file
write.csv(regression_results_entry, file = "results/IEGactivity/250812_entriesinzone_byactivationpercent_bysex_62IEGs.csv", row.names = FALSE)
```


```{r scatters rewards by pct activated with sex}
# Time in Zone
plot_cluster_withsex_time <- function(cluster, regression_results, cluster_ieg_Sex_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Sex_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$TimeinZone,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Time in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 1.5,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sex == "M", "#4478AB", "#ED6677")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(TimeinZone ~ Activated_Percentage + Sex, data = cluster_df)
  fit_male <- lm(TimeinZone ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "M", ])
  fit_female <- lm(TimeinZone ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "F", ])
  
  # Add the fit lines
  abline(fit_all, col = "black", lwd = 2)
  abline(fit_male, col = "#4478AB", lwd = 2)
  abline(fit_female, col = "#ED6677", lwd = 2)
  
  # Calculate R-squared values for each fit
  r2_all <- summary(fit_all)$r.squared
  r2_male <- summary(fit_male)$r.squared
  r2_female <- summary(fit_female)$r.squared
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("All (R² =", round(r2_all, 3), ")"),
    paste("Male (R² =", round(r2_male, 3), ")"),
    paste("Female (R² =", round(r2_female, 3), ")"),
    #paste("Activated % Slope =", round(cluster_results$Slope_Activated_Percentage, 3)),
    #paste("Sex Slope =", round(cluster_results$Slope_Sex, 3)),
    #paste("Intercept =", round(cluster_results$Intercept, 3)),
    paste("Activated % p-value =", format.pval(cluster_results$p_value_activated_percentage, digits = 3)),
    paste("Sex p-value =", format.pval(cluster_results$p_value_sex, digits = 3))
    #paste("Activated % Adj. p-value =", format.pval(cluster_results$adjusted_p_values_activated, digits = 3)),
    #paste("Sex Adj. p-value =", format.pval(cluster_results$adjusted_p_values_sex, digits = 3))
  )
  
  # Add legend
  legend("topright",
         legend = legend_labels,
         col = c("black", "#4478AB", "#ED6677", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 2)),
         cex = 1.2,
         bty = "n"
  )
}

# Get all unique cluster numbers
clusters <- unique(regression_results_time$Cluster)

# Iterate through clusters
for (cluster in clusters) {
  # 
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/TimeinZoneActivatedPercentage_62IEGs_250812/cluster_", cluster, "_plot_bySex.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_withsex_time(cluster, regression_results_time, cluster_ieg_Sex_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'TimeinZoneActivatedPercentage' directory.\n")

# Entries in Zone
plot_cluster_withsex_entry <- function(cluster, regression_results, cluster_ieg_Sex_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Sex_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$EntriesinZone,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Entries in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 1.5,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sex == "M", "#4478AB", "#ED6677")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(EntriesinZone ~ Activated_Percentage + Sex, data = cluster_df)
  fit_male <- lm(EntriesinZone ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "M", ])
  fit_female <- lm(EntriesinZone ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "F", ])
  
  # Add the fit lines
  abline(fit_all, col = "black", lwd = 2)
  abline(fit_male, col = "#4478AB", lwd = 2)
  abline(fit_female, col = "#ED6677", lwd = 2)
  
  # Calculate R-squared values for each fit
  r2_all <- summary(fit_all)$r.squared
  r2_male <- summary(fit_male)$r.squared
  r2_female <- summary(fit_female)$r.squared
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("All (R² =", round(r2_all, 3), ")"),
    paste("Male (R² =", round(r2_male, 3), ")"),
    paste("Female (R² =", round(r2_female, 3), ")"),
    #paste("Activated % Slope =", round(cluster_results$Slope_Activated_Percentage, 3)),
    #paste("Sex Slope =", round(cluster_results$Slope_Sex, 3)),
    #paste("Intercept =", round(cluster_results$Intercept, 3)),
    paste("Activated % p-value =", format.pval(cluster_results$p_value_activated_percentage, digits = 3)),
    paste("Sex p-value =", format.pval(cluster_results$p_value_sex, digits = 3))
    #paste("Activated % Adj. p-value =", format.pval(cluster_results$adjusted_p_values_activated, digits = 3)),
    #paste("Sex Adj. p-value =", format.pval(cluster_results$adjusted_p_values_sex, digits = 3))
  )
  
  # Add legend
  legend("topright",
         legend = legend_labels,
         col = c("black", "#4478AB", "#ED6677", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 2)),
         cex = 1.2,
         bty = "n"
  )
}

# Get all unique cluster numbers
clusters <- unique(regression_results_entry$Cluster)

# Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/EntriesinZoneActivatedPercentage_62IEGs_250812/cluster_", cluster, "_plot_bySex.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_withsex_entry(cluster, regression_results_entry, cluster_ieg_Sex_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'RewardsbyActivatedPercentage' directory.\n")
```

### with Learner 
```{r regression rewards by pct activated and learner}
# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_Learner_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sample = combined_neuronal$sample,
  Sac_Learner = combined_neuronal$Sac_Learner,
  TotalRewards = combined_neuronal$Avg_Total_Rewards,
  TimeinZone = combined_neuronal$Avg_Exp_Zone_Time,
  EntriesinZone = combined_neuronal$Avg_Exp_Zone_Entries,
  ieg_expression
)

# Calculate the percentage of cells expressing at least 1 IEG within each cluster, Sac_Learner, and sample
cluster_ieg_Learner_sample_percentages <- cluster_ieg_Learner_sample_df %>%
  group_by(Cluster, Sac_Learner, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Perform linear regression for each cluster
# Run additional analysis on Time in Zone, Entries in Zone
## Time in Zone
regression_results_time <- cluster_ieg_Learner_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(TimeinZone ~ Activated_Percentage + Sac_Learner, data = .x)),
    summary = map(lm_model, summary),
    coefficients = map(lm_model, coef),
    r_squared = map_dbl(summary, ~ .x$r.squared),
    p_value_activated_percentage = map_dbl(summary, ~ {
      coef_matrix <- .x$coefficients
      if (nrow(coef_matrix) >= 2 && ncol(coef_matrix) >= 4) .x$coefficients[2, 4] else NA_real_
    }),
    p_value_sac_learner = map_dbl(summary, ~ {
      coef_matrix <- .x$coefficients
      if (nrow(coef_matrix) >= 3 && ncol(coef_matrix) >= 4) .x$coefficients[3, 4] else NA_real_
    })
  ) %>%
  dplyr::mutate(
    Slope_Activated_Percentage = map_dbl(coefficients, ~ {
      if (length(.x) >= 2) .x[2] else NA_real_
    }),
    Slope_Sac_Learner = map_dbl(coefficients, ~ {
      if (length(.x) >= 3) .x[3] else NA_real_
    }),
    Intercept = map_dbl(coefficients, ~ .x[1])
  ) %>%
  dplyr::select(
    Cluster, Slope_Activated_Percentage, Slope_Sac_Learner, Intercept,
    p_value_activated_percentage, p_value_sac_learner, r_squared
  )

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values_activated <- p.adjust(regression_results_time$p_value_activated_percentage, method = "BH")
adjusted_p_values_sac_learner <- p.adjust(regression_results_time$p_value_sac_learner, method = "BH")

# Add adjusted p-values to the results data frame
regression_results_time$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results_time$adjusted_p_values_sac_learner <- adjusted_p_values_sac_learner

regression_results_time <- regression_results_time %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results_time)

# Write the results to a CSV file
write.csv(regression_results_time, file = "results/IEGactivity/250812_timeinzone_byactivationpercent_bysaclearner_62IEGs.csv", row.names = FALSE)

## Entries in Zone
regression_results_entry <- cluster_ieg_Learner_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(EntriesinZone ~ Activated_Percentage + Sac_Learner, data = .x)),
    summary = map(lm_model, summary),
    coefficients = map(lm_model, coef),
    r_squared = map_dbl(summary, ~ .x$r.squared),
    p_value_activated_percentage = map_dbl(summary, ~ {
      coef_matrix <- .x$coefficients
      if (nrow(coef_matrix) >= 2 && ncol(coef_matrix) >= 4) .x$coefficients[2, 4] else NA_real_
    }),
    p_value_sac_learner = map_dbl(summary, ~ {
      coef_matrix <- .x$coefficients
      if (nrow(coef_matrix) >= 3 && ncol(coef_matrix) >= 4) .x$coefficients[3, 4] else NA_real_
    })
  ) %>%
  dplyr::mutate(
    Slope_Activated_Percentage = map_dbl(coefficients, ~ {
      if (length(.x) >= 2) .x[2] else NA_real_
    }),
    Slope_Sac_Learner = map_dbl(coefficients, ~ {
      if (length(.x) >= 3) .x[3] else NA_real_
    }),
    Intercept = map_dbl(coefficients, ~ .x[1])
  ) %>%
  dplyr::select(
    Cluster, Slope_Activated_Percentage, Slope_Sac_Learner, Intercept,
    p_value_activated_percentage, p_value_sac_learner, r_squared
  )

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values_activated <- p.adjust(regression_results_entry$p_value_activated_percentage, method = "BH")
adjusted_p_values_sac_learner <- p.adjust(regression_results_entry$p_value_sac_learner, method = "BH")

# Add adjusted p-values to the results data frame
regression_results_entry$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results_entry$adjusted_p_values_sac_learner <- adjusted_p_values_sac_learner

regression_results_entry <- regression_results_entry %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results_entry)

# Write the results to a CSV file
write.csv(regression_results, file = "results/IEGactivity/250812_entriesinzone_byactivationpercent_bysaclearner_62IEGs.csv", row.names = FALSE)
```


```{r scatters rewards by pct activated with learner}
# Time in Zone
plot_cluster_withSac_Learner_time <- function(cluster, regression_results, cluster_ieg_Learner_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Learner_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$TimeinZone,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Time in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 1.5,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sac_Learner == "Learner", "#547B80", "#D1D3D4")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(TimeinZone ~ Activated_Percentage + Sac_Learner, data = cluster_df)
  fit_learner <- lm(TimeinZone ~ Activated_Percentage, data = cluster_df[cluster_df$Sac_Learner == "Learner", ])
  fit_nonlearner <- lm(TimeinZone ~ Activated_Percentage, data = cluster_df[cluster_df$Sac_Learner == "Non_Learner", ])
  
  # Add the fit lines
  abline(fit_all, col = "black", lwd = 2)
  abline(fit_learner, col = "#547B80", lwd = 2)
  abline(fit_nonlearner, col = "#D1D3D4", lwd = 2)
  
  # Calculate R-squared values for each fit
  r2_all <- summary(fit_all)$r.squared
  r2_learner <- summary(fit_learner)$r.squared
  r2_nonleaner <- summary(fit_nonlearner)$r.squared
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("All (R² =", round(r2_all, 3), ")"),
    paste("Achiever (R² =", round(r2_learner, 3), ")"),
    paste("Non-Achiever (R² =", round(r2_nonleaner, 3), ")"),
    #paste("Activated % Slope =", round(cluster_results$Slope_Activated_Percentage, 3)),
    #paste("Sac_Learner Slope =", round(cluster_results$Slope_Sac_Learner, 3)),
    #paste("Intercept =", round(cluster_results$Intercept, 3)),
    paste("Activated % p-value =", format.pval(cluster_results$p_value_activated_percentage, digits = 3)),
    paste("Achiever p-value =", format.pval(cluster_results$p_value_sac_learner, digits = 3))
    #paste("Activated % Adj. p-value =", format.pval(cluster_results$adjusted_p_values_activated, digits = 3)),
    #paste("Achiever Adj. p-value =", format.pval(cluster_results$adjusted_p_values_sac_learner, digits = 3))
  )
  
  # Add legend
  legend("topright",
         legend = legend_labels,
         col = c("black", "#547B80", "#D1D3D4", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 2)),
         cex = 1.2,
         bty = "n"
  )
}

# Get all unique cluster numbers
clusters <- unique(regression_results_time$Cluster)

# Iterate through clusters
for (cluster in clusters) {
  pdf(file = paste0("figures/IEGactivity/TimeinZoneActivatedPercentage_62IEGs_250812/cluster_", cluster, "_plot_bySac_Learner.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_withSac_Learner_time(cluster, regression_results_time, cluster_ieg_Learner_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'TimeinZoneActivatedPercentage' directory.\n")

# Entries in Zone
plot_cluster_withSac_Learner_entry <- function(cluster, regression_results, cluster_ieg_Learner_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Learner_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$EntriesinZone,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Entries in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 1.5,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sac_Learner == "Learner", "#547B80", "#D1D3D4")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(EntriesinZone ~ Activated_Percentage + Sac_Learner, data = cluster_df)
  fit_learner <- lm(EntriesinZone ~ Activated_Percentage, data = cluster_df[cluster_df$Sac_Learner == "Learner", ])
  fit_nonlearner <- lm(EntriesinZone ~ Activated_Percentage, data = cluster_df[cluster_df$Sac_Learner == "Non_Learner", ])
  
  # Add the fit lines
  abline(fit_all, col = "black", lwd = 2)
  abline(fit_learner, col = "#547B80", lwd = 2)
  abline(fit_nonlearner, col = "#D1D3D4", lwd = 2)
  
  # Calculate R-squared values for each fit
  r2_all <- summary(fit_all)$r.squared
  r2_learner <- summary(fit_learner)$r.squared
  r2_nonleaner <- summary(fit_nonlearner)$r.squared
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("All (R² =", round(r2_all, 3), ")"),
    paste("Achiever (R² =", round(r2_learner, 3), ")"),
    paste("Non-Achiever (R² =", round(r2_nonleaner, 3), ")"),
    #paste("Activated % Slope =", round(cluster_results$Slope_Activated_Percentage, 3)),
    #paste("Sac_Learner Slope =", round(cluster_results$Slope_Sac_Learner, 3)),
    #paste("Intercept =", round(cluster_results$Intercept, 3)),
    paste("Activated % p-value =", format.pval(cluster_results$p_value_activated_percentage, digits = 3)),
    paste("Achiever p-value =", format.pval(cluster_results$p_value_sac_learner, digits = 3))
    #paste("Activated % Adj. p-value =", format.pval(cluster_results$adjusted_p_values_activated, digits = 3)),
    #paste("Achiever Adj. p-value =", format.pval(cluster_results$adjusted_p_values_sac_learner, digits = 3))
  )
  
  # Add legend
  legend("topright",
         legend = legend_labels,
         col = c("black", "#547B80", "#D1D3D4", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 2)),
         cex = 1.2,
         bty = "n"
  )
}

# Get all unique cluster numbers
clusters <- unique(regression_results_entry$Cluster)

# Iterate through clusters
for (cluster in clusters) {
  # Hist Gaba, AD Serpinb7 Glut, ADP-MPO Trp73 Glut, LH-MH Glut, SNc-VTA-RAmb Foxa1 Dopa, SCsg Gabrr2 Gaba, TH Prkcd Grin2c Glut_1
  pdf(file = paste0("figures/IEGactivity/EntriesinZoneActivatedPercentage_62IEGs_250812/cluster_", cluster, "_plot_bySac_Learner.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_withSac_Learner_entry(cluster, regression_results_entry, cluster_ieg_Learner_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'EntriesinZoneActivatedPercentage' directory.\n")
```

### Control: by Pct Activity in All Clusters
```{r scatter rewards by pct activity in all neurons}
# First calculate the activation percentage across all neurons for each sample
overall_sample_percentages <- cluster_ieg_sample_df %>%
  group_by(Sample) %>%
  summarise(
    Total_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Total_Cells * 100) %>%
  ungroup()

# Create a plotting function for the control
## Time in zone
plot_control_time <- function(sample_data) {
  # Fit regression line
  fit <- lm(TimeinZone ~ Activated_Percentage, data = sample_data)
  fit_summary <- summary(fit)
  
  # Create the plot
  plot(sample_data$Activated_Percentage, sample_data$TimeinZone,
       main = "All Neurons",
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Time in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       # xlim = c(0, max(sample_data$Activated_Percentage)),
       bty = "l"
  )
  
  # Add the fit line
  abline(fit, col = "black", lwd = 2)
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("R² =", round(fit_summary$r.squared, 3)),
    paste("p-value =", format.pval(fit_summary$coefficients[2, 4], digits = 3)),
    paste("Slope =", round(fit_summary$coefficients[2, 1], 3)),
    paste("Intercept =", round(fit_summary$coefficients[1, 1], 3))
  )
  
  # Add legend
  legend("topleft",
         legend = legend_labels,
         col = "black",
         cex = 1.2,
         bty = "n"
  )
}

# Save the  plot
pdf(
  file = "figures/IEGactivity/TimeinZoneActivatedPercentage_62IEGs_250812/control_all_neurons_plot.pdf",
  width = 7, height = 8
)
plot_control_time(overall_sample_percentages)
dev.off()

# Entries in Zone
plot_control_entry <- function(sample_data) {
  # Fit regression line
  fit <- lm(EntriesinZone ~ Activated_Percentage, data = sample_data)
  fit_summary <- summary(fit)
  
  # Create the plot
  plot(sample_data$Activated_Percentage, sample_data$EntriesinZone,
       main = "All Neurons",
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Entries in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       # xlim = c(0, max(sample_data$Activated_Percentage)),
       bty = "l"
  )
  
  # Add the fit line
  abline(fit, col = "black", lwd = 2)
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("R² =", round(fit_summary$r.squared, 3)),
    paste("p-value =", format.pval(fit_summary$coefficients[2, 4], digits = 3)),
    paste("Slope =", round(fit_summary$coefficients[2, 1], 3)),
    paste("Intercept =", round(fit_summary$coefficients[1, 1], 3))
  )
  
  # Add legend
  legend("topleft",
         legend = legend_labels,
         col = "black",
         cex = 1.2,
         bty = "n"
  )
}

# Save the  plot
pdf(
  file = "figures/IEGactivity/EntriesinZoneActivatedPercentage_62IEGs_250812/control_all_neurons_plot.pdf",
  width = 7, height = 8
)
plot_control_entry(overall_sample_percentages)
dev.off()
```

```{r regression rewards by pct activity in all neurons}
# Time in zone
overall_stats_time <- overall_sample_percentages %>%
  {
    fit <- lm(TimeinZone ~ Activated_Percentage, data = .)
    sum_fit <- summary(fit)
    
    data.frame(
      Analysis = "All_Neurons",
      Slope = coef(fit)[2],
      Intercept = coef(fit)[1],
      p_value = sum_fit$coefficients[2, 4],
      r_squared = sum_fit$r.squared
    )
  }

# Write results to CSV
write.csv(overall_stats_time,
          file = "results/IEGactivity/250812_timeinzone_byactivationpercent_allneurons_62IEGs.csv",
          row.names = FALSE
)

# Entries in Zone
overall_stats_entry <- overall_sample_percentages %>%
  {
    fit <- lm(EntriesinZone ~ Activated_Percentage, data = .)
    sum_fit <- summary(fit)
    
    data.frame(
      Analysis = "All_Neurons",
      Slope = coef(fit)[2],
      Intercept = coef(fit)[1],
      p_value = sum_fit$coefficients[2, 4],
      r_squared = sum_fit$r.squared
    )
  }

# Write results to CSV
write.csv(overall_stats_entry,
          file = "results/IEGactivity/250812_entriesinzone_byactivationpercent_allneurons_62IEGs.csv",
          row.names = FALSE
)
```

```{r scatter rewards by pct activity in all neurons by sex}
# Calculate overall activation by sample and sex
overall_sample_sex_percentages <- cluster_ieg_sample_df %>%
  group_by(Sample, Sex) %>%
  summarise(
    Total_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Total_Cells * 100) %>%
  ungroup()

# Create a plotting function for the control with sex
## Time in Zone
plot_control_withsex_time <- function(sample_data) {
  # Fit lines for all data, male, and female
  fit_all <- lm(TimeinZone ~ Activated_Percentage + Sex, data = sample_data)
  fit_male <- lm(TimeinZone ~ Activated_Percentage, data = sample_data[sample_data$Sex == "M", ])
  fit_female <- lm(TimeinZone ~ Activated_Percentage, data = sample_data[sample_data$Sex == "F", ])
  
  # Create the plot
  plot(sample_data$Activated_Percentage, sample_data$TimeinZone,
       main = "All Neurons",
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Time in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(min(sample_data$Activated_Percentage), max(sample_data$Activated_Percentage)),
       bty = "l",
       col = ifelse(sample_data$Sex == "M", "#4478AB", "#ED6677")
  )
  
  # Add the fit lines
  abline(fit_all, col = "black", lwd = 2)
  abline(fit_male, col = "#4478AB", lwd = 2)
  abline(fit_female, col = "#ED6677", lwd = 2)
  
  # Calculate R-squared values for each fit
  r2_all <- summary(fit_all)$r.squared
  r2_male <- summary(fit_male)$r.squared
  r2_female <- summary(fit_female)$r.squared
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("All (R² =", round(r2_all, 3), ")"),
    paste("Male (R² =", round(r2_male, 3), ")"),
    paste("Female (R² =", round(r2_female, 3), ")"),
    #paste("Activated % Slope =", round(summary(fit_all)$coefficients[2, 1], 3)),
    #paste("Sex Slope =", round(summary(fit_all)$coefficients[3, 1], 3)),
    #paste("Intercept =", round(summary(fit_all)$coefficients[1, 1], 3)),
    paste("Activated % p-value =", format.pval(summary(fit_all)$coefficients[2, 4], digits = 3)),
    paste("Sex p-value =", format.pval(summary(fit_all)$coefficients[3, 4], digits = 3))
  )
  
  # Add legend
  legend("topright",
         legend = legend_labels,
         col = c("black", "#4478AB", "#ED6677", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 2)),
         cex = 1.2,
         bty = "n"
  )
}

# Save the control plot with sex
pdf(
  file = "figures/IEGactivity/TimeinZoneActivatedPercentage_62IEGs_250812/control_all_neurons_plot_bySex.pdf",
  width = 7, height = 8
)
plot_control_withsex_time(overall_sample_sex_percentages)
dev.off()

## Entries in Zone
plot_control_withsex_entry <- function(sample_data) {
  # Fit lines for all data, male, and female
  fit_all <- lm(EntriesinZone ~ Activated_Percentage + Sex, data = sample_data)
  fit_male <- lm(EntriesinZone ~ Activated_Percentage, data = sample_data[sample_data$Sex == "M", ])
  fit_female <- lm(EntriesinZone ~ Activated_Percentage, data = sample_data[sample_data$Sex == "F", ])
  
  # Create the plot
  plot(sample_data$Activated_Percentage, sample_data$EntriesinZone,
       main = "All Neurons",
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Entries in Zone",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(min(sample_data$Activated_Percentage), max(sample_data$Activated_Percentage)),
       bty = "l",
       col = ifelse(sample_data$Sex == "M", "#4478AB", "#ED6677")
  )
  
  # Add the fit lines
  abline(fit_all, col = "black", lwd = 2)
  abline(fit_male, col = "#4478AB", lwd = 2)
  abline(fit_female, col = "#ED6677", lwd = 2)
  
  # Calculate R-squared values for each fit
  r2_all <- summary(fit_all)$r.squared
  r2_male <- summary(fit_male)$r.squared
  r2_female <- summary(fit_female)$r.squared
  
  # Create legend labels with statistics
  legend_labels <- c(
    paste("All (R² =", round(r2_all, 3), ")"),
    paste("Male (R² =", round(r2_male, 3), ")"),
    paste("Female (R² =", round(r2_female, 3), ")"),
    #paste("Activated % Slope =", round(summary(fit_all)$coefficients[2, 1], 3)),
    #paste("Sex Slope =", round(summary(fit_all)$coefficients[3, 1], 3)),
    #paste("Intercept =", round(summary(fit_all)$coefficients[1, 1], 3)),
    paste("Activated % p-value =", format.pval(summary(fit_all)$coefficients[2, 4], digits = 3)),
    paste("Sex p-value =", format.pval(summary(fit_all)$coefficients[3, 4], digits = 3))
  )
  
  # Add legend
  legend("topright",
         legend = legend_labels,
         col = c("black", "#4478AB", "#ED6677", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 2)),
         cex = 1.2,
         bty = "n"
  )
}

# Save the control plot with sex
pdf(
  file = "figures/IEGactivity/EntriesinZoneActivatedPercentage_62IEGs_250812/control_all_neurons_plot_bySex.pdf",
  width = 7, height = 8
)
plot_control_withsex_entry(overall_sample_sex_percentages)
dev.off()
```

```{r regression rewards by pct activity in all neurons by sex}
# Time in Zone
overall_sex_stats_time <- overall_sample_sex_percentages %>%
  {
    # Full model with sex
    fit <- lm(TimeinZone ~ Activated_Percentage + Sex, data = .)
    sum_fit <- summary(fit)
    
    # Sex-specific models
    fit_male <- lm(TimeinZone ~ Activated_Percentage, data = .[.$Sex == "M", ])
    fit_female <- lm(TimeinZone ~ Activated_Percentage, data = .[.$Sex == "F", ])
    
    data.frame(
      Analysis = "All_Neurons",
      Slope_Activated_Percentage = coef(fit)[2],
      Slope_Sex = coef(fit)[3],
      Intercept = coef(fit)[1],
      p_value_activated_percentage = sum_fit$coefficients[2, 4],
      p_value_sex = sum_fit$coefficients[3, 4],
      r_squared_full = sum_fit$r.squared,
      r_squared_male = summary(fit_male)$r.squared,
      r_squared_female = summary(fit_female)$r.squared
    )
  }

# Write results to CSV
write.csv(overall_sex_stats_time,
          file = "results/IEGactivity/250812_timeinzone_byactivationpercent_allneurons_bysex_62IEGs.csv",
          row.names = FALSE
)

# Entries in Zone
overall_sex_stats_entry <- overall_sample_sex_percentages %>%
  {
    # Full model with sex
    fit <- lm(EntriesinZone ~ Activated_Percentage + Sex, data = .)
    sum_fit <- summary(fit)
    
    # Sex-specific models
    fit_male <- lm(EntriesinZone ~ Activated_Percentage, data = .[.$Sex == "M", ])
    fit_female <- lm(EntriesinZone ~ Activated_Percentage, data = .[.$Sex == "F", ])
    
    data.frame(
      Analysis = "All_Neurons",
      Slope_Activated_Percentage = coef(fit)[2],
      Slope_Sex = coef(fit)[3],
      Intercept = coef(fit)[1],
      p_value_activated_percentage = sum_fit$coefficients[2, 4],
      p_value_sex = sum_fit$coefficients[3, 4],
      r_squared_full = sum_fit$r.squared,
      r_squared_male = summary(fit_male)$r.squared,
      r_squared_female = summary(fit_female)$r.squared
    )
  }

# Write results to CSV
write.csv(overall_sex_stats_entry,
          file = "results/IEGactivity/250812_entriesinzone_byactivationpercent_allneurons_bysex_62IEGs.csv",
          row.names = FALSE
)
```
