---
title: "Analysis of IEG Activation in tuberal hypothalamus of MYT1L mice after social operant conditioning"
author: "Jungeun Ji"
date: "2025-07-25"
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

```{r renaming cell_types for saving}
# Saving scatterplots won't work if you have "/" in cluster name
combined_neuronal@meta.data$cell_type <- as.factor(gsub("/", "-", combined_neuronal@meta.data$cell_type))
```

# Use of only core list of IEGs (Fos, Jun, Arc)
```{r all iegs list}
iegs <- c("Fos", "Jun", "Arc")
```

```{r filtering iegs}
# Fetch expression data for all IEGs
ieg_expression <- FetchData(combined_neuronal, vars = iegs)

num_cells_per_ieg <- colSums(ieg_expression > 0)

# Create a data frame with IEG names and the number of cells expressing each IEG
ieg_cell_count_df <- data.frame(IEG = iegs, NumCells = num_cells_per_ieg)

print(ieg_cell_count_df)
# Fos     2016
# Jun      979
# Arc      718

# Filter the IEGs based on the number of cells expressing each IEG
#min_cells <- 300
#max_cells <- 10000
#filtered_iegs <- ieg_cell_count_df$IEG[(ieg_cell_count_df$NumCells >= min_cells) & (ieg_cell_count_df$NumCells <= max_cells)]

# Print the filtered IEG list
cat(paste(filtered_iegs, collapse = '", "'))
```

# Stats on Pct activated by Sample

## by Geno
```{r fishers pct activated by geno}
Idents(combined_neuronal) <- combined_neuronal@meta.data$cell_type

# Fetch expression data for all IEGs
ieg_expression <- FetchData(combined_neuronal, vars = iegs)

# Create a data frame with cluster assignments, genotype, sample, and IEG expression
cluster_ieg_genotype_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Genotype = combined_neuronal$genotype,
  Sample = combined_neuronal$sample,
  ieg_expression
)

# Calculate the percentage of cells expressing at least one IEGs within each cluster, genotype, and sample
cluster_ieg_genotype_sample_percentages <- cluster_ieg_genotype_sample_df %>%
   group_by(Cluster, Genotype, Sample) %>%
   summarise(
     Cluster_Cells = n(),
     Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1)
   ) %>%
   dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
   ungroup()

# Perform Fisher's exact test on activated percentages between genotypes within each cluster
cluster_ieg_genotype_comparison <- cluster_ieg_genotype_sample_percentages %>%
  group_by(Cluster) %>%
  summarise(
    Activated_Percentage_WT_mean = mean(Activated_Percentage[Genotype == "WT"]),
    Activated_Percentage_Het_mean = mean(Activated_Percentage[Genotype == "Het"]),
    Activated_Percentage_Difference = Activated_Percentage_WT_mean - Activated_Percentage_Het_mean,
    p_value = fisher.test(
      x = matrix(
        c(
          sum(Activated_Cells[Genotype == "WT"]),
          sum(Cluster_Cells[Genotype == "WT"] - Activated_Cells[Genotype == "WT"]),
          sum(Activated_Cells[Genotype == "Het"]),
          sum(Cluster_Cells[Genotype == "Het"] - Activated_Cells[Genotype == "Het"])
        ),
        nrow = 2
      )
    )$p.value
  )

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values <- p.adjust(cluster_ieg_genotype_comparison$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_ieg_genotype_comparison$AdjustedPValue <- adjusted_p_values

cluster_ieg_genotype_comparison <- cluster_ieg_genotype_comparison %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_ieg_genotype_comparison)


# Write the results to a CSV file
write.csv(cluster_ieg_genotype_comparison, file = "results/IEGactivity/250725_only3IEGs_activationpercent_bygeno.csv", row.names = FALSE)
```

## by Sex

```{r fishers pct activated by sex}
Idents(combined_neuronal) <- combined_neuronal@meta.data$cell_type

# Fetch expression data for all IEGs
# ieg_expression <- FetchData(combined_neuronal, vars = filtered_iegs)

# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_Sex_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  ieg_expression
)

# Calculate the percentage of cells expressing at least one IEGs within each cluster, sex, and sample
cluster_ieg_Sex_sample_percentages <- cluster_ieg_Sex_sample_df %>%
  group_by(Cluster, Sex, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()


# Perform Fisher's exact test on activated percentages between Sexs within each cluster
cluster_ieg_Sex_comparison <- cluster_ieg_Sex_sample_percentages %>%
  group_by(Cluster) %>%
  summarise(
    Activated_Percentage_M_mean = mean(Activated_Percentage[Sex == "M"]),
    Activated_Percentage_F_mean = mean(Activated_Percentage[Sex == "F"]),
    Activated_Percentage_Difference = Activated_Percentage_M_mean - Activated_Percentage_F_mean,
    p_value = fisher.test(
      x = matrix(
        c(
          sum(Activated_Cells[Sex == "M"]),
          sum(Cluster_Cells[Sex == "M"] - Activated_Cells[Sex == "M"]),
          sum(Activated_Cells[Sex == "F"]),
          sum(Cluster_Cells[Sex == "F"] - Activated_Cells[Sex == "F"])
        ),
        nrow = 2
      )
    )$p.value
  )

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values <- p.adjust(cluster_ieg_Sex_comparison$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_ieg_Sex_comparison$AdjustedPValue <- adjusted_p_values

cluster_ieg_Sex_comparison <- cluster_ieg_Sex_comparison %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_ieg_Sex_comparison)


# # Write the results to a CSV file
write.csv(cluster_ieg_Sex_comparison, file = "results/IEGactivity/250725_only3IEGs_activationpercent_bysex.csv", row.names = FALSE)
```


## by Learner 
```{r fishers pct activated by learner}
Idents(combined_neuronal) <- combined_neuronal@meta.data$cell_type

# Fetch expression data for all IEGs
# ieg_expression <- FetchData(combined_neuronal, vars = filtered_iegs)

# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_Learner_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Learner = combined_neuronal$Sac_Learner,
  Sample = combined_neuronal$sample,
  ieg_expression
)

# Calculate the percentage of cells expressing at least 1 IEGs within each cluster, Sex, and sample
cluster_ieg_Learner_sample_percentages <- cluster_ieg_Learner_sample_df %>%
  group_by(Cluster, Learner, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()


# Perform Fisher's exact test on activated percentages between Sexs within each cluster
cluster_ieg_Learner_comparison <- cluster_ieg_Learner_sample_percentages %>%
  group_by(Cluster) %>%
  summarise(
    Activated_Percentage_Learner_mean = mean(Activated_Percentage[Learner == "Learner"]),
    Activated_Percentage_Non_Learner_mean = mean(Activated_Percentage[Learner == "Non_Learner"]),
    Activated_Percentage_Difference = Activated_Percentage_Learner_mean - Activated_Percentage_Non_Learner_mean,
    p_value = fisher.test(
      x = matrix(
        c(
          sum(Activated_Cells[Learner == "Learner"]),
          sum(Cluster_Cells[Learner == "Learner"] - Activated_Cells[Learner == "Learner"]),
          sum(Activated_Cells[Learner == "Non_Learner"]),
          sum(Cluster_Cells[Learner == "Non_Learner"] - Activated_Cells[Learner == "Non_Learner"])
        ),
        nrow = 2
      )
    )$p.value
  )

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values <- p.adjust(cluster_ieg_Learner_comparison$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_ieg_Learner_comparison$AdjustedPValue <- adjusted_p_values

cluster_ieg_Learner_comparison <- cluster_ieg_Learner_comparison %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_ieg_Learner_comparison)


# # Write the results to a CSV file
write.csv(cluster_ieg_Learner_comparison, file = "results/IEGactivity/250725_only3IEGs_activationpercent_bylearner.csv", row.names = FALSE)
```

# Percentage bar plots

## by Geno
```{r pct activated by geno df}
# Fetch expression data for all IEGs
# ieg_expression <- FetchData(combined_neuronal, vars = filtered_iegs)

# Create a data frame with cluster assignments, Geno, and IEG expression
cluster_ieg_Geno_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sample = combined_neuronal@meta.data$sample,
  Geno = combined_neuronal@meta.data$genotype,
  ieg_expression
)

# Calculate the percentage of cells expressing at least 1 IEG within each cluster and sample
cluster_ieg_Geno_percentages <- cluster_ieg_Geno_df %>%
  group_by(Cluster, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Convert Cluster to a factor with desired order of levels
cluster_ieg_Geno_percentages$Cluster <- factor(cluster_ieg_Geno_percentages$Cluster,
                                               levels = sort(as.character(levels(combined_neuronal)))
)

cluster_ieg_Geno_percentages$Geno <- gsub("^(.*?)_.*$", "\\1", cluster_ieg_Geno_percentages$Sample)
```

```{r pct activated by geno calculations}
# Determine mean average per cluster for Geno
cluster_ieg_Geno_percentages_norm <- cluster_ieg_Geno_percentages %>%
  group_by(Cluster, Geno) %>%
  dplyr::mutate(Mean_Activated_Percentage = mean(Activated_Percentage)) %>%
  dplyr::mutate(SE_Mean = sd(Activated_Percentage) / sqrt(36))

WT_M_1_values <- cluster_ieg_Geno_percentages_norm[cluster_ieg_Geno_percentages_norm$Sample == "WT_M_1", c("Cluster", "Mean_Activated_Percentage")]
cluster_ieg_Geno_percentages_norm <- merge(cluster_ieg_Geno_percentages_norm, WT_M_1_values, by = "Cluster", all.x = TRUE)
colnames(cluster_ieg_Geno_percentages_norm)[colnames(cluster_ieg_Geno_percentages_norm) == "Mean_Activated_Percentage.y"] <- "Mean_Activated_Percentage_M"

cluster_ieg_Geno_percentages_norm <- cluster_ieg_Geno_percentages_norm %>%
  dplyr::mutate(Activated_Percentage_normalizedtoWT = ((Activated_Percentage / Mean_Activated_Percentage_M) - 1))

# Make sure the WT mean is around 0
cluster_ieg_Geno_percentages_norm %>%
  filter(Geno == "WT") %>%
  group_by(Cluster) %>%
  summarize(mean = mean(Activated_Percentage_normalizedtoWT))

cluster_ieg_Geno_percentages_norm <- cluster_ieg_Geno_percentages_norm %>%
  group_by(Cluster, Geno) %>%
  dplyr::mutate(Mean_Activated_Percentage_normalizedtoWT = mean(Activated_Percentage_normalizedtoWT)) %>%
  dplyr::mutate(SE_Mean_Activated_Percentage_normalizedtoWT = sd(Activated_Percentage_normalizedtoWT) / sqrt(36))
```


```{r pct activated by geno plot}
p1 <- ggplot(cluster_ieg_Geno_percentages_norm, aes(x = Mean_Activated_Percentage.x, y = Cluster, fill = Geno)) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(aes(y = Cluster, xmin = Mean_Activated_Percentage.x - SE_Mean, xmax = Mean_Activated_Percentage.x + SE_Mean),
                position = position_dodge(0.7),
                width = 0.2, size = 0.2
  ) +
  geom_point(aes(x = Activated_Percentage, y = Cluster),
             position = position_dodge(width = 0.7),
             color = "black",
             size = 0.01
  ) +
  labs(x = "Proportion", y = "Cluster") +
  theme_classic() +
  scale_fill_manual(values = palette_geno) +
  theme(legend.position = "right")

p2 <- ggplot(cluster_ieg_Geno_percentages_norm, aes(x = Mean_Activated_Percentage_normalizedtoWT, y = Cluster, fill = Geno)) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(aes(y = Cluster, xmin = Mean_Activated_Percentage_normalizedtoWT - SE_Mean_Activated_Percentage_normalizedtoWT, xmax = Mean_Activated_Percentage_normalizedtoWT + SE_Mean_Activated_Percentage_normalizedtoWT),
                position = position_dodge(0.7),
                width = 0.2, size = 0.2
  ) +
  labs(x = "Proportion", y = "") +
  theme_classic() +
  scale_fill_manual(values = palette_geno) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none"
  )


wrap_plots(p1, p2)
ggsave("figures/IEGactivity/250725-only3IEGs_activationbyGeno_neuronal_1-5_5pct_50_1.8.pdf",
       device = pdf,
       width = 12,
       height = 8
)
```

## by Sex
```{r pct activated by sex df}
# Fetch expression data for all IEGs
# ieg_expression <- FetchData(combined_neuronal, vars = filtered_iegs)

# Create a data frame with cluster assignments, sex, and IEG expression
cluster_ieg_sex_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sample = combined_neuronal@meta.data$sample,
  Sex = combined_neuronal@meta.data$sex,
  ieg_expression
)

# Calculate the percentage of cells expressing at least one IEG within each cluster and sample
cluster_ieg_sex_percentages <- cluster_ieg_sex_df %>%
  group_by(Cluster, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Convert Cluster to a factor with desired order of levels
cluster_ieg_sex_percentages$Cluster <- factor(cluster_ieg_sex_percentages$Cluster,
                                              levels = sort(as.character(levels(combined_neuronal)))
)

cluster_ieg_sex_percentages$Sex <- gsub(".*_(\\w)_.*", "\\1", cluster_ieg_sex_percentages$Sample)
```


```{r pct activated by sex calculations}
# Determine mean average per cluster for sex
cluster_ieg_sex_percentages_norm <- cluster_ieg_sex_percentages %>%
  group_by(Cluster, Sex) %>%
  dplyr::mutate(Mean_Activated_Percentage = mean(Activated_Percentage)) %>%
  dplyr::mutate(SE_Mean = sd(Activated_Percentage) / sqrt(36))

WT_M_1_values <- cluster_ieg_sex_percentages_norm[cluster_ieg_sex_percentages_norm$Sample == "WT_M_1", c("Cluster", "Mean_Activated_Percentage")]
cluster_ieg_sex_percentages_norm <- merge(cluster_ieg_sex_percentages_norm, WT_M_1_values, by = "Cluster", all.x = TRUE)
colnames(cluster_ieg_sex_percentages_norm)[colnames(cluster_ieg_sex_percentages_norm) == "Mean_Activated_Percentage.y"] <- "Mean_Activated_Percentage_M"

cluster_ieg_sex_percentages_norm <- cluster_ieg_sex_percentages_norm %>%
  dplyr::mutate(Activated_Percentage_normalizedtoM = ((Activated_Percentage / Mean_Activated_Percentage_M) - 1))

# Make sure the WT mean is around 0
cluster_ieg_sex_percentages_norm %>%
  filter(Sex == "M") %>%
  group_by(Cluster) %>%
  summarize(mean = mean(Activated_Percentage_normalizedtoM))

cluster_ieg_sex_percentages_norm <- cluster_ieg_sex_percentages_norm %>%
  group_by(Cluster, Sex) %>%
  dplyr::mutate(Mean_Activated_Percentage_normalizedtoM = mean(Activated_Percentage_normalizedtoM)) %>%
  dplyr::mutate(SE_Mean_Activated_Percentage_normalizedtoM = sd(Activated_Percentage_normalizedtoM) / sqrt(36))
```


```{r pct activated by geno plot}
p1 <- ggplot(cluster_ieg_sex_percentages_norm, aes(x = Mean_Activated_Percentage.x, y = Cluster, fill = Sex)) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(aes(y = Cluster, xmin = Mean_Activated_Percentage.x - SE_Mean, xmax = Mean_Activated_Percentage.x + SE_Mean),
                position = position_dodge(0.7),
                width = 0.2, size = 0.2
  ) +
  geom_point(aes(x = Activated_Percentage, y = Cluster),
             position = position_dodge(width = 0.7),
             color = "black",
             size = 0.01
  ) +
  labs(x = "Proportion", y = "Cluster") +
  theme_classic() +
  scale_fill_manual(values = palette_sex) +
  theme(legend.position = "right")

p2 <- ggplot(cluster_ieg_sex_percentages_norm, aes(x = Mean_Activated_Percentage_normalizedtoM, y = Cluster, fill = Sex)) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(aes(y = Cluster, xmin = Mean_Activated_Percentage_normalizedtoM - SE_Mean_Activated_Percentage_normalizedtoM, xmax = Mean_Activated_Percentage_normalizedtoM + SE_Mean_Activated_Percentage_normalizedtoM),
                position = position_dodge(0.7),
                width = 0.2, size = 0.2
  ) +
  labs(x = "Proportion", y = "") +
  theme_classic() +
  scale_fill_manual(values = palette_sex) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none"
  )


wrap_plots(p1, p2)
ggsave("figures/IEGactivity/250725-only3IEGs_activationbySex_neuronal_1-5_5pct_50_1.8.pdf",
       device = pdf,
       width = 12,
       height = 8
)
```

```{r}
significant_clusters <- cluster_ieg_Sex_comparison %>%
  filter(p_value < 0.1) %>%
  pull(Cluster) %>%
  as.character()

sig_cluster_ieg_Sex_percentages_norm <- cluster_ieg_sex_percentages_norm %>%
  dplyr::filter(Cluster %in% significant_clusters) %>%
  arrange(match(Cluster, significant_clusters))

sig_cluster_ieg_Sex_percentages_norm$Cluster <-
  factor(sig_cluster_ieg_Sex_percentages_norm$Cluster,
         levels = rev(sort(significant_clusters))
  )

sig_p1 <- ggplot(sig_cluster_ieg_Sex_percentages_norm, aes(x = Mean_Activated_Percentage.x, y = Cluster, fill = Sex)) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(aes(y = Cluster, xmin = Mean_Activated_Percentage.x - SE_Mean, xmax = Mean_Activated_Percentage.x + SE_Mean),
                position = position_dodge(0.7),
                width = 0.2, size = 0.2
  ) +
  labs(x = "Proportion", y = "Cluster") +
  theme_classic() +
  scale_fill_manual(values = palette_sex) +
  geom_point(aes(x = Activated_Percentage, y = Cluster),
             position = position_dodge(width = 0.7),
             color = "black",
             size = 0.01
  ) +
  theme(legend.position = "right")

sig_p1

ggsave("figures/IEGactivity/250725-only3IEGs_activationbySex_neuronal_1-5_5pct_50_1.8_sigonly.pdf",
       device = pdf,
       width = 6.5,
       height = 8
)
```

## by Learner
```{r pct activated by learner df}
# Fetch expression data for all IEGs
# ieg_expression <- FetchData(combined_neuronal, vars = filtered_iegs)

# Create a data frame with cluster assignments, sex, and IEG expression
cluster_ieg_Learner_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sample = combined_neuronal@meta.data$sample,
  Learner = combined_neuronal@meta.data$Sac_Learner,
  ieg_expression
)

# Calculate the percentage of cells expressing at least 1 IEG within each cluster and sample
cluster_ieg_Learner_percentages <- cluster_ieg_Learner_df %>%
  group_by(Cluster, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    Learner = dplyr::first(Learner)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()
```


```{r pct activated by learner calculations}
# Determine mean average per cluster for sex
cluster_ieg_Learner_percentages_norm <- cluster_ieg_Learner_percentages %>%
  group_by(Cluster, Learner) %>%
  dplyr::mutate(Mean_Activated_Percentage = mean(Activated_Percentage)) %>%
  dplyr::mutate(SE_Mean = sd(Activated_Percentage) / sqrt(36))

WT_M_10_values <- cluster_ieg_Learner_percentages_norm[cluster_ieg_Learner_percentages_norm$Sample == "WT_M_10", c("Cluster", "Mean_Activated_Percentage")]
cluster_ieg_Learner_percentages_norm <- merge(cluster_ieg_Learner_percentages_norm, WT_M_10_values, by = "Cluster", all.x = TRUE)
colnames(cluster_ieg_Learner_percentages_norm)[colnames(cluster_ieg_Learner_percentages_norm) == "Mean_Activated_Percentage.y"] <- "Mean_Activated_Percentage_NL"

cluster_ieg_Learner_percentages_norm <- cluster_ieg_Learner_percentages_norm %>%
  dplyr::mutate(Activated_Percentage_normalizedtoNL = ((Activated_Percentage / Mean_Activated_Percentage_NL) - 1))

# Make sure the WT mean is around 0
cluster_ieg_Learner_percentages_norm %>%
  filter(Learner == "Non_Learner") %>%
  group_by(Cluster) %>%
  summarize(mean = mean(Activated_Percentage_normalizedtoNL))

cluster_ieg_Learner_percentages_norm <- cluster_ieg_Learner_percentages_norm %>%
  group_by(Cluster, Learner) %>%
  dplyr::mutate(Mean_Activated_Percentage_normalizedtoNL = mean(Activated_Percentage_normalizedtoNL)) %>%
  dplyr::mutate(SE_Mean_Activated_Percentage_normalizedtoNL = sd(Activated_Percentage_normalizedtoNL) / sqrt(36))
```


```{r pct activated by learner plot}
# Change "Learner" label to "Achiever"
cluster_ieg_Learner_percentages_norm$Achiever = ifelse(cluster_ieg_Learner_percentages_norm$Learner == "Learner", "Achiever", "Non-Achiever")

p1 <- ggplot(cluster_ieg_Learner_percentages_norm, aes(x = Mean_Activated_Percentage.x, y = Cluster, fill = Achiever)) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(aes(y = Cluster, xmin = Mean_Activated_Percentage.x - SE_Mean, xmax = Mean_Activated_Percentage.x + SE_Mean),
                position = position_dodge(0.7),
                width = 0.2, size = 0.2
  ) +
  geom_point(aes(x = Activated_Percentage, y = Cluster),
             position = position_dodge(width = 0.7),
             color = "black",
             size = 0.01
  ) +
  labs(x = "Proportion", y = "Cluster") +
  theme_classic() +
  scale_fill_manual(values = palette_achiever) +
  theme(legend.position = "right")

p2 <- ggplot(cluster_ieg_Learner_percentages_norm, aes(x = Mean_Activated_Percentage_normalizedtoNL, y = Cluster, fill = Achiever)) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(aes(y = Cluster, xmin = Mean_Activated_Percentage_normalizedtoNL - SE_Mean_Activated_Percentage_normalizedtoNL, xmax = Mean_Activated_Percentage_normalizedtoNL + SE_Mean_Activated_Percentage_normalizedtoNL),
                position = position_dodge(0.7),
                width = 0.2, size = 0.2
  ) +
  labs(x = "Proportion", y = "") +
  theme_classic() +
  scale_fill_manual(values = palette_achiever) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none"
  )

wrap_plots(p1, p2)
ggsave("figures/IEGactivity/250725-only3IEGs_activationbyAchiever_neuronal_1-5_5pct_50_1.8.pdf",
       device = pdf,
       width = 12,
       height = 8
)
```


```{r sig clusters only plotr}
significant_clusters <- cluster_ieg_Learner_comparison %>%
  filter(p_value < 0.1) %>%
  pull(Cluster) %>%
  as.character()

sig_cluster_ieg_Learner_percentages_norm <- cluster_ieg_Learner_percentages_norm %>%
  dplyr::filter(Cluster %in% significant_clusters) %>%
  arrange(match(Cluster, significant_clusters))

sig_p1 <- ggplot(sig_cluster_ieg_Learner_percentages_norm, aes(x = Mean_Activated_Percentage.x, y = Cluster, fill = Achiever)) +
  geom_bar(
    stat = "identity", position = "dodge",
    width = 0.7
  ) +
  geom_errorbar(aes(y = Cluster, xmin = Mean_Activated_Percentage.x - SE_Mean, xmax = Mean_Activated_Percentage.x + SE_Mean),
                position = position_dodge(0.7),
                width = 0.2, size = 0.2
  ) +
  labs(x = "Proportion", y = "Cluster") +
  theme_classic() +
  scale_fill_manual(values = palette_achiever) +
  geom_point(aes(x = Activated_Percentage, y = Cluster),
             position = position_dodge(width = 0.7),
             color = "black",
             size = 0.01
  ) +
  theme(legend.position = "right")

sig_p1

ggsave("figures/IEGactivity/250725-only3IEGs_activationbyAchiever_neuronal_1-5_5pct_50_1.8_sigonly.pdf",
       device = pdf,
       width = 6.5,
       height = 8
)
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
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

cluster_regression_results <- cluster_ieg_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(TotalRewards ~ Activated_Percentage, data = .x)),
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
adjusted_p_values <- p.adjust(cluster_regression_results$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_regression_results$AdjustedPValue <- adjusted_p_values
cluster_regression_results <- cluster_regression_results %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_regression_results)

# # Write the results to a CSV file
write.csv(cluster_regression_results, file = "results/IEGactivity/250725_only3IEGs_rewards_byactivationpercent.csv", row.names = FALSE)

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
write.csv(cluster_regression_time_results, file = "results/IEGactivity/250812_only3IEGs_timeinzone_byactivationpercent.csv", row.names = FALSE)

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
write.csv(cluster_regression_entry_results, file = "results/IEGactivity/250812_only3IEGs_entriesinzone_byactivationpercent.csv", row.names = FALSE)
```

```{r scatters rewards/timeinzone/entriesinzone by pct activated}
# Total rewards
plot_cluster <- function(cluster, cluster_regression_results, cluster_ieg_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- cluster_regression_results[cluster_regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$TotalRewards,
       #main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Total Rewards",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l"
  )
  
  
  # Fit line
  fit <- lm(TotalRewards ~ Activated_Percentage, data = cluster_df)
  
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
clusters <- unique(cluster_regression_results$Cluster)

## Create a directory to save the PDFs if it doesn't exist
dir.create("figures/IEGactivity/RewardsbyActivatedPercentage_only3IEGs_250812", showWarnings = FALSE)

## Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/RewardsbyActivatedPercentage_only3IEGs_250812/cluster_", cluster, "_plot.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster(cluster, cluster_regression_results, cluster_ieg_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'RewardsbyActivatedPercentage' directory.\n")

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
dir.create("figures/IEGactivity/TimeinZoneActivatedPercentage_only3IEGs_250812/", showWarnings = FALSE)

## Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/TimeinZoneActivatedPercentage_only3IEGs_250812/cluster_", cluster, "_plot.pdf"), width = 7, height = 8)
  
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
dir.create("figures/IEGactivity/EntriesinZoneActivatedPercentage_only3IEGs_250812/", showWarnings = FALSE)

# Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/EntriesinZoneActivatedPercentage_only3IEGs_250812/cluster_", cluster, "_plot.pdf"), width = 7, height = 8)
  
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
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Perform linear regression for each cluster
regression_results <- cluster_ieg_Sex_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(TotalRewards ~ Activated_Percentage + Sex, data = .x)),
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
adjusted_p_values_activated <- p.adjust(regression_results$p_value_activated_percentage, method = "BH")
adjusted_p_values_sex <- p.adjust(regression_results$p_value_sex, method = "BH")

# Add adjusted p-values to the results data frame
regression_results$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results$adjusted_p_values_sex <- adjusted_p_values_sex

regression_results <- regression_results %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results)

# Write the results to a CSV file
write.csv(regression_results, file = "results/IEGactivity/250725_rewards_byactivationpercent_bysex.csv", row.names = FALSE)

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
write.csv(regression_results, file = "results/IEGactivity/250812_timeinzone_byactivationpercent_bysex_only3IEGs.csv", row.names = FALSE)

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
write.csv(regression_results_entry, file = "results/IEGactivity/250812_entriesinzone_byactivationpercent_bysex_only3IEGs.csv", row.names = FALSE)
```


```{r scatters rewards by pct activated with sex}
# Total rewards
plot_cluster_withsex <- function(cluster, regression_results, cluster_ieg_Sex_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Sex_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$TotalRewards,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Rewards",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 1.5,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sex == "M", "#4478AB", "#ED6677")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(TotalRewards ~ Activated_Percentage + Sex, data = cluster_df)
  fit_male <- lm(TotalRewards ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "M", ])
  fit_female <- lm(TotalRewards ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "F", ])
  
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
clusters <- unique(regression_results$Cluster)

# Iterate through clusters
for (cluster in clusters) {
  # SCsg Gabrr2 Gaba, ADP-MPO Trp73 Glut, SNc-VTA-RAmb Foxa1 Dopa not included because all activated percentage is at 0
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/RewardsbyActivatedPercentage_only3IEGs_250725/cluster_", cluster, "_plot_bySex.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_withsex(cluster, regression_results, cluster_ieg_Sex_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'RewardsbyActivatedPercentage' directory.\n")

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
  # ADP-MPO Trp73 Glut, SCsg Gabrr2 Gaba, SNc-VTA-RAmb Foxa1 Dopa not plotted
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/TimeinZoneActivatedPercentage_only3IEGs_250812/cluster_", cluster, "_plot_bySex.pdf"), width = 7, height = 8)
  
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
  # ADP-MPO Trp73 Glut, SCsg Gabrr2 Gaba, SNc-VTA-RAmb Foxa1 Dopa not plotted
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/EntriesinZoneActivatedPercentage_only3IEGs_250812/cluster_", cluster, "_plot_bySex.pdf"), width = 7, height = 8)
  
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
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Perform linear regression for each cluster
## Total rewards
regression_results <- cluster_ieg_Learner_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(TotalRewards ~ Activated_Percentage + Sac_Learner, data = .x)),
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
adjusted_p_values_activated <- p.adjust(regression_results$p_value_activated_percentage, method = "BH")
adjusted_p_values_sac_learner <- p.adjust(regression_results$p_value_sac_learner, method = "BH")

# Add adjusted p-values to the results data frame
regression_results$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results$adjusted_p_values_sac_learner <- adjusted_p_values_sac_learner

regression_results <- regression_results %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results)

# Write the results to a CSV file
write.csv(regression_results, file = "results/IEGactivity/250725_only3IEGs_rewards_byactivationpercent_bysaclearner.csv", row.names = FALSE)

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
write.csv(regression_results_time, file = "results/IEGactivity/250812_timeinzone_byactivationpercent_bysaclearner_only3IEGs.csv", row.names = FALSE)

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
write.csv(regression_results, file = "results/IEGactivity/250812_entriesinzone_byactivationpercent_bysaclearner_only3IEGs.csv", row.names = FALSE)
```


```{r scatters rewards by pct activated with learner}
# Total rewards
plot_cluster_withSac_Learner <- function(cluster, regression_results, cluster_ieg_Learner_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Learner_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$TotalRewards,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Rewards",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 1.5,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sac_Learner == "Learner", "#547B80", "#D1D3D4")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(TotalRewards ~ Activated_Percentage + Sac_Learner, data = cluster_df)
  fit_learner <- lm(TotalRewards ~ Activated_Percentage, data = cluster_df[cluster_df$Sac_Learner == "Learner", ])
  fit_nonlearner <- lm(TotalRewards ~ Activated_Percentage, data = cluster_df[cluster_df$Sac_Learner == "Non_Learner", ])
  
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
clusters <- unique(regression_results$Cluster)

# Iterate through clusters
for (cluster in clusters) {
  # AD SErpinb7 Glut, SCsg Gabrr2 Gaba, LH-MH Glut, ADP-MPO Trp73 Glut, Hist Gaba, SNc-VTA-RAmb Foxa1 Dopa, TH Prkcd Grin2c Glut_1 not included because at least one of the groups had all 0 activated percentage
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/RewardsbyActivatedPercentage_only3IEGs_250725/cluster_", cluster, "_plot_bySac_Learner.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_withSac_Learner(cluster, regression_results, cluster_ieg_Learner_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'RewardsbyActivatedPercentage' directory.\n")

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
  # AD Serpinb7 Glut, ADP-MPO Trp73 Glut, Hist Gaba, LH-MH Glut.SCsg Gabrr2 Gaba, SNc-VTA-RAmb Foxa1 Dopa, TH Prkcd Grin2c Glut_1 not plotted
  pdf(file = paste0("figures/IEGactivity/TimeinZoneActivatedPercentage_only3IEGs_250812/cluster_", cluster, "_plot_bySac_Learner.pdf"), width = 7, height = 8)
  
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
  pdf(file = paste0("figures/IEGactivity/EntriesinZoneActivatedPercentage_only3IEGs_250812/cluster_", cluster, "_plot_bySac_Learner.pdf"), width = 7, height = 8)
  
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
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Total_Cells * 100) %>%
  ungroup()

# Create a plotting function for the control
## Total rewards
plot_control <- function(sample_data) {
  # Fit regression line
  fit <- lm(TotalRewards ~ Activated_Percentage, data = sample_data)
  fit_summary <- summary(fit)
  
  # Create the plot
  plot(sample_data$Activated_Percentage, sample_data$TotalRewards,
       main = "All Neurons",
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Rewards",
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
  file = "figures/IEGactivity/RewardsbyActivatedPercentage_only3IEGs_250725/control_all_neurons_plot.pdf",
  width = 7, height = 8
)
plot_control(overall_sample_percentages)
dev.off()

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
  file = "figures/IEGactivity/TimeinZoneActivatedPercentage_only3IEGs_250812/control_all_neurons_plot.pdf",
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
  file = "figures/IEGactivity/EntriesinZoneActivatedPercentage_only3IEGs_250812/control_all_neurons_plot.pdf",
  width = 7, height = 8
)
plot_control_entry(overall_sample_percentages)
dev.off()
```

```{r regression rewards by pct activity in all neurons}
# Total rewards
overall_stats <- overall_sample_percentages %>%
  {
    fit <- lm(TotalRewards ~ Activated_Percentage, data = .)
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
write.csv(overall_stats,
          file = "results/IEGactivity/250725_only3IEGs_Rewards_byactivationpercent_allneurons.csv",
          row.names = FALSE
)

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
          file = "results/IEGactivity/250812_timeinzone_byactivationpercent_allneurons_only3IEGs.csv",
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
          file = "results/IEGactivity/250812_entriesinzone_byactivationpercent_allneurons_only3IEGs.csv",
          row.names = FALSE
)
```

```{r scatter rewards by pct activity in all neurons by sex}
# Calculate overall activation by sample and sex
overall_sample_sex_percentages <- cluster_ieg_sample_df %>%
  group_by(Sample, Sex) %>%
  summarise(
    Total_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    TotalRewards = dplyr::first(TotalRewards),
    TimeinZone = dplyr::first(TimeinZone),
    EntriesinZone = dplyr::first(EntriesinZone)
    ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Total_Cells * 100) %>%
  ungroup()

# Create a plotting function for the control with sex
## Total rewards
plot_control_withsex <- function(sample_data) {
  # Fit lines for all data, male, and female
  fit_all <- lm(TotalRewards ~ Activated_Percentage + Sex, data = sample_data)
  fit_male <- lm(TotalRewards ~ Activated_Percentage, data = sample_data[sample_data$Sex == "M", ])
  fit_female <- lm(TotalRewards ~ Activated_Percentage, data = sample_data[sample_data$Sex == "F", ])
  
  # Create the plot
  plot(sample_data$Activated_Percentage, sample_data$TotalRewards,
       main = "All Neurons",
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Rewards",
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
  file = "figures/IEGactivity/RewardsbyActivatedPercentage_only3IEGs_250725/control_all_neurons_plot_bySex.pdf",
  width = 7, height = 8
)
plot_control_withsex(overall_sample_sex_percentages)
dev.off()

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
  file = "figures/IEGactivity/TimeinZoneActivatedPercentage_only3IEGs_250812/control_all_neurons_plot_bySex.pdf",
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
  file = "figures/IEGactivity/EntriesinZoneActivatedPercentage_only3IEGs_250812/control_all_neurons_plot_bySex.pdf",
  width = 7, height = 8
)
plot_control_withsex_entry(overall_sample_sex_percentages)
dev.off()
```

```{r regression rewards by pct activity in all neurons by sex}
# Total rewards
overall_sex_stats <- overall_sample_sex_percentages %>%
  {
    # Full model with sex
    fit <- lm(TotalRewards ~ Activated_Percentage + Sex, data = .)
    sum_fit <- summary(fit)
    
    # Sex-specific models
    fit_male <- lm(TotalRewards ~ Activated_Percentage, data = .[.$Sex == "M", ])
    fit_female <- lm(TotalRewards ~ Activated_Percentage, data = .[.$Sex == "F", ])
    
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
write.csv(overall_sex_stats,
          file = "results/IEGactivity/250725_only3IEGs_Rewards_byactivationpercent_allneurons_bysex.csv",
          row.names = FALSE
)

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
          file = "results/IEGactivity/250812_timeinzone_byactivationpercent_allneurons_bysex_only3IEGs.csv",
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
          file = "results/IEGactivity/250812_entriesinzone_byactivationpercent_allneurons_bysex_only3IEGs.csv",
          row.names = FALSE
)
```


## Plotting B values
```{r plotting beta coefficients rewards}
# First, go back to the nested data and models to get standard errors
standard_errors <- cluster_ieg_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(TotalRewards ~ Activated_Percentage, data = .x)),
    std_error = map_dbl(lm_model, ~ summary(.x)$coefficients[2, 2])
  ) %>%
  dplyr::select(Cluster, std_error)

# Join the standard errors with your existing results
cluster_regression_results <- cluster_regression_results %>%
  left_join(standard_errors, by = "Cluster")

# Now create the plot with error bars
p_beta <- ggplot(
  cluster_regression_results,
  aes(x = Slope, y = reorder(Cluster, Slope))
) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    fill = "#888888"
  ) +
  geom_errorbar(
    aes(
      xmin = Slope - std_error,
      xmax = Slope + std_error
    ),
    position = position_dodge(0.7),
    width = 0.2,
    size = 0.2
  ) +
  labs(
    x = "Beta Coefficient (Slope)",
    y = "Cluster"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color = rep(c("#000000", "#555555"),
                                           length.out = length(unique(cluster_regression_results$Cluster))
    ))
  )

p_beta

ggsave("figures/IEGactivity/intermed-figs/250725_only3IEGs_rewardsbyactivity_betacoeffs.pdf",
       device = pdf,
       width = 6.5,
       height = 8
)
```

```{r selecting sig clusters rewards}
sig_clusters_rewards <- cluster_regression_results[cluster_regression_results$p_value < 0.10, ]

sig_clusters_rewards_names <- sort(as.character(sig_clusters_rewards$Cluster))
print(sig_clusters_rewards)
```

```{r plotting beta coefficients rewards sig only}
sig_clusters_rewards$Cluster <- factor(sig_clusters_rewards$Cluster,
                                       levels = rev(sort(unique(as.character(sig_clusters_rewards$Cluster))))
)

# Now create the plot
p_beta <- ggplot(
  sig_clusters_rewards,
  aes(x = Slope, y = Cluster)
) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    fill = "#888888"
  ) +
  geom_errorbar(
    aes(
      xmin = Slope - std_error,
      xmax = Slope + std_error
    ),
    position = position_dodge(0.7),
    width = 0.2,
    size = 0.2
  ) +
  labs(
    x = "Beta Coefficient (Slope)",
    y = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color = rep(c("#000000", "#555555"),
                                           length.out = length(unique(sig_clusters_rewards))
    ))
  )

p_beta

ggsave("figures/IEGactivity/250725_only3IEGs_rewardsbyactivity_betacoeffs_sigonly.pdf",
       device = pdf,
       width = 6,
       height = 8
)
```

## Sac Rewards

```{r regression sac rewards by pct activated}
# Fetch expression data for all IEGs
# ieg_expression <- FetchData(combined_neuronal, vars = filtered_iegs)

# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  SacRewards = combined_neuronal$Sac_Total_Rewards,
  ieg_expression
)

# Calculate the percentage of cells expressing at least 1 IEG within each cluster, Sex, and sample
cluster_ieg_sample_percentages <- cluster_ieg_sample_df %>%
  group_by(Cluster, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    SacRewards = dplyr::first(SacRewards)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

cluster_regression_results <- cluster_ieg_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(SacRewards ~ Activated_Percentage, data = .x)),
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
adjusted_p_values <- p.adjust(cluster_regression_results$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_regression_results$AdjustedPValue <- adjusted_p_values
cluster_regression_results <- cluster_regression_results %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_regression_results)

# # # Write the results to a CSV file
write.csv(cluster_regression_results, file = "results/IEGactivity/250728_only3IEGs_sac_rewards_byactivationpercent.csv", row.names = FALSE)
```

```{r scatters sac rewards by pct activated}
plot_cluster <- function(cluster, cluster_regression_results, cluster_ieg_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- cluster_regression_results[cluster_regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$SacRewards,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Sac Rewards",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l"
  )
  
  # Fit line
  fit <- lm(SacRewards ~ Activated_Percentage, data = cluster_df)
  
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

# Get all unique cluster numbers
clusters <- unique(cluster_regression_results$Cluster)


# Get all unique cluster numbers
clusters <- unique(cluster_regression_results$Cluster)

# Create a directory to save the PDFs if it doesn't exist
dir.create("figures/IEGactivity/SacRewardsbyActivatedPercentage_only3IEGs_250728", showWarnings = FALSE)

# Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/SacRewardsbyActivatedPercentage_only3IEGs_250728/cluster_", cluster, "_plot.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster(cluster, cluster_regression_results, cluster_ieg_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'SacRewardsbyActivatedPercentage' directory.\n")
```

### with Sex
```{r regression sac rewards by pct activated and sex }
# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_Sex_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  SacRewards = combined_neuronal$Sac_Total_Rewards,
  ieg_expression
)
# Calculate the percentage of cells expressing at least 1 IEG within each cluster, Sex, and sample
cluster_ieg_Sex_sample_percentages <- cluster_ieg_Sex_sample_df %>%
  group_by(Cluster, Sex, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    SacRewards = dplyr::first(SacRewards)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Perform linear regression for each cluster
regression_results <- cluster_ieg_Sex_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(SacRewards ~ Activated_Percentage + Sex, data = .x)),
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
adjusted_p_values_activated <- p.adjust(regression_results$p_value_activated_percentage, method = "BH")
adjusted_p_values_sex <- p.adjust(regression_results$p_value_sex, method = "BH")

# Add adjusted p-values to the results data frame
regression_results$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results$adjusted_p_values_sex <- adjusted_p_values_sex

regression_results <- regression_results %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results)

# # Write the results to a CSV file
write.csv(regression_results, file = "results/IEGactivity/250728_sacrewards_byactivationpercent_bysex.csv", row.names = FALSE)
```

```{r scatters sac rewards by pct activated with sex}
# First, let's create a function for our plotting code
plot_cluster_withsex <- function(cluster, regression_results, cluster_ieg_Sex_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Sex_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$SacRewards,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Sac Rewards",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sex == "M", "#4478AB", "#ED6677")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(SacRewards ~ Activated_Percentage + Sex, data = cluster_df)
  fit_male <- lm(SacRewards ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "M", ])
  fit_female <- lm(SacRewards ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "F", ])
  
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
  legend("topleft",
         legend = legend_labels,
         col = c("black", "#4478AB", "#ED6677", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 2)),
         cex = 1.2,
         bty = "n"
  )
}
# Get all unique cluster numbers
clusters <- unique(regression_results$Cluster)

# Iterate through clusters
for (cluster in clusters[62:65]) {
  # SCsg Gabrr2 Gaba, ADP-MPO Trp73 Glut, SNc-VTA-RAmb Foxa1 Dopa not plotted because at least one of the groups had 0 activated percentage
  # Open a PDF device
  pdf(file = paste0("figures/IEGactivity/SacRewardsbyActivatedPercentage_only3IEGs_250728/cluster_", cluster, "_plot_bySex.pdf"), width = 7, height = 8)
  
  # Create the plot
  plot_cluster_withsex(cluster, regression_results, cluster_ieg_Sex_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'SacRewardsbyActivatedPercentage' directory.\n")
```

### Control: by Pct Activity in All Clusters

```{r scatter sac rewards by pct activity in all neurons}
# First calculate the activation percentage across all neurons for each sample
overall_sample_percentages <- cluster_ieg_sample_df %>%
  group_by(Sample) %>%
  summarise(
    Total_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    SacRewards = dplyr::first(SacRewards)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Total_Cells * 100) %>%
  ungroup()

# Create a plotting function for the control
plot_control <- function(sample_data) {
  # Fit regression line
  fit <- lm(SacRewards ~ Activated_Percentage, data = sample_data)
  fit_summary <- summary(fit)
  
  # Create the plot
  plot(sample_data$Activated_Percentage, sample_data$SacRewards,
       main = "All Neurons",
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Sac Rewards",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(sample_data$Activated_Percentage)),
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
  legend("topright",
         legend = legend_labels,
         col = "black",
         cex = 1.2,
         bty = "n"
  )
}

# Save the  plot
pdf(
  file = "figures/IEGactivity/SacRewardsbyActivatedPercentage_only3IEGs_250728/control_all_neurons_plot.pdf",
  width = 7, height = 8
)
plot_control(overall_sample_percentages)
dev.off()
```

```{r regression sac rewards by pct activity in all neurons}
overall_stats <- overall_sample_percentages %>%
  {
    fit <- lm(SacRewards ~ Activated_Percentage, data = .)
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
write.csv(overall_stats,
          file = "results/IEGactivity/250728_only3IEGs_sacrewards_byactivationpercent_allneurons.csv",
          row.names = FALSE
)
```

```{r scatter sac rewards by pct activity in all neurons by sex}
# Calculate overall activation by sample and sex
overall_sample_sex_percentages <- cluster_ieg_sample_df %>%
  group_by(Sample, Sex) %>%
  summarise(
    Total_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(iegs)) > 0) >= 1),
    SacRewards = dplyr::first(SacRewards)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Total_Cells * 100) %>%
  ungroup()

# Create a plotting function for the control with sex
plot_control_withsex <- function(sample_data) {
  # Fit lines for all data, male, and female
  fit_all <- lm(SacRewards ~ Activated_Percentage + Sex, data = sample_data)
  fit_male <- lm(SacRewards ~ Activated_Percentage, data = sample_data[sample_data$Sex == "M", ])
  fit_female <- lm(SacRewards ~ Activated_Percentage, data = sample_data[sample_data$Sex == "F", ])
  
  # Create the plot
  plot(sample_data$Activated_Percentage, sample_data$SacRewards,
       main = "All Neurons",
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Sac Rewards",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(sample_data$Activated_Percentage)),
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
  legend("topleft",
         legend = legend_labels,
         col = c("black", "#4478AB", "#ED6677", rep("black", 2)),
         lty = c(1, 1, 1, rep(NA, 5)),
         cex = 1.2,
         bty = "n"
  )
}

# Save the control plot with sex
pdf(
  file = "figures/IEGactivity/SacRewardsbyActivatedPercentage_only3IEGs_250728/control_all_neurons_plot_bySex.pdf",
  width = 7, height = 8
)
plot_control_withsex(overall_sample_sex_percentages)
dev.off()
```

```{r get overall sex stats}
overall_sex_stats <- overall_sample_sex_percentages %>%
  {
    # Full model with sex
    fit <- lm(SacRewards ~ Activated_Percentage + Sex, data = .)
    sum_fit <- summary(fit)
    
    # Sex-specific models
    fit_male <- lm(SacRewards ~ Activated_Percentage, data = .[.$Sex == "M", ])
    fit_female <- lm(SacRewards ~ Activated_Percentage, data = .[.$Sex == "F", ])
    
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
write.csv(overall_sex_stats,
          file = "results/IEGactivity/250728_only3IEGs_sacrewards_byactivationpercent_allneurons_bysex.csv",
          row.names = FALSE
)
```

## Distance
```{r regressions distance by pct activated}
# Fetch expression data for all IEGs
# ieg_expression <- FetchData(combined_neuronal, vars = filtered_iegs)

# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  Distance = combined_neuronal$Avg_Test_Distance,
  ieg_expression
)

# Calculate the percentage of cells expressing 2 or more IEGs within each cluster, Sex, and sample
cluster_ieg_sample_percentages <- cluster_ieg_sample_df %>%
  group_by(Cluster, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    Distance = dplyr::first(Distance)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Removing outlier for distance
cluster_ieg_sample_percentages <- cluster_ieg_sample_percentages %>%
  filter(Sample != "Het_F_6")

# Perform linear regression for each cluster
cluster_regression_results <- cluster_ieg_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(Distance ~ Activated_Percentage, data = .x)),
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
adjusted_p_values <- p.adjust(cluster_regression_results$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_regression_results$AdjustedPValue <- adjusted_p_values
cluster_regression_results <- cluster_regression_results %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_regression_results)

# # Write the results to a CSV file
write.csv(cluster_regression_results, file = "../../results/neuronal_analyses/IEGactivity/250303_Distance_byactivationpercent.csv", row.names = FALSE)
```

```{r scatters distance by pct activated}
# First, let's create a function for our plotting code
plot_cluster <- function(cluster, cluster_regression_results, cluster_ieg_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- cluster_regression_results[cluster_regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$Distance,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Total Distance",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l"
  )
  
  # Fit line
  fit <- lm(Distance ~ Activated_Percentage, data = cluster_df)
  
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

# Get all unique cluster numbers
clusters <- unique(cluster_regression_results$Cluster)

# Create a directory to save the PDFs if it doesn't exist
dir.create("../../results/neuronal_analyses/IEGactivity/DistancebyActivatedPercentage_250303", showWarnings = FALSE)

# Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("../../results/neuronal_analyses/IEGactivity/DistancebyActivatedPercentage_250303/cluster_", cluster, "_plot.pdf"), width = 6, height = 8)
  
  # Create the plot
  plot_cluster(cluster, cluster_regression_results, cluster_ieg_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'DistancebyActivatedPercentage' directory.\n")
```

### with Sex
```{r regressions distance by pct activated and sex }
# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_Sex_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  Distance = combined_neuronal$Avg_Test_Distance,
  ieg_expression
)
# Calculate the percentage of cells expressing 2 or more IEGs within each cluster, Sex, and sample
cluster_ieg_Sex_sample_percentages <- cluster_ieg_Sex_sample_df %>%
  group_by(Cluster, Sex, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    Distance = dplyr::first(Distance)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Removing outlier for distance
cluster_ieg_Sex_sample_percentages <- cluster_ieg_Sex_sample_percentages %>%
  filter(Sample != "Het_F_6")

# Perform linear regression for each cluster
regression_results <- cluster_ieg_Sex_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(Distance ~ Activated_Percentage + Sex, data = .x)),
    summary = map(lm_model, summary),
    coefficients = map(lm_model, coef),
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
    p_value_activated_percentage, p_value_sex
  )

# Adjust p-values for multiple comparisons (optional)
adjusted_p_values_activated <- p.adjust(regression_results$p_value_activated_percentage, method = "BH")
adjusted_p_values_sex <- p.adjust(regression_results$p_value_sex, method = "BH")

# Add adjusted p-values to the results data frame
regression_results$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results$adjusted_p_values_sex <- adjusted_p_values_sex

regression_results <- regression_results %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results)

# Write the results to a CSV file
write.csv(regression_results, file = "../../results/neuronal_analyses/IEGactivity/250303_distance_byactivationpercent_bysex.csv", row.names = FALSE)
```

```{r scatters distance by pct activated with sex}
# First, let's create a function for our plotting code
plot_cluster_withsex <- function(cluster, regression_results, cluster_ieg_Sex_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Sex_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$Distance,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Mean Distance",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sex == "M", "#4478AB", "#ED6677")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(Distance ~ Activated_Percentage + Sex, data = cluster_df)
  fit_male <- lm(Distance ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "M", ])
  fit_female <- lm(Distance ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "F", ])
  
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
    paste("Activated % Slope =", round(cluster_results$Slope_Activated_Percentage, 3)),
    paste("Sex Slope =", round(cluster_results$Slope_Sex, 3)),
    paste("Intercept =", round(cluster_results$Intercept, 3)),
    paste("Activated % p-value =", format.pval(cluster_results$p_value_activated_percentage, digits = 3)),
    paste("Sex p-value =", format.pval(cluster_results$p_value_sex, digits = 3)),
    paste("Activated % Adj. p-value =", format.pval(cluster_results$adjusted_p_values_activated, digits = 3)),
    paste("Sex Adj. p-value =", format.pval(cluster_results$adjusted_p_values_sex, digits = 3))
  )
  
  # Add legend
  legend("topleft",
         legend = legend_labels,
         col = c("black", "#4478AB", "#ED6677", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 5)),
         cex = 1.2,
         bty = "n"
  )
}
# Get all unique cluster numbers
clusters <- unique(regression_results$Cluster)

# Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("../../results/neuronal_analyses/IEGactivity/DistancebyActivatedPercentage_250303/cluster_", cluster, "_plot_bySex.pdf"), width = 6, height = 8)
  
  # Create the plot
  plot_cluster_withsex(cluster, regression_results, cluster_ieg_Sex_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'DistancebyActivatedPercentage' directory.\n")
```

## Sac Distance

```{r regressions sac distance by pct activated}
# Fetch expression data for all IEGs
# ieg_expression <- FetchData(combined_neuronal, vars = filtered_iegs)

# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  SacDistance = combined_neuronal$Sac_Test_Distance,
  ieg_expression
)

# Calculate the percentage of cells expressing 2 or more IEGs within each cluster, Sex, and sample
cluster_ieg_sample_percentages <- cluster_ieg_sample_df %>%
  group_by(Cluster, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    SacDistance = dplyr::first(SacDistance)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Removing outlier for distance
cluster_ieg_sample_percentages <- cluster_ieg_sample_percentages %>%
  filter(Sample != "Het_F_6")

cluster_regression_results <- cluster_ieg_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(SacDistance ~ Activated_Percentage, data = .x)),
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
adjusted_p_values <- p.adjust(cluster_regression_results$p_value, method = "BH")

# Add adjusted p-values to the results data frame
cluster_regression_results$AdjustedPValue <- adjusted_p_values
cluster_regression_results <- cluster_regression_results %>%
  arrange(AdjustedPValue)

# Print the updated results
print(cluster_regression_results)

# # # Write the results to a CSV file
write.csv(cluster_regression_results, file = "../../results/neuronal_analyses/IEGactivity/250303_sac_distance_byactivationpercent.csv", row.names = FALSE)
```


```{r scatters sac distance by pct activated}
# First, let's create a function for our plotting code
plot_cluster <- function(cluster, cluster_regression_results, cluster_ieg_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- cluster_regression_results[cluster_regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$SacDistance,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Sac Distance",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l"
  )
  
  # Fit line
  fit <- lm(SacDistance ~ Activated_Percentage, data = cluster_df)
  
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

# Get all unique cluster numbers
clusters <- unique(cluster_regression_results$Cluster)

# Create a directory to save the PDFs if it doesn't exist
dir.create("../../results/neuronal_analyses/IEGactivity/Sac_DistancebyActivatedPercentage_250303", showWarnings = FALSE)

# Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("../../results/neuronal_analyses/IEGactivity/Sac_DistancebyActivatedPercentage_250303/cluster_", cluster, "_plot.pdf"), width = 6, height = 8)
  
  # Create the plot
  plot_cluster(cluster, cluster_regression_results, cluster_ieg_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'Sac_DistancebyActivatedPercentage' directory.\n")
```

### with Sex
```{r regression sac distance by pct activated and sex }
# Create a data frame with cluster assignments, Sex, sample, and IEG expression
cluster_ieg_Sex_sample_df <- data.frame(
  Cluster = Idents(combined_neuronal),
  Sex = combined_neuronal$sex,
  Sample = combined_neuronal$sample,
  SacDistance = combined_neuronal$Sac_Test_Distance,
  ieg_expression
)
# Calculate the percentage of cells expressing 2 or more IEGs within each cluster, Sex, and sample
cluster_ieg_Sex_sample_percentages <- cluster_ieg_Sex_sample_df %>%
  group_by(Cluster, Sex, Sample) %>%
  summarise(
    Cluster_Cells = n(),
    Activated_Cells = sum(rowSums(across(all_of(filtered_iegs)) > 0) >= 2),
    SacDistance = dplyr::first(SacDistance)
  ) %>%
  dplyr::mutate(Activated_Percentage = Activated_Cells / Cluster_Cells * 100) %>%
  ungroup()

# Removing outlier for distance
cluster_ieg_Sex_sample_percentages <- cluster_ieg_Sex_sample_percentages %>%
  filter(Sample != "Het_F_6")

# Perform linear regression for each cluster
regression_results <- cluster_ieg_Sex_sample_percentages %>%
  group_by(Cluster) %>%
  nest() %>%
  dplyr::mutate(
    lm_model = map(data, ~ lm(SacDistance ~ Activated_Percentage + Sex, data = .x)),
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
adjusted_p_values_activated <- p.adjust(regression_results$p_value_activated_percentage, method = "BH")
adjusted_p_values_sex <- p.adjust(regression_results$p_value_sex, method = "BH")

# Add adjusted p-values to the results data frame
regression_results$adjusted_p_values_activated <- adjusted_p_values_activated
regression_results$adjusted_p_values_sex <- adjusted_p_values_sex

regression_results <- regression_results %>%
  arrange(adjusted_p_values_activated)

# Print the updated results
print(regression_results)

# # Write the results to a CSV file
write.csv(regression_results, file = "../../results/neuronal_analyses/IEGactivity/250303_SacDistance_byactivationpercent_bysex.csv", row.names = FALSE)
```

```{r scatters sac distance by pct activated with sex}
# First, let's create a function for our plotting code
plot_cluster_withsex <- function(cluster, regression_results, cluster_ieg_Sex_sample_percentages) {
  # Extract the results for the specific cluster we're plotting
  cluster_results <- regression_results[regression_results$Cluster == cluster, ]
  
  # Filter the data for the specific cluster
  cluster_df <- cluster_ieg_Sex_sample_percentages %>%
    filter(Cluster == cluster)
  
  # Create the plot
  plot(cluster_df$Activated_Percentage, cluster_df$SacDistance,
       main = paste("Cluster", cluster),
       pch = 16, cex = 2,
       xlab = "Activated Percentage", ylab = "Sac Distance",
       cex.lab = 1.5, cex.axis = 1.5,
       cex.main = 2,
       xlim = c(0, max(cluster_df$Activated_Percentage)),
       bty = "l",
       col = ifelse(cluster_df$Sex == "M", "#4478AB", "#ED6677")
  )
  
  # Fit lines for all data, male, and female
  fit_all <- lm(SacDistance ~ Activated_Percentage + Sex, data = cluster_df)
  fit_male <- lm(SacDistance ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "M", ])
  fit_female <- lm(SacDistance ~ Activated_Percentage, data = cluster_df[cluster_df$Sex == "F", ])
  
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
    paste("Activated % Slope =", round(cluster_results$Slope_Activated_Percentage, 3)),
    paste("Sex Slope =", round(cluster_results$Slope_Sex, 3)),
    paste("Intercept =", round(cluster_results$Intercept, 3)),
    paste("Activated % p-value =", format.pval(cluster_results$p_value_activated_percentage, digits = 3)),
    paste("Sex p-value =", format.pval(cluster_results$p_value_sex, digits = 3)),
    paste("Activated % Adj. p-value =", format.pval(cluster_results$adjusted_p_values_activated, digits = 3)),
    paste("Sex Adj. p-value =", format.pval(cluster_results$adjusted_p_values_sex, digits = 3))
  )
  
  # Add legend
  legend("topleft",
         legend = legend_labels,
         col = c("black", "#4478AB", "#ED6677", rep("black", 5)),
         lty = c(1, 1, 1, rep(NA, 5)),
         cex = 1.2,
         bty = "n"
  )
}

# Get all unique cluster numbers
clusters <- unique(regression_results$Cluster)

# Iterate through clusters
for (cluster in clusters) {
  # Open a PDF device
  pdf(file = paste0("../../results/neuronal_analyses/IEGactivity/Sac_DistancebyActivatedPercentage_250303/cluster_", cluster, "_plot_bySex.pdf"), width = 6, height = 8)
  
  # Create the plot
  plot_cluster_withsex(cluster, regression_results, cluster_ieg_Sex_sample_percentages)
  
  # Close the PDF device
  dev.off()
  
  # Print progress
  cat("Saved plot for cluster", cluster, "\n")
}

cat("All plots have been saved in the 'SacDistancebyActivatedPercentage' directory.\n")
```

# Gene Module Analysis

```{r add module scores}
ieg_lst <- list(ieg_62 = filtered_iegs)

combined_neuronal <- AddModuleScore(
  object = combined_neuronal,
  features = ieg_lst,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "IEG_score",
  seed = 1,
  search = FALSE,
  slot = "data"
)
```

```{r set up dataframe}
sample_info <- combined_neuronal@meta.data %>%
  dplyr::select(
    sample,
    rewards = Avg_Total_Rewards,
    sac_rewards = Sac_Total_Rewards,
    distance = Avg_Test_Distance,
    sac_distance = Sac_Test_Distance,
    cell_type,
    genotype,
    sex,
    score_ieg62 = IEG_score1
  )
```

## Rewards

```{r linear regression for IEG module rewards}
model_summary_df <- data.frame()

for (cluster in unique(sample_info$cell_type)) {
  # Filter the data frame for the current cluster
  cluster_df <- sample_info %>%
    filter(cell_type == cluster) %>%
    distinct()
  
  # Fit the linear regression model
  lm_model <- lm(rewards ~ score_ieg62, data = cluster_df)
  
  # Extract the model summary
  model_summary <- summary(lm_model)
  
  # Create a named vector with model summary statistics
  model_summary_vec <- c(
    `(Intercept)` = coef(lm_model)[1],
    `(Intercept) p-value` = coef(model_summary)[1, "Pr(>|t|)"],
    coef(lm_model)[2],
    `Std. Error` = coef(model_summary)[2, "Std. Error"],
    `t_value` = coef(model_summary)[2, "t value"],
    `p_value` = coef(model_summary)[2, "Pr(>|t|)"],
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared,
    f_statistic = model_summary$fstatistic[1],
    f_pvalue = pf(model_summary$fstatistic[1], model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = FALSE)
  )
  
  # Convert the named vector to a data frame
  model_summary_vec_df <- data.frame(t(model_summary_vec))
  
  # Add a column for the cluster
  model_summary_vec_df$Cluster <- cluster
  
  # Append the data frame to the overall result
  model_summary_df <- rbind(model_summary_df, model_summary_vec_df)
}

model_summary_df$adj_p_value <- p.adjust(model_summary_df$p_value, method = "BH")

model_summary_df %>%
  write_csv("../../results/neuronal_analyses/IEGactivity/250303_genemodulescore_results_rewards.csv")
```

```{r plots for rewards}
# filter only significant results
sig_df <- model_summary_df %>% filter(p_value < 0.05)

cts <- sig_df$Cluster %>% unique()

dir.create("../../results/neuronal_analyses/IEGactivity/RewardsbyModuleScore_Figs", showWarnings = FALSE)

for (ct in cts) {
  plot_df <- sample_info %>%
    filter(cell_type == ct)
  
  pp <- ggplot(plot_df, aes(score_ieg62, rewards)) +
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, color = "black", linewidth = 0.75) +
    geom_smooth(method = "lm", aes(color = sex), se = FALSE, fullrange = TRUE, linewidth = 0.75) +
    scale_color_manual(values = palette_sex) +
    scale_x_continuous(expand = expansion(c(0, 0))) +
    labs(title = ct, color = "Sex", x = "Score_IEG62", y = "Rewards") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 16),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14)
    )
  
  print(pp)
  ct_name <- gsub("/", "-", ct)
  ggsave(paste0("../../results/neuronal_analyses/IEGactivity/RewardsbyModuleScore_Figs/250303_linearplot_", ct_name, "_genemodulescore.pdf"), pp, width = 5, height = 3.5)
}
```

## Sac rewards

```{r linear regression for IEG module sac_rewards}
model_summary_df <- data.frame()

for (cluster in unique(sample_info$cell_type)) {
  # Filter the data frame for the current cluster
  cluster_df <- sample_info %>%
    filter(cell_type == cluster) %>%
    distinct()
  
  # Fit the linear regression model
  lm_model <- lm(sac_rewards ~ score_ieg62, data = cluster_df)
  
  # Extract the model summary
  model_summary <- summary(lm_model)
  
  # Create a named vector with model summary statistics
  model_summary_vec <- c(
    `(Intercept)` = coef(lm_model)[1],
    `(Intercept) p-value` = coef(model_summary)[1, "Pr(>|t|)"],
    coef(lm_model)[2],
    `Std. Error` = coef(model_summary)[2, "Std. Error"],
    `t_value` = coef(model_summary)[2, "t value"],
    `p_value` = coef(model_summary)[2, "Pr(>|t|)"],
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared,
    f_statistic = model_summary$fstatistic[1],
    f_pvalue = pf(model_summary$fstatistic[1], model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = FALSE)
  )
  
  # Convert the named vector to a data frame
  model_summary_vec_df <- data.frame(t(model_summary_vec))
  
  # Add a column for the cluster
  model_summary_vec_df$Cluster <- cluster
  
  # Append the data frame to the overall result
  model_summary_df <- rbind(model_summary_df, model_summary_vec_df)
}

model_summary_df$adj_p_value <- p.adjust(model_summary_df$p_value, method = "BH")

model_summary_df %>%
  write_csv("../../results/neuronal_analyses/IEGactivity/250303_genemodulescore_results_sac_rewards.csv")
```

```{r plots for sac_rewards}
# filter only significant results
sig_df <- model_summary_df %>% filter(p_value < 0.05)

cts <- sig_df$Cluster %>% unique()

dir.create("../../results/neuronal_analyses/IEGactivity/sac_rewardsbyModuleScore_Figs", showWarnings = FALSE)

for (ct in cts) {
  plot_df <- sample_info %>%
    filter(cell_type == ct)
  
  pp <- ggplot(plot_df, aes(score_ieg62, sac_rewards)) +
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, color = "black", linewidth = 0.75) +
    geom_smooth(method = "lm", aes(color = sex), se = FALSE, fullrange = TRUE, linewidth = 0.75) +
    scale_color_manual(values = palette_sex) +
    scale_x_continuous(expand = expansion(c(0, 0))) +
    labs(title = ct, color = "Sex", x = "Score_IEG62", y = "sac_rewards") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 16),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14)
    )
  
  print(pp)
  ct_name <- gsub("/", "-", ct)
  ggsave(paste0("../../results/neuronal_analyses/IEGactivity/sac_rewardsbyModuleScore_Figs/250303_linearplot_", ct_name, "_genemodulescore.pdf"), pp, width = 5, height = 3.5)
}
```
