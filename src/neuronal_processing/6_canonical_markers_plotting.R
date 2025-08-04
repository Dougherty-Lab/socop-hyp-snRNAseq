---
title: "Plotting of cananoical marker genes in tuberal hypothalamus of MYT1L mice after social operant conditioning"
author: "Jungeun Ji"
date: "2025-07-28"
output:
  ---
  
# Set up environment
```{r message=FALSE, warning=TRUE}
# Change this to your directory
.libPaths("/ref/jdlab/software/r-envs/simona/r-packages/4.2/")
setwd("/scratch/jdlab/jungeun/")

library(Seurat)

library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
library(circlize)
library(viridis)

set.seed(420)

source("color_utils.R")
```

# Load in data
```{r load in data}
combined_neuronal <- readRDS("data/250228_combined.neuronal_1-5_5pct.rds")
combined_neuronal

# Saving scatterplots won't work if you have "/" in cluster name
combined_neuronal@meta.data$cell_type <- as.factor(gsub("/", "-", combined_neuronal@meta.data$cell_type))

# Subset data to only hypothalamus region
only_hypo = subset(combined_neuronal, subset = anatomical_region == "Hypothalamus")
```

# Finding region markers
```{r check marker genes for each region}
# Set marker genes
regions = c("Amygdala/Hypothalamus", "Hypothalamus", "Midbrain", "Striatum/Pallidum", "Thalamus", "Zona Incerta")

Idents(combined_neuronal) = "anatomical_region"
region_markers = FindAllMarkers(combined_neuronal, verbose = TRUE, only.pos = TRUE)

write.csv(region_markers, file = "results/neuronal_processing/250729_markers_neuronal_anatomical_regions.csv")
region_markers = read.csv("results/neuronal_processing/250729_markers_neuronal_anatomical_regions.csv")

saveRDS(region_markers, file = "results/neuronal_processing/250729_markers_neuronal_anatomical_regions.rds")

region_markers = region_markers %>% mutate(cluster = factor(cluster, levels = regions)) %>% arrange(cluster)
sig_region_markers = region_markers %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.6)
region_marker_genes = sig_region_markers %>% pull(gene) %>% unique()

selected_region_markers = c("Pde1a", #Amygdala/Hypothalamus
                   "Ndst4", #"Ghr", #Hypothalamus
                   "Lmx1b", # Midbrain
                   "Bcl11b", #"Phactr1", #Striatum/Pallidum
                   "Tcf7l2", #"Lef1", #Thalamus
                   "Pou6f2" #Zona Incerta
                   )
```

# Finding cluster markers
```{r findMarkers, message=FALSE, results=FALSE}
Idents(combined_neuronal) = "cell_type"
combined_neuronal_ct_markers_all <- FindAllMarkers(combined_neuronal,
                                                          only.pos = TRUE,
                                                          min.pct = 0.25,
                                                          logfc.threshold = 0.5,
                                                          return.thresh = 0.01,
                                                   verbose = TRUE
)

write.csv(combined_neuronal_ct_markers_all, file = "results/neuronal_processing/250729_markers_neuronal_cell_types.csv")
combined_neuronal_ct_markers_all = read.csv("results/neuronal_processing/250729_markers_neuronal_cell_types.csv")
saveRDS(combined_neuronal_ct_markers_all, file = "results/neuronal_processing/250729_markers_neuronal_cell_types.rds")
```

# Plot heatmap (region markers)
```{r Heatmap, message=FALSE, results=FALSE}
# Extract expression of genes of interest 
subset_df <- FetchData(combined_neuronal, vars = c(region_marker_genes, "cell_type", "anatomical_region"), slot = 'scale.data')
avg_expr <- subset_df %>%
  group_by(cell_type, anatomical_region) %>%
  summarise(across(all_of(region_marker_genes), mean), .groups = "drop")

avg_expr = avg_expr %>% arrange(anatomical_region)

column_regions = factor(avg_expr$anatomical_region)
# avg_expr <- avg_expr %>%
#   mutate(celltype_label = paste(anatomical_region, cell_type, sep = "_"))

# Reshape: rows = genes, columns = cell types
expr_mat <- avg_expr %>%
  select(cell_type, all_of(region_marker_genes)) %>%
  pivot_longer(-cell_type, names_to = "gene", values_to = "avg_expression") %>%
  pivot_wider(names_from = cell_type, values_from = avg_expression) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Annotation for columns
ha_col <- HeatmapAnnotation(
  Region = column_regions,
  col = list(Region = c(region_muted_colors)),
  annotation_legend_param = list(
    title_gp = gpar(fontface = "plain")  
  ),
 show_annotation_name = FALSE
)

row_labels <- ifelse(rownames(expr_mat) %in% selected_region_markers, rownames(expr_mat), "")
min_value <- quantile(expr_mat, 0.0001, na.rm = TRUE)
max_value <- quantile(expr_mat, 0.9999, na.rm = TRUE)

viridis_colors <- viridis(50)

ht <- Heatmap(expr_mat,
              name = "Avg. scaled expr",
              top_annotation = ha_col,
              col = colorRamp2(
                c(min_value + 0.5, 0, max_value - 0.5),
                viridis_colors[c(1, 25, 50)]  # map low, mid, high to viridis
              ),
              cluster_columns = FALSE, 
              cluster_rows = FALSE,
              column_split = column_regions,
              column_title = NULL,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 12),
              row_labels = row_labels,
              heatmap_legend_param = list(
                title_gp = gpar(fontface = "plain")  
              )
)

pdf("figures/neuronal_processing/250731_heatmap_avg_celltype_expr_region_markers.pdf", width = 10, height = 8)
draw(ht,
     heatmap_legend_side = "right",       # put heatmap legend at bottom
     annotation_legend_side = "right",    # put annotation legend below it
     merge_legend = TRUE                  # <--- show them as separate blocks, stacked
)
dev.off()
```

# Plot Violin
```{r VlnPlot, message=FALSE, results=FALSE}
vln_list <- VlnPlot(combined_neuronal,
                    features = selected_region_markers,
                    pt.size = 0.0,
                    cols = region_muted_colors,
                    group.by = 'anatomical_region',
                    combine = FALSE)

for (i in seq_along(vln_list)) {
  vln_list[[i]] <- vln_list[[i]] + 
  theme(legend.position = "none",
        plot.title = element_text(face = "plain", size = 14),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 90, size = 12),
        axis.title.x = element_blank())
    if (i != length(vln_list)) {
      vln_list[[i]] <- vln_list[[i]] + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
}

pdf("figures/neuronal_processing/250731_vlnplot_region_markers_perregion.pdf", width = 7.5, height = 12)
wrap_plots(vln_list, ncol = 1, guides = "collect")
dev.off()
```

# Plot heatmap (only hypothalamus, cell type markers)
```{r Heatmap, message=FALSE, results=FALSE}
selected_cluster_markers = c("Slc17a6", # Glutamatergic
                             "Gad2", # GABAergic
                             #"Tbx3", # ARC
                             "Ptpn3", # PH
                             "Sim1", # PVH
                             "Nr5a1", # VMH
                             #"Slc6a3", # Dopa
                             "Chat" # Cholinergic
                             )

# Create column that groups glutamatergic/GABAergic
only_hypo@meta.data = only_hypo@meta.data %>% mutate(
  glut_gaba = cell_type %>%
    str_split(" ") %>%
    sapply(function(x) str_split(x[length(x)], "_")[[1]][1])
) %>% 
  mutate(glut_gaba = ifelse(cell_type == "Chat GABA", "Gaba", glut_gaba)) %>% 
  mutate(glut_gaba = ifelse(glut_gaba %in% c("Glut", "Gaba"), glut_gaba, "Others"))

# Extract expression of genes of interest
subset_df <- FetchData(only_hypo, vars = c(selected_cluster_markers, "cell_type", "glut_gaba"), slot = 'scale.data') %>% 
  group_by(glut_gaba) %>% arrange(cell_type, .by_group = TRUE)

subset_df$cell_type = droplevels(subset_df$cell_type)

ordered_cell_types <- subset_df %>%
  distinct(cell_type) %>%
  pull(cell_type)

## Average expression per cell type
avg_expr <- subset_df %>%
  group_by(cell_type) %>%
  summarise(across(all_of(selected_cluster_markers), mean), .groups = "drop") %>%
  arrange(factor(cell_type, levels = ordered_cell_types))

## Reshape: rows = genes, columns = cell types
expr_mat <- avg_expr %>%
  select(cell_type, all_of(selected_cluster_markers)) %>%
  pivot_longer(-cell_type, names_to = "gene", values_to = "avg_expression") %>%
  pivot_wider(names_from = cell_type, values_from = avg_expression) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Annotation for columns
## Create palette 
celltype_colors = data.frame(cell_type = levels(combined_neuronal$cell_type), color = palette_70[1:65])
hypo_colors = celltype_colors %>% filter(cell_type %in% only_hypo$cell_type) %>% arrange(factor(cell_type, levels = ordered_cell_types))

glut_gaba_colors <- c(
  Glut = "#3B9AB2",  
  Gaba = "#ed682f",
  Others = "grey"
)

# Create annotation for heatmap
column_info <- subset_df %>% distinct(cell_type, glut_gaba) %>% arrange(factor(cell_type, levels = ordered_cell_types)) %>% select(glut_gaba)

ha_col <- HeatmapAnnotation(
  Glut_GABA = column_info$glut_gaba,
  col = list(
    Glut_GABA = glut_gaba_colors
  ),
  annotation_legend_param = list(
    title_gp = gpar(fontface = "plain")
  ),
  show_annotation_name = FALSE
)

min_value <- quantile(expr_mat, 0.001, na.rm = TRUE)
max_value <- quantile(expr_mat, 0.999, na.rm = TRUE)


ht <- Heatmap(expr_mat,
              name = "Avg. scaled expr",
              top_annotation = ha_col,
              col = colorRamp2(
                c(min_value + 0.5, 0, max_value - 0.5),
                viridis_colors[c(1, 25, 50)]  # map low, mid, high to viridis
              ),
              cluster_columns = FALSE, 
              cluster_rows = FALSE,
              column_split = column_info,
              column_title = NULL,
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 8),
              row_names_gp = gpar(fontsize = 12),
              heatmap_legend_param = list(
                title_gp = gpar(fontface = "plain")  
              )
)

pdf("figures/neuronal_processing/250731_heatmap_onlyhypo_avg_celltype_expr_cluster_markers.pdf", width = 10, height = 4)
draw(ht,
     heatmap_legend_side = "right",       # put heatmap legend at bottom
     annotation_legend_side = "right",    # put annotation legend below it
     merge_legend = TRUE                  # <--- show them as separate blocks, stacked
)
dev.off()
```

# Plot Violin
```{r VlnPlot, message=FALSE, results=FALSE}
only_hypo$cell_type <- factor(only_hypo$cell_type, levels = ordered_cell_types)

vln_list <- VlnPlot(only_hypo,
                    features = selected_cluster_markers,
                    pt.size = 0.0,
                    group.by = 'cell_type',
                    combine = FALSE)

for (i in seq_along(vln_list)) {
  vln_list[[i]] <- vln_list[[i]] + 
    scale_fill_manual(values = rep("steelblue", 31)) +
    theme(legend.position = "none",
          plot.title = element_text(face = "plain", size = 14),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(angle = 90, size = 12),
          axis.title.x = element_blank())
  if (i != length(vln_list)) {
    vln_list[[i]] <- vln_list[[i]] + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
}

pdf("figures/neuronal_processing/250731_vlnplot_onlyhypo_cluster_markers_percelltype.pdf", width = 7.5, height = 12)
wrap_plots(vln_list, ncol = 1, guides = "collect")
dev.off()
```
