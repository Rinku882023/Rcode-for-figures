# First Approach: 
# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(ggrepel)
library(Cairo)
library(pheatmap)

# Load data
chr1_association <- read.delim("chr1_association")

# Arrange miRNAs by genomic position
chr1_association <- chr1_association %>%
  arrange(start) %>%
  mutate(miRNA_order = row_number())  # Assign vertical order to miRNAs

# Arrange metabolites by sub-pathway and alphabetically
chr1_association <- chr1_association %>%
  arrange(Sub.Pathway, Metabolite) %>%
  mutate(metabolite_order = row_number())  # Assign vertical order to metabolites

# Create a new column for association type (positive/negative)
chr1_association <- chr1_association %>%
  mutate(association_type_Endo2 = ifelse(coef.Endo2. >= 0, "positive", "negative"))

# Define color schemes
association_colors <- c("positive" = "#F8766D", "negative" = "#619CFF")
Sub.Pathway_colors <- hue_pal()(length(unique(chr1_association$Sub.Pathway)))

# Plot for Dataset 1
p1 <- ggplot(chr1_association) +
  geom_segment(
    aes(
      x = 1, xend = 2,
      y = miRNA_order, yend = metabolite_order,
      size = coef.Endo2.,
      color = association_type_Endo2
    ),
    alpha = 0.7
  ) +
  geom_text(
    data = chr1_association %>% distinct(miR, .keep_all = TRUE),
    aes(x = 1, y = miRNA_order, label = miR),
    hjust = 1, fontface = "bold"
  ) +
  geom_text(
    data = chr1_association %>% distinct(Metabolite, .keep_all = TRUE),
    aes(x = 2, y = metabolite_order, label = Metabolite, color = Sub.Pathway),
    hjust = 0, fontface = "bold"
  ) +
  scale_color_manual(values = c(association_colors, setNames(Sub.Pathway_colors, unique(chr1_association$Sub.Pathway)))) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(title = "(A) Metabo-Endotype 2", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none"
  ) +
  xlim(0.8, 2.5)

# Create a new column for Endo3 association type
chr1_association <- chr1_association %>%
  mutate(association_type_Endo3 = ifelse(coef.Endo3. >= 0, "positive", "negative"))

# Plot for Dataset 2
p2 <- ggplot(chr1_association) +
  geom_segment(
    aes(
      x = 1, xend = 2,
      y = miRNA_order, yend = metabolite_order,
      size = coef.Endo3.,
      color = association_type_Endo3
    ),
    alpha = 0.7
  ) +
  geom_text(
    data = chr1_association %>% distinct(miR, .keep_all = TRUE),
    aes(x = 1, y = miRNA_order, label = miR),
    hjust = 1, fontface = "bold"
  ) +
  geom_text(
    data = chr1_association %>% distinct(Metabolite, .keep_all = TRUE),
    aes(x = 2, y = metabolite_order, label = Metabolite, color = Sub.Pathway),
    hjust = 0, fontface = "bold"
  ) +
  scale_color_manual(values = c(association_colors, setNames(Sub.Pathway_colors, unique(chr1_association$Sub.Pathway)))) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(title = "(B) Metabo-Endotype 3", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none"
  ) +
  xlim(0.8, 2.5)

# Combine plots
combined_plot <- p1 + p2
print(combined_plot)

# Second Approach: Sort miRNAs by genomic position and plot
position_data <- chr1_association %>%
  distinct(miR, start) %>%
  arrange(start) %>%
  mutate(order = row_number())

data <- chr1_association %>%
  left_join(position_data, by = "miR") %>%
  arrange(order)

ggplot(data, aes(x = reorder(miR, order), y = coef.Endo2.)) +
  geom_point(aes(color = Sub.Pathway), size = 3) +
  geom_text_repel(aes(label = Metabolite), size = 3, fontface = "bold") +
  labs(
    x = "miRNAs (sorted by genomic position)",
    y = "Association Estimate",
    color = "Sub-Pathway"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )

# Third Approach: Heatmap for Endo3
Endo3 <- read.delim("Endo3.txt")
Endo3$coef_Endo3 <- as.numeric(Endo3$coef_Endo3)
matrix_Endo3 <- with(Endo3, tapply(coef_Endo3, list(miRNA, Metabolite), FUN = mean))
matrix_Endo3[is.na(matrix_Endo3)] <- 0
rownames(matrix_Endo3) <- gsub("hsa\\-", "", rownames(matrix_Endo3))
custom_colors <- colorRampPalette(c("blue", "white", "red"))(201)
breaks <- c(seq(min(matrix_Endo3), 0, length.out = 101), seq(0, max(matrix_Endo3), length.out = 101)[-1])
CairoPNG("Endo3_association.png", width = 6000, height = 6000, res = 600)
pheatmap(matrix_Endo3,
         color = custom_colors,
         breaks = breaks,
         na_col = "grey",
         fontsize_row = 3,
         fontsize_col = 3,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
dev.off()
