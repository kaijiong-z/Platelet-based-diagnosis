# =======================================================================
# Fig. 4. The discovery of BC-associated lncRNA biomarkers in the platelet
# =======================================================================

rm(list = ls())

# library package
library(tidyverse)
library(broom)
library(GGally)
library(ggplotify)
library(pheatmap)
library(cowplot)

## load data
load("plt_lncRNA_expr.rdata")
load("group_list.rdata")
load("TEPS_DEG.rdata")

deg_expr <- plt_lncRNA_tpm_expr %>% filter(rownames(.) %in% rownames(sig_DEG))
deg_expr <- log2(deg_expr + 1)

##--- fig. 4a-b. feature selection via linear models ----

##== 1.1 linear regression/ANOVA==

# linear data
lr_data <- deg_expr %>%
  t() %>%
  as.data.frame() %>%
  mutate(
    groups = group_list,
    groups = factor(groups, levels = c("Normal", "Tumor"))
  ) %>%
  relocate(groups)

# long data
lr_long_data <- lr_data %>%
  pivot_longer(-groups, names_to = "gene", values_to = "value")


# ----- ANOVA  -----
anova_res <- lr_long_data %>%
  group_by(gene) %>%
  do({
    g_data <- .

    lm_model <- glm(value ~ groups, data = g_data)

    model <- aov(value ~ groups, data = g_data)

    res <-
      model %>%
      broom::tidy() %>%
      filter(term != "Residuals") %>%
      left_join(
        confint(lm_model) %>%
          as_tibble(rownames = "term") %>%
          dplyr::rename(CIlow = 2, CIhigh = 3) %>%
          right_join(tidy(lm_model), by = "term") %>%
          mutate(term = str_replace(term, "groupsTumor", "groups")) %>%
          filter(term != "(Intercept)") %>%
          dplyr::select(term, CIlow, CIhigh, estimate, estimate_pval = p.value),
        by = "term"
      )

    res
  }) %>%
  group_by(term) %>%
  mutate(
    adj_pval = p.adjust(p.value, method = "BH"),
    estimate.adj_pval = p.adjust(estimate_pval, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(significant = adj_pval < 0.05)


anova_total_res <-
  anova_res %>%
  mutate(
    sign_label = case_when(
      adj_pval < 0.001 ~ "***",
      adj_pval < 0.01 ~ "**",
      adj_pval < 0.05 ~ "*",
      T ~ ""
    )
  )

anova_res_sign <-
  anova_total_res %>%
  filter(significant)

##----- fig4a plot -----
fig4.a <- anova_total_res %>%
  ggplot(aes(gene, statistic, fill = significant)) +
  geom_col() +
  geom_text(
    aes(label = sign_label),
    position = "stack",
    color = "white",
    size = 8,
    hjust = 0.5,
    vjust = 1
  ) +
  #coord_flip() +
  scale_fill_manual(values = c("#075f75", "#b12115")) +
  #scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 18, color = "black"),
    axis.text.x = element_text(
      size = 15,
      color = "black",
      angle = 60,
      hjust = 1
    ),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    #legend.position = "none"
    # panel.spacing = unit(1, "mm"),
    # panel.border = element_blank(),
    plot.margin = unit(c(0, 1, 0, 1), 'lines')
  ) +
  labs(
    x = "",
    y = "Variance",
    title = "Feature selection of individual linear models"
  )

fig4.a

##== 1.2 Boruta feature selection ==

set.seed(12345)
# install.packages('Boruta')
library(Boruta)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(RColorBrewer)

# broruta data
# remove genes that are not significant in ANOVA
br_data <- lr_data %>%
  dplyr::select(
    -(anova_total_res %>%
      filter(!significant) %>%
      filter(!significant) %>%
      pull(gene))
  )


# Perform Boruta search
boruta.train <- Boruta(
  groups ~ .,
  data = br_data,
  pValue = 0.05,
  mcAdj = TRUE,
  maxRuns = 100,
  doTrace = 0
)
print(boruta.train)
table(boruta.train$finalDecision)

# Get significant variables
confirmed.var <- getSelectedAttributes(boruta.train, withTentative = F) ##Tentative=T，选择含Tentative
confirmed.var

boruta.imp <- function(x) {
  imp <- reshape2::melt(x$ImpHistory, na.rm = T)[, -1]
  colnames(imp) <- c("Variable", "Importance")
  imp <- imp[is.finite(imp$Importance), ]

  variableGrp <- data.frame(
    Variable = names(x$finalDecision),
    finalDecision = x$finalDecision
  )

  showGrp <- data.frame(
    Variable = c("shadowMax", "shadowMean", "shadowMin"),
    finalDecision = c("shadowMax", "shadowMean", "shadowMin")
  )

  variableGrp <- rbind(variableGrp, showGrp)

  boruta.variable.imp <- merge(imp, variableGrp, all.x = T)

  sortedVariable <- boruta.variable.imp %>%
    group_by(Variable) %>%
    dplyr::summarise(median = median(Importance)) %>%
    arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)

  boruta.variable.imp$Variable <- factor(
    boruta.variable.imp$Variable,
    levels = sortedVariable
  )

  invisible(boruta.variable.imp)
}

boruta.variable.imp <- boruta.imp(boruta.train)

head(boruta.variable.imp)

# Variable Importance Scores
boruta.variable.imp$finalDecision <- factor(
  boruta.variable.imp$finalDecision,
  levels = c(
    "Confirmed",
    "Tentative",
    "Rejected",
    "shadowMax",
    "shadowMean",
    "shadowMin"
  )
)

##----- fig4b plot -----
fig4.b <- boruta.variable.imp %>%
  ggplot(aes(x = Variable, y = Importance, fill = finalDecision)) +
  geom_boxplot() +
  #scale_fill_brewer(palette="Set1")+
  #scale_fill_manual(values =mycolors)+
  scale_fill_manual(
    values = c(
      "Confirmed" = "#b12115",
      "Tentative" = "#F46D43",
      "shadowMax" = "#075f75",
      "shadowMean" = "#ef9c02",
      "Rejected" = "#2e7ebb",
      "shadowMin" = "#5f5f5f"
    )
  ) +
  theme_bw() + #coord_flip()+
  labs(
    x = "",
    y = "Importance",
    title = "Feature selection of Boruta algorithm"
  ) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 18, color = "black"),
    axis.text.x = element_text(
      size = 15,
      color = "black",
      angle = 60,
      hjust = 1
    ),
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    #legend.position = "none",
    plot.margin = unit(c(0, 1, 0, 1), 'lines')
  )

fig4.b


##----- fig. 4c. Expression levels of the optimal variables in TEPs ---

library(tidyverse)
library(gridExtra)
library(ggpubr)

# optimal variables
data <- lr_data %>%
  select(groups, all_of(confirmed.var)) %>%
  relocate(groups, .after = last_col())

##----- fig.4c plot -----
plot_column_boxplot <- function(data, column) {
  ggplot(data = data, aes(x = groups, y = get(column))) +
    geom_boxplot(width = 0.7) +
    geom_jitter(aes(colour = groups), size = 6, width = 0.2, alpha = 0.9) +
    scale_color_manual(values = c(Normal = "#075f75", Tumor = "#b12115")) +
    # scale_fill_manual(values = c(Normal = "#075f75", Tumor = "#b12115")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    geom_signif(
      comparisons = list(c("Tumor", "Normal")),
      map_signif_level = function(p) sprintf("p = %.2g", p),
      position = "identity",
      test = "wilcox.test",
      textsize = 6,
      tip_length = 0.02,
      margin_top = 0.1,
      vjust = 0
    ) +
    labs(title = column, x = "", y = "Expression(TPM)") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      panel.grid = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 12),
      legend.position = "none",
      plot.margin = unit(c(0, 1, 0, 1), 'lines')
    )
}


myplots_boxplot <- lapply(colnames(data)[1:8], plot_column_boxplot, data = data)


# fig4.c
fig4.c <- plot_grid(
  myplots_boxplot[[1]],
  myplots_boxplot[[2]],
  myplots_boxplot[[3]],
  myplots_boxplot[[4]],
  myplots_boxplot[[5]],
  myplots_boxplot[[6]],
  myplots_boxplot[[7]],
  myplots_boxplot[[8]],
  #labels=LETTERS[1:8],
  ncol = 4,
  label_fontface = "bold",
  label_size = 20,
  scale = 1
) +
  theme(plot.margin = unit(c(1, 1, 1, 1), 'lines'))


# ##-- Correlation between variables --

# correlation_data <- br_data %>% select(confirmed.var)

# ## Create our correlation matrix
# correlation_matrix <-
#   correlation_data %>%
#   na.omit() %>%
#   select(where(is.numeric)) %>%
#   #column_to_rownames("id") %>%
#   as.matrix() %>%
#   cor(use = "pairwise.complete.obs", method = "spearman") ##spearman

# correlation_matrix

# library(GGally)
# library(ggplotify)
# library(pheatmap)

# gene_correlation <- as.ggplot(pheatmap(
#   correlation_matrix,
#   fontsize = 12,
#   #color = inferno(10, direction = -1),
#   color = c(
#     rev(
#       RColorBrewer::brewer.pal(9, "Blues")
#     ),
#     RColorBrewer::brewer.pal(9, "YlOrRd")
#   ),
#   #color = c(colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)),
#   na_col = "white",
#   breaks = seq(-1, 1, length.out = 18),
#   border_color = NA,
#   clustering_method = "ward.D2",
#   #clustering_method = "average",
#   # cluster_rows = F,
#   #cluster_cols = F,
#   display_numbers = T,
#   number_format = "%.2f",
#   number_color = "black",
#   fontsize_number = 10,
#   main = "Spearman correlation",
#   width = 12,
#   height = 12
# )) +
#   theme(
#     panel.grid = element_blank(),
#     plot.margin = unit(c(1, 1, 0, 1), 'lines')
#   )

# gene_correlation

##----- fig. 4d. Heatmap of confirmed markers in TEPs -------

## Heatmap of confirmed markers
n <- deg_expr %>% filter(rownames(.) %in% confirmed.var)
n <- log2(n + 1)
n <- t(scale(t(n)))
n[n > 2] <- 2
n[n < -2] <- -2

## column annotations
colD <- data.frame(group = group_list)
rownames(colD) <- colnames(n)
# annotation colors
ann_colors <- list(
  group = c(Tumor = "#E41A1C", Normal = "#377EB8")
)

## plot
marker_pheatmap <- as.ggplot(pheatmap(
  n,
  cluster_rows = T,
  cluster_cols = F,
  #scale = "row",
  border_color = "white",
  legend = T,
  labCol = T,
  fontsize = 15,
  annotation_legend = T,
  show_colnames = F,
  show_rownames = T,
  annotation_col = colD,
  annotation_colors = ann_colors
)) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = unit(c(1, 1, 1, 1), 'lines')
  )

marker_pheatmap


## final data
final_data <- lr_data %>% select(groups, confirmed.var)

##save data
save(final_data, file = "final_data.rdata")


##----- fig. 4e. Enrichment analysis of lncRNA co-expression genes -------

## load data
load("final_data.rdata")
load("function_data.Rdata")

## load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggsci)
library(RColorBrewer)
library(bruceR)
library(scales)

func_gene <- str_trim(func$mRNA, 'both')

# Convert gene names to gene IDs (ENTREZID)
gene <- bitr(
  func_gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)

## KEGG enrichment analysis
kegg <- enrichKEGG(
  gene = gene$ENTREZID,
  keyType = "kegg",
  organism = 'hsa',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  qvalueCutoff = 0.1
)

dotplot(kegg, showCategory = 20) +
  ggtitle("KEGG Pathway") +
  theme(plot.title = element_text(size = 15, hjust = 0.5))

kegg_results <- kegg@result

# Convert fraction string to decimal
calculate_decimal <- function(str) {
  parts <- str_split(str, "/")[[1]]
  numerator <- as.numeric(parts[1])
  denominator <- as.numeric(parts[2])
  decimal <- round(numerator / denominator, 2)
  return(decimal)
}

kegg_results_p <- kegg_results %>%
  mutate(GeneRatio = map_dbl(GeneRatio, calculate_decimal)) %>%
  filter(p.adjust < 0.05)


# Set the number of colors in the palette
colourCount <- length(unique(kegg_results_p$Description))
getPalette <- colorRampPalette(rev(brewer.pal(9, "Set1")))

# fig. 4e. plot
kegg_results <- kegg_results_p %>%
  mutate(
    Description = factor(
      Description,
      levels = unique(Description)[order(GeneRatio, decreasing = F)]
    )
  ) %>%
  arrange(desc(GeneRatio)) %>%
  slice_head(n = 20) %>% # top 20
  ggplot(aes(
    x = Description,
    y = GeneRatio,
    color = Description,
    size = Count
  )) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  geom_segment(aes(xend = Description, yend = 0), size = 1.5) +
  coord_flip() +
  theme_minimal() +
  scale_color_manual(values = getPalette(colourCount)) +
  labs(
    title = "",
    x = "",
    y = "GeneRatio"
  ) +
  guides(color = "none") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.ticks.length.x = unit(0.02, 'cm'),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    axis.text.y = element_text(size = 15, colour = "black"),
    plot.margin = unit(c(0, 1, 1, 1), 'lines')
  )


kegg_results


# fig. 4d-e combination
fig4.d_e <- plot_grid(
  marker_pheatmap,
  kegg_results,
  labels = LETTERS[4:5],
  ncol = 2,
  label_fontface = "bold",
  label_size = 20,
  scale = 1,
  rel_widths = c(1.4, 1.8)
)


## ------- fig 4 plot grid -------
library(patchwork)

fig4 <- plot_grid(
  fig4.a,
  fig4.b,
  fig4.c,
  fig4.d_e,
  labels = LETTERS[1:3],
  ncol = 1,
  label_fontface = "bold",
  label_size = 20,
  scale = 1,
  rel_heights = c(0.2, 0.25, 0.3, 0.25)
)
fig4
ggsave("../results/p4.plt_expr_kegg1.PDF", width = 17, height = 21)
ggsave(
  '../results/p4.plt_expr_kegg1.jpg',
  width = 17,
  height = 20,
  bg = "white",
  dpi = 400
)
