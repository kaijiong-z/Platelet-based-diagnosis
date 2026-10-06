# =======================================================================
# Fig. 3. the differential expression analysis of platelet RNA-sequencing
# =======================================================================

rm(list = ls())

## load package
pacman::p_load(
  tidyverse,
  limma,
  pheatmap,
  viridis,
  cowplot,
  edgeR,
  gridExtra,
  ggsci,
  data.table
)

##========== step 1: platelet lncRNAs annotation and distribution  ============

##--  groups and annotation information --

##-- PLT groups --
## tumor
plt_tumor_names <- c(
  "A383-03",
  "A383-26",
  "A383-27", ##HER2
  "A383-24",
  "A383-25",
  "A383-05",
  "A383-S01-17", ##LuminalB
  "A383-23",
  "A383-28",
  "A383-S01-16",
  "A383-S01-21", ##luminalA
  "A383-S01-18",
  "A383-S01-19",
  "A383-S01-15" ##TNBC
)

## normal
plt_normal_names <- c(
  "A352-07",
  "A352-08",
  "A383-S01-08",
  "A383-S01-09",
  "A383-S01-10"
)


## group list
group_list <- rep(c('Tumor', 'Normal'), c(14, 5))
group_list <- factor(group_list, levels = c("Normal", "Tumor"))
table(group_list)

## save groups
save(plt_tumor_names, plt_normal_names, group_list, file = "group_list.rdata")


##-- Gene Annotation --

rm(list = ls())

## read gene annotation file
gene_anno <- read.table("../data/gencode.gene.info.v22.tsv", header = T)

## lncRNA gene types
lncRNA <- c(
  '3prime_overlapping_ncRNA',
  'antisense',
  'bidirectional_promoter_lncRNA',
  'lincRNA',
  'macro_lncRNA',
  'non_coding',
  'processed_transcript',
  'sense_intronic',
  'sense_overlapping'
)

## mRNA gene types
mRNA <- c(
  "protein_coding"
)


## lncRNA gene annotation list
lncRNA_anno_list <- gene_anno %>%
  dplyr::select(c(gene_id, gene_name, gene_type)) %>%
  dplyr::filter(gene_type %in% lncRNA)


## mRNA gene annotation list
mRNA_anno_list <- gene_anno %>%
  dplyr::select(c(gene_id, gene_name, gene_type)) %>%
  dplyr::filter(gene_type %in% mRNA)

## save gene anno
save(lncRNA_anno_list, mRNA_anno_list, file = "gene_anno_list.rdata")


# ---- fig. 3a. lncRNA type distribution ----

rm(list = ls())

## load groups data
load("group_list.rdata")
load("gene_anno_list.rdata")

## Read counts expr data
lncRNA_counts <- read_delim("../data/lncRNA_counts.xls", show_col_types = FALSE)

## lncrna percent
gene_percent <- lncRNA_counts %>%
  mutate(
    gene_type = ifelse(
      str_detect(Symbol, pattern = "MSTRG."),
      "newgene",
      "lncRNA"
    )
  ) %>%
  select(gene_type) %>%
  group_by(gene_type) %>% # Variable to be transformed
  dplyr::count() %>%
  ungroup() %>%
  mutate(perc = n / sum(n)) %>%
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

## lncrna percent plot
gene_percent_p <- gene_percent %>%
  ggplot(aes(x = "", y = perc, fill = gene_type)) +
  geom_col(color = "white") +
  geom_label(
    aes(label = labels),
    color = c("white", "white"),
    size = 6,
    position = position_stack(vjust = 0.5),
    show.legend = FALSE
  ) +
  #scale_fill_viridis_d() +
  scale_fill_manual(values = c("#075f75", "#b12115")) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(
    #legend.key.size = unit(0.8, "cm"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )
gene_percent_p


# ---- fig. 3b. lncRNA subgroup distribution ----

## lncrna subgroup percent
plt_lncrna_counts <- lncRNA_counts %>%
  select(
    c("#ID", "Symbol"),
    all_of(plt_tumor_names),
    all_of(plt_normal_names)
  ) %>%
  filter(!str_detect(Symbol, pattern = "MSTRG."))


plt_lncRNA_counts_expr <- lncRNA_anno_list %>%
  separate(gene_id, into = c("gene_id"), sep = "\\.") %>%
  dplyr::inner_join(plt_lncrna_counts, by = c("gene_id" = "#ID")) %>% #与表达谱合并
  mutate(avg = rowMeans(dplyr::select(., 5:ncol(.)))) %>% #计算每个基因平均表达量avg
  arrange(desc(avg)) %>% #根据avg列对表达矩阵降序
  distinct(gene_name, .keep_all = T) %>%
  dplyr::select(-avg)


## lncrna subgroup
lncrna_percent <- plt_lncRNA_counts_expr %>%
  mutate(
    lncRNA_type = factor(
      case_when(
        gene_type == "antisense" ~ "antisense",
        gene_type == "lincRNA" ~ "lincRNA",
        gene_type == "sense_intronic" ~ "sense_intr",
        gene_type == "processed_transcript" ~ "pro_trans",
        TRUE ~ "others"
      ),
      levels = c("lincRNA", "antisense", "sense_intr", "pro_trans", "others")
    )
  ) %>%
  group_by(lncRNA_type) %>%
  count() %>%
  ungroup() %>%
  mutate(perc = n / sum(n)) %>%
  arrange(perc) %>%
  mutate(labels = scales::percent(round(perc, 3)))

lncrna_percent


##lncrna sugroup percent plot
lncrna_percent_p <- lncrna_percent %>%
  ggplot(aes(x = "", y = perc, fill = lncRNA_type)) +
  geom_col(color = "white") +
  geom_label(
    aes(label = labels),
    color = c("white", "white", "white", "white", "white"),
    size = 6,
    position = position_stack(vjust = 0.5),
    show.legend = FALSE
  ) +
  #scale_fill_viridis_d() +
  scale_fill_manual(
    values = c(
      "#E41A1C",
      "#7F3C8D",
      "#11A579",
      "#3969AC",
      "#F2B701"
    )
  ) +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(
    #legend.key.size = unit(0.8, "cm"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )
# guides(fill = guide_legend(ncol = 2,byrow = TRUE,
#                            title.position = "top",
#                            title.hjust = 0.5))

lncrna_percent_p


## save lncRNA counts data
save(plt_lncRNA_counts_expr, file = "plt_lncRNA_counts_expr.rdata")


##========== step 2: Platelet differential expression analysis  ============

## ----- 1. Platelet expression matrix   ------

rm(list = ls())

##load package
library(tidyverse)

## load groups data
load("group_list.rdata")
load("gene_anno_list.rdata")

## PLT counts expr rawdata
lncRNA_counts <- read_delim("../data/lncRNA_counts.xls", show_col_types = FALSE)


plt_lncrna_counts <- lncRNA_counts %>%
  select(
    c("#ID", "Symbol"),
    all_of(plt_tumor_names),
    all_of(plt_normal_names)
  ) %>%
  filter(!str_detect(Symbol, pattern = "MSTRG."))

plt_lncRNA_counts_expr <- lncRNA_anno_list %>%
  separate(gene_id, into = c("gene_id"), sep = "\\.") %>%
  dplyr::inner_join(plt_lncrna_counts, by = c("gene_id" = "#ID")) %>%
  mutate(avg = rowMeans(dplyr::select(., 5:ncol(.)))) %>%
  arrange(desc(avg)) %>%
  distinct(Symbol, .keep_all = T) %>%
  dplyr::select(-c(gene_id, gene_type, gene_name, avg)) %>%
  column_to_rownames("Symbol")

plt_lncRNA_counts_expr[1:6, 1:6]


## PLT fpkm expr data
lncRNA_fpkm <- read_delim(
  "../data/lncRNA_expression.xls",
  show_col_types = FALSE
)

plt_lncRNA_fpkm <- lncRNA_fpkm %>%
  select(
    c("#ID", "Symbol"),
    all_of(plt_tumor_names),
    all_of(plt_normal_names)
  ) %>%
  filter(!str_detect(Symbol, pattern = "MSTRG."))


plt_lncRNA_fpkm_expr <- lncRNA_anno_list %>%
  separate(gene_id, into = c("gene_id"), sep = "\\.") %>%
  dplyr::inner_join(plt_lncRNA_fpkm, by = c("gene_id" = "#ID")) %>%
  mutate(avg = rowMeans(dplyr::select(., 5:ncol(.)))) %>%
  arrange(desc(avg)) %>%
  distinct(Symbol, .keep_all = T) %>%
  dplyr::select(-c(gene_id, gene_type, gene_name, avg)) %>%
  column_to_rownames("Symbol")

plt_lncRNA_fpkm_expr[1:6, 1:6]

## FPKM to TPM
FPKMToTPM <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

plt_lncRNA_tpm_expr <- plt_lncRNA_fpkm_expr %>%
  apply(., 2, FPKMToTPM) %>%
  as.data.frame()

## save expr data
save(
  plt_lncRNA_counts_expr,
  plt_lncRNA_tpm_expr,
  file = "plt_lncRNA_expr.rdata"
)

## ----- 2. Platelet differential analysis   ------

rm(list = ls())

## load package
library(tidyverse)
library(stringr)
library(edgeR)
library(ggrepel)
library(cowplot)

## load data
load("plt_lncRNA_expr.rdata")
load("group_list.rdata")

## filter low-expression genes
expr <- plt_lncRNA_counts_expr %>% filter(rowSums(. > 1) >= 0.50 * ncol(.))


##== Data preprocessing ==

##(1) Construct DGEList object
dgelist <- DGEList(counts = expr, group = group_list)

##(2) Filter low count data
keep <- rowSums(cpm(dgelist) > 1) >= 2

dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist$samples$lib.size <- colSums(dgelist$counts)

##(3) Normalization, using TMM as an example
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

##== Differential analysis ==

##(1) Construct design matrix
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(dgelist)

##(2) Estimate gene expression dispersion
dge <- estimateDisp(dgelist_norm, design)
fit <- glmFit(dge, design)
fit2 <- glmLRT(fit, contrast = c(-1, 1))

##(3) Filter differentially expressed genes
DEG <- topTags(fit2, n = nrow(dge))
DEG <- as.data.frame(DEG)
head(DEG)

##(4) Set thresholds
#FDR<0.05, fold change(FC)>2, i.e., logFC>1
PValue_cutoff <- 0.05
logFC_cutoff <- 1
DEG$change <- as.factor(ifelse(
  DEG$PValue > PValue_cutoff,
  'No',
  ifelse(
    DEG$logFC > logFC_cutoff,
    'Up',
    ifelse(DEG$logFC < -logFC_cutoff, 'Down', 'No')
  )
))

table(DEG$change)

##(5) Differentially expressed results

## all DEG
sig_DEG.name <- rownames(DEG)[DEG$change != "No"]
sig_DEG <- DEG[sig_DEG.name, ]

TEPS_DEG <- DEG

##save data
save(TEPS_DEG, sig_DEG, file = "TEPS_DEG.rdata")

##================ step3: Visualization =========================

# ---- fig. 3e. Volcano plot ----

head(DEG)
DEG$Symbol <- rownames(DEG)

Volcano_plot <- ggplot(data = DEG, aes(x = logFC, y = -log10(PValue))) +
  geom_point(alpha = 1, size = 3, aes(color = change)) +
  scale_color_manual(values = c("#075f75", "grey", "#b12115")) +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by = 2)) +
  # scale_y_continuous(#expand = expansion(add = c(0.5, 0)),
  #                    limits = c(0, 5),
  #                    breaks = seq(0, 5, by = 1)) +
  geom_vline(
    xintercept = c(-logFC_cutoff, logFC_cutoff),
    lty = 4,
    col = "grey",
    lwd = 0.8
  ) +
  geom_hline(
    yintercept = -log10(PValue_cutoff),
    lty = 4,
    col = "grey",
    lwd = 0.8
  ) +
  #scale_x_continuous(limits = c(-3.5, 3.5))+
  labs(
    title = "Diff Expressed lncRNAs",
    x = "log2(FoldChange)",
    y = "-log10(PValue)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    # legend.justification = c(0, 0), legend.position = c(0, 0),
    legend.background = element_blank()
  )

Volcano_plot


# ---- fig. 3f. histogram plot of DE-lncRNAs ----
diff_gene_P <- sig_DEG %>%
  dplyr::count(change) %>% #count DEGs
  mutate(change = factor(change, levels = c("Up", "Down"))) %>%
  ggplot(aes(x = change, y = n, fill = change)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b12115", "#075f75")) +
  # 在条形图上添加文本标注
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  scale_y_continuous(
    #expand = expansion(add = c(1, 0)),
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5)
  ) +
  labs(title = "Number of DE-lncRNAs", x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_blank()
  )


diff_gene_P

# ---- fig. 3c. PCA plot ----

#edgeR recommends using logCPM values for heatmap and PCA plotting.
logCPM <- cpm(dgelist_norm, log = TRUE, prior.count = 2) #logCPM
logCPM[1:4, 1:4]
logCPM_SigDEG_expr <- logCPM[rownames(sig_DEG), ]
logCPM_SigDEG_expr[1:4, 1:4]
dat <- logCPM_SigDEG_expr

## Plot PCA
library(FactoMineR)
library(factoextra)
dat <- as.data.frame(t(dat))
dat.pca <- PCA(dat, graph = FALSE)

pca <- fviz_pca_ind(
  dat.pca,
  geom = c("point"),
  #geom.ind = c("point","text"), # show points only (nbut not "text")
  pointsize = 5,
  pointshape = 21,
  fill.ind = group_list,
  col.ind = group_list, # color by groups
  palette = c("#075f75", "#b12115"),
  addEllipses = TRUE, # Concentration ellipses
  mean.point = F,
  legend.title = "group"
) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), 'lines')
  )

pca

# ---- fig. 3d. heatmap plot ----
library(pheatmap)
library(ggplotify)

n <- logCPM_SigDEG_expr
n <- log2(n + 0.1)
n <- t(scale(t(n)))
n[n > 2] <- 2
n[n < -2] <- -2

## Construct column annotation information
#annotation_col <- data.frame(group = group_list)
annotation_col <- data.frame(
  subtype = factor(c(
    rep(c("HER2"), 3),
    rep(c("LuminalB"), 4),
    rep(c("LuminalA"), 4),
    rep(c("TNBC"), 3),
    rep(c("Normal"), 5)
  )),

  group = factor(c(rep(c("Tumor"), 14), rep(c("Normal"), 5)))
)

rownames(annotation_col) <- colnames(n)
head(annotation_col)


# defined color list for annotation information
ann_colors <- list(
  group = c(Tumor = "#b12115", Normal = "#075f75"),
  subtype = c(
    HER2 = "#2e7ebb",
    LuminalB = "#984EA3",
    LuminalA = "#FF7F00",
    TNBC = "#4DAF4A",
    Normal = "#075f75"
  )
)
head(ann_colors)

# Plot heatmap
library(pheatmap)
library(ggplotify)
pheatmap <- as.ggplot(pheatmap(
  n,
  cluster_rows = T,
  cluster_cols = F,
  scale = "row",
  #border_color="white",
  legend = T,
  labCol = T,
  annotation_legend = T,
  show_colnames = F,
  show_rownames = F,
  annotation_col = annotation_col,
  annotation_colors = ann_colors
))

pheatmap


## correlation heatmap
cor_pheatmap <- as.ggplot(pheatmap(
  cor(n),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = T,
  show_colnames = T
))


##==fig3. all plot==
library(cowplot)

# top row
top_row <- plot_grid(
  gene_percent_p,
  lncrna_percent_p,
  pca,
  labels = c("A", 'B', 'C'),
  ncol = 3,
  align = 'h',
  rel_widths = c(1, 1, 1),
  label_fontface = "bold",
  label_size = 20,
  scale = 1
)

# bottom row
bottom_row <- plot_grid(
  pheatmap,
  Volcano_plot,
  diff_gene_P,
  labels = c('D', 'E', 'F'),
  ncol = 3,
  align = 'h',
  rel_widths = c(1, 1, 1),
  label_fontface = "bold",
  label_size = 20,
  scale = 1
)


lncrna_expr_plot <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  nrow = 2,
  rel_heights = c(1, 1.2),
  label_fontface = "bold",
  label_size = 20,
  scale = 1
) +
  theme(plot.margin = unit(c(1, 1, 1, 1), 'lines'))

print(lncrna_expr_plot)

# save plot
ggsave("../results/P2.lncrna_expr_plot.pdf", width = 17, height = 10)
ggsave(
  "../results/P2.lncrna_expr_plot.jpg",
  width = 17,
  height = 10,
  bg = "white",
  dpi = 400
)
