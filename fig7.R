# =======================================================================
# fig. 7 — lncRNA panel performance in breast cancer subgroups
# ----------------------------------------------------------------------
# fig.s6  Distribution of TNM stage, subtype and BI-RADS

# 7a-b  CA153 / CEA distribution across TNM stage
# 7c    Early stage (I-II): lncRNA panel vs CA153 vs CEA
# 7d    BI-RADS 2-4a
# 7e    Marker-negative patients
# 7f    Molecular subtypes
# =======================================================================

rm(list = ls())

library(tidyverse)
library(mlr3verse)
library(cowplot)
library(pROC)

load("model_data0521.rdata")
load("final_learner.rdata")        # trainset / testset / bmr / best_learner(_id)
load("marker_honest_preds.rdata")  # ca153_honest / cea_honest (from fig6.R)

lgr::get_logger("mlr3")$set_threshold("warn")

# theme set
base_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.background = element_blank()
  )


# ----------------------------------------------------------------------
# fig.s6. Distribution of TNM stage, subtype and BI-RADS
# ----------------------------------------------------------------------
pal5 <- c("#075f75", "#b12115", "#7363ac", "#ef9c02", "#5f5f5f")

pct_bar <- function(data, var, xlab) {
  data |>
    mutate(grp = forcats::fct_na_value_to_level({{ var }}, level = "NA")) |>
    count(grp) |>
    mutate(pct = n / sum(n) * 100) |>
    ggplot(aes(grp, pct, fill = grp)) +
    geom_col() +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_fill_manual(values = pal5) +
    labs(x = xlab, y = "Percentage") +
    base_theme +
    theme(legend.position = "none")
}

# TNM stage percent
figs6.stage <- pct_bar(zlyy_cli_tumor, stage, "TNM stage")

## breast subtype
figs6.subtype <- zlyy_cli_tumor  %>% 
  mutate(subtype = factor(subtype,levels = c("Her2","lumina A","lumina B","TNBC" , "NA"))) %>% 
  pct_bar(subtype,"Subtype groups")

## breast BI-RADS
figs6.birads <- zlyy_cli_tumor  %>% 
  mutate(BI_RADS = case_when(
    BI_RADS %in% c("BI-RADS 2", "BI-RADS 3", "BI-RADS 4", "BI-RADS 4a") ~ "Below 4a",
    BI_RADS %in% c("BI-RADS 4b", "BI-RADS 4c", "BI-RADS 5", "BI-RADS 6") ~ "Above 4a",
    .default = BI_RADS
  )) |>
  pct_bar(BI_RADS, "BI-RADS groups")

fig.s6 <- plot_grid(
  figs6.stage, figs6.subtype, figs6.birads,
  labels = LETTERS[1:3],
  ncol = 3,
  rel_widths = c(1, 1, 1),
  label_fontface = "bold",
  label_size = 20
)

fig.s6

ggsave(
  "../results/fig.s6.clinical_distribution.pdf",
  fig.s6,
  width = 16,
  height = 5
)

ggsave(
  "../results/fig.s6.clinical_distribution.jpg",
  fig.s6,
  width = 16,
  height = 5,
  bg = "white",
  dpi = 400
)



# ----------------------------------------------------------------------
# fig7.The performance of platelet lncRNA panel in BC subgroups
# ----------------------------------------------------------------------

# -- 1.ROC plotting function --
plot_roc_curves <- function(
  preds,
  positive = "tumor",
  colors = NULL,
  title = "",
  theme_fn = theme_bw()
) {
  if (inherits(preds, "PredictionClassif")) preds <- list(Model = preds)
  
  nms <- names(preds)
  neg <- setdiff(levels(preds[[1]]$truth), positive)

  rocs <- map(preds, \(p) {
    roc(
      response = p$truth, predictor = p$prob[, positive],
      levels = c(neg, positive), direction = "<", quiet = TRUE
    )
  })

  auc_lab <- imap_chr(rocs, \(r, nm) {
    v <- as.numeric(ci.auc(r))
    sprintf("%s  %.3f (%.3f-%.3f)", nm, v[2], v[1], v[3])
  })

  curve_df <- imap(rocs, \(r, nm) {
    tibble(Model = nm, fpr = 1 - r$specificities, tpr = r$sensitivities) |>
      arrange(fpr, tpr)
  }) |>
    list_rbind() |>
    mutate(Model = factor(Model, levels = nms))

  if (is.null(colors)) {
    colors <- set_names(scales::hue_pal()(length(nms)), nms)
  } else if (is.null(names(colors))) {
    colors <- set_names(colors, nms)
  }

  ggplot(curve_df, aes(fpr, tpr, color = Model)) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "grey70", linewidth = 1
    ) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(
      values = colors, labels = auc_lab[nms], name = "AUC (95% CI)"
    ) +
    coord_fixed() +
    labs(x = "1 - Specificity", y = "Sensitivity", title = title) +
    theme_fn +
    theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0))
}


# -- 2. Out-of-fold / held-out predictions for every sample --
idx <- which(bmr$aggregate()$learner_id == best_learner_id)
oof_pred <- bmr$resample_result(idx)$prediction()
tst_pred <- best_learner$predict(testset)

lnc_honest <- bind_rows(
  tibble(
    sample_id = trainset$data(cols = "sample_id")$sample_id[
      match(oof_pred$row_ids, trainset$row_ids)
    ],
    truth = oof_pred$truth,
    response = oof_pred$response,
    prob = oof_pred$prob[, "tumor"]
  ),
  tibble(
    sample_id = testset$data(cols = "sample_id")$sample_id[
      match(tst_pred$row_ids, testset$row_ids)
    ],
    truth = tst_pred$truth,
    response = tst_pred$response,
    prob = tst_pred$prob[, "tumor"]
  )
)

stopifnot(!any(is.na(lnc_honest$sample_id)))

lev <- levels(lnc_honest$truth)

# subset a prediction table by sample id
make_pred <- function(honest_tbl, ids) {
  d <- honest_tbl |> filter(sample_id %in% ids)
  prob <- cbind(d$prob, 1 - d$prob)
  colnames(prob) <- c("tumor", setdiff(lev, "tumor"))

  PredictionClassif$new(
    row_ids = seq_len(nrow(d)),
    truth = factor(d$truth, levels = lev),
    response = factor(d$response, levels = lev),
    prob = prob[, lev, drop = FALSE]
  )
}


# -- 3. Subgroup sample ids --

table(zlyy_cli_tumor$subtype, useNA = "ifany")
table(zlyy_cli_tumor$BI_RADS, useNA = "ifany")

subtype_levels <- c("lumina A", "lumina B", "Her2", "TNBC")
normal_ids <- zlyy_cli_normal$sample_id

# each subgroup = qualifying tumor samples + all normal controls
sub_ids <- list(
  `Early stage (I-II)` = zlyy_cli_tumor |>
    filter(!is.na(stage), stage %in% c("I", "II")) |>
    pull(sample_id),

  `BI-RADS 2-4a` = zlyy_cli_tumor |>
    filter(BI_RADS %in% c("BI-RADS 2", "BI-RADS 3", "BI-RADS 4", "BI-RADS 4a")) |>
    pull(sample_id),

  `CA153 negative` = zlyy_cli_tumor |>
    filter(!is.na(CA153), CA153 <= 30) |>
    pull(sample_id),

  `CEA negative` = zlyy_cli_tumor |>
    filter(!is.na(CEA), CEA <= 5) |>
    pull(sample_id),

  `Both negative` = zlyy_cli_tumor |>
    filter(!is.na(CA153), !is.na(CEA), CA153 <= 30, CEA <= 5) |>
    pull(sample_id)
) |>
  c(set_names(
    map(subtype_levels, \(s) zlyy_cli_tumor |> filter(subtype == s) |> pull(sample_id)),
    subtype_levels
  )) |>
  map(\(x) c(x, normal_ids))

sub_preds <- map(sub_ids, \(ids) make_pred(lnc_honest, ids))


# fig7a-b. CA153 / CEA distribution across TNM stage

stage_colors <- c("#075f75", "#b12115", "#7363ac", "#ef9c02")

fig7a.ca153_dist <- zlyy_cli_tumor |>
  filter(!is.na(stage), !is.na(CA153)) |>
  ggplot(aes(stage, CA153)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(aes(colour = stage), size = 3, width = 0.2, alpha = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = stage_colors) +
  labs(title = "CA153 distribution", x = "TNM stage", y = "CA153") +
  base_theme +
  theme(legend.position = "none")

fig7a.ca153_dist

fig7b.cea_dist <- zlyy_cli_tumor |>
  filter(!is.na(stage), !is.na(CEA)) |>
  ggplot(aes(stage, CEA)) +
  geom_boxplot(width = 0.7, outlier.shape = NA) +
  geom_jitter(aes(colour = stage), size = 3, width = 0.2, alpha = 0.8) +
  scale_y_log10() +
  scale_color_manual(values = stage_colors) +
  labs(title = "CEA distribution", x = "TNM stage", y = "CEA") +
  base_theme +
  theme(legend.position = "none")

fig7b.cea_dist


# fig7c. Early stage (I-II): lncRNA panel vs CA153 vs CEA

early_ids <- sub_ids[["Early stage (I-II)"]]

early_pred_list <- list(
  lncRNAs = make_pred(lnc_honest, early_ids),
  CA153 = make_pred(ca153_honest, early_ids),
  CEA = make_pred(cea_honest, early_ids)
)

# the three must share the same patients for paired comparison
map_int(early_pred_list, \(p) length(p$truth))

fig7c.early_roc <- plot_roc_curves(
  early_pred_list,
  colors = c("#b12115", "#f29325", "#075f75"),
  title = "Early stage (I-II)",
  theme_fn = base_theme
)

fig7c.early_roc


# fig7d. BI-RADS 2-4a

fig7d.birads_roc <- plot_roc_curves(
  list(`BI-RADS 2-4a` = sub_preds[["BI-RADS 2-4a"]]),
  colors = "#075f75",
  title = "BI-RADS 2-4a",
  theme_fn = base_theme
)

fig7d.birads_roc


# fig7e. Marker-negative patients

fig7e.negative_roc <- plot_roc_curves(
  sub_preds[c("CA153 negative", "CEA negative", "Both negative")],
  colors = c("#b12115", "#ef9c02", "#075f75"),
  title = "Marker-negative patients",
  theme_fn = base_theme
)

fig7e.negative_roc


# fig7f. Molecular subtypes

fig7f.subtype_roc <- plot_roc_curves(
  sub_preds[subtype_levels],
  colors = c("#b12115", "#ef9c02", "#075f75", "#5f5f5f"),
  title = "Molecular subtype",
  theme_fn = base_theme
)

fig7f.subtype_roc

# ----------------------------------------------------------------------
# fig7 plot grid
# ----------------------------------------------------------------------
fig7 <- plot_grid(
  fig7a.ca153_dist,
  fig7b.cea_dist,
  fig7c.early_roc,
  fig7d.birads_roc,
  fig7e.negative_roc,
  fig7f.subtype_roc,
  labels = LETTERS[1:6],
  ncol = 3,
  rel_widths = c(1, 1, 1),
  label_fontface = "bold",
  label_size = 20,
  scale = 1
) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

fig7

ggsave(
  "../results/p7.subgroup_roc.pdf",
  fig7,
  width = 15,
  height = 9
)

ggsave(
  "../results/p7.subgroup_roc.jpg",
  fig7,
  width = 15,
  height = 9,
  bg = "white",
  dpi = 400
)



