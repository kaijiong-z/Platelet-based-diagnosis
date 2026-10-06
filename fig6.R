# =======================================================================
# Figure 6 — Diagnostic performance of the locked lncRNA model versus
#            CA15-3 and CEA across three cohorts, and baseline
#            characteristics of the study population
# -----------------------------------------------------------------------
# Part A  Discrimination, clinical utility and calibration
#   6a-c  ROC curves (training / testing / external validation)
#   6d-f  Decision curves
#   Figure S3  Calibration curves
#   Table 1   Performance of the three predictors in all three cohorts

# Part B  Baseline characteristics and cohort composition
#   Table s2  Clinical characteristics by cohort (tumour / non-tumour)
#   Figure S4  Cohort composition and confusion matrices
#   Figure S5  lncRNA expression in the external validation cohort
# =======================================================================

rm(list = ls())

library(mlr3verse)
library(mlr3extralearners)
library(mlr3viz)
library(mlr3mbo)
library(future)
library(tidyverse)
library(cowplot)
library(ggsignif)
library(compareGroups)

load("model_data0521.rdata")
load("final_learner.rdata")

lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

source("plot_evaluate_functions.R")  # plot_roc_curves / plot_dca_curves / plot_cal_curves
source("boot_ci_and_compare.R")      # boot_ci_metrics

seed <- 1234
set.seed(seed)

cv_folds <- 5
n_evals <- 50

base_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.background = element_blank()
  )

## ======================================================================
## Part A. Diagnostic performance across the three cohorts
## ======================================================================

# Predictions arising from cross-validation are pooled in fold order, whereas
# predictions from a fitted model follow task order. Sorting by row id puts
# every prediction vector in the same patient order, which is a prerequisite
# for any paired comparison between predictors.
sort_pred <- function(pred) {
  ord <- order(pred$row_ids)
  PredictionClassif$new(
    row_ids = pred$row_ids[ord],
    truth = pred$truth[ord],
    response = pred$response[ord],
    prob = pred$prob[ord, , drop = FALSE]
  )
}

# ----------------------------------------------------------------------
# 1. train/test/validation set lncRNA tasks
# ----------------------------------------------------------------------
# Training/testing lncRNA tasks and id
trainset_id <- trainset$data(cols = c("sample_id", "groups"))
testset_id <- testset$data(cols = c("sample_id", "groups"))

# validation lncRNA task and id
val_task <- as_task_classif(
  validation_pcr_data,
  target = "groups",
  positive = "tumor"
)
val_task$set_col_roles("sample_id", roles = "name")

val_id <- val_task$data(cols = c("sample_id", "groups"))


# ----------------------------------------------------------------------
# 2. Serum markers aligned to each cohort
# ----------------------------------------------------------------------
zlyy_marker_data <- bind_rows(
  zlyy_cli_tumor %>% select(sample_id, groups, CA153, CEA),
  zlyy_cli_normal %>% select(sample_id, groups, CA153, CEA)
)

val_marker_data <- bind_rows(
  validation_cli_tumor %>% select(sample_id, groups, CA153, CEA),
  validation_cli_normal %>% select(sample_id, groups, CA153, CEA)
)

train_marker <- trainset_id %>%
  as_tibble() %>%
  left_join(zlyy_marker_data, by = c("sample_id", "groups"))

test_marker <- testset_id %>%
  as_tibble() %>%
  left_join(zlyy_marker_data, by = c("sample_id", "groups"))

val_marker <- val_id %>%
  as_tibble() %>%
  left_join(val_marker_data, by = c("sample_id", "groups"))


# ----------------------------------------------------------------------
# 3. Single-marker tasks
# ----------------------------------------------------------------------
classif_marker_task <- function(marker_df, marker) {
  df <- marker_df %>%
    dplyr::select(sample_id, groups, all_of(marker))
  
  task <- as_task_classif(df, target = "groups", positive = "tumor")
  task$set_col_roles("sample_id", roles = "name")
  task
}

train_task_ca153 <- classif_marker_task(train_marker, "CA153")
train_task_cea <- classif_marker_task(train_marker, "CEA")

test_task_ca153 <- classif_marker_task(test_marker, "CA153")
test_task_cea <- classif_marker_task(test_marker, "CEA")

val_task_ca153 <- classif_marker_task(val_marker, "CA153")
val_task_cea <- classif_marker_task(val_marker, "CEA")


# ----------------------------------------------------------------------
# 4. Predictions under a common protocol
# ----------------------------------------------------------------------
# preprocessing
data_prep <- po("removeconstants") %>>% po("scale")

# best model (Bayesian optimization + 5-CV)
rf_learner <- auto_tuner(
  tuner = tnr("mbo"),
  learner = lrn("classif.ranger", predict_type = "prob"),
  resampling = rsmp("cv", folds = cv_folds),
  measure = msr("classif.auc"),
  terminator = trm("evals", n_evals = n_evals),
  search_space = ps(
    num.trees = p_int(200, 2000),
    mtry.ratio = p_dbl(0.1, 1.0),
    min.node.size = p_int(3, 30),
    sample.fraction = p_dbl(0.4, 0.9)
  ),
  store_tuning_instance = TRUE
)
rf_learner <- as_learner(data_prep$clone() %>>% rf_learner)

## ---- Training cohort ----
# Out-of-fold predictions of the selected lncRNA model.
idx <- which(bmr$aggregate()$learner_id == best_learner_id)
train_lnc_pred <- sort_pred(bmr$resample_result(idx)$prediction())

set.seed(seed)
train_ca153_pred <- resample(
  task = train_task_ca153,
  learner = rf_learner,
  resampling = rsmp("cv", folds = cv_folds),
  store_models = FALSE
)$prediction() %>%
  sort_pred()

set.seed(seed)
train_cea_pred <- resample(
  task = train_task_cea,
  learner = rf_learner,
  resampling = rsmp("cv", folds = cv_folds),
  store_models = FALSE
)$prediction() %>%
  sort_pred()

## Marker models fitted on the complete training cohort 
final_ca153 <- rf_learner$clone()$train(train_task_ca153)
final_cea <- rf_learner$clone()$train(train_task_cea)

## ---- Testing cohort ----
test_lnc_pred <- sort_pred(best_learner$predict(testset))
test_ca153_pred <- sort_pred(final_ca153$predict(test_task_ca153))
test_cea_pred <- sort_pred(final_cea$predict(test_task_cea))

## ---- External validation cohort ----
val_lnc_pred <- sort_pred(best_learner$predict(val_task))
val_ca153_pred <- sort_pred(final_ca153$predict(val_task_ca153))
val_cea_pred <- sort_pred(final_cea$predict(val_task_cea))

train_pred_list <- list(
  lncRNAs = train_lnc_pred,
  CA153 = train_ca153_pred,
  CEA = train_cea_pred
)
test_pred_list <- list(
  lncRNAs = test_lnc_pred,
  CA153 = test_ca153_pred,
  CEA = test_cea_pred
)
val_pred_list <- list(
  lncRNAs = val_lnc_pred,
  CA153 = val_ca153_pred,
  CEA = val_cea_pred
)


# Patient-level predictions are exported for the subgroup analyses reported in Figure 7.
marker_tbl <- function(pred, task) {
  tibble(
    sample_id = task$data(cols = "sample_id")$sample_id[
      match(pred$row_ids, task$row_ids)
    ],
    truth = pred$truth,
    response = pred$response,
    prob = pred$prob[, "tumor"]
  )
}

ca153_honest <- bind_rows(
  marker_tbl(train_ca153_pred, train_task_ca153),
  marker_tbl(test_ca153_pred, test_task_ca153)
)
cea_honest <- bind_rows(
  marker_tbl(train_cea_pred, train_task_cea),
  marker_tbl(test_cea_pred, test_task_cea)
)

save(ca153_honest, cea_honest, file = "marker_honest_preds.rdata")


# ----------------------------------------------------------------------
# 5. ROC, decision and calibration curves
# ----------------------------------------------------------------------

## ---- Training cohort ----
train_roc_p <- plot_roc_curves(
  train_pred_list,
  type = "overlay",
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)

train_dca_raw <- plot_dca_curves(
  train_pred_list,
  type = "overlay",
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)
train_benefit <- attr(train_dca_raw, "benefit_range")
train_dca_p <- train_dca_raw + theme(legend.position = c(0.76, 0.75))

train_cal_p <- plot_cal_curves(
  train_pred_list,
  n_bins = 4,
  type = "overlay",
  error_bars = FALSE,
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)
train_calibration <- attr(train_cal_p, "calibration_metrics")

## ---- Testing cohort ----
test_roc_p <- plot_roc_curves(
  test_pred_list,
  type = "overlay",
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)

test_dca_raw <- plot_dca_curves(
  test_pred_list,
  type = "overlay",
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)
test_benefit <- attr(test_dca_raw, "benefit_range")
test_dca_p <- test_dca_raw + theme(legend.position = c(0.76, 0.75))

test_cal_p <- plot_cal_curves(
  test_pred_list,
  n_bins = 4,
  type = "overlay",
  error_bars = FALSE,
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)
test_calibration <- attr(test_cal_p, "calibration_metrics")

## ---- External validation cohort ----
val_roc_p <- plot_roc_curves(
  val_pred_list,
  type = "overlay",
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)

val_dca_raw <- plot_dca_curves(
  val_pred_list,
  type = "overlay",
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)
val_benefit <- attr(val_dca_raw, "benefit_range")
val_dca_p <- val_dca_raw + theme(legend.position = c(0.76, 0.75))

val_cal_p <- plot_cal_curves(
  val_pred_list,
  n_bins = 4,
  type = "overlay",
  error_bars = FALSE,
  colors = c("#b12115", "#f29325", "#075f75"),
  theme_set = base_theme
)
val_calibration <- attr(val_cal_p, "calibration_metrics")

# Net-benefit ranges and calibration statistics
benefit_summary <- bind_rows(
  train_benefit %>% mutate(cohort = "Training"),
  test_benefit %>% mutate(cohort = "Testing"),
  val_benefit %>% mutate(cohort = "External validation")
) %>%
  relocate(cohort)

calibration_summary <- bind_rows(
  train_calibration %>% mutate(cohort = "Training"),
  test_calibration %>% mutate(cohort = "Testing"),
  val_calibration %>% mutate(cohort = "External validation")
) %>%
  relocate(cohort)

benefit_summary
calibration_summary


# ----------------------------------------------------------------------
# Figure 6 assembly
# ----------------------------------------------------------------------
fig6 <- plot_grid(
  train_roc_p,
  test_roc_p,
  val_roc_p,
  train_dca_p,
  test_dca_p,
  val_dca_p,
  labels = LETTERS[1:6],
  ncol = 3,
  rel_widths = c(1, 1, 1),
  label_fontface = "bold",
  label_size = 15,
  scale = 1
) +
  theme(plot.margin = unit(c(0, 1, 1, 1), "lines"))

fig6

ggsave("../results/p6.roc_dca_result.pdf", fig6, width = 16, height = 8.5)
ggsave("../results/p6.roc_dca_result.jpg", fig6, width = 16, height = 8.5,
       bg = "white", dpi = 400)


# ----------------------------------------------------------------------
# Figure S3. Calibration curves
# ----------------------------------------------------------------------
fig.s3 <- plot_grid(
  train_cal_p,
  test_cal_p,
  val_cal_p,
  labels = LETTERS[1:3],
  ncol = 3,
  label_fontface = "bold",
  label_size = 15
)

fig.s3

ggsave("../results/fig.s3.calibration.pdf", fig.s3, width = 15, height = 4.5)
ggsave("../results/fig.s3.calibration.jpg", fig.s3, width = 15, height = 4.5,
       bg = "white", dpi = 400)


# ----------------------------------------------------------------------
# Table 1. Performance of the three predictors in all cohorts
# ----------------------------------------------------------------------
create_performance_table <- function(pred_list, boot_seed = 1234) {
  set.seed(boot_seed)
  imap(pred_list, \(pred, name) boot_ci_metrics(pred, name)) %>%
    list_rbind() %>%
    dplyr::filter(metric != "ce") %>%
    mutate(
      value = if_else(
        metric %in% c("auc", "sensitivity", "specificity"),
        label,
        sprintf("%.3f", estimate)
      )
    ) %>%
    dplyr::select(Models, metric, value) %>%
    pivot_wider(names_from = metric, values_from = value)
}

allset_perform_table <- bind_rows(
  create_performance_table(train_pred_list) %>% mutate(dataset = "Training"),
  create_performance_table(test_pred_list) %>% mutate(dataset = "Testing"),
  create_performance_table(val_pred_list) %>% mutate(dataset = "Validation")
) %>%
  relocate(sensitivity,specificity,.before = "acc")

allset_perform_table

print_table(
  allset_perform_table,
  file = "../results/Table1.all_dataset_performance.doc"
)


## ======================================================================
## Part B. Baseline characteristics and cohort composition
## ======================================================================

# ----------------------------------------------------------------------
# 1. Table s2. Clinical characteristics by cohort
# ----------------------------------------------------------------------
train_tumor <- zlyy_cli_tumor %>%
  right_join(
    trainset_id %>% filter(groups == "tumor"),
    by = c("sample_id", "groups")
  ) %>%
  select(age, stage, BI_RADS, subtype, CA153, CEA) %>%
  mutate(across(c(stage, BI_RADS, subtype), ~ replace_na(as.character(.x), "NA"))) %>%
  mutate(group = "Training set")

train_normal <- zlyy_cli_normal %>%
  right_join(
    trainset_id %>% filter(groups != "tumor"),
    by = c("sample_id", "groups")
  ) %>%
  select(age, CA153, CEA) %>%
  mutate(group = "Training set")

test_tumor <- zlyy_cli_tumor %>%
  right_join(
    testset_id %>% filter(groups == "tumor"),
    by = c("sample_id", "groups")
  ) %>%
  select(age, stage, BI_RADS, subtype, CA153, CEA) %>%
  mutate(across(c(stage, BI_RADS, subtype), ~ replace_na(as.character(.x), "NA"))) %>%
  mutate(group = "Testing set")

test_normal <- zlyy_cli_normal %>%
  right_join(
    testset_id %>% filter(groups != "tumor"),
    by = c("sample_id", "groups")
  ) %>%
  select(age, CA153, CEA) %>%
  mutate(group = "Testing set")

validation_tumor <- validation_cli_tumor %>%
  filter(groups == "tumor") %>%
  select(age, stage, BI_RADS, subtype, CA153, CEA) %>%
  mutate(across(c(stage, BI_RADS, subtype), ~ replace_na(as.character(.x), "NA"))) %>%
  mutate(group = "External validation set")

validation_normal <- validation_cli_normal %>%
  filter(groups == "normal") %>%
  select(age, CA153, CEA) %>%
  mutate(group = "External validation set")

group_levels <- c("Training set", "Testing set", "External validation set")

pcr_clin_tumor <- bind_rows(train_tumor, test_tumor, validation_tumor) %>%
  mutate(group = factor(group, levels = group_levels))

pcr_clin_normal <- bind_rows(train_normal, test_normal, validation_normal) %>%
  mutate(group = factor(group, levels = group_levels))

pcr_clin_tumor_table <- descrTable(
  group ~ .,
  data = pcr_clin_tumor,
  method = NA,
  show.all = FALSE,
  p.corrected = TRUE,
  show.p.overall = FALSE
)

export2word(
  pcr_clin_tumor_table,
  file = "../results/TableS2.clinical_tumor.docx"
)

pcr_clin_normal_table <- descrTable(
  group ~ .,
  data = pcr_clin_normal,
  method = NA,
  show.all = FALSE,
  p.corrected = TRUE,
  show.p.overall = FALSE
)

export2word(
  pcr_clin_normal_table,
  file = "../results/TableS2.clinical_normal.docx"
)


# ----------------------------------------------------------------------
# 2. Cohort composition
# ----------------------------------------------------------------------
plot_per_chart <- function(data, fill_colors, title) {
  data %>%
    ggplot(aes(x = "", y = value, fill = group)) +
    geom_bar(stat = "identity", width = 2, color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(x = "", y = "", title = title) +
    geom_text(
      aes(y = value, label = paste0(group, "\n", value)),
      color = "white",
      size = 8,
      position = position_stack(vjust = 0.5)
    ) +
    scale_fill_manual(values = fill_colors) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      legend.position = "none"
    )
}

fill_colors <- c("#007172", "#b12115")

train_df <- trainset$data() %>%
  mutate(groups = ifelse(groups == "tumor", "tumor", "nontumor")) %>%
  count(group = groups, name = "value")

test_df <- testset$data() %>%
  mutate(groups = ifelse(groups == "tumor", "tumor", "nontumor")) %>%
  count(group = groups, name = "value")

validation_df <- bind_rows(
  validation_cli_tumor %>% select(sample_id, groups),
  validation_cli_normal %>% select(sample_id, groups)
) %>%
  mutate(groups = ifelse(groups == "tumor", "tumor", "nontumor")) %>%
  count(group = groups, name = "value")

train_per <- plot_per_chart(train_df, fill_colors, "Training Set")
test_per <- plot_per_chart(test_df, fill_colors, "Testing Set")
validation_per <- plot_per_chart(validation_df, fill_colors, "External Validation Set")


# ----------------------------------------------------------------------
# 3. Confusion matrices of the selected model
# ----------------------------------------------------------------------
plot_confusion_matrix <- function(confusion_data, title = "") {
  ggplot(confusion_data, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = color)) +
    scale_fill_manual(values = c("#007172", "grey100", "grey100", "#b12115")) +
    geom_text(aes(label = percent), size = 8) +
    geom_vline(xintercept = c(0.5, 1.5, 2.5), linewidth = 0.5) +
    geom_hline(yintercept = c(0.5, 1.5, 2.5), linewidth = 0.5) +
    labs(title = title, x = "Actual", y = "Predicted") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_text(size = 18, colour = "black"),
      axis.text = element_text(size = 15, colour = "black"),
      legend.position = "none"
    )
}


confusion_tbl <- function(pred) {
  pred$confusion %>%
    as_tibble() %>%
    mutate(
      response = forcats::fct_recode(response, nontumor = "normal"),
      truth = forcats::fct_recode(truth, nontumor = "normal"),
      Prediction = factor(response, levels = c("nontumor", "tumor")),
      Reference = factor(truth, levels = c("tumor", "nontumor")),
      percent = n,
      color = paste(Prediction, Reference, sep = "_")
    )
}

train_confusion_p <- plot_confusion_matrix(confusion_tbl(train_lnc_pred))
test_confusion_p <- plot_confusion_matrix(confusion_tbl(test_lnc_pred))
val_confusion_p <- plot_confusion_matrix(confusion_tbl(val_lnc_pred))

upper_per <- plot_grid(
  train_per, test_per, validation_per,
  labels = c("A", "C", "E"),
  ncol = 3,
  label_fontface = "bold",
  label_size = 20,
  scale = 1
)

bottom_conf <- plot_grid(
  train_confusion_p, test_confusion_p, val_confusion_p,
  labels = c("B", "D", "F"),
  ncol = 3,
  label_fontface = "bold",
  label_size = 20,
  scale = 1
) +
  theme(plot.margin = unit(c(0, 1, 1, 1), "lines"))

fig.s4 <- plot_grid(upper_per, bottom_conf, nrow = 2)

fig.s4

ggsave("../results/fig.s4.cohort_confusion.pdf", fig.s4, width = 16, height = 8)
ggsave("../results/fig.s4.cohort_confusion.jpg", fig.s4, width = 16, height = 8,
       bg = "white", dpi = 400)


# ----------------------------------------------------------------------
# 4. Figure S5. lncRNA expression in the external validation cohort
# ----------------------------------------------------------------------
plot_validation_boxplot <- function(data, column) {
  ggplot(data, aes(x = groups, y = .data[[column]])) +
    geom_boxplot(width = 0.7, outlier.shape = NA) +
    geom_jitter(aes(colour = groups), size = 2.5, width = 0.2, alpha = 0.9) +
    scale_y_log10() +
    scale_color_manual(values = c("nontumor" = "#075f75", "tumor" = "#b12115")) +
    geom_signif(
      comparisons = list(c("tumor", "nontumor")),
      map_signif_level = function(p) sprintf("p = %.2g", p),
      test = "wilcox.test",
      textsize = 4, tip_length = 0, margin_top = 0.05, vjust = 0
    ) +
    labs(
      title = column,
      x = "",
      y = expression("Expression value (" * 2^{-Delta * Ct} * ")")
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      legend.position = "none"
    )
}

lncrna_vars <- c(
  "FTX", "Z93241.1", "BAIAP2_DT",
  "AC008083.4", "LINC01275", "AL360012.1"
)

validation_box_data <- validation_pcr_data %>%
  mutate(groups = ifelse(groups == "tumor", "tumor", "nontumor"))

validation_boxplots <- map(
  lncrna_vars,
  \(v) plot_validation_boxplot(validation_box_data, v)
)

fig.s5 <- plot_grid(
  plotlist = validation_boxplots,
  labels = LETTERS[1:6],
  ncol = 3,
  label_fontface = "bold",
  label_size = 20,
  scale = 1
)

fig.s5

ggsave("../results/fig.s5.validation_expression.pdf", fig.s5, width = 12, height = 8.5)
ggsave("../results/fig.s5.validation_expression.jpg", fig.s5, width = 12, height = 8.5,
       bg = "white", dpi = 300)
