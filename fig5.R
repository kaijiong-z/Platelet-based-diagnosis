# =======================================================================
# Figure 5 — Platelet lncRNA expression, univariable discrimination and
#            machine-learning model development
# -----------------------------------------------------------------------
# 5a  Expression of six candidate lncRNAs (tumour vs non-tumour)
# 5b  Univariable discrimination of six lncRNAs and two serum markers
# 5c  ROC comparison of eight machine-learning algorithms (nested CV)
# 5d  Classification error across algorithms
#
# Supplementary outputs
#   Table S3  Univariable diagnostic performance
#   Table S4  Model performance with 95% confidence intervals
#   Table S5  Tuned hyper-parameters of each algorithm
#   Figure S1 Performance heat map and AUC forest plot
#   Figure S2 Permutation variable importance of the selected model
# =======================================================================

rm(list = ls())

library(tidyverse)
library(bruceR)
library(pROC)
library(ggsignif)
library(broom)
library(cowplot)

load("model_data0521.rdata")

# theme set
base_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.background = element_blank()
  )

# Candidate lncRNAs and serum tumour markers
lncrna_vars <- c(
  "FTX", "Z93241.1", "BAIAP2_DT",
  "AC008083.4", "LINC01275", "AL360012.1"
)
marker_vars <- c("CA153", "CEA")


# ----------------------------------------------------------------------
# fig5a. Expression of the six candidate lncRNAs
# ----------------------------------------------------------------------
# define plot function
plot_lncrna_boxplot <- function(data, column) {
  # Computing p-values using Wilcoxon rank sum test
  p_value <- wilcox.test(
    data[[column]][data$groups == "tumor"],
    data[[column]][data$groups == "nontumor"]
  )$p.value
  
  ggplot(data, aes(x = groups, y = .data[[column]])) +
    geom_boxplot(width = 0.7, outlier.shape = NA) +
    geom_jitter(aes(colour = groups), size = 2.5, width = 0.2, alpha = 0.9) +
    scale_y_log10(expand = expansion(mult = c(0.02, 0.15))) +
    scale_color_manual(
      values = c("tumor" = "#b12115", "nontumor" = "#075f75")
    ) +
    geom_signif(
      comparisons = list(c("tumor", "nontumor")),
      annotations = ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value)),
      test = "wilcox.test",textsize = 4, tip_length = 0, margin_top = 0.12, vjust = 0
    ) +
    labs(
      title = column,
      x = "",
      y = expression("Expression value (" * 2^{-Delta * Ct} * ")")
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text = element_text(size = 15, angle = 30, hjust = 1),
      axis.title = element_text(size = 15),
      legend.position = "none"
    )
}

data_box <- zlyy_pcr_data %>%
  mutate(
    groups = ifelse(groups == "tumor", "tumor", "nontumor"),
    groups = factor(groups, levels = c("tumor", "nontumor"))
  )

boxplot_list <- map(lncrna_vars, \(v) plot_lncrna_boxplot(data_box, v))

# fig5a
fig5a.lncRNA_pcr_expr <- plot_grid(
  plotlist = boxplot_list,
  ncol = 6,
  label_fontface = "bold",
  label_size = 18,
  scale = 1
)

fig5a.lncRNA_pcr_expr


# ----------------------------------------------------------------------
# fig5b. Univariable discrimination of lncRNAs and serum markers
# ----------------------------------------------------------------------
roc_data <- zlyy_pcr_data %>%
  left_join(
    bind_rows(
      zlyy_cli_tumor %>% select(sample_id, all_of(marker_vars)),
      zlyy_cli_normal %>% select(sample_id, all_of(marker_vars))
    ),
    by = "sample_id"
  ) %>%
  select(-sample_id) %>%
  mutate(groups = ifelse(groups == "tumor", 1, 0))

covs <- c(lncrna_vars, marker_vars)

# Define  AUC calculation function
univariable_fit <- map(set_names(covs), function(v) {
  model <- glm(
    as.formula(paste0("groups ~ ", v)),
    data = roc_data, family = binomial
  )
  
  fitted_prob <- predict(model, roc_data, type = "response")
  
  list(
    model = model,
    roc = pROC::roc(
      roc_data$groups, fitted_prob,
      levels = c(0, 1), direction = "<", quiet = TRUE
    )
  )
})

roc_list <- map(univariable_fit, "roc")

univariable_res <- imap(univariable_fit, function(fit, v) {
  auc_ci <- pROC::ci.auc(fit$roc)
  
  best <- pROC::coords(
    fit$roc,
    x = "best", input = "threshold",
    ret = c("threshold", "sensitivity", "specificity"),
    best.method = "youden", transpose = FALSE
  )
  
  ci_best <- pROC::ci.coords(
    fit$roc,
    x = "best", input = "threshold",
    ret = c("sensitivity", "specificity"),
    best.method = "youden", best.policy = "random",
    boot.n = 500, transpose = FALSE
  )
  
  # Return results as data frame
  tibble(
    variable = v,
    n_analysed = length(fit$roc$cases) + length(fit$roc$controls),
    auc = as.numeric(auc_ci[2]),
    auc_low = as.numeric(auc_ci[1]),
    auc_high = as.numeric(auc_ci[3]),
    threshold = as.numeric(best$threshold),
    sens = as.numeric(best$sensitivity),
    sens_low = ci_best$sensitivity[1, "2.5%"],
    sens_high = ci_best$sensitivity[1, "97.5%"],
    spec = as.numeric(best$specificity),
    spec_low = ci_best$specificity[1, "2.5%"],
    spec_high = ci_best$specificity[1, "97.5%"],
    p_value = tidy(fit$model) %>% filter(term != "(Intercept)") %>% pull(p.value)
  )
}) %>%
  list_rbind() %>%
  arrange(desc(auc))   

univariable_tbl <- univariable_res %>%
  transmute(
    Variable = variable,
    `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", auc, auc_low, auc_high),
    `Sensitivity (95% CI)` = sprintf("%.3f (%.3f-%.3f)", sens, sens_low, sens_high),
    `Specificity (95% CI)` = sprintf("%.3f (%.3f-%.3f)", spec, spec_low, spec_high),
    # `Youden cut-off` = sprintf("%.3f", threshold),
    `P value` = format.pval(p_value, digits = 3, eps = 1e-3)
  )

univariable_tbl

print_table(
  univariable_tbl,
  file = "../results/TableS3.univariable_performance.doc"
)


marker_colors <- set_names(
  c("#075f75", "#7363ac", "#e3550b", "#b12115",
    "#ef9c02", "#377EB8", "#2d9850", "#5f5f5f"),
  covs
)

auc_label <- set_names(
  sprintf(
    "%s = %.3f (%.3f-%.3f)",
    univariable_res$variable, univariable_res$auc,
    univariable_res$auc_low, univariable_res$auc_high
  ),
  univariable_res$variable
)

fig5b.biomarker_roc <- pROC::ggroc(
  roc_list,
  legacy.axes = TRUE,
  linewidth = 1.4
) +
  annotate(
    geom = "segment", x = 0, xend = 1, y = 0, yend = 1,
    linetype = "dashed", linewidth = 1.4, color = "grey50"
  ) +
  scale_color_manual(
    values = marker_colors[names(roc_list)],
    labels = auc_label[names(roc_list)],
    name = "Tumour vs non-tumour"
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  base_theme +
  theme(
    legend.key.size = unit(0.7, "cm"),
    legend.title = element_text(hjust = 0.5),
    legend.position = c(0.71, 0.26)
  )

fig5b.biomarker_roc


# =======================================================================
# fig. 5c-d — Machine-learning model development
# =======================================================================

library(mlr3verse)
library(mlr3extralearners)
library(mlr3viz)
library(mlr3mbo)
library(future)

lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

# ----------------------------------------------------------------------
# Step 0. Data Splitting and Preprocessing
# ----------------------------------------------------------------------
# task
task <- as_task_classif(zlyy_pcr_data, target = "groups", positive = "tumor")
task$set_col_roles("sample_id", roles = "name")

seed <- 1234
set.seed(seed)

# split data
split_task <- partition(task, ratio = 0.70)
trainset <- task$clone()$filter(split_task$train)
testset <- task$clone()$filter(split_task$test)

# Pre-processing 
data_prep <- po("removeconstants") %>>% po("scale")


# ----------------------------------------------------------------------
# Step 1. Learner and search spaces
# ----------------------------------------------------------------------
learners <- list(
  classif.log_reg = lrn("classif.log_reg", predict_type = "prob"),
  classif.kknn = lrn("classif.kknn", predict_type = "prob"),
  classif.rpart = lrn("classif.rpart", predict_type = "prob"),
  classif.svm = lrn(
    "classif.svm",
    type = "C-classification",
    kernel = "radial",
    predict_type = "prob"
  ),
  classif.ranger = lrn(
    "classif.ranger",
    importance = "permutation",
    predict_type = "prob"
  ),
  classif.xgboost = lrn("classif.xgboost", predict_type = "prob"),
  classif.gbm = lrn("classif.gbm", predict_type = "prob"),
  classif.nnet = lrn("classif.nnet", predict_type = "prob")
)


# Search spaces 
search_spaces <- list(
  # classif.log_reg = ps(),
  classif.kknn = ps(k = p_int(1, 10)),
  classif.rpart = ps(cp = p_dbl(0.001, 0.1), minsplit = p_int(10, 40)),
  classif.svm = ps(
    kernel = p_fct(c("radial")),
    cost = p_dbl(0.01, 10),
    gamma = p_dbl(0.001, 1)
  ),
  classif.ranger = ps(
    num.trees = p_int(200, 2000),
    mtry.ratio = p_dbl(0.1, 1.0),
    min.node.size = p_int(3, 30),
    sample.fraction = p_dbl(0.4, 0.9)
  ),
  classif.xgboost = ps(
    nrounds = p_int(50, 500),
    eta = p_dbl(0.01, 0.3),
    max_depth = p_int(1, 10)
  ),
  classif.gbm = ps(
    n.trees = p_int(50, 500),
    shrinkage = p_dbl(0.01, 0.3),
    interaction.depth = p_int(1, 10),
    n.minobsinnode = p_int(5, 25)
  ),
  classif.nnet = ps(
    size = p_int(1, 10),
    decay = p_dbl(0, 0.1)
  )
)


# ----------------------------------------------------------------------
# Step 2. Nested resampling
# ----------------------------------------------------------------------
# Inner loop : 5-fold CV, model-based (Bayesian) optimisation, 50 evaluations,
#              objective = AUC
# Outer loop : 5-fold CV, used only for performance estimation
cv_folds <- 5
n_evals <- 50

# inner_cv
inner_cv <- rsmp("cv", folds = cv_folds)

# AutoTuner
autotuners <- imap(learners, function(l, id) {
  obj <- if (id %in% names(search_spaces)) {
    auto_tuner(
      tuner = tnr("mbo"),
      learner = l$clone(),
      resampling = inner_cv,
      measure = msr("classif.auc"),
      terminator = trm("evals", n_evals = n_evals),
      search_space = search_spaces[[id]],
      store_tuning_instance = TRUE
    )
  } else {
    l$clone()
  }
  
  # Pre-processing
  glrn <- as_learner(data_prep$clone() %>>% obj)
  glrn$id <- id
  glrn
})

future::plan("multisession", workers = 12)
set.seed(seed)

outer_cv <- rsmp("cv", folds = cv_folds)
design <- benchmark_grid(trainset, autotuners, outer_cv)
bmr <- benchmark(design, store_models = TRUE)

future::plan("sequential")

model_name_map <- function(x) {
  dplyr::recode_values(
    x,
    "classif.log_reg" ~ "Logistic",
    "classif.kknn" ~ "KNN",
    "classif.rpart" ~ "Rpart",
    "classif.svm" ~ "SVM",
    "classif.ranger" ~ "Random Forest",
    "classif.xgboost" ~ "XGboost",
    "classif.gbm" ~ "GBM",
    "classif.nnet" ~ "NNET"
  )
}

model_colors <- c(
  "Logistic" = "#ef9c02",
  "KNN" = "#7363ac",
  "Rpart" = "#e3550b",
  "SVM" = "#80BA5A",
  "Random Forest" = "#b12115",
  "XGboost" = "#377EB8",
  "GBM" = "#2d9850",
  "NNET" = "#5f5f5f"
)

measures <- msrs(c(
  "classif.auc", "classif.acc", "classif.bacc", "classif.ce",
  "classif.sensitivity", "classif.specificity",
  "classif.fbeta", "classif.npv", "classif.ppv"
))


# ----------------------------------------------------------------------
# Step 3. Performance with 95% confidence intervals
# ----------------------------------------------------------------------
source("boot_ci_and_compare.R")

set.seed(seed)
ci_results <- map(seq_len(nrow(bmr$aggregate())), function(i) {
  rr <- bmr$resample_result(i)
  boot_ci_metrics(
    pred = rr$prediction(),
    model_name = model_name_map(rr$learner$id),
    positive = "tumor",
    boot_n = 1000
  )
}) %>%
  list_rbind()


# ----------------------------------------------------------------------
# fig5c. ROC comparison of the eight algorithms
# ----------------------------------------------------------------------
roc_curve_df <- map(seq_len(nrow(bmr$aggregate())), function(i) {
  rr <- bmr$resample_result(i)
  pred <- rr$prediction()
  neg <- setdiff(levels(pred$truth), "tumor")
  
  roc_obj <- pROC::roc(
    response = pred$truth,
    predictor = pred$prob[, "tumor"],
    levels = c(neg, "tumor"),
    direction = "<",
    quiet = TRUE
  )
  
  tibble(
    Models = model_name_map(rr$learner$id),
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  ) %>%
    arrange(fpr, tpr)
}) %>%
  list_rbind()

# fig5c
fig5c.model_comparison_roc <- roc_curve_df %>%
  ggplot(aes(fpr, tpr, color = Models)) +
  geom_line(linewidth = 1.4, alpha = 0.9) +
  annotate(
    geom = "segment", x = 0, xend = 1, y = 0, yend = 1,
    linetype = "dashed", linewidth = 1.4, color = "grey"
  ) +
  scale_color_manual(values = model_colors, name = "Models") +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  base_theme +
  theme(legend.key.size = unit(0.7, "cm"), legend.position = c(0.8, 0.3))

fig5c.model_comparison_roc


# ----------------------------------------------------------------------
# fig5d. Classification error across algorithms
# ----------------------------------------------------------------------
fig5d.model_resample_ce <- ci_results %>%
  filter(metric == "ce") %>%
  ggplot(aes(x = estimate, y = Models, color = Models)) +
  geom_point(size = 4) +
  geom_errorbarh(
    aes(xmin = conf_low, xmax = conf_high),
    width = 0.2,   
    linewidth = 1
  ) +
  scale_color_manual(values = model_colors) +
  labs(x = "Classification error", y = "") +
  base_theme +
  theme(legend.position = "none")

fig5d.model_resample_ce


# ----------------------------------------------------------------------
# Figure 5 combination
# ----------------------------------------------------------------------
fig5.b_d <- plot_grid(
  fig5b.biomarker_roc,
  fig5c.model_comparison_roc,
  fig5d.model_resample_ce,
  labels = c("B", "C", "D"),
  align = "h",
  ncol = 3,
  rel_widths = c(1.09, 1.01, 0.9),
  label_fontface = "bold",
  label_size = 20,
  scale = 1
)

fig5 <- plot_grid(
  fig5a.lncRNA_pcr_expr,
  fig5.b_d,
  nrow = 2,
  labels = c("A", ""),
  rel_heights = c(0.9, 1.2),
  label_fontface = "bold",
  label_size = 20,
  scale = 1
) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))

fig5

ggsave(
  "../results/p5.pcr_result_train_result.pdf",
  fig5, width = 17, height = 10, bg = "white"
)

ggsave(
  "../results/p5.pcr_result_train_result.jpg",
  fig5, width = 17, height = 10, bg = "white", dpi = 400
)


# ----------------------------------------------------------------------
# Table S4. Model performance with 95% CI
# ----------------------------------------------------------------------
model_perform_ci <- ci_results %>%
  filter(metric != "ce") %>%
  mutate(
    value = if_else(
      metric %in% c("auc", "sensitivity", "specificity"),
      label,
      sprintf("%.3f", estimate)
    )
  ) %>%
  select(Models, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value) %>% 
  relocate(sensitivity,specificity,.after = "auc")


model_perform_ci

print_table(
  model_perform_ci,
  file = "../results/TableS4.model_performance_95CI.doc"
)


# ----------------------------------------------------------------------
# Step 4. Final model fitting and selection
# ----------------------------------------------------------------------
# best tuned hyper-parameters
get_best_params <- function(m) {
  slot <- grep("\\.tuned$", names(m$model), value = TRUE)
  if (length(slot) == 0) {
    return(NULL)
  }
  m$model[[slot[1]]]$model$learner$param_set$values
}

future::plan("multisession", workers = 12)

final_models <- map(autotuners, function(obj) {
  m <- obj$clone(deep = TRUE)
  m$train(trainset)
  m
})

future::plan("sequential")

# Export best hyperparameters for each model
tuning_params <- imap(final_models, function(m, id) {
  best <- get_best_params(m)
  if (is.null(best)) {
    return(tibble(Model = id, Tuned = "No", Best_Params = "N/A"))
  }
  tibble(
    Model = id,
    Tuned = "Yes",
    Best_Params = imap_chr(best, \(v, p) paste0(p, "=", format(v, digits = 4))) %>%
      paste(collapse = ", ")
  )
}) %>%
  list_rbind()

tuning_params_res <- tuning_params %>%
  mutate(Model = model_name_map(Model))

tuning_params_res

# Table S5
print_table(
  tuning_params_res,
  file = "../results/TableS5.tuned_hyperparameters.doc"
)

# Obtain the best-performing model
best_id <- ci_results %>%
  filter(metric == "auc") %>%
  slice_max(estimate, n = 1) %>%
  pull(Models)

id_lookup <- set_names(names(final_models), model_name_map(names(final_models)))
best_learner_id <- id_lookup[[best_id]]
best_learner_id

best_learner <- final_models[[best_learner_id]]
best_learner

# save data
save(
  seed,
  trainset,
  testset,
  bmr,
  tuning_params,
  best_learner_id,
  best_learner,
  file = "final_learner.rdata"
)


# ----------------------------------------------------------------------
# Figure S1. Performance heatmap and AUC forest plot
# ----------------------------------------------------------------------

# fig s1.a:  Heatmap of Model Performance
model_perform_data <- ci_results %>%
  filter(metric != "ce") %>%
  mutate(metric = substr(toupper(metric), 1, 4))

figs1.a <- model_perform_data %>%
  ggplot(aes(
    x = reorder(metric, estimate),
    y = reorder(Models, estimate),
    fill = estimate
  )) +
  geom_tile(color = "black") +
  geom_text(aes(label = sprintf("%.3f", estimate)), color = "white", size = 5) +
  scale_fill_viridis_c(direction = -1) +
  labs(title = "Model performance (nested 5-fold CV)", x = "Metrics", y = "") +
  base_theme +
  theme(legend.position = "none")

figs1.a

# fig s1.b: AUC forest plot with 95% CIs
figs1.b <- ci_results %>%
  filter(metric == "auc") %>%
  ggplot(aes(x = estimate, y = reorder(Models, estimate), color = Models)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.2, linewidth = 1) +
  scale_color_manual(values = model_colors) +
  labs(x = "AUC (95% CI)", y = "", title = "AUC forest plot") +
  base_theme +
  theme(legend.position = "none")

figs1.b

fig.s1 <- plot_grid(
  figs1.a,
  figs1.b,
  labels = LETTERS[1:2],
  ncol = 2,
  label_fontface = "bold",
  label_size = 20,
  scale = 1,
  rel_widths = c(0.55, 0.45)
)

fig.s1

ggsave("../results/ps1.model_perform_cv.pdf", fig.s1, width = 12, height = 5)
ggsave("../results/ps1.model_perform_cv.jpg", fig.s1, width = 12, height = 5, dpi = 350)


# ----------------------------------------------------------------------
# Figure S2. Permutation variable importance of the selected model
# ----------------------------------------------------------------------
imp <- best_learner$importance()
importance_df <- as.data.table(imp, keep.rownames = "feature")

figs2.rf_importance <- importance_df %>%
  ggplot(aes(x = reorder(feature, imp), y = imp)) +
  geom_col(fill = "#457b9d") +
  coord_flip() +
  labs(x = "Feature", y = "Permutation importance", title = "Random forest") +
  base_theme

figs2.rf_importance

ggsave("../results/ps2.rf_importance.pdf", figs2.rf_importance, width = 6, height = 6)
ggsave("../results/ps2.rf_importance.jpg", figs2.rf_importance, width = 6, height = 6, dpi = 350)
