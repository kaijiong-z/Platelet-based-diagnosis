# ==================================================================
# boot_ci_metrics()  单组性能 + 95% CI
# boot_compare_models()  同一批患者上的头对头配对比较
# ==================================================================
# 依赖：tidyverse、pROC
# ==================================================================

#' 单组预测的性能指标与 95% CI
#'
#' @param pred mlr3 PredictionClassif（建议用 rr$prediction() 的池化折外预测）
#' @param model_name 显示名
#' @param positive 阳性类标签
#' @param boot_n 自助次数；0 表示跳过自助（CI 置 NA）
#' @param alpha 显著性水平，0.05 对应 95% CI
#' @param auc_ci AUC 的 CI 方法："delong"（默认，解析法）或 "bootstrap"
#' @param seed 随机种子。NULL（默认）= 不干预全局随机流；
#'   给定数值时函数内固定种子（仅在需要单独复现某一次调用时使用）
#' @param boot_idx 外部传入的自助索引 list（配对比较用）。给定时忽略
#'   seed 与内部抽样，直接使用这批索引 —— 这样多个模型共享同一批
#'   重抽样患者，差值分布才能反映配对结构
#' @param return_boot 是否把自助矩阵挂在 attr "boot_mat" 上（配对比较需要）

library(tidyverse);library(pROC)
boot_ci_metrics <- function(
  pred,
  model_name,
  positive = "tumor",
  boot_n = 1000,
  alpha = 0.05,
  auc_ci = c("delong", "bootstrap"),
  seed = NULL,
  boot_idx = NULL,
  return_boot = FALSE
) {
  auc_ci <- match.arg(auc_ci)

  stopifnot(inherits(pred, "PredictionClassif"))
  if (is.null(pred$prob)) {
    stop("Prediction 不含概率，请设置 predict_type = 'prob'")
  }
  truth <- pred$truth
  if (!positive %in% levels(truth)) {
    stop(sprintf(
      "positive = '%s' 不在 truth 水平中：%s",
      positive,
      paste(levels(truth), collapse = ", ")
    ))
  }

  prob_pos <- pred$prob[, positive]
  response <- pred$response
  n <- length(truth)
  negative <- setdiff(levels(truth), positive)
  stopifnot(length(negative) == 1L)

  # ---- 单次计算全部指标（含除零守卫，避免 NaN 被 na.rm 静默丢弃）----
  calc_all <- function(idx) {
    t_i <- truth[idx]
    p_i <- prob_pos[idx]
    r_i <- response[idx]

    auc_i <- if (length(unique(t_i)) < 2L) {
      NA_real_
    } else {
      as.numeric(pROC::auc(pROC::roc(
        response = t_i,
        predictor = p_i,
        levels = c(negative, positive),
        direction = "<",
        quiet = TRUE
      )))
    }

    tp <- sum(t_i == positive & r_i == positive)
    tn <- sum(t_i != positive & r_i != positive)
    fp <- sum(t_i != positive & r_i == positive)
    fn <- sum(t_i == positive & r_i != positive)

    sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
    ppv <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    npv <- if (tn + fn > 0) tn / (tn + fn) else NA_real_
    acc <- (tp + tn) / length(t_i)
    f1 <- if (is.na(ppv) || is.na(sens) || ppv + sens == 0) {
      NA_real_
    } else {
      2 * ppv * sens / (ppv + sens)
    }

    c(
      auc = auc_i,
      acc = acc,
      bacc = (sens + spec) / 2,
      ce = 1 - acc,
      sensitivity = sens,
      specificity = spec,
      f1 = f1,
      npv = npv,
      ppv = ppv
    )
  }

  est <- calc_all(seq_len(n))

  # ---- 自助分布 ----
  boot_mat <- NULL
  if (!is.null(boot_idx)) {
    # 外部索引：与其它模型共享同一批重抽样患者
    boot_mat <- vapply(boot_idx, calc_all, numeric(length(est)))
  } else if (boot_n > 0) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    idx_pos <- which(truth == positive)
    idx_neg <- which(truth != positive)
    boot_mat <- vapply(
      seq_len(boot_n),
      function(b) {
        idx <- c(
          sample(idx_pos, length(idx_pos), replace = TRUE),
          sample(idx_neg, length(idx_neg), replace = TRUE)
        )
        calc_all(idx)
      },
      numeric(length(est))
    )
  }

  if (!is.null(boot_mat)) {
    lower <- apply(boot_mat, 1, quantile, probs = alpha / 2, na.rm = TRUE)
    upper <- apply(boot_mat, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  } else {
    lower <- setNames(rep(NA_real_, length(est)), names(est))
    upper <- setNames(rep(NA_real_, length(est)), names(est))
  }

  # ---- AUC 改用 DeLong 解析区间 ----
  if (auc_ci == "delong") {
    ca <- tryCatch(
      {
        ro <- pROC::roc(
          response = truth,
          predictor = prob_pos,
          levels = c(negative, positive),
          direction = "<",
          quiet = TRUE
        )
        as.numeric(pROC::ci.auc(ro, method = "delong", conf.level = 1 - alpha))
      },
      error = function(e) c(NA_real_, NA_real_, NA_real_)
    )
    lower[["auc"]] <- ca[1]
    upper[["auc"]] <- ca[3]
  }

  out <- tibble(
    Models = model_name,
    metric = names(est),
    estimate = as.numeric(est),
    conf_low = as.numeric(lower),
    conf_high = as.numeric(upper)
  ) %>%
    mutate(
      label = if_else(
        is.na(conf_low),
        sprintf("%.3f", estimate),
        sprintf("%.3f (%.3f-%.3f)", estimate, conf_low, conf_high)
      )
    ) 

  if (return_boot) {
    attr(out, "boot_mat") <- boot_mat
  }
  out
}


#' 同一批患者上多个模型的头对头配对比较
#'
#' 每次自助重抽样抽的是"患者"，同一批索引同时用于所有模型，
#' 因此差值分布保留了配对结构。这是正确的比较方式 ——
#' 直接比较两条独立 CI 是否重叠会低估显著性（重叠仍可能显著）。
#'
#' AUC 额外给出 DeLong 配对检验的 p 值（pROC::roc.test），
#' 这是临床文献比较相关 ROC 曲线的标准方法。
#'
#' @param preds named list of PredictionClassif；必须是同一批患者、同一顺序
#' @param reference 作为参照的组名（其余组与它相减）
#' @param positive 阳性类标签
#' @param boot_n 自助次数
#' @param alpha 显著性水平
#' @param seed 随机种子（索引只生成一次，所有模型共享）
#' @return list(metrics = 各组性能表, comparison = 与参照的差值表)
boot_compare_models <- function(
  preds,
  reference,
  positive = "tumor",
  boot_n = 1000,
  alpha = 0.05,
  seed = 42
) {
  stopifnot(is.list(preds), !is.null(names(preds)))
  group_names <- names(preds)
  if (!reference %in% group_names) {
    stop(sprintf("reference = '%s' 不在 preds 的名字中", reference))
  }

  # ---- 配对前提：同一批患者、同一顺序 ----
  truths <- map(preds, ~ .x$truth)
  n_all <- map_int(truths, length)
  if (length(unique(n_all)) > 1L) {
    stop(sprintf(
      "各组样本量不同（%s），无法做配对比较。\n请先取各组共同的完整病例再比较。",
      paste(sprintf("%s=%d", group_names, n_all), collapse = ", ")
    ))
  }
  if (
    !all(map_lgl(
      truths[-1],
      ~ identical(as.character(.x), as.character(truths[[1]]))
    ))
  ) {
    stop("各组的 truth 不一致，说明不是同一批患者或顺序不同，无法配对比较")
  }

  # ---- 索引只生成一次，所有模型共享 ----
  truth1 <- truths[[1]]
  idx_pos <- which(truth1 == positive)
  idx_neg <- which(truth1 != positive)
  set.seed(seed)
  boot_idx <- map(seq_len(boot_n), function(b) {
    c(
      sample(idx_pos, length(idx_pos), replace = TRUE),
      sample(idx_neg, length(idx_neg), replace = TRUE)
    )
  })

  res <- imap(preds, function(p, g) {
    boot_ci_metrics(
      p,
      g,
      positive = positive,
      alpha = alpha,
      boot_idx = boot_idx,
      return_boot = TRUE
    )
  })

  metrics_tbl <- map(res, \(x) x) |> list_rbind()

  # ---- 差值分布：Δ = 各组 - 参照 ----
  ref_mat <- attr(res[[reference]], "boot_mat")
  metric_names <- rownames(ref_mat)

  comparison <- imap(res[setdiff(group_names, reference)], function(r, g) {
    d_mat <- attr(r, "boot_mat") - ref_mat # 逐次自助的差值
    est_diff <- r$estimate - res[[reference]]$estimate

    tibble(
      Comparison = sprintf("%s - %s", g, reference),
      metric = metric_names,
      diff = est_diff,
      diff_low = apply(d_mat, 1, quantile, probs = alpha / 2, na.rm = TRUE),
      diff_high = apply(d_mat, 1, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
    ) |>
      mutate(
        significant = !(diff_low <= 0 & diff_high >= 0), # CI 不含 0
        label = sprintf("%+.3f (%+.3f, %+.3f)", diff, diff_low, diff_high)
      )
  }) |>
    list_rbind()

  # ---- AUC 的 DeLong 配对检验 ----
  neg <- setdiff(levels(truth1), positive)
  roc_of <- function(p) {
    pROC::roc(
      response = p$truth,
      predictor = p$prob[, positive],
      levels = c(neg, positive),
      direction = "<",
      quiet = TRUE
    )
  }
  roc_ref <- roc_of(preds[[reference]])

  delong_p <- imap(preds[setdiff(group_names, reference)], function(p, g) {
    tt <- tryCatch(
      pROC::roc.test(roc_of(p), roc_ref, method = "delong", paired = TRUE),
      error = function(e) NULL
    )
    tibble(
      Comparison = sprintf("%s - %s", g, reference),
      metric = "auc",
      delong_p = if (is.null(tt)) NA_real_ else unname(tt$p.value)
    )
  }) |>
    list_rbind()

  comparison <- comparison |>
    left_join(delong_p, by = c("Comparison", "metric"))

  list(metrics = metrics_tbl, comparison = comparison)
}

# ==================================================================
# 用法
# ==================================================================
## ---- A. 三队列（不同患者，不能配对）：各自算 CI ----
## preds_cohort <- list(
##   Training = bmr$resample_result(idx)$prediction(),
##   Testing = best_learner$predict(testset),
##   `External validation` = best_learner$predict(validation_task)
## )
## set.seed(seed)   # 在 map 外设一次
## cohort_ci <- imap(preds_cohort, \(p, g) boot_ci_metrics(p, g)) |>
##   list_rbind() |>
##   rename(Cohort = Models)
##
## cohort_ci |>
##   select(Cohort, metric, label) |>
##   pivot_wider(names_from = metric, values_from = label)
##
## ---- B. 同一队列三个标志物（同一批患者）：头对头配对比较 ----
## preds_marker <- list(
##   lncrna = pred_lncrna,
##   ca153 = pred_ca153,
##   cea = pred_cea
## )
## cmp <- boot_compare_models(preds_marker, reference = "lncrna")
##
## cmp$metrics       # 各标志物性能 + CI
## cmp$comparison    # Δ 与其 95% CI；significant = CI 不含 0
##                   # auc 行另有 delong_p（DeLong 配对检验 p 值）
##
## # 只看 AUC 的比较
## cmp$comparison |> filter(metric == "auc")
##
## ---- 注意 ----
## 配对比较要求三组是同一批患者、同一顺序。若 CA153/CEA 曾各自
## drop_na() 导致 N 不同，函数会直接报错；需先取三者共同的完整病例，
## 用同一批样本重新生成预测后再比较。
