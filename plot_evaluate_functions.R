# ==================================================================
# 多组临床评价绘图函数：plot_roc_curves() / plot_dca_curves() / plot_cal_curves()
# ==================================================================
# 输入统一为"预测 list"：named list of mlr3 PredictionClassif，
#   名字即分组名（队列/模型/亚组均可，决定图例与面板标题及顺序），如
#   preds <- list(
#     Training  = bmr$resample_result(idx)$prediction(),  # 折外预测
#     Testing   = best_learner$predict(testset),
#     `External validation` = best_learner$predict(validation_task)
#   )
# 两个函数均支持：
#   type = "overlay"  多组叠加一张图
#   type = "panel"    每队列一个面板（patchwork 拼接）
#   type = "both"     返回 list(overlay = ..., panel = ...)
# 依赖：tidyverse, pROC, dcurves, patchwork；ggplot 主题可经 theme 参数传入
# ==================================================================

library(tidyverse)
library(pROC)
library(dcurves)
library(patchwork)
library(cowplot)
# ------------------------------------------------------------------
# 1. ROC
# ------------------------------------------------------------------

#' 多组 ROC 绘图（叠加 / 分面板）
#'
#' @param preds named list of PredictionClassif；名字作分组名
#' @param positive 阳性类标签
#' @param type "overlay" 叠加一张图 / "panel" 每组一个面板
#' @param colors 分组配色；NULL 时自动取内置色板
#' @param panel_ncol 分面板列数（assemble = TRUE 时生效）
#' @param assemble panel 模式下是否 wrap_plots 拼好；FALSE 返回单图 list 供自行拼装
#' @param theme_set 追加到每张图的 ggplot 主题
#' @return ggplot；panel 且 assemble=FALSE 时返回 named list of ggplot
plot_roc_curves <- function(
  preds,
  positive = "tumor",
  type = c("overlay", "panel"),
  colors = NULL,
  panel_ncol = 3,
  assemble = TRUE,
  theme_set = ggplot2::theme_bw()
) {
  type <- match.arg(type)
  stopifnot(is.list(preds), !is.null(names(preds)), all(nzchar(names(preds))))

  group_names <- names(preds)
  if (is.null(colors)) {
    palette_default <- c(
      "#b12115",
      "#377EB8",
      "#2d9850",
      "#ef9c02",
      "#7363ac",
      "#5f5f5f"
    )
    colors <- set_names(palette_default[seq_along(group_names)], group_names)
  }

  # 用户若传入未命名颜色向量，按组名命名，
  # 使 panel 分支的 colors[[group]] 与 overlay 的 scale_color_manual 都可用
  if (is.null(names(colors)) && length(colors) == length(group_names)) {
    colors <- set_names(colors, group_names)
  }

  ## roc 对象：levels 显式指定
  roc_objs <- imap(preds, function(pred, group) {
    neg <- setdiff(levels(pred$truth), positive)
    stopifnot(length(neg) == 1L)
    pROC::roc(
      response = pred$truth,
      predictor = pred$prob[, positive],
      levels = c(neg, positive),
      direction = "<",
      quiet = TRUE
    )
  })

  ## AUC (DeLong 95% CI)
  auc_ci <- imap(roc_objs, function(r, group) {
    v <- as.numeric(pROC::ci.auc(r))
    tibble(Group = group, auc = v[2], ci_low = v[1], ci_high = v[3])
  }) |>
    list_rbind()

  roc_df <- imap(roc_objs, function(r, group) {
    tibble(Group = group, fpr = 1 - r$specificities, tpr = r$sensitivities) |>
      arrange(fpr, tpr)
  }) |>
    list_rbind() |>
    mutate(Group = factor(Group, levels = group_names))

  ## ---- overlay ----
  make_overlay <- function() {
    legend_labels <- auc_ci |>
      mutate(
        lab = sprintf("%s: %.3f (%.3f-%.3f)", Group, auc, ci_low, ci_high)
      ) |>
      pull(lab) |>
      set_names(auc_ci$Group)

    ggplot(roc_df, aes(fpr, tpr, color = Group)) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "grey70",
        linewidth = 1
      ) +
      geom_line(linewidth = 1.3, alpha = 0.9) +
      scale_color_manual(
        values = colors,
        labels = legend_labels,
        name = "AUC (95% CI)"
      ) +
      labs(x = "1 - Specificity", y = "Sensitivity") +
      # coord_fixed() +
      theme_set +
      theme(legend.position = c(0.62, 0.18))
  }

  ## ---- panel ----
  make_panel <- function() {
    plots <- imap(roc_objs, function(r, group) {
      ci <- auc_ci |> filter(Group == group)
      tibble(fpr = 1 - r$specificities, tpr = r$sensitivities) |>
        arrange(fpr, tpr) |>
        ggplot(aes(fpr, tpr)) +
        geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          color = "grey70"
        ) +
        geom_line(linewidth = 1.3, color = colors[[group]]) +
        annotate(
          "text",
          x = 0.50,
          y = 0.12,
          hjust = 0,
          size = 4,
          label = sprintf("AUC %.3f (%.3f-%.3f)", ci$auc, ci$ci_low, ci$ci_high)
        ) +
        labs(x = "1 - Specificity", y = "Sensitivity", title = group) +
        # coord_fixed() +
        theme_set
    })
    if (assemble) {
      patchwork::wrap_plots(plotlist = plots, ncol = panel_ncol)
    } else {
      plots
    }
  }

  switch(
    type,
    overlay = make_overlay(),
    panel = make_panel()
  )
}

##----- example ------

# # pred list
# preds <- list(
#   lncrna = train_lncrna_pred,
#   ca153 = train_ca153_pred,
#   cea = train_cea_pred
# )

# # overlay
# p_roc <- plot_roc_curves(preds, type = "overlay")
#
# # panel
# p_roc_panel <- plot_roc_curves(preds, type = "panel")

# #single plot list
# roc_list <- plot_roc_curves(preds, type = "panel", assemble = FALSE)

# ------------------------------------------------------------------
# 2. DCA
# ------------------------------------------------------------------

#' 多组 DCA 绘图（叠加 / 分面板）
#'
#' @param preds named list of PredictionClassif；名字作分组名
#' @param positive 阳性类标签
#' @param type "overlay" 叠加一张图 / "panel" 每组一个面板
#' @param model_label 模型显示名。NULL（默认）= 逐组使用分组名
#'   （多模型对比场景）；单个字符串 = 所有组共用同一名称
#'   （多队列验证同一模型的场景）；也可传命名向量按组指定
#' @param thresholds DCA 阈值网格
#' @param min_run 连续区段最少阈值点数；短于此的碎段按噪声丢弃（默认 3，即 thresholds 步长 0.01 下约 0.03 宽）
#' @param colors 分组配色（overlay 用）；NULL 自动取色板
#' @param panel_ncol 分面板列数（assemble = TRUE 时生效）
#' @param assemble panel 模式下是否 wrap_plots 拼好；FALSE 返回单图 list 供自行拼装
#' @param theme_set 追加主题
#' @return ggplot（或 panel+assemble=FALSE 时 named list）；附 attr "benefit_range"（每组主区间）与 "benefit_segments"（全部连续区段）
plot_dca_curves <- function(
  preds,
  positive = "tumor",
  type = c("overlay", "panel"),
  model_label = NULL,
  thresholds = seq(0, 0.99, by = 0.01),
  min_run = 3,
  colors = NULL,
  panel_ncol = 3,
  assemble = TRUE,
  theme_set = ggplot2::theme_bw()
) {
  type <- match.arg(type)
  stopifnot(is.list(preds), !is.null(names(preds)), all(nzchar(names(preds))))

  group_names <- names(preds)
  if (is.null(colors)) {
    palette_default <- c(
      "#b12115",
      "#377EB8",
      "#2d9850",
      "#ef9c02",
      "#7363ac",
      "#5f5f5f"
    )
    colors <- set_names(palette_default[seq_along(group_names)], group_names)
  }

  ## 逐组模型显示名
  model_labs <- if (is.null(model_label)) {
    set_names(group_names, group_names) # 多模型场景：面板各用组名
  } else if (length(model_label) == 1L) {
    set_names(rep(model_label, length(group_names)), group_names) # 多队列场景
  } else {
    stopifnot(all(group_names %in% names(model_label)))
    model_label[group_names]
  }

  ## 逐组 DCA
  dca_objs <- imap(preds, function(pred, group) {
    d <- tibble(
      outcome = as.integer(pred$truth == positive),
      prob = pred$prob[, positive]
    )
    dcurves::dca(
      outcome ~ prob,
      data = d,
      thresholds = thresholds,
      label = list(prob = model_labs[[group]])
    )$dca
  })

  ## 获益阈值区间：模型净获益 > max(All, None) 且 > 0。
  ## 满足条件的阈值可能不连续（尾部数据稀疏常产生孤立尖峰），
  ## 因此识别"连续区段"，主区间取最长区段；min_run 个阈值点以下的
  ## 碎段视为噪声丢弃。所有保留区段附在 attr("benefit_segments")。
  benefit_segments <- imap(dca_objs, function(d, group) {
    wide <- d |>
      select(label, threshold, net_benefit) |>
      pivot_wider(names_from = label, values_from = net_benefit) |>
      arrange(threshold) |>
      mutate(
        ok = .data[[model_labs[[group]]]] >
          pmax(`Treat All`, `Treat None`, na.rm = TRUE) &
          .data[[model_labs[[group]]]] > 0,
        ok = replace_na(ok, FALSE),
        run_id = consecutive_id(ok) # 相邻同值成段
      )
    wide |>
      filter(ok) |>
      summarise(
        range_low = min(threshold),
        range_high = max(threshold),
        n_points = n(),
        .by = run_id
      ) |>
      mutate(Group = group, .before = 1) |>
      select(-run_id) |>
      filter(n_points >= min_run) # 丢弃噪声碎段
  }) |>
    list_rbind()

  ## 主区间 = 每组最长的连续区段
  benefit_range <- benefit_segments |>
    slice_max(n_points, n = 1, with_ties = FALSE, by = Group) |>
    select(Group, range_low, range_high)

  ## 保证每组都有一行（可能全被过滤）
  benefit_range <- tibble(Group = group_names) |>
    left_join(benefit_range, by = "Group")

  model_df <- imap(dca_objs, function(d, group) {
    d |>
      filter(!label %in% c("Treat All", "Treat None")) |>
      transmute(Group = group, threshold, net_benefit)
  }) |>
    list_rbind() |>
    mutate(Group = factor(Group, levels = group_names))

  all_df <- imap(dca_objs, function(d, group) {
    d |>
      filter(label == "Treat All") |>
      transmute(Group = group, threshold, net_benefit)
  }) |>
    list_rbind() |>
    mutate(Group = factor(Group, levels = group_names))

  ## 各组 treat-all 是否完全一致（同一队列上比较多个模型时成立）？
  ## 一致则叠画只会显示最后一组的颜色，误导读者——合并为一条灰色虚线。
  all_identical <- all_df |>
    summarise(
      n_distinct_nb = n_distinct(round(net_benefit, 10)),
      .by = threshold
    ) |>
    summarise(all(n_distinct_nb == 1)) |>
    pull()

  y_max <- max(model_df$net_benefit, na.rm = TRUE) * 1.08

  ## ---- overlay ----
  make_overlay <- function() {
    legend_labels <- benefit_range |>
      mutate(
        lab = if_else(
          is.na(range_low),
          paste0(Group, ": none"),
          sprintf("%s: %.2f-%.2f", Group, range_low, range_high)
        )
      ) |>
      pull(lab) |>
      set_names(benefit_range$Group)

    ggplot() +
      geom_hline(
        aes(yintercept = 0, linetype = "None"),
        color = "black",
        linewidth = 1.1
      ) +
      {
        if (all_identical) {
          # 同一队列多模型：一条灰色 treat-all
          geom_line(
            data = all_df |> filter(Group == group_names[1]),
            aes(threshold, net_benefit, linetype = "All"),
            color = "grey55",
            linewidth = 1.1
          )
        } else {
          # 多队列：各组一条同色 treat-all（患病率不同，不能合并）
          geom_line(
            data = all_df,
            aes(threshold, net_benefit, color = Group, linetype = "All"),
            linewidth = 0.7,
            alpha = 0.55
          )
        }
      } +
      geom_line(
        data = model_df,
        aes(threshold, net_benefit, color = Group),
        linewidth = 1.2
      ) +
      scale_color_manual(
        values = colors,
        labels = legend_labels,
        name = "Benefit range"
      ) +
      scale_linetype_manual(
        values = c(All = "dashed", None = "dotted"),
        name = NULL
      ) +
      guides(
        color = guide_legend(order = 1),
        # 隐藏 All / None 参考线的图例（线条仍绘制，仅不在图上显示图例）
        linetype = "none"
        # linetype = guide_legend(
        #   order = 2,
        #   override.aes = list(color = c("grey55", "black"))
        # )
      ) +
      coord_cartesian(ylim = c(-0.05, y_max)) +
      labs(x = "Threshold Probability", y = "Net Benefit") +
      theme_set +
      theme(
        legend.position = c(0.80, 0.75),
        # 两组图例（color 与 linetype）无法真正合并——映射不同；
        # 收紧间距、去掉边距，使其视觉上连为一个分段图例
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt")
      )
  }

  ## ---- panel ----
  make_panel <- function() {
    plots <- imap(dca_objs, function(d, group) {
      rng <- benefit_range |> filter(Group == group)
      m_lab <- model_labs[[group]]
      rng_lab <- if (!is.na(rng$range_low)) {
        sprintf(
          "Benefit range:\n%s: %.2f - %.2f",
          m_lab,
          rng$range_low,
          rng$range_high
        )
      } else {
        "Benefit range: none"
      }

      d |>
        # 注意：recode_values 对未匹配值置 NA（会静默丢掉模型曲线），
        # 这里用 if_else 链做无歧义的选择性替换
        mutate(
          label = as.character(label),
          label = if_else(label == "Treat All", "All", label),
          label = if_else(label == "Treat None", "None", label),
          # 因子水平即图例顺序：模型在前，All/None 在后
          label = factor(label, levels = c(m_lab, "All", "None"))
        ) |>
        ggplot(aes(threshold, net_benefit, color = label, linetype = label)) +
        geom_line(linewidth = 1.1) +
        # scale_color_manual(
        #   values = set_names(
        #     c(colors[[group]], "grey60", "black"),
        #     c(m_lab, "All", "None")
        #   )
        # ) +
        scale_linetype_manual(
          values = set_names(
            c("solid", "dashed", "dotted"),
            c(m_lab, "All", "None")
          )
        ) +
        coord_cartesian(ylim = c(-0.05, y_max)) +
        annotate(
          "text",
          x = 0.02,
          y = y_max * 0.15,
          hjust = 0,
          size = 3.6,
          label = rng_lab
        ) +
        labs(
          x = "Threshold Probability",
          y = "Net Benefit",
          title = group,
          color = NULL,
          linetype = NULL
        ) +
        theme_set +
        theme(legend.position = c(0.8, 0.8))
    })
    if (assemble) {
      patchwork::wrap_plots(plots, ncol = panel_ncol) # 不加 tag，留给调用侧
    } else {
      plots # named list of ggplot，自行拼装
    }
  }

  out <- switch(
    type,
    overlay = make_overlay(),
    panel = make_panel()
  )
  attr(out, "benefit_range") <- benefit_range
  attr(out, "benefit_segments") <- benefit_segments
  out
}


##----- example ------

# # pred list
# preds <- list(
#   lncrna = train_lncrna_pred,
#   ca153 = train_ca153_pred,
#   cea = train_cea_pred
# )

# # overlay
# p_dca <- plot_dca_curves(preds, type = "overlay")
#
# # panel
# p_dca_panel <- plot_dca_curves(preds, type = "panel")

# #single plot list
# dca_list <- plot_dca_curves(preds, type = "panel", assemble = FALSE)

## # 取 DCA 获益区间表（主区间 = 最长连续区段）
## attr(p_dca, "benefit_range")
## # 全部连续区段（诊断尾部噪声段用）
## attr(p_dca, "benefit_segments")

# ------------------------------------------------------------------
# 3. calibration curve
# ------------------------------------------------------------------

#' 多组校准曲线（叠加 / 分面板）
#'
#' 校准曲线，两种样式可选：
#'   method = "bins"  等频分箱的（预测均值, 实际发生率）点（空心圆）直线相连
#'   method = "loess" LOESS 平滑曲线（panel 模式带置信带）
#' 均注记 Brier 分数、校准截距（calibration-in-the-large）与校准斜率。
#' 截距/斜率来自 logistic 再校准模型：
#'   glm(outcome ~ offset(logit_p))  -> 截距（理想 0）
#'   glm(outcome ~ logit_p)          -> 斜率（理想 1；<1 过拟合特征）
#'
#' @param preds named list of PredictionClassif；名字作分组名
#' @param positive 阳性类标签
#' @param type "overlay" 叠加一张图 / "panel" 每组一个面板
#' @param method "bins" 分箱折线（默认）/ "loess" 平滑曲线
#' @param n_bins 分箱箱数（method = "bins" 时生效；等频分箱）
#' @param error_bars 分箱点是否加 95% Wilson 区间误差棒
#' @param show_bins method = "loess" 时是否在平滑曲线上叠加分箱点
#'   （method = "bins" 时分箱点本就是主体，此参数不生效）。
#'   LOESS 本身是连续曲线、不产生离散点，其不确定性由置信带表达；
#'   叠加的分箱点及其 Wilson 误差棒提供经验证据与样本量感
#' @param loess_span LOESS 带宽（method = "loess" 时生效；小样本可调大）
#' @param colors 分组配色；NULL 自动取色板
#' @param panel_ncol 分面板列数（assemble = TRUE 时生效）
#' @param assemble panel 模式下是否 wrap_plots 拼好；FALSE 返回单图 list
#' @param theme_set 追加主题
#' @return ggplot（或 panel+assemble=FALSE 时 named list）；
#'   附 attr "calibration_metrics"（Brier / 截距 / 斜率表）

plot_cal_curves <- function(
  preds,
  positive = "tumor",
  type = c("overlay", "panel"),
  method = c("bins", "loess"),
  n_bins = 5,
  error_bars = TRUE,
  show_bins = TRUE,
  loess_span = 0.75,
  colors = NULL,
  panel_ncol = 3,
  assemble = TRUE,
  theme_set = ggplot2::theme_bw()
) {
  type <- match.arg(type)
  method <- match.arg(method)
  stopifnot(is.list(preds), !is.null(names(preds)), all(nzchar(names(preds))))

  group_names <- names(preds)
  if (is.null(colors)) {
    palette_default <- c(
      "#b12115",
      "#377EB8",
      "#2d9850",
      "#ef9c02",
      "#7363ac",
      "#5f5f5f"
    )
    colors <- set_names(palette_default[seq_along(group_names)], group_names)
  }

  # 用户若传入未命名颜色向量（如 c("#b12115", ...)），按组名命名，
  # 使 panel 分支的 colors[[group]] 与 overlay 的 scale_*_manual 都可用
  if (is.null(names(colors)) && length(colors) == length(group_names)) {
    colors <- set_names(colors, group_names)
  }

  ## 逐组整理 outcome + 预测概率（裁掉 0/1 端点避免 logit 无穷）
  cal_data <- imap(preds, function(pred, group) {
    tibble(
      Group = group,
      outcome = as.integer(pred$truth == positive),
      prob = pmin(pmax(pred$prob[, positive], 1e-6), 1 - 1e-6)
    )
  })

  ## 量化指标：Brier + 截距 + 斜率
  calibration_metrics <- imap(cal_data, function(d, group) {
    logit_p <- qlogis(d$prob)
    fit_int <- glm(outcome ~ offset(logit_p), family = binomial, data = d)
    fit_slp <- glm(outcome ~ logit_p, family = binomial, data = d)
    tibble(
      Group = group,
      brier = mean((d$prob - d$outcome)^2),
      intercept = unname(coef(fit_int)[1]),
      slope = unname(coef(fit_slp)["logit_p"])
    )
  }) |>
    list_rbind()

  ## 分箱观测点（等频分箱：每箱样本量近似相等，小样本下比等宽稳）
  smooth_df <- list_rbind(cal_data) |>
    mutate(Group = factor(Group, levels = group_names))

  ## 每箱实际发生率为二项比例；95% CI 用 Wilson 区间
  ## （小箱样本量下比正态近似稳，且不会越出 [0,1]）
  wilson_ci <- function(x, n, conf = 0.95) {
    z <- qnorm(1 - (1 - conf) / 2)
    p <- x / n
    denom <- 1 + z^2 / n
    center <- (p + z^2 / (2 * n)) / denom
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
    tibble(ci_low = pmax(0, center - half), ci_high = pmin(1, center + half))
  }

  bin_df <- imap(cal_data, function(d, group) {
    d |>
      mutate(bin = ntile(prob, n_bins)) |>
      summarise(
        pred_mean = mean(prob),
        obs_rate = mean(outcome),
        events = sum(outcome),
        n = n(),
        .by = bin
      ) |>
      (\(b) bind_cols(b, wilson_ci(b$events, b$n)))() |>
      mutate(Group = group, .before = 1)
  }) |>
    list_rbind() |>
    mutate(Group = factor(Group, levels = group_names))

  ## ---- overlay ----
  make_overlay <- function() {
    legend_labels <-
      # calibration_metrics |>
      # mutate(lab = sprintf(
      #   "%s: Brier %.3f, Int %.2f, Slope %.2f",
      #   Group, brier, intercept, slope
      # )) |>
      # pull(lab) |>
      set_names(calibration_metrics$Group)

    p <- ggplot() +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "grey60",
        linewidth = 1
      )

    # bins
    p <- if (method == "bins") {
      p <- p +
        geom_line(
          data = bin_df,
          aes(pred_mean, obs_rate, color = Group),
          linewidth = 1.2
        )
      if (error_bars) {
        p <- p +
          geom_errorbar(
            data = bin_df,
            aes(pred_mean, ymin = ci_low, ymax = ci_high, color = Group),
            width = 0.02,
            linewidth = 0.6,
            alpha = 0.8
          )
      }
      p +
        geom_point(
          data = bin_df,
          aes(pred_mean, obs_rate, color = Group),
          shape = 1,
          size = 3,
          stroke = 1.2 # 空心圆
        )
    } else {
      # loess
      p <- p +
        geom_smooth(
          data = smooth_df,
          aes(prob, outcome, color = Group, fill = Group),
          method = "loess",
          span = loess_span,
          se = TRUE,
          alpha = 0.15,
          linewidth = 1.2
        )
      if (show_bins) {
        if (error_bars) {
          p <- p +
            geom_errorbar(
              data = bin_df,
              aes(pred_mean, ymin = ci_low, ymax = ci_high, color = Group),
              width = 0.02,
              linewidth = 0.6,
              alpha = 0.8
            )
        }
        p <- p +
          geom_point(
            data = bin_df,
            aes(pred_mean, obs_rate, color = Group),
            shape = 1,
            size = 3,
            stroke = 1.2
          )
      }
      p
    }

    p <- p +
      scale_color_manual(
        values = colors, #labels = legend_labels,
        name = NULL
      )

    if (method == "loess") {
      p <- p +
        scale_fill_manual(
          values = colors, #labels = legend_labels,
          name = NULL
        )
    }

    p +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "Predicted probability", y = "Actual probability") +
      theme_set +
      theme(legend.position = c(0.68, 0.14))
  }

  ## ---- panel ----
  make_panel <- function() {
    plots <- imap(group_names, function(group, i) {
      met <- calibration_metrics |> filter(Group == group)
      d_b <- bin_df |> filter(Group == group)

      d_s <- smooth_df |> filter(Group == group)

      p <- ggplot() +
        geom_abline(
          slope = 1,
          intercept = 0,
          linetype = "dashed",
          color = "grey60"
        )

      p <- if (method == "bins") {
        p <- p +
          geom_line(
            data = d_b,
            aes(pred_mean, obs_rate),
            color = colors[[group]],
            linewidth = 1.2
          )
        if (error_bars) {
          p <- p +
            geom_errorbar(
              data = d_b,
              aes(pred_mean, ymin = ci_low, ymax = ci_high),
              color = colors[[group]],
              width = 0.02,
              linewidth = 0.6,
              alpha = 0.8
            )
        }
        p +
          geom_point(
            data = d_b,
            aes(pred_mean, obs_rate),
            color = colors[[group]],
            shape = 1,
            size = 3,
            stroke = 1.2 # 空心圆
          )
      } else {
        p <- p +
          geom_smooth(
            data = d_s,
            aes(prob, outcome),
            method = "loess",
            span = loess_span,
            se = TRUE,
            alpha = 0.15,
            color = colors[[group]],
            fill = colors[[group]],
            linewidth = 1.2
          )
        if (show_bins) {
          if (error_bars) {
            p <- p +
              geom_errorbar(
                data = d_b,
                aes(pred_mean, ymin = ci_low, ymax = ci_high),
                color = colors[[group]],
                width = 0.02,
                linewidth = 0.6,
                alpha = 0.8
              )
          }
          p <- p +
            geom_point(
              data = d_b,
              aes(pred_mean, obs_rate),
              color = colors[[group]],
              shape = 1,
              size = 3,
              stroke = 1.2
            )
        }
        p
      }

      p +
        annotate(
          "text",
          x = 0.02,
          y = 0.90,
          hjust = 0,
          size = 3.6,
          label = sprintf(
            "Brier %.3f\nIntercept %.2f\nSlope %.2f",
            met$brier,
            met$intercept,
            met$slope
          )
        ) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(
          x = "Predicted probability",
          y = "Actual probability",
          title = group
        ) +
        theme_set
    }) |>
      set_names(group_names)

    if (assemble) {
      patchwork::wrap_plots(plots, ncol = panel_ncol) # 不加 tag，留给调用侧
    } else {
      plots
    }
  }

  out <- switch(
    type,
    overlay = make_overlay(),
    panel = make_panel()
  )
  attr(out, "calibration_metrics") <- calibration_metrics
  out
}

##----- example ------

# # pred list
# preds <- list(
#   lncrna = train_lncrna_pred,
#   ca153 = train_ca153_pred,
#   cea = train_cea_pred
# )

# # overlay

# # bins + overlay模式（默认）
# p_cal <- plot_cal_curves(preds, type = "overlay")   # 默认bins 分箱折线

# # LOESS 平滑曲线 + 分箱点 + Wilson 误差棒（默认）
# p_cal_sm <- plot_cal_curves(preds, type = "overlay", method = "loess", loess_span = 0.9)

# # 只要平滑曲线，无点
# plot_cal_curves(preds, method = "loess", n_bins = 5, show_bins = FALSE)
#
# # 要点但不要误差棒
# plot_cal_curves(preds, method = "loess", n_bins = 5,error_bars = FALSE)

# # panel

# # bins
# p_cal_panel <- plot_cal_curves(preds, type = "panel", method = "bins", n_bins = 5)
# # loess平滑
# p_cal_panel <- plot_cal_curves(preds, type = "panel", method = "loess", loess_span = 0.9)

# #single plot list
# cal_list <- plot_cal_curves(preds, type = "panel", assemble = FALSE)

## attr(p_cal, "calibration_metrics")   # Brier / 截距 / 斜率表

## # 九面板：ROC + DCA + Calibration 三行
## nine <- patchwork::wrap_plots(
##   c(
##     plot_roc_curves(preds, type = "panel", assemble = FALSE, theme_set = base_theme),
##     plot_dca_curves(preds, type = "panel", assemble = FALSE, theme_set = base_theme),
##     plot_cal_curves(preds, type = "panel", assemble = FALSE, theme_set = base_theme)
##   ),
##   ncol = 3
## ) + patchwork::plot_annotation(tag_levels = "A")
##
