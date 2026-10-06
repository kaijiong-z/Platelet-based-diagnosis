# =====================================================================
# Fig. 2. the association between platelet parameters and breast cancer
# =====================================================================

rm(list = ls())

library(tidyverse)
library(broom)
library(ggpubr)
library(skimr)
library(pROC)
library(ggsignif)
library(ggsci)
library(scales)
library(RColorBrewer)
library(cowplot)

## ======= fig2a.data distribution of participants =======

## load data
load("rawdata.rdata")

# Sample statistics
rawdata_groups_p <- rawdata %>%
  dplyr::group_by(groups) %>%
  dplyr::summarise(count = n()) %>%
  ggplot(aes(x = groups, y = count, fill = groups)) +
  geom_bar(stat = "identity") +
  labs(title = "Numbers of samples", x = "", y = "") +
  scale_fill_manual(values = c("#075f75", "#ef9c02", "#b12115")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.position = "none"
  )

rawdata_groups_p

## age distribution
age_groups_p <- rawdata %>%
  ggplot(aes(x = groups, y = age, fill = groups)) +
  geom_boxplot(width = 0.7) +
  scale_color_manual(values = c("#b12115", "#ef9c02", "#075f75")) +
  scale_fill_manual(values = c("#b12115", "#ef9c02", "#075f75")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  geom_signif(
    comparisons = list(
      c("tumor", "benign"),
      c("benign", "normal"),
      c("tumor", "normal")
    ),
    #map_signif_level =  function(p) sprintf("p = %.2g", p),
    step_increase = 0.2,
    map_signif_level = TRUE,
    position = "identity",
    test = "wilcox.test",
    textsize = 6,
    tip_length = 0.05,
    margin_top = 0.05,
    vjust = 0.2
  ) +
  labs(title = "Distribution of age", x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "none"
    #plot.margin=unit(rep(1,4),'lines')
  )

age_groups_p


# Tumor stage statistics
tumor_stage_p <- tumor_data %>%
  mutate(stage = ifelse(stage == "0", "I", stage)) %>%
  dplyr::group_by(stage) %>%
  dplyr::summarise(count = n()) %>%
  ggplot(aes(x = stage, y = count, fill = stage)) +
  geom_bar(stat = "identity") +
  labs(title = "AJCC8thTNM Stage", x = "", y = "") +
  scale_fill_manual(
    values = c("#075f75", "#b12115", "#7363ac", "#ef9c02", "#5f5f5f")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.position = "none"
  )

tumor_stage_p

# graphic combination
rawdata_plot <- plot_grid(
  rawdata_groups_p,
  age_groups_p,
  tumor_stage_p,
  #labels=LETTERS[1:3],
  ncol = 3,
  label_fontface = "bold",
  label_size = 20,
  scale = 1
)


rawdata_plot


## ======= fig2b.platelet parameters expression  =======

plot_column_boxplot <- function(data, column) {
  ggplot(data = data, aes(x = groups, y = get(column), fill = groups)) +
    geom_boxplot(width = 0.7) +
    #geom_point(aes(colour = groups), size = 4, position = position_jitterdodge(jitter.width = 1))  +
    scale_color_manual(values = c("#b12115", "#075f75")) +
    scale_fill_manual(values = c("#b12115", "#075f75")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    geom_signif(
      comparisons = list(c("tumor", "nontumor")),
      #map_signif_level =  function(p) sprintf("p = %.2g", p),
      map_signif_level = TRUE,
      position = "identity",
      test = "wilcox.test",
      textsize = 6,
      tip_length = 0.02,
      margin_top = 0.05,
      vjust = 0.2
    ) +
    labs(title = column, x = "", y = "") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "none"
      #plot.margin=unit(rep(1,4),'lines')
    )
}


# Define a function to remove outliers, processed by group
remove_outliers_by_group <- function(
  data_frame,
  group_column,
  threshold = 1.5
) {
  unique_groups <- unique(data_frame[[group_column]])

  # Loop through each group
  for (group_value in unique_groups) {
    # Subset data frame containing only observations for the specific group
    group_data <- data_frame[data_frame[[group_column]] == group_value, ]

    # Loop through each column in the grouped data frame
    for (col in names(group_data)) {
      # Process only numeric variables
      if (is.numeric(group_data[[col]])) {
        # Calculate Tukey's IQR (Interquartile Range)
        q1 <- quantile(group_data[[col]], 0.25, na.rm = TRUE)
        q3 <- quantile(group_data[[col]], 0.75, na.rm = TRUE)
        iqr <- q3 - q1

        # Determine the lower and upper bounds for outliers
        lower_bound <- q1 - threshold * iqr
        upper_bound <- q3 + threshold * iqr

        # Replace outliers with NA
        group_data[[col]][
          group_data[[col]] < lower_bound | group_data[[col]] > upper_bound
        ] <- NA
      }
    }

    # Put the processed grouped data frame back into the original data frame
    data_frame[data_frame[[group_column]] == group_value, ] <- group_data
  }

  return(data_frame)
}

# Remove outliers by group
data_without_outliers_by_group <- remove_outliers_by_group(
  rawdata,
  group_column = "groups"
)


data <- data_without_outliers_by_group %>%
  dplyr::select(PLT:PDW, groups) %>%
  mutate(groups = ifelse(groups == "tumor", "tumor", "nontumor")) %>%
  mutate(groups = factor(groups, levels = c("tumor", "nontumor")))


myplots_boxplot <- lapply(colnames(data)[1:3], plot_column_boxplot, data = data)


# graphic combination
feature_expr <- plot_grid(
  myplots_boxplot[[1]],
  myplots_boxplot[[2]],
  myplots_boxplot[[3]],
  #labels=LETTERS[1:8],
  ncol = 3,
  label_fontface = "bold",
  label_size = 20,
  scale = 1
)


feature_expr


##======= fig2c.survival analysis  =======

library(survival)
library(survminer)

## load Survival data
load("BC_surdata.Rdata")
load("rawdata.rdata")

sur_data <- rawdata %>%
  filter(!groups == "benign") %>%
  distinct(id, .keep_all = T) %>%
  inner_join(BC_surdata, by = "id") %>%
  dplyr::select(c("status", "time", "PLT", "MPV", "PDW")) %>%
  dplyr::mutate(status = ifelse(status == "death", 1, 0), time = time / 30)

## Determine the optimal cutoff value----
vars <- sur_data %>% dplyr::select(-c("status", "time")) %>% colnames()
sur.cut <- surv_cutpoint(
  sur_data,
  time = "time",
  event = "status",
  variables = vars
) ## Select cutoff value

# summary(sur.cut)# Output cutoff
# plot(sur.cut, "MPV", palette = "npg")# Plot graph
df <- surv_categorize(sur.cut) ## Convert continuous variables to categorical variables


## Survival curve function
surv_plot_func <- function(df, vars, time, status) {
  results_list <- lapply(vars, \(x) {
    # # Creating a formula as a string
    form <<- paste0("Surv(", time, ", ", status, ") ~ ", x)
    fit <- survfit(as.formula(form), data = df)

    # # Plot the Kaplan-Meier curve using ggsurvplot
    ggsurv <- ggsurvplot(
      fit,
      pval = TRUE,
      conf.int = F,
      pval.coord = c(1, 0.9), # p-value location
      pval.size = 6.5, # adjust p-value size
      risk.table = T, #legend.title = "Groups",
      surv.median.line = "hv",
      xlab = "Time(Month)",
      ylab = "Survival Probability(OS)",
      palette = c("#b12115", "#007172"),
      font.x = 15,
      font.y = 15,
      font.tickslab = 15,
      font.legend = 15,
      xlim = c(0, 80),
      ylim = c(0.85, 1),
      ggtheme = theme_classic()
    )

    # Return the fit and ggsurv as a list
    list(fit = fit, ggsurv = ggsurv)
  })

  # Return the list of results
  return(results_list)
}


res_list <- surv_plot_func(df, vars, "time", "status")

## PLT
fit_PLT <- res_list[[1]]$ggsurv$plot

## MPV
fit_MPV <- res_list[[2]]$ggsurv$plot

## PDW
fit_PDW <- res_list[[3]]$ggsurv$plot


## Combine plots
surv_plot <- plot_grid(
  fit_PLT,
  fit_MPV,
  fit_PDW,
  #labels=LETTERS[1:3],
  ncol = 3,
  label_fontface = "bold",
  label_size = 20,
  scale = 1
)

surv_plot


## Combine plots
expr_surv_plot <- plot_grid(
  rawdata_plot,
  feature_expr,
  surv_plot,
  labels = LETTERS[1:3],
  ncol = 1,
  #align = "v",
  axis = "b",
  label_fontface = "bold",
  label_size = 20,
  rel_widths = c(1, 1.8)
) +
  theme(plot.margin = unit(rep(1, 4), 'lines'))

expr_surv_plot


# fig2 save
ggsave("../results/P1.clin_all_plot.pdf", width = 12, height = 10)
ggsave(
  "../results/P1.clin_all_plot.jpg",
  width = 12,
  height = 10,
  bg = "white",
  dpi = 300
)
