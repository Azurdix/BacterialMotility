# =============================================================================================
# Author: Mateusz Glenszczyk
# Email: mateusz.glenszczyk@gmail.com
# Date: 2026-03-19
# Description: Bacterial-Motility
# =============================================================================================

#==============================================================================================
# My PhD script to analyse bacterial-motility on low-agar content plates. You may first need to read
# about these tests. For this script to work, you need to have another R file with data matrix. 
# My data matrices had a clear structure - diameter (cm) in 12 repetitions x 4 treatments x 2 species. 
# Therefore you need to have according matrices for this to work... accordingly :D
# Warning - I upload this script faster than my publication - So you may want to look for
# Glenszczyk et al., 2026 or 2027 :) Code may have some parts of polish language. Sorr for that.
# Overall it's english version. This code was created to be highly adaptable in case future experiments.
#==============================================================================================

# =========================================
# LIBRARIES
# =========================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(grid)
  library(tibble)  # rownames_to_column()
})

source("/Users/Mateusz/OneDrive/Desktop/motility_data_matrices.R")

if (!exists("M")) stop("Obiekt 'M' nie istnieje po source(). Sprawdź motility_data_matrices.R")

# =========================================
# THINGIES TO TWEAK
# =========================================
TEST <- "wilcox"
ADJ  <- "holm"
ERR_TYPE <- "SE"  # albo "SD"

SPIDER_COLS   <- c(PRD = "#274DF5", PTSD = "#0E7005")
SPIDER_SHAPES <- c(PRD = 21, PTSD = 21)

SPIDER_LABELS <- c(
  PRD  = expression(italic("Pardosa lugubris")),
  PTSD = expression(italic("Parasteatoda tepidariorum"))
)

SAMPLE_LEVELS <- c("CONTROL","ANTIBIO","EGGS","SILK")
SAMPLE_LABELS <- c(
  CONTROL = "CONTROL",
  ANTIBIO = "ANTIBIOTIC",
  EGGS    = "EGGS",
  SILK    = "SILK"
)

# oś Y: dokładne podziały co 0.25 cm do 2.5 cm
Y_BREAKS <- seq(0, 2.75, by = 0.25)

p_to_stars <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

# =========================================
# LONG FORMAT
# =========================================
df <- as.data.frame(M) |>
  rownames_to_column("rep") |>
  mutate(
    strain = case_when(
      str_detect(rep, "^Ecoli_") ~ "Escherichia coli (NCTC 9001)",
      str_detect(rep, "^Paeruginosa_") ~ "Pseudomonas aeruginosa (NCTC 10662)",
      TRUE ~ "Unknown"
    )
  ) |>
  pivot_longer(
    cols = -c(rep, strain),
    names_to = "var",
    values_to = "value"
  ) |>
  mutate(
    spider = ifelse(str_detect(var, "^PRD_"), "PRD", "PTSD"),
    sample = str_remove(var, "^(PRD|PTSD)_"),
    sample = factor(sample, levels = SAMPLE_LEVELS),
    spider = factor(spider, levels = c("PRD","PTSD"))
  )

# =========================================
# PLOT: BOXPLOT + JITTER
# =========================================
dodge_w  <- 0.78
jitter_w <- 0.20

p <- ggplot(df, aes(x = sample, y = value, fill = spider)) +
  geom_boxplot(
    aes(group = interaction(sample, spider)),
    position = position_dodge(width = dodge_w),
    width = 0.62,
    notch = FALSE,
    outlier.shape = NA,
    alpha = 0.90,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_point(
    aes(shape = spider),
    position = position_jitterdodge(
      jitter.width = jitter_w,
      dodge.width  = dodge_w
    ),
    size = 2.6,
    alpha = 0.45,
    colour = "black",
    stroke = 0.35
  ) +
  facet_wrap(~ strain, nrow = 1) +
  scale_x_discrete(labels = SAMPLE_LABELS) +
  scale_y_continuous(
    limits = c(0, 2.75),
    breaks = Y_BREAKS,
    labels = function(x) sprintf("%.2f", x)  # równe formatowanie (np. 0.25, 1.00)
  ) +
  scale_fill_manual(values = SPIDER_COLS, labels = SPIDER_LABELS, drop = FALSE) +
  scale_shape_manual(values = SPIDER_SHAPES, labels = SPIDER_LABELS, drop = FALSE) +
  labs(
    x = NULL,
    y = "Motility Diameter (cm)",
    fill = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.ticks.length = grid::unit(3, "pt"),
    legend.position = "bottom",
    plot.margin = margin(8, 10, 8, 10),
    panel.spacing = grid::unit(8, "pt")
  ) +
  guides(
    fill = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  )

print(p)

# ==================================================================
# COMPARISONS // (ANTIBIO/EGGS/SILK) VS CONTROL in (STRAIN X SPIDER)
# ==================================================================
compare_to_control <- function(d){
  d <- d |> filter(!is.na(value))
  
  ctrl <- d |> filter(sample == "CONTROL") |> pull(value)
  if (length(ctrl) < 1) return(tibble())
  
  others <- d |> filter(sample != "CONTROL") |> distinct(sample) |> pull(sample)
  
  res <- lapply(others, function(smp){
    x <- d |> filter(sample == smp) |> pull(value)
    if (length(x) < 1) return(NULL)
    
    if (TEST == "t") {
      tt <- t.test(x, ctrl)
      tibble(
        comparison = paste0(as.character(smp), " vs CONTROL"),
        test = "t-test",
        statistic = unname(tt$statistic),
        p_value = tt$p.value,
        n_treat = length(x),
        n_ctrl = length(ctrl),
        median_treat = median(x),
        median_ctrl = median(ctrl)
      )
    } else {
      wt <- wilcox.test(x, ctrl, exact = FALSE)
      tibble(
        comparison = paste0(as.character(smp), " vs CONTROL"),
        test = "Wilcoxon rank-sum",
        statistic = unname(wt$statistic),
        p_value = wt$p.value,
        n_treat = length(x),
        n_ctrl = length(ctrl),
        median_treat = median(x),
        median_ctrl = median(ctrl)
      )
    }
  })
  
  bind_rows(res)
}

tests_vs_control <- df |>
  group_by(strain, spider) |>
  group_modify(~ compare_to_control(.x)) |>
  ungroup() |>
  group_by(strain, spider) |>
  mutate(
    p_adj = p.adjust(p_value, method = ADJ),
    sig   = p_to_stars(p_adj)
  ) |>
  ungroup() |>
  arrange(strain, spider, comparison)

tests_vs_control

# =========================================
# TABLE (THE PRETTIER ONE)
# =========================================
tests_table <- tests_vs_control |>
  mutate(
    sample = gsub(" vs CONTROL", "", comparison),
    sample = factor(sample, levels = c("ANTIBIO","EGGS","SILK")),
    sample_label = recode(as.character(sample), !!!SAMPLE_LABELS),
    spider_label = recode(as.character(spider),
                          PRD  = "Pardosa lugubris",
                          PTSD = "Parasteatoda tepidariorum")
  ) |>
  arrange(strain, spider, sample) |>
  transmute(
    strain,
    spider = spider_label,
    comparison = paste0(sample_label, " vs CONTROL"),
    test,
    n_ctrl, n_treat,
    median_ctrl,
    median_treat,
    statistic = round(statistic, 3),
    p_value = signif(p_value, 3),
    p_adj   = signif(p_adj, 3),
    sig
  )

tests_table

# =================================================================
# ADDING SIGNIFICANCE STARS TO THE PLOT, SO ITS EVEN MORE BEAUTIFUL
# =================================================================
stars_df <- df |>
  group_by(strain, spider, sample) |>
  summarise(y = max(value, na.rm = TRUE), .groups = "drop") |>
  left_join(
    tests_vs_control |>
      mutate(sample = gsub(" vs CONTROL", "", comparison)) |>
      select(strain, spider, sample, sig),
    by = c("strain","spider","sample")
  ) |>
  filter(sample != "CONTROL") |>
  mutate(y = pmin(y + 0.15, 2.65))  # żeby nie wychodziło poza limit 2.5

final_plot <- p + geom_text(
  data = stars_df,
  aes(x = sample, y = y, label = sig, group = spider),
  position = position_dodge(width = dodge_w),
  vjust = 0,
  size = 4
)

final_plot

# =========================================
# SAVING IT
# =========================================

ggsave(
  filename = "motility_boxplot_vs_control.tiff",
  plot     = final_plot,
  width    = 12,
  height   = 7,
  dpi      = 600,
  units    = "in",
  compression = "lzw"
)
