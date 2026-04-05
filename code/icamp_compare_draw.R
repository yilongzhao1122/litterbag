library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)

## =========================
## 全局设置（颜色 + 次序）
## =========================
proc_levels <- c("Heterogeneous.Selection", "Homogeneous.Selection",
                 "Dispersal.Limitation", "Homogenizing.Dispersal",
                 "Drift.and.Others")
proc_labels <- c("Heterogeneous\nSelection", "Homogeneous\nSelection",
                 "Dispersal\nLimitation", "Homogenizing\nDispersal",
                 "Drift")

groups_3 <- c("half","normal","double")
group_colors <- c("half"="#f39b7f","normal"="#3c5488","double"="#4dbbd5")

pair_levels <- c("half-normal","normal-double","half-double")
offsets <- c(half=-0.27, normal=0.00, double=+0.27)  # 三根柱的相对位置
pair_height <- c(`half-normal`=0.90, `normal-double`=0.95, `half-double`=1.00)
label_lift  <- 0.02

## ==========================================
## 数据准备函数（5 个过程；三组；两两比较）
##  - 删除 BootSummary 的第13列以避免重复列名问题
##  - Compare 保留“等级标签 + 星号”，不做数值 d 的计算
## ==========================================
prepare_icamp_data <- function(boot_summary_file, compare_file, community){
  # --- BootSummary：读入后先删第13列，再筛选 ---
  boot_raw <- read.csv(boot_summary_file, check.names = FALSE)
  if(ncol(boot_raw) >= 13) boot_raw <- boot_raw[ , -13]  # 删除第13列
  
  boot_summary <- boot_raw %>%
    dplyr::filter(Group %in% groups_3,
                  Process %in% proc_levels) %>%
    dplyr::transmute(
      Group   = factor(Group, levels = groups_3),
      Process = factor(Process, levels = proc_levels, labels = proc_labels),
      Observed = as.numeric(Observed),
      Stdev    = as.numeric(Stdev),
      Community = community
    )
  
  # --- Compare：三种配对；使用 等级标签(Effect.Size) + 星号(P.value) ---
  cmp_raw <- read.csv(compare_file, check.names = FALSE)
  
  cmp_long <- cmp_raw %>%
    dplyr::mutate(
      Pair = dplyr::case_when(
        (Group1=="half"   & Group2=="normal") | (Group1=="normal" & Group2=="half")   ~ "half-normal",
        (Group1=="normal" & Group2=="double") | (Group1=="double" & Group2=="normal") ~ "normal-double",
        (Group1=="half"   & Group2=="double") | (Group1=="double" & Group2=="half")   ~ "half-double",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Pair)) %>%
    # 把 *_Effect.Size / *_P.value 两类列收成长表
    tidyr::pivot_longer(
      cols = dplyr::matches("(_Effect.Size|_P.value)$"),
      names_to = c("ProcessRaw", ".value"),
      names_sep = "_"
    ) %>%
    dplyr::mutate(
      Process = factor(ProcessRaw, levels = proc_levels, labels = proc_labels),
      # 保留 Compare 里写好的“大小标识”（small/medium/large/negligible 等）
      Effect.Size = as.character(Effect.Size),
      # p 值转星号
      P.value = suppressWarnings(as.numeric(P.value)),
      Significance = dplyr::case_when(
        !is.na(P.value) & P.value < 0.001 ~ "***",
        !is.na(P.value) & P.value < 0.01  ~ "**",
        !is.na(P.value) & P.value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Pair = factor(Pair, levels = pair_levels)
    ) %>%
    # 每个 Process × Pair 只留一条（若有重复就取第一条非空标签）
    dplyr::group_by(Pair, Process) %>%
    dplyr::summarise(
      EffectLabel = dplyr::first(Effect.Size[!is.na(Effect.Size)]),
      Significance = dplyr::first(Significance),
      .groups = "drop"
    ) %>%
    # 组装最终用于标注的文本：如 "medium **"
    dplyr::mutate(Label = stringr::str_trim(paste0(EffectLabel, " ", Significance)))
  
  list(boot = boot_summary, cmp = cmp_long)
}

## =======================================================
## 绘图函数：三组柱 + 每个过程三条“括号线”（pairwise）
##  - 括号旁的文字用 “等级标签 + 星号”
##  - 使用你指定的新颜色
## =======================================================
create_process_plot <- function(boot_df, cmp_df, community){
  plot_data <- boot_df %>% dplyr::filter(Community == community)
  
  y_max <- max(plot_data$Observed + plot_data$Stdev, na.rm = TRUE) * 1.25
  
  base_df <- expand.grid(
    Process = levels(plot_data$Process),
    Pair    = pair_levels,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) %>% 
    dplyr::as_tibble() %>%
    dplyr::left_join(cmp_df, by = c("Process","Pair")) %>%
    dplyr::mutate(
      x_center = as.numeric(factor(Process, levels = levels(plot_data$Process))),
      x1 = dplyr::case_when(
        Pair == "half-normal"   ~ x_center + offsets["half"],
        Pair == "normal-double" ~ x_center + offsets["normal"],
        Pair == "half-double"   ~ x_center + offsets["half"]
      ),
      x2 = dplyr::case_when(
        Pair == "half-normal"   ~ x_center + offsets["normal"],
        Pair == "normal-double" ~ x_center + offsets["double"],
        Pair == "half-double"   ~ x_center + offsets["double"]
      ),
      y      = y_max * unname(pair_height[Pair]),
      ytick  = y - 0.01 * y_max,
      label_x = (x1 + x2) / 2,
      label_y = y + label_lift * y_max
    )
  
  ggplot(plot_data, aes(x = Process, y = Observed, fill = Group)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8),
             width = 0.7, color = "black", size = 0.5) +
    geom_errorbar(aes(ymin = pmax(0, Observed - Stdev),
                      ymax = Observed + Stdev),
                  position = position_dodge(width = 0.8),
                  width = 0.25, size = 0.8, color = "black") +
    # 括号线
    geom_segment(data = base_df,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, size = 0.8, color = "black") +
    geom_segment(data = base_df,
                 aes(x = x1, xend = x1, y = ytick, yend = y),
                 inherit.aes = FALSE, size = 0.8, color = "black") +
    geom_segment(data = base_df,
                 aes(x = x2, xend = x2, y = ytick, yend = y),
                 inherit.aes = FALSE, size = 0.8, color = "black") +
    # 等级标签 + 星号（例如 "medium **"）
    geom_text(data = base_df,
              aes(x = label_x, y = label_y, label = Label),
              inherit.aes = FALSE, size = 4.8, fontface = "bold") +
    scale_fill_manual(values = group_colors, breaks = groups_3) +
    labs(
      title = paste("Process Importance in", community),
      subtitle = "Pairwise comparisons: half–normal, normal–double, half–double",
      x = "Ecological Process",
      y = "Relative Importance",
      fill = "Precipitation"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    ) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.1)))
}


# =========================================================
# 3. 加载数据（Bacteria / Fungi）
#    使用你上传的 /mnt/data 路径
# =========================================================
# 细菌
bac <- prepare_icamp_data(
  boot_summary_file = "/home/ZhaoYilong/litterbag/iCAMP2/iCAMP_16S_bin=48/litterbag.16S.iCAMP.BootSummary.csv",
  compare_file      = "/home/ZhaoYilong/litterbag/iCAMP2/iCAMP_16S_bin=48/litterbag.16S.iCAMP.Compare.csv",
  community         = "Bacteria"
)

# 真菌
fun <- prepare_icamp_data(
  boot_summary_file = "/home/ZhaoYilong/litterbag/iCAMP2/iCAMP_ITS_bin=48/litterbag.ITS..iCAMP.BootSummary.csv",
  compare_file      = "/home/ZhaoYilong/litterbag/iCAMP2/iCAMP_ITS_bin=48/litterbag.ITS..iCAMP.Compare.csv",
  community         = "Fungi"
)

# =========================================================
# 4. 生成图
# =========================================================
plot_bac <- create_process_plot(bac$boot, bac$cmp, "Bacteria")
plot_fun <- create_process_plot(fun$boot, fun$cmp, "Fungi")

combined_plot <- ggarrange(
  plot_bac, plot_fun,
  ncol = 2,
  labels = c("A", "B"),
  font.label = list(size = 18, face = "bold"),
  common.legend = TRUE,
  legend = "top"
)

final_plot <- annotate_figure(
  combined_plot,
  top = text_grob(
    "Comparative Analysis of Community Assembly Processes (Pairwise only)",
    face = "bold", size = 18, color = "black"
  ),
  bottom = text_grob(
    paste(
      "Error bars: Standard deviation (Stdev) | ",
      "Effect size labels from iCAMP Compare | ",
      "Significance: *p<0.05, **p<0.01, ***p<0.001"
    ),
    color = "grey40", size = 11, face = "italic"
  )
)

ggsave("Process_Importance_Comparison_3groups_pairwise.png",
       final_plot, width = 16, height = 9, dpi = 300)

# =========================================================
# 5. 可选：Cohen's d 解释图例（保持你的原版）
# =========================================================
effect_legend <- ggplot() +
  annotate("text", x = 0, y = 0.8, label = "Cohen's d Effect Size Interpretation",
           size = 5, fontface = "bold", hjust = 0) +
  annotate("text", x = 0, y = 0.6, label = "Large effect: |d| ≥ 0.8",
           size = 4, hjust = 0, color = "#D55E00") +
  annotate("text", x = 0, y = 0.4, label = "Medium effect: 0.5 ≤ |d| < 0.8",
           size = 4, hjust = 0, color = "#E69F00") +
  annotate("text", x = 0, y = 0.2, label = "Small effect: 0.2 ≤ |d| < 0.5",
           size = 4, hjust = 0, color = "#56B4E9") +
  annotate("text", x = 0, y = 0.0, label = "Negligible effect: |d| < 0.2",
           size = 4, hjust = 0, color = "#999999") +
  theme_void() + xlim(0,1) + ylim(0,1)

final_with_legend <- ggarrange(
  final_plot, effect_legend,
  ncol = 1,
  heights = c(5, 0.5)
)

# ggsave("Process_Importance_Comparison_3groups_pairwise_withLegend.png",
#        final_with_legend, width = 16, height = 10, dpi = 300)

# =========================================================
# Stochasticity：数据准备（仅 Process == "Stochasticity"）
# =========================================================
prepare_stochasticity_data <- function(boot_summary_file, compare_file, community){
  boot <- read.csv(boot_summary_file, check.names = FALSE) %>%
    { .[, -13] } %>%
    filter(Group %in% groups_3, Process == "Stochasticity") %>%
    transmute(
      Group = factor(Group, levels = groups_3),
      Observed = as.numeric(Observed),
      Stdev    = as.numeric(Stdev),
      Community = community
    )
  
  cmp_raw <- read.csv(compare_file, check.names = FALSE)
  
  cmp <- cmp_raw %>%
    mutate(
      Pair = case_when(
        (Group1=="half"   & Group2=="normal") | (Group1=="normal" & Group2=="half")   ~ "half-normal",
        (Group1=="normal" & Group2=="double") | (Group1=="double" & Group2=="normal") ~ "normal-double",
        (Group1=="half"   & Group2=="double") | (Group1=="double" & Group2=="half")   ~ "half-double",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Pair)) %>%
    select(Pair, Stochasticity_Effect.Size, Stochasticity_P.value) %>%
    mutate(
      Effect.Size = suppressWarnings(as.numeric(Stochasticity_Effect.Size)),
      P.value     = suppressWarnings(as.numeric(Stochasticity_P.value)),
      Significance = case_when(
        !is.na(P.value) & P.value < 0.001 ~ "***",
        !is.na(P.value) & P.value < 0.01  ~ "**",
        !is.na(P.value) & P.value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      EffectLabel = ifelse(is.na(Effect.Size), "", sprintf("%.2f", Effect.Size)),
      Pair = factor(Pair, levels = pair_levels)
    ) %>%
    group_by(Pair) %>%
    summarise(EffectLabel = paste(unique(EffectLabel), collapse = ","),
              Significance = paste(unique(Significance), collapse = ""),
              .groups="drop")
  
  list(boot=boot, cmp=cmp)
}

# 载入
bac_stoch <- prepare_stochasticity_data(
  "/home/ZhaoYilong/litterbag/iCAMP2/iCAMP_16S_bin=48/litterbag.16S.iCAMP.BootSummary.csv",
  "/home/ZhaoYilong/litterbag/iCAMP2/iCAMP_16S_bin=48/litterbag.16S.iCAMP.Compare.csv",
  "Bacteria"
)
fun_stoch <- prepare_stochasticity_data(
  "/home/ZhaoYilong/litterbag/iCAMP2/iCAMP_ITS_bin=48/litterbag.ITS..iCAMP.BootSummary.csv",
  "/home/ZhaoYilong/litterbag/iCAMP2/iCAMP_ITS_bin=48/litterbag.ITS..iCAMP.Compare.csv",
  "Fungi"
)

# 绘图函数（单群落）
create_single_stochasticity_plot <- function(boot_cmp_list, community_title=""){
  data <- boot_cmp_list$boot
  cmp  <- boot_cmp_list$cmp
  
  y_max <- max(data$Observed + data$Stdev, na.rm = TRUE) * 1.25
  
  # 三条配对括号（x 在 1 个类目内的左中右）
  pair_df <- tibble(
    Pair = factor(pair_levels, levels = pair_levels),
    x1 = c(1 + offsets["half"], 1 + offsets["normal"], 1 + offsets["half"]),
    x2 = c(1 + offsets["normal"], 1 + offsets["double"], 1 + offsets["double"]),
    y  = y_max * unname(pair_height[pair_levels]),
    ytick = y - 0.01 * y_max,
    label_x = (x1 + x2) / 2,
    label_y = y + label_lift * y_max
  ) %>% left_join(cmp %>% mutate(Pair=factor(Pair, levels=pair_levels)), by="Pair") %>%
    mutate(Label = str_trim(paste0(EffectLabel, " ", Significance)))
  
  ggplot(data, aes(x = factor(1), y = Observed, fill = Group)) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.7, color = "black", size = 0.5) +
    geom_errorbar(aes(ymin = pmax(0, Observed - Stdev), ymax = Observed + Stdev),
                  position = position_dodge(width = 0.8),
                  width = 0.18, size = 0.8, color = "black") +
    # 阈值线（50%）
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = 0.5, y = 0.52, label = "Random/Deterministic Threshold (50%)",
             color = "red", size = 4, hjust = 0) +
    # 括号
    geom_segment(data = pair_df, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, color = "black", size = 0.8) +
    geom_segment(data = pair_df, aes(x = x1, xend = x1, y = ytick, yend = y),
                 inherit.aes = FALSE, color = "black", size = 0.8) +
    geom_segment(data = pair_df, aes(x = x2, xend = x2, y = ytick, yend = y),
                 inherit.aes = FALSE, color = "black", size = 0.8) +
    geom_text(data = pair_df, aes(x = label_x, y = label_y, label = Label),
              inherit.aes = FALSE, size = 5, fontface = "bold") +
    scale_x_discrete(labels = "Stochasticity") +
    scale_fill_manual(values = group_colors, breaks = groups_3) +
    labs(title = paste(community_title, "Stochasticity"),
         x = NULL, y = "Relative Importance", fill = "Treatment") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      legend.position = "top"
    ) +
    scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.1)))
}

bac_plot <- create_single_stochasticity_plot(bac_stoch, "Bacterial")
fun_plot <- create_single_stochasticity_plot(fun_stoch, "Fungal")
# ggsave("Stochasticity_3groups_pairwise_Bacteria.png", bac_plot, width=6, height=5, dpi=300)
# ggsave("Stochasticity_3groups_pairwise_Fungi.png",    fun_plot, width=6, height=5, dpi=300)

