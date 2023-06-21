# Attach requirement packages and setting WD -----------------------------------
packages_names <-
  c(
    "tidyverse",
    "data.table",
    "readxl",
    "reshape2",
    "rstudioapi",
    "caTools",
    "car",
    "quantmod",
    "MASS",
    "corrplot",
    "janitor",
    "nortest",
    "survminer",
    "survival",
    "broom",
    "forestplot"
  )

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate


# Reading data
merged_data <-
  fread("Merged annotated data AURKA, KRAS, TP53, EGFR.csv") %>%
  as_tibble()

# 1 Question -------------------------------------------------------------------
merged_data %>%
  filter(!is.na(AURKA_cna)) %>%
  mutate(AURKA_cna = as.factor(AURKA_cna)) %>%
  ggplot(aes(x = AURKA_CNA_log2, fill = AURKA_cna)) +
  geom_density(alpha = .6) +
  theme_minimal()

merged_data %>%
  filter(!is.na(AURKA_cna)) %>%
  mutate(AURKA_cna = as.factor(AURKA_cna)) %>%
  ggplot(aes(y = AURKA_CNA_log2, x = AURKA_cna, fill = AURKA_cna)) +
  geom_violin(trim = F) +
  theme_minimal()

merged_data %>%
  filter(!is.na(AURKA_cna), !is.na(AURKA_CNA_log2)) %>%
  mutate(AURKA_cna = as.factor(AURKA_cna)) %>%
  group_by(AURKA_cna) %>%
  reframe(
    min = min(AURKA_CNA_log2),
    max = max(AURKA_CNA_log2)
  )

# 2 Question ------------
plot_violin_chart <- function(df, 
                              ylab = "Density",
                              stat_test = "wilcox.test") {
  p <- df %>% 
    ggviolin(
      y = "value",
      x = "variable",
      fill = "variable",
      add = "boxplot",
      add.params = list(fill = "white"),
      palette = c("#00AFBB", "#E7B800")
    ) +
    stat_compare_means(
      aes(group = variable),
      comparisons = list(c("a", "b")),
      method = {stat_test}
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      x = "Value",
      y = {ylab},
      title = {stat_test}
    )
  return(p)
}

generate_random_data <- function(n = 100,
                                 mean_1 = 10,
                                 mean_2 = 11,
                                 sd_1 = 2,
                                 sd_2 = 2,
                                 z_trans = F) {
  df <- tibble(
    a = rnorm(n = {n}, mean = {mean_1}, sd = {sd_1}),
    b = rnorm(n = {n}, mean = {mean_2}, sd = {sd_2})
  ) %>% 
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>% 
    mutate(variable = as.factor(variable))

  return(df)
}

compare_chart <- function(n = 100,
                          mean_1 = 10,
                          mean_2 = 11,
                          sd_1 = 2,
                          sd_2 = 2,
                          stat_test = "wilcox.test") {
  
  df_plain = generate_random_data(n = n,
                                  mean_1 = mean_1,
                                  mean_2 = mean_2,
                                  sd_1 = sd_1,
                                  sd_2 = sd_2)
  
  df_z_trans = df_plain %>% 
    mutate(value = log(value, base = 2),
           value = scale(value))
  
  chart_plain = plot_violin_chart(df_plain,
                                  ylab = "Density",
                                  stat_test = {stat_test})
  
  chart_z_trans = plot_violin_chart(df_z_trans,
                                  ylab = "Density, log2, Z-scored",
                                  stat_test = {stat_test})
  
  combined_chart = ggarrange(chart_plain, chart_z_trans)
  
  return(combined_chart)
}

compare_chart(stat_test = "t.test")
compare_chart(stat_test = "wilcox.test")

# 3 Question ------------
# 4 Question ------------
processed_rna_data_tp53 <- merged_data %>%
  filter(EGFR != KRAS) %>%
  select(EGFR, TP53, KRAS, AURKA_CNA_log2) %>%
  mutate(
    EGFR_only = ifelse(EGFR == "ALT" & KRAS == "WT", "ALT", "WT"),
    KRAS_only = ifelse(KRAS == "ALT" & EGFR == "WT", "ALT", "WT"),
    across(
      !AURKA_CNA_log2,
      .fns = function(x) {
        as.factor(x)
      }
    )
  ) %>%
  pivot_longer(
    cols = contains("only"),
    names_to = "subgroup",
    values_to = "alteration"
  ) %>%
  filter(!alteration == "WT",
         !is.na(AURKA_CNA_log2)) %>% 
  mutate(TP53 = paste("TP53", TP53, sep = "_"),
         subgroup = paste(subgroup, TP53, sep = "_")) %>% 
  select(-EGFR, -KRAS, -alteration, -TP53) %>% 
  mutate(subgroup = case_when(subgroup == "KRAS_only_TP53_WT" ~ "K",
                              subgroup == "KRAS_only_TP53_ALT" ~ "KT",
                              subgroup == "EGFR_only_TP53_WT" ~ "E",
                              subgroup == "EGFR_only_TP53_ALT" ~ "ET")) %>% 
  mutate(subgroup = as.factor(subgroup))

processed_rna_data_tp53 %>%
  ggviolin(
    y = "AURKA_CNA_log2",
    x = "subgroup",
    fill = "subgroup",
    add = "boxplot",
    add.params = list(fill = "white"),
    palette = c("#EEC900", "#B0E2FF",
                "#EEC900", "#B0E2FF")
  ) +
  theme(legend.position = "none") +
  stat_compare_means(
    aes(group = variable),
    comparisons = list(c("E", "ET"), c("K", "KT")),
    method = "wilcox.test",
    label = "p.signif"
  )
