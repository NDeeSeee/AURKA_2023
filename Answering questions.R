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
df_test <- tibble(
  a = rnorm(100, mean = 10, sd = 2),
  b = rnorm(100, mean = 11, sd = 1.5)
) %>% 
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>% 
  mutate(variable = as.factor(variable))

plot_violin_chart <- function(df, ylab = "Density") {
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
      method = "t.test"
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      x = "Value",
      y = {ylab},
      title = "Violin plots with mean comparison"
    )
  return(p)
  }

plot_violin_chart(df_test)

df_test %>% 
  mutate(value = log(value, base = 2),
         value = scale(value)) %>% 
  plot_violin_chart(ylab = "Density, log2, Z-scored")

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
