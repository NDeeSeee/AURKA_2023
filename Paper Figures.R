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

# Figure XA --------------------------------------------------------------------
processed_rna_data <- merged_data %>%
  filter(EGFR != KRAS) %>%
  select(AURKA_rna_exp, EGFR, KRAS) %>%
  pivot_longer(
    cols = c("EGFR", "KRAS"),
    names_to = "gene",
    values_to = "impact"
  ) %>%
  filter(impact != "WT") %>%
  select(-impact) %>%
  mutate(gene = as.factor(gene))

processed_rna_data %>%
  ggviolin(
    y = "AURKA_rna_exp",
    x = "gene",
    fill = "gene",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test",
                     label.y = 4.5,
                     label.x = 1.3)

ggsave(
  "Paper Figures/Fig XA.png",
  dpi = 400,
  height = 5,
  width = 5,
  units = "in"
)

# Statistics
EGFR_only_AURKA_rna_exp <-
  filter(processed_rna_data, gene == "EGFR")$AURKA_rna_exp
KRAS_only_AURKA_rna_exp <-
  filter(processed_rna_data, gene == "KRAS")$AURKA_rna_exp

# Check for normality
ad.test(EGFR_only_AURKA_rna_exp)
ad.test(KRAS_only_AURKA_rna_exp)

# Compute variance for t.test
var(EGFR_only_AURKA_rna_exp)
# 0.9
var(KRAS_only_AURKA_rna_exp)
# 1

ks.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)
t.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)


# Figure XB --------------------------------------------------------------------
processed_cna_data <- merged_data %>%
  filter(EGFR != KRAS, !is.na(AURKA_CNA_log2)) %>%
  select(AURKA_CNA_log2, EGFR, KRAS) %>%
  pivot_longer(
    cols = c("EGFR", "KRAS"),
    names_to = "gene",
    values_to = "impact"
  ) %>%
  filter(impact != "WT") %>%
  select(-impact) %>%
  mutate(gene = as.factor(gene))

processed_cna_data %>%
  ggviolin(
    y = "AURKA_CNA_log2",
    x = "gene",
    fill = "gene",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test",
                     label.y = 4.5,
                     label.x = 1.3)

ggsave(
  "Paper Figures/Fig XB.png",
  dpi = 400,
  height = 5,
  width = 5,
  units = "in"
)

# Statistics
EGFR_only_AURKA_CNA_log2 <-
  filter(processed_cna_data, gene == "EGFR")$AURKA_CNA_log2
KRAS_only_AURKA_CNA_log2 <-
  filter(processed_cna_data, gene == "KRAS")$AURKA_CNA_log2

# Check for normality
ad.test(EGFR_only_AURKA_CNA_log2)
ad.test(KRAS_only_AURKA_CNA_log2)

# Compute variance for t.test
var(EGFR_only_AURKA_CNA_log2)
# 0.276
var(KRAS_only_AURKA_CNA_log2)
# 0.164

ks.test(EGFR_only_AURKA_CNA_log2, KRAS_only_AURKA_CNA_log2)
t.test(EGFR_only_AURKA_CNA_log2, KRAS_only_AURKA_CNA_log2)


# Figure XC --------------------------------------------------------------------
processed_rna_data_tp53 <- merged_data %>%
  filter(EGFR != KRAS) %>%
  select(EGFR, TP53, KRAS, AURKA_rna_exp) %>%
  mutate(
    EGFR_only = ifelse(EGFR == "ALT" & KRAS == "WT", "ALT", "WT"),
    KRAS_only = ifelse(KRAS == "ALT" & EGFR == "WT", "ALT", "WT"),
    across(
      !AURKA_rna_exp,
      .fns = function(x) {
        as.factor(x)
      }
    )
  ) %>% 
  pivot_longer(cols = contains("only"),
                     names_to = "subgroup",
                     values_to = "alteration") %>%
  filter(!alteration == "WT")

processed_rna_data_tp53 %>%
  ggplot(aes(
    y = AURKA_rna_exp,
    fill = subgroup,
    x = subgroup,
    col = TP53
  )) +
  geom_violin(show.legend = F) +
  geom_boxplot(width = .2, position = position_dodge(.9)) +
  scale_color_manual(values = c("ALT" = "gray80", "WT" = "black")) +
  theme_classic()

my_comparisons <- list(c("EGFR_only", "KRAS_only"))

processed_rna_data_tp53 %>%
  ggviolin(
    y = "AURKA_rna_exp",
    x = "subgroup",
    col = "TP53",
    fill = "subgroup",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
  theme(legend.position = "none") +
  # stat_compare_means(
  #   method = "wilcox.test",
  #   label.y = 4.5,
  #   label.x = 1.3,
  #   paired = T
  # ) +
  scale_color_manual(values = c("ALT" = "gray80", "WT" = "black"))

ggsave(
  "Paper Figures/Fig XC.png",
  dpi = 400,
  height = 5,
  width = 5,
  units = "in"
)

# Statistics
EGFR_only_TP53_mut_AURKA_rna_exp <-
  filter(processed_rna_data_tp53, EGFR == "ALT", TP53 == "ALT")$AURKA_rna_exp
EGFR_only_TP53_wt_AURKA_rna_exp <-
  filter(processed_rna_data_tp53, EGFR == "ALT", TP53 == "WT")$AURKA_rna_exp

KRAS_only_TP53_mut_AURKA_rna_exp <-
  filter(processed_rna_data_tp53, KRAS == "ALT", TP53 == "ALT")$AURKA_rna_exp
KRAS_only_TP53_wt_AURKA_rna_exp <-
  filter(processed_rna_data_tp53, KRAS == "ALT", TP53 == "WT")$AURKA_rna_exp


# Check for normality
ad.test(EGFR_only_TP53_mut_AURKA_rna_exp)
ad.test(EGFR_only_TP53_wt_AURKA_rna_exp)
ad.test(KRAS_only_TP53_mut_AURKA_rna_exp)
ad.test(KRAS_only_TP53_wt_AURKA_rna_exp)

# Compute variance for t.test
var(EGFR_only_TP53_mut_AURKA_rna_exp)
# 0.88
var(EGFR_only_TP53_wt_AURKA_rna_exp)
# 0.71
var(KRAS_only_TP53_mut_AURKA_rna_exp)
# 0.73
var(KRAS_only_TP53_wt_AURKA_rna_exp)
# 1.03

ks.test(EGFR_only_TP53_mut_AURKA_rna_exp,
        EGFR_only_TP53_wt_AURKA_rna_exp)
ks.test(EGFR_only_TP53_wt_AURKA_rna_exp,
        KRAS_only_TP53_wt_AURKA_rna_exp)
ks.test(KRAS_only_TP53_mut_AURKA_rna_exp,
        KRAS_only_TP53_wt_AURKA_rna_exp)
ks.test(EGFR_only_TP53_mut_AURKA_rna_exp,
        KRAS_only_TP53_mut_AURKA_rna_exp)

t.test(EGFR_only_TP53_mut_AURKA_rna_exp,
       EGFR_only_TP53_wt_AURKA_rna_exp)
t.test(EGFR_only_TP53_wt_AURKA_rna_exp,
       KRAS_only_TP53_wt_AURKA_rna_exp)
t.test(KRAS_only_TP53_mut_AURKA_rna_exp,
       KRAS_only_TP53_wt_AURKA_rna_exp)
t.test(EGFR_only_TP53_mut_AURKA_rna_exp,
       KRAS_only_TP53_mut_AURKA_rna_exp)


# Figure XD -----------------------------------------------------------------------
processed_rna_os_data <- merged_data %>%
  filter(EGFR != KRAS) %>%
  select(AURKA_rna_exp, EGFR, KRAS, TP53, os_months, os_status) %>%
  pivot_longer(
    cols = c("EGFR", "KRAS"),
    names_to = "gene",
    values_to = "impact"
  ) %>%
  filter(impact != "WT", !is.na(os_months)) %>%
  select(-impact) %>%
  mutate(os_status = as.numeric(str_remove(os_status, ":.+")))

survfit(
  Surv(time = os_months, os_status) ~ AURKA_rna_exp,
  data = filter(processed_rna_os_data, gene == "EGFR")
)


get_cox_statistics <- function(input_data) {
  coxph(Surv(time = os_months, os_status) ~ AURKA_rna_exp, data = input_data) %>%
    tidy(conf.int = T) %>%
    return()
}

filter(processed_rna_os_data, gene == "EGFR", TP53 == "ALT") %>%
  get_cox_statistics()
filter(processed_rna_os_data, gene == "EGFR", TP53 == "WT") %>%
  get_cox_statistics()
filter(processed_rna_os_data, gene == "KRAS", TP53 == "ALT") %>%
  get_cox_statistics()
filter(processed_rna_os_data, gene == "KRAS", TP53 == "WT") %>%
  get_cox_statistics()


# surv_data_EGFR <-
#   survival::Surv(time = filter(processed_rna_os_data, gene == "EGFR")$os_months, filter(processed_rna_os_data, gene == "EGFR")$os_status) ~ filter(processed_rna_os_data, gene == "EGFR")$AURKA_rna_exp)
#
# surv_data_KRAS <-
#   survival::Surv(time = os_months, os_status) ~ AURKA_rna_exp, data = filter(processed_rna_os_data, gene == "KRAS"))

egfr_data_fit <-
  surv_fit(Surv(time = os_months, os_status) ~ TP53, data = egfr_data)
kras_data_fit <-
  surv_fit(Surv(time = os_months, os_status) ~ TP53, data = kras_data)

ggsurvplot(
  egfr_data_fit,
  # survfit object with calculated statistics.
  data = egfr_data,
  risk.table = T,
  # data used to fit survival curves.
  pval = T,
  # show p-value of log-rank test.
  conf.int = TRUE,
  # show confidence intervals for
  # point estimates of survival curves.
  xlim = c(0, 150),
  # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",
  # customize X axis label.
  break.time.by = 50,
  # break X axis in time intervals by 500.
  ggtheme = theme_light() # customize plot and risk table with a theme.
  # in legend of risk table
)

ggsurvplot(
  kras_data_fit,
  # survfit object with calculated statistics.
  data = egfr_data,
  risk.table = T,
  # data used to fit survival curves.
  pval = T,
  # show p-value of log-rank test.
  conf.int = TRUE,
  # show confidence intervals for
  # point estimates of survival curves.
  xlim = c(0, 150),
  # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in months",
  # customize X axis label.
  break.time.by = 50,
  # break X axis in time intervals by 500.
  ggtheme = theme_light() # customize plot and risk table with a theme.
  # in legend of risk table
)


# Figure XA EGFR version -------------------------------------------------------
processed_rna_data <- merged_data %>%
  filter(EGFR != KRAS) %>%
  select(EGFR_rna_exp, EGFR, KRAS) %>%
  pivot_longer(
    cols = c("EGFR", "KRAS"),
    names_to = "gene",
    values_to = "impact"
  ) %>%
  filter(impact != "WT") %>%
  select(-impact) %>%
  mutate(gene = as.factor(gene))

processed_rna_data %>%
  ggviolin(
    y = "EGFR_rna_exp",
    x = "gene",
    fill = "gene",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test",
                     label.y = 4.5,
                     label.x = 1.3)

ggsave(
  "Paper Figures/Fig XA EGFR.png",
  dpi = 400,
  height = 5,
  width = 5,
  units = "in"
)

# Statistics
EGFR_only_EGFR_rna_exp <-
  filter(processed_rna_data, gene == "EGFR")$EGFR_rna_exp
KRAS_only_EGFR_rna_exp <-
  filter(processed_rna_data, gene == "KRAS")$EGFR_rna_exp

# Check for normality
ad.test(EGFR_only_EGFR_rna_exp)
ad.test(KRAS_only_EGFR_rna_exp)

# Compute variance for t.test
var(EGFR_only_EGFR_rna_exp)
# 0.9
var(KRAS_only_EGFR_rna_exp)
# 1

ks.test(EGFR_only_EGFR_rna_exp, KRAS_only_EGFR_rna_exp)
t.test(EGFR_only_EGFR_rna_exp, KRAS_only_EGFR_rna_exp)

# Figure XA EGFR protein level -------------------------------------------------
processed_protein_data <- merged_data %>%
  filter(EGFR != KRAS) %>%
  select(EGFR_protein_exp, EGFR, KRAS) %>%
  pivot_longer(
    cols = c("EGFR", "KRAS"),
    names_to = "gene",
    values_to = "impact"
  ) %>%
  filter(impact != "WT") %>%
  select(-impact) %>%
  mutate(gene = as.factor(gene))

processed_protein_data %>%
  ggviolin(
    y = "EGFR_protein_exp",
    x = "gene",
    fill = "gene",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test",
                     label.y = 4.5,
                     label.x = 1.3)

ggsave(
  "Paper Figures/Fig XA EGFR protein level.png",
  dpi = 400,
  height = 5,
  width = 5,
  units = "in"
)

# Statistics
EGFR_only_EGFR_protein_exp <-
  filter(processed_protein_data, gene == "EGFR")$EGFR_protein_exp
KRAS_only_EGFR_protein_exp <-
  filter(processed_protein_data, gene == "KRAS")$EGFR_protein_exp

# Check for normality
ad.test(EGFR_only_EGFR_protein_exp)
ad.test(KRAS_only_EGFR_protein_exp)

# Compute variance for t.test
var(EGFR_only_EGFR_protein_exp, na.rm = T)
# 0.9
var(KRAS_only_EGFR_protein_exp, na.rm = T)
# 1

ks.test(EGFR_only_EGFR_protein_exp, KRAS_only_EGFR_protein_exp)
t.test(EGFR_only_EGFR_protein_exp, KRAS_only_EGFR_protein_exp)

processed_cna_data <- merged_data %>%
  filter(EGFR != KRAS, !is.na(AURKA_CNA_log2)) %>%
  select(AURKA_CNA_log2, EGFR, KRAS) %>%
  pivot_longer(
    cols = c("EGFR", "KRAS"),
    names_to = "gene",
    values_to = "impact"
  ) %>%
  filter(impact != "WT") %>%
  select(-impact) %>%
  mutate(gene = as.factor(gene))

processed_cna_data %>%
  ggviolin(
    y = "AURKA_CNA_log2",
    x = "gene",
    fill = "gene",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox.test",
                     label.y = 4.5,
                     label.x = 1.3)

ggsave(
  "Paper Figures/Fig XB.png",
  dpi = 400,
  height = 5,
  width = 5,
  units = "in"
)

# Statistics
EGFR_only_AURKA_CNA_log2 <-
  filter(processed_cna_data, gene == "EGFR")$AURKA_CNA_log2
KRAS_only_AURKA_CNA_log2 <-
  filter(processed_cna_data, gene == "KRAS")$AURKA_CNA_log2

# Check for normality
ad.test(EGFR_only_AURKA_CNA_log2)
ad.test(KRAS_only_AURKA_CNA_log2)

# Compute variance for t.test
var(EGFR_only_AURKA_CNA_log2)
# 0.276
var(KRAS_only_AURKA_CNA_log2)
# 0.164

ks.test(EGFR_only_AURKA_CNA_log2, KRAS_only_AURKA_CNA_log2)
t.test(EGFR_only_AURKA_CNA_log2, KRAS_only_AURKA_CNA_log2)
# Figure XB alternative --------------------------------------------------------
processed_cna_data_discrete <- merged_data %>%
  filter(EGFR != KRAS, !is.na(AURKA_cna)) %>%
  select(AURKA_cna, EGFR, KRAS) %>%
  pivot_longer(
    cols = c("EGFR", "KRAS"),
    names_to = "gene",
    values_to = "impact"
  ) %>%
  filter(impact != "WT") %>%
  select(-impact) %>%
  mutate(gene = as.factor(gene),
         AURKA_cna = as.factor(AURKA_cna)) %>%
  group_by(AURKA_cna, gene) %>%
  summarise(count = n()) %>%
  mutate(
    AURKA = fct_recode(
      AURKA_cna,
      "Shallow Deletion" = "-1",
      "Diploid" = "0",
      "Gain" = "1",
      "Amplification" = "2"
    )
  ) %>%
  ungroup() %>%
  group_by(gene) %>%
  reframe(proportion = count / sum(count),
          AURKA = AURKA,
          count = count)

processed_cna_data_discrete %>%
  mutate(AURKA = fct_relevel(AURKA, "Amplification", "Gain", "Diploid", "Shallow Deletion")) %>%
  ggplot(aes(
    x = gene,
    y = proportion,
    fill = AURKA,
    label = count
  )) +
  geom_col(position = "fill",
           width = .5,
           colour = "black") +
  theme_minimal() +
  xlab("Genomic Variation") +
  ylab("Relative Frequency") +
  geom_text(position = position_stack(vjust = 0.5), size = 4) +
  coord_flip() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    legend.text = element_text(size = 15, colour = "black")
  ) +
  guides(fill = guide_legend(reverse = TRUE))

ggsave(
  "Paper Figures/Fig XB alternative.png",
  dpi = 400,
  height = 5,
  width = 10,
  units = "in"
)

# Statistics
EGFR_only_AURKA_CNA <-
  filter(processed_cna_data_discrete, gene == "EGFR")$count
KRAS_only_AURKA_CNA <-
  filter(processed_cna_data_discrete, gene == "KRAS")$count

fisher.test(matrix(
  c(EGFR_only_AURKA_CNA, KRAS_only_AURKA_CNA),
  byrow = T,
  nrow = 2
))

# Run a proportions test for each pair of proportions
for (i in 1:4) {
  test_result <-
    prop.test(c(EGFR_only_AURKA_CNA[i], KRAS_only_AURKA_CNA[i]),
              c(sum(EGFR_only_AURKA_CNA), sum(KRAS_only_AURKA_CNA)))
  print(test_result)
}

# Figure XB alternative TP53 ---------------------------------------------------
processed_cna_data_discrete_tp53 <- merged_data %>%
  filter(EGFR != KRAS, !is.na(AURKA_cna)) %>%
  select(AURKA_cna, EGFR, KRAS, TP53) %>%
  pivot_longer(
    cols = c("EGFR", "KRAS"),
    names_to = "gene",
    values_to = "impact"
  ) %>%
  filter(impact != "WT") %>%
  select(-impact) %>%
  mutate(
    gene = paste("TP53", TP53, gene, "MUT", sep = "_"),
    gene = as.factor(gene),
    AURKA_cna = as.factor(AURKA_cna)
  ) %>%
  select(-TP53) %>%
  group_by(AURKA_cna, gene) %>%
  summarise(count = n()) %>%
  mutate(
    AURKA = fct_recode(
      AURKA_cna,
      "Shallow Deletion" = "-1",
      "Diploid" = "0",
      "Gain" = "1",
      "Amplification" = "2"
    )
  ) %>%
  ungroup() %>%
  group_by(gene) %>%
  reframe(proportion = count / sum(count),
          AURKA = AURKA,
          count = count)

processed_cna_data_discrete_tp53 %>%
  mutate(AURKA = fct_relevel(AURKA, "Amplification", "Gain", "Diploid", "Shallow Deletion")) %>%
  ggplot(aes(
    x = gene,
    y = proportion,
    fill = AURKA,
    label = count
  )) +
  geom_col(position = "fill",
           width = .5,
           colour = "black") +
  theme_minimal() +
  xlab("Genomic Variation") +
  ylab("Relative Frequency") +
  geom_text(position = position_stack(vjust = 0.5), size = 4) +
  coord_flip() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 15, colour = "black"),
    legend.text = element_text(size = 15, colour = "black")
  ) +
  guides(fill = guide_legend(reverse = TRUE))

ggsave(
  "Paper Figures/Fig XB alternative with TP53.png",
  dpi = 400,
  height = 5,
  width = 13,
  units = "in"
)

# Statistics
EGFR_only_AURKA_CNA <-
  filter(processed_cna_data_discrete_tp53, gene == "EGFR")$count
KRAS_only_AURKA_CNA <-
  filter(processed_cna_data_discrete_tp53, gene == "KRAS")$count

fisher.test(matrix(
  c(EGFR_only_AURKA_CNA, KRAS_only_AURKA_CNA),
  byrow = T,
  nrow = 2
))

# Run a proportions test for each pair of proportions
# test_df <- processed_cna_data_discrete_tp53 %>%
#   select(-proportion) %>%
#   pivot_wider(names_from = gene, values_from = count)
#
# for (i in 1:4) {
#   for (j in 1:3) {
#     test_result <- prop.test(c(as.numeric(test_df[i, j + 1]), as.numeric(test_df[i, j + 2])),
#                              c(sum(test_df[, j + 1]), sum(test_df[, j + 2])))
#     if (test_result$p.value < 0.05) {
#       print(test_df[i, j])
#       print(test_df[i, j + 1])
#       print(test_result)
#     }
#   }
# }



# NEW way
statistics <- processed_cna_data_discrete_tp53 %>%
  mutate(AURKA = str_replace(AURKA, "Amplification", "Gain")) %>%
  mutate(AURKA = as.factor(AURKA)) %>%
  group_by(gene, AURKA) %>%
  reframe(gene = gene, count = sum(count)) %>%
  distinct() %>%
  group_by(gene) %>%
  reframe(sum = sum(count),
          count = count,
          AURKA = AURKA) %>%
  group_by(AURKA) %>%
  filter(AURKA == "Gain")

i_ext <- c()
for (i in 1:3) {
  i_ext <- c(1:4)[which(c(1:4) > i)]
  for (k in i_ext) {
    print(c(i, k))
    print(prop.test(
      c(statistics$count[i], statistics$count[k]),
      c(statistics$sum[i], statistics$sum[k])
    )$p.value)
  }
}

# Figure XF --------------------------------------------------------------------
processed_rna_data_tp53 <- merged_data %>%
  filter(EGFR != KRAS) %>%
  select(EGFR, TP53, KRAS, EGFR_rna_exp) %>%
  mutate(
    EGFR_only = ifelse(EGFR == "ALT" & KRAS == "WT", "ALT", "WT"),
    KRAS_only = ifelse(KRAS == "ALT" & EGFR == "WT", "ALT", "WT"),
    across(
      !EGFR_rna_exp,
      .fns = function(x) {
        as.factor(x)
      }
    )
  ) %>% 
  pivot_longer(cols = contains("only"),
               names_to = "subgroup",
               values_to = "alteration") %>%
  filter(!alteration == "WT") 

processed_rna_data_tp53 %>%
  ggplot(aes(
    y = EGFR_rna_exp,
    fill = subgroup,
    x = subgroup,
    col = TP53
  )) +
  geom_violin(show.legend = F) +
  geom_boxplot(width = .2, position = position_dodge(.9)) +
  scale_color_manual(values = c("ALT" = "gray80", "WT" = "black")) +
  theme_classic()

my_comparisons <- list(c("EGFR_only", "KRAS_only"))

processed_rna_data_tp53 %>%
  ggviolin(
    y = "EGFR_rna_exp",
    fill = "subgroup",
    x = "subgroup",
    col = "TP53",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
  theme(legend.position = "none") +
  # stat_compare_means(method = "wilcox.test",
  #                    label.y = 4.5,
  #                    label.x = 1.3) +
  scale_color_manual(values = c("ALT" = "gray80", "WT" = "black"))


ggsave(
  "Paper Figures/Fig XF.png",
  dpi = 400,
  height = 5,
  width = 5,
  units = "in"
)


# Statistics
EGFR_only_TP53_mut_EGFR_rna_exp <-
  filter(processed_rna_data_tp53, EGFR == "ALT", TP53 == "ALT")$EGFR_rna_exp
EGFR_only_TP53_wt_EGFR_rna_exp <-
  filter(processed_rna_data_tp53, EGFR == "ALT", TP53 == "WT")$EGFR_rna_exp

KRAS_only_TP53_mut_EGFR_rna_exp <-
  filter(processed_rna_data_tp53, KRAS == "ALT", TP53 == "ALT")$EGFR_rna_exp
KRAS_only_TP53_wt_EGFR_rna_exp <-
  filter(processed_rna_data_tp53, KRAS == "ALT", TP53 == "WT")$EGFR_rna_exp


# Check for normality
ad.test(EGFR_only_TP53_mut_EGFR_rna_exp)
ad.test(EGFR_only_TP53_wt_EGFR_rna_exp)
ad.test(KRAS_only_TP53_mut_EGFR_rna_exp)
ad.test(KRAS_only_TP53_wt_EGFR_rna_exp)

# Compute variance for t.test
var(EGFR_only_TP53_mut_EGFR_rna_exp)
# 0.88
var(EGFR_only_TP53_wt_EGFR_rna_exp)
# 0.71
var(KRAS_only_TP53_mut_EGFR_rna_exp)
# 0.73
var(KRAS_only_TP53_wt_EGFR_rna_exp)
# 1.03

ks.test(EGFR_only_TP53_mut_EGFR_rna_exp,
        EGFR_only_TP53_wt_EGFR_rna_exp)
ks.test(EGFR_only_TP53_wt_EGFR_rna_exp,
        KRAS_only_TP53_wt_EGFR_rna_exp)
ks.test(KRAS_only_TP53_mut_EGFR_rna_exp,
        KRAS_only_TP53_wt_EGFR_rna_exp)
ks.test(EGFR_only_TP53_mut_EGFR_rna_exp,
        KRAS_only_TP53_mut_EGFR_rna_exp)

t.test(EGFR_only_TP53_mut_EGFR_rna_exp,
       EGFR_only_TP53_wt_EGFR_rna_exp)
t.test(EGFR_only_TP53_wt_EGFR_rna_exp,
       KRAS_only_TP53_wt_EGFR_rna_exp)
t.test(KRAS_only_TP53_mut_EGFR_rna_exp,
       KRAS_only_TP53_wt_EGFR_rna_exp)
t.test(EGFR_only_TP53_mut_EGFR_rna_exp,
       KRAS_only_TP53_mut_EGFR_rna_exp)
