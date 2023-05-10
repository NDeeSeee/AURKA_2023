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
  )

processed_rna_data_tp53 %>%
  pivot_longer(cols = contains("only"),
               names_to = "subgroup",
               values_to = "alteration") %>%
  filter(!alteration == "WT") %>%
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
  pivot_longer(cols = contains("only"),
               names_to = "subgroup",
               values_to = "alteration") %>%
  filter(!alteration == "WT") %>%
  ggviolin(
    y = "AURKA_rna_exp",
    x = "subgroup",
    col = "TP53",
    fill = "subgroup",
    add = "boxplot",
    add.params = list(fill = "white")
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label.y = 4.5,
    label.x = 1.3,
    paired = T
  )

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

egfr_data <- filter(processed_rna_os_data, gene == "EGFR") %>%
  select(-gene)
kras_data <- filter(processed_rna_os_data, gene == "KRAS") %>%
  select(-gene)

cox_model_egfr <-
  coxph(Surv(time = os_months, os_status) ~ AURKA_rna_exp, data = egfr_data)
cox_model_kras <-
  coxph(Surv(time = os_months, os_status) ~ AURKA_rna_exp, data = kras_data)

visualise_cox_model <- function(cox_model_data) {
  # Assuming cox_model is your fitted model
  # Use broom to tidy the results
  tidy_cox <- tidy(cox_model_data, conf.int = T)
  
  # Create data frame for forest plot
  fp_data <- tidy_cox %>%
    select(term, estimate, conf.low, conf.high) %>%
    mutate(term = as.character(term)) %>%
    rename(
      `Variables` = term,
      `Hazard Ratio` = estimate,
      `Lower CI` = conf.low,
      `Upper CI` = conf.high
    ) %>%
    as.data.frame() %>%
    clean_names()
  
  return(fp_data)
}

cox_model_egfr <- visualise_cox_model(cox_model_egfr)
cox_model_kras <- visualise_cox_model(cox_model_kras)

cox_model_egfr %>%
  mutate(variables = "EGFR") %>%
  bind_rows(mutate(cox_model_kras, variables = "KRAS")) %>%
  as_tibble() %>%
  mutate(variables = as.factor(variables)) %>%
  ggplot(aes(x = variables, y = hazard_ratio, col = variables)) +
  geom_hline(yintercept = 0) +
  geom_point(position = position_dodge(0.05),
             size = 2.5,
             show.legend = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = .2,
    position = position_dodge(0.05),
    size = 1,
    show.legend = F
  ) +
  coord_cartesian(ylim = c(-3, 4.5)) +
  theme_classic()


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
