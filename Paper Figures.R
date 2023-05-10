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
    "nortest"
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
  ggplot(aes(y = AURKA_rna_exp, x = gene, fill = gene)) +
  geom_violin(show.legend = F) +
  geom_boxplot(width = .2, show.legend = F) +
  theme_classic()

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
  ggplot(aes(y = AURKA_CNA_log2, x = gene, fill = gene)) +
  geom_violin(show.legend = F) +
  geom_boxplot(width = .2, show.legend = F) +
  theme_classic()

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
  theme_minimal()

ggsave("Paper Figures/Fig XC.png", dpi = 400, height = 5, width = 5, units = "in")

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

ks.test(EGFR_only_TP53_mut_AURKA_rna_exp, EGFR_only_TP53_wt_AURKA_rna_exp)
ks.test(EGFR_only_TP53_wt_AURKA_rna_exp, KRAS_only_TP53_wt_AURKA_rna_exp)
ks.test(KRAS_only_TP53_mut_AURKA_rna_exp, KRAS_only_TP53_wt_AURKA_rna_exp)
ks.test(EGFR_only_TP53_mut_AURKA_rna_exp, KRAS_only_TP53_mut_AURKA_rna_exp)

t.test(EGFR_only_TP53_mut_AURKA_rna_exp, EGFR_only_TP53_wt_AURKA_rna_exp)
t.test(EGFR_only_TP53_wt_AURKA_rna_exp, KRAS_only_TP53_wt_AURKA_rna_exp)
t.test(KRAS_only_TP53_mut_AURKA_rna_exp, KRAS_only_TP53_wt_AURKA_rna_exp)
t.test(EGFR_only_TP53_mut_AURKA_rna_exp, KRAS_only_TP53_mut_AURKA_rna_exp)
