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
merged_data <- fread("Merged annotated data AURKA, KRAS, TP53, EGFR.csv") %>% 
  as_tibble()

# Figure XA --------------------------------------------------------------------
merged_data %>% 
  filter(EGFR != KRAS) %>% 
  select(AURKA_rna_exp, EGFR, KRAS) %>% 
  pivot_longer(cols = c("EGFR", "KRAS"), names_to = "gene", values_to = "impact") %>% 
  filter(impact != "WT") %>% 
  select(-impact) %>% 
  mutate(gene = as.factor(gene)) %>%
  ggplot(aes(y = AURKA_rna_exp, x = gene, fill = gene)) +
  geom_violin(show.legend = F) +
  geom_boxplot(width = .2, show.legend = F) +
  theme_classic()

ggsave("Paper Figures/Fig XA.png", dpi = 400, height = 5, width = 5, units = "in")

# Statistics
EGFR_only_AURKA_rna_exp = filter(test, gene == "EGFR")$AURKA_rna_exp
KRAS_only_AURKA_rna_exp = filter(test, gene == "KRAS")$AURKA_rna_exp

# Check for normality
nortest::ad.test(EGFR_only_AURKA_rna_exp)
nortest::ad.test(KRAS_only_AURKA_rna_exp)

# Compute variance for t.test
var(EGFR_only_AURKA_rna_exp)
# 0.9
var(KRAS_only_AURKA_rna_exp)
# 1

ks.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)
t.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)

