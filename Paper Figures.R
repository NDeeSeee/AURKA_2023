# AURKA paper figures

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

ggsave("Fig XA.png", dpi = 400, height = 5, width = 5, units = "in")

# Statistics
EGFR_only_AURKA_rna_exp = filter(test, gene == "EGFR")$AURKA_rna_exp
KRAS_only_AURKA_rna_exp = filter(test, gene == "KRAS")$AURKA_rna_exp

# Check for normality
install.packages("nortest")
nortest::ad.test(EGFR_only_AURKA_rna_exp)
nortest::ad.test(KRAS_only_AURKA_rna_exp)

# Compute variance for t.test
var(EGFR_only_AURKA_rna_exp)
# 0.9
var(KRAS_only_AURKA_rna_exp)
# 1

ks.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)
t.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)

