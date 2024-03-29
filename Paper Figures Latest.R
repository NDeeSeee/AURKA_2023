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

# Figure A -------------------------------------------------------
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
  ylab("EGFR mRNA") +
  xlab("") +
  scale_x_discrete(labels = c("EGFR" = "EGFR MUT", "KRAS" = "KRAS MUT"))
# stat_compare_means(method = "wilcox.test",
#                    label.y = 4.5,
#                    label.x = 1.3)

ggsave(
  "Paper Figures/Fig 1A.png",
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
# 1
var(KRAS_only_EGFR_rna_exp)
# 0.64

# Compute median
median(EGFR_only_EGFR_rna_exp)
# 0.53
median(KRAS_only_EGFR_rna_exp)
# -0.14

# Compute means and n
mean(EGFR_only_EGFR_rna_exp)
length(EGFR_only_EGFR_rna_exp)
mean(KRAS_only_EGFR_rna_exp)
length(KRAS_only_EGFR_rna_exp)

summary(EGFR_only_EGFR_rna_exp)
summary(KRAS_only_EGFR_rna_exp)


ks.test(EGFR_only_EGFR_rna_exp, KRAS_only_EGFR_rna_exp)
t.test(EGFR_only_EGFR_rna_exp, KRAS_only_EGFR_rna_exp)
wilcox.test(EGFR_only_EGFR_rna_exp, KRAS_only_EGFR_rna_exp)


# Figure B -------------------------------------------------
processed_protein_data <- merged_data %>%
  filter(EGFR != KRAS, study_id != "luad_cptac_2020") %>%
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
  ylab("EGFR protein") +
  xlab("") +
  scale_x_discrete(labels = c("EGFR" = "EGFR MUT", "KRAS" = "KRAS MUT"))
# stat_compare_means(method = "wilcox.test",
#                    label.y = 4.5,
#                    label.x = 1.3)

ggsave(
  "Paper Figures/Fig 1B.png",
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
# 9.424e-05
ad.test(KRAS_only_EGFR_protein_exp)
# 0.0008593

# Compute variance for t.test
var(EGFR_only_EGFR_protein_exp, na.rm = T)
# 2.35
var(KRAS_only_EGFR_protein_exp, na.rm = T)
# 0.57

# Compute medians
median(EGFR_only_EGFR_protein_exp, na.rm = T)
# 0.48
median(KRAS_only_EGFR_protein_exp, na.rm = T)
# -0.3

mean(EGFR_only_EGFR_protein_exp, na.rm = T)
length(EGFR_only_EGFR_protein_exp)
mean(KRAS_only_EGFR_protein_exp, na.rm = T)
length(KRAS_only_EGFR_protein_exp)

summary(EGFR_only_EGFR_protein_exp)
summary(KRAS_only_EGFR_protein_exp)

ks.test(EGFR_only_EGFR_protein_exp, KRAS_only_EGFR_protein_exp)
t.test(EGFR_only_EGFR_protein_exp, KRAS_only_EGFR_protein_exp)
wilcox.test(EGFR_only_EGFR_protein_exp, KRAS_only_EGFR_protein_exp)

# Figure C --------------------------------------------------------------------
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
  ylab("AURKA mRNA") +
  xlab("") +
  scale_x_discrete(labels = c("EGFR" = "EGFR MUT", "KRAS" = "KRAS MUT"))
# stat_compare_means(method = "wilcox.test",
#                    label.y = 4.5,
#                    label.x = 1.3)

ggsave(
  "Paper Figures/Fig 1C.png",
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
# 0.1635
ad.test(KRAS_only_AURKA_rna_exp)
# 0.03213

# Compute variance for t.test
var(EGFR_only_AURKA_rna_exp)
# 0.9
var(KRAS_only_AURKA_rna_exp)
# 1

# Medians
median(EGFR_only_AURKA_rna_exp)
# -0.16
median(KRAS_only_AURKA_rna_exp)
# 0.18

mean(EGFR_only_AURKA_rna_exp, na.rm = T)
length(EGFR_only_AURKA_rna_exp)
mean(KRAS_only_AURKA_rna_exp, na.rm = T)
length(KRAS_only_AURKA_rna_exp)

summary(EGFR_only_AURKA_rna_exp)
summary(KRAS_only_AURKA_rna_exp)

ks.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)
t.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)
wilcox.test(EGFR_only_AURKA_rna_exp, KRAS_only_AURKA_rna_exp)

# Figure D ---------------------------------------------------
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
  geom_text(position = position_stack(vjust = 0.5), size = 6) +
  coord_flip() +
  scale_fill_brewer(palette = "Pastel1") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 20, colour = "black"),
    legend.text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE))

ggsave(
  "Paper Figures/Fig 1D.png",
  dpi = 400,
  height = 5,
  width = 13,
  units = "in"
)

# Statistics
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


# Figure E --------------------------------------------------------------------
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

# processed_rna_data_tp53 %>%
#   ggplot(aes(
#     y = AURKA_rna_exp,
#     fill = subgroup,
#     x = subgroup,
#     col = TP53
#   )) +
#   geom_violin(show.legend = F) +
#   geom_boxplot(width = .2, position = position_dodge(.9)) +
#   scale_color_manual(values = c("ALT" = "gray80", "WT" = "black")) +
#   theme_classic()

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
  scale_color_manual(values = c("ALT" = "gray80", "WT" = "black")) +
  ylab("AURKA mRNA") +
  xlab("") +
  scale_x_discrete(labels = c("EGFR_only" = "EGFR MUT", "KRAS_only" = "KRAS MUT"))

ggsave(
  "Paper Figures/Fig 1E.png",
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

wilcox.test(EGFR_only_TP53_mut_AURKA_rna_exp,
       EGFR_only_TP53_wt_AURKA_rna_exp)
wilcox.test(EGFR_only_TP53_wt_AURKA_rna_exp,
       KRAS_only_TP53_wt_AURKA_rna_exp)
wilcox.test(KRAS_only_TP53_mut_AURKA_rna_exp,
       KRAS_only_TP53_wt_AURKA_rna_exp)
wilcox.test(EGFR_only_TP53_mut_AURKA_rna_exp,
       KRAS_only_TP53_mut_AURKA_rna_exp)

mean(EGFR_only_TP53_mut_AURKA_rna_exp, na.rm = T)
length(EGFR_only_TP53_mut_AURKA_rna_exp)
mean(EGFR_only_TP53_wt_AURKA_rna_exp, na.rm = T)
length(EGFR_only_TP53_wt_AURKA_rna_exp)
mean(KRAS_only_TP53_mut_AURKA_rna_exp, na.rm = T)
length(KRAS_only_TP53_mut_AURKA_rna_exp)
mean(KRAS_only_TP53_wt_AURKA_rna_exp, na.rm = T)
length(KRAS_only_TP53_wt_AURKA_rna_exp)

summary(EGFR_only_TP53_mut_AURKA_rna_exp)
summary(EGFR_only_TP53_wt_AURKA_rna_exp)
summary(KRAS_only_TP53_mut_AURKA_rna_exp)
summary(KRAS_only_TP53_wt_AURKA_rna_exp)


# Figure F --------------------------------------------------------------------
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

# processed_rna_data_tp53 %>%
#   ggplot(aes(
#     y = EGFR_rna_exp,
#     fill = subgroup,
#     x = subgroup,
#     col = TP53
#   )) +
#   geom_violin(show.legend = F) +
#   geom_boxplot(width = .2, position = position_dodge(.9)) +
#   scale_color_manual(values = c("ALT" = "gray80", "WT" = "black")) +
#   theme_classic()

# my_comparisons <- list(c("EGFR_only", "KRAS_only"))

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
  scale_color_manual(values = c("ALT" = "gray80", "WT" = "black")) +
  ylab("EGFR mRNA") +
  xlab("") +
  scale_x_discrete(labels = c("EGFR_only" = "EGFR MUT", "KRAS_only" = "KRAS MUT"))


ggsave(
  "Paper Figures/Fig 1F.png",
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
# 1.4
var(EGFR_only_TP53_wt_EGFR_rna_exp)
# 0.67
var(KRAS_only_TP53_mut_EGFR_rna_exp)
# 0.79
var(KRAS_only_TP53_wt_EGFR_rna_exp)
# 0.51


median(EGFR_only_TP53_mut_EGFR_rna_exp)
# 0.72
median(EGFR_only_TP53_wt_EGFR_rna_exp)
# 0.37
median(KRAS_only_TP53_mut_EGFR_rna_exp)
# 0.05
median(KRAS_only_TP53_wt_EGFR_rna_exp)
# -0.31

ks.test(EGFR_only_TP53_mut_EGFR_rna_exp,
        EGFR_only_TP53_wt_EGFR_rna_exp)
ks.test(EGFR_only_TP53_wt_EGFR_rna_exp,
        KRAS_only_TP53_wt_EGFR_rna_exp)
ks.test(KRAS_only_TP53_mut_EGFR_rna_exp,
        KRAS_only_TP53_wt_EGFR_rna_exp)
ks.test(EGFR_only_TP53_mut_EGFR_rna_exp,
        KRAS_only_TP53_mut_EGFR_rna_exp)

wilcox.test(EGFR_only_TP53_mut_EGFR_rna_exp,
       EGFR_only_TP53_wt_EGFR_rna_exp)
wilcox.test(EGFR_only_TP53_wt_EGFR_rna_exp,
       KRAS_only_TP53_wt_EGFR_rna_exp)
wilcox.test(KRAS_only_TP53_mut_EGFR_rna_exp,
       KRAS_only_TP53_wt_EGFR_rna_exp)
wilcox.test(EGFR_only_TP53_mut_EGFR_rna_exp,
       KRAS_only_TP53_mut_EGFR_rna_exp)

mean(EGFR_only_TP53_mut_EGFR_rna_exp, na.rm = T)
length(EGFR_only_TP53_mut_EGFR_rna_exp)
mean(EGFR_only_TP53_wt_EGFR_rna_exp, na.rm = T)
length(EGFR_only_TP53_wt_EGFR_rna_exp)
mean(KRAS_only_TP53_mut_EGFR_rna_exp, na.rm = T)
length(KRAS_only_TP53_mut_EGFR_rna_exp)
mean(KRAS_only_TP53_wt_EGFR_rna_exp, na.rm = T)
length(KRAS_only_TP53_wt_EGFR_rna_exp)



summary(EGFR_only_TP53_mut_EGFR_rna_exp)
summary(EGFR_only_TP53_wt_EGFR_rna_exp)
summary(KRAS_only_TP53_mut_EGFR_rna_exp)
summary(KRAS_only_TP53_wt_EGFR_rna_exp)
