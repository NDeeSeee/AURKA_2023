# Comparing EGFR and AURKA expression ------------------------------------------
# 4 samples only for BOTH (KRAS+EGFR) mutation

KRAS_wt_EGFR_wt <-
  fread("~/Downloads/mRNA expression KRAS wt EGFR wt.txt") %>%
  select(-STUDY_ID) %>%
  mutate(status = "no one")

KRAS_wt_EGFR_mut <-
  fread("~/Downloads/mRNA expression KRAS wt EGFR mut.txt") %>%
  select(-STUDY_ID) %>%
  mutate(status = "EGFR only")

KRAS_mut_EGFR_wt <-
  fread("~/Downloads/mRNA expression KRAS mut EGFR wt.txt") %>%
  select(-STUDY_ID) %>%
  mutate(status = "KRAS only")


all_data <-
  rbind(KRAS_wt_EGFR_wt, KRAS_wt_EGFR_mut, KRAS_mut_EGFR_wt)

ggplot(all_data, aes(x = as.factor(status), y = AURKA)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, color = "red")

ggplot(all_data, aes(x = as.factor(status), y = EGFR)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, color = "red")

# KRAS vs non-KRAS
KRAS_mut <-
  fread(
    "~/Downloads/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM) KRAS mut only.txt"
  ) %>%
  select(-STUDY_ID) %>%
  mutate(status = "KRAS mut")

KRAS_wt <-
  fread(
    "~/Downloads/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM) KRAS wt.txt"
  ) %>%
  select(-STUDY_ID) %>%
  mutate(status = "KRAS wt")

KRAS_data <- rbind(KRAS_mut, KRAS_wt)

ggplot(KRAS_data, aes(x = as.factor(status), y = AURKA)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, color = "red")

ggplot(KRAS_data, aes(x = as.factor(status), y = EGFR)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, color = "red")
