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
    "janitor"
  )

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename <- dplyr::rename
select <- dplyr::select
filter <- dplyr::filter
group_by <- dplyr::group_by
mutate <- dplyr::mutate

# Functions --------------------------------------------------------------------
parse_cbioportal_data <-
  function(sample_data,
           patient_data,
           alterations_data,
           rna_seq_data,
           cna_data,
           custom_columns) {
    custom_columns <- custom_columns
    
    sample_data <- fread(sample_data,
                         skip = 4,
                         na.strings = c("", "[Not Available]")) %>%
      as_tibble()
    
    sample_data <- sample_data %>%
      select(names(sample_data)[which(names(sample_data) %in% names(custom_columns))])
    
    names(sample_data) <- custom_columns[names(sample_data)]
    
    alterations_data <-
      fread(alterations_data, na.strings = c("", "[Not Available]")) %>%
      as_tibble() %>%
      mutate(sample_id = `Sample ID`, study_id = `Study ID`) %>%
      select(sample_id, study_id, TP53, EGFR, KRAS, AURKA) %>%
      mutate(across(
        .cols = matches("TP53|EGFR|KRAS|AURKA"),
        .fns = function(x) {
          ifelse(x == "no alteration", "WT", "ALT")
        }
      ))
    
    rna_seq_data <-
      fread(rna_seq_data, na.strings = c("", "[Not Available]")) %>%
      as_tibble() %>%
      filter(!if_any(
        .cols = everything(),
        .fns = function(x) {
          is.na(x)
        }
      )) %>%
      rename(
        sample_id = SAMPLE_ID,
        study_id = STUDY_ID,
        TP53_rna_exp = TP53,
        EGFR_rna_exp = EGFR,
        KRAS_rna_exp = KRAS,
        AURKA_rna_exp = AURKA
      ) %>%
      select(-study_id)
    
    cna_data <-
      fread(cna_data, na.strings = c("", "[Not Available]")) %>%
      as_tibble() %>%
      filter(!if_any(
        .cols = everything(),
        .fns = function(x) {
          is.na(x)
        }
      )) %>%
      rename(
        sample_id = SAMPLE_ID,
        study_id = STUDY_ID,
        TP53_cna = TP53,
        EGFR_cna = EGFR,
        KRAS_cna = KRAS,
        AURKA_cna = AURKA
      ) %>%
      select(-study_id)
    
    # protein_exp_data <- fread(protein_exp_data) %>%
    #   as_tibble() %>%
    #   filter(!if_all(
    #     .cols = c("TP53", "EGFR", "KRAS", "AURKA"),
    #     .fns = function(x)
    #       is.na(x)
    #   ))  %>%
    #   rename(sample_id = SAMPLE_ID, study_id = STUDY_ID,
    #          TP53_protein_exp = TP53,
    #          EGFR_protein_exp = EGFR,
    #          KRAS_protein_exp = KRAS,
    #          AURKA_protein_exp = AURKA) %>%
    #   select(-study_id)
    
    clinical_data <- fread(patient_data,
                           skip = 4,
                           na.strings = c("", "[Not Available]")) %>%
      as_tibble()
    
    clinical_data <- clinical_data %>%
      select(names(clinical_data)[which(names(clinical_data) %in% names(custom_columns))])
    
    names(clinical_data) <- custom_columns[names(clinical_data)]
    
    clinical_data <- clinical_data %>%
      left_join(sample_data) %>%
      left_join(alterations_data) %>%
      left_join(rna_seq_data) %>%
      left_join(cna_data) %>%
      # left_join(protein_exp_data)
      mutate(
        across(
          .cols = matches("exp|tmb|age|months|smoking_pack_years"),
          .fns = function(x) {
            as.numeric(x)
          }
        ),
        across(
          .cols = matches("sex"),
          .fns = function(x) {
            tolower(x)
          }
        ),
        across(
          .cols = matches("smoking_history"),
          .fns = function(x) {
            as.character(x)
          }
        )
      )
    
    return(clinical_data)
  }

custom_columns <- c(
  "PATIENT_ID" = "patient_id",
  "SAMPLE_ID" = "sample_id",
  "TMB_NONSYNONYMOUS" = "tmb",
  "TOBACCO_SMOKING_HISTORY_INDICATOR" = "smoking_history",
  "SMOKING_HISTORY" = "smoking_history",
  "AJCC_PATHOLOGIC_TUMOR_STAGE" = "stg",
  "STAGE" = "stg",
  "AGE" = "age",
  "SEX" = "sex",
  "RACE" = "race",
  "RADIATION_THERAPY" = "rad_therapy",
  "OS_STATUS" = "os_status",
  "OS_MONTHS" = "os_months",
  "CHEMOTHERAPY" = "chem_therapy",
  "SMOKING_STATUS" = "smoking_status",
  "SMOKING_PACK_YEARS" = "smoking_pack_years",
  "PACK_YEARS_SMOKED" = "smoking_pack_years"
)


OncoSG_data <- parse_cbioportal_data(
  "OncoSG lung adenocarcinoma data/data_clinical_sample.txt",
  "OncoSG lung adenocarcinoma data/data_clinical_patient.txt",
  "OncoSG lung adenocarcinoma data/alterations_across_samples.tsv",
  "OncoSG lung adenocarcinoma data/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM)-2.txt",
  "OncoSG lung adenocarcinoma data/cna.txt",
  custom_columns
)

CPTAC_data <- parse_cbioportal_data(
  "CPTAC lung adenocarcinoma data/data_clinical_sample.txt",
  "CPTAC lung adenocarcinoma data/data_clinical_patient.txt",
  "CPTAC lung adenocarcinoma data/alterations_across_samples.tsv",
  "CPTAC lung adenocarcinoma data/Z-scores of mRNA expression (RPKM, log2 transformed).txt",
  "CPTAC lung adenocarcinoma data/cna.txt",
  custom_columns
)


TCGA_data_extended <- parse_cbioportal_data(
  "TCGA pancancer lung adenocarcinoma data/data_bcr_clinical_data_sample.txt",
  "TCGA pancancer lung adenocarcinoma data/data_bcr_clinical_data_patient.txt",
  "TCGA pancancer lung adenocarcinoma data/alterations_across_samples.tsv",
  "TCGA pancancer lung adenocarcinoma data/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt",
  "TCGA pancancer lung adenocarcinoma data/cna.txt",
  custom_columns
)


# Add hypoxia and cna data specifically for TCGA
hyp_score_filenames <-
  c(
    "Winter_Hypoxia_Score.txt",
    "Ragnum_Hypoxia_Score.txt",
    "Buffa_Hypoxia_Score.txt"
  )
rel_hyp_score_filenames <-
  paste("TCGA pancancer lung adenocarcinoma data",
        hyp_score_filenames,
        sep = "/")

TCGA_hypoxia_scores <- lapply(rel_hyp_score_filenames, fread) %>%
  reduce(full_join) %>%
  clean_names() %>%
  as_tibble()

TCGA_cna_log2_scores <-
  fread("TCGA pancancer lung adenocarcinoma data/Log2 copy-number values.txt") %>%
  clean_names() %>%
  as_tibble() %>%
  pivot_longer(cols = !c(study_id, sample_id)) %>%
  mutate(name = paste(toupper(name), "CNA_log2", sep = "_")) %>%
  pivot_wider(names_from = name, values_from = value)

TCGA_data <- parse_cbioportal_data(
  "TCGA pancancer lung adenocarcinoma data/data_clinical_sample.txt",
  "TCGA pancancer lung adenocarcinoma data/data_clinical_patient.txt",
  "TCGA pancancer lung adenocarcinoma data/alterations_across_samples.tsv",
  "TCGA pancancer lung adenocarcinoma data/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt",
  "TCGA pancancer lung adenocarcinoma data/cna.txt",
  custom_columns
) %>%
  left_join(y = select(
    TCGA_data_extended,
    sample_id,
    smoking_history,
    smoking_pack_years
  )) %>%
  left_join(TCGA_hypoxia_scores) %>%
  left_join(select(TCGA_cna_log2_scores, sample_id, AURKA_CNA_log2))

merged_data <- bind_rows(OncoSG_data, CPTAC_data) %>%
  bind_rows(TCGA_data) %>%
  mutate(
    stg = str_remove(stg, "STAGE "),
    stg = str_remove(stg, "[A|B]"),
    stg = case_when(
      stg == "I" ~ "1",
      stg == "II" ~ "2",
      stg == "III" ~ "3",
      stg == "IV" ~ "4",
      str_detect(stg, "1") ~ "1",
      str_detect(stg, "2") ~ "2",
      str_detect(stg, "3") ~ "3"
    ),
    stg = as.numeric(stg),
    across(
      .cols = everything(),
      .fns = function(x) {
        ifelse(x == "NA", NA, x)
      }
    ),
    smoking_status = case_when(
      smoking_status == "No" ~ "Non-Smoker",
      smoking_status == "Yes" ~ "Smoker",
      smoking_status == "Non-Smoker" ~ "Non-Smoker",
      smoking_status == "Smoker" ~ "Smoker"
    ),
    across(
      c("AURKA", "TP53", "KRAS", "EGFR"),
      .fns = function(x) {
        factor(x, levels = c("WT", "ALT"))
      }
    ),
    across(
      c("rad_therapy", "smoking_status", "stg", "sex", "chem_therapy"),
      .fns = function(x) {
        as.factor(x)
      }
    ),
    across(
      contains("cna"),
      .fns = function(x) {
        as.factor(ifelse(x == "NP", NA, x))
      }
    ),
    study_id = as.factor(study_id)
  ) %>% # ,
  #   across(
  #     c(
  #       "age",
  #       "TP53_rna_exp",
  #       "KRAS_rna_exp",
  #       "EGFR_rna_exp",
  #       "smoking_pack_years",
  #       "tmb"
  #     ),
  #     .fns = function(x) {
  #       factor(
#         ifelse(x > median(x, na.rm = T), "upper_median", "lower_median"),
#         levels = c("lower_median", "upper_median")
#       )
#     }
#   )
# ) %>%
select(
  age,
  stg,
  sex,
  study_id,
  smoking_status,
  smoking_pack_years,
  smoking_history,
  tmb,
  TP53,
  EGFR,
  KRAS,
  AURKA,
  AURKA_rna_exp,
  TP53_rna_exp,
  KRAS_rna_exp,
  EGFR_rna_exp,
  AURKA_cna,
  AURKA_CNA_log2,
  winter_hypoxia_score,
  buffa_hypoxia_score,
  ragnum_hypoxia_score
) %>%
  filter(!is.na(AURKA_rna_exp)) %>%
  mutate(
    smoking_status = case_when(
      !is.na(smoking_status) ~ smoking_status,
      is.na(smoking_status) & smoking_pack_years > 0 ~ "Smoker",
      is.na(smoking_status) & smoking_history == "1" ~ "Non-Smoker",
      is.na(smoking_status) &
        smoking_history %in% c("2", "3", "4", "5") ~ "Smoker"
    ),
    smoking_status = as.factor(smoking_status),
    AURKA_CNA_log2 = as.numeric(as.character(AURKA_CNA_log2)),
  ) %>%
  select(-smoking_history, -smoking_pack_years)

write.csv(merged_data,
          "Merged annotated data AURKA, KRAS, TP53, EGFR.csv",
          row.names = F)

merged_data <- fread("Merged annotated data AURKA, KRAS, TP53, EGFR.csv") %>% 
  as_tibble()


# corr function ----------------------------------------------------------------
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# Correlations -----------------------------------------------------------------
num_cor_data <- merged_data %>%
  # filter(!if_any(.cols = everything(), .fns = is.na)) %>%
  select(where(is.numeric), -AURKA_rna_exp)

p.mat <- cor.mtest(num_cor_data)
# Scatter plots of numeric variables
pairs(num_cor_data)
corrplot(
  cor(num_cor_data),
  diag = F,
  p.mat = p.mat,
  sig.level = 0.005
)

cat_cor_data <-
  merged_data %>%
  # filter(!if_any(.cols = everything(), .fns = is.na)) %>%
  select(where(is.factor), -AURKA_rna_exp, -stg) %>%
  select(!matches("cna|study"))

f_test_odds_c <- c()
f_test_lower_c <- c()
f_test_upper_c <- c()
f_test_p_value_c <- c()
i_c <- c()
j_c <- c()
both_check <- c()
for (i in names(cat_cor_data)) {
  for (j in names(cat_cor_data)) {
    if (i != j) {
      both <- paste0(sort(c(i, j)), collapse = "")
      if (!both %in% both_check) {
        both_check <- c(both_check, both)
        f_test <-
          fisher.test(table(select(all_of(
            cat_cor_data
          ), i, j)))
        f_test_odds_c <- c(f_test_odds_c, f_test$estimate)
        f_test_lower_c <- c(f_test_lower_c, f_test$conf.int[1])
        f_test_upper_c <- c(f_test_upper_c, f_test$conf.int[2])
        f_test_p_value_c <- c(f_test_p_value_c, f_test$p.value)
        i_c <- c(i_c, i)
        j_c <- c(j_c, j)
      }
    }
  }
}

cat_cor_data_frame <- tibble(
  first_cat = i_c,
  second_cat = j_c,
  odds_ratio = f_test_odds_c,
  upper_conf = f_test_lower_c,
  lower_conf = f_test_upper_c,
  p_value = f_test_p_value_c
) %>%
  mutate(
    odds_ratio = log2(odds_ratio),
    upper_conf = log2(upper_conf),
    lower_conf = log2(lower_conf)
  ) %>%
  mutate(
    both_cat = paste(first_cat, second_cat, sep = "-"),
    both_cat = as.factor(both_cat)
  )

cat_cor_data_frame %>%
  mutate(signif = ifelse(sign(upper_conf * lower_conf) == 1, "Signif", "Non-signif"),
         signif = as.factor(signif)) %>%
  ggplot(aes(x = both_cat,
             y = odds_ratio,
             col = signif)) +
  geom_hline(yintercept = 0, col = "gray35") +
  geom_point(size = 3, position = position_dodge(0.4)) +
  geom_errorbar(
    aes(ymin = lower_conf, ymax = upper_conf),
    width = .2,
    position = position_dodge(0.4),
    linewidth = 1
  ) +
  theme_minimal() +
  xlab("") +
  ylab("log2 odds ratio") +
  theme(legend.title = element_blank(), text = element_text(angle = 25))


# Density chart
merged_data %>%
  pivot_longer(cols = matches("exp"),
               names_to = "hugo_symbol",
               values_to = "mrna_exp_zscore") %>%
  mutate(
    hugo_symbol = str_remove(hugo_symbol, "_rna_exp"),
    hugo_symbol = as.factor(hugo_symbol)
  ) %>%
  ggplot(aes(x = mrna_exp_zscore, fill = hugo_symbol)) +
  geom_density(alpha = .6) +
  xlab("mRNA expression z-scored") +
  # coord_cartesian(xlim = c(-3, 3)) +
  theme_minimal()

merged_data %>%
  ggplot(aes(x = AURKA_rna_exp, y = KRAS_rna_exp)) +
  geom_point()

merged_data %>%
  group_by(KRAS, TP53, EGFR) %>%
  summarise(
    AURKA_rna_exp_median = median(AURKA_rna_exp),
    AURKA_rna_exp_mean = mean(AURKA_rna_exp),
    AURKA_rna_exp_sd = sd(AURKA_rna_exp)
  ) %>%
  write.csv("Basic stats AURKA mRNA expression level.csv", row.names = F)

# Plain regression model -------------------------------------------------------
# Mixed

# Dummy Coding Using Regression
merged_data <- merged_data %>%
  mutate(across(
    matches("hypoxia"),
    .fns = function(x)
      as.numeric(x)
  ))

merged_data_simplified <- merged_data %>%
  select(!matches("hypoxia|log2")) %>%
  filter(if_all(
    .cols = everything(),
    .fns = function(x)
      !is.na(x)
  ))

plain_model_simplified <- summary(lm(AURKA_rna_exp ~ .,
                                     data = merged_data_simplified))

plain_model_tcga <- summary(lm(AURKA_rna_exp ~ .,
                               data = select(
                                 filter(merged_data, study_id == "luad_tcga_pan_can_atlas_2018"),
                                 -study_id
                               )))

# Difference Coding Using Regression
k_stage <- nlevels(merged_data$stg)
contrast_matrix <-
  diag(k_stage - 1) - matrix(1, k_stage - 1, k_stage - 1, byrow = TRUE)

merged_data_modified <- merged_data %>%
  mutate(stg_modified = contrast_matrix)



# TCGA SHOULD BE ANALYSED SEPARATELY

plain_model <- summary(lm(AURKA_rna_exp ~ .,
                          data = filter(
                            merged_data,
                            if_all(
                              .cols = everything(),
                              .fns = function(x)
                                ! is.na(x)
                            )
                          )))

plain_model <- summary(lm(AURKA_rna_exp ~ .,
                          data = select(merged_data, !matches("hypoxia"))))

plain_model <- summary(lm(AURKA_rna_exp ~ .,
                          data = filter(
                            merged_data, if_any(
                              matches("hypoxia"),
                              .fns = function(x)
                                is.na(x)
                            )
                          )))

plain_model <- summary(lm(AURKA_rna_exp ~ .,
                          data = filter(
                            merged_data, if_any(
                              .cols = everything(),
                              .fns = function(x)
                                is.na(x)
                            )
                          )))

# KRAS only
plain_model_kras <- summary(lm(AURKA_rna_exp ~ .,
                               data = select(
                                 filter(merged_data, KRAS == "ALT", EGFR == "WT"),
                                 -KRAS,
                                 -EGFR,
                                 -smoking_status
                               )))

# EGFR only
plain_model_egfr <- summary(lm(AURKA_rna_exp ~ .,
                               data = select(
                                 filter(merged_data, KRAS == "WT", EGFR == "ALT"),
                                 -KRAS,
                                 -EGFR,
                                 -smoking_status
                               )))



male_exp <- filter(merged_data, sex == "male")$AURKA_rna_exp
female_exp <- filter(merged_data, sex == "female")$AURKA_rna_exp
t.test(male_exp, female_exp)

# VIF --------------------------------------------------------------------------
vif(plain_model)



# Regression model with interaction between variables --------------------------
# Mixed
mixed_model <-
  summary(
    lm(
      AURKA_rna_exp ~ age * stg * sex * tmb * TP53 * EGFR * KRAS * AURKA * TP53_rna_exp * KRAS_rna_exp * EGFR_rna_exp,
      data = merged_data
    )
  )

mixed_model_coeffs <-
  as_tibble(mixed_model$coefficients, .name_repair = "universal") %>%
  mutate(intercept = attributes(mixed_model$coefficients)$dimnames[[1]]) %>%
  rename(
    p_value = Pr...t..,
    sd = Std..Error,
    estimate = Estimate,
    t_value = t.value
  )

mixed_model_coeffs %>%
  filter(p_value < 0.01)

# KRAS only, simplified
mixed_model_kras <-
  summary(
    lm(
      AURKA_rna_exp ~ age + stg + sex + tmb + TP53 + AURKA + TP53_rna_exp * KRAS_rna_exp + EGFR_rna_exp,
      data = select(
        filter(merged_data, KRAS == "ALT", EGFR == "WT", !is.na(age), !is.na(stg)),
        -KRAS,
        -EGFR
      )
    )
  )


mixed_model_coeffs <-
  as_tibble(mixed_model$coefficients, .name_repair = "universal") %>%
  mutate(intercept = attributes(mixed_model$coefficients)$dimnames[[1]]) %>%
  rename(
    p_value = Pr...t..,
    sd = Std..Error,
    estimate = Estimate,
    t_value = t.value
  )

mixed_model_coeffs %>%
  filter(p_value < 0.01)

# EGFR only
mixed_model_egfr <-
  summary(
    lm(
      AURKA_rna_exp ~ age + stg + sex + tmb + TP53 + AURKA + TP53_rna_exp * KRAS_rna_exp + EGFR_rna_exp,
      data = select(filter(merged_data, KRAS == "WT", EGFR == "ALT"), -KRAS, -EGFR)
    )
  )

mixed_model_coeffs <-
  as_tibble(mixed_model$coefficients, .name_repair = "universal") %>%
  mutate(intercept = attributes(mixed_model$coefficients)$dimnames[[1]]) %>%
  rename(
    p_value = Pr...t..,
    sd = Std..Error,
    estimate = Estimate,
    t_value = t.value
  )

mixed_model_coeffs %>%
  filter(p_value < 0.01)





# Analysis from 04.05.2023 -----------------------------------------------------
sample_data <-
  fread(
    "TCGA pancancer lung adenocarcinoma data/data_clinical_sample.txt",
    skip = 4,
    na.strings = ""
  ) %>%
  as_tibble() %>%
  rename(sample_id = SAMPLE_ID,
         patient_id = PATIENT_ID,
         tmb = TMB_NONSYNONYMOUS) %>%
  select(sample_id, patient_id, tmb)

clinical_data <-
  fread(
    "TCGA pancancer lung adenocarcinoma data/data_clinical_patient.txt",
    skip = 4,
    na.strings = ""
  ) %>%
  as_tibble() %>%
  rename(
    patient_id = PATIENT_ID,
    stage = AJCC_PATHOLOGIC_TUMOR_STAGE,
    age = AGE,
    sex = SEX,
    race = RACE,
    rad_therapy = RADIATION_THERAPY,
    os_status = OS_STATUS,
    os_months = OS_MONTHS
  ) %>%
  select(patient_id, stage, age, race, rad_therapy, os_status, os_months) %>%
  left_join(sample_data, by = "patient_id") %>%
  mutate(stage = str_remove(stage, "STAGE "))

rna_seq_data <-
  fread(
    "TCGA pancancer lung adenocarcinoma data/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt"
  ) %>%
  as_tibble() %>%
  filter(!if_any(
    .cols = everything(),
    .fns = function(x) {
      is.na(x)
    }
  )) %>%
  rename(
    sample_id = SAMPLE_ID,
    study_id = STUDY_ID,
    TP53_rna_exp = TP53,
    EGFR_rna_exp = EGFR,
    KRAS_rna_exp = KRAS,
    AURKA_rna_exp = AURKA
  ) %>%
  select(-study_id)

protein_exp_data <-
  fread("TCGA pancancer lung adenocarcinoma data/Protein expression z-scores (RPPA).txt") %>%
  as_tibble() %>%
  filter(!if_all(
    .cols = c("TP53", "EGFR", "KRAS", "AURKA"),
    .fns = function(x) {
      is.na(x)
    }
  )) %>%
  rename(
    sample_id = SAMPLE_ID,
    study_id = STUDY_ID,
    TP53_protein_exp = TP53,
    EGFR_protein_exp = EGFR,
    KRAS_protein_exp = KRAS,
    AURKA_protein_exp = AURKA
  ) %>%
  select(-study_id)


alterations_data <-
  fread("TCGA pancancer lung adenocarcinoma data/alterations_across_samples.tsv") %>%
  as_tibble() %>%
  mutate(sample_id = `Sample ID`, study_id = `Study ID`) %>%
  select(sample_id, study_id, TP53, EGFR, KRAS, AURKA) %>%
  mutate(across(
    .cols = matches("TP53|EGFR|KRAS|AURKA"),
    .fns = function(x) {
      ifelse(x == "no alteration", "WT", "ALT")
    }
  ))

annotated_data <- alterations_data %>%
  left_join(rna_seq_data, by = "sample_id") %>%
  left_join(protein_exp_data, by = "sample_id") %>%
  filter(!if_all(
    .cols = matches("exp"),
    .fns = function(x) {
      is.na(x)
    }
  ))

# mRNA expression RSEM ---------------------------------------------------------
# Across all genes
# Z-score density charts
rna_seq_data %>%
  pivot_longer(cols = matches("exp"),
               names_to = "hugo_symbol",
               values_to = "mrna_exp_zscore") %>%
  mutate(hugo_symbol = as.factor(hugo_symbol)) %>%
  ggplot(aes(x = mrna_exp_zscore, fill = hugo_symbol)) +
  geom_density(alpha = .7) +
  # coord_cartesian(xlim = c(-3, 3)) +
  theme_minimal()

# Z-score line chart
rna_seq_data %>%
  mutate(y_axis = 1:nrow(.)) %>%
  pivot_longer(cols = matches("exp"),
               names_to = "hugo_symbol",
               values_to = "mrna_exp_zscore") %>%
  group_by(hugo_symbol) %>%
  reframe(mrna_exp_zscore = sort(mrna_exp_zscore),
          y_axis = y_axis) %>%
  ggplot(aes(y = y_axis, x = mrna_exp_zscore, col = hugo_symbol)) +
  geom_line() +
  theme_minimal()

# Z-score violin plots
annotated_data %>%
  select(contains("rna")) %>%
  pivot_longer(cols = everything(),
               names_to = "hugo_symbol",
               values_to = "mrna_exp_zscore") %>%
  mutate(
    hugo_symbol = str_remove_all(hugo_symbol, "_rna_exp"),
    hugo_symbol = as.factor(hugo_symbol)
  ) %>%
  ggplot(aes(y = mrna_exp_zscore, x = hugo_symbol, fill = hugo_symbol)) +
  geom_violin() +
  geom_boxplot(width = .2) +
  theme_minimal()

annotated_data %>%
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


# Protein expression RPPA ------------------------------------------------------
# Z-score density charts
protein_exp_data %>%
  pivot_longer(cols = matches("exp"),
               names_to = "hugo_symbol",
               values_to = "prot_exp_zscore") %>%
  mutate(hugo_symbol = as.factor(hugo_symbol)) %>%
  ggplot(aes(x = prot_exp_zscore, fill = hugo_symbol)) +
  geom_density(alpha = .7) +
  # coord_cartesian(xlim = c(-3, 3)) +
  theme_minimal()

# Z-score line chart
protein_exp_data %>%
  mutate(y_axis = 1:nrow(.)) %>%
  pivot_longer(cols = matches("exp"),
               names_to = "hugo_symbol",
               values_to = "prot_exp_zscore") %>%
  group_by(hugo_symbol) %>%
  reframe(prot_exp_zscore = sort(prot_exp_zscore, na.last = F),
          y_axis = y_axis) %>%
  ggplot(aes(y = y_axis, x = prot_exp_zscore, col = hugo_symbol)) +
  geom_line() +
  theme_minimal()
