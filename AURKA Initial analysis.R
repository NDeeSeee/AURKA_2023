# Attach requirement packages and setting WD -----------------------------------
packages_names <-
  c("tidyverse", "data.table", "readxl", "reshape2", "rstudioapi")

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
  custom_columns
)

CPTAC_data <- parse_cbioportal_data(
  "CPTAC lung adenocarcinoma data/data_clinical_sample.txt",
  "CPTAC lung adenocarcinoma data/data_clinical_patient.txt",
  "CPTAC lung adenocarcinoma data/alterations_across_samples.tsv",
  "CPTAC lung adenocarcinoma data/Z-scores of mRNA expression (RPKM, log2 transformed).txt",
  custom_columns
)


TCGA_data_extended <- parse_cbioportal_data(
  "TCGA pancancer lung adenocarcinoma data/data_bcr_clinical_data_sample.txt",
  "TCGA pancancer lung adenocarcinoma data/data_bcr_clinical_data_patient.txt",
  "TCGA pancancer lung adenocarcinoma data/alterations_across_samples.tsv",
  "TCGA pancancer lung adenocarcinoma data/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt",
  custom_columns
)

TCGA_data <- parse_cbioportal_data(
  "TCGA pancancer lung adenocarcinoma data/data_clinical_sample.txt",
  "TCGA pancancer lung adenocarcinoma data/data_clinical_patient.txt",
  "TCGA pancancer lung adenocarcinoma data/alterations_across_samples.tsv",
  "TCGA pancancer lung adenocarcinoma data/mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt",
  custom_columns
) %>%
  left_join(y = select(
    TCGA_data_extended,
    sample_id,
    smoking_history,
    smoking_pack_years
  ))


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
      smoking_status == "No" ~ "0",
      smoking_status == "Yes" ~ "1",
      smoking_status == "Non-Smoker" ~ "0",
      smoking_status == "Smoker" ~ "1"
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
      c(
        "age",
        "TP53_rna_exp",
        "KRAS_rna_exp",
        "EGFR_rna_exp",
        "smoking_pack_years",
        "tmb"
      ),
      .fns = function(x) {
        factor(
          ifelse(x > median(x, na.rm = T), "upper_median", "lower_median"),
          levels = c("lower_median", "upper_median")
        )
      }
    )
  ) %>%
  select(
    age,
    sex,
    smoking_status,
    smoking_pack_years,
    chem_therapy,
    rad_therapy,
    tmb,
    TP53,
    EGFR,
    KRAS,
    AURKA,
    AURKA_rna_exp,
    TP53_rna_exp,
    KRAS_rna_exp,
    EGFR_rna_exp
  ) %>%
  filter(!is.na(AURKA_rna_exp))

write.csv(merged_data,
          "Merged annotated data AURKA, KRAS, TP53, EGFR.csv",
          row.names = F)

lm(AURKA_rna_exp ~ c(sex, age), data = merged_data)

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
