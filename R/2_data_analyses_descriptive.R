library(tidyverse)
source("global/fun.R")
source("global/obj.R")

df_epi_coded_eng <- read.csv("inputs/df_epi_cap_cleaned.csv") %>% 
  dplyr::select(
    -hos_admission_date,
    -hos_admission_time,
    -hos_abx_start_date,
    -hos_abx_start_time,
    -contains("radiology_")
  ) %>% 
  dplyr::left_join(
    read.csv("inputs/df_lab1_cap_cleaned.csv") %>% 
      dplyr::transmute(
        id = id,
        
        # redefine workLab results to systemic & carriage 
        # workLab_urine = convert_bool(workLab_urine_binax_now),
        workLab_systemic = case_when(
          workLab_urine_binax_now == 1L |
            (workLab_blood_bactec_result == "1" & 
               workLab_culture_result_blood_maldi_tof_confirmed_1 == "Streptococcus pneumoniae") ~ 1L,
          TRUE ~ 0L
          ),
        workLab_carriage = case_when(
          workLab_np_result == "1" |
            workLab_sputum == "1" ~ 1L,
          TRUE ~ 0L
          ),
        
        workLab_overall = case_when(
          workLab_systemic == "1" |
            workLab_carriage == "1" ~ 1L,
          TRUE ~ 0L
          )
      )
    ,
    by = "id"
  ) %>%
  # mutate some numeric columns to categorical
  dplyr::mutate(
    hos_com_curb65 = as.character(hos_com_curb65)
  ) %>% 
  glimpse()

num_val <- c(
  "age",
  "hos_symptom_onset",
  "hos_bmi",
  "epi_smoking_day_n",
  "hos_vacc_covid",
  
  # "hos_com_psi_class", # value in range of I to IV
  # "hos_com_curb65", # value in range of 2 to 5
  "hos_com_cci_score",
  "hos_los"
)


# ptest for categorical data ###################################################
ptest_matrix_all <- sumstats::generate_or_chisq_report(df_input = df_epi_coded_eng %>% 
                                                         dplyr::select(-id,
                                                                       -all_of(num_val),
                                                                       -contains("hos_abx_")),
                                                       binary_disease = "workLab_overall")

ptest_matrix_table_report <- sumstats::chisq_fisher_to_2D_table(result = ptest_matrix_all) %>%
  dplyr::select(-significance) %>% 
  dplyr::mutate(
    final_p = case_when(
      (preferred_usage == "Chi-squared" &
         chisq_p < 0.01) ~ "< 0.01",
      (preferred_usage == "Chi-squared" &
         (chisq_p >= 0.01 & chisq_p < 0.05)) ~ "< 0.05",
      (preferred_usage == "Chi-squared" &
         (chisq_p >= 0.05)) ~ as.character(round(chisq_p, 2)),
      
      (preferred_usage == "Fisher's Exact" &
         fisher_p < 0.01) ~ "< 0.01",
      (preferred_usage == "Fisher's Exact" &
         (fisher_p >= 0.01 & fisher_p < 0.05)) ~ "< 0.05",
      (preferred_usage == "Fisher's Exact" &
         (fisher_p >= 0.05)) ~ as.character(round(fisher_p, 2))
    ),
    significance = ifelse(final_p < 0.05, "occur", "no")
  ) %>% 
  glimpse()

# compile all assessment
df_assessed <- df_epi_coded_eng %>% 
  dplyr::select(-id,
                -all_of(num_val)
                ) %>% 
  dplyr::mutate(across(everything(), as.character)) %>% 
  tidyr::pivot_longer(cols = everything(),
                      names_to = "variable", values_to = "value") %>% 
  dplyr::group_by(variable, value) %>% 
  dplyr::summarise(count_all = n(), .groups = "drop") %>% 
  dplyr::mutate(
    percent = round(count_all/nrow(df_epi_coded_eng)*100, 1),
    report_assessed_all = paste0(count_all, " (", percent, "%)")
  ) %>% 
  dplyr::select(-count_all,
                -percent) %>% 
  # view() %>%
  glimpse()

# compiled pneumo positive
df_compiled_positivePneumo <- dplyr::left_join(
  df_epi_coded_eng %>% 
    dplyr::select(-id,
                  -all_of(num_val)
    ) %>% 
    dplyr::filter(workLab_overall == 1) %>% 
    dplyr::mutate(across(everything(), as.character)) %>%
    tidyr::pivot_longer(cols = everything(),
                        names_to = "variable", values_to = "value") %>%
    dplyr::group_by(variable) %>% # variable only to make this total of positives
    dplyr::summarise(count_all = n(), .groups = "drop")
  ,
  df_epi_coded_eng %>% 
    dplyr::select(-id,
                  -all_of(num_val)
    ) %>% 
    dplyr::mutate(across(-workLab_overall, as.character)) %>%
    tidyr::pivot_longer(cols = -workLab_overall, 
                        names_to = "variable", 
                        values_to = "value"
                        ) %>% 
    dplyr::group_by(variable, value, workLab_overall) %>% 
    dplyr::summarise(count_positivePneumo = n(), .groups = "drop") %>% 
    dplyr::filter(workLab_overall == 1) %>% 
    dplyr::select(-workLab_overall)
  ,
  by = c("variable")
) %>% 
  dplyr::mutate(
    count_positivePneumo = case_when(
      is.na(count_positivePneumo) ~ 0,
      TRUE ~ count_positivePneumo
    ),
    percent = round(count_positivePneumo/count_all*100, 1),
    report_positive_all = paste0(count_positivePneumo, " (", percent, "%)")
  ) %>% 
  dplyr::select(-count_all,
                -count_positivePneumo,
                -percent) %>% 
  # view() %>% 
  glimpse()


# ptest for numerical data #####################################################
ptest_cont_all <- df_epi_coded_eng %>% 
  dplyr::select(all_of(num_val),
                workLab_overall
                ) %>% 
  tidyr::pivot_longer(
    cols = all_of(num_val),
    names_to = "variable",
    values_to = "value"
  ) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    wilcox_p = tryCatch(
      wilcox.test(value ~ workLab_overall, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    final_p = case_when(
      wilcox_p < 0.01 ~ "< 0.01",
      wilcox_p >= 0.01 & wilcox_p < 0.05 ~ "< 0.05",
      wilcox_p >= 0.05 ~ as.character(round(wilcox_p, 2))
    ),
    significance = ifelse(final_p < 0.05, "occur", "no")
  ) %>% 
  glimpse()

# compile all assessment
df_assessed_num <- df_epi_coded_eng %>% 
  dplyr::select(all_of(num_val)) %>% 
  dplyr::summarise(across(everything(), list(
    "Median (IQR)" = ~sprintf("%.1f (%.1f - %.1f)",
                              median(.x, na.rm = TRUE),
                              quantile(.x, 0.25, na.rm = TRUE),
                              quantile(.x, 0.75, na.rm = TRUE)),
    Unknown = ~paste0(sum(is.na(.x)), " (",
                      round(mean(is.na(.x)) * 100, 1), "%)")
  ))
  ) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = c("variable", "value"),
    names_sep = "_(?=[^_]+$)",
    values_to = "report_assessed_all"
  ) %>%
  # view() %>%
  glimpse()

# compiled pneumo positive
df_compiled_positivePneumo_num <- df_epi_coded_eng %>% 
  dplyr::filter(workLab_overall == "1") %>% 
  dplyr::select(all_of(num_val)) %>% 
  dplyr::summarise(across(everything(), list(
    "Median (IQR)" = ~sprintf("%.1f (%.1f - %.1f)",
                              median(.x, na.rm = TRUE),
                              quantile(.x, 0.25, na.rm = TRUE),
                              quantile(.x, 0.75, na.rm = TRUE)),
    Unknown = ~paste0(sum(is.na(.x)), " (",
                      round(mean(is.na(.x)) * 100, 1), "%)")
  ))
  ) %>%
  tidyr::pivot_longer(
    everything(),
    names_to = c("variable", "value"),
    names_sep = "_(?=[^_]+$)",
    values_to = "report_positive_all"
  ) %>%
  # view() %>%
  glimpse()


# combine report ###############################################################
compile_all_report_with_pValues <- dplyr::bind_rows(
  dplyr::left_join(
    df_assessed,
    df_compiled_positivePneumo,
    by = c("variable", "value"),
    relationship = "many-to-many"
  ) %>% 
    dplyr::left_join(
      ptest_matrix_table_report,
      by = c("variable")
    )
  ,
  dplyr::left_join(
    df_assessed_num,
    df_compiled_positivePneumo_num,
    by = c("variable", "value"),
    relationship = "many-to-many"
  ) %>% 
    dplyr::left_join(
      ptest_cont_all,
      by = c("variable")
    )
) %>% 
  
  dplyr::distinct(variable, value, .keep_all = T) %>% 
  dplyr::mutate(
    value = case_when(
      is.na(value) ~ "Unknown",
      value == "0" ~ "No",
      value == "1" ~ "Yes",
      TRUE ~ value
    ),
    varClass = varClass_map[variable],
    arrangement = as.numeric(recode(variable, !!!var_arrangement)),
    
    variable = recode(variable, !!!var_map),
    
    # correct value for smoking
    value = case_when(
      variable == "Smoking intensity (cig/day)" &
        value == "Unknown" ~ "Not smoking",
      TRUE ~ value
    ),
    preferred_usage = ifelse(is.na(preferred_usage), "Wilcoxon Rank Sum",
                             preferred_usage)
  ) %>%
  dplyr::relocate(varClass) %>% 
  arrange(arrangement) %>%
  dplyr::filter(!(value == "Unknown" &
                  report_assessed_all == "0 (0%)" &
                  report_positive_all == "0 (0%)")
                ) %>%
  # omit varclass Antibiotics
  dplyr::filter(!varClass %in% c("Antibiotics", "Antifungal", "Antiviral")
                ) %>% 
  # view() %>%
  glimpse()

write.csv(compile_all_report_with_pValues %>% 
            dplyr::select(-arrangement),
          "outputs/out_1allpositives_1all_tables.csv",
          row.names = F)


# Additional data for workLab results ##########################################
workLab_list <- c(
  "workLab_overall",
  "workLab_carriage",
  "workLab_systemic"
)

# moron R dimension requirement
compile_list <- list()

# pre-compute overall positive counts BEFORE the loop
df_overall_positive_carriage <- df_epi_coded_eng %>%
  dplyr::filter(.data[["workLab_overall"]] == 1) %>%
  dplyr::select(-id, -all_of(num_val)) %>%
  dplyr::mutate(across(everything(), as.character)) %>%
  tidyr::pivot_longer(-"workLab_overall",
                      names_to = "variable",
                      values_to = "value") %>%
  dplyr::count(variable, value, name = "count_overall_pos") %>%  # denominator
  glimpse()

df_overall_positive_num_n <- df_epi_coded_eng %>%
  dplyr::filter(.data[["workLab_overall"]] == 1) %>%
  dplyr::select(all_of(num_val)) %>%
  nrow()  # numeric denominator

for(i in workLab_list){
  
  # categorical data
  ptest_matrix_all <- sumstats::generate_or_chisq_report(
    df_input = df_epi_coded_eng %>% 
      dplyr::select(-id,
                    -all_of(num_val),
                    contains("hos_abx_")),
    binary_disease = i
  )
  
  ptest_matrix_table_report <- sumstats::chisq_fisher_to_2D_table(
    result = ptest_matrix_all
  ) %>%
    dplyr::mutate(
      raw_p = dplyr::case_when(
        preferred_usage == "Chi-squared" ~ chisq_p,
        preferred_usage == "Fisher's Exact" ~ fisher_p
      ),
      final_p = case_when(
        raw_p < 0.01 ~ "< 0.01",
        raw_p < 0.05 ~ "< 0.05",
        TRUE ~ as.character(round(raw_p, 2))
      ),
      significance = ifelse(raw_p < 0.05, "occur", "no")
    ) %>%
    dplyr::select(variable, final_p, significance) %>%
    dplyr::rename(
      !!paste0("p_", i) := final_p,
      !!paste0("sig_", i) := significance
    )
  
  # assessed
  df_assessed <- df_epi_coded_eng %>% 
    dplyr::select(-id, -all_of(num_val)) %>% 
    dplyr::mutate(across(everything(), as.character)) %>% 
    tidyr::pivot_longer(everything(),
                        names_to = "variable",
                        values_to = "value") %>% 
    dplyr::count(variable, value, name = "count_all") %>% 
    dplyr::mutate(
      percent = round(count_all / nrow(df_epi_coded_eng) * 100, 1),
      report_assessed_all = paste0(count_all, " (", percent, "%)")
    ) %>% 
    dplyr::select(variable, value, count_all, report_assessed_all)  # Keep count_all for joining
  
  # determine denominator based on which dataset
  # workLab_overall = count per variable/overall positives
  # workLab_carriage/systemic= count per variable/workLab_overall positives
  
  if(i == "workLab_overall") {
    denom_carriage <- df_assessed %>% 
      dplyr::select(variable, value, count_all) %>%
      dplyr::rename(denom = count_all)
  } else {
    # workLab_overall positive counts as denominator
    denom_carriage <- df_overall_positive_carriage %>%
      dplyr::rename(denom = count_overall_pos)
  }
  
  # positives
  df_positive_carriage <- df_epi_coded_eng %>%
    dplyr::filter(.data[[i]] == 1) %>%
    dplyr::select(-id, -all_of(num_val)) %>%
    dplyr::mutate(across(everything(), as.character)) %>%
    tidyr::pivot_longer(-all_of(i),
                        names_to = "variable",
                        values_to = "value") %>%
    dplyr::count(variable, value, name = "count_pos") %>%
    dplyr::right_join(df_assessed,
                      by = c("variable", "value")
    ) %>%
    dplyr::left_join(denom_carriage,
                     by = c("variable", "value")
    ) %>%
    dplyr::mutate(
      count_pos = ifelse(is.na(count_pos), 0, count_pos),
      denom = ifelse(is.na(denom), 1, denom),  # Avoid division by zero
      percent = round(count_pos / denom * 100, 1),  # FIXED: divide by correct denom
      !!paste0("report_positive_", i) :=
        paste0(count_pos, " (", percent, "%)")
    ) %>%
    dplyr::select(variable, value, starts_with("report_positive_"))
  
  # p-values
  ptest_cont_all <- df_epi_coded_eng %>% 
    dplyr::select(all_of(num_val), all_of(i)) %>% 
    tidyr::pivot_longer(all_of(num_val),
                        names_to = "variable",
                        values_to = "value") %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      raw_p = tryCatch(
        wilcox.test(value ~ .data[[i]], exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      !!paste0("p_", i) := case_when(
        raw_p < 0.01 ~ "< 0.01",
        raw_p < 0.05 ~ "< 0.05",
        TRUE ~ as.character(round(raw_p, 2))
      )
    ) %>%
    dplyr::select(-raw_p)
  
  # numeric data
  df_assessed_num <- df_epi_coded_eng %>% 
    dplyr::select(all_of(num_val)) %>% 
    dplyr::summarise(across(everything(), list(
      "Median (IQR)" = ~sprintf("%.1f (%.1f - %.1f)",
                                median(.x, na.rm = TRUE),
                                quantile(.x, 0.25, na.rm = TRUE),
                                quantile(.x, 0.75, na.rm = TRUE)),
      Unknown = ~paste0(sum(is.na(.x)), " (",
                        round(mean(is.na(.x)) * 100, 1), "%)")
    ))) %>%
    tidyr::pivot_longer(
      everything(),
      names_to = c("variable", "value"),
      names_sep = "_(?=[^_]+$)",
      values_to = "report_assessed_all"
    )
  
  # numeric denominator
  if(i == "workLab_overall") {
    denom_num <- df_epi_coded_eng %>%
      dplyr::filter(.data[["workLab_overall"]] == 1) %>%
      nrow()
  } else {
    denom_num <- df_overall_positive_num_n  # workLab_overall n as denominator
  }
  
  # positives
  df_positive_num <- df_epi_coded_eng %>% 
    dplyr::filter(.data[[i]] == 1) %>% 
    dplyr::select(all_of(num_val)) %>% 
    dplyr::summarise(across(everything(), list(
      "Median (IQR)" = ~sprintf("%.1f (%.1f - %.1f)",
                                median(.x, na.rm = TRUE),
                                quantile(.x, 0.25, na.rm = TRUE),
                                quantile(.x, 0.75, na.rm = TRUE)),
      Unknown = ~paste0(
        sum(is.na(.x)), " (",
        round(sum(is.na(.x)) / denom_num * 100, 1),  # FIXED: divide by correct denom
        "%)"
      )
    ))) %>%
    tidyr::pivot_longer(
      everything(),
      names_to = c("variable", "value"),
      names_sep = "_(?=[^_]+$)",
      values_to = paste0("report_positive_", i)
    )
  
  # combine
  compile_i <- bind_rows(
    dplyr::left_join(df_assessed %>% select(-count_all),  # Drop helper column
                     df_positive_carriage,
                     by = c("variable", "value")) %>%
      dplyr::left_join(ptest_matrix_table_report, by = "variable")
    ,
    dplyr::left_join(df_assessed_num, df_positive_num,
                     by = c("variable", "value")) %>%
      dplyr::left_join(ptest_cont_all, by = "variable")
  )
  
  compile_list[[i]] <- compile_i
}

additional_workLab <- Reduce(function(x, y){
  dplyr::left_join(x, y,
                   by = c("variable", "value", "report_assessed_all"))
}, compile_list) %>%
  dplyr::mutate(
    value = case_when(
      is.na(value) ~ "Unknown",
      value == "0" ~ "No",
      value == "1" ~ "Yes",
      TRUE ~ value
    ),
    varClass = varClass_map[variable],
    arrangement = as.numeric(recode(variable, !!!var_arrangement)),
    
    variable = recode(variable, !!!var_map),
    
    # correct value for smoking
    value = case_when(
      variable == "Smoking intensity (cig/day)" &
        value == "Unknown" ~ "Not smoking",
      TRUE ~ value
    )
  ) %>%
  dplyr::relocate(varClass) %>% 
  arrange(arrangement) %>%
  dplyr::filter(!varClass %in% c("Antibiotics", "Antifungal", "Antiviral")
  ) %>% 
  # view() %>%
  glimpse()

write.csv(additional_workLab %>% 
            dplyr::select(-arrangement),
          "outputs/out_additional_each_workLab_numeric.csv",
          row.names = F)



# Additional data for workLab results (kappa agreement) ########################
test_workLab <- read.csv("inputs/df_epi_cap_cleaned.csv") %>% 
  dplyr::select(
    -hos_admission_date,
    -hos_admission_time,
    -hos_abx_start_date,
    -hos_abx_start_time,
    -contains("radiology_")
  ) %>% 
  dplyr::left_join(
    read.csv("inputs/df_lab1_cap_cleaned.csv")
    ,
    by = "id"
  ) %>% 
  # dplyr::select(contains("workLab")) %>% 
  dplyr::transmute(
    NP = workLab_np_result,
    sputum = as.numeric(workLab_sputum),
    blood = case_when(
      (workLab_blood_bactec_result == "1" & 
         workLab_culture_result_blood_maldi_tof_confirmed_1 == "Streptococcus pneumoniae") ~ 1L,
      is.na(workLab_blood_bactec_result) ~ NA,
      TRUE ~ 0L
    ),
    urine = workLab_urine_binax_now
  ) %>% 
  # slice(1:100) %>% # test model
  glimpse()

test_workLab[, -1] <- lapply(test_workLab[, -1], as.numeric)


library(psych)
ppa_npa <- function(tab) {
  a <- tab[2,2]
  b <- tab[2,1]
  c <- tab[1,2]
  d <- tab[1,1]
  
  ppa <- if ((a + c) == 0) NA else a / (a + c)
  npa <- if ((b + d) == 0) NA else d / (b + d)
  
  return(c(PPA = ppa, NPA = npa))
}

analyze_pair <- function(data, var1, var2, ref = NULL){
  
  # Remove NA pairwise
  sub <- data[!is.na(data[[var1]]) & !is.na(data[[var2]]), ]
  
  # 2×2 table
  tab <- table(sub[[var1]], sub[[var2]])
  
  # Ensure full 2×2
  tab <- as.matrix(tab)
  tab <- tab[match(c("0","1"), rownames(tab)), 
             match(c("0","1"), colnames(tab))]
  tab[is.na(tab)] <- 0
  
  # Rename
  rownames(tab) <- paste(var1, c("-", "+"))
  colnames(tab) <- paste(var2, c("-", "+"))
  
  # Kappa
  kappa <- psych::cohen.kappa(tab)$kappa
  
  # Agreement
  agreement <- sum(diag(tab)) / sum(tab)
  
  # Operating characteristics (if reference provided)
  sens <- spec <- ppv <- npv <- NA
  
  if (!is.null(ref) && ref == var1) {
    a <- tab[2,2]; b <- tab[2,1]; c <- tab[1,2]; d <- tab[1,1]
    sens <- a / (a + c)
    spec <- d / (b + d)
    ppv  <- a / (a + b)
    npv  <- d / (c + d)
  }
  
  # PPA-NPA calculation
  ppa <- ppa_npa(tab)["PPA"]
  npa <- ppa_npa(tab)["NPA"]
  
  
  return(list(
    table = tab,
    kappa = kappa,
    agreement = agreement,
    sensitivity = sens,
    specificity = spec,
    ppv = ppv,
    npv = npv,
    ppa = ppa,
    npa = npa,
    n = sum(tab)
  ))
}



pairs <- list(
  c("NP", "sputum"),
  c("NP", "blood"),
  c("NP", "urine"),
  c("sputum", "blood"),
  c("sputum", "urine"),
  c("blood", "urine")
)

results <- lapply(pairs, function(p) {
  analyze_pair(test_workLab, p[1], p[2], ref = "NP")  # set NP as reference if appropriate
})

names(results) <- sapply(pairs, paste, collapse = "_vs_")


summary_test_workLab <- data.frame(
  "Sample Comparison" = names(results),
  N = sapply(results, function(x) x$n),
  Agreement = sapply(results, function(x) x$agreement),
  Kappa = sapply(results, function(x) x$kappa),
  Sensitivity = sapply(results, function(x) x$sensitivity),
  Specificity = sapply(results, function(x) x$specificity),
  PPA = sapply(results, function(x) x$ppa),
  NPA = sapply(results, function(x) x$npa)
)

summary_test_workLab

write.csv(summary_test_workLab,
          "outputs/out_additional_each_workLab_sample_agreement.csv",
          row.names = F)

