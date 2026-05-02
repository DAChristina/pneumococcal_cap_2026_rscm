library(tidyverse)
library(logistf)
source("global/fun.R")
source("global/obj.R")

df_epi_coded_eng_pos <- read.csv("inputs/df_epi_cap_cleaned.csv") %>% 
  dplyr::select(
    -age,
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
  dplyr::mutate(
    # change numerical data to categorical data
    hos_symptom_onset = as.factor(
      ifelse(hos_symptom_onset <= 2, "Early (<3days)", "Late (>=3days)")
      ),
    hos_bmi = case_when(
      hos_bmi < 18.5 ~ "Underweight",
      hos_bmi >= 18.5 & hos_bmi < 25 ~ "Normal",
      hos_bmi >= 25 & hos_bmi < 30 ~ "Overweight", 
      hos_bmi >= 30 ~ "Obese",
      TRUE ~ NA_character_
    ),
    hos_bmi = factor(hos_bmi,
                     levels = c("Underweight", "Normal",
                                "Overweight", "Obese")
                     ),
    
    epi_smoking_day_n = as.factor(
      ifelse(epi_smoking_day_n <= 20, "A pack", "More")
      ),
    hos_com_curb65 = factor(ifelse(is.na(hos_com_curb65), "Unknown",
                                   ifelse(hos_com_curb65 >= 3, "High (4-5)",
                                          "Low (1-3)")),
                            levels = c("Low (1-3)", "High (4-5)", "Unknown")
      ),
    
    hos_com_cci_score = case_when(
      hos_com_cci_score <= 2 ~ "Low (0-2)",
      hos_com_cci_score %in% 3:4 ~ "Moderate (3-4)",
      hos_com_cci_score >= 5 ~ "High (≥5)",
      TRUE ~ NA_character_
    ),
    hos_com_cci_score = factor(hos_com_cci_score, 
                          levels = c("Low (0-2)", "Moderate (3-4)",
                                     "High (≥5)")
                          ),
    hos_los = case_when(
      hos_los < 10 ~ "Short Stay (<10 days)",     # Systemic???
      hos_los >= 10 ~ "Prolonged Stay (≥10 days)", # Carriage???
      TRUE ~ NA_character_
    ),
    hos_los = factor(hos_los,
                     levels = c("Short Stay (<10 days)",
                                "Prolonged Stay (≥10 days)")
                     ),
    
    # correct factor type
    ageGroup = factor(ageGroup,
                      levels = c("18-64 (Adult)", "65+ (Elderly)")
                      ),
    sex = factor(sex,
                 levels = c("Female", "Male")
                 ),
    epi_smoking_history = factor(epi_smoking_history,
                                 levels = c("Never smoked",
                                            "Previously smoking",
                                            "Currently smoking",
                                            "Unknown")
                                 ),
    hos_vacc_covid = as.factor(hos_vacc_covid),
    
    across(contains(c("_bool_", "_died", "rad_")), 
           ~ case_when(
             .x == 1 ~ "Yes",
             .x == 0 ~ "No",
             is.na(.x) ~ "Unknown",
             TRUE ~ "Unknown"
           ) %>% factor(levels = c("No", "Yes", "Unknown"))
      ),
    
  ) %>% 
  dplyr::select(
    -contains(c("hos_abx_",
              "hos_treat_bool_",
              "hos_therapy_bool_"))
    # -hos_los # due to severe missing data (> 35% overall, n = 7 in positives)
  ) %>% 
  dplyr::filter(workLab_overall == 1) %>% 
  glimpse()


# Firth logistic: What predicts systemic vs carriage? ##########################
separation_results <- scan_perfect_separation_all(df_epi_coded_eng_pos,
                                              "workLab_systemic") %>% 
  dplyr::filter(perfect_sep == TRUE) %>% 
  # view() %>%
  glimpse()

get_columns <- df_epi_coded_eng_pos %>%
  dplyr::select(
    -all_of(separation_results$variable),
    -id,
    -workLab_systemic,
    -contains(c("workLab_"))
  ) %>%
  colnames()

# univar first for selecting variables
univar_results <- map_dfr(get_columns, ~{
  formula_str <- paste("workLab_systemic ~", .x)
  
  model <- logistf(as.formula(formula_str), 
                   data = df_epi_coded_eng_pos, 
                   firth = TRUE)
  
  coef_table <- data.frame(
    value = names(model$coefficients),
    estimate = model$coefficients,
    p.value = model$prob,
    conf.low = model$ci.lower,
    conf.high = model$ci.upper,
    stringsAsFactors = FALSE
  )
  
  # Filter out intercept and add variable name
  coef_table %>%
    filter(value != "(Intercept)") %>%
    mutate(
      variable = .x,
      OR = exp(estimate),
      OR_conf.low = exp(conf.low),
      OR_conf.high = exp(conf.high),
      OR_report = paste0(round(OR, 2), " (",
                         round(OR_conf.low, 2), "-",
                         round(OR_conf.high, 2), ")")
    ) %>%
    select(variable, value, estimate, OR_report, OR, p.value)
})

print_univar <- univar_results %>%
  dplyr::mutate(
    p.report = case_when(
      p.value < 0.01 ~ "< 0.01",
      p.value < 0.05 ~ "< 0.05",
      TRUE ~ as.character(round(p.value, 2))
    )
  ) %>% 
  arrange(p.value) %>% 
  glimpse()


# run multivar
get_columns2 <- c(print_univar$variable[1:2]
                   # print_univar$variable[5],
                   # print_univar$variable[8]
                  )
# get_columns2
# "hos_symptom_onset" "hos_com_curb65"    "hos_los"
# "hos_com_bool_cad" 
# "hos_com_cci_score"

formula <- as.formula(paste("workLab_systemic", "~",
                            paste(get_columns2, collapse = " + ")))

model_primary <- logistf(
  formula,  # Max 2-3 predictors
  data = df_epi_coded_eng_pos,
  firth = TRUE,
  pl = TRUE  # Profile likelihood CIs
)

# result check
summary(model_primary)

print_multivar <- data.frame(
  value = names(model_primary$coefficients),
  OR = exp(model_primary$coefficients),
  OR_report = paste0(
    round(exp(model_primary$coefficients), 2), " (", 
    round(exp(model_primary$ci.lower), 2), "–", 
    round(exp(model_primary$ci.upper), 2), ")"
  ),
  p.value = model_primary$prob
) %>% 
  dplyr::mutate(
    p.report = case_when(
      p.value < 0.01 ~ "< 0.01",
      p.value < 0.05 ~ "< 0.05",
      TRUE ~ as.character(round(p.value, 2))
      )
  ) %>% 
  glimpse()

combine_report <- dplyr::left_join(
  print_univar %>% 
    dplyr::select(-estimate) %>% 
    dplyr::rename(
      "Univariable OR" = OR,
      "Univariable OR (95% CI)" = OR_report,
      "Univariable p" = p.value,
      "Univariable p report" = p.report
    )
  ,
  print_multivar %>% 
    dplyr::rename(
      "Multivariable OR" = OR,
      "Multivariable OR (95% CI)" = OR_report,
      "Multivariable p" = p.value,
      "Multivariable p report" = p.report
    )
  ,
  by = "value"
) %>% 
  dplyr::mutate(
    value = stringr::str_remove(value, variable),
    
    varClass = varClass_map[variable],
    arrangement = as.numeric(recode(variable, !!!var_arrangement)),
    variable = recode(variable, !!!var_map)
    ) %>% 
  dplyr::relocate(varClass) %>% 
  arrange(arrangement) %>%
  glimpse()

write.csv(combine_report %>% 
            dplyr::select(-arrangement),
          "outputs/out_1allpositives_2OR.csv",
          row.names = F)

# Further testing ##############################################################
# 1. Likelihood ratio test
# Test hos_symptom_onset
model_without_onset <- logistf(
  workLab_systemic ~ hos_died,
  data = df_epi_coded_eng_pos,
  firth = TRUE
)

lr_onset <- 2 * (model_primary$loglik - model_without_onset$loglik)
p_onset <- 1 - pchisq(lr_onset, df = 1)

# Test hos_died
model_without_curb65 <- logistf(
  workLab_systemic ~ hos_symptom_onset,
  data = df_epi_coded_eng_pos,
  firth = TRUE
)

lr_curb65 <- 2 * (model_primary$loglik - model_without_curb65$loglik)
p_curb65 <- 1 - pchisq(lr_curb65, df = 1)

cat("LR test - hos_symptom_onset: p =", round(p_onset, 4), "\n") # adding Onset significantly improve the model
cat("LR test - hos_died: p =", round(p_curb65, 4), "\n")

# 2. Diagnostic plots
# Get model data (cleaner with 2 variables)
model_primary <- df_epi_coded_eng_pos %>%
  filter(complete.cases(select(., workLab_systemic, hos_symptom_onset, hos_died))) %>%
  mutate(
    predicted_prob = predict(model_primary, type = "response"),
    deviance_resid = residuals(model_primary, type = "deviance")
  )

# Key diagnostic plot: Predicted probabilities by actual outcome
diagnostic_plot <- ggplot(model_primary,
                          aes(x = factor(workLab_systemic),
                              y = predicted_prob)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2) +
  labs(title = "Model Diagnostic: Predicted vs Actual",
       x = "Actual Outcome (0=Carriage, 1=Systemic)",
       y = "Predicted Probability",
       caption = paste("n =", nrow(model_primary))) +
  theme_bw()

print(diagnostic_plot)

# Simple separation check
separation_check <- any(model_primary$predicted_prob %in% c(0, 1))
cat("Perfect separation detected:", separation_check, "\n")
cat("Predicted probability range:", 
    round(min(model_primary$predicted_prob), 3), "to", 
    round(max(model_primary$predicted_prob), 3), "\n")

