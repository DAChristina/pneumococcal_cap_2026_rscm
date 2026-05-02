# raw data were given to me at 24/2/2026.
# adjustment needed for merged columns and tables.
# ID inconsistency occurs.

# epiData n = 360; eligible only 354 (355?) 6 ID duplicated = 349
# epiData: Data streptococcus_Dr. Gurmeet_data cleaning_dc.xlsx
# labData: dr. Gurmeet Penelitian CAP Laboratory Data 11022024 igaiw_wt_dc.xlsx
#          Data streptococcus_Dr. Gurmeet_data cleaning_dc.xlsx SHEET analysis pos spn
# genData: supplementary genome_dc.xlsx
# radiology data was given to me on 1/4/2026

library(tidyverse)
source("global/obj.R")
source("global/fun.R")

# data cleaning
dat <- readxl::read_excel("raw_data/Data streptococcus_Dr. Gurmeet_data cleaning_dc.xlsx") %>% 
  janitor::clean_names() %>% 
  
  # omit personal identity & old IDs
  dplyr::select(-no,
                -new_no,
                -nama,
                -id, -id_lab, -revised_lab_id) %>% 
  dplyr::filter(!stringr::str_detect(dc_id, "invalid")) %>%

  # classify columns based on data source
  dplyr::rename_with(~ mask_raw[.x], .cols = all_of(names(mask_raw))) %>%
  dplyr::select(-contains(c("workLab_",
                            "workWGS_"))
                ) %>%

  # fix boolean values (validated on 27/2/2026; focused on yes)
  dplyr::mutate(
    across(where(is.character), tolower),
    across(where(is.character), ~ case_when(
      . %in% c("v", "ya", "iya") ~ "1",
      . %in% c("-", "--",
               "tidak", "tidak ada", "tidak diketahui",
               "N/A") ~ "0", # omit | is.na(.) because too much data missing
      TRUE ~ as.character(.)
      )),

    hos_vacc_covid = ifelse(is.na(hos_vacc_covid), 0,
                            hos_vacc_covid)
                ) %>%
  glimpse()


# Detect duplicated IDs
report_duplicated_id <- detoxdats::report_duplicate(df = dat,
                                                    column_id = id,
                                                    other_iden = hos_admission_date) %>% 
  # view() %>% 
  glimpse()

# test unique values
report_unique_values <- detoxdats::report_uniqval(df = dat) %>%
  # view() %>% 
  glimpse()


# pending
df_epi_clean <- dat %>% 
  dplyr::transmute(
    id = id,
    age = as.numeric(age),
    ageCat = ifelse(age >= 65, "65+ (Elderly)", "18-64 (Adult)"),
    sex = ifelse(sex == "p", "Female",
                 "Male"),
    hos_symptom_onset = as.numeric(hos_symptom_onset),
    hos_admission_date = as.Date(as.numeric(hos_admission_date),
                                 origin = "1899-12-30"),
    hos_admission_time =
      hms::as_hms(as.numeric(hos_admission_time) %% 86400),
    # hos_admission_time = as.POSIXct(
    #   as.numeric(hos_admission_time) * 86400,
    #   origin = "1899-12-30",
    #   tz = "UTC"
    # ),
    hos_bmi = as.numeric(hos_bmi),
    
    # too much missing values
    # epi_smoking = epi_smoking,
    epi_smoking_history = case_when(
      str_detect(epi_smoking, "berhenti") ~ "Previously smoking", 
      str_detect(epi_smoking, c("ya|/hari|bung|bat")) |
        epi_smoking == "1" |
        epi_smoking == "5" |
        epi_smoking == "44685" ~ "Currently smoking",
      epi_smoking == "0" ~ "Never smoked",
      TRUE ~ "Unknown"
    ),
    epi_smoking_day_n = ifelse(
      current_smoking_bungkus_20batang_perhari != 0, as.numeric(current_smoking_bungkus_20batang_perhari),
      NA_integer_
    ),
    
    # antibiotic treatment
    hos_abx_start_date = as.Date(as.numeric(hos_abx_start_date),
                                 origin = "1899-12-30"),
    hos_abx_start_time = 
      hms::as_hms(as.numeric(hos_admission_time) %% 86400),
    # hos_abx_start_time = as.POSIXct(
    #   as.numeric(hos_abx_start_time) * 86400,
    #   origin = "1899-12-30",
    #   tz = "UTC"
    # ),
    hos_abx_type = hos_abx_type, # will be cleaned later
    hos_abx_bool_90d = convert_bool(hos_abx_90d),
    # hos_abx_90d = hos_abx_90d,
    
    # convert to bool
    # Treatments
    hos_treat_bool_dm = convert_bool(hos_treat_dm),
    # hos_treat_dm = hos_treat_dm, # value check
    hos_treat_bool_htn = convert_bool(hos_treat_htn),
    # hos_treat_htn = hos_treat_htn # value check
    hos_treat_bool_renal = convert_bool(hos_treat_renal),
    
    # Therapies 
    hos_therapy_bool_tracheostomy = convert_bool(hos_therapy_tracheostomy),
    hos_therapy_bool_crrt = convert_bool(hos_therapy_crrt),
    hos_therapy_bool_hemodialysis = convert_bool(hos_therapy_hemodialysis),
    hos_therapy_bool_ecmo = convert_bool(hos_therapy_ecmo),
    
    # mechanical Ventilation
    hos_mv_bool_invasive = convert_bool(hos_mv_invasive),
    hos_mv_bool_noninvasive = convert_bool(hos_mv_noninvasive),
    
    # Vaccines
    hos_vacc_covid = hos_vacc_covid,
    # hos_vacc_bool_noncovid = convert_bool(hos_vacc_noncovid), all none
    
    # Comorbidities  
    hos_com_bool_diabetes = convert_bool(hos_com_diabetes),
    hos_com_bool_hypertension = convert_bool(hos_com_hypertension),
    hos_com_bool_ckd = convert_bool(hos_com_ckd),
    hos_com_bool_cad = convert_bool(hos_com_cad),
    hos_com_bool_cvd = convert_bool(hos_com_cvd),
    hos_com_bool_copd = convert_bool(hos_com_copd),
    hos_com_bool_hiv = convert_bool(hos_com_hiv),
    hos_com_bool_malignancy = convert_bool(hos_com_malignancy),
    # hos_com_malignancy = hos_com_malignancy, # sus & curiga included as 1
    hos_com_bool_autoimmune = convert_bool(hos_com_autoimmune),
    
    # Mortality
    # hos_survived = as.integer(hos_survived),
    hos_died = convert_bool(hos_died),
    
    # scores
    hos_com_psi_class = case_when(
      hos_com_psi_class == 1 ~ "I",
      hos_com_psi_class == 2 ~ "II",
      hos_com_psi_class == 3 ~ "III",
      hos_com_psi_class == 4 ~ "IV",
      hos_com_psi_class == 5 ~ "V",
      hos_com_psi_class == 6 ~ "VI",
      hos_com_psi_class == 7 ~ "VII",
      TRUE ~ "Unknown"
      ),
    hos_com_curb65 = as.numeric(hos_com_curb65),
    hos_com_cci_score = as.numeric(hos_com_cci_score),
    hos_los = as.numeric(hos_los)
    
  ) %>% 
  # view() %>% 
  glimpse()


# clean antibiotics name #######################################################
# list of all other antibiotics
name_abx <- c(
  # --- core antibiotics ---
  "ceftriaxone" = "ceftriaxone",
  "cefriaxone" = "ceftriaxone",
  "ceftriaxne" = "ceftriaxone",
  "cefri" = "ceftriaxone",
  
  "levofloxacin" = "levofloxacin",
  "levofloxacine" = "levofloxacin",
  "leveofloxacine" = "levofloxacin",
  "levoflox" = "levofloxacin",
  "lefofloxacin" = "levofloxacin",
  "lefoloxacin" = "levofloxacin",
  "levolfloxacin" = "levofloxacin",
  "levo" = "levofloxacin",
  
  "azitromicin" = "azithromycin",
  "azitromycine" = "azithromycin",
  "azitoromycine" = "azithromycin",
  "azitro" = "azithromycin",
  "azithromycin" = "azithromycin",
  "azitromicicn" = "azithromycin",
  
  "meropenem" = "meropenem",
  
  "metronidazol" = "metronidazole",
  "metronidazole" = "metronidazole",
  
  "cefepim" = "cefepime",
  "cefepime" = "cefepime",
  "cefipime" = "cefepime",
  
  "amikasin" = "amikacin",
  "amikacin" = "amikacin",
  
  "gentamycine" = "gentamicin",
  "gentamycin" = "gentamicin",
  "gentamycin" = "gentamicin",
  
  "cefazolin" = "cefazolin",
  
  "fosmicin" = "fosfomycin",
  
  "cotrimoxazole" = "cotrimoxazole",
  
  # --- combinations / beta-lactam ---
  "ampisilin" = "ampicillin_sulbactam",
  "ampicillin" = "ampicillin_sulbactam",
  "ampicillin sulbacatam" = "ampicillin_sulbactam",
  "ampicillin-sulbactam" = "ampicillin_sulbactam",
  "ampicilin sulbactam" = "ampicillin_sulbactam",
  "ampicillin - sulbactam" = "ampicillin_sulbactam",
  "ampisilin + sulbaktam" = "ampicillin_sulbactam",
  
  "cefoperazone+sulbac" = "cefoperazone_sulbactam",
  "cefoperazone" = "cefoperazone_sulbactam",
  
  "piperacilin-tazobactam" = "piperacillin_tazobactam",
  
  # --- brand/local names ---
  "vicillin" = "ampicillin_sulbactam",
  "viccilin" = "ampicillin_sulbactam",
  "viccillin" = "ampicillin_sulbactam",
  "ampicillin - sulbactam" = "ampicillin_sulbactam",
  "ampicilin sulbactam" = "ampicillin_sulbactam",
  
  "tamicil" = "amoxicillin",
  "picyn" = "ampicillin_sulbactam",
  
  # --- others ---
  "remdesivir" = "remdesivir"
)

weird_abx <- c("oral", "injeksi", "dan",
               "oat", "sx", "po", "1", "2", "x",
               "gr", "mg", "ml",
               "sulbactam", "sulbaktam", "sulbacatam", "sulbac",
               "picyn"
)

cleaned_other_abx <- df_epi_clean %>%
  dplyr::rowwise() %>% # crucial!
  dplyr::mutate(hos_abx_type_cleaned = detoxdats::detox_multitext(text = hos_abx_type,
                                                                       dict = name_abx,
                                                                       exclude = weird_abx)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    id,
    hos_abx_type,
    hos_abx_type_cleaned
  ) %>%
  # view() %>%
  glimpse()

all_unique_abx <- cleaned_other_abx %>%
  dplyr::pull(hos_abx_type_cleaned) %>%
  stringr::str_split(", ") %>%
  unlist() %>%
  unique() %>%
  sort() %>% 
  # view() %>% 
  glimpse()

df_epi_clean2 <- df_epi_clean %>% 
  dplyr::left_join(
    cleaned_other_abx %>%
      dplyr::select(-hos_abx_type,
      ) %>%
      dplyr::filter(!is.na(hos_abx_type_cleaned)) %>%
      tidyr::separate_rows(hos_abx_type_cleaned, sep = ", ") %>%
      dplyr::mutate(value = as.numeric(1)) %>%
      tidyr::pivot_wider(
        id_cols = id,
        names_from = hos_abx_type_cleaned,
        names_prefix = "hos_abx_",
        values_from = value,
        values_fill = 0,
        values_fn = length
      ) %>%
      dplyr::mutate(
        across(where(is.list), ~ as.numeric(.)),
        # across(where(is.list), ~ ifelse(map_lgl(., is.null), 0, .)),
        across(where(is.numeric), ~ replace_na(., 0))
      )
    ,
    by = "id"
  ) %>% 
  dplyr::select(-matches("gr|mg|ml-|750|x1|-"),
                -hos_abx_type,
                -hos_abx_iv
                ) %>% 
  dplyr::rename(ageGroup = ageCat) %>% 
  glimpse()

# combine epidata with radiology
df_epi_clean3 <- df_epi_clean2 %>% 
  dplyr::left_join(
    readxl::read_excel("raw_data/Form Data Penelitian S. Pneumonia (Excel) edit 2_dc.xlsx",
                       sheet = "Radiologi 355") %>% 
      dplyr::select(-No)
    ,
    by = c("id" = "dc_id")
  ) %>% 
  dplyr::mutate(
    rad_assessment = ifelse(radiology_ro_thorax == "1" |
                              radiology_ct_scan_notes == "1", 1,
                            0
    ),
    rad_infiltrate = ifelse(
      stringr::str_detect(radiology_ro_thorax_notes,
                          regex("infiltrat", ignore_case = TRUE)) |
        stringr::str_detect(radiology_ct_scan_notes,
                            regex("infiltrat", ignore_case = TRUE)),
      1, 0
    ),
    rad_infiltrate = ifelse(
      rad_assessment == 1 & is.na(rad_infiltrate), 0,rad_infiltrate
    ),
    
    rad_consolidation = ifelse(
      stringr::str_detect(radiology_ro_thorax_notes,
                          regex("konsolidasi", ignore_case = TRUE)) |
        stringr::str_detect(radiology_ct_scan_notes,
                            regex("konsolidasi", ignore_case = TRUE)),
      1, 0
    ),
    rad_consolidation = ifelse(
      rad_assessment == 1 & is.na(rad_consolidation), 0, rad_consolidation
    ),
    
    rad_pleural_effusion = ifelse(
      stringr::str_detect(radiology_ro_thorax_notes,
                          regex("efusi", ignore_case = TRUE)) |
        stringr::str_detect(radiology_ct_scan_notes,
                            regex("efusi", ignore_case = TRUE)),
      1, 0
    ),
    rad_pleural_effusion = ifelse(
      rad_assessment == 1 & is.na(rad_pleural_effusion), 0, rad_pleural_effusion
    ),
    
  ) %>% 
  dplyr::select(
    -name_etc, -NAMA, -ifelse
  ) %>% 
  # view() %>% 
  glimpse()


write.csv(df_epi_clean3, "inputs/df_epi_cap_cleaned.csv",
          row.names = F)
