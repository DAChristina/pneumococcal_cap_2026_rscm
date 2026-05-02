library(tidyverse)
source("global/obj.R")
source("global/fun.R")

# Data cleaning process for workLab ############################################

labWork1 <- readxl::read_excel("raw_data/dr. Gurmeet Penelitian CAP Laboratory Data 11022024 igaiw_wt_dc.xlsx") %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(dc_excuded = ifelse(is.na(dc_excuded), "included",
                                    "excluded")
                ) %>% 
  dplyr::filter(dc_excuded == "included") %>%
  
  # omit personal identity & old IDs
  dplyr::select(-no,
                -sample_id,
                -dc_excuded,
                -np_optochin_susceptibility, # equivalent to np
                -np_bile_solubility, # all NA
                ) %>%

  # fix boolean values (validated on 27/2/2026; focused on yes)
  dplyr::mutate(
    across(where(is.character), ~ case_when(
      . %in% c("+") ~ "1",
      . %in% c("-") ~ "0",
      TRUE ~ as.character(.)
    )),
    
    # re-adjust some values to bool
    # np_optochin_susceptibility = ifelse(np_optochin_susceptibility == "S"),
    # np_bile_solubility = ifelse(np_bile_solubility),
  ) %>%
  
  # re-adjust column names
  dplyr::rename_all(~ paste0("workLab_", .)) %>%
  dplyr::rename(
    id = workLab_dc_id
  ) %>% 
  glimpse()

# Detect duplicated IDs
report_duplicated_id <- detoxdats::report_duplicate(df = labWork1,
                                                    column_id = id,
                                                    other_iden = workLab_date_of_collection) %>% 
  # view() %>% 
  glimpse()

# test unique values
report_unique_values <- detoxdats::report_uniqval(df = labWork1) %>%
  # view() %>% 
  glimpse()


write.csv(labWork1, "inputs/df_lab1_cap_cleaned.csv",
          row.names = F)


# focused on verified pneumoPositives ##########################################
labWork2 <- readxl::read_excel("raw_data/data streptococcus_Dr. Gurmeet_data cleaning_dc.xlsx",
                               sheet = "analysis pos spn") %>% 
  janitor::clean_names() %>% 
  
  # fix boolean values (validated on 27/2/2026; focused on yes)
  dplyr::mutate(
    across(where(is.character), ~ case_when(
      . %in% c("pos") ~ "1",
      . %in% c("neg") ~ "0",
      TRUE ~ as.character(.)
    )),
    
    source = ifelse(is.na(cg_mlst), "unknown", "NP")
  ) %>% 
  
  # test epi; data compatible with ID.
  # dplyr::left_join(
  #   readxl::read_excel("raw_data/Data streptococcus_Dr. Gurmeet_data cleaning_dc.xlsx") %>% 
  #     janitor::clean_names()
  #   ,
  #   by = c("dc_id")
  # ) %>% 
  # dplyr::select(
  #   #contains(c(".x", ".y")),
  #               nama.x, nama.y) %>% 
  # view() %>% 
  dplyr::select(dc_id,
                source,
                revised_serotype,
                mslt,
                cg_mlst,
                serotype,
                contains(c("_mic", "interp"))
                ) %>%
  dplyr::rename_all(~ paste0("workLab_", .)) %>%
  dplyr::rename(id = workLab_dc_id,
                workLab_revised_mlst = workLab_mslt,
                ) %>% 
  
  # adjust AMR class & MDR
  dplyr::mutate(
    workLab_erythromycin_interp_class = ifelse(
      workLab_erythromycin_interp == "R" |
        workLab_azithromycin_interp == "R",
      "Resistant",
      "Not found"
    ),
    
    workLab_meropenem_interp_class = ifelse(
      workLab_meropenem_interp == "R" |
        workLab_ertapenem_interp == "R",
      "Resistant",
      "Not found"
    ),
    
    workLab_penicillins_interp_class = ifelse(
      workLab_penicillin_interp == "R" |
        workLab_amoxicillin_clavulanic_acid_2_1_interp == "R",
      "Resistant",
      "Not found"
    ),
    
    workLab_cephalosporins_interp_class = ifelse(
      workLab_cefuroxime_interp == "R" |
        workLab_ceftriaxone_interp == "R" |
        workLab_cefotaxime_interp == "R" |
        workLab_cefepime_interp == "R",
      "Resistant",
      "Not found"
    ),
    
    workLab_clindamycin_interp_class = ifelse(
      workLab_clindamycin_interp == "R",
      "Resistant",
      "Not found"
    ),
    
    workLab_chloramphenicol_interp_class = ifelse(
      workLab_chloramphenicol_interp == "R",
      "Resistant",
      "Not found"
    ),
    
    workLab_tetracycline_interp_class = ifelse(
      workLab_tetracycline_interp == "R",
      "Resistant",
      "Not found"
    ),
    
    workLab_antifolates_interp_class = case_when(
      workLab_trimethoprim_sulfamethoxazole_interp == "R" ~ "Resistant",
      workLab_trimethoprim_sulfamethoxazole_interp == "I" ~ "Intermediate",
      TRUE ~ "Not found"
    ),
    
    workLab_MDR_flag_lab = ifelse(
      rowSums(across(starts_with("workLab_") & ends_with("_interp_class"),
                     ~ . == "Resistant")) >= 3,
      "MDR",
      "Not found"
    )
    
  ) %>% 
  glimpse()


write.csv(labWork2, "inputs/df_lab2_cap_positives_cleaned.csv",
          row.names = F)

