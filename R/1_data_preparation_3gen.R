library(tidyverse)
source("global/obj.R")
source("global/fun.R")

# Data cleaning process for genData ############################################

df_gen_all <- readxl::read_excel("raw_data/supplementary genome_dc.xlsx") %>% 
  janitor::clean_names() %>% 
  dplyr::left_join(
    read.table("outputs/result_poppunk/qfile_filtered_19to23.txt") %>% 
      dplyr::mutate(specimen_id = V1,
                    workPoppunk_qc = "pass_qc") %>% 
      dplyr::select(specimen_id, workPoppunk_qc)
    ,
    by = c("no_isolat" = "specimen_id")
    
  ) %>% 
  
  # re-adjust column names
  dplyr::rename_all(~ paste0("workWGS_", .)) %>% 
  dplyr::rename(
    id = workWGS_dc_id,
    workPoppunk_qc = workWGS_workPoppunk_qc,
    
    workWGS_AMR_pbp1a = workWGS_pbp1a,
    workWGS_AMR_pbp2b = workWGS_pbp2b,
    workWGS_AMR_pbp2x = workWGS_pbp2x,
    workWGS_AMR_chloramphenicol = workWGS_chloramphenicol,
    workWGS_AMR_clindamycin = workWGS_clindamycin,
    workWGS_AMR_erythromycin = workWGS_erythromycin,
    workWGS_AMR_fluoroquinolones = workWGS_fluoroquinolones,
    workWGS_AMR_kanamycin = workWGS_kanamycin,
    workWGS_AMR_linezolid = workWGS_linezolid,
    workWGS_AMR_tetracycline = workWGS_tetracycline,
    workWGS_AMR_trimethoprim = workWGS_trimethoprim,
    workWGS_AMR_sulfamethoxazole = workWGS_sulfamethoxazole,
    workWGS_AMR_cotrimoxazole = workWGS_co_trimoxazole,
    workWGS_AMR_amoxicillin = workWGS_amoxicillin,
    workWGS_AMR_ceftriaxone = workWGS_ceftriaxone,
    workWGS_AMR_cefotaxime = workWGS_cefotaxime,
    workWGS_AMR_cefuroxime = workWGS_cefuroxime,
    workWGS_AMR_meropenem = workWGS_meropenem,
    workWGS_AMR_penicillin = workWGS_penicillin
  ) %>% 
  
  # readjust names for antibiotics
  dplyr::mutate(
    id = tolower(id),
    
    
    across(
    .cols = contains("AMR_"),
    .fns = ~ case_when(
      str_detect(.x, "cat_pC194") ~ "R (cat_pC194)",
      str_detect(.x, "ermB") ~ "R (ermB)",
      str_detect(.x, "ermB;mefA_10") ~ "R (ermB;mefA_10)",
      
      str_detect(.x, "tetM_2") ~ "R (tetM_2)",
      str_detect(.x, "tetM_12") ~ "R (tetM_12)",
      str_detect(.x, "tetM_13") ~ "R (tetM_13)",
      str_detect(.x, "tetM_8") ~ "R (tetM_8)",
      
      str_detect(.x, "folA_I100L") ~ "R (folA_I100L)",
      str_detect(.x, "folP_aa_insert_57-70") ~ "R (folP_57-70)",
      str_detect(.x, "folA_I100L; folP_aa_insert_57-70") ~ "R (folA_I100L; folP_aa_insert_57-70)",
      
      # str_detect(tolower(.x), "resistant") ~ "R",
      # str_detect(tolower(.x), "sensitive") ~ "S", 
      # str_detect(tolower(.x), "intermediate") ~ "I",
      
      tolower(.x) %in% c("sensitive", "s") ~ "S",
      tolower(.x) %in% c("resistant", "r") ~ "R",
      tolower(.x) %in% c("intermediate", "i") ~ "I",
      tolower(.x) %in% c("none", "nf", "nfnf", "", "null", "null/null", "nf/nf", "nfnull/null") | 
        is.na(.x) ~ "NF",
      
      tolower(.x) %in% c("sensitive / intermediate", "s/i") ~ "S/I",
      tolower(.x) %in% c("sensitive / resistant", "s/r") ~ "S/R",
      tolower(.x) %in% c("intermediate / resistant", "i/r") ~ "I/R",
      tolower(.x) %in% c("sensitive / sensitive", "s/s") ~ "S/S",
      
      TRUE ~ .x
      ))
    ) %>% 
  # grouping
  dplyr::mutate(
    workWGS_AMR_ceftriaxone_nonMeningitis = case_when(
      is.na(workWGS_AMR_ceftriaxone) ~ NA_character_,
      workWGS_AMR_ceftriaxone == "NF" ~ "NF",
      TRUE ~ sub("/.*", "", workWGS_AMR_ceftriaxone)
    ),
    workWGS_AMR_ceftriaxone_meningitis = case_when(
      is.na(workWGS_AMR_ceftriaxone) ~ NA_character_,
      workWGS_AMR_ceftriaxone == "NF" ~ "NF",
      TRUE ~ sub(".*/", "", workWGS_AMR_ceftriaxone)
    ),
    
    workWGS_AMR_cefotaxime_nonMeningitis = case_when(
      is.na(workWGS_AMR_cefotaxime) ~ NA_character_,
      workWGS_AMR_cefotaxime == "NF" ~ "NF",
      TRUE ~ sub("/.*", "", workWGS_AMR_cefotaxime)
    ),
    workWGS_AMR_cefotaxime_meningitis = case_when(
      is.na(workWGS_AMR_cefotaxime) ~ NA_character_,
      workWGS_AMR_cefotaxime == "NF" ~ "NF",
      TRUE ~ sub(".*/", "", workWGS_AMR_cefotaxime)
    ),
    
    workWGS_AMR_penicillin_nonMeningitis = case_when(
      is.na(workWGS_AMR_penicillin) ~ NA_character_,
      workWGS_AMR_penicillin == "NF" ~ "NF",
      TRUE ~ sub("/.*", "", workWGS_AMR_penicillin)
    ),
    workWGS_AMR_penicillin_meningitis = case_when(
      is.na(workWGS_AMR_penicillin) ~ NA_character_,
      workWGS_AMR_penicillin == "NF" ~ "NF",
      TRUE ~ sub(".*/", "", workWGS_AMR_penicillin)
    ),
    
    workWGS_AMR_class_cephalosporins = case_when(
      workWGS_AMR_ceftriaxone_nonMeningitis == "R" | workWGS_AMR_cefotaxime_nonMeningitis == "R" | workWGS_AMR_cefuroxime == "R" ~ "R",
      workWGS_AMR_ceftriaxone_nonMeningitis == "I" | workWGS_AMR_cefotaxime_nonMeningitis == "I" | workWGS_AMR_cefuroxime == "I" ~ "I",
      workWGS_AMR_ceftriaxone_nonMeningitis == "S" | workWGS_AMR_cefotaxime_nonMeningitis == "S" | workWGS_AMR_cefuroxime == "S" ~ "S",
      TRUE ~ "NF"
    ),
    workWGS_AMR_class_penicillins = case_when(
      workWGS_AMR_penicillin_nonMeningitis == "R" | workWGS_AMR_amoxicillin == "R" ~ "R",
      workWGS_AMR_penicillin_nonMeningitis == "I" | workWGS_AMR_amoxicillin == "I" ~ "I",
      workWGS_AMR_penicillin_nonMeningitis == "S" | workWGS_AMR_amoxicillin == "S" ~ "S",
      TRUE ~ "NF"
    ),
    workWGS_AMR_class_antifolates = case_when( # technically including workWGS_AMR_cotrimoxazole
      workWGS_AMR_trimethoprim == "R (folA_I100L)" & workWGS_AMR_sulfamethoxazole == "R (folP_57-70)" ~ "R (folA_I100L & folP_57-70)",
      workWGS_AMR_trimethoprim == "NF" & workWGS_AMR_sulfamethoxazole == "R (folP_57-70)" ~ "R (folP_57-70)",
      workWGS_AMR_trimethoprim == "R (folA_I100L)" & workWGS_AMR_sulfamethoxazole == "NF" ~ "R (folA_I100L)",
      TRUE ~ "NF"
    ),
    # define MDR flag
    workWGS_AMR_logic_class_chloramphenicol = str_detect(workWGS_AMR_chloramphenicol, "^R"),
    workWGS_AMR_logic_class_clindamycin = str_detect(workWGS_AMR_clindamycin, "^R"),
    workWGS_AMR_logic_class_erythromycin = str_detect(workWGS_AMR_erythromycin, "^R"),
    workWGS_AMR_logic_class_fluoroquinolones = str_detect(workWGS_AMR_fluoroquinolones, "^R"),
    workWGS_AMR_logic_class_kanamycin = str_detect(workWGS_AMR_kanamycin, "^R"),
    workWGS_AMR_logic_class_linezolid = str_detect(workWGS_AMR_linezolid, "^R"),
    workWGS_AMR_logic_class_tetracycline = str_detect(workWGS_AMR_tetracycline, "^R"),
    workWGS_AMR_logic_class_meropenem = str_detect(workWGS_AMR_meropenem, "^R"),
    
    workWGS_AMR_logic_class_cephalosporins = str_detect(workWGS_AMR_class_cephalosporins, "^R"),
    workWGS_AMR_logic_class_penicillins = str_detect(workWGS_AMR_class_penicillins, "^R"),    
    workWGS_AMR_logic_class_antifolates = str_detect(workWGS_AMR_class_antifolates, "^R"),
    
    workWGS_AMR_logic_class_counts = rowSums(across(starts_with("workWGS_AMR_logic_class_")), na.rm = TRUE),
    workWGS_AMR_MDR_flag = case_when(
      workWGS_AMR_logic_class_counts >= 3 ~ "MDR",
      workWGS_AMR_logic_class_counts >= 0 ~ "non-MDR",
      # workWGS_AMR_logic_class_counts == 0 ~ "non-AMR",
      TRUE ~ NA_character_
    )
  ) %>% 
  
  # ignore useless columns
  dplyr::select(
    -workWGS_no,
    -workWGS_on_paper
    
  ) %>% 
  glimpse()

# Detect duplicated IDs not required due to multiple values
report_duplicated_id <- detoxdats::report_duplicate(df = df_gen_all,
                                                    column_id = id,
                                                    other_iden = workWGS_no_isolat) %>% 
  # view() %>% 
  glimpse()

# test unique values
report_unique_values <- detoxdats::report_uniqval(df = df_gen_all) %>%
  # view() %>% 
  glimpse()

write.csv(df_gen_all, "inputs/df_gen_cap_positives_cleaned.csv",
          row.names = F)


# test select data for tree ####################################################
tree_select <- df_gen_all %>% 
  dplyr::select(workWGS_no_isolat, 
                id,
                workWGS_dc_source,
                workPoppunk_qc) %>% 
  dplyr::filter(workPoppunk_qc == "pass_qc") %>% 
  dplyr::arrange(workWGS_dc_source) %>%
  dplyr::distinct(id, .keep_all = TRUE) %>%
  # view() %>%
  glimpse()

