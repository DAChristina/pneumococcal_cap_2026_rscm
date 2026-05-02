library(tidyverse)
library(ComplexUpset)
library(ggsankey)
source("global/fun.R")
source("global/obj.R")

# test comorbidities logic count
logic <- read.csv("inputs/df_epi_cap_cleaned.csv") %>% 
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
  dplyr::select(id, contains(c("_com_bool", "workLab"))) %>% 
  dplyr::mutate(
    hos_com_bool_count = rowSums(
      across(contains("_com_bool_"), ~ as.numeric(.x == 1)), 
      na.rm = TRUE
    ),
  ) %>% 
  glimpse()

df_signs <- read.csv("inputs/df_epi_cap_cleaned.csv") %>%
  dplyr::select(
    id,
    age,
    ageGroup,
    sex,
    contains(c("hos_treat_",
               "hos_therapy_",
               "hos_mv_",
               "hos_vacc_",
               "hos_com_",
               "rad_"
               )),
    -hos_com_curb65,
    -hos_com_cci_score,
    -hos_vacc_covid,
    -contains(c("class", "treat", "mv", "therapy"))
  ) %>% 
  
  dplyr::mutate(
    across(5:ncol(.), as.numeric)
    
    # readjust column naming for publication
    
    ) %>%
  as.data.frame() %>% 
  
  # combine with labWork
  dplyr::left_join(
    read.csv("inputs/df_lab1_cap_cleaned.csv") %>% 
      dplyr::transmute(
        id = id,
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
    across(where(is.integer), as.double),
    Test = case_when(
      workLab_systemic == 0 & workLab_carriage == 0 ~ "Negative",
      workLab_systemic == 1 & workLab_carriage == 1 ~ "Positive (systemic and carriage)",
      workLab_systemic == 1 & workLab_carriage == 0 ~ "Positive (systemic)",
      workLab_systemic == 0 & workLab_carriage == 1 ~ "Positive (carriage)",
      TRUE ~ "Negative"
      ),
    ) %>% 
  dplyr::rename(!!!col_pub) %>%
  # dplyr::filter(Test != "Negative") %>% 
  glimpse()

png(file = "pictures/epi_combinations_comorbidities.png",
    width = 29, height = 15, unit = "cm", res = 600)
target <- df_signs %>% 
  # non-laboratory results
  dplyr::select(
    -contains("Positive"),
    -age, -sex, -ageGroup,
    -Test,
    -`Chest radiography`,
    -Infiltrates,
    -Consolidation,
    -`Pleural effusion`
  ) %>% 
  colnames()

ComplexUpset::upset(
  df_signs,
  intersect = target,
  # encode_sets=FALSE, # for annotate() to select the set by name disable encoding
  base_annotations = list(
    "Test result" = intersection_size(
      counts = FALSE,
      mapping = aes(fill=Test),
      text = list(
        check_overlap=TRUE,
        vjust = -0.1,
        hjust = -0.1
      ))), # +
  # scale_fill_manual(values=col_map)),
  annotations = list(
    "Age\n(%)" = ggplot(mapping=aes(fill=ageGroup)) +
      geom_bar(stat='count', position='fill') +
      # scale_fill_manual(values=col_map) +
      scale_y_continuous(labels=scales::percent_format())
  ),
  # queries=list(
  #   upset_query(
  #     intersect=c('pneumonia', 'pleural_effusion'),
  #     color='red',
  #     fill='red',
  #     only_components=c('intersections_matrix', 'Intersection size')
  #   ),
  # ),
  wrap = F,
  mode = "distinct",
  width_ratio=0.1,
  min_degree=1, # min combination
  set_sizes=upset_set_size(mapping = aes(),
                           geom = geom_bar(aes(fill=`Overall Positive`,
                                               x=group),
                                           width = 0.6),
                           position = "left",
                           filter_intersections = T) +
  # scale_fill_manual(values=col_map) +
  geom_text(aes(label=..count..),
            hjust=1.1, size = 3,
            stat='count') +
  expand_limits(y=200) +
  theme(legend.position = "none",
        axis.text.x=element_text(angle=90))
)
dev.off()


