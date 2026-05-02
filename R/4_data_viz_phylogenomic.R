library(tidyverse)
library(ggtree)
library(ggtreeExtra)
source("global/fun.R")
source("global/obj.R")

df_epi_gen_pneumo <- read.csv("inputs/df_gen_cap_positives_cleaned.csv") %>% 
  dplyr::filter(workPoppunk_qc == "pass_qc") %>%
  dplyr::rename(epi_label = id) %>%
  dplyr::mutate(workWGS_no_isolat = gsub("^Streptococcus_pneumoniae_", "",
                                         workWGS_no_isolat),
                workWGS_serotype = factor(workWGS_serotype,
                                                 levels = c(
                                                   # VT
                                                   # "3", "6A/6B", "6A/6B/6C/6D", "serogroup 6",
                                                   # "14", "17F", "18C/18B", "19A", "19F", "23F",
                                                   "1", "3", "4", "5", "7F",
                                                   "6A", "6B", "9V", "14", "18C",
                                                   "19A", "19F", "23F",
                                                   # NVT
                                                   # "7C", "10A", "10B", "11A/11D", "13", "15A", "15B/15C",
                                                   # "16F", "19B", "20", "23A", "23B", "24F/24B", "25F/25A/38",
                                                   # "28F/28A", "31", "34", "35A/35C/42", "35B/35D", "35F",
                                                   # "37", "39", "mixed serogroups",
                                                   "serogroup 6", "6C", "7C",
                                                   "10", "10A", "10B", "11A", "13",
                                                   "15A", "15B", "15C","15B/15C", "16F",
                                                   "17F", "18A", "18B", "19B", "20", "20B",
                                                   "21", "23A", "23B", "23B1",
                                                   "24F", "24B/C/F", "24B/24C/24F", "serogroup 24",
                                                   "25B", "25F",
                                                   "28A", "31", "33B", "33G",
                                                   "34", "35A", "35B", "35C", "35F", "37",
                                                   "37F", "38", "39",
                                                   "NT")),
                # test individual serotype
                # workWGS_serotype = ifelse(workWGS_serotype == "NT", "NT", "others"),
                
                workWGS_genome_length = as.numeric(workWGS_genome_length),
                
                workWGS_gpsc_strain = ifelse(workWGS_gpsc_strain == "Not assigned", "not assigned",
                                             workWGS_gpsc_strain),
                
                serotype_classification_PCV13_final_decision = case_when(
                  workWGS_serotype %in% c("1", "3", "4", "5", "7F",
                                                 "6A", "6B", "9V", "14", "18C",
                                                 "19A", "19F", "23F") ~ "VT",
                  workWGS_serotype == "NT" ~ "NT",
                  is.na(workWGS_serotype) ~ NA,
                  TRUE ~ "NVT"
                ),
                
                serotype_classification_PCV13_final_decision = factor(serotype_classification_PCV13_final_decision,
                                                                      levels = c("VT", "NVT", "NT")),
                
                serotype_classification_PCV15_final_decision = case_when(
                  workWGS_serotype %in% c("1", "3", "4", "5", "7F",
                                                 "6A", "6B", "9V", "14", "18C",
                                                 "19A", "19F", "23F",
                                                 "22F", "33F") ~ "VT",
                  workWGS_serotype == "NT" ~ "NT",
                  is.na(workWGS_serotype) ~ NA,
                  TRUE ~ "NVT"
                ),
                
                workWGS_dc_source = case_when(
                  workWGS_dc_source == "darah" ~ "Blood",
                  workWGS_dc_source == "NP" ~ "Nasopharyngeal Swab",
                  workWGS_dc_source == "sputum" ~ "Sputum",
                  TRUE ~ workWGS_dc_source
                ),
                
                workWGS_dc_source = factor(workWGS_dc_source,
                                           levels = c("Nasopharyngeal Swab",
                                                      "Sputum",
                                                      "Blood")),
                
                
                # reset logic label for AMR-MDR viz
                workWGS_AMR_logic_class_chloramphenicol = case_when(
                  workWGS_AMR_chloramphenicol != "NF" ~ stringr::str_extract(workWGS_AMR_chloramphenicol,
                                                                                   "(?<=R \\().*(?=\\))"),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_chloramphenicol = case_when(
                  workWGS_AMR_logic_class_chloramphenicol == "cat_pC194" ~ "cat (pC194)",
                  workWGS_AMR_logic_class_chloramphenicol == "cat_q" ~ "catQ",
                  TRUE ~ workWGS_AMR_logic_class_chloramphenicol
                ),
                workWGS_AMR_logic_class_chloramphenicol = factor(workWGS_AMR_logic_class_chloramphenicol),
                
                workWGS_AMR_logic_class_clindamycin = case_when(
                  workWGS_AMR_clindamycin != "NF" ~ stringr::str_extract(workWGS_AMR_clindamycin,
                                                                               "(?<=R \\().*(?=\\))"),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_clindamycin = factor(workWGS_AMR_logic_class_clindamycin),
                
                workWGS_AMR_logic_class_erythromycin = case_when(
                  workWGS_AMR_erythromycin != "NF" ~ stringr::str_extract(workWGS_AMR_erythromycin,
                      "(?<=R \\().*(?=\\))"
                    ),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_erythromycin = case_when(
                  workWGS_AMR_logic_class_erythromycin == "mefA_10" ~ "mefA",
                  TRUE ~ workWGS_AMR_logic_class_erythromycin
                ),
                workWGS_AMR_logic_class_erythromycin = factor(workWGS_AMR_logic_class_erythromycin),
                
                workWGS_AMR_logic_class_fluoroquinolones = case_when(
                  workWGS_AMR_fluoroquinolones != "NF" ~ stringr::str_extract(workWGS_AMR_fluoroquinolones,
                      "(?<=R \\().*(?=\\))"
                    ),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_fluoroquinolones = case_when(
                  workWGS_AMR_logic_class_fluoroquinolones == "parC_S79Y" |
                    workWGS_AMR_logic_class_fluoroquinolones == "parC_D83N" ~ "parC",
                  TRUE ~ workWGS_AMR_logic_class_fluoroquinolones
                ),
                workWGS_AMR_logic_class_fluoroquinolones = factor(workWGS_AMR_logic_class_fluoroquinolones),
                
                workWGS_AMR_logic_class_tetracycline = case_when(
                  workWGS_AMR_tetracycline != "NF" ~ stringr::str_extract(workWGS_AMR_tetracycline,
                      "(?<=R \\().*(?=\\))"
                    ),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_tetracycline = case_when(
                  stringr::str_detect(workWGS_AMR_logic_class_tetracycline, "tetM") ~ "tet(M)",
                  stringr::str_detect(workWGS_AMR_logic_class_tetracycline, "tet_") ~ "tet",
                  workWGS_AMR_logic_class_tetracycline == "tetK" ~ "tet(K)",
                  TRUE ~ workWGS_AMR_logic_class_tetracycline
                ),
                workWGS_AMR_logic_class_tetracycline = factor(workWGS_AMR_logic_class_tetracycline),
                
                workWGS_AMR_logic_class_antifolates = case_when(
                  workWGS_AMR_class_antifolates != "NF" ~ stringr::str_extract(workWGS_AMR_class_antifolates,
                      "(?<=R \\().*(?=\\))"
                    ),
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_antifolates = case_when(
                  workWGS_AMR_logic_class_antifolates == "folA_I100L" ~ "folA",
                  workWGS_AMR_logic_class_antifolates == "folP_57-70" ~ "folP",
                  workWGS_AMR_logic_class_antifolates == "folA_I100L & folP_57-70" ~ "folA & folP",
                  TRUE ~ workWGS_AMR_logic_class_antifolates
                ),
                workWGS_AMR_logic_class_antifolates = factor(workWGS_AMR_logic_class_antifolates),
                
                # SIR format
                workWGS_AMR_logic_class_cephalosporins = case_when(
                  workWGS_AMR_logic_class_cephalosporins ~ "  Resistant",
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_cephalosporins = factor(workWGS_AMR_logic_class_cephalosporins, levels = c(" Not found", "  Resistant")),

                workWGS_AMR_logic_class_penicillins = case_when(
                  workWGS_AMR_logic_class_penicillins ~ "  Resistant",
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_penicillins = factor(workWGS_AMR_logic_class_penicillins, levels = c(" Not found", "  Resistant")),
                
                workWGS_AMR_logic_class_meropenem = case_when(
                  workWGS_AMR_logic_class_meropenem ~ "  Resistant",
                  TRUE ~ " Not found"
                ),
                workWGS_AMR_logic_class_meropenem = factor(workWGS_AMR_logic_class_meropenem, levels = c(" Not found", "  Resistant")),
                
                workWGS_AMR_MDR_flag = ifelse(workWGS_AMR_MDR_flag == "non-MDR", " Not found", " MDR"),
                workWGS_AMR_MDR_flag = factor(workWGS_AMR_MDR_flag,
                                              levels = c(" Not found", " MDR")),

                ) %>% 
  
  # combine with labWork data (phenotypic AMR tests)
  dplyr::left_join(
    read.csv("inputs/df_lab2_cap_positives_cleaned.csv") %>%
      dplyr::mutate(workLab_source = case_when(
        workLab_source == "darah" ~ "Blood",
        workLab_source == "NP" ~ "Nasopharyngeal Swab",
        workLab_source == "sputum" ~ "Sputum",
        TRUE ~ workLab_source
      ),
      
      workLab_source = factor(workLab_source,
                              levels = c("Nasopharyngeal Swab",
                                         "Sputum",
                                         "Blood")
                              ),

      across(
        ends_with("_class") | ends_with("_flag_lab"),
        ~ case_when(
          . == "Resistant" ~ "  Resistant",
          . == "Intermediate" ~ "  Intermediate",
          . == "Not found"  ~ " Not found",
          . == "MDR"        ~ " MDR",
          # is.na(.) ~ " Not tested",
          TRUE ~ .
        )
      ),
      
      # across(everything(), as.factor)
      
      ) %>% 
      dplyr::select(
        id,
        contains(c("_interp_class", "_source", "_MDR_"))
      )
    ,
    by = c("epi_label" = "id", "workWGS_dc_source" = "workLab_source")
  ) %>%
  
  # re-adjust not-tested data in labWork
  dplyr::mutate(
    across(
      ends_with("_class") | ends_with("_flag_lab"),
      ~ case_when(is.na(.) ~ " Not tested",
                  TRUE ~ .
      )
    ),
    
    across(everything(), as.factor)
  ) %>% 
  glimpse()

# Detect duplicated IDs not required due to multiple values
report_duplicated_id <- detoxdats::report_duplicate(df = df_epi_gen_pneumo,
                                                    column_id = workWGS_no_isolat,
                                                    other_iden = epi_label) %>% 
  # view() %>% 
  glimpse()

tre_raxml <- ape::read.tree("outputs/result_raxml_from_panaroo_1all/RAxML_bestTree.1_output_tree")
tre_raxml$tip.label <- gsub("^Streptococcus_pneumoniae_", "", tre_raxml$tip.label)
tre_raxml$tip.label <- gsub(".contigs_velvet$", "", tre_raxml$tip.label)

# focused on raxml tree: rearrange label coz' ggtree link label to row_names
df_epi_gen_pneumo <- 
  dplyr::left_join(
    data.frame(tre_raxml$tip.label),
    df_epi_gen_pneumo,
    by = c("tre_raxml.tip.label" = "workWGS_no_isolat")
  )

rownames(df_epi_gen_pneumo) <- tre_raxml$tip.label

# test node
ggtree(tre_raxml) + 
  geom_tiplab(size = 2) +
  geom_label2(aes(subset=!isTip, label=node),
              size=2, color="darkred", alpha=0.5)


show_raxml <- ggtree(tre_raxml,
                  layout = "rectangular",
                  # open.angle=30,
                  size=0.25
                  # aes(colour=Clade),
) %<+% 
  df_epi_gen_pneumo +
  # geom_tiplab(size = 2) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm"),
    # plot.margin = margin(0, 0, 20, 0, "mm") # additional margin for AMR-MDR flag
  ) +
  geom_tiplab(size = 2) #+
  # geom_label2(aes(subset=!isTip, label=node),
  #             size=2, color="darkred", alpha=0.5)
show_raxml

# gen tree #####################################################################
tree_gen_raxml <- show_raxml %<+%
  df_epi_gen_pneumo +
  # vaccine classification
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=df_epi_gen_pneumo$serotype_classification_PCV13_final_decision),
    width=0.001,
    offset=0.21
  ) +
  scale_fill_manual(
    name = "PCV13 Serotype Coverage",
    values = c("indianred3", "skyblue2", "deepskyblue3"), # c(col_map),
    breaks = c("VT", "NVT", "NT"),
    labels = c("VT", "NVT", "NT"),
    guide=guide_legend(keywidth=0.3, keyheight=0.3,
                       ncol=2, order=1)
  ) +
  theme(
    legend.title=element_text(size=12), 
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm"),
    # plot.margin = margin(0, 0, 20, 0, "mm") # additional margin for AMR-MDR flag
  ) +
  # serotype
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=df_epi_gen_pneumo$workWGS_serotype),
    width=0.001,
    offset=0.05
  ) +
  scale_fill_viridis_d(
    name = "Serotype",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 4, order = 2)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm"),
    # plot.margin = margin(0, 0, 20, 0, "mm") # additional margin for AMR-MDR flag
  ) +
  # GPSC
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=df_epi_gen_pneumo$workWGS_gpsc_strain), #gpsc_dominant_filter),
    width=0.001,
    offset=0.05
  ) +
  scale_fill_viridis_d(
    name = "GPSC",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 4, order = 3)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm"),
    # plot.margin = margin(0, 0, 20, 0, "mm") # additional margin for AMR-MDR flag
  ) +
  # Source
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=df_epi_gen_pneumo$workWGS_dc_source), #gpsc_dominant_filter),
    width=0.001,
    offset=0.05
  ) +
  scale_fill_viridis_d(
    name = "Sample Source",
    option = "C",
    direction = -1,
    guide = guide_legend(keywidth = 0.3, keyheight = 0.3,
                         ncol = 1, order = 4)
  ) +
  theme(
    legend.title=element_text(size=12),
    legend.text=element_text(size=9),
    legend.spacing.y = unit(0.02, "cm"),
    # plot.margin = margin(0, 0, 20, 0, "mm") # additional margin for AMR-MDR flag
  )

# png("pictures/phylo_raxml_1epiTree.png",
#     width = 30, height = 20, units = "cm", res = 800)
tree_gen_raxml
# dev.off()


# AMR tree ver2 ################################################################
filtered_df <- df_epi_gen_pneumo %>% 
  dplyr::select(
    # serotype_classification_PCV13_final_decision,
    # workWGS_serotype,
    # tre_raxml.tip.label,
    contains("workWGS_AMR_logic_class_"),
    -workWGS_AMR_logic_class_counts,
    -workWGS_AMR_logic_class_fluoroquinolones, # 0 Resistant
    -workWGS_AMR_logic_class_kanamycin, # 0 Resistant
    -workWGS_AMR_logic_class_linezolid, # 0 Resistant
    workWGS_AMR_MDR_flag,
    contains("_interp_class"),
    workLab_MDR_flag_lab
  ) %>% 
  dplyr::rename_with(
    ~ sub("workWGS_AMR_logic_class_|workWGS_AMR_|gene_present_absent_|workLab_", "", .)
  ) %>% 
  dplyr::transmute(
    "MEM (g)" = meropenem,
    "MEM (p)" = meropenem_interp_class,
    
    "PEN (g)" = penicillins,
    "PEN (p)" = penicillins_interp_class,
    
    "Ceph (g)" = cephalosporins,
    "Ceph (p)" = cephalosporins_interp_class,
    
    "CHL (g)" = chloramphenicol,
    "CHL (p)" = chloramphenicol_interp_class,
    
    # the big three
    "TET (g)" = tetracycline,
    "TET (p)" = tetracycline_interp_class,
    
    "FOL (g)" = antifolates,
    "FOL (p)" = antifolates_interp_class,
    
    "MDR (g)" = MDR_flag,
    "MDR (p)" = MDR_flag_lab
    
  ) %>% 
  glimpse()

rownames(filtered_df) <- tre_raxml$tip.label

all_labels <- unique(unlist(lapply(filtered_df, as.character)))
manual <- c(" Not found" = "palegoldenrod",
            " Not tested" = "lightgoldenrodyellow")

# "lightgoldenrod",
# "lightgoldenrodyellow",
others <- setdiff(all_labels, names(manual))

auto_col <- scales::hue_pal()(length(others))
names(auto_col) <- others
final_col <- c(manual, auto_col)


png("pictures/phylo_raxml_2AMR_all.png",
    width = 24, height = 12, units = "cm", res = 800)
library(ggnewscale)
p2 <- tree_gen_raxml + ggnewscale::new_scale_fill()
ggtree::gheatmap(p2, filtered_df,
                 offset=0.005, width=0.7, font.size=1.7, 
                 colnames_angle=-30, hjust=0) +
  scale_fill_manual(values = final_col,
                    na.value = "white",
                    drop = F,
                    # na.translate = FALSE,
                    name = "Antimicrobial Resistant\n(g) = Genomic Prediction\n(p) = Phenotypic Test",
                    guide = guide_legend(nrow = 5),
                    labels = function(x) parse(text = parse_italiced_legends(x))
                    
  ) +
  # readjust legends
  theme(
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )
dev.off()



