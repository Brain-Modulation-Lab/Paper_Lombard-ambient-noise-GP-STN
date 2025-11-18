
library(tidyverse)
library(glue)
library(readr)
library(magrittr)
library(jsonlite)
library(RColorBrewer)
library(lme4)
library(ggpubr)

  
  PATH_ANALYSIS <- "/Volumes/Nexus4/DBS/groupanalyses/task-lombard/20230524-subctx-lfp-group-PLB"
  
  # p <- r"(Y:\DBS\participants.tsv)"
  # p <- gsub("\\\\", "/", p)
  
  setwd(PATH_ANALYSIS)
  PATH_DATA <- '/Volumes/Nexus4/DBS/derivatives'
  TASK <- 'lombard'
  SESSION <- 'intraop'
  
  PATH_FIG = glue('{PATH_ANALYSIS}/fig')
  
  SUBJECTS_META = read_tsv("/Volumes/Nexus4/DBS/participants.tsv") %>%
    mutate(target = str_remove(dbs_target, "_.*")) %>%
    relocate(target, .after=dbs_target)
  
  
  COLORS <- read_tsv("../20230602-metadata-and-groupdata-PLB/lombard_colors.tsv")
  color_map_conditions <- COLORS %>% filter(category == 'condition') 
  color_map_targets <- COLORS %>% filter(category == 'target')
  color_map_disease <- COLORS %>% filter(category == 'disease') 
  # scale_color_manual(values = color_map_targets$hex, breaks=color_map_targets$name) +
  
  
  subject_ids = c(1001:1020)
  
  
  # A couple of helper functions
  summarize_all_stats <- function(.data) {
    dplyr::summarise(.data, across(where(is.numeric), .fns = 
                                     list(mean = mean,
                                          stdev = sd,
                                          median = median,
                                          q25 = ~quantile(., 0.25),
                                          q75 = ~quantile(., 0.75),
                                          max = max,
                                          min = min,
                                          n = length)))
  }
  
  
  # load a bunch of tables
  annot_name = 'produced-sentences-acoustics'
  sentences <-
    tibble(subject=paste0('DM',subject_ids)) %>%
    mutate(path_annot=glue("{PATH_DATA}/sub-{subject}/annot/sub-{subject}_ses-{SESSION}_task-{TASK}_annot-{annot_name}.tsv")) %>%
    mutate(exists_annot=file.exists(path_annot)) %>%
    filter(exists_annot) %>% 
    pull(path_annot) %>% 
    map_df(read_tsv) %>%
    rename(any_of(c("subject_id" = "sub"))) %>%
    mutate(sentence_on = onset, 
           sentence_off = onset + duration, 
           noise_type = as.factor(noise_type)) %>%
    group_by(subject_id, run_id, trial_id, noise_type)
  # sentences$noise_type <- recode(sentences$noise_type, '0' = "QUIET",  '1' = "LMBRD")
  
  sentences_save <-
    sentences %>%
    group_by(subject_id, run_id, noise_type) %>%
    arrange(intensity_praat, .by_group=TRUE) %>%
    mutate(intensity_praat_rank = row_number()) %>% 
    add_count() %>%
    mutate(intensity_praat_prct = intensity_praat_rank/n) %>%
    relocate(subject_id, run_id, trial_id, noise_type, intensity_praat, intensity_praat_prct) 
  fname = glue('A02b_produced-sentences-acoustics-all-subj-DM1001-DM1020')
  # ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=700*2.5, height=600*2.5, units='px')
  write_tsv(sentences_save, glue('{fname}.tsv'))




  