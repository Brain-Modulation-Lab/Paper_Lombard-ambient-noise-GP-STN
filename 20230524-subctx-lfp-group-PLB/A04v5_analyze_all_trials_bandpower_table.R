# A04v5 need clean start of figures
# A04v5 need clean start of figures--LFPs weren't log(10) in the last run of A04v4

# TODO upgrade pipeline to have an epochSTIMULUS in addition to epochSPEECH
# - co-fitting the time model and the noise model in LFP
# - run full R pipeline (permutation testing) for auditory electrode subset (on turbo) 
# - correlate loudness level with BGA in cortical electrodes which show a NOISExSPEECH effect
# - take subjects with at least one NOISExSPEECH contact and run dPCA on them each individually, or pull them all together
# - we had two main hypotheses about global state changes induced by the lombard effect after discovering that there were few electrode-level changes
# - from Dichter 2018 pitch control... correlate pitch () and intensity() to see if we observe a dLMC

# COMPLETE
# - plot power across all electrodes in periods with NOISE vs periods of QUIET... is the effect a global state-change?
# - co-fitting the time model and the noise model in LFP


# 2025 10 03 goals LFP: A04v5_00_trial-beh-variability_tlock-timewarp_t0.pdf
# produced current onset time densities F2 



# Setup ---- 
# library(signal)
# detach("package:signal", unload = TRUE)


library(tidyverse)
library(glue)
library(readr)
library(magrittr)
library(jsonlite)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(ggpubr)
library(plotly)
library(broom)
library(arrow)
library(testit)
library(boot)
# install.packages("lmerTest")
# install.packages("plotly")
# install.packages("broom")


if (Sys.info()['nodename']=="MGB032340") {       PATH_ROOT = "/Volumes/Nexus4/DBS"
} else if (Sys.info()['nodename']=="NSSBML01") { PATH_ROOT = "Y:/DBS"}

PATH_DATA <- glue('{PATH_ROOT}/derivatives')

TASK <- 'lombard'
SESSION <- 'intraop'

MODALITY <- "SU" # SU or LFP, or LFP-CTX
SUBCTX_LFP_REGION <- 'BG' # BG for SU/subcort LFP, 'BG-epochSTIMULUS' for special config
CTX_FLAG <- FALSE
CTX_LFP_REGION <- "" # pos-AC, pre-CG-post-CG, DLPFC-IFG-POC

RM_OUTLIERS_FLAG <- TRUE

if (MODALITY=="LFP") {
  PATH_ANALYSIS <- glue("{PATH_ROOT}/groupanalyses/task-lombard/20230524-subctx-lfp-group-PLB")
  if (CTX_FLAG) {
    PATH_FIG <- glue('{PATH_ANALYSIS}/fig/channel_summary_ctx_{CTX_LFP_REGION}')
    PATH_DATA_INPUT <-  glue('{PATH_ANALYSIS}/data/tlock-audio-on_chs-{CTX_LFP_REGION}')
  } else {
    PATH_FIG <- glue('{PATH_ANALYSIS}/fig/channel_summary_{SUBCTX_LFP_REGION}')
    if (SUBCTX_LFP_REGION=='BG') {PATH_FIG <- glue('{PATH_ANALYSIS}/fig/channel_summary')}
    
    
    PATH_DATA_INPUT <-  glue('{PATH_ANALYSIS}/data/tlock-audio-on_chs-{SUBCTX_LFP_REGION}')
    if (SUBCTX_LFP_REGION=='BG' | SUBCTX_LFP_REGION=='BG-epochSTIMULUS') {PATH_DATA_INPUT <-  glue('{PATH_ANALYSIS}/data/tlock-audio-on_chs-macro-ecogL-ecogL-ecogL-audio-envaudio')}
  }
} else {
  PATH_ANALYSIS <- glue("{PATH_ROOT}/groupanalyses/task-lombard/20231122-firing-rate_rep/latane") # single unit} 
  PATH_FIG <- glue('{PATH_ANALYSIS}/fig/channel_summary')
  if (SUBCTX_LFP_REGION=='BG-epochSTIMULUS') {PATH_FIG <- glue('{PATH_ANALYSIS}/fig/channel_summary_{SUBCTX_LFP_REGION}')}
  PATH_DATA_INPUT <-  glue('/Volumes/Nexus4/DBS/groupanalyses/task-lombard/20231122-firing-rate_rep/results/output/rate/loghash_20250219140301/sentence')
}
setwd(PATH_ANALYSIS)
# {PATH_ANALYSIS}/data/{path_data_subfolder}



SUBJECTS_META <- read_tsv(glue("{PATH_ROOT}/participants.tsv")) %>%
  mutate(target = str_remove(dbs_target, "_.*")) %>%
  relocate(target, .after=dbs_target)

COLORS <- read_tsv(glue("{PATH_ROOT}/groupanalyses/task-lombard/20230602-metadata-and-groupdata-PLB/lombard_colors.tsv"))
color_map_conditions <- COLORS %>% filter(category == 'condition') 
color_map_targets <- COLORS %>% filter(category == 'target')
color_map_disease <- COLORS %>% filter(category == 'disease') 

color_map_targets_simple <- color_map_targets %>% mutate(name = case_when(name=="GPi" ~ "GP",
                                                                       name=="VIM" ~ "other",
                                                                       .default = name))

COLORS <- setNames(COLORS$hex, COLORS$name)


# theme_set(theme_grey(base_size = 10))
# theme_set(theme_grey(base_size = 18, axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))) 
theme_set(theme_bw())

subject_ids = c(1001:1039)


# A couple of helper functions
summarize_all_stats <- function(.data) {
  dplyr::summarise(.data, across(where(is.numeric), .fns = 
                                   list(mean = mean,
                                        stdev = sd,
                                        median = median,
                                        q25 = ~quantile(., 0.25),
                                        q75 = ~quantile(., 0.75),
                                        q05 = ~quantile(., 0.05),
                                        q95 = ~quantile(., 0.95),
                                        max = max,
                                        min = min,
                                        n = length,
                                        sum = sum)))
}

# testing: get_nt_nc_ns_counts(tmp_plot2, "")
get_nt_nc_ns_counts <- function(df, groupby) {
  df %>%
    group_by(!!sym(groupby)) %>%
    # ungroup() %>%
    select(subject_id, run_id, channel, trial_id) %>% 
    distinct() %>% 
    summarize(nc = n_distinct(interaction(subject_id, run_id, channel)), 
              ns = n_distinct(interaction(subject_id)), 
              nr = n_distinct(interaction(subject_id, run_id)), 
              nt = n_distinct(interaction(subject_id, run_id, trial_id)))
}

plot_events <- function() {
  list(
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8),
    geom_vline(xintercept = 2, color = "black", linetype = "dashed", linewidth = 0.8),
    geom_vline(xintercept = 3, color = "black", linetype = "dashed", linewidth = 0.8)
    # annotate("text", x = 0, y = Inf, label = "STIM on", vjust = 2, hjust = -0.2, size = 3),
    # annotate("text", x = 2, y = Inf, label = "STIM off", vjust = 2, hjust = -0.2, size = 3),
    # annotate("text", x = 3, y = Inf, label = "SPEECH on", vjust = 2, hjust = -0.2, size = 3)
  )
}

vars_channel_id <- c("subject_id", "run_id", "channel")
vars_trial_id <- c("subject_id", "run_id", "trial_id")

# Load tables ----
electrodes = read_tsv(glue('{PATH_ROOT}/groupanalyses/task-lombard/20230328-subctx-ctx-group-coverage-PLB/A02b_electrodes-all-subj-DM1001-DM1039-with-lombard-run-id.tsv'))
if (Sys.info()['nodename']=="MGB032340") {      hcplabels2regions <- read_tsv('/Volumes/Nexus/Resources/HCPMMP1-Labeling-Atlas/HCPMMP1-labels2areas.tsv') %>% rename(HCPMMP1_label_1=label) 
} else if (Sys.info()['nodename']=="NSSBML01") {hcplabels2regions <- read_tsv('Z:/Resources/HCPMMP1-Labeling-Atlas/HCPMMP1-labels2areas.tsv') %>% rename(HCPMMP1_label_1=label)}

electrodes <-
  electrodes %>%
  left_join(hcplabels2regions, by=join_by(HCPMMP1_label_1)) %>%
  filter(str_detect(name, 'macro|micro')) %>%
  left_join(SUBJECTS_META %>% select(subject_id, target)) %>%
  mutate(region_clean = case_when(str_detect(DISTAL_label_1, 'STN') ~ 'STN',
                                  str_detect(DISTAL_label_1, 'GPi') ~ 'GPi',
                                  str_detect(DISTAL_label_1, 'GPe') ~ 'GPe',
                                  str_detect(DISTAL_label_1, 'Striatum') ~ 'striatum',
                                  str_detect(DISTAL_label_1, 'medullary_lamina') ~ 'pallidal lamina',
                                  str_detect(DISTAL_label_1, 'Substantia_nigra') ~ 'SNR',
                                  str_detect(DISTAL_label_1, 'incerta') ~ 'ZI',
                                  str_detect(DISTAL_label_1, 'capsule') &  target=='STN'~ 'STN-capsule',
                                  str_detect(DISTAL_label_1, 'capsule') &  target=='GPi'~ 'GPi-capsule',
                                  str_detect(DISTAL_label_1, 'capsule') &  target=='VIM'~ 'VIM-capsule',
                                  str_detect(DISTAL_label_1, 'capsule') ~ 'capsule',
                                  str_detect(DISTAL_label_1, 'Thal') ~ 'thal',
                                  .default = glue('target_{target}'))) %>%
  group_by(region_clean) %>% mutate(region_clean_n = n()) %>% ungroup() %>%
  relocate(c('target', 'region_clean'), .after='DISTAL_label_1') %>%
  mutate(region_clean = factor(region_clean, levels = c('GPe', 'GPi', 'STN-capsule', 'STN'))) 


# electrodes <-
#   electrodes %>%
#   left_join(SUBJECTS_META %>% select(subject_id, dbs_target), by=join_by(subject_id)) %>%
#   mutate(target=case_when(str_detect(dbs_target, regex('STN', ignore_case = TRUE)) ~ 'STN', 
#                           str_detect(dbs_target, regex('GPi', ignore_case = TRUE)) ~ 'GPi'))


                 #                                    str_detect(dbs_target, 'GPi') ~ 'GPi', 


# Load acoustics tables to summarize
annot_name = 'produced-sentences-acoustics'
sentences <-
  tibble(subject_id_tmp=paste0('DM',subject_ids)) %>%
  mutate(path_annot=glue("{PATH_DATA}/sub-{subject_id_tmp}/annot/sub-{subject_id_tmp}_ses-{SESSION}_task-{TASK}_annot-{annot_name}.tsv")) %>%
  mutate(exists_annot=file.exists(path_annot)) %>%
  filter(exists_annot) %>%
  mutate(t = map(pull(., path_annot), function(x) read_tsv(x))) %>%
  unnest(t) %>%
  rename(subject_id = subject_id_tmp) %>% # unnest the tibbles on each row
  mutate(sentence_on = onset) %>%
  mutate(sentence_off = onset + duration) %>%
  group_by(subject_id, run_id, trial_id)



# create a time variable within-run
# this is used for modeling random effects of time on BGA in each trial
sentences_timenorm <-
  sentences %>%
  group_by(subject_id, run_id) %>%
  filter(row_number()==1) %>%
  mutate(t0_run = onset) %>%
  select(t0_run)
sentences <-
  sentences %>%
  left_join(sentences_timenorm, by=join_by(subject_id, run_id)) %>%
  mutate(time_withinrun = onset-t0_run)

# fname = "/Volumes/Nexus4/DBS/groupanalyses/task-lombard/20230524-subctx-lfp-group-PLB/data/tlock-audio-on_chs-macro-ecogL-ecogL-ecogL-audio-envaudio/sub-DM1019_ses-intraop_task-lombard_tlock-audio-on_TFR"
# eg <- read_parquet(glue('{fname.parquet}'))
# write_tsv(eg, glue('{fname}.tsv'))


# -------------- for LFPs
if (MODALITY=="LFP") {
# # A03b_all-trials-bandpower_20240626.tsv, macro and ecog
# A03b_all-trials-bandpower_20240620.tsv, macro only
# eeg = read_tsv('A03b_all-trials-bandpower_20240620.tsv') # A03b_channel_summary.tsv

# ecogelec_exclude <- # doing this for computational/space reasons
#   tibble(elec_num = c(seq(101, 161, 2), seq(201, 261, 2))) %>%
#   mutate(channel = glue("ecog_L{elec_num}")) %>%
#   pull(channel)

eeg <-
  tibble(subject_id_tmp=paste0('DM',subject_ids)) %>%
  mutate(path_annot=glue("{PATH_DATA_INPUT}/sub-{subject_id_tmp}_ses-intraop_task-lombard_tlock-audio-on_TFR.parquet")) %>%
  mutate(exists_annot=file.exists(path_annot)) %>%
  filter(exists_annot) %>%
  mutate(t = map(pull(., path_annot), function(x) read_parquet(x))) %>%
  unnest(t) %>%
  select(-path_annot, -exists_annot, -subject_id_tmp) # %>%
  # select(-one_of(ecogelec_exclude))


 # --------------  for SINGLE UNITS 
} else if (MODALITY=="SU") {
process_spike_parquet <- function(df) {
  # df <- eeg$data[[4]]
  newcols <- sub("(DM[A-Za-z0-9]+_)", "", colnames(df))
  newcols <- sub("(run[0-9]+_)", "", (newcols))
  colnames(df) <- newcols
  
  channel <- df %>% select(matches("micro")) %>% colnames(.) 
  channel <- 
    tibble(channel) %>%
    mutate(channel_orig = channel) %>%
    separate(channel, into = c("micro", "track", "unitID")) %>%
    group_by(micro, track) %>%
    mutate(index = row_number(), 
           channel = glue('{micro}_{track}_{index}'))
  
  channel_notephys <- df %>% select(!matches("micro")) %>% colnames(.) 
  channel_notephys <- tibble(channel_notephys) %>%
    rename(channel = channel_notephys) %>%
    mutate(channel_orig = channel)
  
  lookup <- 
    bind_rows(channel, channel_notephys) %>%
    mutate(subject_id = df$subject_id[[1]], 
           run_id =     df$run_id[[1]])
  
  
  # Rename the columns in the data frame
  df <- df %>%
      rename_with(~ lookup$channel[match(., lookup$channel_orig)], .cols = everything())
  
  df <- 
    df %>% 
    select(where(~ !all(is.na(.))))
  
  list(data = df, channel_map = lookup)
}  

# data_folder <- '/Volumes/Nexus4/DBS/groupanalyses/task-lombard/20231122-firing-rate_rep/results/output/rate/loghash_20250219140301/sentence'
# stop('update the folder location//debug as paths have changed')
file <- list.files(PATH_DATA_INPUT, pattern =  "*.parquet")
eeg <- 
  tibble(file) %>%
  # filter(row_number() %in% c(7, 8,9, 10 )) %>% # for testing
  mutate(subject_id = sub("sub-([A-Za-z0-9]+)_.*", "\\1", file),
         run_id =     sub(".*run-(\\d+)_.*", "\\1", file)) %>%
  mutate(path = glue('{PATH_DATA_INPUT}/{file}')) %>%
  mutate(data = map(pull(., path), function(x) read_parquet(x)),
         data_and_chanmap = map(data, process_spike_parquet), 
         data = map(data_and_chanmap, "data"),
         chanmap = map(data_and_chanmap, "channel_map"))
# aa <- eeg$data[[6]] %>% arrange(time_t0)

# testing
channel_map <- 
  eeg %>%
  select(chanmap) %>%
  unnest(chanmap)

eeg <-  
  eeg %>%
  select(data) %>%
  unnest(data) %>%
  relocate(!matches("micro")) # Move non-matching columns first

# ------------- end LFP/SU 
} 
  
eeg <- eeg %>%  
  # rename(subject_id = subject_id_tmp) %>% # unnest the tibbles on each row
  group_by(subject_id) %>%
  mutate(noise_type=case_when(beh_noise_type==0 ~ 'QUIET', 
                              beh_noise_type==1 ~ 'NOISE')) %>%
  mutate(noise_type=factor(noise_type, levels=c("QUIET", "NOISE"))) %>%
  mutate(epoch     =case_when(beh_itg==1 ~ 'BSLN', 
                              beh_sentence==1 ~ 'SPEECH', 
                              beh_audio==1 ~ 'STIMULUS'), 
         ) %>%
  group_by(subject_id, run_id, trial_id) %>%
  mutate(epoch = case_when(time_t0<0 & epoch!="BSLN" ~ NA,
                           time_t0>=-0.25 & time_t0<=0.25 & is.na(epoch) ~ "BSLN", 
                           time_t0>=0.25 & time_t0<=3.5 & is.na(epoch) ~ "RT", 
                           time_t0>=3.5 & epoch!="SPEECH" ~ NA, 
                           .default = epoch)) %>%
  mutate(epoch     = factor(epoch, levels=c("PREBSLN", "BSLN", "STIMULUS", "RT", "SPEECH", "POSTSPEECH"))) %>%
  mutate(noise_typeXepoch = interaction(noise_type, epoch)) %>%
  filter(!is.na(epoch))

# DM1022-unit422
# aa <- eeg %>% group_by(subject_id, run_id) %>% nest() %>% ungroup() %>% filter(subject_id == "DM1022")
# aa <- aa$data[[2]]  %>% select(where(~ !all(is.na(.))))

# Extend data table with further annotations ----
t0_epochwise <- # define start of epochs
  eeg %>%
  group_by(subject_id, run_id, trial_id, epoch) %>%
  arrange(time_t0, .by_group = TRUE) %>%
  filter(row_number()==1) %>%
  select(time_t0) %>%
  group_by(subject_id, run_id, trial_id) %>%
  pivot_wider(values_from = "time_t0", names_from = "epoch", names_prefix = "time_t0_") %>%
  arrange(.by_group = TRUE)

t1_epochwise <- # define end of epochs
  eeg %>%
  group_by(subject_id, run_id, trial_id, epoch) %>%
  arrange(time_t0, .by_group = TRUE) %>%
  filter(row_number()==n()) %>%
  select(time_t0) %>%
  group_by(subject_id, run_id, trial_id) %>%
  pivot_wider(values_from = "time_t0", names_from = "epoch", names_prefix = "time_t1_") %>%
  arrange(.by_group = TRUE)  
  
eeg <-
  eeg %>%
  group_by(subject_id, run_id, trial_id, epoch) %>%
  arrange(time_t0, .by_group = TRUE) %>%  # Ensure data is sorted by time within each trial and epoch
  mutate(
    epoch_duration = max(time_t0) - min(time_t0),  # Calculate the duration of each epoch within the trial
    timeperc_t0_epoch = (time_t0 - min(time_t0)) / epoch_duration,  # Calculate percentage time within epoch
    timewarp_t0_epoch = case_when(epoch=="BSLN" ~     timeperc_t0_epoch * 1, 
                                  epoch=="STIMULUS" ~ timeperc_t0_epoch * 2, 
                                  epoch=="RT" ~       timeperc_t0_epoch * 1,
                                  epoch=="SPEECH" ~   timeperc_t0_epoch * 2.3), 
    timewarp_t0       = case_when(epoch=="BSLN" ~     timewarp_t0_epoch - 1, 
                                  epoch=="STIMULUS" ~ timewarp_t0_epoch + 0, 
                                  epoch=="RT" ~       timewarp_t0_epoch + 2,
                                  epoch=="SPEECH" ~   timewarp_t0_epoch + 3), 
  ) %>%
  left_join(t0_epochwise) %>%
  left_join(t1_epochwise) %>%
  mutate(time_t0_SPEECH = time_t0 - time_t0_SPEECH, 
         time_t1_SPEECH = time_t0 - time_t1_SPEECH, 
         time_t0_STIMULUS = time_t0 - time_t0_STIMULUS, 
         time_t1_STIMULUS = time_t0 - time_t1_STIMULUS,
         time_t0_RT = time_t0 - time_t0_RT,
         time_t1_RT = time_t0 - time_t1_RT)

loudness <- # define start of epochs
  eeg %>%
  group_by(subject_id, run_id, trial_id) %>%
  filter((time_t0_SPEECH>0) & (time_t0_SPEECH<=1)) %>%
  summarise(acousticspectrum_intensity_trialwise = mean(acousticspectrum_intensity, na.rm=TRUE)) %>%
  group_by(subject_id, run_id) %>%
  mutate(perc = percent_rank(acousticspectrum_intensity_trialwise), 
         loudness_lvl = cut(perc, breaks=3, labels=c("low3rd", "mid3rd", "loud3rd"))) %>%
  group_by(subject_id, run_id, trial_id) %>%
  select(loudness_lvl)


# # sanity check loudness level coding 
# loudness %>%
#   group_by(subject_id, run_id) %>% nest() %>% filter(row_number() %in% sample.int(40, 10)) %>%
#   unnest(data) %>%
#   ggplot(aes(x=trial_id, y=acousticspectrum_intensity_trialwise, color=loudness_lvl)) +
#   geom_point() +
#   facet_grid(subject_id+run_id~., scales='free')

FCR <- 
  sentences %>%
  group_by(subject_id, run_id) %>%
  filter(!all(is.na(FCR))) %>%
  mutate(perc = percent_rank(FCR), 
         fcr_lvl = cut(perc, breaks=3, labels=c("low3rdFCR", "mid3rdFCR", "loud3rdFCR"))) %>%
  group_by(subject_id, run_id, trial_id) %>%
  select(FCR, fcr_lvl)
  
eeg <-
  eeg %>%
  group_by(subject_id, run_id, trial_id) %>%
  left_join(loudness) %>%
  left_join(FCR %>% group_by(subject_id, run_id, trial_id) %>% select(FCR, fcr_lvl))


# plot a few trials' behavioral annotation to ensure they are correct
# randsamp <- sample.int(3000, 10)
randsamp <- c(2240 , 440 , 158 ,1569 ,1621 ,2125, 2506  ,302 ,1799 , 1164)
tlock <- "time_t0" # time_t0 or timewarp_t0
eeg %>%  
  group_by(subject_id, run_id, trial_id) %>%
  nest() %>% 
  ungroup() %>%
  filter(row_number() %in% randsamp) %>% 
  unnest(data) %>%
  relocate(subject_id, run_id, trial_id, time_t0) %>%
  ggplot(aes(!!sym(tlock), epoch, color=epoch)) +
  geom_point() + 
  facet_grid(subject_id+run_id+trial_id~.)
ggsave(filename=glue('{PATH_FIG}/A04v5_00_trial-beh-variability_tlock-{tlock}.pdf'), dpi=300, width=5, height=6, units='in')


# summarize behavioral data
eeg %>% 
  # unnest(data) %>%
  mutate(STIMULUS_dur = time_t1_STIMULUS - time_t0_STIMULUS, 
         BSLN_dur     = time_t1_BSLN - time_t0_BSLN, 
         RT_dur     = time_t1_RT - time_t0_RT, 
         SPEECH_dur   = time_t1_SPEECH - time_t0_SPEECH) %>%
  group_by(subject_id, run_id) %>%
  summarize(STIMULUS_dur = mean(STIMULUS_dur), 
            BSLN_dur = mean(BSLN_dur),
            RT_dur = mean(RT_dur, na.rm=TRUE),
            SPEECH_dur = mean(SPEECH_dur, na.rm=TRUE)) %>% # print here for subject-wise
  ungroup() %>%
  summarize(STIMULUS_dur_m = mean(STIMULUS_dur,na.rm=TRUE), STIMULUS_dur_sd = sd(STIMULUS_dur,na.rm=TRUE),
            BSLN_dur_m = mean(BSLN_dur, na.rm=TRUE),        BSLN_dur_sd = sd(BSLN_dur,na.rm=TRUE),
            RT_dur_m = mean(RT_dur, na.rm=TRUE),            RT_dur_sd = sd(RT_dur,na.rm=TRUE),
            SPEECH_dur_m = mean(SPEECH_dur, na.rm=TRUE) ,   SPEECH_dur_sd = sd(SPEECH_dur,na.rm=TRUE))
# >> STIMULUS_dur_m STIMULUS_dur_sd BSLN_dur_m BSLN_dur_sd RT_dur_m RT_dur_sd SPEECH_dur_m SPEECH_dur_sd
#    1.98          0.0121       1.01      0.0151    0.961     0.264        -2.29         0.378
# Define behavioral epochs
EPOCHS_MEAN_TIMES <- data.frame(
  epoch = c("BSLN", "STIMULUS", "RT", "SPEECH"),
  start = c(-1.0, 0, 2, 3),
  end =   c(0, 2, 3, 5.3)
)

  

# ensure behavioral coding is correct--epochs are mutually exclusive 
assert(nrow(eeg %>%
                mutate(tot = beh_itg + beh_audio + beh_sentence) %>%
                filter(tot > 2))
       ==0)

eeg %>% group_by(subject_id, run_id) %>% select() %>% distinct()

# # old version--needed to create more efficient version fo 
# eeg_long <- eeg %>%
#   pivot_longer(cols = matches("macro_L|dbs_L|ecog_L|micro"), 
#                names_to = "channel", 
#                values_to = "value", 
#                values_drop_na = TRUE) %>%
#   {if (MODALITY=="LFP") {mutate(., value = log10(value))} else {.} } %>%
#   filter( !(is.na(value) | is.infinite(value)) ) %>% # this removes a lot of rows
#   filter( !(trial_id==0) & !(trial_id==1) & !(is.na(trial_id)) )
#   # distinct(subject_id, run_id, time_GTC, .keep_all = TRUE) %>%

# df <- eeg_nest_tmp$data[[9]]
not_all_na <- function(x) any(!is.na(x))
not_any_na <- function(x) all(!is.na(x))
pivot_long_efficient <- function(df) {
  # browser()
  df_notallna <-
    df %>% 
    select(where(not_all_na))
  
  df_tmp <- df_notallna %>%
    select(matches("macro_L|dbs_L|ecog_L|micro"))
  if (dim(df_tmp)[1]==0 | dim(df_tmp)[2]==0) { 
    # browser()
    return(df_notallna) 
    }
  
  # %>% # some columns are all NaN
  df_notallna %>%
    pivot_longer(cols = matches("macro_L|dbs_L|ecog_L|micro"), 
               names_to = "channel", 
               values_to = "value", 
               values_drop_na = TRUE) %>%
    filter( !(is.na(value) | is.infinite(value)) ) %>% # this removes a lot of rows
    filter( !(trial_id==0) & !(trial_id==1) & !(is.na(trial_id)) )
}

eeg_long <-
# df <- tmp$data[[5]]
# tmp <- 
  eeg %>% 
  group_by(subject_id,  run_id) %>%
  nest() %>%
  mutate(npts = map_int(data, \(df) nrow(df))) %>%
  filter(npts>1000) %>% # some SU dataframes have like have of a trial and the rest NaN
  {if (CTX_FLAG & CTX_LFP_REGION=="pre-CG-post-CG") {
      mutate(., subject_id_num = as.numeric(str_remove(subject_id, "^DM"))) %>%
      filter(subject_id_num < 1020)
  } else {.} } %>% 
  mutate(data = map(data, \(df) pivot_long_efficient(df) )) %>%
  unnest(data) %>%
  {if (MODALITY=="LFP") {mutate(., value = log10(value))} else {.} } 
  
  
# eeg_long %>%
#   group_by(subject_id, run_id, channel) %>%
#   arrange(.by_group = TRUE) %>%
#   nest()
  

eeg_long_summaryBSLN <- 
  eeg_long %>% 
  filter(epoch=="BSLN") %>% 
  group_by(subject_id, run_id, channel) %>%
  summarize(value_BSLN_run = mean(value, na.rm = TRUE))
  

# aa <- eeg_long %>% filter(subject_id=="DM1022") %>% group_by(subject_id, run_id) %>%
#   select(where(~ !all(is.na(.))))


# Subset data table for relevant analyses ---- 
eeg_nest <-
  eeg_long %>% 
  # trial-wise baseline value
  group_by(subject_id, run_id, channel, trial_id) %>%
  # left_join(eeg_long_summaryBSLN) %>% # we baseline later in fitting the time linear model
  # filter( !(is.na(value_BSLN_run) | is.infinite(value_BSLN_run)) ) %>% # this removes a lot of rows
  # mutate(value = value - value_BSLN_run) %>% # effectively baselining the data
  # grab t0 in each trial with a left join, create new column that is speech-locked
  group_by(subject_id, run_id, trial_id) %>%
  # mutate(time_t0_SPEECH  = time_t0 - time_t0_SPEECH) %>%
  # nest each channel
  group_by(subject_id, run_id, channel) %>%
  nest() 

eeg_nest_test <- eeg_nest %>%
  # filter(subject_id %in% c("DM1005", "DM1019", "DM1021")) %>%
  filter(run_id != 0) %>%
  mutate(chantype = ifelse(channel %in% c("macro_Ll", "macro_Lc", "macro_Lp", "macro_Lp"), "monopolar", 'bipolar')) %>%
  {if (str_detect(SUBCTX_LFP_REGION, 'bipolar')) { filter(., chantype=='bipolar')} else {.}}  %>%
  filter(!grepl("CA", channel)) %>% # select filter out common average 
  {if (CTX_FLAG) {.} else {filter(., !str_detect(channel, 'ecog'))}} 

  # filter(subject_id=="DM1012" & channel=="ecog_L163") 
  # filter(subject_id=="DM1003" & run_id==4 & channel=="macro_Lc")


# Optionally, remove outliers -
# df <- eeg_nest$data[[18]] # 18=== 2 DM1006     ecog_L147 
rm_outliers <- function(df) {
  df_stats <- df %>%
    group_by(trial_id) %>%
    select(value) %>%
    summarize_all_stats()
  
  # Calculate IQR thresholds for value_mean
  iqr_mean <- IQR(df_stats$value_mean)
  q1_mean <- quantile(df_stats$value_mean, 0.25)
  q3_mean <- quantile(df_stats$value_mean, 0.75)
  
  lower_bound <- q1_mean - 1.5 * iqr_mean
  upper_bound <- q3_mean + 1.5 * iqr_mean
  
  df_filtered <- df_stats %>%
    filter(
      between(value_mean, lower_bound, upper_bound)
    )
  
  df_cleaned <- df %>%
    filter(trial_id %in% df_filtered$trial_id)
  df_cleaned
}
if (RM_OUTLIERS_FLAG) {
  print('Removing outliers in eeg_nest_test')
  eeg_nest_test <-
    eeg_nest_test %>%
    mutate(data = map(data, \(df) rm_outliers(df)))
}

# # testing
# df %>%
#   rm_outliers(.) %>%
#   group_by(trial_id) %>%
#   select(value) %>%
#   summarize_all_stats()

# cohort statistics 

# summarize counts of patients



if (!CTX_FLAG) { # this section is too computationally intense for large datasets
cohort_tabulation <- 
  get_nt_nc_ns_counts(eeg_nest_test %>% unnest(), "subject_id") %>%
  left_join(SUBJECTS_META %>% select(subject_id, age, sex, diagnosis, target), 
            by=join_by(subject_id))
fname <- "recordkeeping-and-tabulation_no-qualitycontrol"
write_tsv(cohort_tabulation, file=glue('{PATH_FIG}/A04v5_00_{fname}_{MODALITY}.tsv'))
fname <- "recordkeeping-and-tabulation_no-qualitycontrol"
write_tsv(cohort_tabulation %>% 
            summarize_all_stats() %>% 
            pivot_longer(everything()), 
          file=glue('{PATH_FIG}/A04v5_00_{fname}_{MODALITY}-stats.tsv'))
cohort_tabulation %>%
  group_by(target) %>%
  summarize(ns = sum(ns), nc = sum(nc), nr = sum(nr), nt = sum(nt))
# >> for LFP
# target    ns    nc    nr    nt
# GPi       11    51    18  1048
# STN       12    61    22  1375
}



if (MODALITY=="SU") {
  eeg_nest_test <-
    eeg_nest_test %>%
    left_join(channel_map %>% select(subject_id, run_id, channel, channel_orig), 
              by = join_by(subject_id, run_id, channel))
}

# aa <- eeg_nest$data[[20]]
# DM1022-unit422
# aa <- eeg_nest %>% group_by(subject_id, run_id) %>% nest() %>% ungroup() %>% filter(subject_id == "DM1022")
# aa <- eeg_nest %>% ungroup() %>% filter(subject_id == "DM1022")
# aa <- aa$data[[1]]  %>% select(where(~ !all(is.na(.))))


# Plot time-resolved warped behavioral data ----
tlock <- 'time_t0' # timewarp_t0, time_t0
beh_subj <-
  eeg_nest_test %>%
  mutate(data = map(data, \(df) df %>% select(noise_type, acousticspectrum_intensity, timewarp_t0, time_t0))) %>%
  group_by(subject_id) %>%
  filter(row_number()==1) %>%
  unnest(data) %>%
  pivot_longer(cols = c("timewarp_t0", "time_t0"), values_to = "time", names_to = "time_type") %>%
  mutate(time = floor(time *50)/50) %>% 
  # mutate(time = cut(time, 200)) %>%
  group_by(subject_id, run_id, time_type, noise_type, time) %>%
  summarize(acousticspectrum_intensity = mean(acousticspectrum_intensity, na.rm=T)) %>%
  filter(time_type==tlock)

beh_all <-
  beh_subj %>%
  group_by(time_type, noise_type, time) %>%
  summarize(acousticspectrum_intensity = mean(acousticspectrum_intensity, na.rm=T)) %>%
  mutate(acousticspectrum_intensity = signal::sgolayfilt(acousticspectrum_intensity, p = 4, n = 51)) # p = polynomial order, n = window size (odd number)

    
ggplot(data=beh_subj, aes(x=time, y=acousticspectrum_intensity, color=noise_type)) + 
  # geom_line(aes(group=interaction(noise_type, subject_id, run_id)), alpha=0.2) +
    # geom_line(alpha=0.3) + 
    
  geom_line(data=beh_all, alpha=1, size=2) +
  facet_grid(time_type~., scales='free_x') + 
  {if (tlock=="timewarp_t0") plot_events() } +  
  scale_color_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) + 
  xlim(c(-1, 6)) + 
  ylim(c(-0.005, 0.4))

ggsave(filename=glue('{PATH_FIG}/A04v5_00_01_trial-warping_tlock-{tlock}.pdf'), dpi=300, width=3.5, height=2, units='in')

  
# Fit a linear model with just global time and intercept ---- 
fit_lm_time <- function(df) {
  # if () {browser()}
  # browser()
  df <- df %>% mutate(time_GTC2 = time_GTC^2)
  
  
  
  # model <- lm(value ~ 1 + time_GTC + time_GTC2, data = df)
  # # summary <- tidy(model)
  # # significant <- summary %>% filter(p.value < 0.05)
  # 
  # # PLB 2025 01 21 this is a bit of a cheap trick to compare different approaches 
  # # in v2, we fit time and then fit conditions...
  # # in v3, we want to fit time and conditions jointly....
  # 
  # df$value_residualtime <- residuals(model)
  # # df$value_residualtime <- df$value # cAREFUL THIS REMOVES THE EFFECT OF THE nonstationarity
  
  
  
  # --------- update 2025 05 26
  # Fit full model including noise_type
  # browser()
  model <- lm(value ~ time_GTC + time_GTC2 + noise_type, data = df)
  
  # Predict the contribution from time only (keep noise_type fixed)
  # Use model.matrix() to build design matrix and zero out noise_type terms
  X <- model.matrix(model)
  
  # Identify which columns in the design matrix correspond to time terms
  time_terms <- grepl("time_GTC", colnames(X))
  time_terms[1] <- TRUE  # keep intercept
  
  # Zero out noise_type columns in the matrix
  X_partial <- X
  X_partial[, !time_terms] <- 0
  
  # Multiply by full model coefficients
  pred_time_only <- X_partial %*% coef(model)
  # browser() # make sure this is a single vector not matrix
  
  # Subtract time-predicted values to get residuals of time
  df$value_residualtime <- df$value - as.vector(pred_time_only)
  # --------- update 2025 05 26
  
  

  # shift all data such that BSLN==0
  # 2025 04 24 added noise_type=="QUIET" to baseline QUIET only 
  df$value_residualtime = df$value_residualtime - mean(df %>% filter(epoch=="BSLN" & noise_type=="QUIET") %>% pull(value_residualtime), na.rm = TRUE)

  list(model = model, 
       data = df)
       # stat_noise_type = summary$statistic[summary$term=="noise_type"], # %>%filter(term=noise_type) %>% pull(statistic))
       # stat_timewin = summary$statistic[summary$term=="timewin"])
}

# apply linear model to each channel, collapse time series into two points per trial
lmtime_results_channelwise <- eeg_nest_test %>%
  # filter(subject_id=="DM1023") %>% # for testing purposes
  # option 1: data are summarized by trial, usually for statistics 
  mutate(data = map(data, \(df) df %>% filter(epoch=="BSLN" | 
                                                (epoch=="STIMULUS" & time_t0_STIMULUS<1) |
                                                (epoch=="SPEECH" & time_t0_SPEECH<1)) %>% # | (epoch=="SPEECH" & time_t0_SPEECH<1) | (epoch=="STIMULUS" & time_t0_STIMULUS<1)) %>%
                            group_by(trial_id, epoch, noise_type) %>%
                            summarize(value     =mean(value),
                                      time_GTC  =mean(time_GTC)))) %>% # Filter data for epochs
  # option 2: data are time-resolved, ~10 ms resolution
  # mutate(data = map(data, ~filter(.x, epoch=="BSLN" | (epoch=="SPEECH" & time_t0_SPEECH<1)))) %>% # Filter data for epochs
  # option 3: data are time-resolved, ~10 ms resolution, no filtering for epoch, usually for plotting
  # mutate(data = data) %>%
  # options end:---
  # mutate(model_and_data = map(data, fit_lm_time)) %>%
  mutate(model_and_data = map(data, possibly(fit_lm_time, otherwise = 'error'))) %>% # update 20250502
  mutate(model_time = map(model_and_data, "model"),
         data = map(model_and_data, "data")) %>%
         # stat_noise_type = map_dbl(model_and_data, "stat_noise_type"), 
         # stat_timewin = map_dbl(model_and_data, "stat_timewin")) %>%
  select(-model_and_data)

filtered_out <- 
  lmtime_results_channelwise %>%
  filter(map_lgl(data, is.null)) # there are a few channels that fail in the ecog data because all data are NAN
print(filtered_out)  

lmtime_results_channelwise <-
  lmtime_results_channelwise %>%
  anti_join(filtered_out)


# apply linear model to each channel, time-resolved for visualization
lmtime_results_channelwise_timeresolved <- eeg_nest_test %>%
  # filter(subject_id=="DM1023") %>% # for testing purposes
  # option 1: data are summarized by trial, usually for statistics 
  # mutate(data = map(data, \(df) df %>% filter(epoch=="BSLN" | (epoch=="SPEECH" & time_t0_SPEECH<1)) %>%
  #                     group_by(trial_id, epoch, noise_type) %>%
  #                     summarize(value     =mean(value),
  #                               time_GTC  =mean(time_GTC)))) %>% # Filter data for epochs
  # option 2: data are time-resolved, ~10 ms resolution
  # mutate(data = map(data, ~filter(.x, epoch=="BSLN" | (epoch=="SPEECH" & time_t0_SPEECH<1)))) %>% # Filter data for epochs
  # option 3: data are time-resolved, ~10 ms resolution, no filtering for epoch, usually for plotting
  mutate(data = data) %>%
  # options end:---
  # filter(subject_id == "DM1002" & run_id == 3) %>%  # & channel == "ecog_L121") %>% # DEV only
  mutate(model_and_data = map(data, possibly(fit_lm_time, otherwise = 'error'))) %>%
  mutate(model_time = map(model_and_data, "model"),
         data = map(model_and_data, "data")) %>%
  # stat_noise_type = map_dbl(model_and_data, "stat_noise_type"), 
  # stat_timewin = map_dbl(model_and_data, "stat_timewin")) %>%
  select(-model_and_data) 

filtered_out <- 
  lmtime_results_channelwise_timeresolved %>%
  filter(map_lgl(data, is.null)) # there are a few channels that fail in the ecog data because all data are NAN
print(filtered_out)  

lmtime_results_channelwise_timeresolved <-
  lmtime_results_channelwise_timeresolved %>%
  anti_join(filtered_out)




# quality control -- verify that we have enough data for statistics
qc <-
  lmtime_results_channelwise %>%
  mutate(qc = map(data, \(df) df %>%
      group_by(interaction(epoch, noise_type)) %>%
      summarize(n_trial = n_distinct(trial_id)))) %>%
  select(-data, -model_time) %>%
  unnest(qc) %>%
  filter(n_trial < 10)
# 27 SU were removed 
qc %>% select() %>% distinct() 
  

# apply quality control to both lmtime tables 
lmtime_results_channelwise_timeresolved <- 
  lmtime_results_channelwise_timeresolved %>%
  anti_join(qc %>% select() %>% distinct(), by = vars_channel_id)

lmtime_results_channelwise <- 
  lmtime_results_channelwise %>%
  anti_join(qc %>% select() %>% distinct(), by = vars_channel_id)

# create summary for recordkeeping and reporting in manuscript
cohort_tabulation <-
  lmtime_results_channelwise %>%
  left_join(SUBJECTS_META %>% select(subject_id, target, diagnosis)) %>%
  {if (MODALITY=="SU") {separate_wider_delim(., channel, "_", names = c("micmac", "track", "unitidx"), cols_remove = FALSE) %>%
      mutate(channel_SUtrunc = glue("{micmac}_{track}"))} else . } %>%
  select(-data, -model_time) %>%
  ungroup() 

fname <- "recordkeeping-and-tabulation-after-qualitycontrol"
cohort_tabulation %>%
  group_by(subject_id) %>%
  summarize(n_run = n_distinct(run_id), 
            n_chan = n(), 
            target = first(target), 
            diagnosis = first(diagnosis)) %>%
  write_tsv(., file=glue('{PATH_FIG}/A04v5_00_{fname}_{MODALITY}.tsv'))

# sum %>%
#   summarize(n_chan = n(),
#             n_subj = n_distinct(subject_id),
#             avg_chanpersubj = n_chan / n_subj) 
# sum %>%
#   group_by(subject_id) %>%
#   summarize(n_chan = n()) %>%
#   ungroup() %>%
#   summarize(n_chan_mean = mean(n_chan),
#             n_chan_med = median(n_chan),
#             n_chan_max = max(n_chan),
#             n_chan_mmin = min(n_chan))
#   
# # for SU, we need to identify number of micro_l, micro_p, micro_c channels
# sum %>%
#   group_by(subject_id, channel_SUtrunc) %>%
#   summarize(n_chan = n()) %>%
#   ungroup() %>%
#   summarize(n_chanSUtrunc_mean = mean(n_chan), 
#             n_chanSUtrunc_med = median(n_chan),
#             n_chanSUtrunc_max = max(n_chan),
#             n_chanSUtrunc_min = min(n_chan))
  
  


# Visualize: Plot gamma over trials or within condition  ----

# plot eeg power across task epochs
df <- lmtime_results_channelwise$data[[3]]
df %>%
  # mutate(time_bin = floor(time_t0*100)) %>%
  # mutate(audioR_s = (audioR_s - mean(audioR_s)) / sd(audioR_s)) %>%
  # filter(epoch=="BSLN") %>%
  rename(eeg_res=value_residualtime, 
         eeg=value) %>%
  # mutate(time_t0=time_t0_SPEECH) %>%
  # filter(epoch %in% c("BSLN", "SPEECH")) %>%
  group_by(trial_id, epoch, noise_type) %>%
  summarize(eeg_res = mean(eeg_res), 
            eeg = mean(eeg), 
            time_GTC = mean(time_GTC)) %>%
  # mutate(eeg_diff = eeg_res - eeg)
  # select(trial_id, noise_type, eeg, eeg_res, time_t0, epoch) %>%
  # filter(trial_id < 10) %>%
  pivot_longer(cols = c("eeg_res"), names_to = "series") %>% #   pivot_longer(cols = c("audioR_s", "eeg")) %>%
  # group_by(time_bin, epoch) %>%
  # summarise(valuebin = mean(valuelog)) %>%
  # summarise(valuebin = mean(value_residualtime), intensity = mean(audioR_s)) %>%
  # mutate(intensity_logit = log(intensity / (1 - intensity))) %>%
  # pivot_longer(cols = c("valuebin", "intensity"), names) %>%
  ggplot(aes(x=epoch, y=value, color=factor(noise_type))) +
  # geom_line() +
  facet_grid(~epoch) %>%
  # stat_summary(geom = "ribbon", fun.data = "mean_cl_normal", aes(fill = factor(name)), alpha=0.3, linewidth=0) + 
  # stat_summary_bin(geom = "ribbon", fun.data = "mean_cl_normal", aes(fill = factor(noise_type)), alpha=0.3, linewidth=0, binwidth=0.100) +
  # stat_summary_bin(geom = "line", fun.data = "mean_cl_normal",   aes(color = factor(noise_type)), binwidth=0.100)
  geom_boxplot() +
  geom_jitter(aes(alpha=0.3)) +
  # geom_line()  +
  # geom_point()
  # geom_smooth(span = )
  # geom_bin2d(bins = 60)
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  # scale_fill_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  # ylim(-0.25, 0.25) + 
  # xlim(-1, 3)
# scale_alpha_manual(values = c("speech" = 0.9, "basel" = 0.4))+ 
# scale_shape_manual(values = c("speech" = 16, "basel" = 1)) +

df %>%
  filter(trial_id > 3 & trial_id <= 7) %>%
  group_by(trial_id, noise_type) %>%
  # mutate(value_norm = (value - value_BSLN)/value_BSLN) %>%
  # summarize(value = mean(value), value_BSLN = mean(value_BSLN), time = mean(time_GTC)) %>%
  ggplot(aes(x=time_GTC, y=time_t0_SPEECH)) + 
  stat_summary(geom = "point", fun.data = "mean_cl_normal", aes(color = factor(trial_id)), alpha=0.3, linewidth=0, binwidth=0.3) 
  # geom_histogram()
  # geom_point() + 
  scale_fill_gradient()
  geom_abline() + 
  coord_fixed()
  
summary(lmtime_results_channelwise$model_time[[1]])


# Validate SU data with visualizations ---- 

# 3 DM1015     micro_m_1  micro_m_unit1057... small negative deflection
# data[[21]] ... 1 DM1011     micro_c_1  micro_c_unit1343 ... channel has small fluctations but no clear activity
# data[43] ... 2 DM1022     micro_c_1 micro_c_unit422 ... channel has speech onset activity, starts ramping in the PREP window
# data  34 ....      2 DM1019     micro_m_1 <tibble [16,349 × 31]> micro_m_unit2099 <lm>      
# data  56 ....      3 DM1033     micro_l_2 <tibble [13,320 × 31]> micro_l_unit3217 <lm>      
chid = 'DM1022_2_micro_c_1'
df <- lmtime_results_channelwise_timeresolved$data[[43]]

# # pull from either 
# df <-  lmtime_results_channelwise_timeresolved %>% # eeg_nest_test for raw values or lmtime_results_channelwise
#   # DM1022-unit422, DM1014-unit867, DM1002-unit765, DM1019 run 2 micro_p_2, 
#   # DM1003 micro_l_4 unit1249, decreasing type at stimulus onset
#   filter(str_detect(channel_orig, 'unit867')) %>% 
#   pull(data) %>%
#   .[[1]] 

# Plot FR in the baseline over the course of the run
df %>%
  filter(epoch=="BSLN") %>% 
  group_by(trial_id, noise_type) %>%
  summarize(eeg = mean(value_residualtime), time_GTC = mean(time_GTC)) %>%

  ggplot(aes(x=time_GTC, y=eeg, color=noise_type)) + 
  geom_point() + 
  scale_color_manual(values = COLORS)

# visualizet the same but with more samples s
shades <- df %>%
  mutate(block_id = ceiling((trial_id)/10)) %>%
  group_by(noise_type, block_id) %>%
  summarise(start = min(time_GTC)-49290, end = max(time_GTC)-49290) 


df %>%
  mutate(timebin = cut(time_GTC, breaks = 150)) %>%
  group_by(timebin, noise_type) %>%
  summarize(eeg = mean(value_residualtime), 
            time_GTC = round(mean(time_GTC))-49290) %>%
  
  ggplot(aes(x=time_GTC, y=eeg, color=noise_type)) + 
  geom_point() + 
  scale_color_manual(values = COLORS) +   
  geom_rect(data = shades, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = noise_type), 
            alpha = 0.2, inherit.aes = FALSE) + 
  scale_fill_manual(values = COLORS)
ggsave(filename=glue('{PATH_FIG}/A04v5_z01_{chid}_eeg_over_time.pdf'), dpi=300, width=4, height=2, units='in')


# df <-
#   df %>% 
#   select(-value) %>%
#   rename(pwr=value_residualtime) %>% 
#   pivot_longer(cols = c("pwr", 'predicted'), values_to = c("value"), names_to = "series")


# chid = df$channel_id[[1]]
df %>%
  filter(epoch %in% c("BSLN", "SPEECH")) %>%
  
  ggplot(aes(x=trial_id, y=value, color=noise_type, shape=interaction(epoch, series))) + 
  geom_point() + 
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) + 
  scale_shape_manual(values=c(2, 17, 1, 16)) + 
  geom_rect(data = shades, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = noise_type), 
            alpha = 0.2, inherit.aes = FALSE) +
  scale_fill_manual(values = COLORS) +
  labs(title=glue('{chid}')) + 
  theme_bw() 




# Plot average FR time-locked to speech onset
# aa <-
  df %>%
  mutate(time_t0 = timewarp_t0) %>% # time_t0 = time_t0... or time_t0 = time_t0_SPEECH
  mutate(time = floor(time_t0*20)/20) %>%
  mutate(value = value_residualtime) %>%
  # mutate(noise_type = sample(noise_type)) %>% # split-half testing
  group_by(time, noise_type) %>%
  summarize(eeg_mean = mean(value), 
            eeg_sem = sd(value) / sqrt(n()), 
            n_dpoints = n(), 
            n_trls = n_distinct(trial_id)) %>%
  filter(time>-1 & time<5) %>%
  
  ggplot(aes(time)) +
  geom_ribbon(aes(ymin = eeg_mean-eeg_sem, ymax = eeg_mean+eeg_sem, alpha = 0.5, fill=noise_type)) +
  geom_line(aes(y=eeg_mean, color=noise_type)) + 
  geom_vline(xintercept = 0) 
  
  
  

  




# Plot long-term gamma trends over the course of the session   ---- 
aa <- 
    eeg_nest_test %>%
    mutate(data = map(data, ~ .x %>% mutate(time_since_session_start = time_GTC - min(time_GTC, na.rm = TRUE)))) %>%
    unnest(data) %>%
    mutate(
      time_bin = cut(time_since_session_start,
                     breaks = seq(0, ceiling(max(time_since_session_start, na.rm = TRUE)), by = 10),
                     include.lowest = TRUE,
                     right = FALSE,
                     labels = seq(10, ceiling(max(time_since_session_start, na.rm = TRUE)), by = 10))
    ) %>%
    group_by(subject_id, run_id, channel, time_bin) %>%
    summarize(
      value = mean(value, na.rm = TRUE),
      time = mean(time_since_session_start, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(time_bin) %>%
    summarize(
      value_mean = mean(value, na.rm = TRUE),
      value_sem = sd(value, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(time_bin = as.numeric(as.character(time_bin)))
  
    
aa %>%
  ggplot(aes(x = time_bin, y = value_mean)) +
    geom_ribbon(aes(ymin = value_mean - value_sem, ymax = value_mean + value_sem),
                fill = "lightblue", alpha = 0.4) +
    geom_line(color = "blue", size = 1) +
    labs(
      x = "Time (s)",
      y = "Signal (mean ± SEM)",
      title = "Average Neural Activity Across Sessions"
    ) +
    theme_minimal() 
    

aa %>%
  ggplot(aes(x = time_bin, y = n)) +
  # geom_ribbon(aes(ymin = value_mean - value_sem, ymax = value_mean + value_sem),
  #             fill = "lightblue", alpha = 0.4) +
  geom_line(color = "blue", size = 1) +
  labs(
    x = "Time (s)",
    y = "Signal (mean ± SEM)",
    title = "Average Neural Activity Across Sessions"
  ) +
  theme_minimal() 

    
# Case study of DM1015 ---- 


tmp1 <-
  lmtime_results_channelwise_timeresolved %>%
  select(-model_time) %>%
  filter(subject_id=="DM1015", run_id==2, channel=="macro_Lc") %>%
  unnest(data) %>%
  relocate(time_t0_SPEECH, acousticspectrum_intensity, timewarp_t0) %>%
  filter(time_t0_SPEECH>-0.2 & time_t0_SPEECH<2) # grab speech window only
  
corLag<-ccf(tmp1$value_residualtime, tmp1$time_t0_SPEECH, lag.max=50)

tibble(acf=drop(corLag$acf), lag=drop(corLag$lag)/100) %>%
  ggplot(aes(x=lag, y=acf)) + 
  geom_point() + 
  geom_hline(aes(yintercept=0))
  

lag_sec <- -0.1  # -100 ms
bin_width <- 0.1
lag_bins <- round(lag_sec / bin_width) 

tmp3 <- 
  tmp2 %>%
  filter(time_t0_SPEECH>-0.2 & time_t0_SPEECH<2) %>%
  group_by(timewarp_t0) %>%
  mutate(
    time_bin = cut(timewarp_t0,
                   breaks = seq(-1, 5.3, by = bin_width),
                   include.lowest = TRUE,
                   right = FALSE), 
                   # labels = seq(-1, 5.3, by = 0.050))
  ) %>%
  group_by(channel, run_id, trial_id, time_bin, loudness_lvl) %>%
  summarize(time_bin = mean(timewarp_t0),
            eeg = mean(value_residualtime), 
            intensity = mean(acousticspectrum_intensity)) %>%
  group_by(channel, run_id, trial_id) %>%
  arrange(time_bin) %>%
  mutate(eeg_shifted = dplyr::lead(eeg, n = abs(lag_bins)))
  # filter(time_bin > 2.8)
  

tmp3 %>%
  filter(time_bin > 2.8) %>%
  ggplot(aes(x=intensity, y=eeg_shifted, color=time_bin)) + 
  geom_point() + 
  stat_cor(method="spearman") + 
  facet_grid(run_id~channel) +
  coord_fixed(ratio = 1) + 
  geom_abline(a=1, b=1)

ggsave(filename=glue('{PATH_FIG}/A04v5_s1_DM1015-case-study_scatter_binwidth-{bin_width}_shift-{lag_bins}.pdf'), dpi=300, width=6, height=4, units='in', create.dir = TRUE)




# Fit conditions LM: here we take the residuals from the time LM ----
fit_lm_conditions <- function(df, is_quality_control=FALSE) {
  # df <- lmtime_results_channelwise$data[[39]] ℹ In group 27: `subject_id = "DM1015"`, `run_id = 3`, `channel = "micro_m_2"`.
  # browser()
  
  
  # we have to massage the df to ensure that trials don't have more than 2 samples
  # --this interferes with the blockwise shuffle algorithm
  # also ensure here that we have enough trials for statistics
  ngrp_trl <- df %>% group_by(trial_id) %>% summarize(ngrp_trl = n())
  df <- df %>% 
    left_join(ngrp_trl, by = join_by(trial_id)) # %>% 
    # filter(!(ngrp_trl==3 & row_number()==1 & epoch=="BSLN")) %>%
    # filter(!(ngrp_trl==1))
  ngrp_trlepoch <- df %>% group_by(trial_id, epoch) %>% summarize(ngrp_trlepoch = n())
  df <- df %>% 
    left_join(ngrp_trlepoch, by = join_by(trial_id, epoch)) %>% 
    filter(ngrp_trlepoch==1) %>%
    select(-ngrp_trl, -ngrp_trlepoch)
  # now quality control--do we have enough data?
  qc <- df %>% group_by(interaction(noise_type, epoch)) %>% summarize(n = n()) # quality control
  if (is_quality_control & (nrow(qc)<4 | any(qc$n <8))) {print(qc); return(NA)}

  
  
  # explicitly force no intercept with -1 
  # update 2024 12 29 this wasn't allowing BSLN to be used as the base factor, so I allow intercept
  model <- lm(value_residualtime ~ 1 + time_GTC + time_GTC2 + noise_type + epoch + noise_type:epoch, data = df)
  summary <- tidy(model) # this will return the beta coeffs
  aov_summary <- tidy(anova(model))
  
  # model <- lm(value_residualtime ~ 1 + noise_type + epoch + interaction(epoch, noise_type), data = df)
  # model <- aov(value_residualtime ~ 1 + noise_type + epoch + noise_type:epoch, data = df)
  # aov_summary <- tidy(model) # this will return the f-stat
    # rownames_to_column(aov_summary, var = "rowname")
  
  significant <- summary %>% filter(p.value < 0.05)
  df$residual_conditions <- residuals(model)
  # list(model = model, 
  #      data = df, 
  #      stat_noise_typeNOISE = summary$statistic[summary$term=="noise_typeNOISE"], # %>%filter(term=noise_type) %>% pull(statistic))
  #      stat_epochSPEECH = summary$statistic[summary$term=="epochSPEECH"],
  #      stat_epochSTIMULUS = summary$statistic[summary$term=="epochSTIMULUS"], 
  #      stat_noise_typeNOISExepochSPEECH = summary$statistic[summary$term=="noise_typeNOISE:epochSPEECH"])
  
  # anova/f-stat version
  list(stat_noise_typeNOISE = aov_summary$statistic[aov_summary$term=="noise_type"], # %>%filter(term=noise_type) %>% pull(statistic))
       coeff_noise_typeNOISE = summary$estimate[summary$term=="noise_typeNOISE"],
       stat_epochSPEECH = aov_summary$statistic[aov_summary$term=="epoch"],
       coeff_epochSPEECH = summary$estimate[summary$term=="epochSPEECH"],
       coeff_epochSTIMULUS = summary$estimate[summary$term=="epochSTIMULUS"],
       model = model,
       data = df,
       # stat_epochSTIMULUS = aov_summary$statistic[aov_summary$term=="epoch"], 
       stat_noise_typeNOISExepochSPEECH = aov_summary$statistic[aov_summary$term=="noise_type:epoch"], 
       coeff_noise_typeNOISExepochSPEECH = summary$estimate[summary$term=="noise_typeNOISE:epochSPEECH"])
  
  # # tidy/linear model/beta coeff version 
  # list(stat_noise_typeNOISE = summary$statistic[summary$term=="noise_typeNOISE"], # %>%filter(term=noise_type) %>% pull(statistic))
  #      stat_epochSPEECH = summary$statistic[summary$term=="epochSPEECH"],
  #      stat_epochSTIMULUS = summary$statistic[summary$term=="epochSTIMULUS"], 
  #      stat_noise_typeNOISExepochSPEECH = summary$statistic[summary$term=="noise_typeNOISE:epochSPEECH"])
}

# df <- lmtime_results_channelwise$data[[3]] %>% filter(length(m[[1]])==1)
lmcondition_results_channelwise <- lmtime_results_channelwise %>%
  mutate(
    m = map(data, \(d) fit_lm_conditions(d, is_quality_control=TRUE))  # Apply the function
  )
# testing: lmcondition_results_channelwise$model_conditions[[3]]

# these models failed--usually due to low trial count
lmcondition_null <- lmcondition_results_channelwise %>%
  filter(length(m[[1]])==1)  


lmcondition_results_channelwise <- 
  lmcondition_results_channelwise %>%
  anti_join(lmcondition_null) %>% # filter for only good units
  mutate(model_conditions = map(m, "model"), # OR this for visualizations
         data = map(m, "data")) %>%
  mutate(m = map(m, \(l) l[!(names(l) %in% c("model", "data"))])) %>% unnest_wider(m) # this for permutation testing
  

lmcondition_results_channelwise %>%
  arrange(desc(stat_epochSPEECH))

# Plot individual channels over trials, predicted vs actual from conditions LM ---- 
plot_single_channel_over_trials <- function(df, mdl, chid) {
  # df <- lmcondition_results_channelwise$data[[3]]
  # mdl <- lmcondition_results_channelwise$model_conditions[[3]]
  
  df$predicted <- predict(mdl, df)
  df$residual <- residuals(mdl)
  
  df <-
    df %>% 
    # select(-value) %>% # plot value or value_residualtime
    rename(pwr=value) %>% 
    pivot_longer(cols = c("pwr", 'predicted'), values_to = c("value"), names_to = "series")
  
  shades <- df %>%
    mutate(block_id = ceiling((trial_id)/10)) %>%
    group_by(noise_type, block_id) %>%
    summarise(start = min(trial_id), end = max(trial_id)) 
  
  # chid = df$channel_id[[1]]
  df %>%
  filter(epoch %in% c("BSLN"),  # c("BSLN", "SPEECH")
         series %in% c("pwr")) %>%
    
  ggplot(aes(x=trial_id, y=value, color=noise_type, shape=interaction(epoch, series))) + 
    geom_point() + 
    scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) + 
    scale_shape_manual(values=c(2, 17, 1, 16)) + 
    geom_rect(data = shades, 
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = noise_type), 
              alpha = 0.2, inherit.aes = FALSE) +
    scale_fill_manual(values = COLORS) +
    labs(title=glue('{chid}')) + 
    theme_bw() 
  # {PATH_FIG}/A04v5_04_channelwise-over-trials/{chid}_real-value-vs-predicted_BSLN-and-SPEECH.pdf'
  ggsave(filename=glue('{PATH_FIG}/A04v5_04_channelwise-over-trials/{chid}_real-value_BSLN.pdf'), dpi=300, width=5, height=3, units='in', create.dir = TRUE)
  
  df
}
# run the above function for each channel
lmcondition_results_channelwise %>%
  mutate(chid = glue('{subject_id}_run_{run_id}_{channel}')) %>%
  mutate(result = pmap(list(data, model_conditions, chid), plot_single_channel_over_trials))



# aa <- lmtime_results_channelwise %>% mutate(data = map(data, ~filter(.x, epoch %in% c("BSLN", "SPEECH"))))
# df <- aa$data[[1]]
# 
# aa <- lmtime_results_channelwise %>%
#   mutate(model_and_data = map(data, fit_lm_conditions)) 

# summary(lmcondition_results_channelwise$model_conditions[[1]])


# Permute and fit the condition linear models ----
permute_and_fit <- function(df, permute_type, predictor_of_interest) {
  # df <- lmcondition_results_channelwise$data[[27]]     group 27: `subject_id = "DM1015"`, `run_id = 3`, `channel = "micro_m_2"`.
  # predictor_of_interest <- "noise_type:epoch"

  # this is a circle-shift appproach
  shift_vector <- function(x, n) {
    c(tail(x, -n), head(x, n))
  }
  
  # Function to remove a predictor
  remove_predictor <- function(df, mdl_full, predictor_of_interest) {
    # Create a formula without the specified predictor
    # predictors <- names(all.vars(mdl_full))[-1]
    orig_formula <- formula(mdl_full)
    term_labels <- attr(terms(orig_formula), "term.labels")
    # terms <- all.vars(orig_formula)[-1]  # The first element is intercept, so we exclude it
    # term <- terms[2]
    
    # have to account for interaction terms... when we remove one variable, 
    # we also need to remove any interaction terms that invovle the variable
    predictors_to_remove <- term_labels[str_detect(term_labels, predictor_of_interest)]
    if (length(predictors_to_remove) > 1) {
      remove_string <- paste(predictors_to_remove[1], predictors_to_remove[2], sep = "-")
    } else {
      remove_string <- predictors_to_remove
    }
    new_formula <- update.formula(orig_formula, paste(". ~ . -", remove_string))
    
    # formula_without_predictor <- as.formula(paste("value ~", 
    #                                               paste(setdiff(, predictor_of_interest), collapse = " + ")))
    
    # Fit the new model
    new_model <- lm(new_formula, data = df)
  }
  

  # # alternatively, shuffle noise_type variable block-wise
  # shuffle_blockwise <- function(df) {
  #   df <- df %>% 
  #     mutate(block_id = floor((trial_id-1) / 10)) %>%
  #     ungroup()
  #   blocks_shuffled <- df %>%
  #     select(noise_type, block_id) %>% 
  #     distinct() %>% 
  #     mutate(block_id_shuffle = sample(block_id)) %>%
  #     select(-noise_type) 
  #   
  #   df <-
  #     df %>%
  #     group_by(block_id) %>%
  #     left_join(blocks_shuffled, by = join_by(block_id)) %>%
  #     # turn noise_type_shuffle into noise_type
  #     relocate(value_shuffle) %>% 
  #     select(-value) %>% 
  #     rename(value=value_shuffle) 
  #   df 
  # }
  
  shuffle_trialid_blockwise <- function(df) {
    # requires df$trial_id
    # adds column trial_id_shuff
    
    trls <- unique(df %>% pull(trial_id))
    # block <- ceiling(trls / 10) 
    # trl_within_block <- mod(trls, 10)
    
    tens <- rep(sample(seq(0,9)*10), each=10)
    ones <- rep((seq(0,9)), 10)
    df_shuff <-
      tibble(tens, ones) %>% 
      group_by(tens) %>% 
      mutate(ones = sample(ones), 
             trial_id_shuff = tens + ones) %>%
      filter(trial_id_shuff %in% trls) %>% 
      add_column(trial_id = trls) %>%
      ungroup() %>%
      select(-tens, -ones)
    
    df_shuff
  } # df_shuff <- shuffle_blockwise(df)
  
  
  # print(df)
  if (permute_type == "permute_trials") {
    df <- df %>% mutate(value = sample(value))
    
  } else if (permute_type == "permute_blocks") {
    df <- shuffle_blockwise(df) # assumes noise_type
    
  } else if (permute_type == "circshift_trials") {
    shift_n <- sample(1:nrow(df), 1)
    df <- df %>%
      mutate(value = shift_vector(value, shift_n))
    
  } else if (permute_type == "permute_residuals") {
    # create reduced/null models
    # mdl_full <- lm(value_residualtime ~ 1 + epoch + noise_type + noise_typeXepoch, data = df)
    mdl_full <- lm(value_residualtime ~ 1 + noise_type + epoch + noise_type:epoch, data = df)
    
    # r_squared_values <- sapply(predictors, function(p) remove_predictor(mdl_full, p))
    mdl_null <- remove_predictor(df, mdl_full, predictor_of_interest)
    
    null_residuals <- residuals(mdl_null)
    null_fitted <- fitted.values(mdl_null)
    
    random_null_residuals <- sample(null_residuals, length(null_residuals), replace = FALSE)
    df$value_residualtime <- null_fitted + random_null_residuals
    
  } else if (permute_type == "permute_residuals_blockwise") {
    # create reduced/null models
    # mdl_full <- lm(value_residualtime ~ 1 + epoch + noise_type + noise_typeXepoch, data = df)
    # value_residualtime ~ 1 + value_BSLN + annot_time_GTC + annot_time_GTC2 + noise_type + epoch + noise_type:epoch
    mdl_full <- lm(value_residualtime ~ 1 + time_GTC + time_GTC2 + noise_type + epoch + noise_type:epoch, data = df)
    
    # mdl_full <- lm(value_residualtime ~ 1 + noise_type + epoch + noise_type:epoch, data = df)
    
    # r_squared_values <- sapply(predictors, function(p) remove_predictor(mdl_full, p))
    mdl_null <- remove_predictor(df, mdl_full, predictor_of_interest)
    
    df$null_residuals <- residuals(mdl_null)
    df$null_fitted <- fitted.values(mdl_null)
    
    if (predictor_of_interest == "noise_type") {
      
      # shuffle trial_id blockwise
      df_trialid_shuff <- 
        df %>%
        group_by(trial_id, epoch) %>%
        select(trial_id, epoch, null_residuals) %>% # needs to be updated if there is also a time component
        left_join(shuffle_trialid_blockwise(df %>% ungroup() %>% select(trial_id)), 
                  join_by(trial_id)) %>%
        group_by(trial_id_shuff, epoch) %>%
        select(-trial_id) %>%
        rename(random_null_residuals=null_residuals, 
               trial_id=trial_id_shuff) # rename variables for join_by
      # validation
      # ngrp <- df %>% group_by(trial_id, epoch) %>% summarize(ngrp = n())
      # ngrp <- df_trialid_shuff %>% group_by(trial_id, epoch) %>% summarize(ngrp = n())
      
      
      # left-join df with the shuffled data
      df <-
        df %>%
        group_by(trial_id, epoch) %>%
        left_join(df_trialid_shuff, 
                  join_by(trial_id, epoch))
      
      # df$random_null_residuals <- df_trialid_shuff$null_residuals
      
    } else if (predictor_of_interest == "epoch") {
      # random_null_residuals <- sample(null_residuals, length(null_residuals), replace = FALSE)
      df$random_null_residuals <- sample(df$null_residuals, replace = FALSE)
      
    } else {
      stop("Invalid predictor variable")
    }
    
    df$value_residualtime <- df$null_fitted + df$random_null_residuals
    
    
  } else {
    stop("Invalid shuffle_type")
  }
  
  # is_quality_control FALSE means we don't want to filter out data with low trial counts
  # we already did this when we weren't permutation testing
  fit_lm_conditions(df, is_quality_control=FALSE) 
}


get_percentile <- function(distr, val) {
  fn<-ecdf(distr)
  fn(val)
}

# Permutation test
set.seed(123)
num_permutations <- 500
permute_type <- "permute_residuals" # permute_trials, permute_blocks, circshift_trials, permute_residuals, permute_residuals_blockwise
# predictors <- c("epochSPEECH", "noise_typeNOISE", "noise_typeNOISE:epochSPEECH")
# predictors <- attr(terms(formula(lmcondition_results_channelwise$model_conditions[[1]])), "term.labels")
var_strs <-           c("noise_type",      "epoch", "noise_type:epoch") #      "noise_type:epoch") # predictor of interest... "stat_noise_type" or "stat_timewin" or "stat_intensity_praat"
var_strs_coeffname <- c("noise_typeNOISE", "epochSPEECH", "noise_typeNOISExepochSPEECH") # "noise_typeNOISExepochSPEECH")
for (i in seq_along(var_strs)) {
  # Access the current var_str and corresponding coeffname by index
  var_str <- var_strs[i]
  var_str_coeffname <- var_strs_coeffname[i]
  
  # debugging: aa<-permuted_lmcondition_results_channelwise$stat_permute_epochSPEECH[[19]]
  # df <- lmcondition_results_channelwise$data[[3]]
  permuted_lmcondition_results_channelwise <- 
    replicate(num_permutations, {
      permuted_results <- 
        lmcondition_results_channelwise %>%
        # ungroup() %>% filter(row_number() %in% c(50, 123, 189)) %>% group_by(subject_id, run_id, channel) %>% # TESTING
        mutate(model_and_data = map(data, \(x) permute_and_fit(x, permute_type, var_str))) %>%
        mutate(!!glue("stat_{var_str_coeffname}") := map_dbl(model_and_data, glue("stat_{var_str_coeffname}"))) %>% 
        select(!!glue("stat_{var_str_coeffname}"))
      # sum(permuted_results)
    }, simplify = FALSE)

  permuted_lmcondition_results_channelwise <- 
    bind_rows(permuted_lmcondition_results_channelwise) %>%
    nest() %>%
    rename(!!glue("stat_permute_{var_str_coeffname}") := data)
  
  lmcondition_results_channelwise <-
    lmcondition_results_channelwise %>%
    select(-any_of(!!(glue("stat_permute_{var_str_coeffname}")))) %>% # make sure the left join gets the updated values
    left_join(permuted_lmcondition_results_channelwise) %>%
    # rowwise() %>% 
    # mutate(stat_noise_type_pval = mean(abs(data_permute$stat_noise_type) >= abs(stat_noise_type))) %>%
    # mutate(stat_timewin_pval = mean(abs(data_permute$stat_timewin) >= abs(stat_timewin))) 
    mutate(!!glue("stat_{var_str_coeffname}_prctl") := map2_dbl(!!sym(glue("stat_permute_{var_str_coeffname}")), 
                                                                !!sym(glue("stat_{var_str_coeffname}")), 
                                                                \(x, y) sum((y) > (x[[1]])) / length(x[[1]])), # \(x, y) get_percentile(x[[1]], y)), 
                                                                # \(x, y) sum(abs(y) > abs(x[[1]])) / length(x[[1]])), # \(x, y) get_percentile(x[[1]], y)), 
           !!glue("sig_{var_str_coeffname}") := !!sym(glue("stat_{var_str_coeffname}_prctl"))>0.95) 
          
  # tail(lmcondition_results_channelwise %>% select(c(1:4, 13:18)), n = 10)
  # !!sym(glue("stat_permute_{var_str}"))
  # !!sym(glue("stat_{var_str_coeffname}"))
}

# testing for robustness... want to see if two different runs with the same params returns the same results
# permute_type <- glue('{permute_type}-iter02') 

lmcondition_results_channelwise <-
  lmcondition_results_channelwise %>%
  mutate(sig = case_when(sig_noise_typeNOISE & sig_epochSPEECH ~ 'both', 
                     sig_noise_typeNOISE ~ 'noise_type', 
                     sig_epochSPEECH ~ 'epoch',
                     TRUE~'none')) 

fname <- 'channelwise_lm_timeres-trial-speech1s_stat-fstat_mdl-time-joint-fit'
write_tsv(lmcondition_results_channelwise, file=glue('{PATH_FIG}/A04v5_03_{fname}_permute-{permute_type}_channelwise_{format(Sys.time(), "%Y%m%d%H%M")}.tsv'))
# reload an old permutation test
# for LFP
# lmcondition_results_channelwise <- read_tsv(file=glue('{PATH_FIG}/A04v3_03_{fname}_permute-{permute_type}_channelwise.tsv'))
#  lmcondition_results_channelwise <- read_tsv(glue('{PATH_FIG}/A04v5_03_channelwise_lm_timeres-trial-speech1s_stat-fstat_mdl-time-joint-fit_permute-permute_residuals_channelwise_202504141618.tsv'))
# for SU
# lmcondition_results_channelwise <- read_tsv(glue('{PATH_FIG}/A04v5_03_channelwise_lm_timeres-trial-speech1s_stat-fstat_mdl-time-joint-fit_permute-permute_residuals_blockwise_channelwise_202503102040.tsv'))
# lmcondition_results_channelwise <- read_tsv(glue('{PATH_FIG}/A04v5_03_channelwise_lm_timeres-trial-speech1s_stat-fstat_mdl-time-joint-fit_permute-permute_residuals_channelwise_202504141520.tsv'))



# Stand-in/dummy/fake rather than running actual permutation test... ----
# summary(lmcondition_results_channelwise$m[[1]])
permute_type <- "dummy"
lmcondition_results_channelwise <-
  lmcondition_results_channelwise %>% 
  ungroup() %>%
  arrange(desc(stat_noise_typeNOISExepochSPEECH)) %>%
  # mutate(sig_noise_typeNOISExepochSPEECH = if_else(row_number()<30, TRUE, FALSE)) %>%
  mutate(sig_noise_typeNOISExepochSPEECH = stat_noise_typeNOISExepochSPEECH>5)  %>%
  
    
  arrange(desc(stat_noise_typeNOISE)) %>%
  # mutate(sig_noise_typeNOISE = if_else(row_number()<30, TRUE, FALSE)) %>%
  mutate(sig_noise_typeNOISE = stat_noise_typeNOISE>5)  %>%
  
  arrange(desc(stat_epochSPEECH)) %>%
  # mutate(sig_epochSPEECH = if_else(row_number()<30, TRUE, FALSE)) %>%
  mutate(sig_epochSPEECH = stat_epochSPEECH>5)  %>%
  
  group_by(subject_id, run_id, channel) %>%
  
  mutate(sig = case_when(sig_noise_typeNOISE & sig_epochSPEECH ~ 'both', 
                         sig_noise_typeNOISE ~ 'noise_type', 
                         sig_epochSPEECH ~ 'epoch',
                         TRUE~'none')) 
# %>% select(tail(names(.), 5))




# Plot channelwise permutation null distributions and actual values ----
n_obs <- nrow(lmcondition_results_channelwise) 
plot_top_n <- 20 # number of channels to plot,
var_str_coeffname <- "epochSPEECH" # noise_typeNOISE", "epochSPEECH", "noise_typeNOISExepochSPEECH
lmcondition_results_channelwise_subset <- 
  lmcondition_results_channelwise %>% 
  ungroup() %>% 
  arrange(desc(abs(across(glue("stat_{var_str_coeffname}"))))) %>% 
  slice_head(n=plot_top_n) %>%
  mutate(channel_id = paste(subject_id, channel, run_id, sep = "_")) %>%
  # select(-all_of(vars_channel_id)) %>%
  ungroup() %>%
  mutate(facet_order = fct_inorder(channel_id)) %>%
  group_by(facet_order)

ggplot() + 
  geom_histogram(data = lmcondition_results_channelwise_subset %>% 
                        select(all_of(c(glue("stat_permute_{var_str_coeffname}")))) %>% 
                        unnest(), 
                 aes(x = !!sym(glue("stat_{var_str_coeffname}")), y = ..density..), 
                 fill = "grey") + 
  # geom+(data = anova_results_channelwise_subset, 
  #               aes(y = 0, xmin = stat_noise_type_q05, xmax = stat_noise_type_q95), height = 0.1, color = "red") +
  geom_vline(data = lmcondition_results_channelwise_subset, aes(xintercept=!!sym(glue("stat_{var_str_coeffname}"))) , color="red") + 
  facet_wrap(facet_order~., scales = "free") + 
  theme_minimal() +
  theme(axis.text.y = element_blank())
# coord_cartesian(xlim = c(-40, 40), ylim = c(0, 0.2))
fname <- "null-dist-vs-actual"
ggsave(file=glue('{PATH_FIG}/A04v5_03_{fname}_permute-{permute_type}_var-{var_str_coeffname}_channelwise-noecog.pdf'), dpi=300, width=10, height=8,  units='in') # 



# print summary table
lmcondition_results_channelwise %>% 
  ungroup() %>%
  group_by(sig) %>%
  summarise(n = n())
  # summarise(sig_noise_typeNOISE = sum(sig_noise_typeNOISE), 
  #           sig_epochSPEECH = sum(sig_epochSPEECH),
  #           sig_noise_typeNOISExepochSPEECH = sum(sig_noise_typeNOISExepochSPEECH), 
  #           n = n())

axis_limits <- if (MODALITY == "LFP" & !CTX_FLAG) { list(xlim(0, 10), ylim(0, 50)) 
        } else if (MODALITY == "LFP" & CTX_FLAG) {  list(xlim(0, 100), ylim(0, 100)) 
        } else if (MODALITY == "SU" )             { list(xlim(0, 20), ylim(0, 100)) }
                   

# plot scatter of f-stats
lmcondition_results_channelwise <-
  lmcondition_results_channelwise %>%
  separate(channel, into = c('type', 'track', 'unit_idx'), remove = FALSE) %>%
  mutate(modality = MODALITY,
         channel_join = if_else(modality=="LFP", channel, glue('{type}_L{track}'))) %>%
  # select(subject_id, run_id, channel_join) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, target, region_clean),
            by = join_by(subject_id, run_id==lombard_run_id, channel_join==name))

lmcondition_results_channelwise %>%
  mutate(
    sig = case_when(
      sig == "noise_type" ~ "NOISE",
      sig == "epoch" ~ "SPEECH",
      sig == "both" ~ "NOISE & SPEECH",
      .default = sig
    )
  ) %>%
  ggplot(aes(
    x = stat_noise_typeNOISE,
    y = stat_epochSPEECH,
    shape = target,
    fill = sig
  )) +
  geom_point(alpha = 0.7, size = 4, color = "black", stroke = 0.5) +
  scale_fill_manual(
    values = c(
      "NOISE & SPEECH" = "blue",
      "NOISE" = COLORS[["NOISE"]],
      "none" = "lightgrey",
      "SPEECH" = "black"
    ),
    breaks = c("NOISE & SPEECH", "NOISE", "none", "SPEECH")
  ) +
  scale_shape_manual(
    values = c("STN" = 21, "GPi" = 24)
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, color = "black")),
    shape = guide_legend(override.aes = list(fill = "grey"))
  ) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  # {if(MODALITY=="LFP") xlim(0, 15) +  ylim(0, 20)} +
  # else if (MODALITY=="LFP") {xlim(0, 15) + ylim(0, 20)}} +
  # xlim(0, 15) + ylim(0, 20) +   # for BGA fstat
  # ylim(0, 100) for SU-FR, ylim(0, 20) for BGA
  # xlim(-10, 10) + ylim(-22, 30) + # for SU-FR betacoeff
  # xlim(-.1, 10) + ylim(-.1, 30) + # for LFP-FR betacoeff
  labs(
    title = glue("All channels: permutation ({permute_type}) test percentile"),
    x = "stat_noise_typeNOISE",
    y = "stat_epochSPEECH",
    shape = "Target",
    fill = "Significance"
  ) +
  # axis_limits +
  theme_bw()

fname <- 'epoch-vs-noise-type-scatter'
ggsave(file=glue('{PATH_FIG}/A04v5_03-02v2_{fname}_permute-{permute_type}_stat-fstat.pdf'), dpi=300, width=6, height=4, units='in')



# fisher exact test to determine if the frequencies are significant
lmcondition_results_channelwise %>%
mutate(
  sig_noise_typeNOISE = recode_factor(as.character(sig_noise_typeNOISE), `FALSE` = "non-NOISE", `TRUE` = "NOISE"),
  sig_noise_typeNOISE = fct_relevel(sig_noise_typeNOISE, "non-NOISE", "NOISE"),  # Set level order
  sig_epochSPEECH = recode_factor(as.character(sig_epochSPEECH), `FALSE` = "non-SPEECH", `TRUE` = "SPEECH"),
  sig_epochSPEECH = fct_relevel(sig_epochSPEECH, "non-SPEECH", "SPEECH")  # Set level order
)

counts <-
  lmcondition_results_channelwise %>%
  # mutate(sig_epochSPEECH = recode_factor("non-SPEECH", "SPEECH",), 
  #        sig_noise_typeNOISE = if_else(!sig_noise_typeNOISE, "non-NOISE", "NOISE")) 
  mutate(
    sig_noise_typeNOISE = recode_factor(as.character(sig_noise_typeNOISE), `FALSE` = "non-NOISE", `TRUE` = "NOISE"),
    sig_noise_typeNOISE = fct_relevel(sig_noise_typeNOISE, "non-NOISE", "NOISE"),  # Set level order
    sig_epochSPEECH = recode_factor(as.character(sig_epochSPEECH), `FALSE` = "non-SPEECH", `TRUE` = "SPEECH"),
    sig_epochSPEECH = fct_relevel(sig_epochSPEECH, "SPEECH", "non-SPEECH")  # Set level order
  ) %>%
  # mutate(i = interaction(sig_noise_typeNOISE, sig_epochSPEECH )) %>%
  group_by(sig_noise_typeNOISE, sig_epochSPEECH) %>%
  arrange(.by_group = T) %>%
  summarise(n = n())  # %>%
  # pivot_wider(names_from = sig_noise_typeNOISE, values_from = n) %>% 
  # column_to_rownames("sig_epochSPEECH")

# # evidently this converts the binary vectors to a contingency table
# contingency_table <- table(lmcondition_results_channelwise$sig_noise_typeNOISE,
#                            lmcondition_results_channelwise$sig_epochSPEECH)


# install.packages('vcd')
# install.packages('Barnard')
library("vcd")
library("Barnard")

xt <- xtabs(n ~ sig_epochSPEECH + sig_noise_typeNOISE, counts)
mosaic(xt, shade=T, colorize = T, , labeling = labeling_values,
       gp = gpar(fill=matrix(c("black", "lightgrey", "blue", COLORS[["NOISE"]]), 2, 2), 
                 alpha=0.7))
fname <- glue('mosaic-channeltypes') # statistically significant only 
dev.print(pdf, glue('{PATH_FIG}/A04v5_03-03_{fname}_permute-{permute_type}_stat-fstat.pdf'), width=2.5, height=2.5)
# ggsave(file=glue('{PATH_FIG}/A04v5_03-03_{fname}.pdf'), dpi=300, width=3, height=3, units='in') # height=25 for multi-channel plots
counts$fisherp <- fisher.test(xt)$p.value
# counts$barnardp <- barnard.test(44, 6, 37, 3, dp = 0.001, pooled = TRUE)$p.value[[1]] # LFP - GP&STN
# counts$barnardp <- barnard.test(15, 23, 11, 6, dp = 0.001, pooled = TRUE)$p.value[[1]] # SU
# counts$barnardp <- barnard.test(8, 2, 26, 4, dp = 0.001, pooled = TRUE)$p.value[[1]] # LFP - VIM

write_tsv(counts, file=glue('{PATH_FIG}/A04v5_03-03_{fname}-counts_permute-{permute_type}_stat-fstat.tsv'))







# Test whether speech channels are spatially aggregated ---- 
# install.packages("rdist"); library(rdist)
channel_subset_name <- "sig_epochSPEECH-STN" # sig_epochSPEECH-INC, sig_noise_typeNOISE-INC, active-INC

# prepare plotting table
data_all <-
  lmcondition_results_channelwise %>%
  separate(channel, into = c('type', 'track', 'unit_idx'), remove = FALSE) %>%
  mutate(modality = MODALITY,
         channel_join = if_else(modality=="LFP", channel, glue('{type}_L{track}'))) %>%
  # select(data, channel_join) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, target, region_clean, mni_x, mni_y, mni_z),
            by = join_by(subject_id, run_id==lombard_run_id, channel_join==name)) %>%
  # left_join(SUBJECTS_META %>% select(subject_id, target),
  #           by = join_by(subject_id)) %>%
  # mutate(region_clean = target) %>%
  mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) %>%
  mutate(channel_id = factor(channel_id, levels=channel_levels$channel_id)) %>%
  
  # lmtime_results_channelwise_slttest contains info from t-test INC or DEC type channel
  left_join(lmtime_results_channelwise_slttest %>% group_by(across(all_of(vars_channel_id))) %>% select(-data, -model_time)) %>% 
  ungroup() %>%
  filter(!is.na(mni_x))
  
data_all <-
  data_all %>%
  {if      (str_detect(channel_subset_name, "-INC")) {filter(., resp_type=="INC")} 
    else if (str_detect(channel_subset_name, "-DEC")) {filter(., resp_type=="DEC")}
    else . } %>%
  
  {if      (str_detect(channel_subset_name, "STN")) {filter(., target=="STN")} 
    else if (str_detect(channel_subset_name, "GPi")) {filter(., target=="GPi")}
    else . }


data_subset <- 
  data_all %>%
  # filter(row_number() %in% sample.int(100, 15)) %>%
  {if (str_detect(channel_subset_name, "active")) {filter(., sig!="none")} 
    else if (str_detect(channel_subset_name, "sig_epochSPEECH")) {filter(., sig_epochSPEECH)}
    else if (str_detect(channel_subset_name, "sig_noise_typeNOISE")) {filter(., sig_noise_typeNOISE)}
    else if (channel_subset_name=="all") {.}
    else { stop(glue('channel_subset_name: "{channel_subset_name}" not recognized'))} } # get significant channels only
  



# Define the bootstrap function to compute mean 
library(rdist)
boot_mean <- function(df, indices) {
  sample_data <- df[indices, ]  # Resample the data
  pos <- as.matrix(sample_data %>% select(mni_x, mni_y, mni_z))
  mean_A <- mean(rdist(pos), na.rm = TRUE)
  return(mean_A )  # Compute the difference in means
}
# # Perform bootstrap resampling
boot_result_all <- boot(data_all, statistic = boot_mean, R = 1000)
boot_result_subset <- boot(data_subset, statistic = boot_mean, R = 1000) # turn this into a single value
subset_dist <- boot_mean(data_subset, seq_len(nrow(data_subset))) # actual distance amongst speech-active 
percentile <- mean(as.double(boot_result_all$t <= subset_dist)) 


tibble(group = "all", dist = boot_result_all$t) %>%
  # rbind(. , tibble(group = "STN", dist = boot_result_subset$t) ) %>% # ADD or REMOVE the 'actual' distribution
  ggplot(aes(x=dist, fill=group)) +
  geom_histogram(aes(y=after_stat(density))) +
  # geom_density(alpha=0.5) +
  geom_vline(xintercept = subset_dist, color='red', size=2) +
  annotate(
    "text",
    x = subset_dist,
    y = Inf,            # set y to 90% of the histogram height
    label = glue("p={round(percentile, 3)}"),
    angle = 90,                 # rotate text vertically
    color = "black",
    hjust = 1,                  # left-align
    vjust = 2,
    size = 4
  ) +
  scale_fill_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) +
  scale_color_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) +
  xlab('Distance [mm]')
fname <- glue('spatial-aggregation_permutation-{permute_type}_susbet-{channel_subset_name}_onehist')
# fname <- glue('spatial-aggregation_permutation-{permute_type}_susbet-{channel_subset_name}_twohist') 
# fname <- glue('spatial-aggregation_permutation-{permute_type}_susbet-{channel_subset_name}_twodensity') 
ggsave(file=glue('{PATH_FIG}/A04v5_12_{fname}.pdf'), dpi=300, width=2.5, height=1.5, units='in') # height=25 for multi-channel plots


data_all %>% 
  mutate(group='all') %>%
  rbind(data_subset %>% 
        mutate(group='subset')) %>%
  group_by(group) %>%
  select(mni_x, mni_y, mni_z) %>%
  summarize_all_stats() %>%
  # pivot_longer(cols = everything()) %>%
  write_tsv(., file=glue('{PATH_FIG}/A04v5_12_{fname}_data_all_stats.tsv'))



# v0 Fit sliding window t-tests for each channel ---- 
# df <- lmtime_results_channelwise$data[[27]] # 93=DM1033 macro_Lc
winsize <- 0.250 # window size of data for the t-test
stride <- 0.1

channel_subset_name <- "sig_epochSPEECH"
num_permutations <- 500
alpha = 0.05

# ----- option 1
# tlock <- 'time_t0_SPEECH'
# t_centers <- (seq(-2, 2, by = stride))
# ----- option 2
tlock <- 'timewarp_t0'
t_centers <- (seq(2, 5, by = stride))


apply_window_ttest_v0 <- function(t_center, df) {
  # commented in 20250422  ↓
  df2 <- 
    df %>%
    mutate(epoch = if_else(epoch!="BSLN", "SPEECH", epoch)) %>%
    # filter((timewarp_t0>-0.75 & timewarp_t0<=0) | (timewarp_t0>(t_center-winsize) & timewarp_t0<=t_center)) %>%
    {if (tlock=="time_t0_SPEECH")  { filter(., (epoch=="BSLN" & time_t0>-0.75 & time_t0<=0) | (time_t0_SPEECH>(t_center-winsize) & time_t0_SPEECH<=t_center))
    } else if (tlock=="timewarp_t0") { filter(., (timewarp_t0>-0.75 & timewarp_t0<=0) | (timewarp_t0>(t_center-winsize) & timewarp_t0<=t_center))
    } else assert(0, "tlock not recognized.")} %>%
    group_by(epoch) # %>% summarize(n = n(), max = max(timewarp_t0), min = min(timewarp_t0)) 
  # commented out 20250422 ---- ↓
  # df2 <-
  #   df %>%
  #   mutate(epoch = if_else(epoch!="BSLN", "SPEECH", epoch)) %>%
  #   filter((timewarp_t0>-0.75 & timewarp_t0<0) |
  #            (timewarp_t0>(t_center-winsize) & timewarp_t0<t_center)) %>%
  #   group_by(epoch) # %>% summarize(n = n(), max = max(timewarp_t0), min = min(timewarp_t0))
  # commented out 20250422 ----
  
  # inspect df with
  # df2 %>%
  #   group_by(trial_id, epoch) %>%
  #   summarize(value_residualtime = mean(value_residualtime, na.rm=T), 
  #             time_t0_SPEECH = mean(time_t0_SPEECH, na.rm=T)) %>%
  #   ggplot(aes(x = time_t0_SPEECH, y= value_residualtime, color=epoch, group=trial_id)) +
  #   geom_jitter() + 
  #   geom_line()
  
  td <- tidy(lm(value_residualtime ~ epoch, data = df2))
  tibble(p = td$p.value[[2]], beta = td$estimate[[2]])
  # td$p.value[[2]] * sign(td$estimate[[2]])
  # tibble(td$estimate[[2]]
}

# t_center <- results$t_centers[[93]]
# apply_window_ttest(t_center, df)

find_first_sig_v0 <- function(df) {
  # browser()
  results <-
    tibble(t_centers) %>%
    mutate(d = map(t_centers, \(x) apply_window_ttest_v0(x, df))) %>%
    unnest_wider(d) %>%
    mutate(sig = p<(0.05)) %>%
    filter(sig & beta>0) %>%
    {
      if (tlock=="timewarp_t0")   {        filter(., t_centers>2.5) # 2.5 in previous versions
      } else if (tlock=="time_t0_SPEECH") {  filter(., t_centers>-2)
      } else assert(0, "tlock not recognized.")
    } %>%
    mutate(ismaxamplitude = (beta==max(beta)))
  
  # r2 %>%
  #   # mutate(p = log(p)) %>%
  #   ggplot(aes(t_centers, beta, color=sig)) +
  #   geom_point() +
  #   geom_hline(yintercept = 0)
  
  if (nrow(results)==0) {
    earliest_time <- NA
    peak_time <- NA
    peak_amplitude <- NA
  } else {
    earliest_time <- min(results %>% pull(t_centers))
    peak_time <- results %>% filter(ismaxamplitude) %>% pull(t_centers)
    peak_amplitude <- results %>% filter(ismaxamplitude) %>% pull(beta)
  }
  
  tibble(earliest_time, peak_time, peak_amplitude)
}


lmtime_results_channelwise2 <- 
  lmtime_results_channelwise_timeresolved %>%
  # filter(row_number() %in% sample.int(100, 4)) %>%
  mutate(results = map(data, find_first_sig_v0)) %>%
  unnest_wider(results)

lmtime_results_channelwise3 <- 
  lmtime_results_channelwise2 %>% 
  select(-data, -model_time) %>%
  mutate(earliest_time = earliest_time, 
         peak_time = peak_time) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, region_clean), 
            by = join_by(subject_id, run_id==lombard_run_id, channel==name)) %>%
  mutate(region_clean_label = case_when(str_detect(region_clean, "GP") ~ "GP",
                                        str_detect(region_clean, "STN") ~ "STN",
                                        .default = NA)) %>%
  semi_join(lmcondition_results_channelwise %>%
              filter(!!sym(channel_subset_name)),
            by = join_by(subject_id, run_id, channel)) # get significant channels only

# determine what our sample size should be 
lmtime_results_channelwise3 %>%
  group_by(region_clean_label) %>%
  summarize(n = n())


#  added 20250421  plot boxplot+jitter of points 
var_of_interest <- "earliest_time" # earliest_time, peak_time, peak_amplitude, 
lmtime_results_channelwise3 %>%
  filter(!is.na(region_clean_label)) %>%
  ggplot(aes(x=region_clean_label, y=!!sym(var_of_interest), color=region_clean_label)) + 
  geom_boxplot(outliers=F) +
  geom_jitter(alpha=0.8, width = 0.1) + 
  stat_compare_means() +
  scale_color_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) + 
  scale_y_continuous(expand = expansion(mult=c(0.1, 0.1)))
# tmp break ---
fname <- glue('boxplot_tlock-{tlock}_sig-{var_of_interest}-ttest-{channel_subset_name}_STN-and-GP_winsize-{winsize}-casual') # statistically significant only 
ggsave(file=glue('{PATH_FIG}/A04v5_10b_v0_{fname}.pdf'), dpi=300, width=3.5, height=2, units='in') # width=4, height=2 for 



# v1 Fit sliding window t-tests for each channel ---- 
# df <- lmtime_results_channelwise$data[[27]] # 93=DM1033 macro_Lc
winsize <- 0.250 # window size of data for the t-test
stride <- 0.1

# ----- option 1
tlock <- 'time_t0_SPEECH'
t_centers <- (seq(-2, 2, by = stride))
# # ----- option 2
# tlock <- 'timewarp_t0'
# t_centers <- (seq(2, 5, by = stride))

num_permutations <- 500
alpha = 0.01
min_cluster_size <- 5

apply_window_ttest <- function(t_center, df) {
  df2 <- 
    df %>%
    mutate(epoch = if_else(epoch!="BSLN", "SPEECH", epoch)) %>%
    # filter((timewarp_t0>-0.75 & timewarp_t0<=0) | (timewarp_t0>(t_center-winsize) & timewarp_t0<=t_center)) %>%
    {if (tlock=="time_t0_SPEECH")  { filter(., (epoch=="BSLN" & time_t0>-0.75 & time_t0<=0) | (time_t0_SPEECH>(t_center-winsize) & time_t0_SPEECH<=t_center))
    } else if (tlock=="timewarp_t0") { filter(., (timewarp_t0>-0.75 & timewarp_t0<=0) | (timewarp_t0>(t_center-winsize) & timewarp_t0<=t_center))
    } else assert(0, "tlock not recognized.")} %>%
    group_by(epoch) # %>% summarize(n = n(), max = max(timewarp_t0), min = min(timewarp_t0)) 
  
  # testing/debugging
  # df2 %>% group_by(trial_id, epoch) %>% summarize(n = n())
  # 
  # df2 %>%
  #   ggplot(aes(x = value_residualtime, color=epoch)) +
  #   geom_histogram()
  
  
  # # Define the bootstrap function to compute mean difference
  # boot_mean_diff <- function(df, indices) {
  #   sample_data <- df[indices, ]  # Resample the data
  #   mean_A <- mean(sample_data$value_residualtime[sample_data$epoch == "SPEECH"], na.rm = TRUE)
  #   mean_B <- mean(sample_data$value_residualtime[sample_data$epoch == "BSLN"], na.rm = TRUE)
  #   return(mean_A - mean_B)  # Compute the difference in means
  # }
  # # Perform bootstrap resampling
  # boot_result <- boot(df2, statistic = boot_mean_diff, R = 500)
  # # # Calculate 95% confidence interval
  # boot_ci <- boot.ci(boot_result, type = "perc")  # Percentile-based CI
  # p_value <- mean(sign(boot_result$t0)*(boot_result$t) <= 0)
  # beta <- boot_result$t0
  # lbound <- boot_ci$percent[4]
  # ubound <- boot_ci$percent[5]
  # tibble(p = p_value, beta, lbound, ubound)
  
  # t-tests
  td <- tidy(lm(value_residualtime ~ epoch, data = df2))
  tibble(p = td$p.value[[2]], beta = td$estimate[[2]])
  
  # td$p.value[[2]] * sign(td$estimate[[2]])
  # tibble(td$estimate[[2]]
  
}

# t_center <- results$t_centers[[93]]
# apply_window_ttest(t_center, df)

find_first_sig <- function(df) {
  
  print("Running sliding window stat tests.")
  results_all <- 
    tibble(t_centers) %>%
    mutate(d = map(t_centers, \(x) apply_window_ttest(x, df))) %>%
    unnest_wider(d) %>%
    mutate(p_adj = p.adjust(p, method="fdr")) %>% # Benjamini-Hochberg//FDR
    mutate(sig = as.double(p<alpha) * sign(beta))
  
  run_lengths <- rle(results_all$sig)
  # results_all$cluster_size <- rep(run_lengths$lengths, run_lengths$lengths)
  # results_all$cluster_val <- rep(run_lengths$values, run_lengths$lengths)
  # 
  results_all <- 
    results_all %>%
    mutate(cluster_size = rep(run_lengths$lengths, run_lengths$lengths), 
           cluster_val = rep(run_lengths$values, run_lengths$lengths) ) %>%
    mutate(sig = sig * as.double(cluster_size>min_cluster_size)) 
  
  
  results <- 
    results_all %>%
    filter(sig!=0) %>%  # for LFPs--require beta>0?
    mutate(ismaxamplitude = (abs(beta)==max(abs(beta)))) 
  
  # r2 %>%
  #   # mutate(p = log(p)) %>%
  #   ggplot(aes(t_centers, beta, color=sig)) + 
  #   geom_point() +
  #   geom_hline(yintercept = 0)
  
  if (nrow(results)==0) {
    earliest_time <- NA
    earliest_time_amplitude <- NA
    peak_time <- NA
    peak_amplitude <- NA
    
    earliest_time_INC <- NA
    peak_time_INC <- NA
    peak_amplitude_INC <- NA
  } else {
    earliest_time <- min(results %>% pull(t_centers))
    earliest_time_amplitude <- results$beta[results$t_centers==earliest_time] 
    peak_time <- results %>% filter(ismaxamplitude) %>% pull(t_centers)
    peak_amplitude <- results %>% filter(ismaxamplitude) %>% pull(beta)
    
    earliest_time_INC <- min(results %>% filter(beta>0) %>%                  pull(t_centers))
    peak_time_INC <-         results %>% filter(beta>0 & ismaxamplitude) %>% pull(t_centers)
    peak_amplitude_INC <-   results %>% filter(beta>0 & ismaxamplitude) %>% pull(beta)
  }
  
  # have to do some shenanigans because tibble(A, B) returns 0x2 table if EITHER 
  # A, B is numeric(0)
  safe_scalar <- function(x) {
    if (length(x) == 0) return(NA_real_)
    else return(x)
  }
  results_summary <- 
    tibble(
      earliest_time = safe_scalar(earliest_time),
      earliest_time_amplitude = safe_scalar(earliest_time_amplitude),
      peak_time = safe_scalar(peak_time),
      peak_amplitude = safe_scalar(peak_amplitude),
      earliest_time_INC = safe_scalar(earliest_time_INC),
      peak_time_INC = safe_scalar(peak_time_INC),
      peak_amplitude_INC = safe_scalar(peak_amplitude_INC)
    )
  # results_summary <- list(earliest_time, earliest_time_amplitude,  peak_time, peak_amplitude, 
  #                           earliest_time_INC, peak_time_INC, peak_amplitude_INC)
  # browser()
  
  list("results_summary"=results_summary, "slwin_results"=results_all)
  
}

# lmtime_results_channelwise$data[[43]] %>%
#   mutate(time = floor(timewarp_t0 * 50) / 50) %>%
#   group_by(time) %>%
#   summarize(m = mean(value), 
#             sem = sd(value) / sqrt(n())) %>%
#   ggplot(aes(x = time, y = m)) + 
#   geom_line()
# df <- lmtime_results_channelwise$data[[14]]


# df <- lmtime_results_channelwise_timeresolved$data[[5]]
# lmtime_results_channelwise_slttest[19:21,]
lmtime_results_channelwise_slttest  <- # slttest=sliding t-test
  lmtime_results_channelwise_timeresolved %>%
  ungroup() %>%
  # filter(row_number() %in% sample.int(100, 15)) %>%
  # filter(row_number() %in% c(5, 42, 30)) %>%
  mutate(chid = glue('{subject_id}_run_{run_id}_{channel}')) %>%
  # filter(chid=="DM1023_run_3_macro_Lc") %>% # testing
  mutate(results = map(data, find_first_sig)) %>%
  mutate(results_summary = map(results, "results_summary"), 
         slwin_results = map(results, "slwin_results")) %>%
  unnest_wider(results_summary)
# save tsv since the above computation takes so long 
fname <- glue('A04v5_11_channelwise_lm_timeres_sliding-window-wrt-bsln_test-TTEST_tlock-{tlock}_winsize-{winsize}-causal')  # statistically significant only 
write_tsv(lmtime_results_channelwise_slttest, file=glue('{PATH_FIG}/{fname}_{format(Sys.time(), "%Y%m%d%H%M")}.tsv'))
# reload an old permutation test
# lmtime_results_channelwise_slttest <- read_tsv(glue('{PATH_FIG}/{fname}.tsv'))

lmtime_results_channelwise_slttest  <- # slttest=sliding t-test
  lmtime_results_channelwise_slttest %>%
  mutate(resp_type = case_when(sign(peak_amplitude)==1 ~ 'INC', 
                               sign(peak_amplitude)==-1 ~ 'DEC', 
                               is.na(peak_amplitude) ~ 'ns'), 
         resp_type = fct_relevel(factor(resp_type), "ns", "DEC", "INC")) 


# Prepare data for time-resolved plotting ----
# df <- lmtime_results_channelwise$data[[1]]
time_breaks <- seq(0, 5, length.out = 500)
# timewarp_t0 c(-1, 5.3), time_t0_SPEECH c(-1.2, 1), time_t0_STIMULUS c(-0.5, 1)
# time_t0_RT
tlock <- "timewarp_t0" 
twin <- switch(tlock, timewarp_t0=c(-1, 5.3), 
               time_t0_RT=c(-0.5, 1), 
               time_t0_SPEECH=c(-1.2, 1), 
               time_t0_STIMULUS=c(-0.5, 1))


# sig_epochSPEECH-INC, sig_noise_typeNOISE-INC, active-INC, 
# sig-noise_type, sig-epoch
channel_subset_name <- "active" # active or one of the sig-* above 
channel_subset_polarity <- "" # "", "INC", "DEC" 

# channel_levels  <- lmcondition_results_channelwise %>%
#   # mutate(stat_epochSPEECH = map(m, "stat_epochSPEECH"))
#   mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) %>%
#   mutate(direction = sign(coeff_epochSPEECH)) %>%
#   arrange((direction), (stat_epochSPEECH))
channel_levels <-
  lmtime_results_channelwise_slttest %>%
  arrange(peak_amplitude) %>%
  mutate(channel_id = glue('{subject_id}_{run_id}_{channel}'))
# channel_levels  <- lmtime_results_channelwise %>% 
# mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) # %>%
# arrange((stat_epochSPEECH))

# prepare plotting table
tmp_plot <- 
  lmtime_results_channelwise_timeresolved %>% 
  separate(channel, into = c('type', 'track', 'unit_idx'), remove = FALSE) %>%
  mutate(modality = MODALITY,
         channel_join = if_else(modality=="LFP", channel, glue('{type}_L{track}'))) %>%
  select(data, channel_join) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, target, region_clean),
            by = join_by(subject_id, run_id==lombard_run_id, channel_join==name)) %>%
  # left_join(SUBJECTS_META %>% select(subject_id, target),
  #           by = join_by(subject_id)) %>%
  # mutate(region_clean = target) %>%
  mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) %>%
  mutate(channel_id = factor(channel_id, levels=channel_levels$channel_id)) %>%
  
  # lmtime_results_channelwise_slttest contains info from t-test INC or DEC type channel
  left_join(lmtime_results_channelwise_slttest %>% group_by(across(all_of(vars_channel_id))) %>% select(-data, -model_time)) %>% 
  
  # CHANGED 20240425
  ungroup() %>%
  # filter(row_number() %in% sample.int(100, 15)) %>%
  {if (str_detect(channel_subset_name, "active")) 
  { semi_join(., lmcondition_results_channelwise %>% filter(sig!="none"), by = join_by(subject_id, run_id, channel)) }
    
    else if (channel_subset_name=="sig_epochSPEECH")
    { semi_join(., lmcondition_results_channelwise %>%  filter(sig_epochSPEECH), by = join_by(subject_id, run_id, channel)) }
    
    else if (channel_subset_name=="sig_noise_typeNOISE")
    { semi_join(., lmcondition_results_channelwise %>% filter(sig_noise_typeNOISE), by = join_by(subject_id, run_id, channel)) }
    
    else if (channel_subset_name=="sig_noise_typeNOISExepochSPEECH")
    { semi_join(., lmcondition_results_channelwise %>% filter(sig_noise_typeNOISExepochSPEECH), by = join_by(subject_id, run_id, channel)) }
    
    else if (channel_subset_name=="sig-epoch")
    { semi_join(., lmcondition_results_channelwise %>% filter(sig=="epoch"), by = join_by(subject_id, run_id, channel)) }
    
    else if (channel_subset_name=="sig-noise_type")
    { semi_join(., lmcondition_results_channelwise %>% filter(sig=="noise_type"), by = join_by(subject_id, run_id, channel)) }
    
    else if (channel_subset_name=="sig-both")
    { semi_join(., lmcondition_results_channelwise %>% filter(sig=="both"), by = join_by(subject_id, run_id, channel)) }
    
    else if (channel_subset_name=="sig-none")
    { semi_join(., lmcondition_results_channelwise %>% filter(sig=="none"), by = join_by(subject_id, run_id, channel)) }
    
    else if (channel_subset_name=="all") {.}
    else { stop(glue('channel_subset_name: "{channel_subset_name}" not recognized'))} } %>% # get significant channels only
  {if      (str_detect(channel_subset_polarity, "INC")) {filter(., resp_type=="INC")} 
    else if (str_detect(channel_subset_polarity, "DEC")) {filter(., resp_type=="DEC")}
    else . } %>%
  unnest(data) %>%
  # mutate(noise_type = lag(noise_type, 5, default = last(noise_type))) %>% # shuffle
  # filter(subject_id!="DM1019") %>%
  group_by(subject_id, run_id, channel, trial_id) %>%
  # mutate(eeg = (value - value_BSLN) / value_BSLN) %>%
  mutate(eeg = value_residualtime) %>% # + value_BSLN_run) %>% # TESTING for SU... want to look at absolute firing rate
  mutate(time = floor(!!sym(tlock) *50)/50) %>% 
  filter(!is.na(loudness_lvl)) %>%
  # mutate(time = cut(timewarp_t0_epoch, breaks=time_breaks, labels=time_breaks)) %>% # cut(timewarp_t0_epoch, breaks=100, labels=FALSE)) %>%
  mutate(eeg_smoothed = signal::sgolayfilt(eeg, p = 3, n = 51)) %>%  # p = polynomial order, n = window size (odd number)
  # mutate(eeg_smoothed = eeg) %>%  # p = polynomial order, n = window size (odd number)
  # filter(time>-1 & time<3) %>%
  mutate(clrby = noise_type) %>% # clrby = loudness_lvl
  group_by(subject_id, run_id, channel, time) 



tmp_plot2 <- 
  tmp_plot %>%
  # --- op1
  # mutate(region_clean = "region-all") # no facet/subgroup
  # --- op2
  mutate(region_clean = target) %>% # very coarse assignment, comment for more precise
  mutate(region_clean = case_when(str_detect(region_clean, 'STN') ~ 'STN',
                                  str_detect(region_clean, 'GP') ~ 'GP',
                                  str_detect(region_clean, 'VIM') ~ 'VIM',
                                  .default = 'other')) %>% relocate(region_clean) %>%
  filter(region_clean!="other")
  # --- op2

# # inspect baseline firing rates 
# tmp_plot2 %>%
#   group_by(subject_id, run_id, channel, region_clean) %>%
#   summarize(value_BSLN_run = mean(value_BSLN_run)) %>%
#   ggplot(aes(x=value_BSLN_run, color=region_clean)) + 
#   geom_histogram(aes(y=after_stat(density)))


# Tile plot//heat map of all speech-responsive electrodes  ----
# SU-FR
if (MODALITY=="SU"){
  scale_max = 20
  scale_min = -20
} else if (MODALITY=="LFP") {
  if (CTX_FLAG) {
    scale_max = 0.3
    scale_min = -0.2
  } else {
    scale_max = 0.1
    scale_min = -0.025
  }

}

tmp_plot %>%
  group_by(channel_id, time) %>%
  summarise(
    eeg_mean = mean(eeg_smoothed, na.rm = TRUE),
    eeg_sem = sd(eeg_smoothed, na.rm = TRUE) / sqrt(n())  # Standard error of the mean
  ) %>%
  mutate(eeg_mean = case_when(eeg_mean>scale_max~scale_max, eeg_mean<scale_min~scale_min, .default=eeg_mean)) %>%
  ggplot(aes(x=time, y=channel_id, fill=eeg_mean)) + 
  geom_raster(interpolate = TRUE) + 
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "#FF0000", limits=c(scale_min, scale_max)) + 
  {if (tlock=="time_t0_SPEECH") xlim(-2.5, 2.5)} + # for t
  coord_cartesian(expand = FALSE) + 
  plot_events() 
fname <- glue('channelwise-heatmap-eeg_chan-subset-{channel_subset_name}_tlock-{tlock}') # statistically significant only 
ggsave(file=glue('{PATH_FIG}/A04v5_07_{fname}_shortheight.pdf'), dpi=300, width=6, height=2.5, units='in') # w=6 h=2.5 for LFP active


# Plot over trialtime... two curves (color-by=noise_type), averaged first within channel (grand average or facet channel) ----
colorby <- "loudness_lvl" # or loudness_lvl or noise_type or fcr_lvl
facetby <- "region_all" # region_clean, region_all or channel
tmp_plot3 <-
  tmp_plot2 %>%
  filter(!is.na(!!sym(colorby))) %>%
  # mutate(region_clean = "region-all") %>%
  { if (facetby == "region_all") { mutate(., region_clean = "region-all")
    } else {.} } %>%
  # left join with a summary of the number of channels within each group
  mutate(noise_type = !!sym(colorby)) %>% # CAREFUL when this is uncommented...must change filename below manually
  left_join(. , get_nt_nc_ns_counts(., "region_clean") %>%
              mutate(region_clean_label = glue('{region_clean} ({nt}, {nc}, {ns})')), #region_clean, " (", n, ")")) %>%
            # select(-n), 
            by = join_by(region_clean)) %>% # %>% relocate(region_clean_label)
  # mutate(facet_label = paste0(region_clean, " (", n, ")")) %>%
  # mutate(noise_type = sample(noise_type)) %>% # DANGER shuffle noise type
  group_by(subject_id, run_id, channel, time, noise_type, region_clean_label) %>%
  summarise( # first summarize--across trials within each channel
    eeg_mean = mean(eeg_smoothed, na.rm = TRUE),
    eeg_sem = sd(eeg_smoothed, na.rm = TRUE) / sqrt(n())  # Standard error of the mean
  ) 

tmp_plot3 %>% 
  
  # comment out section for channel-wise plots
  group_by(time, noise_type, region_clean_label) %>% # arrange(.by_group=TRUE) %>%
  summarise( # second summarize--across channels
    eeg_mean_grp = mean(eeg_mean, na.rm = TRUE),
    eeg_sem_grp = sd(eeg_mean, na.rm = TRUE) / sqrt(n()) # Standard error of the mean
  ) %>% rename(eeg_mean = eeg_mean_grp, eeg_sem = eeg_sem_grp) %>%
  
  ggplot(aes(x=time, color=noise_type)) +
  # geom_rect(data = EPOCHS_MEAN_TIMES, aes(xmin = start+0.1, xmax = end-0.1, ymin = -Inf, ymax = Inf, inherit.aes = FALSE, fill = 'grey'), alpha = 0.2) +  # Add shaded regions
  #  + 
  geom_ribbon(aes(ymin = eeg_mean - eeg_sem, ymax = eeg_mean + eeg_sem, fill=noise_type), alpha = 0.3, linewidth=0) +  # Shaded SE ribbon
  geom_line(aes(y=eeg_mean)) +
  
  {if (facetby=="channel") facet_grid(subject_id+run_id+channel~., scales = 'free')} +
  {if (facetby=="region_clean" || facetby=="region_all" ) facet_grid(region_clean_label~., scales = 'fixed')} +

  # stat_summary(geom = "ribbon", fun.data = "mean_cl_normal", aes(fill = factor(noise_type)), alpha=0.3, linewidth=0) +
  scale_color_manual(values = COLORS) +
  scale_fill_manual(values = COLORS) + 
  # scale_fill_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  geom_hline(yintercept = 0) + 
  # {if (tlock=="time_t0_SPEECH") lims(x=c(-2.5, 2.5), y=c(0, 15))} +
  # {if (MODALITY=="LFP" & facetby!="channel") lims(y=c(-.025/2, 0.08))} + # BG-LFP-BGA
  # {if (MODALITY=="LFP") lims(y=c(-.02, 0.09))} + # tmp change for VIM data
  
  {if (SUBCTX_LFP_REGION=="macro-thetaalpha") lims(y=c(-0.15, 0.1))} +
  {if (CTX_FLAG) scale_y_continuous(expand = expansion(-0.2, 0.2))  } +
  
  {if (MODALITY=="SU" & facetby!="channel") lims(y=c(-3, 10))} +
  
  plot_events() +
  coord_cartesian(expand = FALSE) +
  
  ggtitle(glue('subset {sym(channel_subset_name)}')) + 
  theme_bw()
# tmp ---
fname <- glue('grandavg-eeg_chan-subset-{channel_subset_name}_tlock-{tlock}_color-by-{colorby}_facet-{facetby}') # statistically significant only 
ggsave(file=glue('{PATH_FIG}/A04v5_05_{fname}.pdf'), dpi=300, width=5, height=2.5, units='in') # w=5 h=2.5 for 1 panel, w=5 h=4 for 2 panels, h=25 for multichannel
  

# statistical testing on the above data: want to compare loudness_lvl loud and soft 
# around speech onset
assert('code doesnt currently support tlocks other than timewarp_t0', tlock=="timewarp_t0")
time_windows <- list(
  neg0p9to0  = c(-0.9, 0),
  pos2p5to3 = c(2.5, 3),
  pos3to3p5   = c(3, 3.5),
  pos3p5to4   = c(3.5, 4)
)

for (window_name in names(time_windows)) {
  t1 <- time_windows[[window_name]][1]
  t2 <- time_windows[[window_name]][2]
  
  tmp_plot3 %>%
    filter(time>t1 & time<t2) %>%
    # filter(noise_type!="mid3rd") %>%
    filter(subject_id!="DM1015") %>% # testing/debugging
    group_by(subject_id, run_id, channel, noise_type, region_clean_label) %>%
    summarize(eeg_mean = mean(eeg_mean, na.rm = TRUE)) %>%
    ggplot() + 
    aes(x=noise_type, y=eeg_mean, color=noise_type) + 
    geom_boxplot(outliers = F) +
    geom_point(position = position_jitter(seed = 3, width=0.1), alpha=0.5) +
    scale_color_manual(values = COLORS) +
    # geom_point() +
    facet_grid(.~region_clean_label) + 
    geom_line(aes(group=interaction(subject_id, run_id, channel)), alpha=0.5) + 
    scale_y_continuous(expand = expansion(mult = c(0.1,0.15))) +
    {if (colorby=="noise_type")   stat_compare_means(method = 'wilcox.test', comparisons=list(c('QUIET', 'NOISE')), paired = TRUE, size=3)} +
    {if (colorby=="loudness_lvl")   stat_compare_means(method = 'wilcox.test', comparisons=list(c('low3rd', 'mid3rd'), 
                                                                                                c('mid3rd', 'loud3rd'), 
                                                                                                c('low3rd', 'loud3rd')), paired = TRUE, size=3)} +
    {if (colorby=="fcr_lvl")   stat_compare_means(method = 'wilcox.test', comparisons=list(c('low3rdFCR', 'mid3rdFCR'), 
                                                                                                c('mid3rdFCR', 'loud3rdFCR'), 
                                                                                                c('low3rdFCR', 'loud3rdFCR')), paired = TRUE, size=3)} + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file=glue('{PATH_FIG}/A04v5_05_{fname}_summary_{window_name}-mean_noDM1015.pdf'), dpi=300, width=4, height=3, units='in') # w=2 h=3 for one panel, w=4 h=2 for 2 panels
}

tmp_plot3 %>%
  filter(time>3 & time<3.5) %>%
  # filter(noise_type!="mid3rd") %>%
  group_by(subject_id, run_id, channel, noise_type, region_clean_label) %>%
  summarize(eeg_mean = mean(eeg_mean, na.rm = TRUE)) %>%
  pivot_wider(values_from = eeg_mean, names_from = noise_type) 



# Plot grand average, color-by=region, facet-by=resp_type  ---- 

# facet_by <- "all"
plot_single_channel_flag <- FALSE  # or FALSE

tmp_plot3 <- 
  tmp_plot2 %>%
  mutate(facet_by = resp_type) %>% # if no split by region
  left_join(. , get_nt_nc_ns_counts(., "facet_by") %>%
              mutate(facet_by_label = glue('{facet_by} ({nt}, {nc}, {ns})')), #region_clean, " (", n, ")")) %>%
            # select(-n), 
            by = join_by(facet_by)) %>%
    # mutate(facet_by_label = fct_relevel(facet_by_label, rev)) %>%
    mutate(
      facet_by_label = fct_relevel(
        facet_by_label,
        function(lvls) {
          # non‑INC first, then INC so it appears on top in facet_grid
          c(lvls[!grepl("^INC", lvls)], lvls[grepl("^INC", lvls)])
        }
      )
    ) %>%
  filter(!str_detect(facet_by_label, 'ns')) %>% # for resp_type==ns
  # mutate(region_clean_label = case_when(region_clean=="GPi" ~ "GP",
  #                                       is.na(region_clean) ~ "other",
  #                                       .default = region_clean)) %>%
  # filter(region_clean!="other") %>%
  # mutate(facet_by_label = factor(facet_by_label))
  group_by(subject_id, run_id, channel, time, region_clean, facet_by_label) %>%
  summarise(
    eeg_mean = mean(eeg_smoothed, na.rm = TRUE),
    eeg_sem = sd(eeg_smoothed, na.rm = TRUE) / sqrt(n())  # Standard error of the mean
  )



# tmp_plot3 %>%
#   mutate(facet_by_label = fct_infreq(facet_by_label)) %>%
#   
#   group_by(time, region_clean, facet_by_label) %>% # arrange(.by_group=TRUE) %>%
#   summarise(
#     eeg_mean_grp = mean(eeg_mean, na.rm = TRUE),
#     eeg_sem_grp = sd(eeg_mean, na.rm = TRUE) / sqrt(n()) # Standard error of the mean
#   ) %>% rename(eeg_mean = eeg_mean_grp, eeg_sem = eeg_sem_grp) %>%
#   # filter(facet_label != "all (100)") %>%
#   
#   ggplot(aes(x=time, color=region_clean)) +
#   # geom_rect(data = EPOCHS_MEAN_TIMES, aes(xmin = start+0.1, xmax = end-0.1, ymin = -Inf, ymax = Inf, inherit.aes = FALSE, fill = 'grey'), alpha = 0.2) +  # Add shaded regions
#   #  + 
#   geom_ribbon(aes(ymin = eeg_mean - eeg_sem, ymax = eeg_mean + eeg_sem, fill=region_clean), alpha = 0.3, linewidth=0) +  # Shaded SE ribbon
#   geom_line(aes(y=eeg_mean)) +
#   # facet_grid(subject_id+run_id+channel~., scales = 'free') +
#   # facet_grid(resp_type~., scales = 'free') +
#   facet_grid(facet_by_label~., scales = 'fixed') +
#   # stat_summary(geom =f "ribbon", fun.data = "mean_cl_normal", aes(fill = factor(noise_type)), alpha=0.3, linewidth=0) +
#   scale_color_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) +
#   scale_fill_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) +
#   geom_hline(yintercept = 0) + 
#   plot_events() +
#   coord_cartesian(expand = FALSE) +
#   # ylim(-0.025, 0.05) +
#   {if (tlock=="time_t0_SPEECH" & MODALITY=="LFP") list(xlim(-1, 1), ylim(-0.1, 0.1))} + # + ylim(-0.1, 0.1)} + # for speech, tlock <- "time_t0_SPEECH
#   ggtitle(glue('subset {sym(channel_subset_name)}')) + 
#   theme_bw()
tmp_plot3 %>%
  mutate(facet_by_label = fct_infreq(facet_by_label)) %>%
  ggplot(aes(x = time, color = region_clean)) +
  
  {if (plot_single_channel_flag) {
    list(
      geom_line(aes(y = eeg_mean, group = interaction(subject_id, run_id, channel)), alpha = 0.6, linewidth = 0.7)
    )
  } else {
    tmp_plot3_summary <- tmp_plot3 %>%
      group_by(time, region_clean, facet_by_label) %>%
      summarise(
        eeg_mean_grp = mean(eeg_mean, na.rm = TRUE),
        eeg_sem_grp = sd(eeg_mean, na.rm = TRUE) / sqrt(n()) # Standard error of the mean
      ) %>% rename(eeg_mean = eeg_mean_grp, eeg_sem = eeg_sem_grp)
      # filter(facet_label != "all (100)") %>%
    #   group_by(time, region_clean, facet_by_label) %>% # arrange(.by_group=TRUE) %>%
    #   summarise(

    
    list(
      geom_ribbon(data = tmp_plot3_summary,
                  aes(ymin = eeg_mean - eeg_sem, ymax = eeg_mean + eeg_sem, fill = region_clean),
                  alpha = 0.3, linewidth = 0),
      geom_line(data = tmp_plot3_summary, aes(y = eeg_mean))
    )
  }} +
  
  facet_grid(facet_by_label ~ ., scales = 'fixed') +
  scale_color_manual(values = color_map_targets_simple$hex, breaks = color_map_targets_simple$name) +
  scale_fill_manual(values = color_map_targets_simple$hex, breaks = color_map_targets_simple$name) +
  geom_hline(yintercept = 0) +
  plot_events() +
  coord_cartesian(expand = FALSE) +
  {
    if (tlock %in% c("time_t0_SPEECH", "time_t0_RT") & MODALITY == "LFP") {
      list(ylim(-0.015, 0.03))
    } else if (tlock %in% c("time_t0_SPEECH", "time_t0_RT") & MODALITY == "SU") {
      list(ylim(-6, 15))
    }
  } +
  xlim(twin[1], twin[2]) + 
  ggtitle(glue('subset {sym(channel_subset_name)}')) +
  theme_bw() -> g
g
fname <- glue('A04v5_08_grandavg-eeg_chan-subset-{channel_subset_name}_tlock-{tlock}_colorby-region-clean_facet-by-resp-type') # statistically significant only 
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=4, height=4, units='in') # w=4 h=4 for two facet, 


tmp_plot3 %>%
  group_by(region_clean) %>% #   group_by(resp_type, region_clean_label) %>%
  select(subject_id, run_id, channel) %>% distinct() %>% 
  summarize(n_chan = n(), 
            n_subj = n_distinct(subject_id)) %>%
  write_tsv(., file=glue('{PATH_FIG}/{fname}_summary.tsv'))

# does STN activate before GP? // are there differences around speech onset
num_boot <- 100
alpha <- 0.05

bootstrap_diff <- function(df, nboot = num_boot) {
  replicate(nboot, {
    stn <- df %>% filter(region_clean == "STN") %>% pull(value) %>% sample(replace = TRUE)
    gp  <- df %>% filter(region_clean == "GP")  %>% pull(value) %>% sample(replace = TRUE)
    
    mean(gp, na.rm = TRUE) - mean(stn, na.rm = TRUE)
  })
}

df_boot_input <- tmp_plot3 %>%
  filter(region_clean %in% c("STN", "GP"), 
         time > twin[1], time < twin[2]) %>%
  group_by(time, subject_id, run_id, channel, region_clean) %>%
  summarize(value = mean(eeg_mean, na.rm = TRUE), .groups = "drop")

boot_results <- df_boot_input %>%
  group_by(time) %>%
  nest() %>%
  mutate(
    boots = map(data, bootstrap_diff),
    boot_mean = map_dbl(boots, mean),
    ci_lower = map_dbl(boots, ~ quantile(.x, alpha / 2)),
    ci_upper = map_dbl(boots, ~ quantile(.x, 1 - alpha / 2))
  ) %>%
  select(time, boot_mean, ci_lower, ci_upper)

boot_results <- boot_results %>%
  mutate(sig = ci_lower > 0 | ci_upper < 0)

ggplot(boot_results, aes(x = time, y = boot_mean)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = COLORS["GPi"], alpha = 0.3) +
  geom_line(color = COLORS["GPi"], linewidth = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_segment(data = boot_results %>% filter(sig),
               aes(x = time, xend = time, y = -0.005/2, yend = 0.005/2),  # adjust y range as needed
               color = "purple", linewidth = 1.5) +
  labs(
    title = "Bootstrapped STN vs GP Difference Over Time",
    x = "Time (s)",
    y = "Mean(GP - STN)"
  ) + 
  coord_cartesian(expand = FALSE) +
  plot_events() + 
  xlim(twin[1], twin[2]) +
  theme_bw()
fname <- glue('A04v5_08b_grandavg-eeg_chan-subset-{channel_subset_name}_tlock-{tlock}_colorby-region-clean_facet-by-resp-type_difference') # statistically significant only 
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=2.5, height=2, units='in') # w=4 h=4 for two facet



# plot grand average NOISE vs QUIET difference after taking mean within channel ----- 
tmp_plot %>%
  group_by(subject_id, run_id, channel, epoch, time, noise_type) %>%
  summarise(
    eeg_mean = mean(eeg_smoothed, na.rm = TRUE),
    eeg_sem = sd(eeg_smoothed, na.rm = TRUE) / sqrt(n())  # Standard error of the mean
  ) %>% 
  group_by(time, epoch, noise_type) %>% # arrange(.by_group=TRUE) %>%
  # summarize(n = n())
  summarise(
    eeg_mean_grp = mean(eeg_mean, na.rm = TRUE),
    eeg_sem_grp = sd(eeg_mean, na.rm = TRUE) / sqrt(n()),  # Standard error of the mean
    n = n(), 
  ) %>% rename(eeg_mean = eeg_mean_grp, eeg_sem = eeg_sem_grp) %>%
  pivot_wider(values_from = eeg_mean, names_from = noise_type) %>%
  mutate(eeg_mean_NOISEminusQUIET = NOISE - QUIET) %>%
  group_by(epoch, time) %>%
  summarise(
    eeg_mean = mean(eeg_mean_NOISEminusQUIET, na.rm = TRUE),
    eeg_sem = sd(eeg_mean_NOISEminusQUIET, na.rm = TRUE) / sqrt(n())  # Standard error of the mean
  ) %>%
  
  ggplot(aes(x=time, color=noise_type)) +
  # geom_rect(data = EPOCHS_MEAN_TIMES, aes(xmin = start+0.1, xmax = end-0.1, ymin = -Inf, ymax = Inf, inherit.aes = FALSE, fill = 'grey'), alpha = 0.2) +  # Add shaded regions
  #  + 
  geom_ribbon(aes(ymin = eeg_mean - eeg_sem, ymax = eeg_mean + eeg_sem, fill=noise_type), alpha = 0.3, linewidth=0) +  # Shaded SE ribbon
  geom_line(aes(y=eeg_mean)) +
  # facet_grid(subject_id+run_id+channel~., scales = 'free') +
  facet_grid(.~., scales = 'free_y') +
  # stat_summary(geom = "ribbon", fun.data = "mean_cl_normal", aes(fill = factor(noise_type)), alpha=0.3, linewidth=0) +
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  scale_fill_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  geom_hline(yintercept = 0) + 
  coord_cartesian(expand = FALSE) + 
  ylim(-0.025, 0.05) +
  plot_events() +
  # xlim(-1, 2.5) + 
  theme_bw()
# fname <- 'grandavg-eeg_split-by-noise-type_task-difference-responsive-only_facet-epoch_no-lmtime-valueBSLN_channelwise' # statistically significant only 
fname <- 'grandavg-eeg_split-by-noise-type_task-difference-responsive-only_lmtime-no-valueBSLN_channelwise' # statistically significant only 
ggsave(file=glue('{PATH_FIG}/A04v5_06_{fname}.pdf'), dpi=300, width=6, height=3, units='in') # height=25 for multi-channel plots


# Plot boxplots for different speech epochs to visualize LME ---- 
channel_subset_name <- "sig_epochSPEECH"

channel_levels  <- lmcondition_results_channelwise %>% 
  mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) %>%
  arrange((stat_epochSPEECH))


tmp_plot <- 
  lmtime_results_channelwise %>% 
  select(data) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, region_clean), 
            by = join_by(subject_id, run_id==lombard_run_id, channel==name)) %>%
  mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) %>%
  group_by(channel_id) %>%
  ungroup() %>%
  # filter(row_number() %in% sample.int(100, 15)) %>%
  semi_join(lmcondition_results_channelwise,
            # filter(!!sym(channel_subset_name) & !str_detect(channel, 'ecog')),
            by = join_by(subject_id, run_id, channel)) %>% # get significant channels only
  unnest(data) %>%
  # mutate(noise_type = lag(noise_type, 5, default = last(noise_type))) %>% # shuffle
  # filter(subject_id!="DM1019") %>%
  group_by(subject_id, run_id, channel, trial_id) %>%
  # mutate(eeg = (value - value_BSLN) / value_BSLN) %>%
  mutate(eeg = value_residualtime) %>%
  mutate(time = floor(timewarp_t0*50)/50) %>% 
  # filter(!is.na(loudness_lvl)) %>%
  mutate(clrby = noise_type) %>%
  group_by(subject_id, run_id, channel, time) %>%
  # mutate(channel_id = factor(channel_id, levels=channel_levels$channel_id)) %>% 
  
  tmp_plot2 <-
  tmp_plot %>%
  mutate(epoch_factor = case_when(epoch=="RT" ~ "PREP", 
                                  epoch=="BSLN" ~ "BSLN", 
                                  epoch=="SPEECH" & timeperc_t0_epoch<=0.33 ~ "BEG",
                                  epoch=="SPEECH" & timeperc_t0_epoch>0.33 & timeperc_t0_epoch<=0.66 ~ "MID",
                                  epoch=="SPEECH" & timeperc_t0_epoch>0.66 ~ "END", 
                                  .default = NA)) %>%
  mutate(epoch_factor = factor(epoch_factor, levels = c("BSLN", "PREP", "BEG", "MID", "END"))) %>%
  mutate(region_clean_label = case_when(str_detect(region_clean, "GP") ~ "GP",
                                        str_detect(region_clean, "STN") ~ "STN",
                                        .default = NA)) %>%
  filter(!is.na(epoch_factor)) %>%
  filter(!is.na(region_clean_label)) %>%
  group_by(channel_id, epoch_factor, region_clean_label) %>%
  # mutate(eeg = log10(eeg)) %>%
  summarize(eeg = mean(eeg, na.rm = TRUE)) 

tmp_plot2 %>%    
  ggplot(aes(x=epoch_factor, y=eeg, color=region_clean_label)) + 
  # geom_line(aes(group=channel_id, alpha=0.2)) + 
  geom_boxplot(aes(group=interaction(epoch_factor, region_clean_label)), outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(), alpha=0.7) + # position = position_jitterdodge()
  theme_bw() + 
  scale_color_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) +
  scale_fill_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) 
fname <- glue('boxplots-across-speech-epochs_subset-{channel_subset_name}_split-by-region_channelwise') # statistically significant only 
ggsave(file=glue('{PATH_FIG}/A04v5_09_{fname}.pdf'), dpi=300, width=6, height=3, units='in') # height=25 for multi-channel plots
write_tsv(tmp_plot2, file=glue('{PATH_FIG}/A04v5_09_{fname}.tsv'))


tmp_plot3 <-
  tmp_plot2 %>% 
  pivot_wider(names_from=epoch_factor, values_from = eeg)
t.test(tmp_plot3$PRE, tmp_plot3$BSLN, paired = TRUE) # t = -0.44278, df = 29, p-value = 0.6612
t.test(tmp_plot3$BEG, tmp_plot3$BSLN, paired = TRUE) # t = 6.1622, df = 29, p-value = 1.022e-06
t.test(tmp_plot3$MID, tmp_plot3$BSLN, paired = TRUE) # t = 4.4308, df = 29, p-value = 0.0001229
t.test(tmp_plot3$END, tmp_plot3$BSLN, paired = TRUE) # t = 3.9364, df = 29, p-value = 0.0004754





# Note PLB 20250116... DM1019 for some reason has unreasonably high values 
# taking the max over trials..
# # Groups:   subject_id, run_id [25]
# subject_id run_id channel  eeg_smoothed
# <chr>       <dbl> <chr>           <dbl>
#   1 DM1019          3 macro_Lc      2888.  
# 2 DM1019          3 macro_Lp       246.  
# 3 DM1019          3 macro_Lm        97.1 
# 4 DM1019          2 macro_Lc        91.7 
# 5 DM1019          4 macro_Lp        23.7 
# 6 DM1027          4 macro_Lp         9.91
# 7 DM1034          2 macro_Lc         4.66
# 8 DM1028          4 macro_Ll         3.71
# 9 DM1019          2 macro_Lp         3.47
# 10 DM1028          4 macro_Lp         3.39



df %>%
  mutate(time_bin = floor(time_t0*100)) %>%
  mutate(audioR_s = (audioR_s - mean(audioR_s)) / sd(audioR_s)) %>%
  # filter(epoch=="BSLN") %>%
  rename(eeg_res=value_residualtime, 
         eeg=value) %>%
  # mutate(time_t0=time_t0_SPEECH) %>%
  # filter(epoch %in% c("BSLN", "SPEECH")) %>%
  group_by(trial_id, epoch, noise_type) %>%
  summarize(eeg_res = mean(eeg_res), 
            eeg = mean(eeg), 
            time_GTC = mean(time_GTC)) %>%
  # mutate(eeg_diff = eeg_res - eeg)
  # select(trial_id, noise_type, eeg, eeg_res, time_t0, epoch) %>%
  # filter(trial_id < 10) %>%
  pivot_longer(cols = c("eeg_res"), names_to = "series") %>% #   pivot_longer(cols = c("audioR_s", "eeg")) %>%
  # group_by(time_bin, epoch) %>%
  # summarise(valuebin = mean(valuelog)) %>%
  # summarise(valuebin = mean(value_residualtime), intensity = mean(audioR_s)) %>%
  # mutate(intensity_logit = log(intensity / (1 - intensity))) %>%
  # pivot_longer(cols = c("valuebin", "intensity"), names) %>%
  ggplot(aes(x=epoch, y=value, color=factor(noise_type))) +
  # geom_line() +
  # facet_wrap(~trial_id)
  stat_summary(geom = "ribbon", fun.data = "mean_cl_normal", aes(fill = factor(name)), alpha=0.3, linewidth=0) +
  # stat_summary_bin(geom = "ribbon", fun.data = "mean_cl_normal", aes(fill = factor(noise_type)), alpha=0.3, linewidth=0, binwidth=0.100) +
  # stat_summary_bin(geom = "line", fun.data = "mean_cl_normal",   aes(color = factor(noise_type)), binwidth=0.100)
  # geom_boxplot() +
  # geom_jitter(aes(alpha=0.3)) +
  # geom_line()  +
  # geom_point()
  # geom_smooth(span = )
  # geom_bin2d(bins = 60)
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  # scale_fill_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  # ylim(-0.25, 0.25) + 
  # xlim(-1, 3)
  
  
  
  
  #   
  # df %>%
  #   group_by(time_t0, noise_type) %>%
  #   summarize(
  #     mean_eeg = mean(residual_time, na.rm = TRUE),
  #     se_eeg = sd(residual_time, na.rm = TRUE) / sqrt(n())  # Standard error calculation
  #   ) %>%
  #   ggplot(aes(x=time_t0, group=noise_type)) +
  #   geom_ribbon(aes(ymin = mean_eeg - se_eeg, ymax = mean_eeg + se_eeg), fill = "red", alpha = 0.3) +  # Shaded SE ribbon
  #   geom_smooth(span=0.1, aes(y = mean_eeg), color = "blue", size = 1)   # Mean trace
  # 
  #   geom_smooth(span=0.05)
  # 
  # 
  
  
  
  # print subjects from the DM10XX series which are NOT present in the dataset 
  tibble(subject_id=paste0('DM',subject_ids)) %>%
  anti_join(trials %>% ungroup() %>% distinct(subject_id))
# DM1001	research aborted
# DM1004	research aborted
# DM1007	ET/SMSL
# DM1009	?? too few trials? TFR fieldtrip exists
# DM1010	no lombard
# DM1014	?? TFR fieldtrip exists, only 17 good trials
# DM1016	research cancelled
# DM1024	ET/SMSL
# DM1025	ET/SMSL
# DM1026	ET/SMSL
# DM1030	research aborted
# DM1031	research aborted
# DM1032	ET/SMSL
# DM1034	?? TFR fieldtrip exists, but all labeled bad channel
# DM1035	ET/SMSL
# DM1036	no data collected--FINE
# DM1038	ET/SMSL
# DM1039	ET/SMSL


trials <-
  trials %>%
  group_by(subject_id, run_id, channel_id) %>%
  arrange(subject_id, run_id, channel_id, trial_id, desc(timestamp))

# remove duplicate rows
duplicates <- trials[trials %>% select(subject_id, run_id, channel_id, trial_id) %>% duplicated(), ]
if (nrow(duplicates)>0) {
  warning(paste("There are duplicate rows in the data."))
  print(duplicates)
  
  trials <-trials %>% group_by(subject_id, run_id, channel_id, trial_id) %>% filter(row_number()==1)
}

# massage trials, pivot longer
trials <-
  trials %>%
  mutate(ref_type = ifelse(str_detect(channel, '-'), 'bipolar' , 'monopolar')) %>%
  filter(!str_detect(channel, 'CA')) %>%
  left_join(SUBJECTS_META %>% select(subject_id, dbs_target), by=join_by(subject_id)) %>%
  mutate(target=case_when(str_detect(dbs_target, regex('STN', ignore_case = TRUE)) ~ 'STN', 
                          str_detect(dbs_target, regex('GPi', ignore_case = TRUE)) ~ 'GPi')) %>%
  mutate(noise_type=case_when(noise_type==0 ~ 'quiet', 
                              noise_type==1 ~ 'lombard')) %>%
  select(where(negate(is.numeric)), noise_type, trial_id, run_id,
         contains('mean')) %>%
  pivot_longer(cols = beta_basel_mean:last_col(), 
               names_to=c('band', 'timewin', 'valtype'), 
               names_sep='_', 
               values_to = "value") %>%
  mutate(noise_type=factor(noise_type, levels=c("quiet", "lombard")), 
         timewin=factor(timewin, levels=c("basel", "prespeech", "speech")))


trials <-
  trials %>%
  left_join(sentences %>% select(subject_id, run_id, trial_id, time_withinrun, intensity_praat), by=vars_trial_id)

# trials %>% filter(is.na(intensity_praat)) %>% ungroup() %>% distinct(subject_id)



# Verify sliding window with plots of individual electrodes ----

# fname <- glue('A04v5_08_grandavg-eeg_chan-subset-{channel_subset_name}_colorby-region-clean_facet-by-resp-type') # statistically significant only 
# ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=4, height=4, units='in') # w=4 h=4 for two facet
plot_single_channel_over_trialtime <- function(df, df_slwindow, chid) {
  # df <- lmcondition_results_channelwise$data[[3]]
  # mdl <- lmcondition_results_channelwise$model_conditions[[3]]
  
  df_slwindow <- 
    df_slwindow %>%
    mutate(y_constant=0, 
           region_clean = "na", 
           sig = factor(sig))
  
  df <-
    df %>%
    mutate(eeg = value_residualtime) %>%
    mutate(time = floor(!!sym(tlock) *50)/50) %>% 
    mutate(region_clean = "na") %>%
    filter(!is.na(loudness_lvl)) %>%
    group_by(region_clean, time) %>% # arrange(.by_group=TRUE) %>%
    # mutate(time = cut(timewarp_t0_epoch, breaks=time_breaks, labels=time_breaks)) %>% # cut(timewarp_t0_epoch, breaks=100, labels=FALSE)) %>%
    # mutate(eeg_smoothed = signal::sgolayfilt(eeg, p = 3, n = 51)) %>%  # p = polynomial order, n = window size (odd number)
    mutate(eeg_smoothed = eeg) %>%
    summarise(
      eeg_mean = mean(eeg_smoothed, na.rm = TRUE),
      eeg_sem = sd(eeg_smoothed, na.rm = TRUE) / sqrt(n())  # Standard error of the mean
    )
  
  ggplot(data=df, aes(x=time, color=region_clean)) +
    # geom_rect(data = EPOCHS_MEAN_TIMES, aes(xmin = start+0.1, xmax = end-0.1, ymin = -Inf, ymax = Inf, inherit.aes = FALSE, fill = 'grey'), alpha = 0.2) +  # Add shaded regions
    #  + 
    geom_ribbon(aes(ymin = eeg_mean - eeg_sem, ymax = eeg_mean + eeg_sem, fill=region_clean), alpha = 0.3, linewidth=0) +  # Shaded SE ribbon
    geom_line(aes(y=eeg_mean)) +
    # facet_grid(subject_id+run_id+channel~., scales = 'free') +
    # facet_grid(resp_type~., scales = 'free') +
    # facet_grid(facet_by_label~., scales = 'fixed') +
    # stat_summary(geom = "ribbon", fun.data = "mean_cl_normal", aes(fill = factor(noise_type)), alpha=0.3, linewidth=0) +
    scale_color_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) +
    scale_fill_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) +
    geom_hline(yintercept = 0) + 
    coord_cartesian(expand = FALSE) + 
    geom_point(data=df_slwindow %>% filter(sig!=0), aes(x=t_centers, y=y_constant, shape=sig), color="red", size=3) +
    # ylim(-0.025, 0.05) +
    # xlim(-1, 2.5) + 
    ggtitle(glue('{chid}')) + 
    theme_bw()
  
  ggsave(filename=glue('{PATH_FIG}/A04v5_11b_channelwise-over-trialtime-slwindow/{chid}_bga-trace.pdf'), dpi=300, width=5, height=3, units='in')
  
  df
}

# run the above function for each channel
lmtime_results_channelwise_timeresolved %>%
  left_join(lmtime_results_channelwise_slttest %>% select(-data, -model_time), 
            by=vars_channel_id) %>%
  mutate(chid = glue('{subject_id}_run_{run_id}_{channel}')) %>%
  # filter(chid=="DM1023_run_3_macro_Lc") %>%
  mutate(result = pmap(list(data, slwin_results, chid), plot_single_channel_over_trialtime))


  
  
# Plot sliding window statistics ---- 
channel_subset_name <- "active"
resp_type_filter <- c("INC") # , "DEC", "ns") # "INC", "DEC", "ns"
var_of_interest <- "peak_time" # earliest_time, peak_time, peak_amplitude, 
num_permutations <- 1000 

lmtime_results_channelwise3 <- 
  lmtime_results_channelwise_slttest %>%
  # filter(row_number() %in% sample.int(100, 15)) %>%
  {if (str_detect(channel_subset_name, "active")) 
  { semi_join(., lmcondition_results_channelwise %>% 
                filter(sig!="none"),
              by = join_by(subject_id, run_id, channel)) }
    else if (str_detect(channel_subset_name, "sig_epochSPEECH")) 
    { semi_join(., lmcondition_results_channelwise %>% 
                  filter(sig_epochSPEECH),
                by = join_by(subject_id, run_id, channel)) }
    else if (str_detect(channel_subset_name, "sig_noise_typeNOISE"))
    { semi_join(., lmcondition_results_channelwise %>% 
                  filter(sig_noise_typeNOISE),
                by = join_by(subject_id, run_id, channel)) }
    else if (channel_subset_name=="all") {.}
    else { stop(glue('channel_subset_name: "{channel_subset_name}" not recognized'))} } %>% # get significant channels only
  # {if      (str_detect(channel_subset_name, "-INC")) {filter(., resp_type=="INC")} 
  #   else if (str_detect(channel_subset_name, "-DEC")) {filter(., resp_type=="DEC")}
  #   else . } %>%
  # unnest(data) %>%
  select(-data, -model_time) %>%
  filter(str_detect(resp_type, paste(resp_type_filter, collapse="|"))) %>%
  # filter(earliest_time_amplitude>0) %>% # TESTING
  # mutate(earliest_time_INC = earliest_time_INC -3, # careful--minus 3 is only for time-warped data
  #        peak_time =     peak_time - 3 ) %>% # 
  separate(channel, into = c('type', 'track', 'unit_idx'), remove = FALSE) %>%
  mutate(modality = MODALITY,
         channel_join = if_else(modality=="LFP", channel, glue('{type}_L{track}'))) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, target, region_clean),
            by = join_by(subject_id, run_id==lombard_run_id, channel_join==name)) %>%
  mutate(region_clean = target) %>% # coarsest version
  mutate(region_clean_label = case_when(str_detect(region_clean, "GP") ~ "GP",
                                        str_detect(region_clean, "STN") ~ "STN",
                                        .default = NA)) %>% # less coarse version--but leaves out some data
  arrange(earliest_time_INC) 
  # filter(!(region_clean_label=="GP" & earliest_time_INC < -0.5))
# semi_join(lmcondition_results_channelwise %>%
#             filter(!str_detect(channel, 'ecog')),
#           # filter(!!sym(channel_subset_name) & !str_detect(channel, 'ecog')),
#           by = vars_channel_id) # get significant channels only

# determine what our sample size should be in the permutation tests
lmtime_results_channelwise3 %>%
  group_by(region_clean_label) %>%
  summarize(n = n(), 
            earliest_time_INC_mean = mean(earliest_time_INC, na.rm=T), 
            earliest_time_INC_med = median(earliest_time_INC, na.rm=T))

lmtime_results_channelwise3 %>%
  filter(!is.na(region_clean_label)) %>%
  ggplot(aes(x=region_clean_label, y=earliest_time_INC, color=region_clean_label)) + 
  geom_boxplot(outliers=F) +
  geom_jitter(alpha=0.8, width = 0.1) + 
  stat_compare_means() +
  # ylabel() + 
  scale_color_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) 

# x <-  lmtime_results_channelwise3 %>% filter(region_clean_label=="GP") %>% pull(earliest_time_INC) 
# y <-  lmtime_results_channelwise3 %>% filter(region_clean_label=="STN") %>% pull(earliest_time_INC) 
# ks.test(x, y) 

# leveneTest(y, group)
# install.packages("car")
# library(car)
# result = leveneTest(earliest_time_INC ~ region_clean_label, lmtime_results_channelwise3)
# print(result)
  
# fname <- glue('sig-{var_of_interest}-bootstrap-{channel_subset_name}-{paste(resp_type_filter, collapse="-")}_STN-and-GP_winsize-{winsize}-casual') # statistically significant only 
# ggsave(file=glue('{PATH_FIG}/A04v5_10_{fname}.pdf'), dpi=300, width=3, height=2, units='in') # height=25 for multi-channel plots
# write_tsv(onset_boot_summary, file=glue('{PATH_FIG}/A04v5_10b_{fname}.tsv'))


# to-do: consider doing this with the boot() package 
nsamp <- 20 # should be set according number of channels in each group
onset_boot <- 
  replicate(num_permutations, {
    permuted_results <- 
      lmtime_results_channelwise3 %>%
      mutate(value = !!sym(var_of_interest)) %>%
      group_by(region_clean_label) %>%
      slice_sample(n=nsamp, replace=TRUE) %>% 
      summarize(value = mean(value, na.rm = TRUE)) %>%
      mutate(sample_ID = sample.int(1e9, 1))
    # sum(permuted_results)
  }, simplify = FALSE)
onset_boot <- bind_rows(onset_boot)

onset_boot_summary <-
  onset_boot %>%
  group_by(region_clean_label) %>%
  summarise(lowerbound = quantile(value, alpha/2),
            upperbound = quantile(value, 1- alpha/2),
            mean = mean(value),
            n = n())
# region_clean_label lowerbound upperbound    mean     n
# <chr>                   <dbl>      <dbl>   <dbl> <int>
#   1 GP                     -0.181     0.0344 -0.0722  1000
# 2 STN                    -0.360    -0.174  -0.273   1000
# 3 NA                     -0.404    -0.217  -0.312   1000

onset_boot %>%
  # filter(!is.na(region_clean_label)) %>%
  filter(region_clean_label %in% c("STN", "GP")) %>% # "STN", 
  ggplot(aes(x=value, fill=region_clean_label)) +
  geom_histogram(aes(y=after_stat(density)), alpha=0.9) + # aes(y=after_stat(density)) + stat = 'density'
  # geom_density(alpha=0.5) +
  scale_fill_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) + 
  geom_vline(xintercept = 0) + 
  # theme(legend.position="none")
  theme_bw() 
fname <- glue('sig-{var_of_interest}-bootstrap_tlock-{tlock}_{channel_subset_name}-{paste(resp_type_filter, collapse="-")}_STN-and-GP_winsize-{winsize}-casual') # statistically significant only 
ggsave(file=glue('{PATH_FIG}/A04v5_10_{fname}.pdf'), dpi=300, width=4, height=1.7, units='in') # height=25 for multi-channel plots
write_tsv(onset_boot_summary, file=glue('{PATH_FIG}/A04v5_10_{fname}.tsv'))


# now the same, but for the difference between GP and STN onset  
onset_boot_diff <-
  onset_boot %>% 
  pivot_wider(names_from = region_clean_label, values_from = value) %>%
  mutate(value = GP-STN,
         region_clean_label = "GP")

onset_boot_diff_summary <- 
  onset_boot_diff %>%
  summarise(lowerbound = quantile(value, alpha/2), 
            upperbound = quantile(value, 1- alpha/2), 
            mean = mean(value)) 
# lowerbound upperbound  mean
# 1     0.0996      0.309 0.205

onset_boot_diff %>%
  ggplot(aes(x=value, fill = region_clean_label)) +
  geom_histogram(aes(y=after_stat(density))) + 
  geom_vline(xintercept = 0) + 
  scale_fill_manual(values = color_map_targets_simple$hex, breaks=color_map_targets_simple$name) + 
  theme_bw()
  # theme(legend.position="none")
fname <- glue('sig-{var_of_interest}-bootstrap_tlock-{tlock}-{channel_subset_name}_GPSTN-diff_winsize-{winsize}-casual') # statistically significant only 
ggsave(file=glue('{PATH_FIG}/A04v5_10_{fname}.pdf'), dpi=300, width=4, height=1.7, units='in') # height=25 for multi-channel plots
write_tsv(onset_boot_diff_summary, file=glue('{PATH_FIG}/A04v5_10_{fname}.tsv'))

# Convert data to tensor format for matlab/dpca ----
# Define constants
conditions <- c("QUIET", "NOISE")
# loudness_lvls <- c("low3rd", "mid3rd", "loud3rd")  # New dimension
# loudness_lvls <- c("lowhalf", "loudhalf")  # New dimension
loudness_lvls <- c("all")

recode_loudness <- function(df) {
  loudness_trialwise <-
    df %>%
    group_by(trial_id) %>%
    filter((time_t0_SPEECH>0) & (time_t0_SPEECH<=1)) %>%
    summarise(acousticspectrum_intensity_trialwise = mean(acousticspectrum_intensity, na.rm=TRUE)) %>%
    mutate(perc = percent_rank(acousticspectrum_intensity_trialwise), 
           loudness_lvl3 = cut(perc, breaks=3, labels=c("low3rd", "mid3rd", "loud3rd")),
           loudness_lvl2 = cut(perc, breaks=2, labels=c("lowhalf", "loudhalf")), 
           loudness_lvl1 = "all") %>%
    select(trial_id, loudness_lvl3, loudness_lvl2, loudness_lvl1)
  
  df %>%
    left_join(loudness_trialwise, by = join_by(trial_id))
}



nTpt <- 100
maxTrialNum <- 40

# Step 1: Extract and bin
df_binned <- 
  lmtime_results_channelwise_timeresolved %>%
  mutate(data = map(data, recode_loudness)) %>%
  mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) %>%
  ungroup() %>%
  unnest(data) %>%
  
  { if (paste(loudness_lvls, collapse='-')=="low3rd-mid3rd-loud3rd") mutate(., loudness_lvl = loudness_lvl3) 
   else if (paste(loudness_lvls, collapse='-')=="lowhalf-loudhalf") mutate(., loudness_lvl = loudness_lvl2) 
   else if (loudness_lvls=="all") mutate(., loudness_lvl = loudness_lvl1) 
   else stop() } %>%
      
  select(channel_id, trial_id, noise_type, epoch, loudness_lvl,
         timewarp_t0, value_residualtime) %>%
  filter(
    noise_type %in% conditions,
    loudness_lvl %in% loudness_lvls,
    timewarp_t0 > -1, timewarp_t0 < 5.3
  ) %>%
  group_by(channel_id, trial_id) %>%
  mutate(time_bin = cut(timewarp_t0, breaks = nTpt, labels = FALSE)) %>%
  group_by(channel_id, noise_type, loudness_lvl, trial_id, time_bin) %>%
  summarise(
    value = mean(value_residualtime, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Assign trial numbers and select relevant columns
df_trials <- df_binned %>%
  mutate(time_bin = as.integer(time_bin)) %>%
  group_by(channel_id, noise_type, loudness_lvl) %>%
  mutate(trial_num = dense_rank(trial_id)) %>%
  filter(trial_num <= maxTrialNum) %>%
  ungroup() %>%
  select(channel_id, noise_type, loudness_lvl, trial_num, time_bin, value)

# # testing-- do we have trials in every loudness level? 
# aa <- df_trials %>% group_by(channel_id, loudness_lvl, noise_type) %>% select(trial_num) %>% distinct() %>% summarize(n = n())
# aa %>%
#   group_by(loudness_lvl, noise_type) %>% 
#   summarize(n = mean(n)) %>%
#   ggplot(aes(x = interaction(loudness_lvl, noise_type), y = n)) + 
#   geom_col() + 
#   # facet_grid(.~channel_id)

# Step 3: Pad missing trials and time bins
df_complete <- df_trials %>%
  complete(
    channel_id,
    noise_type = conditions,
    loudness_lvl = loudness_lvls,
    trial_num = 1:maxTrialNum,
    time_bin = 1:nTpt,
    fill = list(value = NA_real_)
  )

# Step 4: Initialize tensor
channel_ids <- unique(df_complete$channel_id)
N <- length(channel_ids)
S <- length(conditions)
D <- length(loudness_lvls)
T <- nTpt
R <- maxTrialNum

tensor_5d <- array(NA_real_, dim = c(N, S, D, T, R),
                   dimnames = list(channel_ids, conditions, loudness_lvls, NULL, NULL))
tensor_5d_ntrials <- array(NA_real_, dim = c(N, S, D),
                   dimnames = list(channel_ids, conditions, loudness_lvls))

# Step 5: Fill tensor
for (n in seq_along(channel_ids)) {
  ch <- channel_ids[n]
  ch_data <- df_complete %>% filter(channel_id == ch)
  
  for (s in seq_along(conditions)) {
    cond <- conditions[s]
    
    for (d in seq_along(loudness_lvls)) {
      loud <- loudness_lvls[d]
      
      cond_data <- ch_data %>%
        filter(noise_type == cond, loudness_lvl == loud) %>%
        arrange(trial_num, time_bin)
      
      mat <- matrix(NA_real_, nrow = T, ncol = R)
      
      if (nrow(cond_data) > 0) {
        filled_data <- cond_data %>%
          filter(trial_num <= R, time_bin <= T) %>%
          mutate(row = time_bin, col = trial_num)
        
        # browser()
        ntrials <- filled_data %>% group_by(noise_type, loudness_lvl, trial_num) %>% summarize(missing = all(is.na(value))) %>%
          summarise(ntrials = sum(!missing)) %>% pull(ntrials)
        
        mat[cbind(filled_data$row, filled_data$col)] <- filled_data$value
      } else {
        # browser()
      }
      
      tensor_5d[n, s, d, , ] <- mat
      tensor_5d_ntrials[n, s, d] <- ntrials
    }
  }
}



meta <-
  lmtime_results_channelwise_timeresolved %>% 
  separate(channel, into = c('type', 'track', 'unit_idx'), remove = FALSE) %>%
  mutate(modality = MODALITY,
         channel_join = if_else(modality=="LFP", channel, glue('{type}_L{track}'))) %>%
  select(channel_join) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, target, region_clean),
            by = join_by(subject_id, run_id==lombard_run_id, channel_join==name)) %>%
  # left_join(SUBJECTS_META %>% select(subject_id, target),
  #           by = join_by(subject_id)) %>%
  # mutate(region_clean = target) %>%
  mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) 


library(hdf5r)
t <- format(Sys.time(), "%Y%m%d%H%M")
loudness_str <- paste(loudness_lvls, collapse='-')
if (MODALITY=="SU") {file <- H5File$new(glue("{PATH_ANALYSIS}/data/dpca_tensor_timewarped_cofit-lmtime-noisetype_loudness-{loudness_str}_{t}.h5"), mode = "w")
} else              {file <- H5File$new(glue("{PATH_DATA_INPUT}/dpca_tensor_timewarped_cofit-lmtime-noisetype_loudness-{loudness_str}_{t}.h5"), mode = "w")}
file[["X"]] <- tensor_5d
file[["X_ntrials"]] <- tensor_5d_ntrials

for (col_name in names(meta)) {
  file[[col_name]] <- meta[[col_name]]
}
file$close_all()


## Plot respiration signal, lombard vs quiet ----
# Load resp signals
resp <-
  tibble(subject_id_tmp=paste0('DM',subject_ids)) %>%
  mutate(path_annot=glue("{PATH_ANALYSIS}/data/resp/sub-{subject_id_tmp}_ses-intraop_task-lombard_resp.tsv")) %>%
  mutate(exists_annot=file.exists(path_annot)) %>%
  filter(exists_annot) %>%
  mutate(t = map(pull(., path_annot), function(x) read_tsv(x))) %>%
  unnest(t) %>%
  select(-path_annot, -exists_annot, -subject_id_tmp) %>%
  # rename(subject_id = subject_id_tmp) %>% # unnest the tibbles on each row
  group_by(subject_id) %>%
  mutate(noise_type=case_when(noise_type==0 ~ 'QUIET',
                              noise_type==1 ~ 'NOISE')) %>%
  mutate(noise_type=factor(noise_type, levels=c("QUIET", "NOISE")))

resp_plot <-
  resp %>% # filter(subject_id=="DM1033") %>% filter()
  # filter(any_of(c('rms', 'min', 'max', 'range')) >4) %>% (resp_value < 4 & resp_value >-4) %>%)
  pivot_longer(cols=c('rms', 'min', 'max', 'range'), names_to = 'resp_name', values_to = 'resp_value') %>%
  filter(resp_name=='range' & resp_value<7)


# by SUBJECT 
resp_plot %>%
  ggplot(aes(x=subject_id,
             y=resp_value, 
             color=noise_type,
             group=interaction(subject_id, noise_type))) +
  geom_jitter(position=position_dodge(width=0.75,preserve="total"), alpha=0.4) + 
  geom_boxplot(position=position_dodge(width=0.75,preserve="total"),outlier.alpha=0, show.legend = FALSE) + 
  stat_compare_means(
    method = "wilcox.test",  # Specify the method, e.g., t.test
    comparisons = c("NOISE", "QUIET"),  # Adjust based on your groups
    label = "p.signif",  # Display significance as stars
    hide.ns = FALSE  # Hide non-significant comparisons
  ) + 
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  # scale_alpha_manual(values = c("speech" = 0.9, "basel" = 0.4))+ 
  # scale_shape_manual(values = c("speech" = 16, "basel" = 1)) +
  facet_grid(resp_name~., scales = 'free_y')+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

fname = glue('A04v5_13_resp-all-subj-all-metrics')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=10, height=4, units='in')
# write_tsv(trials_tmp, glue('{PATH_FIG}/{fname}.tsv'))

# averaged within-subject 
resp_plot %>%
  group_by(subject_id, noise_type) %>%
  summarize(
    resp_mean = mean(resp_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(noise_type = factor(noise_type, levels = c("QUIET", "NOISE"))) %>%
  ggplot(aes(x = noise_type, y = resp_mean, color = noise_type)) +
  geom_boxplot() + 
  geom_line(aes(group = subject_id), linewidth = 0.8, alpha = 0.6) +
  geom_point(aes(color = noise_type), size = 3, position = position_jitter(width=0.1, seed=)) + 
  stat_compare_means(comparisons = list(c('QUIET', 'NOISE')), method = "wilcox.test", paired = TRUE) +
  scale_color_manual(values = color_map_conditions$hex, breaks = color_map_conditions$name) +
  labs(
    title = "Respiration Range by Condition",
    x = "Condition",
    y = "Mean Respiration Range [norm]"
  ) +
  theme_minimal()
fname = glue('A04v5_13_resp-cohort')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=4, height=4, units='in')





# Old convert code ----
conditions <- c("QUIET", "NOISE")

N <- nrow(lmtime_results_channelwise_timeresolved)
S <- length(conditions) # number of conditions
D <- 1 # 
nTpt <- 100  # number of time points per trial
maxTrialNum <- 40 # 80 lombard trials / 2


# Step 1: Extract and bin
df_binned <- lmtime_results_channelwise_timeresolved %>%
  mutate(channel_id=glue('{subject_id}_{run_id}_{channel}')) %>%
  ungroup() %>%
  unnest(data) %>%
  select(channel_id, trial_id, noise_type, epoch, loudness_lvl,
         timewarp_t0, value_residualtime) %>%
  filter(
    noise_type %in% c("NOISE", "QUIET"),         # only keep desired conditions
    timewarp_t0 > -1, timewarp_t0 < 5.3          # ensure time is in valid range
  ) %>%
  group_by(channel_id, trial_id) %>%
  mutate(time_bin = cut(timewarp_t0, breaks = nTpt, labels = FALSE)) %>%
  group_by(channel_id, noise_type, trial_id, time_bin) %>%
  summarise(
    value = mean(value_residualtime, na.rm = TRUE),
    timewarp_t0 = mean(timewarp_t0, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Assign trial_num (for padding)
df_trials <- df_binned %>%
  mutate(time_bin = as.integer(time_bin)) %>%
  group_by(channel_id, noise_type) %>%
  mutate(trial_num = dense_rank(trial_id)) %>%
  filter(trial_num <= maxTrialNum) %>%
  ungroup() %>%
  select(channel_id, noise_type, trial_num, time_bin, value)


# Step 2: Pad missing trials/time bins with NA
df_complete <- df_trials %>%
  complete(
    channel_id,
    noise_type = conditions,
    trial_num = 1:maxTrialNum,
    time_bin = 1:nTpt,
    fill = list(value = NA_real_)
  )

# Step 3: Prepare array dimensions
channel_ids <- unique(df_complete$channel_id)
N <- length(channel_ids)
S <- S
D <- 1
T <- nTpt
R <- maxTrialNum

# Initialize empty array
tensor_5d <- array(NA_real_, dim = c(N, S, D, T, R),
                   dimnames = list(channel_ids, conditions, NULL, NULL, NULL))

# Step 4: Fill in the 5D tensor
for (n in seq_along(channel_ids)) {
  ch <- channel_ids[n]
  ch_data <- df_complete %>% filter(channel_id == ch)
  
  for (s in seq_along(conditions)) {
    cond <- conditions[s]
    cond_data <- ch_data %>%
      filter(noise_type == cond) %>%
      arrange(trial_num, time_bin)
    
    if (nrow(cond_data) == T * R) {
      mat <- matrix(cond_data$value, nrow = T, ncol = R)
      tensor_5d[n, s, 1, , ] <- mat
    } else {
      # Partial data — pad manually if needed
      browser()
      mat <- matrix(NA_real_, nrow = T, ncol = R)
      filled_data <- cond_data %>%
        filter(trial_num <= R, time_bin <= T) %>%
        mutate(row = time_bin, col = trial_num)
      
      mat[cbind(filled_data$row, filled_data$col)] <- filled_data$value
      tensor_5d[n, s, 1, , ] <- mat
    }
  }
}


meta <-
  lmtime_results_channelwise_timeresolved %>% 
  separate(channel, into = c('type', 'track', 'unit_idx'), remove = FALSE) %>%
  mutate(modality = MODALITY,
         channel_join = if_else(modality=="LFP", channel, glue('{type}_L{track}'))) %>%
  select(channel_join) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, target, region_clean),
            by = join_by(subject_id, run_id==lombard_run_id, channel_join==name)) %>%
  # left_join(SUBJECTS_META %>% select(subject_id, target),
  #           by = join_by(subject_id)) %>%
  # mutate(region_clean = target) %>%
  mutate(channel_id = glue('{subject_id}_{run_id}_{channel}')) 


# d = list(X = tensor_5d, X_meta = meta)
# # install.packages("R.matlab")
# library(R.matlab)
# 
# writeMat(glue("{PATH_ANALYSIS}/data/dpca_tensor_timewarped_withmeta_{t}.mat"), d = d )
# # writeMat(glue("{PATH_ANALYSIS}/data/dpca_tensor_timewarped_{t}.mat"), X = tensor_5d )

library(hdf5r)
t <- format(Sys.time(), "%Y%m%d%H%M")
# # Write tensor and metadata
# h5write(tensor_5d, h5file, "X")
# h5write(meta$channel_id, h5file, "channel_id")
# h5write(meta$subject_id, h5file, "subject_id")
# h5write(meta$run_id, h5file, "run_id")
# h5write(meta$channel, h5file, "channel")


file <- H5File$new(glue("{PATH_ANALYSIS}/data/dpca_tensor_withmeta_{t}.h5"), mode = "w")
file[["X"]] <- tensor_5d
# for (i in seq_along(names(meta))) {
#   file[[i]] <- meta$i
#   
# }
for (col_name in names(meta)) {
  file[[col_name]] <- meta[[col_name]]
}
file$close_all()






# OLD CODE for conversion Step 1: Clean and bin time ---- 
# df_binned <- lmtime_results_channelwise_timeresolved %>%
#   ungroup() %>%
#   # filter(row_number() %in% c(20, 35, 55, 77)) %>%  # testing
#   unnest(data) %>%
#   select(subject_id, run_id, channel, trial_id, noise_type, epoch, loudness_lvl,
#          timewarp_t0, value_residualtime) %>%
#   filter(noise_type %in% c("NOISE", "QUIET"), timewarp_t0>-1 & timewarp_t0<5.3) %>%
#   
#   group_by(subject_id, run_id, channel, trial_id) %>%
#   mutate(time_bin = cut(timewarp_t0, breaks = nTpt, labels = FALSE)) %>%
#   group_by(subject_id, run_id, channel, noise_type, trial_id, time_bin) %>%
#   summarise(value = mean(value_residualtime, na.rm = TRUE), 
#             timewarp_t0 = mean(timewarp_t0, na.rm = TRUE), 
#                         .groups = "drop")
# 
# # Step 2: Assign trial number within each group (up to 40)
# df_trials <- df_binned %>%
#   group_by(subject_id, run_id, channel, noise_type, trial_id) %>%
#   mutate(time_bin = as.integer(time_bin)) %>%
#   group_by(subject_id, run_id, channel, noise_type) %>%
#   mutate(trial_num = dense_rank(trial_id)) %>%  # Assign 1, 2, 3... to trial IDs
#   filter(trial_num <= maxTrialNum) %>%
#   ungroup() %>%
#   select(subject_id, run_id, channel, noise_type, trial_num, time_bin, value = value_residualtime)
# 
# 
# # Step 3: Pad with NAs up to maxTrialNum
# df_complete <- df_trials %>%
#   complete(
#     subject_id, run_id, channel, noise_type,
#     trial_num = 1:maxTrialNum,
#     time_bin = 1:nTpt,
#     fill = list(value = NA_real_)
#   )

# # Unnest and select relevant columns
# df_long <- lmtime_results_channelwise_timeresolved %>%
#   ungroup() %>%
#   filter(row_number() %in% c(20, 35, 55, 77)) %>%  # testing
#   unnest(data) %>%
#   select(subject_id, run_id, channel, trial_id, noise_type, epoch, loudness_lvl,
#          timewarp_t0, value_residualtime) %>%
#   filter(noise_type %in% c("NOISE", "QUIET"), timewarp_t0>-1 & timewarp_t0<5.3)
# 
# df_long <- df_long %>%
#   group_by(subject_id, run_id, channel, trial_id) %>%
#   mutate(time_bin = cut(timewarp_t0, breaks = nTpt, labels = FALSE)) %>%
#   group_by(subject_id, run_id, channel, noise_type, trial_id, time_bin) %>%
#   summarize(value = mean(value_residualtime, na.rm = TRUE), 
#             timewarp_t0 = mean(timewarp_t0, na.rm = TRUE), 
#             .groups = "drop")
# 
# # Ensure exactly 40 trials per condition
# balanced <- df_long %>%
#   group_by(channel, noise_type) %>%
#   filter(n() >= nTpt * maxTrialNum) %>%
#   slice_head(n = nTpt * maxTrialNum)
# 
# # Convert labels to numeric indices
# trial_indexed <- balanced %>%
#   group_by(channel, noise_type) %>%
#   mutate(trial_num = (row_number() - 1) %/% nTr + 1) %>%
#   filter(trial_num <= maxTrialNum)
# channels <- unique(trial_indexed$channel)

# Create a tensor: N × S × D × T × maxTrialNum


# Initialize array
tensor <- array(NA, dim = c(N, S, D, nTpt, maxTrialNum),
                dimnames = list(channels, c("QUIET", "NOISE"), NULL, NULL, NULL))

# Fill in the array
for (n in seq_along(channels)) {
  ch <- channels[n]
  for (s in c("QUIET", "NOISE")) {
    data_slice <- trial_indexed %>%
      filter(channel == ch, noise_type == s) %>%
      arrange(trial_num, time_bin)
    
    values_matrix <- matrix(data_slice$value, nrow = nTpt, ncol = maxTrialNum)
    tensor[n, s, 1, , ] <- values_matrix
  }
}


# -------- PLACE HOLDER, OLD CODE BELOW ----------





















# EDA: plot distribution of beta or gamma in each channel  ------------
ref_str <- "monopolar"
band_str <- "bga"
trials <- 
  trials %>%
  filter(ref_type==ref_str) %>%
  filter(timewin != 'prespeech') %>%
  filter(band == band_str) %>%
  mutate(condition = interaction(noise_type, timewin)) %>%
  mutate(value = (value)) %>% # value = log10(value)
  mutate(valuelog10 = log10(value)) 

  # filter(subject_id %in% c("DM1002", "DM1019")) %>%
  # filter(channel == 'macro_Lc') %>%
  # group_by(channel_id) %>%
  
# calculate statistics in each channel in each condition 
channel_stats <-  
  trials_tmp %>%
  ungroup() %>%
  select(channel_id, condition, value) %>%
  group_by(across(where(negate(is.numeric)))) %>%
  summarize_all_stats() %>% 
  filter(condition=="quiet.basel") %>% 
  group_by(channel_id) %>% 
  select(value_mean, value_stdev)

# # z-score trials with respect to the quiet baseline mean and std
# trials_tmp <-
#   trials_tmp %>%
#     left_join(channel_stats, by=join_by(channel_id)) %>%
#     mutate(value_normed=(value - value_mean)/value_stdev) %>%
#     select(-value_mean, -value_stdev)
# 
# # Calculate t-values between baseline and other conditions for each channel
# t_values <- 
#   trials_tmp %>%
#   # filter(condition != "quiet.basel") %>% # Exclude baseline for comparison
#   group_by(channel_id, run_id, timewin) %>%  # Group by channel and condition
#   do(tidy(t.test(value ~ timewin, 
#                  data = bind_rows(filter(trials_tmp, condition == "quiet.basel"), .),
#                  paired = TRUE)))  # Paired t-test comparing each condition to baseline


# Calculate t-values between baseline and other conditions for each channel, v2 with t-tests or wilcoxon tests
ttest_speech <-
  trials_tmp %>%
  # mutate(logvalue=(value)) %>% # TAKE THE LOG of the gamma power values
  group_by(subject_id, run_id, channel) %>%
  # pivot_wider(names_from = timewin, values_from = value) %>%
  # do(tidy(t.test(.$speech, .$basel, alternative='greater', paired = TRUE))) %>% # wilcoxon
  mutate(timewin = factor(timewin, levels = c("speech", "basel"))) %>%  # Ensures baseline comes first
  do(tidy(t.test(valuelog10 ~ timewin, data=., alternative='greater', paired = FALSE))) %>% # wilcoxon, alternative=less for beta, greater for BGA
  arrange(p.value) %>%
  left_join(electrodes, 
              # select(subject_id, lombard_run_id, name,  mni_x, mni_y, mni_z), 
            by=join_by(subject_id, run_id==lombard_run_id, channel==name))

ttest_speech %>% # print channels that are in STN-capsule and have strong gamma activation
  filter(region_clean=='STN-capsule') %>%
  arrange(desc(statistic)) %>%
  select(1:7)
#   subject_id run_id channel  estimate estimate1 estimate2 statistic
# <chr>       <dbl> <chr>       <dbl>     <dbl>     <dbl>     <dbl>
#   1 DM1033          3 macro_Ll   0.0472    -0.731    -0.778      7.58
# 2 DM1033          3 macro_Lc   0.0399    -0.736    -0.776      6.00
# 3 DM1027          3 macro_Lp   0.0176    -0.795    -0.812      3.28
# 4 DM1029          2 macro_Lp   0.0286    -0.637    -0.666      3.27
# 5 DM1029          3 macro_Lc   0.0239    -0.666    -0.690      2.85
# 6 DM1029          3 macro_Lp   0.0272    -0.697    -0.724      2.82
# 7 DM1003          3 macro_Ll   0.0186    -0.862    -0.880      2.46

fname = glue('A04b_speech-vs-basel-ttest_band-{band_str}_ch-type-{ref_str}')
write_tsv(ttest_speech, glue('{PATH_ANALYSIS}/{fname}.tsv'))

# Join t-values back to the data for plotting or further analysis
trials_tmp <- trials_tmp %>%
  left_join(t_values, by = c("channel_id", "condition"))


trials_tmp %>% # filter(subject_id=="DM1033") %>% filter()
  ggplot(aes(x=channel, 
             y=value_normed, 
             color=noise_type,
             shape=timewin, 
             group=interaction(channel, condition))) +
  geom_jitter(position=position_dodge(width=0.75,preserve="total")) + 
  geom_boxplot(position=position_dodge(width=0.75,preserve="total"),outlier.alpha=0, show.legend = FALSE) + 
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) +
  # scale_alpha_manual(values = c("speech" = 0.9, "basel" = 0.4))+ 
  scale_shape_manual(values = c("speech" = 16, "basel" = 1)) +
  facet_wrap(target~subject_id+run_id, scales='free_x', ncol=4)+
  theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

fname = glue('A04_01_all-timewins_band-{band_str}_ch-type-{ref_str}_log10')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=10, height=24, units='in')
write_tsv(trials_tmp, glue('{PATH_FIG}/{fname}.tsv'))

trials_tmp %>%
  # filter(timewin=='speech') %>%
  
  ggplot() +
  aes(x=intensity_praat, 
      y=value_normed, 
      color=noise_type, 
      shape=timewin
      ) + 
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_point(alpha=0.5) + 
  facet_wrap(~channel_id, scales='free', ncol=5) +
  scale_shape_manual(values = c("speech" = 16, "basel" = 1)) +
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) + 
  theme_minimal()
  # geom_jitter(position=position_dodge(width=0.75,preserve="total")) + 
  # geom_boxplot(position=position_dodge(width=0.75,preserve="total"),outlier.alpha=0) + 
  # # scale_color_manual(values = color_map_targets$hex, breaks=color_map_targets$name) +
  # facet_wrap(target~subject_id+run_id, scales='free_x', ncol=4)+
  # theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))

fname = glue('A04_02_bandpw-vs-volume-by-channel_band-{band_str}_ch-type-{ref_str}_freex')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=10, height=24, units='in')
# write_tsv(trials_tmp, glue('{PATH_FIG}/{fname}.tsv'))






## Filter and prepare table for statistics ----
ref_str <- "monopolar"
band_str <- "beta"
target_str <- "GPi-STN"
elec_type_str <- "ecog-macro"
is_data_logged <- TRUE

trials_tmp <-
  trials %>% 
  filter(ref_type==ref_str) %>%
  filter(timewin != 'prespeech') %>%
  filter(band == band_str) %>%
  filter(str_detect(target_str, target)) %>%  # GPi, STN, or GPi-STN
  filter(str_detect(channel, "ecog") | str_detect(channel, "macro")) %>%
  # filter(str_detect(elec_type_str, channel)) %>%  # GPi, STN, or GPi-STN
  {if (is_data_logged) mutate(., value = log10(value)) else .}

fname = glue('anova-by-channel__band-{band_str}_electype-{elec_type_str}_islogged-{is_data_logged}')

## Run-channel-wise ANOVAs with on the 'design matrix': baseline vs speech, quiet vs lombard ----
fit_anova <- function(df) {
  # model <- aov(value ~ 1 + noise_type*timewin + time_withinrun+time_withinrun2, data = df)
  model <- lm(value ~ 1 + noise_type*timewin + time_withinrun+time_withinrun2, data = df)
  summary <- tidy(model)
  significant <- summary %>% filter(p.value < 0.05)
  df$residuals <- residuals(model)
  list(model = model, 
       data = df, 
       stat_noise_type = summary$statistic[summary$term=="noise_type"], # %>%filter(term=noise_type) %>% pull(statistic))
       stat_timewin = summary$statistic[summary$term=="timewin"])
}

    





# Original ANOVA results
anova_results_channelwise <- trials_tmp %>%
  mutate(time_withinrun = as.numeric(time_withinrun),
         time_withinrun2 = time_withinrun^2) %>%
  group_by(subject_id, run_id, channel) %>%
  nest() %>%
  mutate(model_and_data = map(data, fit_anova)) %>%
  mutate(model = map(model_and_data, "model"),
         data = map(model_and_data, "data"),
         stat_noise_type = map_dbl(model_and_data, "stat_noise_type"), 
         stat_timewin = map_dbl(model_and_data, "stat_timewin")) %>%
  select(-model_and_data)

# Permutation test
set.seed(123)
num_permutations <- 100
permute_type <- "permute_trials" # permute_trials, permute_blocks, circshift_trials, permute_residuals
var_str <- "stat_timewin" # predictor of interest... "stat_noise_type" or "stat_timewin" or "stat_intensity_praat"
permuted_anova_results_channelwise <- 
  replicate(num_permutations, {
    permuted_results <- 
      trials_tmp %>%
      mutate(time_withinrun = as.numeric(time_withinrun),
             time_withinrun2 = time_withinrun^2) %>%
      group_by(subject_id, run_id, channel) %>%
      nest() %>%
      mutate(model_and_data = map(data, \(x) permute_and_fit(x, permute_type, var_str))) %>%
      mutate(stat_noise_type = map_dbl(model_and_data, "stat_noise_type"), 
             stat_timewin = map_dbl(model_and_data, "stat_timewin")) %>%
      select(-data, -model_and_data)
    # sum(permuted_results)
  }, simplify = FALSE)
permuted_anova_results_channelwise <- bind_rows(permuted_anova_results_channelwise) 


# ¡  nest(data_permute = c(stat_noise_type, stat_timewin))
 # unnest(data)
# started 4:06 pm 
# permuted_anova_results_channelwise 

get_percentile <- function(distr, val) {
  fn<-ecdf(distr)
  fn(val)
}
anova_results_channelwise <-
  anova_results_channelwise %>%
  select(-any_of("data_permute")) %>% # make sure the left join gets the updated values
  left_join(permuted_anova_results_channelwise %>% nest(data_permute = c(stat_noise_type, stat_timewin))) %>%
  # rowwise() %>% 
  # mutate(stat_noise_type_pval = mean(abs(data_permute$stat_noise_type) >= abs(stat_noise_type))) %>%
  # mutate(stat_timewin_pval = mean(abs(data_permute$stat_timewin) >= abs(stat_timewin))) 
  mutate(stat_noise_type_prctl = map2_dbl(data_permute, stat_noise_type, \(x, y) get_percentile(x$stat_noise_type, y))) %>%
  mutate(stat_timewin_prctl = map2_dbl(data_permute, stat_timewin, \(x, y) get_percentile(x$stat_timewin, y))) %>%
  mutate(sig_noise_type=stat_noise_type_prctl>0.95,
         sig_timewin=stat_timewin_prctl>0.95) %>%
  mutate(sig=case_when(sig_noise_type&&sig_timewin ~ 'both', 
                       sig_noise_type ~ 'noise_type', 
                       sig_timewin ~ 'timewin',
                       TRUE~'none')) 

write_tsv(anova_results_channelwise, file=glue('{PATH_FIG}/A04_03_{fname}_fstat_permute-{permute_type}_var-{var_str}_channelwise.tsv'))

# print summary table
anova_results_channelwise %>% 
  ungroup() %>%
  summarise(n_sig_noise_type = sum(sig_noise_type), 
            n_sig_timewin = sum(sig_timewin), 
            n = n())



# Plot CHANNELWISE null distributions and actual values
n_obs <- nrow(anova_results_channelwise) 
plot_top_n <- 20 # number of channels to plot,
anova_results_channelwise_subset <- 
  anova_results_channelwise %>% 
  ungroup() %>% 
  arrange(desc(across(glue(var_str, 'prctl', .sep = '_')))) %>%
  slice_head(n=plot_top_n) %>%
  mutate(channel_id = paste(subject_id, channel, run_id, sep = "_")) %>%
  select(-all_of(vars_channel_id)) %>%
  ungroup() %>%
  mutate(facet_order = fct_inorder(channel_id)) %>%
  group_by(facet_order)
ggplot() + 
  geom_histogram(data = anova_results_channelwise_subset %>% select(data_permute) %>% unnest(), aes(x = !!sym(var_str), y = ..density..), 
                 fill = "grey") + 
  # geom+(data = anova_results_channelwise_subset, 
  #               aes(y = 0, xmin = stat_noise_type_q05, xmax = stat_noise_type_q95), height = 0.1, color = "red") +
  geom_vline(data = anova_results_channelwise_subset, aes(xintercept=!!sym(var_str)) , color="red") + 
  facet_wrap(facet_order~., scales = "free") + 
  theme_minimal() +
  theme(axis.text.y = element_blank())
  
  # coord_cartesian(xlim = c(-40, 40), ylim = c(0, 0.2))
ggsave(file=glue('{PATH_FIG}/A04_03_{fname}_fstat-permute-{permute_type}_var-{var_str}_channelwise.pdf'), dpi=300, width=10, height=8,  units='in') # 

# scatterplot of timewin vs noise_type p-values based on f-statistic permutations
anova_results_channelwise %>%
  # ungroup() %>%
  # filter(term %in% c("noise_type", "timewin")) %>%

  ggplot(aes(x=stat_noise_type_prctl, y=stat_timewin_prctl, fill=sig)) +
  # ggplot(aes(x=stat_noise_type_prctl, y=stat_timewin_prctl)) +
  geom_point(shape=21, color='black', alpha=0.7, size=5) + 
  # stat_cor(method = "spearman") +
  # scale_shape_manual(values=c(1, 16)) +
  scale_fill_manual(values=c('blue', 'black','lightgrey',  'red')) +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title=glue('All channels: permutation ({permute_type}) test percentile ')) + 
  coord_fixed(ratio = 1) +  
  theme_bw()
ggsave(file=glue('{PATH_FIG}/{fname}_fstat_permute-{permute_type}_var-{var_str}_channelwise-scatter-pval.pdf'), dpi=300, width=6, height=5, units='in')


# is the overlap between noise_type and timewin modulations significant? 
n_sig_noise_type <- sum(anova_results_channelwise$sig_noise_type)
n_sig_timewin <- sum(anova_results_channelwise$sig_timewin)

# Permutation test
set.seed(123)
num_permutations <- 1000
n_overlap <- 
  replicate(num_permutations, {
    anova_results_channelwise %>%
      filter(sig_noise_type | sig_timewin) %>%
      ungroup() %>%
      slice_sample(n=n_sig_noise_type+n_sig_timewin, replace = TRUE) %>%
      summarise(n_both = sum(sig_noise_type & sig_timewin)) %>%
      pull(n_both)
    # sum(permuted_results)
  }, simplify = TRUE) 
n_overlap <- tibble(n_overlap=n_overlap) %>%
  group_by(n_overlap) %>% summarise(count = n())

ggplot(n_overlap, aes(x=n_overlap, y=count)) + 
  geom_col() + 
  scale_x_continuous(breaks = seq(0, 10, by = 1)) + 
  geom_vline(xintercept=(anova_results_channelwise %>%
                     filter(sig_noise_type & sig_timewin) %>% nrow()), 
             color='red', size=2) + 
  theme_bw() + 
  labs(title=glue('Noise and speech modulations: permutation test'))
  
  


n <- 85
p_both <- 0.04982701
x <- seq(0, 10, by = 0.01)
df <- tibble(x=x, d=dbinom(x, n, p_both))
idx <- which(df$d>0.95)[1]
ggplot(df, aes(x=x, y=d, width=0.9)) + 
  geom_col() + 
  # geom_hline(yintercept = 0.95, color='grey', size=2) + 
  geom_vline(xintercept = df$x[idx], color='grey', size=2) +
  geom_vline(xintercept = 4, color='red', size=2) + 
  scale_x_continuous(breaks = seq(0, 10, by = 1)) + 
  theme_bw()

# Can also establish independence between speech and noise modulations with a fisher exact 
# detect whether we observe independence between rows and columns. 
modulations <- matrix(c(4, 6, 32, 43), 
                      nrow = 2,
                      dimnames = list(Speech=c("yes", "no"),
                                      Noise=c("yes", "no")))
fisher.test(modulations, alternative = "greater")
# >> 	Fisher's Exact Test for Count Data
# 
# data:  modulations
# p-value = 0.6869
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.2201788       Inf
# sample estimates:
# odds ratio 
#  0.8970075
# .... so we FAIL TO REJECT the null hypothesis


# plot individual channels over trials ---- 
plot_single_channel_over_trials <- function(df, mdl) {

df$predicted <- predict(mdl, df)
df$residual <- residuals(mdl)

df <-
  df %>% 
  rename(pwr=value) %>% 
  pivot_longer(cols = c("pwr", 'predicted'), values_to = c("value"), names_to = "series")

shades <- df %>%
  mutate(block_id = floor((trial_id-1)/10)) %>%
  group_by(noise_type, block_id) %>%
  summarise(start = min(trial_id), end = max(trial_id)) 

chid = df$channel_id[[1]]
ggplot(df, aes(x=trial_id, y=value, color=noise_type, shape=interaction(timewin, series))) + 
  geom_point() + 
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) + 
  scale_shape_manual(values=c(2, 17, 1, 16)) + 
  geom_rect(data = shades, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = noise_type), 
            alpha = 0.2, inherit.aes = FALSE) +
  scale_fill_manual(values = color_map_conditions$hex, breaks = color_map_conditions$name) +
  labs(title=glue('{chid}')) + 
  theme_bw() 
ggsave(file=glue('{PATH_FIG}/channelwise/{chid}_real-vs-predicted.pdf'), dpi=300, width=8, height=5, units='in')

df
}

anova_results_channelwise %>%
  mutate(result = pmap(list(data, model), plot_single_channel_over_trials))



# Compare different types of permutations: scatterplot of fstat permute1 vs permute2 
fname_base <- 
var_str <- "stat_timewin" # "stat_noise_type" or "stat_timewin"

permute_type1 <- "permute_trials" # permute_trials, permute_blocks, circshift_trials, permute_residuals
anova_results_channelwise_1 <- read_tsv(glue("{PATH_FIG}/anova-by-channel__band-bga_electype-ecog-macro_islogged-TRUE_fstat_permute-{permute_type1}_var-{var_str}_channelwise.tsv")) 

permute_type2 <- "permute_residuals" # permute_trials, permute_blocks, circshift_trials, permute_residuals
anova_results_channelwise_2 <- read_tsv(glue("{PATH_FIG}/anova-by-channel__band-bga_electype-ecog-macro_islogged-TRUE_fstat_permute-{permute_type2}_var-{var_str}_channelwise.tsv")) 

var_str <- "stat_timewin_prctl" # "stat_noise_type_prctl" or "stat_timewin_prctl"
anova_results_channelwise_all <- 
  bind_rows(anova_results_channelwise_1 %>% mutate(permute_type=permute_type1), 
            anova_results_channelwise_2 %>% mutate(permute_type=permute_type2)) %>%
  group_by(across(vars_channel_id)) %>%
  pivot_wider(id_cols = all_of(vars_channel_id), values_from = !!sym(var_str), names_from = permute_type, names_prefix = glue("{var_str}_"))
  
anova_results_channelwise_all %>%
  ggplot(aes(x=!!sym(glue("{var_str}_{permute_type1}")), 
             y=!!sym(glue("{var_str}_{permute_type2}")), 
             color=subject_id)) + 
  geom_abline() +
  geom_point() + 
  theme_bw()
ggsave(file=glue('{PATH_FIG}/{fname}_fstat_permute-{permute_type1}-vs-{permute_type2}_var-{var_str}_channelwise-scatter-pval.pdf'), dpi=300, width=6, height=4, units='in')







# Plot histograms of pooled data, across all channels
var_str <- "stat_timewin" # "stat_noise_type" or "stat_timewin"
perc_95 <- quantile(permuted_anova_results_channelwise[[var_str]], 0.95)
n_obs <- nrow(anova_results_channelwise)
ggplot() +
  geom_histogram(data = permuted_anova_results_channelwise, aes(x = !!sym(var_str), y = (..count../sum(..count..))*n_obs), 
                 fill = "grey", alpha = 0.9, binwidth = 1, center = 0.5) + 
  geom_dotplot(data = anova_results_channelwise, aes(x = !!sym(var_str), y = ..count..), 
               binwidth = 1, center = 0.5,
               fill = "blue", alpha = 0.3, method = "histodot", ) + 
  geom_vline(xintercept = perc_95, color = "red", linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("Actual" = "blue", "Permuted" = "grey")) +
  theme_minimal() +
  labs(title = "Histograms of Actual and Permuted Data",
       x = "f-stat",
       y = "Count",
       fill = "Type") + 
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 35))
ggsave(file=glue('{PATH_FIG}/{fname}_fstat-permute-circshift-var-{var_str}.pdf'), dpi=300, width=5, height=5,  units='in') # 







# get the residuals by unnesting, plot 
anova_results_channelwise %>%
  unnest(cols=data) %>%
  ungroup() %>%
  # group_by(subject_id, run_id, channel) %>%

  ggplot(aes(sample=resid, color=channel)) +
  geom_point(stat = "qq", alpha = 0.5) +  # Use geom_point with stat="qq" and set alpha
  stat_qq_line() +
  labs(x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  facet_wrap(target~subject_id+run_id, scales='free', ncol=4)
ggsave(file=glue('{PATH_FIG}/{fname}_residuals.pdf'), dpi=300, width=10, height=24, units='in')


# get the ANOVA model and summarize with tidy
anova_results <-
  anova_results_channelwise %>%
  mutate(model_tidy = map(model, tidy)) %>%
  unnest(cols=model_tidy)
write_tsv(anova_results, glue('{PATH_FIG}/{fname}.tsv'))

aa <- anova_results %>% filter(p.value<0.05 & term=="noise_type") %>% arrange(desc(statistic))

# Show correlation between f-statistics in the terms
anova_results_wide <-
  anova_results %>%
  # ungroup() %>%
  filter(term %in% c("noise_type", "timewin")) %>%
  mutate(sig=p.value<0.05) %>%
  select(term, statistic, sig) %>%
  pivot_wider(names_from=term, values_from=c("statistic", "sig")) %>%
  mutate(sig=case_when(sig_noise_type&&sig_timewin ~ 'both', 
                       sig_noise_type ~ 'noise_type', 
                       sig_timewin ~ 'timewin',
                       TRUE~'none'))


 anova_results_wide %>%
  ggplot(aes(x=statistic_noise_type, y=statistic_timewin, fill=sig)) +
  # ggplot(aes(x=statistic_noise_type, y=statistic_timewin)) + 
  geom_point(shape=21, color='black', alpha=0.7) + 
  # stat_cor(method = "spearman") +
  # scale_shape_manual(values=c(1, 16)) +
  scale_fill_manual(values=c('blue', 'black','grey',  'red')) +
  xlim(0, 40) +
  ylim(0, 40) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title='f-stat') + 
  coord_fixed(ratio = 1) +  
  theme_bw()
ggsave(file=glue('{PATH_FIG}/{fname}_fstat-timewin-vs-noise-type.pdf'), dpi=300, width=4, height=3, units='in')

# inspect individual data points 
aa <-
  anova_results_wide %>%
  group_by(subject_id, run_id) %>%
  arrange(desc(statistic_noise_type), desc(statistic_timewin), .by_group = TRUE) %>% 
  filter(row_number() %in% c(1, 2, 3)) %>%
  left_join(electrodes %>% select(subject_id, lombard_run_id, name, area), by=join_by(subject_id, run_id==lombard_run_id, channel==name)) %>%
  filter(area %in% c('PMC', 'SMC'))



anova_summ <-
  anova_results %>%
  filter(term!='Residuals') %>%
  group_by(term) %>%
  mutate(n_channels = n()) %>%
  filter(p.value<0.05) %>%
  # group_by(term) %>% mutate(n_subj=n()) %>%
  arrange(term, subject_id, (p.value)) %>%
  # group_by(term, subject_id, run_id, channel) %>% add_count(name="n_channel") %>%
  group_by(term) %>%
  summarise(perc_sig = n()/mean(n_channels)*100,
            n_chan=n(), 
            n=mean(n_channels),
            n_subj=n_distinct(subject_id))
anova_summ 

fname = glue('{fname}_summary')
write_tsv(anova_summ, glue('{PATH_FIG}/{fname}.tsv'))
  

## Try group-level LME on 'design matrix': baseline vs speech, quiet vs lombard ---- ----
speech_modulated_channels <-
  anova_results %>% 
  filter(term=='timewin' & p.value < 0.05) %>%
  group_by_all()





model <- 
  trials_tmp %>% 
  semi_join(ttest_speech, by=c("subject_id", "run_id", "channel")) %>% # filter for only channels that were sig different during speech
  # filter(value>-4 & value<4)%>%
  lmer(value ~ 1 + noise_type * timewin + (1 | channel_id), data = .)

summary(model)

confint(model, method="Wald")

# Residuals
resid <- residuals(model)

# Create diagnostic plots
plot_residuals <- function(residuals) {
  ggplot() +
    geom_histogram(aes(x = residuals, y = ..density..), bins = 20, fill = "skyblue", color = "black") +
    geom_density(aes(x = residuals), color = "blue") +
    labs(x = "Residuals", y = "Density") +
    theme_minimal()
}

# Plot residuals
plot_residuals(resid)



## Run-channel-wise linear model/correlation on bandpower vs speech volume ----
lm_pw_vs_vol_results <-
  trials_tmp %>%
  group_by(subject_id, run_id, channel) %>%
  mutate(time_withinrun=as.numeric(time_withinrun), 
         time_withinrun2=time_withinrun**2) %>% # calculate square of trial number 
  mutate(logvalue=(value)) %>% # TAKE THE LOG of the gamma power values
  do(tidy(lm(logvalue ~ 1 + intensity_praat, data = .)))  


lm_pw_vs_vol_summ <-
  lm_pw_vs_vol_results %>%
  summarise(n_sig = sum(p.value<0.05),
            n = n(),
            perc_sig = n_sig/n,
            n_sig_pos = sum(estimate<0 & p.value<0.05),
            n_sig_neg = sum(estimate>0 & p.value<0.05),
            n_subj= n_distinct(subject_id))
lm_pw_vs_vol_summ

lm_pw_vs_vol_results %>%
  filter(term=='intensity_praat') %>% # intensity_praat, (Intercept)
  # group_by(term) %>%
  arrange(term, p.value) %>%
  # filter(p.value<0.05 & estimate<0) %>%
  # group_by(term) %>% mutate(n_subj=n()) %>%
  # arrange(term, subject_id) %>%
  # arrange(term, subject_id, (p.value)) %>%
  # group_by(term, subject_id, run_id, channel) %>% add_count(name="n_channel") %>%
  group_by(term) 

 

fname = glue('pw-vs-vol-channelwise_chs-{ref_str}_band-{band_str}_target-{target_str}_summary')
write_tsv(lm_pw_vs_vol_summ, glue('{PATH_FIG}/{fname}.tsv'))




## Try group-level LME on power vs volume ---- ----
model <- 
  trials_tmp %>% 
  # semi_join(ttest_speech, by=c("subject_id", "run_id", "channel")) %>% # filter for only channels that were sig different during speech
  # filter(value>-4 & value<4)%>%
  lmer(value_normed ~ 1 + intensity_praat + (1 | channel_id), data = .)

(summary(model))

# -- FROM CHATGPT --
fixed_effects <- summary(model)$coefficients
fixed_effects_df <- as.data.frame(fixed_effects)
fixed_effects_df$Effect <- rownames(fixed_effects_df)
fixed_effects_df <- fixed_effects_df %>%
  select(Effect, Estimate, `Std. Error`, df, `t value`, `Pr(>|t|)`)

# Extract random effects summary
random_effects <- VarCorr(model)
random_effects_df <- as.data.frame(random_effects)
random_effects_df <- random_effects_df %>%
  select(grp, var1, vcov) %>%
  mutate(Std.Dev. = sqrt(vcov))

# Combine fixed and random effects into one data frame
summary_df <- bind_rows(
  fixed_effects_df %>% mutate(Type = "Fixed"),
  random_effects_df %>% mutate(Type = "Random")
)
# -- FROM CHATGPT --



fname = glue('pw-vs-vol-lme_chs-{ref_str}_band-{band_str}_target-{target_str}')
write_tsv(summary_df, glue('{PATH_FIG}/{fname}.tsv'))

trials_tmp %>%
  ungroup() %>%
  summarise(n_subj = n_distinct(subject_id), 
            n_runs = n_distinct(subject_id, run_id), 
            n_chan = n_distinct(subject_id, run_id, channel))





# PLOT SUMMARY statistics, each point represents a subject
resp_plot_summ <- 
  resp_plot %>%
  group_by(subject_id, noise_type) %>%
  summarize(resp_range_mean = mean(resp_value), 
            resp_range_median = median(resp_value))
resp_plot_summ %>%
  ggplot() + 
  aes(x=noise_type, y=resp_range_mean, color=noise_type, group=subject_id) +
  geom_line() +
  geom_point(alpha=0.7) + 
  scale_color_manual(values = color_map_conditions$hex, breaks=color_map_conditions$name) + 
  theme_minimal()

fname = glue('A04_05b_resp-all-subj-summary-range-absolute')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=3, height=4, units='in')


mdl <- lmer(range ~ 1 + noise_type + (1|subject_id), resp)
confint(mdl, method="Wald")



## EXTRA CODE ----













# ## find some outliers
# aa <-
#   channels_tmp %>%
#   group_by(channel_id) %>%
#   arrange(desc(value_normed), .by_group = TRUE) %>%
#   slice_head(n = 3) %>%
#   filter(subject_id %in% c("DM1003", "DM1019", "DM1033", "DM1015")) %>%
#   filter(subject_id %in% c("DM1033"))  %>%
#   ungroup() %>%
#   group_by(subject_id, run_id) %>%
#   arrange(desc(value_normed), .by_group = TRUE)





## write file out for plotting electrodes ---- 
channels_out <- 
  channels %>%
  separate(channel, c("channel_first", NA), '-', remove=FALSE) %>%
  filter(channel_first!='macro_L_CA') %>%
  left_join(electrodes, by=c("channel_first"="name", "subject_id", "run_id"="lombard_run_id"))
fname = glue('A04b_channels_summary_with_electrode_info.tsv')
write_tsv(channels_out, glue('{fname}'))


# Save a table for best channel from each run 
max_bga_speech <- 
  channels %>% 
  group_by(subject_id, run_id) %>% 
  filter(bga_speech_mean == max(bga_speech_mean)) %>%
  distinct(., .by_group=TRUE, .keep_all = TRUE) # it can happen that two rows have the exact same max()
fname = glue('A04b_channel_summary_max_bga_speech.tsv')
write_tsv(max_bga_speech, glue('{fname}'))


# Summarize data for a power analysis
tmp <-
  channels %>%
  filter(ref_type=='monopolar', target=='GPi') %>%
  group_by(subject_id, run_id) %>%
  select(channel, ref_type, bga_speech_mean)
  

ggplot(tmp) + 
  aes(x = subject_id, y = bga_speech_mean) + 
  geom_jitter() 

tmp <- 
  tmp %>%
  ungroup() %>%
  select(bga_speech_mean) %>%
  summarize_all_stats() %>%
  pivot_longer(cols=everything())

median(tmp$bga_speech_mean) / std(tmp$bga_speech_mean_stdev)


tmp %>%
group_by(subject_id) %>%
select(bga_speech_mean) %>%
summarize_all_stats() %>%
summarise(mean_elecs_per_subj = mean(bga_speech_mean_n))





# ggplot(channels_tmp, aes(x=condition, y=value_normed, color=condition)) +
#   geom_jitter() + 
#   # scale_color_manual(values = color_map_targets$hex, breaks=color_map_targets$name) +
#   facet_grid(~interaction(subject_id, run_id, channel), scales='free')
# 
# ggplot(channels_tmp, aes(x=run_id, y=value_normed, color=condition, 
#                          group=interaction(run_id, condition))) +
#   geom_boxplot(position=position_dodge(width=0.75,preserve="total"),outlier.alpha=0) + 
#   geom_jitter(position=position_dodge(width=0.75,preserve="total")) + 
#   # scale_color_manual(values = color_map_targets$hex, breaks=color_map_targets$name) +
#   facet_wrap(~subject_id+channel, scales='free', ncol=8)+
#   theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0))
# ggave()



# Plot STN vs GPi beta desynchornization  ------------
channels_tmp <- 
  channels %>%
  filter(ref_type=='monopolar') %>%
  select(where(negate(is.numeric)), noise_type, trial_id, 
         contains('mean')) %>%
  pivot_longer(cols = beta_basel_mean:last_col(), 
               names_to=c('band', 'timewin', 'valtype'), 
               names_sep='_', 
               values_to = "value") %>%
  filter(timewin != 'prespeech') %>%
  rename(mean_beta=value)
ggplot(channels_tmp, aes(x=target, y=mean_beta, color=target)) +
  geom_jitter() + 
  stat_compare_means() +
  scale_color_manual(values = color_map_targets$hex, breaks=color_map_targets$name) +
  facet_wrap(vars(timewin), ncol=2)
fname = glue('beta-desynch-prespeech-speech-monopolar')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=700*2.5, height=600*2.5, units='px')
write_tsv(channels_tmp, glue('{PATH_FIG}/{fname}.tsv'))
# ggplotly(p)






# Calculated summary table for percent signifiant sites  ------------

channels_perc_sig  <- 
  channels %>%
  filter(ref_type == 'monopolar') %>%
  # group_by(ref_type, subject_id) %>%
  # group_by(ref_type) %>%
  group_by(target) %>%
  # ungroup() %>%
  summarize_all_stats() %>%
  select(where(negate(is.numeric)),
        contains('sig_mean'), 
        contains('sig_sum'), 
        contains('sig_n')) %>%
  pivot_longer(cols = where(is.numeric), 
               names_to=c('band', 'timewin', 'sig', 'statistic'), 
               names_sep='_', 
               values_to = "value") %>%
  pivot_wider( names_from=statistic, 
               values_from=value)  %>%
  mutate(lab = str_c(as.character(sum), '/', as.character(n))) %>%
  filter(timewin != 'basel') %>%
  arrange(band, timewin)
fname = glue('channel-summary-summary-monopolar')
# ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300)
write_tsv(channels_perc_sig, glue('{PATH_FIG}/{fname}.tsv'))


ggplot(channels_perc_sig, aes(x=target, y=mean, fill=target)) + 
  # geom_col(position = 'dodge') +
  geom_bar(position = 'dodge', width = 0.75, colour="black", stat="identity") + 
  # geom_point() + 
  # geom_line() +
  scale_fill_manual(values = color_map_targets$hex, breaks=color_map_targets$name) +
  facet_grid(vars(band), vars(timewin)) + 
  geom_text(aes(label = lab), vjust = -0.5) +
  ylim(0, 0.8) + 
  ylab('% Significant')
fname = glue('signif-by-target')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300, width=650*2.5, height=500*2.5, units='px')
write_tsv(channels_perc_sig, glue('{PATH_FIG}/{fname}.tsv'))





ggplot(channels_perc_sig, aes(x=subject_id, y=sig_mean, fill=ref_type)) + 
  # geom_col(position = 'dodge') +
  geom_bar(position = 'dodge', width = 0.75,colour="black", stat="identity") + 
  # geom_point() + 
  # geom_line() + 
  facet_grid(vars(band), vars(timewin)) + 
  ylab('% Significant')
fname = glue('signif-monopolar-vs-bipolar-by-subject')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300)
write_tsv(channels_perc_sig, glue('{PATH_FIG}/{fname}.tsv'))



ggplot(channels_perc_sig, aes(x=ref_type, y=sig_mean, group=subject_id, color=ref_type)) + 
  # geom_col(position = 'dodge') +
  # geom_bar(position = 'dodge',width = 0.75,  colour="black", stat="identity") + 
  geom_point() +
  geom_line() +
  facet_grid(vars(band), vars(timewin)) + 
  ylab('% Significant')
fname = glue('signif-monopolar-vs-bipolar-by-subject-collapsed')
ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300)
write_tsv(channels_perc_sig, glue('{PATH_FIG}/{fname}.tsv'))


# Display ------------------------













# 
# # load all trials
# annot_name = 'trials'
# tmp <-
#   tibble(subject=paste0('DM',subject_ids)) %>%
#   mutate(path_annot=glue("{PATH_DATA}/sub-{subject}/annot/sub-{subject}_ses-{SESSION}_task-{TASK}_annot-{annot_name}.tsv")) %>%
#   mutate(exists_annot=file.exists(path_annot)) %>%
#   filter(exists_annot)
# 
# files <- tmp$path_annot
# names(files) <- tmp$subject
# trials <- files %>% 
#   map_df(read_tsv, .id = "subject") %>% 
#   mutate(subject_id = subject) %>%
#   mutate(audio_on = audio_onset) %>%
#   mutate(audio_off = audio_onset + audio_duration)
# 
# 
# # load all sentences
# annot_name = 'produced-sentences-acoustics'
# sentences <-
#   tibble(subject=paste0('DM',subject_ids)) %>%
#   mutate(path_annot=glue("{PATH_DATA}/sub-{subject}/annot/sub-{subject}_ses-{SESSION}_task-{TASK}_annot-{annot_name}.tsv")) %>%
#   mutate(exists_annot=file.exists(path_annot)) %>%
#   filter(exists_annot) %>% 
#   pull(path_annot) %>% 
#   map_df(read_tsv, .id = "subject") %>% 
#   mutate(subject_id = sub) %>%
#   mutate(sentence_on = onset) %>%
#   mutate(sentence_off = onset + duration) %>%
#   mutate(noise_type = as.factor(noise_type)) %>%
#   group_by(subject_id, run_id, trial_id, noise_type)
# sentences$noise_type <- recode(sentences$noise_type, '0' = "QUIET",  '1' = "LMBRD")
# 
# # load all vowels
# annot_name = 'produced-vowels-acoustics'
# vowels <-
#   tibble(subject=paste0('DM',subject_ids)) %>%
#   mutate(path_annot=glue("{PATH_DATA}/sub-{subject}/annot/sub-{subject}_ses-{SESSION}_task-{TASK}_annot-{annot_name}.tsv")) %>%
#   mutate(exists_annot=file.exists(path_annot)) %>%
#   filter(exists_annot) %>% 
#   pull(path_annot) %>% 
#   map_df(read_tsv, .id = "subject") %>% 
#   mutate(subject_id = sub) %>%
#   mutate(sentence_on = onset) %>%
#   mutate(sentence_off = onset + duration) %>%
#   mutate(noise_type = as.factor(noise_type)) 
# vowels$noise_type <- recode(vowels$noise_type, '0' = "QUIET",  '1' = "LMBRD")
# 
# 
# voi <- c('IY', 'UW', 'AA', 'AE')
# # fcr = (mean(u.F2) + mean(a.F2) + mean(i.F1) + mean(u.F1)) / (mean(i.F2) + mean(a.F1));
# 
# # % create FCR measure for every sentence
# FCR <- vowels %>%
#   mutate(phoneme = str_replace_all(phoneme, "[:digit:]", "")) %>%
#   group_by(subject_id, run_id, trial_id, phoneme, noise_type) %>%
#   summarise(F1 = mean(F1), F2 = mean(F2)) %>%
#   filter(phoneme %in% voi) %>%
#   pivot_wider(names_from = phoneme, values_from = c("F1", "F2")) %>%
#   mutate(FCR = (F2_UW + F2_AA + F1_IY + F1_UW)/(F2_IY + F1_AA)) %>%
#   group_by(subject_id, run_id, trial_id, noise_type)
# 
# sentences <- vowels %>%
#   group_by(subject_id, run_id, trial_id, noise_type) %>%
#   summarise_at(vars(intensity_rms:intensity_praat, duration), ~ mean(.x, na.rm = TRUE))
# 
# sentences <- left_join(sentences, FCR)
# sentences <- left_join(sentences, trials %>% select(subject_id, run_id, trial_id, block_id))
# 
# 
# 
# 
# # plot distributions for each individual 
# color_by <- sym("intensity_rms") # F0, FCR, intensity_rms, intensity_praat, duration
# ggplot(sentences, aes(x=noise_type, y=!!color_by, fill=noise_type, colours())) + 
#   geom_jitter() + 
#   geom_boxplot(width=0.5) + scale_fill_manual(values=color_map) + 
#   stat_compare_means() +
#   facet_wrap(~subject_id, scale="free_y", ncol = 5) + 
#   ggtitle(glue('{rlang::as_string(color_by)}-sentences-per-subject'))
# 
# fname = glue('{rlang::as_string(color_by)}-sentences-per-subject')
# # ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), width=10, height=6, dpi=300)
# # ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300)
# dev.print(pdf, glue('{PATH_FIG}/{fname}.pdf'))
# 
# 
# # single plot with summary stats for each patient in 2-d scatter
# sentences %>%
#   group_by(subject_id, noise_type) %>%
#   summarize(m = mean(!!color_by, na.rm=TRUE),
#             std = sd(!!color_by, na.rm=TRUE)) %>%
#   pivot_wider(names_from = noise_type, values_from = c("m", "std")) %>%
#   ggplot() + 
#   theme(aspect.ratio=1) + 
#   geom_abline() + 
#   geom_errorbar( mapping = aes(x = m_QUIET, y = m_LMBRD,
#                                ymin = m_LMBRD - 1.96*std_LMBRD, 
#                                ymax = m_LMBRD + 1.96*std_LMBRD), 
#                  width = 0) +
#   geom_errorbarh( mapping = aes(x = m_QUIET, y = m_LMBRD,
#                                 xmin = m_QUIET - 1.96*std_QUIET, 
#                                 xmax = m_QUIET + 1.96*std_QUIET), 
#                   height = 0) + 
#   geom_point(aes(x = m_QUIET, 
#                  y = m_LMBRD,
#                  fill = subject_id),
#              color = "black", shape = 22, size = 5,
#              alpha = 0.7, show.legend = TRUE) + 
#   ggtitle(glue('{rlang::as_string(color_by)}-sentences-scatter')) + 
#   xlab('QUIET') + 
#   ylab('LOMBARD')
# 
# 
# stats_by_subject <- compare_means(c(intensity_praat, FCR, F0, duration) ~ noise_type, sentences, group.by = c("subject_id")) %>%
#   mutate(sig = p < 0.05) %>%
#   rename(metric=.y.)
# 
# # Heatmap 
# ggplot(stats_by_subject, aes(metric, subject_id, fill=sig)) + 
#   geom_tile() + 
#   theme(aspect.ratio=1) + 
#   scale_fill_manual(values=color_map) 
# 
# fname = glue('all-metrics-heatmap-wilcoxon')
# # ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), width=6.45, height=5.75, dpi=300)
# dev.print(pdf, glue('{PATH_FIG}/{fname}.pdf'))
# 
# 
# 
# # scale out even further
# stats_by_subject_summ <- stats_by_subject %>%
#   group_by(metric, sig) %>%
#   summarise(n = n()) %>%
#   arrange(desc(sig), desc(n))
#   
# ggplot(stats_by_subject_summ, aes(fill=sig, x=metric, y=n)) + 
#   geom_bar(position="stack", stat="identity") + 
#   scale_fill_manual(values=color_map) + 
#   ylab('# participants') + 
#   xlab("") +
#   ggtitle('Summary of Lombard vs Quiet within-participant')
# 
# fname = glue('all-metrics-barplot-wilcoxon')
# # ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), width=6.45, height=5.75, dpi=300)
# dev.print(pdf, glue('{PATH_FIG}/{fname}.pdf'))
# 
# 
# 
# 
# 
# 
# color_by <- sym("FCR")
# ggplot(vowels_grouped, aes(x=noise_type, y=!!color_by, fill=noise_type, colours())) + 
#   geom_jitter() + 
#   geom_boxplot(width=0.5) + scale_fill_manual(values=color_map) + 
#   stat_compare_means() +
#   facet_wrap(~subject_id, scale="free_y", ncol = 5)
# 
# fname = glue('{rlang::as_string(annot_name)}-by-condition-{rlang::as_string(color_by)}')
# ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300)
# 
# 
# # compare_means(FCR ~ noise_type, data = vowels_grouped, group.by = as.vector(groups(vowels_grouped), 'character'))
# 
# 
# 
# 
# ### plot by 
# df_summ <-
#   sentences %>% 
#   group_by(subject_id, run_id, block_id) %>% 
#   mutate(within_block_trial_id = row_number()) %>% 
#   ungroup() %>%
#   select(subject_id, within_block_trial_id, noise_type, intensity_praat) %>%
#   group_by(subject_id, within_block_trial_id, noise_type) %>%
#   drop_na() %>%
#   summarize_all_stats()
# 
# m = sym("intensity_praat_mean")
# s = sym("intensity_praat_stdev")
# ggplot(df_summ, aes(x=within_block_trial_id, y=!!m, group=noise_type, color=noise_type)) +
#   geom_line() + 
#   geom_point() +
#   scale_x_continuous(breaks=c(1:10)) + 
#   geom_errorbar(aes(ymin=!!m-!!s, ymax=!!m+!!s), width=.2,
#                 position=position_dodge(0.05)) + 
#   facet_wrap(vars(subject_id), ncol=5) 
# fname = glue('within-block-trials--by-subject-lombard-vs-quiet-{rlang::as_string(m)}')
# ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), dpi=300)
# write_tsv(df_summ, glue('{PATH_FIG}/{fname}.tsv'))
# 
# 
# 
# 
# 
# 
# 
# # # combine trials and sentences 
# # trials <- trials %>% 
# #   left_join(sentences, by=c('subject_id', 'run_id', 'trial_id')) %>%
# #   left_join(SUBJECTS_META, by=c('subject_id'))
# 
# 
# voi = c('audio_on', 'audio_off', 'go_time', 'sentence_on', 'sentence_off')
# 
# 
# # tlock to audio_on
# tmp = trials
# tmp[voi] = tmp[voi] - tmp$audio_on
# trials_tlock_audio_on = pivot_longer(tmp, cols=voi, names_to='event', values_to='time') %>%
#   filter(event!='audio_on') %>% 
#   mutate(tlock='audio_on')
# 
# 
# tmp = trials
# tmp[voi] = tmp[voi] - tmp$audio_off
# trials_tlock_audio_off = pivot_longer(tmp, cols=voi, names_to='event', values_to='time') %>%
#   filter(event!='audio_off') %>% 
#   mutate(tlock='audio_off')
# 
# tmp = trials
# tmp[voi] = tmp[voi] - tmp$go_time
# trials_tlock_go_time = pivot_longer(tmp, cols=voi, names_to='event', values_to='time') %>%
#   filter(event!='go_time') %>% 
#   mutate(tlock='go_time')
# 
# 
# tmp = trials
# tmp[voi] = tmp[voi] - tmp$sentence_on
# trials_tlock_sentence_on = pivot_longer(tmp, cols=voi, names_to='event', values_to='time') %>%
#   filter(event!='sentence_on') %>% 
#   mutate(tlock='sentence_on')
# 
# 
# all = bind_rows(trials_tlock_audio_on, trials_tlock_audio_off, trials_tlock_go_time, trials_tlock_sentence_on) 
# all = mutate(all, target=case_when(str_detect(dbs_target, 'STN') ~ 'STN', 
#                                    str_detect(dbs_target, 'GPi') ~ 'GPi', 
#                                    str_detect(dbs_target, 'VIM') ~ 'VIM'))
# 
# 
# 
# # Plot type: points 
# tlock_curr = 'sentence_on'
# color_by <- sym("target")
# 
# all %>% 
#   filter(tlock==tlock_curr) %>%
#   filter(event!='audio_on') %>%
#   ggplot() +
#   aes(x=time, y=!!color_by, color=!!color_by) +
#   #geom_density(color = NA, alpha = 0.4, adjust=1.5) + 
#   geom_jitter() + 
#   geom_boxplot() + 
#   geom_vline(xintercept = 0) + 
#   facet_grid(factor(event, voi) ~ ., scales="free") + 
#   xlim(c(-2, 7)) +
#   theme_bw() + 
#   xlab(glue('Time [s] wrt {tlock_curr}'))
# 
# 
# fname = glue('events-by-{rlang::as_string(color_by)}_plot-type-pointwise_tlock-{tlock_curr}')
# ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), width=10, height=6, dpi=300)
# 
# 
# 
# # Plot type: density 
# tlock_curr = 'sentence_on'
# color_by <- sym("target")
# 
# all %>% 
#   filter(tlock==tlock_curr) %>%
#   filter(event!='audio_on') %>%
#   
#   ggplot() +
#   aes(x=time, color=!!color_by, fill=!!color_by) +
#   geom_density(color = NA, alpha = 0.4, adjust=1.5) + 
#   geom_vline(xintercept = 0) + 
#   facet_grid(factor(event, voi) ~ ., scales="free") + 
#   xlim(c(-2, 7)) +
#   theme_bw() + 
#   xlab(glue('Time [s] wrt {tlock_curr}'))
# 
# 
# fname = glue('events-by-{rlang::as_string(color_by)}_plot-type-density_tlock-{tlock_curr}')
# ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), width=10, height=6, dpi=300)
# 
# 
# 
# 
# 
# # Plot type: per patient 
# tlock_curr = 'audio_off'
# color_by <- sym("target")
# 
# 
# all %>% arrange(target) %>% 
#   mutate(subject_id=factor(subject_id)) %>%
#   filter(tlock==tlock_curr) %>%
#   filter(event!='audio_on') %>% 
#   
#   ggplot() +
#   aes(x=time, y=(subject_id), color=!!color_by) +
#   geom_jitter() + 
#   geom_boxplot() + 
#   geom_vline(xintercept = 0) + 
#   facet_grid(factor(event, voi) ~ ., scales="fixed") + 
#   xlim(c(-2, 7)) +
#   theme_bw() + 
#   xlab(glue('Time [s] wrt {tlock_curr}'))
# 
# 
# fname = glue('events-by-{rlang::as_string(color_by)}_plot-type-pointwise-per-patient_tlock-{tlock_curr}')
# ggsave(file=glue('{PATH_FIG}/{fname}.pdf'), width=10, height=6, dpi=300)









