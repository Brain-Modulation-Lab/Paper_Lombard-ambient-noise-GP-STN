library(tidyverse)
library(glue)

SUBJECT <- 'DM1046'
SESSION <- 'intraop'
TASK <- 'ECoGeStim'

PATH_DATASET <- '/Volumes/Nexus4/DBS'
#PATH_DATASET <- 'Y:/DBS'
PATH_DER <- file.path(PATH_DATASET, 'derivatives')
PATH_DER_SUB <- file.path(PATH_DER, paste0('sub-',SUBJECT))  
PATH_PREPROC <- file.path(PATH_DER_SUB, 'preproc')
PATH_ANNOT <- file.path(PATH_DER_SUB, 'annot')
PATH_SRC <- file.path(PATH_DATASET, 'sourcedata')
PATH_SRC_SUB <- file.path(PATH_SRC, paste0('sub-',SUBJECT)) 
PATHS_TASK <- file.path(PATH_SRC_SUB,c('ses-training','ses-preop','ses-intraop'),'task')
PATH_TASK_INTRAOP <- file.path(PATH_SRC_SUB,'ses-intraop','task')

setwd(PATH_PREPROC)

read_bids_tsv <- function(file,...){
  tab <- read_tsv(file,...)  
  fname <- basename(tools::file_path_sans_ext(file))
  sub <- str_match(fname,'sub-([a-zA-Z0-9]+)_')[2] 
  ses <- str_match(fname,'ses-([a-zA-Z0-9]+)_')[2] 
  task <- str_match(fname,'task-([a-zA-Z0-9-]+)_')[2] 
  run <- as.numeric(str_match(fname,'run-([0-9]+)_')[2])
  return(tab %>% mutate(sub=sub,ses=ses,task=task,run=run))
}

events <- dir(PATH_ANNOT, pattern = paste0("^.*_task-",TASK,"_.*_events-sync\\.tsv$"), full.names = TRUE) %>%
  map_df(read_bids_tsv) %>%
  mutate(value8bit = value %% 2^9) %>%
  filter(!event_code %in% c("run id","clock hour","clock minute","clock second", "clock millisecond","Trigger test"))

task_runs <- read_bids_tsv(glue("{PATH_ANNOT}/sub-{SUBJECT}_runs.tsv")) %>%
  filter(session==SESSION) %>%
  filter(task==TASK) %>%
  select(onset, duration, everything())

blocks <- dir(PATH_TASK_INTRAOP, pattern = paste0("^.*_task-",TASK,"_.*_trials\\.tsv$"), full.names = TRUE) %>%
  map_df(read_bids_tsv) %>%
  rename(run_id=run, block_id=trial_id)

pulses <- dir(PATH_ANNOT, pattern = paste0("^.*_task-",TASK,"*_stimulation-pulses\\.tsv$"), full.names = TRUE) %>%
  map_df(read_bids_tsv) %>%
  select(-run) %>%
  cross_join(task_runs %>% select (onset,duration,run) %>% rename(run_onset = onset, run_duration=duration)) %>%
  filter(onset >= run_onset, onset <= run_onset + run_duration) %>%
  group_by(run) %>%
  mutate(trial_id=dense_rank(sample)) %>%
  group_by(run, trial_id) %>%
  mutate(stimchannel_idx = row_number())

bipolar_pulses <- pulses %>%
  pivot_wider(id_cols = c(onset, duration, type, sample, file_id, filename, sub, ses, task, run, trial_id), 
              values_from = stimchannel, names_from = stimchannel_idx, names_prefix = "stimchannel") %>%
  rename(pulse_onset = onset, run_id = run)

# decoding triggers
TRIG_ITG = 4;
TRIG_TRIAL = 8;
TRIG_KEY = 16;
TRIG_STIM = 32;

events_trig <- events %>%
  rowwise() %>%
  mutate(trig_itg = as.integer(intToBits(value))[log2(TRIG_ITG)+1]) %>%
  mutate(trig_trial = as.integer(intToBits(value))[log2(TRIG_TRIAL)+1]) %>%
  mutate(trig_key = as.integer(intToBits(value))[log2(TRIG_KEY)+1]) %>%
  mutate(trig_stim = as.integer(intToBits(value))[log2(TRIG_STIM)+1]) %>%
  ungroup()

events_trig_trial<-events_trig %>%
  filter(event_code != 'Key Pressed') %>%
  rename(run_id=run) %>%
  group_by(run_id) %>%
  mutate(trial_id = cumsum(trig_stim)) %>%
  mutate(block_id = cumsum(trig_itg)+1) %>%
  filter(!is.na(trial_id)) %>%
  filter(event_code != "ITG") %>%
  filter(event_code != 'End Message') 

time_col_names<-tibble(
  event_code=c(   "STIM",    "ERNA", "Zero"),     
  time_col_name=c("stim_trig_time","erna_trig_time","stim_trig_time"))

annot_trials <- events_trig_trial %>%
  left_join(time_col_names) %>%
  select(sub,ses,task,run_id,block_id,trial_id,trial_type,time_col_name,onset) %>%
  group_by(sub,ses,task,run_id,block_id,trial_id,trial_type,time_col_name) %>%
  filter(n()==1) %>% 
  pivot_wider(names_from=time_col_name,values_from=onset) %>%
  ungroup() %>%
  mutate(stim_trig_next = lead(stim_trig_time)) %>%
  filter(!is.na(erna_trig_time)) %>%
  left_join( bipolar_pulses %>% 
      select(pulse_onset, stimchannel1, stimchannel2, sub, ses, task, run_id, trial_id)) 

  
#checking nominal vs actual stim time
annot_trials %>%
  ggplot()+
  aes(x=stim_trig_time - pulse_onset, fill=factor(run_id)) +
  geom_histogram(binwidth=0.0001)+
  labs(x="mismatch trigger to stimulation pulse (s)")
ggsave(glue("figures/sub-{SUBJECT}_ses-{SESSION}_task-{TASK}_proto-A05_plot-stim-mismatch.png"))

annot_trials_clean <- annot_trials %>%
  filter(!is.na(pulse_onset)) %>%
  mutate(onset=pulse_onset, duration = stim_trig_next - onset) %>%
  select(onset, duration, sub, ses, task, run_id, trial_id, everything())

#writing annotation table
annot_trials_clean %>%
  write_tsv(file.path(PATH_ANNOT,glue("sub-{SUBJECT}_ses-{SESSION}_task-{TASK}_annot-trials.tsv")))
