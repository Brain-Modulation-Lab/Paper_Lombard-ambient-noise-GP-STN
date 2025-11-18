# AM

library(tidyverse)
library(testit)
library(glue)

# Defining paths
SUBJECT <- 'DM####'
SESSION <- 'intraop'
TASK <- 'smsl'
RUN <- '' # added by AM for DM1045 because there were some aborted smsl runs; can leave empty if there's only one run

if(Sys.info()[[4]]=="MSI" || Sys.info()[[4]]=="677-GUE-WL-0010") {
  PATH_DATASET <- 'Y:/DBS'  ### drive mapped to AM laptop
} else {
  PATH_DATASET <- '/Volumes/Nexus4/DBS'
}
  
PATH_DER <- file.path(PATH_DATASET, 'derivatives')
PATH_DER_SUB <- file.path(PATH_DER, paste0('sub-',SUBJECT))  
PATH_PREPROC <- file.path(PATH_DER_SUB, 'preproc')
PATH_ANNOT <- file.path(PATH_DER_SUB, 'annot')
PATH_SRC <- file.path(PATH_DATASET, 'sourcedata')
PATH_SRC_SUB <- file.path(PATH_SRC, paste0('sub-',SUBJECT)) 
PATHS_TASK <- file.path(PATH_SRC_SUB,c('ses-training','ses-preop','ses-intraop'),'task')
PATH_TASK_INTRAOP <- file.path(PATH_SRC_SUB,'ses-intraop','task')

setwd(PATH_PREPROC)

# defining utility functions to load tsv tables
read_bids_tsv <- function(file,...){
  tab <- read_tsv(file,...)  
  fname <- basename(tools::file_path_sans_ext(file))
  sub <- str_match(fname,'sub-([a-zA-Z0-9]+)_')[2] 
  ses <- str_match(fname,'ses-([a-zA-Z0-9]+)_')[2] 
  task <- str_match(fname,'task-([a-zA-Z0-9]+)_')[2] 
  run <- as.numeric(str_match(fname,'run-([0-9]+)_')[2])
  return(tab %>% mutate(sub=sub,ses=ses,task=task,run=run))
}

read_events_tsv <- function(...){
  events <- read_bids_tsv(...)
 return(events[(max(which(events$event_code=='Trigger test'))+1):nrow(events),])
  ## in the case of 1008 intraop, Trigger test events have already been cut off in sync file, so start at first Key Press instead
  # return(events[(min(which(events$event_code=='Key Press'))):nrow(events),])
}

# loading events tables
events <-  dir(PATH_ANNOT, pattern = paste("^.*_ses-intraop_task-SMSL_run-",RUN,".*_events-sync\\.tsv$",sep=""), full.names = TRUE) %>%
  map_df(read_events_tsv)

# loading trials tables
trials <- dir(PATH_TASK_INTRAOP, pattern = "^.*_ses-intraop_task-SMSL_.*_spec-Test\\.tsv$", full.names = TRUE) %>%
  map_df(read_bids_tsv) %>%
  rename(run_id=run) %>%
  group_by(sub,ses,task,run_id) %>%
  mutate(trial_id=1:n())

# decoding triggers
TRIG_FIXCROSS = 1
TRIG_AUDITORY = 2
TRIG_VISUAL = 4
TRIG_GO = 8
TRIG_STOP = 16
TRIG_KEYPRESS = 32
TRIG_ITG = 64
TRIG_ESC = 128
TRIG_FLIPSYNC = 2^15

events_trig <- events %>%
  rowwise() %>%
  mutate(trig_fixcross = as.integer(intToBits(value))[log2(TRIG_FIXCROSS)+1]) %>%
  mutate(trig_auditory = as.integer(intToBits(value))[log2(TRIG_AUDITORY)+1]) %>%
  mutate(trig_visual = as.integer(intToBits(value))[log2(TRIG_VISUAL)+1]) %>%
  mutate(trig_go = as.integer(intToBits(value))[log2(TRIG_GO)+1]) %>%
  mutate(trig_stop = as.integer(intToBits(value))[log2(TRIG_STOP)+1]) %>%
  mutate(trig_keypress = as.integer(intToBits(value))[log2(TRIG_KEYPRESS)+1]) %>%
  mutate(trig_itg = as.integer(intToBits(value))[log2(TRIG_ITG)+1]) %>%
  mutate(trig_esc = as.integer(intToBits(value))[log2(TRIG_ESC)+1]) %>%
  mutate(trig_flipsync = as.integer(intToBits(value))[log2(TRIG_FLIPSYNC)+1]) %>%
  ungroup()

##### Defining stim-audio-noise annotation table

events_trig_trial<-events_trig %>%
  filter(event_code != 'Melcome message') %>%
  rename(run_id=run) %>%
  group_by(run_id) %>%
  mutate(trig_visual_onset = (trig_visual - lag(trig_visual,default=0))>0) %>%
  mutate(trig_visual_onset = ifelse(value>0,trig_visual_onset,0)) %>%
  mutate(trial_id = cumsum(trig_visual_onset)) %>%
  filter(!is.na(trial_id)) %>%
  mutate(to_remove = event_code == "Break message") %>% #removing break messages events 
  mutate(to_remove2 = event_code == "Key Press" & lag(event_code,1) ==  "Break message") %>%
  filter(!to_remove) %>%
  filter(!to_remove2) %>%
  mutate(to_remove = event_code == "Go Visual") %>% #removing Go visual events 
  mutate(to_remove2 = event_code == "Fixation Cross" & (lag(event_code,1) ==  "Go Visual" | lag(event_code,2) ==  "Go Visual")) %>%
  filter(!to_remove) %>%
  filter(!to_remove2) %>%
  mutate(to_remove = event_code == "Stop Visual") %>% #removing Stop visual events 
  mutate(to_remove2 = event_code == "Fixation Cross" & (lag(event_code,1) ==  "Stop Visual" | lag(event_code,2) ==  "Stop Visual"| lag(event_code,2) ==  "Stop Visual")) %>%
  filter(!to_remove) %>%
  filter(!to_remove2) %>%
  select(-to_remove,-to_remove2)
  
time_col_names<-tibble(
  event_code=c("Key Press","Fixation Cross","Visual Stim Onset","Audio Stim Onset","Audio Stim Offset","Go Beep","Stop Buzzer"),     
  time_col_name=c("keypress_time","visual_offset","visual_onset"   ,"audio_onset"     ,"audio_offset", "audio_go_onset","audio_stop_onset"))


#checking mapping of audio files
mismatched_audio_trials <- events_trig_trial %>%
  filter(event_code == "Audio Stim Onset") %>% 
  select(sub,ses,task,run_id,trial_id,stim_file) %>%
  rowwise() %>%
  mutate(stim_word = str_match(stim_file,"\\d{2}_([a-z]*).wav")[2]) %>%
  left_join(trials) %>%
  filter(stim_word != str_to_lower(word)) %>%
  select(sub,ses,task,run_id,trial_id,stim_file) 
assert(nrow(mismatched_audio_trials)==0)

#checking mapping of text files
mismatched_text_trials <- events_trig_trial %>%
  filter(event_code == "Visual Stim Onset") %>% 
  select(sub,ses,task,run_id,trial_id,stim_file) %>%
  left_join(trials) %>%
  filter(stim_file != word) %>%
  select(sub,ses,task,run_id,trial_id,stim_file) 
assert(nrow(mismatched_text_trials)==0)


#creating annot table
annot_trials <- events_trig_trial %>%
  left_join(time_col_names) %>%
  filter(!is.na(time_col_name)) %>%
  select(sub,ses,task,run_id,trial_id,time_col_name,onset) %>%
  group_by(sub,ses,task,run_id,trial_id) %>%
  pivot_wider(names_from=time_col_name,values_from=onset) %>%
  ungroup() %>%
  left_join(trials) %>%
  mutate(audio_go_offset = audio_go_onset + 0.05) %>%
  mutate(audio_stop_offset = audio_stop_onset + 0.05) %>%
  mutate(itg_starts = keypress_time, itg_ends = lead(visual_onset,1)) %>%
  mutate(onset = visual_onset, duration = keypress_time - visual_onset) %>%
  mutate(stim_id = stimnum, stim_condition = stim_conditions) %>%
  mutate(block_id = ((trial_id-1)%/%48)+1) %>%
  select(onset,duration,run_id,block_id,trial_id,stim_id,word,stim_condition,is_stoptrial,
         stop_latency_ms,visual_onset,audio_onset,audio_offset,visual_offset,
         audio_go_onset,audio_go_offset,audio_stop_onset,audio_stop_offset,keypress_time,
         itg_starts,itg_ends)

#checking nominal audio vs actual audio duration
mismatched_audio_duration <- events_trig_trial %>%
  filter(event_code == "Audio Stim Onset") %>% 
  select(sub,ses,task,run_id,trial_id,duration) %>%
  rename(audio_duration_nominal=duration) %>%
  left_join(annot_trials) %>% 
  mutate(audio_duration = audio_offset - audio_onset) %>%
  mutate(audio_duration_mismatch = audio_duration - audio_duration_nominal)
assert(max(abs(mismatched_audio_duration$audio_duration_mismatch))<1)

mismatched_audio_duration %>%
  ggplot()+
  aes(x=audio_duration_mismatch*1000,fill=factor(run_id)) +
  geom_histogram(binwidth=0.1)+
  labs(x="mismatch audio duration (ms)")
ggsave(glue("figures/sub-{SUBJECT}_ses-{SESSION}_task-{TASK}_proto-A05_plot-nominal-audio-duration-mismatch.png"))

#writing annotation table
annot_fname <- glue("sub-{SUBJECT}_ses-{SESSION}_task-{TASK}_annot-trials.tsv")
annot_trials %>%
  write_tsv(file.path(PATH_ANNOT,annot_fname))







