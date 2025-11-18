
library(tidyverse)
library(testit)
library(glue)

# Defining paths
SUBJECT <- 'DM####'
SESSION <- 'intraop'
TASK <- 'lombard'

PATH_DATASET <- '/Volumes/Nexus4/DBS'
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
}

# loading events tables
events <- dir(PATH_ANNOT, pattern = "^.*_task-lombard_.*_events-sync\\.tsv$", full.names = TRUE) %>%
  map_df(read_events_tsv)

# loading trials tables
trials <- dir(PATH_TASK_INTRAOP, pattern = "^.*_task-lombard_.*_trials\\.tsv$", full.names = TRUE) %>%
  map_df(read_bids_tsv) %>%
  rename(run_id=run)


# THIS PART IS TASK-SPECIFIC!!!! For example smsl task will have a different section compared to lombard task because the protocol is different. Change this portion accordingly. To avoid any mistakes, compare the scripts sub-DM1007_ses-intraop_task-lombard_proto-A05_triggers-to-annot_20210729 and sub-DM1007_ses-intraop_task-smsl_proto-A05_triggers-to-annot_20210729 to check the differences in this section. If you are working with another task, you need to change it accordingly.


# decoding triggers
TRIG_FIXCROSS = 1
TRIG_AUDITORY = 2
TRIG_VISUAL = 4
TRIG_GO = 8
TRIG_KEY = 16
TRIG_VOLUME = 32
TRIG_UP = 64
TRIG_NOISE = 128

events_trig <- events %>%
  rowwise() %>%
  mutate(trig_fixcross = as.integer(intToBits(value))[log2(TRIG_FIXCROSS)+1]) %>%
  mutate(trig_auditory = as.integer(intToBits(value))[log2(TRIG_AUDITORY)+1]) %>%
  mutate(trig_visual = as.integer(intToBits(value))[log2(TRIG_VISUAL)+1]) %>%
  mutate(trig_go = as.integer(intToBits(value))[log2(TRIG_GO)+1]) %>%
  mutate(trig_key = as.integer(intToBits(value))[log2(TRIG_KEY)+1]) %>%
  mutate(trig_up = as.integer(intToBits(value))[log2(TRIG_UP)+1]) %>%
  mutate(trig_noise = as.integer(intToBits(value))[log2(TRIG_NOISE)+1]) %>%
  ungroup()

##### Defining stim-audio-noise annotation table

events_trig_trial<-events_trig %>%
  filter(event_code != 'Noise Start') %>%
  filter(event_code != 'Welcome Message') %>%
  rename(run_id=run) %>%
  group_by(run_id) %>%
  mutate(trig_visual_onset = (trig_visual - lag(trig_visual,default=0))>0) %>%
  mutate(trig_visual_onset = ifelse(value>0,trig_visual_onset,0)) %>% ### AM changed last ifelse arg from NA to 0; if there are NAs, cumsum function doesn’t work properly
  mutate(trial_id = cumsum(trig_visual_onset)) %>%
  filter(!is.na(trial_id))
events_trig_trial$event_code[events_trig_trial$event_code=="Key Pressed"]<-"Key Press"

time_col_names<-tibble(
  event_code=c("Key Press","Fixation Cross","Fading Noise In","Fading Noise Complete","Sentence Visual Onset","Sentence Audio","End Sentence Audio","Go","Escape","Fading Noise Out"),     
  time_col_name=c("keypress_time","itg_onset","fading_in_starts","fading_complete"  ,"onset"                ,"audio_onset"  ,"audio_ends",   "go_time","escaped","fading_out_starts"))


#checking mapping of audio files
mismatched_audio_trials <- events_trig_trial %>%
  filter(event_code == "Sentence Audio") %>% 
  select(sub,ses,task,run_id,trial_id,stim_file) %>%
  left_join(trials) %>%
  filter(stim_file != file_audio) %>%
  select(sub,ses,task,run_id,trial_id,stim_file,file_audio) 
assert(nrow(mismatched_audio_trials)==0)

#checking mapping of text files
mismatched_text_trials <- events_trig_trial %>%
  filter(event_code == "Sentence Visual Onset") %>% 
  select(sub,ses,task,run_id,trial_id,stim_file) %>%
  left_join(trials) %>%
  filter(stim_file != file_text) %>%
  select(sub,ses,task,run_id,trial_id,stim_file,file_text) 
assert(nrow(mismatched_text_trials)==0)


#creating annot table
annot_trials <- events_trig_trial %>%
  left_join(time_col_names) %>%
  group_by(sub,ses,task,run_id,trial_id) %>%
  mutate(trig_noise=trig_noise[1]) %>%
  ungroup() %>%
  select(sub,ses,task,run_id,trial_id,trig_noise,time_col_name,onset) %>%
  drop_na(time_col_name) %>% # added LB 2023 111 09 because there were a some 'NaN' events
  group_by(sub,ses,task,run_id,trial_id,trig_noise) %>%
  pivot_wider(names_from=time_col_name,values_from=onset) %>%
  ungroup() %>%
  mutate(itg_noise_type = if_else(is.na(fading_complete),
                                  if_else(trig_noise==1,"babble","none"),
                                  if_else(is.na(fading_out_starts),"fading_in","fading_out")
  )
  ) %>%
  left_join(trials) %>%
  group_by(sub,ses,task,run_id) %>%
  mutate(itg_ends = lead(onset,1)) %>%
  mutate(itg_ends = if_else(is.na(itg_ends),fading_complete,itg_ends)) %>%
  ungroup()


#checking mapping of noise conditions
mismatched_noise_trials <- annot_trials %>%
  filter(trig_noise != noise_type)
assert(nrow(mismatched_noise_trials)==0)

#cleaning up 
annot_trials_clean <- annot_trials %>%
  rename(sentence_text=sentence) %>%
  mutate(audio_volume=1) %>%
  mutate(duration=itg_onset - onset) %>%
  mutate(audio_duration = audio_ends - audio_onset) %>%
  mutate(itg_duration = itg_ends - itg_onset) %>%
  select(onset,duration,run_id,block_id,trial_id,sentence_id,sentence_text,noise_type,file_audio,
         audio_onset,audio_duration,go_time,keypress_time,itg_onset,itg_duration,itg_noise_type,audio_volume)


#checking nominal audio vs actual audio duration
mismatched_audio_duration <- events_trig_trial %>%
  filter(event_code == "Sentence Audio") %>% 
  select(sub,ses,task,run_id,trial_id,duration) %>%
  rename(audio_duration_nominal=duration) %>%
  left_join(annot_trials_clean) %>% 
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
annot_trials_clean %>%
  write_tsv(file.path(PATH_ANNOT,annot_fname))
