library(tidyverse)
library(testit)
library(glue)
library(assertthat)

# Defining paths
SUBJECT <- 'DM10XX'
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
PATH_AEC <- file.path(PATH_DER_SUB,'aec')
PATH_P2FA <- file.path(PATH_DER_SUB,'p2fa')
PATH_P2FA_TEXT <- file.path(PATH_P2FA,'text')
PATH_P2FA_FA <- file.path(PATH_P2FA,'forced-alignment')
P2FA_RATE <- 11025


# Loading trials annotation table
basename <- glue('sub-{SUBJECT}_ses-{SESSION}_task-{TASK}')

trials_fpath <- glue('{PATH_ANNOT}/{basename}_annot-trials.tsv')
assert_that(file.exists(trials_fpath))
trials<-read_tsv(trials_fpath)

sync_fpath <- glue('{PATH_ANNOT}/sub-{SUBJECT}_sync.tsv')
assert_that(file.exists(sync_fpath))
sync<-read_tsv(sync_fpath)

# checking preproc directories 
assert_that(file.exists(PATH_AEC))
if (!file.exists(PATH_P2FA)) dir.create(PATH_P2FA)
if (!file.exists(PATH_P2FA_TEXT)) dir.create(PATH_P2FA_TEXT)
if (!file.exists(PATH_P2FA_FA)) dir.create(PATH_P2FA_FA)


# looping through runs
run_ids <- unique(trials$run_id)
for (i in 1:length(run_ids)){
  
  #check if aec wav file exists
  wav_fname <- glue('{basename}_run-{formatC(run_ids[i],width=2,flag="0")}_recording-directionalmicaec_physio.wav')
  wav_fpath <- glue('{PATH_AEC}/{wav_fname}')
  wav_fname_out <- glue('{basename}_run-{formatC(run_ids[i],width=2,flag="0")}_recording-directionalmicaec-{P2FA_RATE}Hz_physio.wav')
  wav_fpath_out <- glue('{PATH_P2FA}/{wav_fname_out}')
  
  if (!file.exists(wav_fpath)){
    warning(glue('{wav_fname} not found'))
    next    
  }  
  
  #resampling and moving wav file
  system(glue('sox {wav_fpath} -r {P2FA_RATE} {wav_fpath_out}'))
  
  #getting t0 for this audio file from sync table
  t0 <- sync %>%
    filter(name == wav_fname) %>%
    with(t1 + (1-s1) * (t2-t1)/(s2-s1)) 
  assert_that(length(t0)==1)
  
  trials_run <- trials %>%
    filter(run_id==run_ids[i]) %>% 
    filter(!is.na(onset)) %>% 
    mutate(exe_starts = audio_onset + audio_duration - t0 - 1) %>% #audio offset in ~wav time -1s
    mutate(exe_ends = keypress_time - t0 + 1) %>% #keypress time +1s in ~wav time
    mutate(text_fpath = glue("{PATH_P2FA_TEXT}/{basename}_run-{formatC(run_id,width=2,flag='0')}_trial-{formatC(trial_id,width=2,flag='0')}_stim.txt")) %>%
    mutate(textgrid_fpath = glue("{PATH_P2FA_FA}/{basename}_run-{formatC(run_id,width=2,flag='0')}_trial-{formatC(trial_id,width=2,flag='0')}_p2fa.TextGrid")) %>%
    mutate(cmd=glue("python2 /usr/local/p2fa/align.py -s {exe_starts} -e {exe_ends} {wav_fpath_out} {text_fpath} {textgrid_fpath}"))
  
  #creating text files
  trials_run %>%
    rowwise() %>%
    group_walk(~cat(.x$sentence_text,file=.x$text_fpath,sep='',append=FALSE))
  
  #running the forced aligner
  trials_run %>%
    rowwise() %>%
    group_walk(~ system(.x$cmd))
  
}

