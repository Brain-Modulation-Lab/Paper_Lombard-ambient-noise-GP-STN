library(tidyverse)
library(testit)
library(glue)
library(assertthat)
library(readtextgrid)

# Defining paths
SUBJECT <- 'DM####'
SESSION <- 'intraop'
TASK <- 'lombard'

FLANKING_GAP <- 0.0 #time in seconds flanking the identified phrase to feed to force aligner

PATH_DATASET <- '/Volumes/Nexus4/DBS'
PATH_SRC <- file.path(PATH_DATASET, 'sourcedata')
PATH_SRC_SUB <- file.path(PATH_SRC, paste0('sub-',SUBJECT)) 
PATHS_TASK <- file.path(PATH_SRC_SUB,c('ses-training','ses-preop','ses-intraop'),'task')
PATH_TASK_INTRAOP <- file.path(PATH_SRC_SUB,'ses-intraop','task')
PATH_DER <- file.path(PATH_DATASET, 'derivatives')
PATH_DER_SUB <- file.path(PATH_DER, paste0('sub-',SUBJECT))  
PATH_PREPROC <- file.path(PATH_DER_SUB, 'preproc')
PATH_ANNOT <- file.path(PATH_DER_SUB, 'annot')
PATH_AEC <- file.path(PATH_DER_SUB,'aec')
PATH_P2FA <- file.path(PATH_DER_SUB,'p2fa')
PATH_P2FA_TEXT2 <- file.path(PATH_P2FA,'text-rerun')
PATH_P2FA_FA2 <- file.path(PATH_P2FA,'forced-alignment-rerun')
PATH_PHONETIC_CODING <- file.path(PATH_DER_SUB,'phonetic-coding')
P2FA_RATE <- 11025

setwd(PATH_PREPROC)

# Loading trials annotation table
basename <- glue('sub-{SUBJECT}_ses-{SESSION}_task-{TASK}')

trials_fpath <- glue('{PATH_ANNOT}/{basename}_annot-trials.tsv')
assert_that(file.exists(trials_fpath))
trials<-read_tsv(trials_fpath)

sync_fpath <- glue('{PATH_ANNOT}/sub-{SUBJECT}_ses-{SESSION}_sync.tsv')
assert_that(file.exists(sync_fpath))
sync<-read_tsv(sync_fpath)

# checking preproc directories 
assert_that(file.exists(PATH_AEC))
if (!file.exists(PATH_P2FA)) error(glue('{PATH_P2FA} should exists'))
if (!file.exists(PATH_PHONETIC_CODING)) error(glue('{PATH_PHONETIC_CODING} should exists'))
if (!file.exists(PATH_P2FA_TEXT2)) dir.create(PATH_P2FA_TEXT2)
if (!file.exists(PATH_P2FA_FA2)) dir.create(PATH_P2FA_FA2)


# loading names of TextGrid files
phonetic_coding_files <- 
  tibble(path = list.dirs(PATH_PHONETIC_CODING,full.names=FALSE,recursive=FALSE)) %>%
  filter(str_detect(path,'sub-[\\w\\d]+_ses-[\\w]+_task-lombard_run-[\\d]+')) %>%
  rowwise() %>%
  group_modify(~file.info(list.files(file.path(PATH_PHONETIC_CODING,.x$path),pattern=glob2rx('*.TextGrid'),full.names=TRUE))) %>%
  rownames_to_column('path') %>%
  filter(size>0 & isdir==FALSE) %>%
  select(path)

# Reading textgrid files
phonetic_coding_files <- phonetic_coding_files$path
names(phonetic_coding_files) <- phonetic_coding_files
phonetic_coding <- phonetic_coding_files %>% 
  map_df(read_textgrid, .id = 'path') %>%
  mutate(path=dirname(path)) %>%
  mutate(subject =str_match(file,'sub-([\\w\\d]+)_ses-[\\w]+_task-[\\w]+_run-[\\d]+_trial-[\\d]+_p2fa_[\\w]+')[,2]) %>%
  mutate(session =str_match(file,'sub-[\\w\\d]+_ses-([\\w]+)_task-[\\w]+_run-[\\d]+_trial-[\\d]+_p2fa_[\\w]+')[,2]) %>%
  mutate(task    =str_match(file,'sub-[\\w\\d]+_ses-[\\w]+_task-([\\w]+)_run-[\\d]+_trial-[\\d]+_p2fa_[\\w]+')[,2]) %>%
  mutate(run_id  =str_match(file,'sub-[\\w\\d]+_ses-[\\w]+_task-[\\w]+_run-([\\d]+)_trial-[\\d]+_p2fa_[\\w]+')[,2]) %>%
  mutate(trial_id=str_match(file,'sub-[\\w\\d]+_ses-[\\w]+_task-[\\w]+_run-[\\d]+_trial-([\\d]+)_p2fa_[\\w]+')[,2]) %>%
  mutate(type    =str_match(file,'sub-[\\w\\d]+_ses-[\\w]+_task-[\\w]+_run-[\\d]+_trial-[\\d]+_p2fa_([\\w]+)')[,2]) %>%
  mutate(run_id=as.numeric(run_id),trial_id=as.numeric(trial_id)) %>%
  relocate(subject,session,task,run_id,trial_id) %>%
  relocate(path,file,.after = last_col())
  
phonetic_coding_rerun <- phonetic_coding %>%
  filter(type=="rerun") %>%
  filter(tier_num==3) %>%
  filter(text!='sp') %>%
  group_by(subject,session,run_id,trial_id) %>%
  mutate(phrase_id=1:n(),.after=trial_id)

# looping through runs
run_ids <- unique(phonetic_coding_rerun$run_id)
for (i in run_ids){
  
  #check if aec wav file exists
  wav_fname_out <- glue('{basename}_run-{formatC(i,width=2,flag="0")}_recording-directionalmicaec-{P2FA_RATE}Hz_physio.wav')
  wav_fpath_out <- glue('{PATH_P2FA}/{wav_fname_out}')
  
  if (!file.exists(wav_fpath_out)){
    warning(glue('{wav_fname_out} not found'))
    next    
  }  
  
  phrase_rerun <- phonetic_coding_rerun %>%
    filter(run_id==i) %>%
    mutate(starts = xmin - FLANKING_GAP, ends = xmax + FLANKING_GAP) %>%
    mutate(text_fpath = glue("{PATH_P2FA_TEXT2}/{basename}_run-{formatC(run_id,width=2,flag='0')}_trial-{formatC(trial_id,width=2,flag='0')}_phrase-{formatC(phrase_id,width=2,flag='0')}.txt")) %>%
    mutate(textgrid_fpath = glue("{PATH_P2FA_FA2}/{basename}_run-{formatC(run_id,width=2,flag='0')}_trial-{formatC(trial_id,width=2,flag='0')}_phrase-{formatC(phrase_id,width=2,flag='0')}_p2fa.TextGrid")) %>%
    mutate(cmd=glue("conda run -n p2fa python2 /usr/local/p2fa/align.py -s {starts} -e {ends} {wav_fpath_out} {text_fpath} {textgrid_fpath}"))
  
  #phrase_rerun  <- phrase_rerun %>% filter(trial_id==40)
  
  #creating text files
  phrase_rerun %>%
    rowwise() %>%
    group_walk(~cat(.x$text,file=.x$text_fpath,sep='',append=FALSE))
  
  #running the forced aligner
  phrase_rerun %>%
    rowwise() %>%
    group_walk(~ system(.x$cmd))

}

                     
                           




