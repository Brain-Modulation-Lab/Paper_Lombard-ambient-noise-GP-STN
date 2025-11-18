library(tidyverse)
library(testit)
library(glue)
library(assertthat)
library(readtextgrid)

# Defining paths
SUBJECT <- 'DM1038'
SESSION <- 'intraop'
TASK <- 'lombard'

PATH_DATASET <- '/Volumes/Nexus4/DBS'
PATH_DER <- file.path(PATH_DATASET, 'derivatives')
PATH_DER_SUB <- file.path(PATH_DER, paste0('sub-',SUBJECT))  
PATH_PREPROC <- file.path(PATH_DER_SUB, 'preproc')
PATH_ANNOT <- file.path(PATH_DER_SUB, 'annot')
PATH_AEC <- file.path(PATH_DER_SUB,'aec')
PATH_PHONETIC_CODING <- file.path(PATH_DER_SUB,'whisperx')

# Loading trials annotation table
basename <- glue('sub-{SUBJECT}_ses-{SESSION}_task-{TASK}')

trials_fpath <- glue('{PATH_ANNOT}/{basename}_annot-trials.tsv')
assert_that(file.exists(trials_fpath))
trials<-read_tsv(trials_fpath)

sync_fpath <- glue('{PATH_ANNOT}/sub-{SUBJECT}_ses-intraop_sync.tsv')
assert_that(file.exists(sync_fpath))
sync<-read_tsv(sync_fpath)

# specific fix for DM1038... there should only be one row per filename
sync <-
  sync %>% 
  arrange(name) %>%
  group_by(name) %>%
  filter(row_number()==1)
  

# checking preproc directories 
if (!file.exists(PATH_PHONETIC_CODING)) error(glue('{PATH_PHONETIC_CODING} should exists'))



# loading names of TextGrid files
phonetic_coding_files <- 
  tibble(path = list.files(PATH_PHONETIC_CODING,full.names=TRUE,recursive=FALSE)) %>%
  filter(str_detect(path,glue(glob2rx('*{TASK}*directionalmicaec_physio_whisper-all.TextGrid'))))
phonetic_coding_files <- phonetic_coding_files$path
names(phonetic_coding_files) <- phonetic_coding_files

assert_that(length(phonetic_coding_files)>0)


# defining utility functions to load tsv textgrids
read_bids_textgrid <- function(file,...){
  tab <- read_textgrid(file,...)  
  fname <- basename(tools::file_path_sans_ext(file))
  subject_id <- str_match(fname,'sub-([a-zA-Z0-9]+)_')[2] 
  session_id <- str_match(fname,'ses-([a-zA-Z0-9]+)_')[2] 
  task_id <- str_match(fname,'task-([a-zA-Z0-9]+)_')[2] 
  run_id <- as.numeric(str_match(fname,'run-([0-9]+)_')[2])
  trial_id <- as.numeric(str_match(fname,'trial-([0-9]+)_')[2])
  phrase_id <- as.numeric(str_match(fname,'phrase-([0-9]+)_')[2])
  type <- str_match(fname,'p2fa_([a-zA-Z0-9]+)')[2]
  return(tab %>% mutate(subject_id=subject_id,session_id=session_id,task_id=task_id,run_id=run_id,trial_id=trial_id,phrase_id=phrase_id,type=type))
}

# Reading textgrid files
phonetic_coding <- phonetic_coding_files %>% 
  map_df(read_bids_textgrid, .id = 'path') %>%
  mutate(name_wav = str_replace(file, '(?<=physio).*$', '.wav'))

#adding format options to glue from https://stackoverflow.com/questions/63731098/zero-padding-with-glue
sprintf_transformer <- function(text, envir) {
  m <- regexpr(":.+$", text)
  if (m != -1) {
    format <- substring(regmatches(text, m), 2)
    regmatches(text, m) <- ""
    res <- eval(parse(text = text, keep.source = FALSE), envir)
    do.call(sprintf, list(glue("%{format}"), res))
  } else {
    eval(parse(text = text, keep.source = FALSE), envir)
  }
}

glue_fmt <- function(..., .envir = parent.frame()) {
  glue(..., .transformer = sprintf_transformer, .envir = .envir)
}

# transforming to global time coordinates
phonetic_coding_sync <- phonetic_coding %>%
  mutate(name = name_wav) %>%
  left_join(sync_curr %>% select(name,s1,t1,s2,t2,Fs), join_by(xmin_s > s1, xmax_s < s2, name==name)) %>%
  mutate(onset = (xmin*Fs - s1)*(t2-t1)/(s2-s1) + t1) %>%
  mutate(duration = xmax - xmin) %>%
  select(onset,duration,subject_id,session_id,task_id,run_id,trial_id,phrase_id,tier_name,text,annotation_num)
phonetic_coding_sync$phrase_id[is.na(phonetic_coding_sync$phrase_id)] <- 0

# #creating produced phonemes table
# phonetic_coding_phonemes <- phonetic_coding_sync %>%
#   filter(tier_name=='phone') %>%
#   rename(phoneme=text,phoneme_id=annotation_num) %>%
#   select(-tier_name) %>%
#   write_tsv(glue('{PATH_ANNOT}/sub-{SUBJECT}_ses-{SESSION}_task-{TASK}_annot-produced-phonemes.tsv'))

# #creating produced words table
# phonetic_coding_words <- phonetic_coding_sync %>%
#   filter(tier_name=='word') %>%
#   rename(word=text,word_id=annotation_num)%>%
#   select(-tier_name) %>%
#   write_tsv(glue('{PATH_ANNOT}/sub-{SUBJECT}_ses-{SESSION}_task-{TASK}_annot-produced-words.tsv'))

#creating produced sentences table
phonetic_coding_sentences <- phonetic_coding_sync %>%
  filter(tier_name=='utterance') %>%
  filter(text!="") %>%
  rename(utterance=text,utterance_id=annotation_num)%>%
  select(-tier_name) %>%
  write_tsv(glue('{PATH_ANNOT}/sub-{SUBJECT}_ses-{SESSION}_task-{TASK}_annot-produced-utterances-whisperx.tsv'))

# #creating produced sentence table
# phonetic_coding_sentence <- phonetic_coding_words %>%
#   filter(word!='sp') %>%
#   mutate(offset = onset + duration) %>%
#   group_by(sub,ses,task,run_id,trial_id) %>%
#   summarise(onset=min(onset),duration=max(offset) - min(onset), 
#             sentence=paste(word, collapse=' '),
#             complete_type = paste(unique(complete_type),collapse='')) %>%
#   relocate(onset,duration) %>%
#   write_tsv(glue('{PATH_ANNOT}/sub-{SUBJECT}_ses-{SESSION}_task-{TASK}_annot-produced-sentences.tsv'))