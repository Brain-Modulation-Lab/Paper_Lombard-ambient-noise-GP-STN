library(tidyverse)

# Defining paths
SUBJECT <-'DM####'
PATH_DATASET <- 'Y:/DBS'

PATH_DER <- file.path(PATH_DATASET, 'derivatives')
PATH_DER_SUB <- file.path(PATH_DER, paste0('sub-',SUBJECT))  
PATH_PREPROC <- file.path(PATH_DER_SUB, 'preproc')
PATH_ANNOT <- file.path(PATH_DER_SUB, 'annot')
PATH_SRC <- file.path(PATH_DATASET, 'sourcedata')
PATH_SRC_SUB <- file.path(PATH_SRC, paste0('sub-',SUBJECT)) 
PATHS_TASK <- file.path(PATH_SRC_SUB,c('ses-training','ses-preop','ses-intraop'),'task')

setwd(PATH_PREPROC)

tasks <- 
  tibble(path = PATHS_TASK) %>%
  mutate(path_exists=file.exists(path)) %>%
  filter(path_exists) %>%
  rowwise() %>%
  group_modify(~file.info(list.files(.x$path,pattern=glob2rx('*_events.tsv'),full.names=TRUE))) %>%
  rownames_to_column('file_name') 

hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0")
        ,formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0")
        ,formatC(t %% 60, width = 2, format = "d", flag = "0")
        ,sep = ":")}

tasks_file_name <- tasks$file_name
names(tasks_file_name) <- tasks$file_name
events_files <- tasks_file_name %>% 
  map_df(read_tsv, .id = 'tasks_file_name') %>%
  group_by(tasks_file_name) %>%
  summarise(duration= max(onset) - min(onset),onset = min(onset)) %>%
  ungroup() %>%
  mutate(time_starts = hms(onset), time_ends=hms(onset + duration)) %>%
  mutate(path = dirname(tasks_file_name), filename = basename(tasks_file_name)) %>%
  select(onset,duration,path,filename,time_starts,time_ends)

write_tsv(events_files,paste0(PATH_ANNOT,'/sub-',SUBJECT,'_events-files.tsv'))
