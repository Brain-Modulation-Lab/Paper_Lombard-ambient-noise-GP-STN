library(tidyverse)

# Defining paths
SUBJECT <-'DM####'
PATH_DATASET <- '/Volumes/Nexus4/DBS' # Y:/DBS
library(tidyverse)

PATH_DER <- file.path(PATH_DATASET, 'derivatives')
PATH_DER_SUB <- file.path(PATH_DER, paste0('sub-',SUBJECT))  
PATH_PREPROC <- file.path(PATH_DER_SUB, 'preproc')
PATH_ANNOT <- file.path(PATH_DER_SUB, 'annot')
PATH_SRC <- file.path(PATH_DATASET, 'sourcedata')
PATH_SRC_SUB <- file.path(PATH_SRC, paste0('sub-',SUBJECT)) 
PATHS_TASK <- file.path(PATH_SRC_SUB,c('ses-training','ses-preop','ses-intraop'),'task')

setwd(PATH_PREPROC)

# Custom session ordering for arranging rows in tables
SESSION_ORDER <- c("training", "preop", "intraop")

# how many days before surgery was each session
# utilized to set GTC
SESSION_DAY <- tibble(session=     c("training", "preop", "intraop"), 
                      session_day     =c(-XX,         0,       0)) 


# Load task files and compile into events files table
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

write_tsv(events_files, paste0(PATH_ANNOT,'/sub-',SUBJECT,'_events-files.tsv'))



# Create runs & sessions table ---- 

runs <- 
  events_files %>%
  
  # Extract session, task, run from filename using regex
  mutate(
    session = str_extract(filename, "ses-[a-zA-Z]+") %>% str_remove("ses-"),
    task    = str_extract(filename, "task-[a-zA-Z]+") %>% str_remove("task-"),
    run     = str_extract(filename, "run-\\d+") %>% str_remove("run-") %>% as.integer()
  ) %>%
  
  left_join(SESSION_DAY, by="session") %>% 
  mutate(
    onset = onset + 86400*session_day
  ) %>%
  
  # Add required manual columns, filled with 'n/a'
  mutate(
    audio       = "n/a",
    trellis     = "datafile000X",
    MER_side    = "n/a",
    MER_track   = "n/a",
    MER_depth   = "n/a",
    DBS_left    = "n/a",
    DBS_right   = "n/a",
    comment     = "n/a"
  ) %>%
  
  # Set session factor level order
  mutate(session = factor(session, levels = SESSION_ORDER, ordered = TRUE)) %>%
  
  # Arrange by session, then task, then run
  arrange(session, (onset)) %>%
  
  # Select final columns in specified order
  select(
    onset, duration, session, task, run, session_day, time_starts, time_ends,
    audio, trellis, MER_side, MER_track, MER_depth,
    DBS_left, DBS_right, comment
  )
write_tsv(runs, paste0(PATH_ANNOT,'/sub-',SUBJECT,'_runs-orig.tsv'))


sessions <-
  runs %>%
  mutate(offset = onset+duration) %>%
  group_by(session, session_day) %>% 
  summarize(onset = min(onset), 
            offset = max(offset), 
            time_starts = min(time_starts), 
            time_ends = max(time_ends)) %>%
  mutate(duration = offset - onset) %>%
  relocate(onset, duration, offset, session)
write_tsv(sessions, paste0(PATH_ANNOT,'/sub-',SUBJECT,'_sessions-orig.tsv'))




