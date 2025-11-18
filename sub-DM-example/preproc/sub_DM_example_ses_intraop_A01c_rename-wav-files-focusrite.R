library(dplyr)
library(tidyverse)
library(glue)
library(tuneR)

SUBJECT <- 'DM10XX'
SESSION <- 'intraop'
#PATH_DATASET <- 'Y:/DBS'
PATH_DATASET <- '/Volumes/Nexus4/DBS/'
PATH_ANNOT <- glue("{PATH_DATASET}/derivatives/sub-{SUBJECT}/annot")
PATH_PREPROC <- glue("{PATH_DATASET}/derivatives/sub-{SUBJECT}/preproc")
PATH_FIGURES <- glue("{PATH_PREPROC}/figures")
PATH_AUDIO <- glue("{PATH_DATASET}/sourcedata/sub-{SUBJECT}/ses-{SESSION}/audio")
PATH_FOCUSRITE <- glue("{PATH_AUDIO}/focusrite")

# MANUALLY CREATE MAPPING between channel id/idx and the recording name ----  
# This will RARELY change. If all goes well in the OR, they should stay as-is
CHANNELS_TABLE <-
  tribble(
    ~focusrite_file_type, ~channel_id, ~recording,
    ".wav",                1,           "recording-directionalmic_physio",
    ".wav",                2,           "recording-ambientmic_physio",
    ".wav",                3,           "events",
    ".wav",                4,           "recording-speaker_stim"
  )
CHANNELS_TABLE


# Loop over focusrite files and write them to disc with new names ----

# Check that second tsv, called zoom_channels table, has columns zoom_file_type, channel_id, recording
if(!all(c("focusrite_file_type", "channel_id", "recording") %in% colnames(CHANNELS_TABLE))){
  stop("The zoom channels table does not have all the required columns.")
}

for(f in dir(PATH_FOCUSRITE,pattern=glob2rx("*.wav"),full.names = FALSE)){
  # Open a second nested loop for each row of the zoom_channels table,
  for(j in 1:nrow(CHANNELS_TABLE)){
    
    # create a filename for current_wav as sub-<subject_id>_ses-<session_id>_task-<task-id>_run-<run_id>_recording-<recording>.wav
    new_filename <- str_replace(f, "(.*)(\\.wav)$", str_c("\\1_", CHANNELS_TABLE$recording[j], "\\2"))
    print(new_filename)
    
    # using sox to convert to 16-bit signed-integer
    system(glue("sox {PATH_FOCUSRITE}/{f} -c 1 -b 16 -e signed-integer {PATH_AUDIO}/{new_filename} remix {CHANNELS_TABLE$channel_id[j]}"))
    
    # load the created wav file,
    wav_data <- readWave(glue("{PATH_AUDIO}/{new_filename}"))
    
    # Calculate the number of samples that are saturating the range for current_wav
    saturation_count <- sum(abs(wav_data@left) >= 2^(wav_data@bit-1))
    
    # plot an oscillogram of current_wav and save it in a figures folder within the pwd
    png(glue("{PATH_FIGURES}/{new_filename}.png"))
    plot(wav_data)
    title(glue("{saturation_count} samples ({100*signif(saturation_count/length(wav_data),2)}%) clipped" ))
    dev.off()
  }
}

