%%%% this script cuts the audio recording from intrap SMSL into a separate audio file for every trial
%%%% .... these audio files will be individually scored for sequencing/epenthesis errors
% trials will be cut starting at audio_go_offset and ending at keypress_time
%%% data on resulting trial audio files stored in PATH_TRIAL_AUDIO/###_audiofiles.tsv

% note: occasionally keypress_time is recorded incorrectly, as earlier than it should have been (or marked NaN)...
% ... in this case, in the output table, keypress_time will be replaced by NaN and....
% ... trial offset time and duration will be determined by corrected_trial_end_prestim


% Andrew Meier

clear

%% manual editing of specific times
% use this table to modify the start/stop times of specific trials
%   ... for use on specific trials when sub starts early or trigger was sent at wrong time
%%%% indicate which trials to modify with 'inds'
%%%% for each ind, use 'sec' to specify how much to move start times [col 1] and stop times [col 2]
%%%%%%%% negative values move timepoints earlier, positive move later
trials_to_modify_inds = []; 
    trials_to_modify_sec = []; 

%% parameters and paths
%%% audio cutting parameters
% wait the following number of seconds after beep GO signal to start response window
% ... this is mostly to eliminate the echo that is still audible after recorded beep offset time
% ... 0.2 can be helpful, but in some subjects (eg 1024) it will cut off early vocal responses
clip_seconds_postbeep = 0; 
% add the following number of seconds after keypress to the duration of STOP trial audio...
% ... this extra time is intended to make it easier to determine whether the utterance was completed or interrupted
stoptrial_post_keypress_extension_sec = 1; 
% for trials with negative or NaN duration, set trial end as this many seconds before 
% the following trial's visual_onset
corrected_trial_end_prestim = 1; % visual onset is programmed to happen 1s after keypress
%%%% if erroneous trial is the last trial, use the following value in seconds as trial duration
last_trial_duration_correction = 2.5; 
%  if trialdur is longer than the following length in seconds, cut it down to this length
%   ...... also give a warning and open it in praat
max_trialdur = 6; 

%%% Defining paths
SUBJECT='DM####';
SESSION = 'intraop';
TASK = 'smsl'; 
RUN = '##'; % see Y:\DBS\derivatives\sub-DM####\annot\sub-DM####_runs.tsv... look at duration and comment columns
PATH_DATASET = 'Y:\DBS';
open_problematic_trials_in_praat = 1; % automatically load trials with unexpected durations in praat to check
start_offset = 25; % optional - set the bounds for session audio outside of trials' boundaries
end_offset = 5; % optional - set the bounds for session audio outside of trials' boundaries

PATH_DER = [PATH_DATASET filesep 'derivatives'];
PATH_DER_SUB = [PATH_DER filesep 'sub-' SUBJECT];  
PATH_PREPROC = [PATH_DER_SUB filesep 'preproc'];
PATH_ANNOT = [PATH_DER_SUB filesep 'annot'];
PATH_AEC = [PATH_DER_SUB filesep 'aec']; 
PATH_SCORING = [PATH_DER_SUB filesep 'analysis' filesep 'task-', TASK, '_scoring'];
PATH_ANALYSIS = [PATH_DER_SUB filesep 'analysis'];
PATH_TRIAL_AUDIO = [PATH_ANALYSIS filesep 'task-', TASK, '_trial-audio'];
PATH_TRIAL_AUDIO_SESSION_GO = [PATH_TRIAL_AUDIO, filesep, 'ses-', SESSION, '_go-trials'];
PATH_TRIAL_AUDIO_SESSION_STOP = [PATH_TRIAL_AUDIO, filesep, 'ses-', SESSION, '_stop-trials']; 

PATH_SRC = [PATH_DATASET filesep 'sourcedata'];
PATH_SRC_SUB = [PATH_SRC filesep 'sub-' SUBJECT];  
PATH_SRC_SESS = [PATH_SRC_SUB filesep 'ses-' SESSION]; 
PATH_AUDIO = [PATH_SRC_SESS filesep 'audio']; 
PATH_TASK = [PATH_SRC_SESS, filesep, 'task'];
sync_file = [PATH_ANNOT filesep 'sub-' SUBJECT, '_ses-', SESSION,  '_sync-audio-all.tsv']; 
landmarks_file = [PATH_ANNOT filesep 'sub-' SUBJECT, '_ses-', SESSION,  '_annot-audio-landmarks.tsv']; 

%%%% use either AEC or raw audio
audiofile_full_session = [PATH_AUDIO filesep 'sub-' SUBJECT, '_ses-', SESSION, '_task-' TASK, '_run-', RUN, '_recording-directionalmic_physio.wav']; %pre-aec


%% load audio and timing info
% make trial audio directories
system(['mkdir ' PATH_TRIAL_AUDIO]);
system(['mkdir ' PATH_TRIAL_AUDIO_SESSION_GO]);
system(['mkdir ' PATH_TRIAL_AUDIO_SESSION_STOP]); 

% load and organize audio and trial data
trials = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_task-', TASK, '_annot-trials.tsv']);

% sometimes (e.g. DM1008 smsl training) a 'trial zero' gets added to the trial list... remove this trial and correct other indices
if trials.trial_id(1) == 0 
    trials = trials(2:end, :); % remove trial zero
    trials.id = trials.id + 1; % update other indices
end
    
% add vars to trial table
ntrials = height(trials);
cellcol = cell(ntrials,1);
nancol = NaN(ntrials,1); 
trials = [trials, table(cellcol,cellcol,nancol, 'VariableNames',...
                         {'dir','wavname','local_file_id'} )];

audiofiles = trials(:,{'run_id','block_id','trial_id','stim_id','word','local_file_id','stim_condition','visual_onset','audio_go_offset','keypress_time','is_stoptrial'});
audiofiles.starts = audiofiles.audio_go_offset + clip_seconds_postbeep; % clip beginning of trial, which may have beep echo
audiofiles.ends = audiofiles.keypress_time;
audiofiles = bml_annot_table(audiofiles);
audiofiles.id = []; 
[~, audiofilename] = fileparts(audiofile_full_session);

%%% get sample vs. audio time info
audinfo = audioinfo(audiofile_full_session); 
[fpath, fname, fext] = fileparts(audiofile_full_session); 
if exist(sync_file,'file') % if protocol A04 syncing has been done
    sync = bml_annot_read_tsv(sync_file);
    roi = sync(contains(sync.name, audiofilename,'IgnoreCase',true),:); % ignore case because sometimes 'SMSL' gets changed to 'smsl'
    audiofiles.starts_file_relative = audiofiles.starts - roi.starts(1); % trial start time relative to full recording beginning
    audiofiles = movevars(audiofiles, 'starts_file_relative', 'Before', 'run_id');
elseif ~exist(sync_file,'file') % if protocol A04 syncing has not yet been done
    roi = table;
    roi.name{1} = [fname, fext]; 
    roi.folder{1} = PATH_AUDIO;
    roi.chantype{1} = 'directionalmic';
    roi.filetype{1} = 'audio.wav';
    roi.Fs(1) = audinfo.SampleRate; 
    roi.nSamples(1) = audinfo.TotalSamples; 
    % get approximate landmark sample using time and sample rate
    landmarks = bml_annot_read_tsv(landmarks_file);
    roi.s1(1) = 1; 
    roi.t1(1) = landmarks.dtc - landmarks.ftc; % approx. audio file start time in global time coords
    roi.s2(1) = landmarks.ftc(1) * audinfo.SampleRate; % approx. sample at ftc
    roi.t2(1) = landmarks.dtc(1); % global time coords 
end

% cut audio trials
%%%% load session audio
% this section can be run without having done syncing in protocol A04
cfg = [];
cfg.roi = roi; 
cfg.roi.starts(1) = min(trials.starts) - start_offset; % set the bounds for session audio outside of trials' boundaries
cfg.roi.ends(1) = max(trials.ends) + end_offset; % set the bounds for session audio outside of trials' boundaries
session_aud = bml_load_continuous(cfg);
check_this_trial = 0; % initialize; if this var == 1, open the trial in praat

for itrial = 1:ntrials
    if audiofiles.is_stoptrial(itrial)
        audiofiles.dir{itrial} = PATH_TRIAL_AUDIO_SESSION_STOP; 
        post_keypress_extension_sec = stoptrial_post_keypress_extension_sec; % add buffer to stop trials for scoring
    elseif ~audiofiles.is_stoptrial(itrial)
        audiofiles.dir{itrial} = PATH_TRIAL_AUDIO_SESSION_GO; 
        post_keypress_extension_sec = 0; 
    end
    audiofiles.filename{itrial} = ['trial', num2str(itrial,'%03.f'), '_', trials.word{itrial}, '.wav'];
    trialdur = audiofiles.duration(itrial);
    
    if isnan(trialdur) || trialdur <= 0 % if erroneous keypress_time
        audiofiles.keypress_time(itrial) = NaN; % erase erroneous time
        % for trials with negative or NaN duration, set trial end as X many seconds before 
        % the following trial's visual_onset
        if itrial < ntrials % if not last trial
            new_endtime = audiofiles.visual_onset(itrial+1) - corrected_trial_end_prestim; 
        % if erroneous trial is the last trial, use X seconds as trial duration
        elseif itrial == ntrials
            new_endtime = audiofiles.starts(itrial) + last_trial_duration_correction;
        end
        audiofiles.ends(itrial) = new_endtime; clear new_endtime
        audiofiles.duration(itrial) = audiofiles.ends(itrial) - audiofiles.starts(itrial); % correct the duration
        check_this_trial = 1; 
        warning([audiofiles.filename{itrial}, ' duration was originally erroneously labeled as ',...
            num2str(trialdur), ' seconds. Replacing with ', num2str(audiofiles.duration(itrial)),...
            'sec. See resulting audio file in Praat.'])
    elseif trialdur > max_trialdur
        audiofiles.ends(itrial) = audiofiles.starts(itrial) + max_trialdur; % cut trial duration
        audiofiles.duration(itrial) = max_trialdur; % correct the duration
        check_this_trial = 1; 
        warning([audiofiles.filename{itrial}, ' duration is unexpectedly long (',...
            num2str(trialdur), ' seconds). Cutting down to ', num2str(max_trialdur), ' seconds. ', ...
            'See resulting audio file in Praat.'])
    end

    % apply manual edits to trialtimes
    if ~isempty(trials_to_modify_inds) && any(itrial == trials_to_modify_inds)
        matchrow = itrial == trials_to_modify_inds;
        start_time_edit = trials_to_modify_sec(matchrow, 1);
        end_time_edit = trials_to_modify_sec(matchrow, 2);
        audiofiles.starts(itrial) = audiofiles.starts(itrial) + start_time_edit; 
        audiofiles.ends(itrial) = audiofiles.ends(itrial) + end_time_edit; 
    end

    cfg=[];
    cfg.epoch=audiofiles(itrial,:);
    cfg.epoch.ends(1) = cfg.epoch.ends(1) + post_keypress_extension_sec; % stoptrial buffer
    cfg.epoch.duration(1) = cfg.epoch.duration(1) + post_keypress_extension_sec; % stoptrial buffer
    thistrial = bml_redefinetrial(cfg,session_aud);
    fs = thistrial.fsample; 
    trialaud = thistrial.trial{1};
    % make sure columns are channels, rows are data points
    [~, longdim] = max(size(trialaud));         [~, shortdim] = min(size(trialaud)); 
    trialaud = permute(trialaud, [longdim, shortdim]); 

    audiowrite([audiofiles.dir{itrial} filesep audiofiles.filename{itrial}], trialaud, fs)
    if check_this_trial && open_problematic_trials_in_praat % load created audio file in praat
        cmd = ['praat --open ',  [audiofiles.dir{itrial} filesep audiofiles.filename{itrial}], ' &'];
        system(cmd);
        check_this_trial = 0; 
    end
end
audiofiles = movevars(audiofiles,{'trial_id','filename','local_file_id','word',},'Before','starts'); 

% separate into GO and STOP trials, then save all 3 file lists
audiofiles_gotrials = audiofiles(~logical(audiofiles.is_stoptrial),:); 
    audiofiles_gotrials.local_file_id = [1:height(audiofiles_gotrials)]'; % idx of .wav file within its folder
audiofiles_stoptrials = audiofiles(logical(audiofiles.is_stoptrial),:); 
    audiofiles_stoptrials.local_file_id = [1:height(audiofiles_stoptrials)]'; % idx of .wav file within its folder

file_prepend = ['sub-' SUBJECT, '_ses-', SESSION, '_task-' TASK, '_']; 
writetable(audiofiles, [PATH_TRIAL_AUDIO, filesep, file_prepend, 'audiofiles.tsv'], 'FileType','text', 'Delimiter','\t')
writetable(audiofiles_gotrials, [PATH_TRIAL_AUDIO, filesep, file_prepend, 'audiofiles_gotrials.tsv'], 'FileType','text', 'Delimiter','\t')
writetable(audiofiles_stoptrials, [PATH_TRIAL_AUDIO, filesep, file_prepend, 'audiofiles_stoptrials.tsv'], 'FileType','text', 'Delimiter','\t')

