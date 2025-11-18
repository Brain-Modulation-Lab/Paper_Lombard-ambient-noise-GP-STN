
%% loading toolbox
bml_defaults
ft_defaults
format long
clear

%% Defining paths
SUBJECT='DM1032';
SESSION='intraop';
TASKS = {'sync', 'lombard', 'smsl'}; 
%PATH_DATASET = '/Volumes/Nexus4/DBS';
PATH_DATASET = 'Y:\DBS';

PATH_DER = [PATH_DATASET filesep 'derivatives'];
PATH_DER_SUB = [PATH_DER filesep 'sub-' SUBJECT];  
PATH_PREPROC = [PATH_DER_SUB filesep 'preproc'];
PATH_ANNOT = [PATH_DER_SUB filesep 'annot'];

PATH_SRC = [PATH_DATASET filesep 'sourcedata'];
PATH_SRC_SUB = [PATH_SRC filesep 'sub-' SUBJECT];   
PATH_RIPPLE = [PATH_SRC_SUB filesep 'ses-intraop' filesep 'ripple'];
PATH_ALHAOMEGA = [PATH_SRC_SUB filesep 'ses-intraop' filesep 'alphaomega']; %path to neuroomega's mpx files.
PATH_ALHAOMEGA_MAT = [PATH_ALHAOMEGA]; %path to neuroomega's converted mat files.
PATH_AUDIO_DER = [PATH_DER_SUB filesep 'aec']; %path to acoustic echo cancelling folder
PATHS_AUDIO_SRC = strcat(PATH_SRC_SUB,filesep,{'ses-training';'ses-preop';'ses-intraop'},filesep,'audio');
PATHS_TASK = strcat(PATH_SRC_SUB,filesep,{'ses-training';'ses-preop';'ses-intraop'},filesep,'task');

% loading annotation tables
events_files = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_events-files.tsv']);
events_files.path = strrep(events_files.path,'/Volumes/Nexus4','Y:');
session = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_sessions.tsv']);
runs = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_runs.tsv']);

cd(PATH_PREPROC)

%% Loading header information
cfg=[];
cfg.path=PATH_RIPPLE;
cfg.pattern='*.nev';
cfg.filetype='trellis.nev';
info_nev=bml_info_raw(cfg);
info_nev.basename = strrep(info_nev.name,'.nev','');

cfg=[];
cfg.path=PATH_RIPPLE;
cfg.pattern='*.ns5';
cfg.filetype='trellis.ns5';
info_ns5=bml_info_raw(cfg);
info_ns5.basename = strrep(info_ns5.name,'.ns5','');

cfg=[];
cfg.path=PATH_RIPPLE;
cfg.pattern='*.nf3';
cfg.filetype='trellis.nf3';
info_nf3=bml_info_raw(cfg);
info_nf3.chantype(:) = {'hires'};
info_nf3.basename = strrep(info_nf3.name,'.nf3','');

cfg=[];
cfg.path=PATH_RIPPLE;
cfg.pattern='*.nf6';
cfg.filetype='trellis.nf6';
info_nf6=bml_info_raw(cfg);
info_nf6.chantype(:) = {'hifreq'};
info_nf6.basename = strrep(info_nf6.name,'.nf6','');

% cfg=[];
% cfg.mpx_path=PATH_ALHAOMEGA;
% cfg.path=PATH_ALHAOMEGA_MAT;
% cfg.filetype='neuroomega.mat';
% cfg.chantype='analog';
% info_neuroomega=bml_neuroomega_info_raw(cfg);
% info_neuroomega.basename = strrep(info_neuroomega.name,'.mat','');

info_wav_src_renamed=[];
for i=1:length(PATHS_AUDIO_SRC)
    cfg=[];
    cfg.path=PATHS_AUDIO_SRC{i};
    cfg.filetype='audio.wav';
    cfg.pattern='sub-*.wav';
    info_wav_src_renamed=bml_annot_rowbind(info_wav_src_renamed,bml_info_raw(cfg));
end
info_wav_src_renamed.basename = strrep(info_wav_src_renamed.name,'.wav','');

cfg=[];
cfg.path=PATH_AUDIO_DER;
cfg.pattern='*.wav';
cfg.filetype='audio.wav';
info_wav=bml_info_raw(cfg);
info_wav.basename = strrep(info_wav.name,'.wav','');
info_wav = info_wav(~ismember(info_wav.name,info_wav_src_renamed.name),:);
info_wav = bml_annot_rowbind(info_wav_src_renamed,info_wav);

%Saving files information table with OS times
%info_orig = bml_annot_rowbind(info_ns5,info_nf3,info_nf6,info_nev,info_neuroomega,info_wav,info_wav_src);
info_orig = bml_annot_rowbind(info_ns5,info_nf3,info_nf6,info_nev,info_wav,info_wav_src_renamed);
bml_annot_write_tsv(info_orig, [PATH_ANNOT filesep 'sub-' SUBJECT '_data-files-orig.tsv']);

%% ALIGNING TRELLIS FILES TO EVENTS BY SCRIPTS (NO TIME WARPING)

%loading all events from event files
assert(~isempty(events_files))
task_events = table();
for j = 1:height(events_files)
    events_tmp = bml_annot_read_tsv(fullfile(events_files.path{j},events_files.filename{j})); 
    if isnumeric(events_tmp.stim_file) && all(isnan(events_tmp.stim_file)); 
        events_tmp.stim_file = repmat({''}, height(events_tmp), 1); 
    end
    events_tmp.filename(:)=events_files.filename(j);
    task_events = bml_annot_rowbind(task_events, events_tmp);
end

%looping through nev files
info_nev.delta_t(:) = nan;
info_nev.min_cost(:) = nan;
info_nev.n(:) = nan;
for i=1:height(info_nev)
    %selecting events files that correspond to this nev file. Doing this by
    %coarse time
    
    %loading events from nev fil
    cfg=[];
    cfg.roi = bml_roi_table(info_nev(i,:));
    nev = bml_annot_read_event_nev(cfg);

    %matching events by similarity
    cfg=[];
    cfg.timetol=0.001;
    [idxs_task_events, idxs_nev] = bml_sync_match_events(cfg, task_events, nev);   

    %Aligning blackrock to presenation task to get GTC. No timewarping.
    cfg=[];
    cfg.scan = 100;
    [slave_delta_t, min_cost, warpfactor] = bml_timealign_annot(cfg, task_events(idxs_task_events,:), nev);
    
    info_nev.delta_t(i)=slave_delta_t;
    info_nev.min_cost(i)=min_cost;
    info_nev.n(i)=height(nev);
        
    fprintf('%f sec, min cost = %f\n', slave_delta_t, min_cost)

end

%adjusting ripple files
info_nev.delta_t(ismissing(info_nev.delta_t)) = 0;
cfg.keys = 'basename';
cinfo_ns5 = bml_annot_left_join(cfg,info_ns5,info_nev(:,{'basename','delta_t'}));
cinfo_nf3 = bml_annot_left_join(cfg,info_nf3,info_nev(:,{'basename','delta_t'}));
cinfo_nf6 = bml_annot_left_join(cfg,info_nf6,info_nev(:,{'basename','delta_t'}));
cinfo_nev = info_nev;
cinfo_nev.starts = cinfo_nev.starts + cinfo_nev.delta_t;
cinfo_nev.ends = cinfo_nev.ends + cinfo_nev.delta_t;
cinfo_ns5.starts = cinfo_ns5.starts + cinfo_ns5.delta_t;
cinfo_ns5.ends = cinfo_ns5.ends + cinfo_ns5.delta_t;
cinfo_nf3.starts = cinfo_nf3.starts + cinfo_nf3.delta_t;
cinfo_nf3.ends = cinfo_nf3.ends + cinfo_nf3.delta_t;
cinfo_nf6.starts = cinfo_nf6.starts + cinfo_nf6.delta_t;
cinfo_nf6.ends = cinfo_nf6.ends + cinfo_nf6.delta_t;
roi_nev = bml_roi_table(cinfo_nev,'trellis');
roi_ns5 = bml_roi_table(cinfo_ns5,'trellis');
roi_nf3 = bml_roi_table(cinfo_nf3,'trellis');
roi_nf6 = bml_roi_table(cinfo_nf6,'trellis');
roi_ripple = bml_annot_rowbind(roi_nev, roi_ns5, roi_nf3, roi_nf6);
%ripple time after solid translation to GTC is master time
sync_ripple = roi_ripple; 
bml_annot_write_tsv(sync_ripple, [PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION  '_sync-ripple.tsv'])

% %% DETECT  STRETCHES CONSTANT DEPTH IN ALPHA OMEGA
% cfg=[];
% cfg.criterion = @(x) length(unique(x.depth))==1 && abs((max(x.ends)-min(x.starts))-sum(x.duration))<.001;
% neuro_cons_depth = bml_annot_consolidate(cfg,roi_neuroomega);
% ao_session = neuro_cons_depth(neuro_cons_depth.duration > 30,:);
% ao_session = ao_session(:,{'starts','ends','depth'});
% bml_annot_write_tsv(ao_session,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_stable-MER.tsv']);

%% LOAD MASTER EVENTS
%loading master events from all nev files
cfg=[];
cfg.roi = roi_nev;
master_events = bml_annot_read_event_nev(cfg);

% SYNCHRONIZING EVENT FILES TO RIPPLE EVENTS
% Will inherit exact timing from ripple. looping through events_files
events_files.delta_t(:)=nan;
events_files.min_cost(:)=nan;
events_files.warpfactor(:)=nan;
events_files.n(:)=nan;
for i=1:height(events_files)
    %selecting nev files that correspond to this events_files. 
    nev_events_files = bml_annot_filter(roi_nev,bml_annot_extend(events_files(i,:),events_files.duration(i) .* 1e-3));
    
    if isempty(nev_events_files); continue; end

    task_events = bml_annot_read_tsv(fullfile(events_files.path{i},events_files.filename{i}));
  
    %matching task_events to master_events
    cfg=[];
    [idxs_master_events,idxs_task_events]=bml_sync_match_events(cfg,master_events,task_events);    
    
    %Aligning events to ripple
    cfg=[];
    cfg.scan = 1;
    cfg.timewarp = true;
    [slave_delta_t, min_cost, warpfactor] = bml_timealign_annot(cfg, master_events(idxs_master_events,:), task_events(idxs_task_events,:));   
    events_files.delta_t(i)=slave_delta_t;
    events_files.min_cost(i)=min_cost;
    events_files.warpfactor(i)=warpfactor;
    events_files.n(i)=height(task_events(idxs_task_events,:));;
    
    %applying timealignment to events as fallback strategy
    tbar = mean(task_events.starts);
    task_events.starts_corrected = (task_events.starts - tbar) .* warpfactor + tbar + slave_delta_t;
        
    %getting time from master_events according to alignemnt (when available)
    task_events.starts = bml_map(1:height(task_events),idxs_task_events,master_events.starts(idxs_master_events),nan);

    %using fallback times when nev times not available
    task_events.starts(isnan(task_events.starts))=task_events.starts_corrected(isnan(task_events.starts));

    bml_annot_write_tsv(task_events,[PATH_ANNOT filesep strrep(events_files.filename{i},'.tsv','-sync.tsv')]);
    fprintf('sync %s, delta_t=%f, min cost=%f, warpfactor=%f \n',events_files.filename{i}, slave_delta_t, min_cost, warpfactor)

end

%% SYNCHRONIZE ZOOM WITH EVENTS USING DIGITAL TRIGGERS

% converting ripple events to output sync trigger values
master_trig = master_events;
master_trig.tmp = [2^15;abs(diff(master_trig.value))];
master_trig = master_trig(master_trig.tmp >= 2^9,:);
master_trig.value = (master_trig.value > 2^9) .* 1.0;

% doing coarse alignment wav file by file (can't assume chronological order in wav_events_files) 
info_wav_events = info_wav_src_renamed(contains(info_wav_src_renamed.name,'_events.wav') & contains(info_wav_src_renamed.name,TASKS),:);
info_wav_events.mean_sim(:)=nan;
info_wav_events.mean_delta_t(:)=nan;
info_wav_events.median_delta_t(:)=nan;
info_wav_events.std_mismatch(:)=nan;
info_wav_events.mad_mismatch(:)=nan;
for i=1:height(info_wav_events)
    wav_events = bml_annot_read_event_wav(info_wav_events(i,:));
    [idxs_master_trig,idxs_wav_events, mean_sim]=bml_sync_match_events([],master_trig,wav_events); 
    dtv = master_trig.starts(idxs_master_trig)-wav_events.starts(idxs_wav_events);
    info_wav_events.mean_sim(i) = mean_sim;
    info_wav_events.mean_delta_t(i) = mean(dtv);
    info_wav_events.median_delta_t(i) = median(dtv);
    info_wav_events.std_mismatch(i) = std(dtv);
    info_wav_events.mad_mismatch(i) = mad(dtv,1);
end
cinfo_wav_events = info_wav_events; 
cinfo_wav_events.starts = cinfo_wav_events.starts + cinfo_wav_events.median_delta_t;
cinfo_wav_events.ends = cinfo_wav_events.starts + cinfo_wav_events.duration;


% SYNCHRONIZE ZOOM WITH EVENTS USING DIGITAL TRIGGERS
cfg =[];
cfg.master_events = master_trig;
cfg.roi = bml_roi_table(cinfo_wav_events);
cfg.roi = bml_annot_intersect('x',cfg.roi,roi_nev); 
cfg.roi = cfg.roi(cfg.roi.duration > 30,:); 
cfg.timewarp = true;
cfg.diagnostic_plot = true;
cfg.timetol = 0.001;
cfg.coarsetol = 0.001;
cfg.restrict_master_by = 'file_id';
sync_wav = bml_sync_audio_event(cfg);
bml_annot_write_tsv(sync_wav,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_sync-audio.tsv']);

% Transfering synchronizatino to other wav files
sync_audio_cons = sync_wav(contains(sync_wav.name,SESSION) & contains(sync_wav.name,TASKS),:);
cinfo_wav = info_wav(contains(info_wav.name,SESSION) & contains(info_wav.name,TASKS),:);

% matching files by basename
tmp = regexp(sync_audio_cons.name,'^(sub-[\w]*_ses-[\w]*_task-[\w]*_run-[0-9]*)_[\w-\.]*$','tokens','once','forceCellOutput');
sync_audio_cons.basename  = [tmp{:}]';
tmp = regexp(cinfo_wav.name,'^(sub-[\w]*_ses-[\w]*_task-[\w]*_run-[0-9]*)_[\w-\.]*$','tokens','once','forceCellOutput');
cinfo_wav.basename  = [tmp{:}]';
tmp=regexp(cinfo_wav.name,'^sub-[\w]*_ses-[\w]*_task-[\w]*_run-[0-9]*_(recording-)?([a-zA-Z0-9-]*)(_[\w]*)?\.wav$','tokens','once');
cinfo_wav.chantype=cellfun(@(x) x(2),tmp);

% transfering sync to other wav files
sync_audio_all = table();
for i=1:height(sync_audio_cons)
    sel_info_wav = cinfo_wav(strcmp(cinfo_wav.basename,sync_audio_cons.basename(i)),:);
    sync_audio_i = repmat(sync_audio_cons(i,:),[height(sel_info_wav),1]);
    sync_audio_i.name = sel_info_wav.name;
    sync_audio_i.folder = sel_info_wav.folder;
    sync_audio_i.nSamples = sel_info_wav.nSamples;
    sync_audio_i.chantype = sel_info_wav.chantype;
    sync_audio_all = bml_annot_rowbind(sync_audio_all,sync_audio_i);
end
bml_annot_write_tsv(sync_audio_all,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_sync-audio-all.tsv']);

% % SYNCHRONIZE NEUROOMEGA FILE WITH EVENTS USING DIGITAL TRIGGERS
% cfg = [];
% cfg.master_events = master_events;
% cfg.roi = roi_neuroomega;
% cfg.timewarp = true;
% cfg.diagnostic_plot = true;
% cfg.timetol = 0.0003;
% cfg.coarsetol = 0.001;
% sync_ao = bml_sync_neuroomega_event(cfg);
% 
% bml_annot_write_tsv(sync_ao,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_sync-alphaomega.tsv']);

%COMBINING SYNC ENTRIES INTO UNIQUE TABLE
%sync = bml_annot_rowbind(sync_audio_all, sync_ripple, sync_ao);
sync = bml_annot_rowbind(sync_audio_all, sync_ripple);
bml_annot_write_tsv(sync,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_sync.tsv']);
