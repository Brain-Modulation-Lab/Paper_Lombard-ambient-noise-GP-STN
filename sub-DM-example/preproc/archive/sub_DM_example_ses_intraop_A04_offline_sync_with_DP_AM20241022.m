%% loading toolbox
bml_defaults
ft_defaults
format long
clear

%% Defining paths
SUBJECT='DM####';
SESSION='intraop';
TASKS = {'sync', 'sent_onset_stim', 'ERNA', 'SMSL'};
if ispc
    PATH_DATASET = 'Y:\DBS';
elseif ismac
    PATH_DATASET = '/Volumes/Nexus4/DBS';
end
%PARAM_AO_CONST_DEPTH_DURATION = 30; % min duration required for AlphaOmega constant depth epochs
PARAM_MIN_N_EVENTS = 10; %min number of events required to attempt sync

PATH_DER = [PATH_DATASET filesep 'derivatives'];
PATH_DER_SUB = [PATH_DER filesep 'sub-' SUBJECT];  
PATH_PREPROC = [PATH_DER_SUB filesep 'preproc'];
PATH_ANNOT = [PATH_DER_SUB filesep 'annot'];

PATH_SRC = [PATH_DATASET filesep 'sourcedata'];
PATH_SRC_SUB = [PATH_SRC filesep 'sub-' SUBJECT];   
PATH_RIPPLE = [PATH_SRC_SUB filesep 'ses-intraop' filesep 'ripple'];
%PATH_ALHAOMEGA = [PATH_SRC_SUB filesep 'ses-intraop' filesep 'alphaomega']; %path to neuroomega's mpx files.
%PATH_ALHAOMEGA_MAT = [PATH_ALHAOMEGA]; %path to neuroomega's converted mat files.
%PATH_AUDIO_DER = [PATH_DER_SUB filesep 'denoised-audio']; %path to acoustic echo cancelling folder
PATHS_AUDIO_SRC = strcat(PATH_SRC_SUB,filesep,{'ses-training';'ses-preop';'ses-intraop'},filesep,'audio');
PATHS_TASK = strcat(PATH_SRC_SUB,filesep,{'ses-training';'ses-preop';'ses-intraop'},filesep,'task');

%% loading annotation tables
events_files = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_events-files.tsv']);
events_files = events_files(contains(events_files.filename,'ses-intraop'),:);
if ispc
    events_files.path = strrep(events_files.path,'/Volumes/Nexus4','Y:');
elseif ismac
    events_files.path = strrep(events_files.path,'Y:','/Volumes/Nexus4');
end
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

% % --- COMMENT OUT IF SUBJECT DOESN'T HAVE ALPHAOMEGA -----
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

% cfg=[];
% cfg.path=PATH_AUDIO_DER;
% cfg.pattern='*.wav';
% cfg.filetype='audio.wav';
% info_wav=bml_info_raw(cfg);
% info_wav.basename = strrep(info_wav.name,'.wav','');
% info_wav = info_wav(~ismember(info_wav.name,info_wav_src_renamed.name),:);
% info_wav = bml_annot_rowbind(info_wav_src_renamed,info_wav);
info_wav = info_wav_src_renamed;


%Saving files information table with OS times
info_orig = bml_annot_rowbind(info_ns5,info_nf3,info_nf6,info_nev, ...
                                info_wav,info_wav_src_renamed);
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
    %loading events from nev 
    cfg=[];
    cfg.roi = bml_roi_table(info_nev(i,:));
    nev = bml_annot_read_event_nev(cfg);

    %matching events by similarity
    cfg=[];
    cfg.timetol=0.001;
    cfg.diagnostic_plot = true;
    [idxs_task_events, idxs_nev] = bml_sync_match_events(cfg, task_events, nev);   

    slave_delta_t = median(task_events.starts(idxs_task_events) - nev.starts(idxs_nev),'omitnan');
    %%Aligning ripple to presenation task to get GTC. No timewarping.
    %cfg=[];
    %cfg.scan = 100;
    %[slave_delta_t, min_cost, warpfactor] = bml_timealign_annot(cfg, task_events(idxs_task_events,:), nev);
    
    info_nev.delta_t(i)=slave_delta_t;
    info_nev.n(i)=height(nev);
end
info_nev.delta_t(ismissing(info_nev.delta_t)) = 0;

%adjusting ripple files
cfg=[];
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
sync_ripple = bml_annot_rowbind(roi_nev, roi_ns5, roi_nf3, roi_nf6);
bml_annot_write_tsv(sync_ripple, [PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION  '_sync-ripple.tsv'])

% %% DETECT  STRETCHES CONSTANT DEPTH IN ALPHA OMEGA
% cfg=[];
% cfg.criterion = @(x) length(unique(x.depth))==1 && abs((max(x.ends)-min(x.starts))-sum(x.duration))<.001;
% neuro_cons_depth = bml_annot_consolidate(cfg,info_neuroomega);
% ao_session = neuro_cons_depth(neuro_cons_depth.duration > PARAM_AO_CONST_DEPTH_DURATION,:);
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
events_files.warpfactor(:)=nan;
events_files.n(:)=nan;
for i=1:height(events_files)

    %verifyng that there are matching nev events for file
    %(prevents attempting to match preop events to nev files)
    nev_events_files = bml_annot_filter(roi_nev,events_files(i,:));
    if isempty(nev_events_files); continue; end
   
    task_events = bml_annot_read_tsv(fullfile(events_files.path{i},events_files.filename{i}));
  
    %matching task_events to master_events
    cfg=[];
    cfg.timetol=0.003;
    cfg.diagnostic_plot = true;
    cfg.weight_onset = 0.1;
    [idxs_master_events,idxs_task_events, mean_sim, simv]=bml_sync_match_events2(cfg,master_events,task_events);    

    disconts_idx = [1, find(diff(idxs_master_events) > 3), length(idxs_master_events)];
    [~, cont_map_idx_idx] = max(diff(disconts_idx));
    start_idx = disconts_idx(cont_map_idx_idx) + 1;
    end_idx = disconts_idx(cont_map_idx_idx+1) - 1;
    
    %linear warping as fallback strategy
    tbar = mean(task_events.starts(idxs_task_events),'omitnan');
    p = polyfit(task_events.starts(idxs_task_events(start_idx:end_idx)) - tbar, master_events.starts(idxs_master_events(start_idx:end_idx)) - tbar,1);
    %figure; plot(task_events.starts(idxs_task_events(start_idx:end_idx)) - tbar, master_events.starts(idxs_master_events(start_idx:end_idx)) - tbar)
    
    events_files.warpfactor(i) = p(1);
    events_files.delta_t(i) = p(2);
    events_files.n(i)=height(task_events(idxs_task_events,:));;

    %applying timealignment to events as fallback strategy
    task_events.starts_corrected = (task_events.starts - tbar) .* p(1) + tbar +  p(2);
        
    %getting time from master_events according to alignemnt (when available)
    task_events.starts = bml_map(1:height(task_events),idxs_task_events,master_events.starts(idxs_master_events),nan);

    %using fallback times when nev times not available
    task_events.starts(isnan(task_events.starts))=task_events.starts_corrected(isnan(task_events.starts));

    bml_annot_write_tsv(task_events,[PATH_ANNOT filesep strrep(events_files.filename{i},'.tsv','-sync.tsv')]);
    fprintf('sync %s, delta_t=%f, mean_sim=%f, warpfactor=%f \n',events_files.filename{i}, p(2), mean_sim,  p(1))

end

%% SYNCHRONIZE AUDIO WITH EVENTS USING DIGITAL TRIGGERS

% converting ripple events to output sync trigger values
master_trig = master_events;
master_trig.tmp = [2^15;abs(diff(master_trig.value))];
master_trig = master_trig(master_trig.tmp >= 2^9,:);
master_trig.value = (master_trig.value > 2^9) .* 1.0;

% doing event based alignment
info_wav_events = info_wav_src_renamed(contains(info_wav_src_renamed.name,'_events.wav') & contains(info_wav_src_renamed.name,TASKS),:);
%info_wav_events = info_wav_src_renamed(contains(info_wav_src_renamed.name,'_events.wav'),:);
sync_wav = table();
%hF = gobjects(0,1);

for i=1:height(info_wav_events)
    cfg1=[];
    cfg1.roi = info_wav_events(i,:);
    cfg1.flip_polarity = 0;
    wav_events = bml_annot_read_event_wav(cfg1);

    %matching task_events to master_events
    cfg=[];
    cfg.timetol=0.003;
    cfg.diagnostic_plot = true;
    [idxs_trig_events,idxs_wav_events, mean_sim, sim]=bml_sync_match_events(cfg,master_trig,wav_events);   

    disconts_idx = [1, find(diff(idxs_trig_events) > 3), length(idxs_trig_events)];
    [~, cont_map_idx_idx] = max(diff(disconts_idx));
    start_idx = disconts_idx(cont_map_idx_idx) + 1;
    end_idx = disconts_idx(cont_map_idx_idx+1) - 1;

    idxs_trig_events = idxs_trig_events(start_idx:end_idx);
    idxs_wav_events = idxs_wav_events(start_idx:end_idx);
    %idxs_trig_events = idxs_trig_events(sim > 0.9);
    %idxs_wav_events = idxs_wav_events(sim > 0.9);
    
    %linear warping
    %tbar = mean(wav_events.starts(idxs_wav_events),'omitnan');
    %p = polyfit(wav_events.starts(idxs_wav_events) - tbar, master_trig.starts(idxs_trig_events) - tbar,1);

%     cfg=[];
%     cfg.timetol=0.001;
%     cfg.diagnostic_plot = 1;
%     cfg.plot_title = strrep(info_wav_events.name(i),'_','\_');
%     [i_sync_wav, i_hF] = bml_sync_event(cfg,master_trig,wav_events);

    pst = polyfit(wav_events.sample(idxs_wav_events), master_trig.starts(idxs_trig_events),1);
    ptt = polyfit(wav_events.starts(idxs_wav_events),master_trig.starts(idxs_trig_events),1);
    i_sync_wav = table();
    i_sync_wav.starts = master_trig.starts(idxs_trig_events(1));
    i_sync_wav.ends = master_trig.starts(idxs_trig_events(end));    
    i_sync_wav.name = height(wav_events(idxs_wav_events,:));;
    i_sync_wav.s1 = wav_events.sample(idxs_wav_events(1));
    i_sync_wav.s2 = wav_events.sample(idxs_wav_events(end));
    i_sync_wav.t1 = i_sync_wav.s1 * pst(1) + pst(2);
    i_sync_wav.t2 = i_sync_wav.s2 * pst(1) + pst(2);
    i_sync_wav.warpfactor = ptt(1);
    i_sync_wav.delta_t = ptt(2);
    i_sync_wav.name = info_wav_events.name(i);
    i_sync_wav.folder = info_wav_events.folder(i);
    i_sync_wav.chantype = info_wav_events.chantype(i);
    i_sync_wav.filetype = info_wav_events.filetype(i);
    i_sync_wav.Fs = info_wav_events.Fs(i);
    i_sync_wav.nSamples = info_wav_events.nSamples(i);

    sync_wav = bml_annot_rowbind(sync_wav, bml_roi_confluence(i_sync_wav));
%    hF = [hF; i_hF];
end

% savefig(hF, ['./figures' filesep 'sub-' SUBJECT '_ses-intraop_proto-A04-fine-sync-audio_' ...
%     char(string(datestr(now, 'YYYYmmDDMM'))) '.fig']);

% Transfering synchronization to other wav files
sync_wav = sync_wav(sync_wav.duration > 1,:);
sync_audio_cons = sync_wav(contains(sync_wav.name,SESSION) & contains(sync_wav.name,TASKS),:);
cinfo_wav = info_wav(contains(info_wav.name,SESSION) & contains(info_wav.name,TASKS),:);

% matching files by basename
tmp = regexp(sync_audio_cons.name,'^(sub-[\w]*_ses-[\w]*_task-[\w-]*_run-[0-9]*)_[\w-\.]*$','tokens','once','forceCellOutput');
sync_audio_cons.basename  = [tmp{:}]';
tmp = regexp(cinfo_wav.name,'^(sub-[\w]*_ses-[\w]*_task-[\w-]*_run-[0-9]*)_[\w-\.]*$','tokens','once','forceCellOutput');
cinfo_wav.basename  = [tmp{:}]';
tmp=regexp(cinfo_wav.name,'^sub-[\w]*_ses-[\w]*_task-[\w-]*_run-[0-9]*_(recording-)?([a-zA-Z0-9-]*)(_[\w]*)?\.wav$','tokens','once');
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

% %% SYNCHRONIZE NEUROOMEGA FILE WITH EVENTS USING DIGITAL TRIGGERS
% 
% % synchronizing by matching events
% info_neuroomega_events = info_neuroomega;
% sync_ao = table();
% hF = gobjects(0,1);
% for i=1:height(info_neuroomega_events)
% 
%     i_slave_events = bml_read_event(info_neuroomega_events(i,:));
%     if isempty(i_slave_events) || length(i_slave_events) < PARAM_MIN_N_EVENTS
%         continue
%     end
% 
%     cfg1=[]; cfg1.Fs = 1; %already in seconds
%     cfg1.starts = info_neuroomega_events.starts(i);
%     i_slave_events = bml_event2annot(cfg1,i_slave_events);
%     i_slave_events.value = strcmp(i_slave_events.type,'CDIG_IN_1_Up');
%     i_slave_events.sample = round(i_slave_events.sample .* 44000);
% 
%     cfg=[];
%     cfg.timetol=0.001;
%     cfg.diagnostic_plot = 1;
%     cfg.plot_title = strrep(info_neuroomega_events.name(i),'_','\_');
%     [i_sync_ao, i_hF] = bml_sync_event(cfg,master_trig,i_slave_events);
% 
%     i_sync_ao.name(:) = info_neuroomega_events.name(i);
%     i_sync_ao.folder(:) = info_neuroomega_events.folder(i);
%     i_sync_ao.chantype(:) = {'events'};
%     i_sync_ao.filetype(:) = info_neuroomega_events.filetype(i);
%     i_sync_ao.Fs(:) = 44000;
%     i_sync_ao.nSamples(:) = info_neuroomega_events.nSamples(i) .* 44000 ./ info_neuroomega_events.Fs(i);
% 
%     sync_ao = bml_annot_rowbind(sync_ao, i_sync_ao);
%     hF = [hF; i_hF];
% 
% end
% 
% savefig(hF, ['./figures' filesep 'sub-' SUBJECT '_ses-intraop_proto-A04-fine-sync-alphaomega_' ...
%     char(string(datestr(now, 'YYYYmmDDMM'))) '.fig']);
% bml_annot_write_tsv(sync_ao,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_sync-alphaomega.tsv']);

%% COMBINING SYNC ENTRIES INTO UNIQUE TABLE
%sync = bml_annot_rowbind(sync_audio_all, sync_ripple, sync_ao);
sync = bml_annot_rowbind(sync_audio_all, sync_ripple);
bml_annot_write_tsv(sync,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_sync.tsv']);


