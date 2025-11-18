-%Loading paths

ft_defaults()
bml_defaults()
format long

%Defining paths 

SUBJECT='DM1046';
SESSION = 'intraop';
TASK = 'ERNA';
%PATH_DATASET = '/Volumes/Nexus4/DBS';
PATH_DATASET = 'Y:\DBS';

PATH_DER = [PATH_DATASET filesep 'derivatives'];
PATH_DER_SUB = [PATH_DER filesep 'sub-' SUBJECT];  
PATH_PREPROC = [PATH_DER_SUB filesep 'preproc'];
PATH_ANNOT = [PATH_DER_SUB filesep 'annot'];
PATH_FIELDTRIP = [PATH_DER_SUB filesep 'fieldtrip'];

PATH_SRC = [PATH_DATASET filesep 'sourcedata'];
PATH_SRC_SUB = [PATH_SRC filesep 'sub-' SUBJECT];  
PATH_RIPPLE = [PATH_SRC_SUB filesep 'ses-intraop' filesep 'ripple'];
PATH_ALHAOMEGA = [PATH_SRC_SUB filesep 'ses-intraop' filesep 'alphaomega']; %path to neuroomega's mpx files.
PATH_ALHAOMEGA_MAT = [PATH_ALHAOMEGA]; %path to neuroomega's converted mat files.
PATH_AUDIO_DER = [PATH_DER_SUB filesep 'aec']; %path to acoustic echo cancelling folder
PATHS_AUDIO_SRC = strcat(PATH_SRC_SUB,filesep,{'ses-training';'ses-preop';'ses-intraop'},filesep,'audio');
PATHS_TASK = strcat(PATH_SRC_SUB,filesep,{'ses-training';'ses-preop';'ses-intraop'},filesep,'task');

% loading annotation tables
sync = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_sync.tsv']);
% sync.folder = strrep(sync.folder,'Y:\DBS',PATH_DATASET);
% sync.folder = strrep(sync.folder,'\',filesep);
sessions = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_sessions.tsv']);
runs = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_runs.tsv']);
%rchannels = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_channels.tsv']);
%stable_MER = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_stable-MER.tsv']);

epoch = bml_annot_extend(runs(strcmp(runs.task,TASK) & strcmp(runs.session,SESSION),:),5,5);
epoch = epoch(~contains(epoch.comment,'test'),:) 

%Saving stimulation events as annots
cfg=[];
cfg.roi = sync(sync.filetype=="trellis.nev",:);
stim_events = bml_annot_read_stim_nev(cfg);
stim_events = bml_annot_filter(stim_events,epoch);
bml_annot_write_tsv(stim_events,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_task-' TASK '_stimulation-pulses.tsv']);

cfg=[];
cfg.criterion = @(x) (abs(max(x.ends)-min(x.starts))<0.1);
stim_events_cons = bml_annot_consolidate(cfg,stim_events);
bml_annot_write_tsv(stim_events_cons,[PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_task-' TASK  '_stimulation-bursts.tsv']);


cfg=[];
cfg.y = 'file_id'; 
figure();
bml_annot_plot(cfg,stim_events_cons)
