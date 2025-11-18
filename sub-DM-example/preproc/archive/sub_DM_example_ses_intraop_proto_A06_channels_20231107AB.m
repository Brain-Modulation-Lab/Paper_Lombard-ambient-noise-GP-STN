%% Loading paths
ft_defaults
bml_defaults
format long
 
%% Defining paths
SUBJECT='DM1026';
SESSION = 'intraop';
PATH_DATASET = '/Volumes/Nexus4/DBS';
 
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
sync = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_ses-intraop_sync.tsv'])
sessions = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_sessions.tsv'])
runs = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_runs.tsv'])
%stable_MER = bml_annot_read_tsv([PATH_ANNOT filesep 'sub-' SUBJECT '_ses-intraop_stable-MER.tsv'])
 
sync.folder = strrep(sync.folder,'Y:\DBS',PATH_DATASET);
sync.folder = strrep(sync.folder,'\','/');
 
%% Loading header information for nf6 files
cfg=[];
cfg.roi=sync(sync.filetype=="trellis.nf6",:);
cfg.roi=cfg.roi(1,:);
hdr=bml_read_header(cfg);
 
nrow = numel(hdr.label);
 
ch_nf6 = table();
ch_nf6.name = repmat({''},[nrow,1]);
ch_nf6.type = repmat({'ECOG'},[nrow,1]);
ch_nf6.units = deblank({hdr.orig.hdr.analoginfo.Units}');
ch_nf6.filetype = repmat({'trellis.nf6'},[nrow,1]);
ch_nf6.channel = hdr.label;
ch_nf6.connector = repmat([0],[nrow,1]);
ch_nf6.port = repmat({''},[nrow,1]);
ch_nf6.description = repmat({''},[nrow,1]);
ch_nf6.sampling_freq = [hdr.orig.hdr.analoginfo.SampleRate]';
ch_nf6.low_cutoff = [hdr.orig.hdr.analoginfo.LowFreqCorner]';
ch_nf6.low_cutoff(ch_nf6.low_cutoff > ch_nf6.sampling_freq) = 0;
ch_nf6.high_cutoff = [hdr.orig.hdr.analoginfo.HighFreqCorner]';
ch_nf6.notch = repmat([0],[nrow,1]);
ch_nf6.status = repmat({''},[nrow,1]);
ch_nf6.status_description = repmat({''},[nrow,1]);
 
 
%% Loading header information for ns5 files
cfg=[];
cfg.roi=sync(sync.filetype=="trellis.ns5",:);
cfg.roi=cfg.roi(1,:);
cfg.chantype = 'analog';
hdr=bml_read_header(cfg);
 
ch_ns5 = table();
ch_ns5.channel = hdr.label;
ch_ns5.units = deblank({hdr.orig.ElectrodesInfo.AnalogUnits}');
ch_ns5.sampling_freq = repmat(hdr.Fs,[height(ch_ns5),1]);
ch_ns5.low_cutoff = [hdr.orig.ElectrodesInfo.LowFreqCorner]';
ch_ns5.high_cutoff = [hdr.orig.ElectrodesInfo.HighFreqCorner]';
ch_ns5.filetype = repmat({'trellis.ns5'},[height(ch_ns5),1]);
 
% %% Loading header information for NeuroOmega files
% cfg=[];
% cfg.roi=sync(sync.filetype=="neuroomega.mat",:);
% cfg.roi.chantype=[];
% [~,fidx] = max(cfg.roi.duration);
% cfg.roi=cfg.roi(fidx ,:);
% hdr=bml_read_header(cfg);
%  
% ch_no = table();
% ch_no.channel = hdr.orig.chaninfo.channel;
% ch_no.sampling_freq = hdr.orig.chaninfo.Fs;
% ch_no.chantype = hdr.orig.chaninfo.chantype;
% ch_no.filetype = repmat({'neuroomega.mat'},[height(ch_no),1]);
%  
 
%% Loading header information audio files
 
cfg=[];
cfg.criterion = @(x) true;
cfg.groupby = 'chantype';
tmp = bml_annot_consolidate(cfg,sync(sync.filetype=="audio.wav",:));
 
ch_audio = table();
ch_audio.channel = tmp.chantype;
ch_audio.sampling_freq = tmp.Fs;
ch_audio.chantype = tmp.chantype;
ch_audio.type = repmat({'AUDIO'},[height(ch_audio),1]);
ch_audio.filetype = repmat({'audio.wav'},[height(ch_audio),1]);
ch_audio.units = repmat({'unknown'},[height(ch_audio),1]);
 
 
%% joining all channels
 
%ch=bml_annot_rowbind(ch_nf6,ch_ns5,ch_no,ch_audio);
ch=bml_annot_rowbind(ch_nf6,ch_ns5,ch_audio);

ch_fname = [PATH_ANNOT filesep 'sub-' SUBJECT '_ses-' SESSION '_channels-orig.tsv'];
writetable(ch,ch_fname,'delimiter','\t','FileType','text');
 
%manually complete names and active interval for each electrode
%rename file removing the '-orig' string and channels-orig archive  


% checking channel content manually
cfg=[]; %ns5
cfg.roi=sync(strcmp(sync.name,'sub-DM1026_ses-intraop_datafile0003.ns5'),:); 
%cfg.channel='analog'; %stimulus audio
ns5_analog = bml_load_continuous(cfg);
bml_praat(ns5_analog);
 
% cfg=[]; %no
% cfg.roi=sync(strcmp(sync.name,'LT1D-0.202F0002.mat'),:); 
% cfg.channel='CANALOG_IN_1'; %stimulus audio
% no_analog = bml_load_continuous(cfg);
% bml_praat(no_analog);
