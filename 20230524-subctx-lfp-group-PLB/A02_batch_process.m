close all; clear all; 

PROTOCOL_PATH = 'Y:\DBS\groupanalyses\task-lombard\20230524-subctx-lfp-group-PLB';
PROTOCOL_FUNCTION = 'A01_prepare_subject_TFR';
% PROTOCOL_TABLE = 'Subjects.tsv';
exe_daytime = datestr(now,'yyyymmdd_HHMM');
LOG_PATH = [PROTOCOL_PATH filesep 'A02_batch-' PROTOCOL_FUNCTION '_' exe_daytime]; 
diary([LOG_PATH '.log']); 

SKIP_OK = false; %Should previous protocols run by this script successfully be skipped
FORCE = true; %Archive all previous versions of the script and run current
%overrides any manual modification
%SUBJECTS = {'DBS3001'};

PATH_DATA = 'Y:\DBS';
cd(PROTOCOL_PATH)

% subject_table = readtable(PROTOCOL_TABLE,'Delimiter','\t','FileType','Text');
% participants overview
subject_table = readtable("Y:\DBS\participants.tsv",'FileType', 'delimitedtext', 'Delimiter', '\t');
SESSION = 'intraop'; 
TASK = 'lombard'; 


fprintf('=== Running protocol %s ===\n',PROTOCOL_FUNCTION)
if FORCE; fprintf('Forced run, overwritting any manual change.\n'); end
if SKIP_OK; fprintf('Skipping previously successfully executed protocols.\n'); end

faillog = table();


%%  2025 04 25 run subject x channel specific  
electrode = bml_annot_read_tsv('../20230328-subctx-ctx-group-coverage-PLB/A02b_electrodes-all-subj-DM1001-DM1039-with-lombard-run-id.tsv'); 
% electrode = bml_annot_read_tsv('A02b_electrodes-all-subj-DM1001-DM1039.tsv'); 
% elecs = bml_annot_read_tsv('..\20230524-subctx-lfp-group-PLB/A04b_channels_summary_with_electrode_info.tsv'); % get this table from anotoher analysis

% load atlas with groupings of HCP
labels2areas = readtable('Z:/Resources/HCPMMP1-Labeling-Atlas/HCPMMP1-labels2areas.tsv','FileType','text','Delimiter', '\t'); 
labels2areas.color = hex2rgb(labels2areas.color); 

cfg=[];cfg.keys='label'; 
electrode.label = electrode.HCPMMP1_label_1; 
electrode = bml_annot_left_join(cfg, electrode, labels2areas); 

ELEC_TYPE = 'ECOG';  % SEEG for macro, MICRO for micro

elec = electrode; 
elec = elec(elec.type==string(ELEC_TYPE), :); 
% elec = elec(randi([1 height(elec)], 1000, 1), :);  % for testing only 

tblstats2 = grpstats(elec,"area"); 


% % op1 posterior auditory cortex -----
% select channels that are in posterior auditory cortex, inspected manually
% for y==-14 estimate
% CHANNELS_SUBSET_NAME = 'pos-AC';
% elec = elec(elec.mni_y<-14, :); 
% elec = elec(ismember(elec.area, {'AAC', 'TPOJ'}), :); 

% % op2 pre and post CG ---------------------
% CHANNELS_SUBSET_NAME = 'pre-CG-post-CG';
% elec = elec(ismember(elec.area, {'PMC', 'SMC', 'POC'}), :); 

% % op3 DLPFC-IFG-POC ---------------------
% CHANNELS_SUBSET_NAME = 'DLPFC-IFG-POC';
% elec = elec(ismember(elec.area, {'DLPFC', 'IFG', 'POC'}), :); 

% % op4 PMC ---------------------
% CHANNELS_SUBSET_NAME = 'PMC';
% elec = elec(ismember(elec.area, {'PMC'}), :); 
% 

% % op4 IPC + TPOJ + AAC + LTC ---------------------
% CHANNELS_SUBSET_NAME = 'IPC-TPOJ-AAC-LTC';
% elec = elec(ismember(elec.area, {'IPC', 'TPOJ', 'AAC', 'LTC'}), :); 

% % op4 SMC ---------------------
% CHANNELS_SUBSET_NAME = 'SMC';
% elec = elec(ismember(elec.area, {'SMC'}), :); 



% op4 IPC + TPOJ + AAC + LTC ---------------------
CHANNELS_SUBSET_NAME = 'ecog-1001-1039';
elec = elec; 


% FOR ECOG analyses... nest channel names for each subject
T = elec; 
% Assuming your table is called T
% Convert from cell arrays of strings to character vectors if needed
T.subject_id = string(T.subject_id);
T.name = string(T.name);
% Group by subject_id
[G, subject_list] = findgroups(T.subject_id);
% Collect 'name' values into nested cell arrays of character vectors
nested_names = splitapply(@(names) ({cellstr(names)}), T.name, G);
nested_names = cellfun(@unique, nested_names, 'UniformOutput',false); 
% Create the nested output table
subject_table = table(subject_list, nested_names, ...
                'VariableNames', {'subject_id', 'channels'});







% % op3 ------ all VIM lombard dbs channels
% soi = {'DM1007','DM1025', 'DM1026', 'DM1032', 'DM1035', 'DM1038', 'DM1039'}; 
% subject_table = subject_table(ismember(subject_table.subject_id, soi), :); 
% CHANNELS = {{'dbs_*', 'audio*', 'envaudio*'}};
% CHANNELS_SUBSET_NAME = 'VIM-DBS-bipolar-TFR'; 
% subject_table.channels(:) = CHANNELS; % inspect only DBS channels
% CFG_TFR = 'tfr-withinscript';

% % op4 ------ theta-alpha extraction for all channels
% CHANNELS = {'macro_*'};
% CHANNELS_SUBSET_NAME = 'macro-thetaalpha'; 
% subject_table.channels(:) = CHANNELS; % inspect only DBS channels
% 
% CFG_TFR = [];
% CFG_TFR.id = '010_thetaalpha_mtmconvol_band-4-12'; 
% CFG_TFR.method = 'ft-mtmconvol';  % ft-dpss, bandpass-hilbert, bandpass-bank-hilbert
% CFG_TFR.fbank_collapse_method = 'none'; % none, pca, median, mean, zscore
% CFG_TFR.fbank_withinband_norm = 'none'; % none, mean, median, zscore
% CFG_TFR.fbank_n = 1; % pca, median, mean
% CFG_TFR.freq_lowhigh = [4 12];
% CFG_TFR.fs_out = 100;


% % % op5 ------ beta extraction for all channels
% CHANNELS = {'macro_*'};
% CHANNELS_SUBSET_NAME = 'macro-beta'; 
% subject_table.channels(:) = CHANNELS; % inspect only DBS channels
% 
% CFG_TFR = [];
% CFG_TFR.id = '020_beta_mtmconvol_band-12-30'; 
% CFG_TFR.method = 'ft-mtmconvol';  % ft-dpss, bandpass-hilbert, bandpass-bank-hilbert
% CFG_TFR.fbank_collapse_method = 'none'; % none, pca, median, mean, zscore
% CFG_TFR.fbank_withinband_norm = 'none'; % none, mean, median, zscore
% CFG_TFR.fbank_n = 1; % pca, median, mean
% CFG_TFR.freq_lowhigh = [12 30];
% CFG_TFR.fs_out = 100;

%% Loop over subjects
for i=1:height(subject_table)
    
    SUBJECT = subject_table.subject_id{i};
    CHANNELS = subject_table.channels{i}; 
    subject_num = str2num(SUBJECT(end-2:end));
    % if subject_num < 20; continue; end
%     if ~ismember(subject_num, [3, 12, 33]); continue; end

    timetol = 0.01;
    fprintf('Running protocol: %s.', SUBJECT)

    %running protocol
    try
        proto = str2func(PROTOCOL_FUNCTION);
        proto(SUBJECT, CHANNELS, CHANNELS_SUBSET_NAME, CFG_TFR);
        fprintf('OK\n');

    catch err
        fprintf('FAILED: %s\n',err.message)
        nr = [];
        nr.subject = string(SUBJECT);
        nr.session = string(SESSION);
        nr.task = string(TASK);
        nr.identifier = string(err.identifier);
        nr.error = string(err.message);
        nr.timestamp = datetime('now');
        if isempty(faillog)
            faillog = struct2table(nr);
        else
            faillog = [faillog; struct2table(nr)];
        end
    end
end

faillog_fname  = [LOG_PATH '.tsv'];  
% faillog_orig = bml_annot_read_tsv(faillog_fname); 
% faillog = [faillog_orig; faillog]; 
writetable(faillog, faillog_fname, 'FileType','text','Delimiter','\t'); 
