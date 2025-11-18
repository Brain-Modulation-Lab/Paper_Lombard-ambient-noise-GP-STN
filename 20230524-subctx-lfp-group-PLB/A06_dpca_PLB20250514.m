%% Basic setup 
clc;
clear all;
close all;
tStart = tic;

%initialize
bml_defaults;
ft_defaults;
prettify_plots;


%% Config
if ispc
    PATH_SERVER = 'Y:';
elseif ismac
    PATH_SERVER = '/Volumes/Nexus4';
end

MODALITY = 'LFP';           % 'LFP' or 'SU'
CHANNEL_SUBSET = 'DLPFC-IFG-POC';      % e.g., 'BG', 'PMC', 'IPC-TPOJ-AAC-LTC', 
CTX_FLAG = true;           % Is this cortex-specific analysis?
loudness_lvls = 'all'; % all, lowhalf-loudhalf
TASK = 'lombard';
level = 'sentence';

% Create log_hash
log_hash = ['202505141139' '-' MODALITY '-' CHANNEL_SUBSET '_loudness-' loudness_lvls]




% MODALITY = 'LFP'; % LFP, SU
% CHANNEL_SUBSET = 'BG'; % BG, IPC-TPOJ-AAC-LTC, 
% CTX_FLAG = false; 
% log_hash = ['202505141139' '-' MODALITY '-' CHANNEL_SUBSET];
% 
% % RMDIM3 = true; % this will replace the loudness level with dummy data   
% % log_hash = ['202505141139' '-' CHANNEL_SUBSET '_RMDIM3-' num2str(RMDIM3)];
% 
% 
% 
% % % LFP--basal ganglia
% loudness_lvls = 'lowhalf-loudhalf'; % 'lowhalf-loudhalf'
% log_hash = [log_hash '_loudness-' loudness_lvls]
% PATH_DATA_INPUT = ['Y:\DBS\groupanalyses\task-lombard\20230524-subctx-lfp-group-PLB\data\' ... 
%                    'tlock-audio-on_chs-macro-ecogL-ecogL-ecogL-audio-envaudio\dpca_tensor_timewarped_loudness-lowhalf-loudhalf_202505141221.h5']; 

% % LFP--cortex PMC
% PATH_DATA_INPUT = "Y:/DBS/groupanalyses/task-lombard/20230524-subctx-lfp-group-PLB/data/tlock-audio-on_chs-PMC/dpca_tensor_timewarped_withmeta_202505071709.h5";

% % SU
% loudness_lvls = 'lowhalf-loudhalf'; % 'lowhalf-loudhalf'
% PATH_DATA_INPUT = 'Y:\DBS/groupanalyses/task-lombard/20231122-firing-rate_rep/latane/data/dpca_tensor_timewarped_loudness-lowhalf-loudhalf_202505141138.h5';
% log_hash = [log_hash '_loudness-' loudness_lvls]
%% Set paths v1  
% % which task?
% TASK = 'lombard';
% level = 'sentence';
% 
% 
% PATH_DBS = fullfile(PATH_SERVER,'DBS');
% 
% if MODALITY=="SU"
%     PATH_DERIVATIVES = fullfile(PATH_DBS,'derivatives');
%     PATH_GROUPANALYSES = fullfile(PATH_DBS,'groupanalyses');
%     PATH_ANALYSIS = fullfile(PATH_GROUPANALYSES,'task-lombard','20231122-firing-rate_rep', 'latane');
%     PATH_ANNOTOUTPUT = fullfile(PATH_ANALYSIS,'annot');
%     PATH_RESULTSOUTPUT = fullfile(PATH_ANALYSIS,'results','output');
%     PATH_RESULTSLOG = fullfile(PATH_ANALYSIS,'results','log');
%     PATH_FIGUREOUTPUT = fullfile(PATH_ANALYSIS,'results','figures');
%     PATH_AUDIOSTIMULUS_ANNOT = fullfile(PATH_GROUPANALYSES,'task-lombard','20210818-phonetic-coding-stimulus-audio','annot');
%     PATH_RESOURCES = fullfile(PATH_ANALYSIS,'resources');
%     PATH_FIGUREOUTPUT_FIRINGRATE = fullfile(PATH_FIGUREOUTPUT,'rate',['loghash_',log_hash],level);
%     PATH_RESULTSOUTPUT_FIRINGRATE = fullfile(PATH_RESULTSOUTPUT,'rate',['loghash_',log_hash],level);
% 
%     PATH_RESULTSDPCAOUTPUT = fullfile(PATH_ANALYSIS,'decoding-results','output','rate',['loghash_',log_hash],level,'dPCA');
% 
% elseif MODALITY=="LFP"
%     PATH_DERIVATIVES = fullfile(PATH_DBS,'derivatives');
%     PATH_GROUPANALYSES = fullfile(PATH_DBS,'groupanalyses');
%     PATH_ANALYSIS = fullfile(PATH_GROUPANALYSES,'task-lombard','20230524-subctx-lfp-group-PLB');
%     PATH_ANNOTOUTPUT = fullfile(PATH_ANALYSIS,'annot');
%     PATH_RESULTSOUTPUT = fullfile(PATH_ANALYSIS,'results','output');
%     PATH_RESULTSLOG = fullfile(PATH_ANALYSIS,'results','log');
%     PATH_FIGUREOUTPUT = fullfile(PATH_ANALYSIS, 'channel_summary', 'fig', 'A06_dpca');
%     PATH_AUDIOSTIMULUS_ANNOT = fullfile(PATH_GROUPANALYSES,'task-lombard','20210818-phonetic-coding-stimulus-audio','annot');
%     PATH_RESOURCES = fullfile(PATH_ANALYSIS,'resources');
% end


%% Set paths v2

% Set Base Paths
PATH_SERVER = 'Y:'; % Or '/Volumes/Nexus4' if on Mac
PATH_DBS = fullfile(PATH_SERVER, 'DBS');
PATH_GROUPANALYSES = fullfile(PATH_DBS, 'groupanalyses');
PATH_DERIVATIVES = fullfile(PATH_DBS, 'derivatives');

% Modality-Specific Path
switch MODALITY
    case 'SU'
        PATH_ANALYSIS = fullfile(PATH_GROUPANALYSES, 'task-lombard', '20231122-firing-rate_rep', 'latane');
        PATH_DATA_INPUT = fullfile(PATH_ANALYSIS, 'data', 'dpca_tensor_timewarped_loudness-lowhalf-loudhalf_202505141138.h5');
         ['dpca_tensor_timewarped_loudness-' loudness_lvls '_202505141346.h5']
        PATH_FIGUREOUTPUT = fullfile(PATH_ANALYSIS, 'results', 'figures');
        PATH_RESULTSOUTPUT = fullfile(PATH_ANALYSIS, 'results', 'output');
        PATH_RESULTSLOG = fullfile(PATH_ANALYSIS, 'results', 'log');
        PATH_ANNOTOUTPUT = fullfile(PATH_ANALYSIS, 'annot');
        PATH_AUDIOSTIMULUS_ANNOT = fullfile(PATH_GROUPANALYSES, 'task-lombard', '20210818-phonetic-coding-stimulus-audio', 'annot');
        PATH_RESOURCES = fullfile(PATH_ANALYSIS, 'resources');

        PATH_FIGUREOUTPUT_FIRINGRATE = fullfile(PATH_FIGUREOUTPUT, 'rate', ['loghash_' log_hash], level);
        PATH_RESULTSOUTPUT_FIRINGRATE = fullfile(PATH_RESULTSOUTPUT, 'rate', ['loghash_' log_hash], level);
        PATH_RESULTSDPCAOUTPUT = fullfile(PATH_ANALYSIS, 'decoding-results', 'output', 'rate', ['loghash_' log_hash], level, 'dPCA');

    case 'LFP'
        PATH_ANALYSIS = fullfile(PATH_GROUPANALYSES, 'task-lombard', '20230524-subctx-lfp-group-PLB');
        PATH_ANNOTOUTPUT = fullfile(PATH_ANALYSIS, 'annot');
        PATH_RESULTSOUTPUT = fullfile(PATH_ANALYSIS, 'results', 'output');
        PATH_RESULTSLOG = fullfile(PATH_ANALYSIS, 'results', 'log');
        PATH_AUDIOSTIMULUS_ANNOT = fullfile(PATH_GROUPANALYSES, 'task-lombard', '20210818-phonetic-coding-stimulus-audio', 'annot');
        PATH_RESOURCES = fullfile(PATH_ANALYSIS, 'resources');

        % Handle channel subsets and CTX flag
        if CTX_FLAG
            % Cortex-specific path
            PATH_FIGUREOUTPUT = fullfile(PATH_ANALYSIS, 'fig', ['channel_summary_ctx_' CHANNEL_SUBSET]);
            PATH_DATA_INPUT = fullfile(PATH_ANALYSIS, 'data', ['tlock-audio-on_chs-' CHANNEL_SUBSET], ... 
                                        ['dpca_tensor_timewarped_loudness-' loudness_lvls '_202505191757.h5']); 
            % _202505141346.h5 for IPC-TPOJ-AAC
            % '_202505142333.h5' for PMC
             % _loudness-lowhalf-loudhalf_202505191757.h5 for DLPFC
             % _loudness-all_202505191634.h5 for DLPFC

        else
            % Subcortical (BG) or named group
            if strcmp(CHANNEL_SUBSET, 'BG')
                PATH_FIGUREOUTPUT = fullfile(PATH_ANALYSIS, 'fig', 'channel_summary');
                PATH_DATA_INPUT = fullfile(PATH_ANALYSIS, 'data', 'tlock-audio-on_chs-macro-ecogL-ecogL-ecogL-audio-envaudio', ...
                    ['dpca_tensor_timewarped_loudness-' loudness_lvls '_202505141310.h5']);

                
            else
                PATH_FIGUREOUTPUT = fullfile(PATH_ANALYSIS, 'fig', ['channel_summary_' CHANNEL_SUBSET]);
                PATH_DATA_INPUT = fullfile(PATH_ANALYSIS, 'data', ['tlock-audio-on_chs-' CHANNEL_SUBSET]);
            end
         
        end
        
end
assert(isfile(PATH_DATA_INPUT), 'PATH_DATA_INPUT not found')

% add toolbox
PATH_TOOLBOX = fullfile(PATH_SERVER,'Users','MV1019','Tools_analysis');
addpath(genpath(fullfile(PATH_TOOLBOX,'utility')))
addpath(genpath(fullfile(PATH_TOOLBOX,'spks_unit')))
% addpath(fullfile(PATH_TOOLBOX,'fieldtrip-20231220'))
addpath(genpath(fullfile(PATH_TOOLBOX,'bml')))
addpath(genpath(fullfile(PATH_TOOLBOX,'spc')))
addpath(genpath(fullfile(PATH_TOOLBOX,'dPCA-master')))


% % load raster data
% fprintf('Loading raster data [log_hash = %s] \n',log_hash)
% load(fullfile(PATH_RESULTSOUTPUT,['rasterData_',log_hash '.mat']),'rasterData')
% nUnits =  numel(rasterData.sentence);



CLRS = bml_annot_read_tsv(fullfile(PATH_GROUPANALYSES, 'task-lombard', '20230602-metadata-and-groupdata-PLB', 'lombard_colors.tsv'));
tmp = []; 
for ic = 1:height(CLRS)
tmp.(CLRS.name{ic}) = hex2rgb(CLRS.hex{ic}); 
end
CLRS = tmp; 

%% Setup: parameters
MIN_TRIALS = 5;

% combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
% decodingClasses = {[1 1; 2 2], [1 2; 1 2], [], [1 2; 3 4]};
% %decodingClasses = {[1; 2], [1 ;2], [], [1 ;2]};
% margNames = {'Lombard',    'Auditory|Speech',  'Condition-independent', 'Lombard/Auditory|Speech Interaction'};
% margColours = [CLRS.NOISE; CLRS.control;   CLRS.control; CLRS.loud3rd];


% 1 NOISE (1 2)
% 2 Loudness (1 2)
% 3 Time (1)
% 4 NOISE x Loudness
combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
decodingClasses = {[1 1; 2 2], [1 2; 1 2], [], [1 2; 3 4]};
%decodingClasses = {[1; 2], [1 ;2], [], [1 ;2]};
margNames = {'NOISE',    'Loudness',  'Independent', 'NOISE x Loudness'};
margColours = [CLRS.NOISE; CLRS.loud3rd;   CLRS.control; [121, 54, 209]/255];


% margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256; Colors = [[253 174 107]/255; [49 130 189]/255];
regularize = true;
ifSimultaneousRecording = false;
maxTrialNum = 80; %we need that a session of lombard is 80 trials
NUM_COMPS = 10;
NUM_SUBSAMPLES = 100;
NUM_SHUFFLES = 100;


% cfg1 = [];
% cfg1.lock = level;
% [~,Evts,Evts_table] = define_events(cfg1);
% nEvts = numel(Evts);
% 
% 
% trlTimes_snt_audio = [];
% trlTimes_snt_prod = [];
% 
% for unit_i = 1 : numel(rasterData.sentence)
%     trlTimes_snt_audio = [trlTimes_snt_audio; rasterData.(level)(unit_i).trialTimes.stimulus];
%     trlTimes_snt_prod = [trlTimes_snt_prod; rasterData.(level)(unit_i).trialTimes.speech];
% end
% 
% Evt_audio = [median(trlTimes_snt_prod.cue_offset - trlTimes_snt_prod.prod_onset) 0 median(trlTimes_snt_prod.prod_offset - trlTimes_snt_prod.prod_onset)];
% Evt_prod = [median(trlTimes_snt_audio.audio_offset - trlTimes_snt_audio.audio_onset) 0 median(trlTimes_snt_prod.prod_onset - trlTimes_snt_prod.cue_onset)];
% timeEvents =[0 Evt_audio(3) Evt_prod(3)];
%% v0 build population matrix Xaverage e X
% need to build XAverage: N x S x D x T
%                      X: N x S x D x T x maxTrialNum
% N is the number of neurons/LFP channels
% S is the number of task conditions (2: nonoise and lombard)
% D is the number of part of the task (2: listening and speaking)
% T is the number of time-points (note that all the trials should have the
% same length in time!)
%
% cc = 0;
% 
% % loop over units
% 
% Xavg = [];
% X = [];
% trialNum_all = [];
% for unit_i = 1 : nUnits
%     % get current unit
%     rasterData_unit =  rasterData.(level)(unit_i);
%     fprintf('Getting firing rate in Unit %d (%d - %d) [%1.f%%] \n',rasterData_unit.unitID,unit_i,nUnits,unit_i/nUnits*100)
% 
%     if rasterData_unit.nTrials.nonoise >= MIN_TRIALS && rasterData_unit.nTrials.lombard >= MIN_TRIALS
%         cc = cc + 1;
% 
%         %if r2data_snt_audio.rlombard_mean(cc) >= 0.05 || r2data_snt_prod.rlombard_mean(cc) >= 0.05 % enough encoding?
% 
%         FR_auditory = rasterData_unit.FR(1).data;
%         FR_speaking = rasterData_unit.FR(2).data;
%         time = rasterData_unit.FR(1).time;
% 
%         % get lombard trials
%         lombard_trials_auditory = rasterData_unit.trialTypes.stimulus.lombard;
%         lombard_trials_speech = rasterData_unit.trialTypes.speech.lombard;
% 
%         FR_auditory_nonoise = FR_auditory(~lombard_trials_auditory,:); FR_auditory_nonoise(isnan(sum(FR_auditory_nonoise,2)) | sum(FR_auditory_nonoise,2) == 0,:) = [];
%         FR_speaking_nonoise = FR_speaking(~lombard_trials_speech,:); FR_speaking_nonoise(isnan(sum(FR_speaking_nonoise,2)) | sum(FR_speaking_nonoise,2) == 0,:) = [];
%         FR_auditory_noise = FR_auditory(lombard_trials_auditory,:); FR_auditory_noise(isnan(sum(FR_auditory_noise,2)) | sum(FR_auditory_noise,2) == 0,:) = [];
%         FR_speaking_noise = FR_speaking(lombard_trials_speech,:); FR_speaking_noise(isnan(sum(FR_speaking_noise,2)) | sum(FR_speaking_noise,2) == 0,:) = [];
% 
%         FR_neuron = nan(1,2,2,numel(time),maxTrialNum/2);
%         FR_neuron(1,1,1,:,1 : size(FR_auditory_nonoise,1)) = FR_auditory_nonoise'; % no noise auditory
%         FR_neuron(1,1,2,:,1 : size(FR_speaking_nonoise,1)) = FR_speaking_nonoise'; % no noise speaking
%         FR_neuron(1,2,1,:,1 : size(FR_auditory_noise,1)) = FR_auditory_noise'; % noise auditory
%         FR_neuron(1,2,2,:,1 : size(FR_speaking_noise,1)) = FR_speaking_noise'; % noise speaking
% 
%         trialNum = zeros(1,2,2);
%         trialNum(1,1,1) = sum(all(~isnan(FR_auditory(~lombard_trials_auditory,:)),2));
%         trialNum(1,1,2) = sum(all(~isnan(FR_speaking(~lombard_trials_speech,:)),2));
%         trialNum(1,2,1) = sum(all(~isnan(FR_auditory(lombard_trials_auditory,:)),2));
%         trialNum(1,2,2) = sum(all(~isnan(FR_speaking(lombard_trials_speech,:)),2));
% 
%         % average
%         FR_auditory_nonoise = mean(FR_auditory(~lombard_trials_auditory,:),1,'omitnan');
%         FR_auditory_lombard = mean(FR_auditory(lombard_trials_auditory,:),1,'omitnan');
%         FR_speaking_nonoise = mean(FR_speaking(~lombard_trials_speech,:),1,'omitnan');
%         FR_speaking_lombard = mean(FR_speaking(lombard_trials_speech,:),1,'omitnan');
% 
%         % cocnatenate first conditions here
%         FRavg_neuron = nan(1,2,2,numel(time));
%         FRavg_neuron(1,1,1,:) = FR_auditory_nonoise; % no noise auditory
%         FRavg_neuron(1,1,2,:) = FR_speaking_nonoise; % no noise speaking
%         FRavg_neuron(1,2,1,:) = FR_auditory_lombard; % noise auditory
%         FRavg_neuron(1,2,2,:) = FR_speaking_lombard; % noise speaking
% 
%         % concatenaet across neurons
%         Xavg = cat(1,Xavg,FRavg_neuron);
%         X = cat(1,X,FR_neuron);
%         trialNum_all = cat(1,trialNum_all,trialNum);
%         % check
%         for s = 1 : 2
%             for d = 1 : 2
%                 assert(isempty(find(isnan(X(cc,s,d,:,1:trialNum_all(cc,s,d))), 1)), 'Something is wrong!')
%             end
%         end
% 
%         % end
%     else
%         warning('No enough trials for unit-%d \n',rasterData_unit.unitID)
%     end
% 
% end



%% v1 Load data preprocessed in R


% h5info(PATH_DATA_INPUT)
X = h5read(PATH_DATA_INPUT, '/X');
X_ntrials = h5read(PATH_DATA_INPUT, '/X_ntrials');
Xmeta = []; 
Xmeta.channel_id = h5read(PATH_DATA_INPUT, '/channel_id');


% % trialNumBAD_all = sum(sum(isnan(X), 4)>50, 5);  % trial unavailable if more than half of its time points are nan
% % trialNumBAD_all = sum(any(isnan(X), 4), 5); % trial unavailable all/all timepoints are nan
% missingtrials =  all(isnan(X), 4); %  squeeze(missingtrials(1, 1, 1, :,:))
% trialNumBAD_all = sum(missingtrials, 5); 
% trialNum_all = 40 - trialNumBAD_all; 

% for trials that just have some missing values, we want to fill them in
X = fillmissing(X, "previous", 4); %  aa = squeeze(Xfillin(1, 1, 1, :,:))
X = fillmissing(X, "next", 4); %  aa = squeeze(Xfillin(1, 1, 1, :,:))

% trialNum_all = trialNum_all(:, :, 1); % this is supposed to be neurons x conditions x 1


% some counts are so low we can't perform cross validation
% FILTER FOR neurons/sites that have more than 5 trials
keepidxs = all(all(X_ntrials >= MIN_TRIALS, 2) , 3); 
X = X(keepidxs, :, :, :, :); 
X_ntrials = X_ntrials(keepidxs, :, :); 

Xsz = size(X)
Xavg  = mean(X, 5, 'omitnan'); % average over trials
size(Xavg)

time = linspace(-1, 5.3, 100); 
timeEvents = [0 2 3]; 

if Xsz(3)==1
    % create a dummy dimension with noise
    warning('CREATING DUMMY DIMENSION')
    noise =  rand(Xsz) * mean(X(:),'omitnan') / 100; 
    X(:, :, 2, :, :) = X(:, :, 1, :, :) + noise;
%     X_ntrials(:, 1, 1) = X_ntrials(:, 1); 
    X_ntrials = repmat(X_ntrials, [1 1 2]); 
    Xsz = size(X)
    size(X_ntrials)
end


trialNum_all = X_ntrials; 
%% Sanity check on loaded data
ich = 40; % 40 DM1015, micro_Lc, we know this to be active 
x = squeeze(X(ich, :, 1, :, :)); 
x = mean(x, 3, 'omitnan'); 

figure; plot(x')

%% Data check
for n = 1:size(X,1)
    for s = 1:size(X,2)
        for d = 1:size(X,3)
            assert(isempty(find(isnan(X(n,s,d,:,1:trialNum_all(n,s,d))), 1)), 'Something is wrong!')
        end
    end
end

aa = squeeze(X(n,s,d,:,1:trialNum_all(n,s,d))); 
%% run dPCA
% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; herewe ignore the 1st dimension
% (neurons), i.e. we have the following parameters:
%    1 - NOISE/QUIET
%    2 - low3rd/mid3rd/loud3rd
%    3 - time
% There are three pairwise interactions:
%    [1 3] - lombard/time interaction

%    [2 3] - auditory|task/time interaction
%    [1 2] - lombard/auditory|task interaction
% And one three-way interaction:
%    [1 2 3] - rest

if ~regularize
    % try without regularization
    [W,V,whichMarg] = dpca(Xavg, NUM_COMPS, ...
        'combinedParams', combinedParams);


    explVar = dpca_explainedVariance(Xavg, W, V, ...
        'combinedParams', combinedParams);

    dpca_plot(Xavg, W, V, @dpca_plot_default, ...
        'explainedVar', explVar, ...
        'marginalizationNames', margNames, ...
        'marginalizationColours', margColours, ...
        'whichMarg', whichMarg,                 ...
        'time', time,                        ...
        'timeEvents', timeEvents,               ...
        'timeMarginalization', 3, ...
        'legendSubplot', 16);
 
else
    mkdir(PATH_FIGUREOUTPUT); 
    if ~isfile(fullfile(PATH_FIGUREOUTPUT,'optimalLambdas.mat'))
        optimalLambda = dpca_optimizeLambda(Xavg, X, trialNum_all, ...
            'combinedParams', combinedParams, ...
            'simultaneous', ifSimultaneousRecording, ...
            'numRep', NUM_SUBSAMPLES, ...  % increase this number to ~10 for better accuracy
            'filename', fullfile(PATH_FIGUREOUTPUT,'A06_optimalLambdas.mat'));
    else
        load(fullfile(PATH_FIGUREOUTPUT,'A06_optimalLambdas.mat'))
    end

    [W,V,whichMarg] = dpca(Xavg, NUM_COMPS, ...
        'combinedParams', combinedParams, ...
        'lambda', optimalLambda);

    explVar = dpca_explainedVariance(Xavg, W, V, ...
        'combinedParams', combinedParams);

    dpca_plot(Xavg, W, V, @dpca_plot_default, ...
        'explainedVar', explVar, ...
        'marginalizationNames', margNames, ...
        'marginalizationColours', margColours, ...
        'whichMarg', whichMarg,                 ...
        'time', time,                        ...
        'timeEvents', timeEvents,               ...
        'timeMarginalization', 3,           ...
        'legendSubplot', 16);
end

exportgraphics(gcf(),fullfile(PATH_FIGUREOUTPUT,'A06_dpcaplot.pdf'),'ContentType','vector')

dPCA = struct();
dPCA.combinedParams = combinedParams;
dPCA.optimalLambda = optimalLambda;
dPCA.whichMarg = whichMarg;
dPCA.explVar = explVar;
dPCA.decoder = W;
dPCA.encoder = V;

save(fullfile(PATH_FIGUREOUTPUT,'A06_dPCA_decomposition.mat'),'dPCA')
% ld = load(fullfile(PATH_FIGUREOUTPUT,'A06_dPCA_decomposition.mat'))

%% Run decoding
 
accuracy_filename = fullfile(PATH_FIGUREOUTPUT,'A06_classification_accuracy.mat');

if ~isfile(fullfile(PATH_FIGUREOUTPUT,'A06_classification_accuracy.mat'))

    % split auditory and decoding
    accuracy = dpca_classificationAccuracy(Xavg, X, trialNum_all, ...
        'lambda', optimalLambda, ...
        'combinedParams', combinedParams, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 100, ...        % increase to 100
        'decodingClasses', decodingClasses, ...
        'filename', accuracy_filename);
    size(accuracy)

    accuracyShuffle = dpca_classificationShuffled(X, trialNum_all, ...
        'lambda', optimalLambda, ...
        'combinedParams', combinedParams, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 100, ...        % increase to 100
        'numShuffles', 100, ...
        'decodingClasses', decodingClasses, ... 
        'filename', accuracy_filename);
    size(accuracyShuffle)


%     % split auditory and decoding
%     accuracy = dpca_classificationAccuracy_separated(Xavg, X, trialNum_all, ...
%         'lambda', optimalLambda, ...
%         'combinedParams', combinedParams, ...
%         'simultaneous', ifSimultaneousRecording, ...
%         'numRep', NUM_SUBSAMPLES, ...        % increase to 100
%         'decodingClasses', decodingClasses, ...
%         'filename', accuracy_filename);
%     accuracyShuffle = dpca_classificationShuffled_separated(X, trialNum_all, ...
%         'lambda', optimalLambda, ...
%         'combinedParams', combinedParams, ...
%         'simultaneous', ifSimultaneousRecording, ...
%         'numRep', NUM_SUBSAMPLES, ...        % increase to 100
%         'numShuffles', NUM_SHUFFLES, ...
%         'decodingClasses', decodingClasses, ...% increase to 100 (takes a lot of time)
%         'filename', accuracy_filename);
%     %dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
%     %dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
    
    %%
    componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
    save(accuracy_filename, 'componentsSignif', '-append')


else
    load(fullfile(PATH_FIGUREOUTPUT,'A06_classification_accuracy.mat'))

end
%% plot pca with signif
dpca_plot(Xavg, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'componentsSignif', componentsSignif, ...
    'legendSubplot', 16);

set(gcf,'Position',[100 100 2000 1000])
exportgraphics(gcf(),fullfile(PATH_FIGUREOUTPUT,'A06_dpcaplot_with-sig.pdf'),'ContentType','vector')



%% Get relevant component
%which_componentsSignif = any(componentsSignif,2)';
NeuronEncoding_Lombard = V(:,find(whichMarg == 1,1)); % get the first one
NeuronEncoding_Task = V(:,find(whichMarg == 2 ,1)); % get the first one



%% plot decoding
% load(accuracy_filename);


figure('Position',[200 200 600 300])
% tiledlayout(1,2)
% nexttile(1)
plot(time,squeeze(accuracy(1,1,:,1)),'linewidth',3,'color','k')
hold on
plotShaded(time,[prctile(squeeze(accuracyShuffle(1,:,:,1)),5,2)'; mean(squeeze(accuracyShuffle(1,:,:,1)),2,'omitnan')'; prctile(squeeze(accuracyShuffle(1,:,:,1)),95,2)'],[.7 .7 .7],'-',2);
% xline([median(trlTimes_snt_audio.audio_offset - trlTimes_snt_audio.audio_onset) 0 median(trlTimes_snt_prod.prod_onset - trlTimes_snt_prod.cue_onset)],'r')
xline(timeEvents)

xlabel("Time [s]")
ylabel('Decoder accuracy [a.u.]')
box off
xlim([-1 5])
ylim([0.3 1])
title(CHANNEL_SUBSET)

exportgraphics(gcf(),fullfile(PATH_FIGUREOUTPUT,'A06_02_decoding.pdf'),'ContentType','vector')


% nexttile(2)
% plot(time,squeeze(accuracy(1,1,:,2)),'linewidth',3,'color','k')
% hold on
% plotShaded(time,[prctile(squeeze(accuracyShuffle(1,:,:,2)),5,2)'; mean(squeeze(accuracyShuffle(1,:,:,2)),2,'omitnan')'; prctile(squeeze(accuracyShuffle(1,:,:,2)),95,2)'],[.7 .7 .7],'-',2);
% % xline([median(trlTimes_snt_prod.cue_offset - trlTimes_snt_prod.prod_onset) 0 median(trlTimes_snt_prod.prod_offset - trlTimes_snt_prod.prod_onset)],'r')
% xline(timeEvents)
% xlabel("time [s]")
% ylabel('Decoder accuracy [a.u.]')
% box off
% xlim([-1 4])
% ylim([0.3 1])
% title('speech window')






%% Get angle between components and distribution of encoder weights
DotProduct = dot(NeuronEncoding_Lombard,NeuronEncoding_Task);
% figure
% scatter(NeuronEncoding_Lombard,NeuronEncoding_Task)
[tauK, psp] = corr(NeuronEncoding_Lombard,NeuronEncoding_Task, 'type', 'Kendall');
CosTheta = max(min(DotProduct/(norm(NeuronEncoding_Lombard)*norm(NeuronEncoding_Task)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));
isSign = abs(DotProduct) > 3.3/sqrt(size(NeuronEncoding_Task,1)) & psp<0.001;
figure("Position",[200 200 900 300])
tiledlayout(1,3)
nexttile(1)
scatter(NeuronEncoding_Lombard,NeuronEncoding_Task,15,'filled')
xlabel('dPCA weight_{Lombard}')
ylabel('dPCA weight_{Task}')
title(['\theta = ', sprintf('%1.2f',ThetaInDegrees), '° \tau = ', sprintf('%1.2f',tauK), ' (p = ',sprintf('%1.2f',psp),')']);%1.2f (p = %1.2f)',,tauK,psp),'Interpreter','latex')
box off
nexttile(2)
pd = fitdist(NeuronEncoding_Lombard,'Normal ');
x = -.5 : .01 : .5;
y = pdf(pd,x);
plot(x,y,'LineWidth',2,'color','b')
xlabel('dPCA weight_{Lombard}')
ylabel('pdf')
xlim([-.5 .5])
box off
nexttile(3)
pd = fitdist(NeuronEncoding_Task,'Normal ');
x = -.5 : .01 : .5;
y = pdf(pd,x);
plot(x,y,'LineWidth',2,'color','b')
xlabel('dPCA weight_{Task}')
ylabel('pdf')
xlim([-.5 .5])
box off
%% Does encoding Lombard related to lombard modulation with task id at the single enruon level?


% get sentence model with task identity
cfg = [];
cfg.log_hash = log_hash;
cfg.lock = 'sentence';
cfg.event = 'all';
cfg.model = 'lmb_task';
cfg.path = PATH_RESULTSOUTPUT_FIRINGRATE;
cfg.time = [0 3]; % greater than 0 (try only [0 2.5] s for audio and [0 3] s for speech [otherwise [0 inf] for both)
Model_Sentence = loadLMComponents(cfg);

% get sentence model without task identity
cfg = [];
cfg.log_hash = log_hash;
cfg.lock = 'sentence';
cfg.event = 'prod_onset';
cfg.model = 'lmb';
cfg.path = PATH_RESULTSOUTPUT_FIRINGRATE;
cfg.time = [0 3]; % greater than 0 (try only [0 2.5] s for audio and [0 3] s for speech [otherwise [0 inf] for both)
Model_Sentence_prod = loadLMComponents(cfg);
cfg.event = 'audio_onset';
cfg.time = [0 2.5]; 
Model_Sentence_audio = loadLMComponents(cfg);

% significant and average among these
cfg = [];
cfg.roi = {'time','post'};
r2data_snt_audio = getRvalues_from_ModelsComponents(cfg,Model_Sentence_audio.lmb);
% get word R values
r2data_snt_prod = getRvalues_from_ModelsComponents(cfg,Model_Sentence_prod.lmb);


% significant and average among these
cfg = [];
cfg.roi = {'time'};
r2data_snt = getRvalues_from_ModelsComponents(cfg,Model_Sentence.lmb_task  );


figure("Position",[200 200 600 300])
tiledlayout(1,2)
nexttile(1)
scatter(r2data_snt_audio.rlombard_max, r2data_snt_prod.rlombard_max,abs(NeuronEncoding_Lombard)*500,abs(NeuronEncoding_Lombard),'filled');
xlabel('\DeltaR^{2} Lombard [auditory]')
ylabel('\DeltaR^{2} Lombard [speech]')
% [Rho,pRho] = compare_stat_scatterplot(r2data_snt_audio.rlombard_post,r2data_snt_audio.rlombard_max,'corr');
% title(sprintf("sentence auditory model (R = %1.2f, p = %1.3f)",Rho,pRho))

colormap(linspecer)
xlim([0 .6])
ylim([0 .6])
box off
nexttile(2)
scatter(r2data_snt.rlombard_max, abs(NeuronEncoding_Lombard),15,[.7 .7 .7],'filled');
hold on
% overlay task invariant and task-variant
scatter(r2data_snt.rlombard_max(r2data_snt.rlombard_combsign), abs(NeuronEncoding_Lombard(r2data_snt.rlombard_combsign)),15,'r','filled');
scatter(r2data_snt.rlombard_max(~r2data_snt.rlombard_combsign & (r2data_snt_audio.rfull_sign | r2data_snt_prod.rfull_sign )), abs(NeuronEncoding_Lombard(~r2data_snt.rlombard_combsign & (r2data_snt_audio.rfull_sign | r2data_snt_prod.rfull_sign ))),15,'b','filled');

%scatter(r2data_snt.rlombard_max, abs(NeuronEncoding_Lombard),[.7 .7 .7],'filled');

xlabel('\DeltaR^{2} Lombard [auditory|speech]')
ylabel('dPCA weight_{Lombard}')
[Rho,pRho] = compare_stat_scatterplot(r2data_snt.rlombard_max,abs(NeuronEncoding_Lombard),'corr');
title(sprintf("(R = %1.2f, p = %1.3f)",Rho,pRho))

%% get all angles

DotProduct = dot(NeuronEncoding_Lombard,r2data_snt_audio.rlombard_max);
% scatter(NeuronEncoding_Lombard,NeuronEncoding_Task)
[tauK, psp] = corr(NeuronEncoding_Lombard,r2data_snt_audio.rlombard_max, 'type', 'Kendall');
CosTheta = max(min(DotProduct/(norm(NeuronEncoding_Lombard)*norm(r2data_snt_audio.rlombard_max)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));
disp(ThetaInDegrees)

DotProduct = dot(NeuronEncoding_Lombard,r2data_snt_prod.rlombard_max);
% scatter(NeuronEncoding_Lombard,NeuronEncoding_Task)
[tauK, psp] = corr(NeuronEncoding_Lombard,r2data_snt_prod.rlombard_max, 'type', 'Kendall');
CosTheta = max(min(DotProduct/(norm(NeuronEncoding_Lombard)*norm(r2data_snt_prod.rlombard_max)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));
disp(ThetaInDegrees)

DotProduct = dot(r2data_snt_prod.rlombard_max,r2data_snt_audio.rlombard_max);
% scatter(NeuronEncoding_Lombard,NeuronEncoding_Task)
[tauK, psp] = corr(r2data_snt_prod.rlombard_max,r2data_snt_audio.rlombard_max, 'type', 'Kendall');
CosTheta = max(min(DotProduct/(norm(r2data_snt_prod.rlombard_max)*norm(r2data_snt_audio.rlombard_max)),1),-1);
ThetaInDegrees = real(acosd(CosTheta));
disp(ThetaInDegrees)

