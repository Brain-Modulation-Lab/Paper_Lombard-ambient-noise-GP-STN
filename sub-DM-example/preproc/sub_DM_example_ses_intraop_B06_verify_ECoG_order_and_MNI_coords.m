

clear;
close all; 

%% Setup 

bml_defaults

SUBJECT = 'DM1033'; 
SESSION = 'intraop'; 

if ispc
    PATH_DBSROOT = 'Y:\DBS'; 
else    
    PATH_DBSROOT = '/Volumes/Nexus4/DBS'; 
end
PATH_ANALYSIS = [PATH_DBSROOT filesep 'groupanalyses/ecog-localizations/20240327-trace-native-space-electrodetsv-PLB'];
cd(PATH_ANALYSIS); 

PATH_FIG         = [PATH_ANALYSIS filesep 'fig'];
PATH_NEWDATA     = [PATH_ANALYSIS filesep 'data'];

PATH_DER         = [PATH_DBSROOT filesep 'derivatives' ];
PATH_DER_SUB     = [PATH_DER filesep 'sub-' SUBJECT ];  
PATH_PREPROC     = [PATH_DER_SUB filesep 'preproc'];
PATH_ANNOT       = [PATH_DER_SUB filesep 'annot'];
PATH_FIELDTRIP   = [PATH_DER_SUB filesep 'fieldtrip'];
PATH_LEADDBS     = [PATH_DER_SUB filesep 'leaddbs'];
PATH_ECOGLOC     = [PATH_DER_SUB filesep 'ecogloc'];
PATH_FREESURFER  = [PATH_DER_SUB filesep 'freesurfer' filesep 'freesurfer'];

addpath('utils\'); 

PATH_MNI_BRAIN = 'C:\Program Files\LeadDBS_Classic\leaddbs\templates\space/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexHiRes.mat';

channels = bml_annot_read_tsv([PATH_ANNOT filesep bml_bids_basefname(SUBJECT, SESSION) '_channels.tsv']); 

COORD_VARNAMES = {'native_x', 'native_y', 'native_z', 'mni_x', 'mni_y', 'mni_z'};

%% Load ecog localizations  

filename = '_eleclocalizer-all-landmarks-final';
ext = '.tsv'; 
p = [PATH_DER_SUB filesep 'ecogloc' filesep bml_bids_basefname(SUBJECT) filename ext]; 
if exist(p, 'file') 

ecogloc = readtable(p, 'FileType','text', 'Delimiter', '\t') ;
ecogloc = ecogloc(ecogloc.coordframe=="recon"&startsWith(ecogloc.name, 'ecog'), :); 

% PLB 20250606 I frankly don't remember how the _matchleaddbs columns are
% generated, so I'm going to take the solutions from 
% ecogloc.native_x = ecogloc.sol_pos_1_matchleaddbs;
% ecogloc.native_y = ecogloc.sol_pos_2_matchleaddbs; 
% ecogloc.native_z = ecogloc.sol_pos_3_matchleaddbs; 
ecogloc.native_x = ecogloc.sol_pos_1;
ecogloc.native_y = ecogloc.sol_pos_2; 
ecogloc.native_z = ecogloc.sol_pos_3; 

ecogloc = ecogloc(:, [{'name'} COORD_VARNAMES{1:3}]); 


channels_ecog = channels(startsWith(channels.name, 'ecog'), {'starts', 'duration', 'name', 'type'});
% channels_ecog = renamevars(channels_ecog, 'starts', 'onset');

% subject has 8 unlocalized macro electrodes in second script 
if ismember(SUBJECT, {'DM1023'})
    channels_ecog = channels_ecog(~contains(channels_ecog.name, 'ecog_L2'), :);
end

assert(height(ecogloc)==height(channels_ecog)); 


% Take new native coordinates from the ecogloc pipeline
ecogloc = join(channels_ecog, ecogloc, ...
    "Keys", "name", "RightVariables", {'native_x', 'native_y', 'native_z'});

end
%% Load MER/LFP localizations  

filename = '_MER-coords';
ext = '.tsv'; 
p = [PATH_LEADDBS filesep 'annot' filesep bml_bids_basefname(SUBJECT, SESSION) filename ext]; 
% merloc = readtable(p, 'FileType','text', 'Delimiter', '\t') ;
if exist(p, 'file') 

merloc = bml_annot_read_tsv(p);

% build electrode name for every row
merloc.name = repmat({''}, [height(merloc), 1]);
for ir = 1:height(merloc)
    name = ''; 
    type = ''; 
    switch merloc.type{ir}
        case 'MER recording'
            name = [name 'micro_']; 
            type = 'MICRO'; 
        case 'LFP recording'
            name = [name 'macro_']; 
            type = 'SEEG'; 
        otherwise 
            error('type not recognized')
    end
    
    sidevar = find(startsWith(merloc.Properties.VariableNames, 'side'), 1); 
    sidevar = merloc.Properties.VariableNames{sidevar}; 
    disp(['Column name detected to indicate side: ', sidevar])
    switch merloc.(sidevar){ir}
        case 'left'
            name = [name 'L']; 
        case 'right'
            name = [name 'R']; 
        otherwise 
            error('type not recognized')
    end

    name = [name merloc.tract{ir}(1)]; 
    merloc.name{ir} = name; 
    merloc.type{ir} = type; 
end

merloc = merloc(:, [{'starts', 'duration', 'name', 'type'} COORD_VARNAMES]); 

end
%% Load DBS lead locations   
reco = load([PATH_LEADDBS filesep 'ea_reconstruction.mat']);
dbsloc = [reco.reco.native.coords_mm{2} reco.reco.mni.coords_mm{2}; 
          reco.reco.native.coords_mm{1} reco.reco.mni.coords_mm{1}]; 
dbsloc = array2table(dbsloc); 
dbsloc.Properties.VariableNames = COORD_VARNAMES;
writetable(dbsloc, [PATH_LEADDBS filesep 'ea_reconstruction_coords.tsv'], 'FileType','text', 'Delimiter','\t'); 


channels_dbs = channels(startsWith(channels.name, 'dbs'), {'starts', 'duration', 'name', 'type'});

% only left dbs channels are present
if ismember(SUBJECT, {'DM1003', 'DM1029', 'DM1032', 'DM1033'})
    dbsloc = dbsloc(1:height(channels_dbs), :); 
    
% no dbs channels are present in channels.tsv because we never record them
elseif ismember(SUBJECT, {'DM1006', 'DM1009', 'DM1011', 'DM1014', 'DM1015', 'DM1017', 'DM1018'})
    dbsloc(:, :) = []; 

    % subject has different left and right constructions
elseif SUBJECT=="DM1010" 
    reco_L = load([PATH_LEADDBS filesep 'ea_reconstruction_left.mat']);
    reco_R = load([PATH_LEADDBS filesep 'ea_reconstruction_right.mat']);
    dbsloc = [reco_L.reco.native.coords_mm{2} reco_L.reco.mni.coords_mm{2}; 
              reco_R.reco.native.coords_mm{1} reco_R.reco.mni.coords_mm{1}]; 
    dbsloc = array2table(dbsloc); 
    dbsloc.Properties.VariableNames = COORD_VARNAMES;
    writetable(dbsloc, [PATH_LEADDBS filesep 'ea_reconstruction_coords.tsv'], 'FileType','text', 'Delimiter','\t'); 
end

assert(height(channels_dbs)==height(dbsloc), 'DBS channels and DBS electrodes have different number of rows');
dbsloc = [channels_dbs dbsloc];
% dbsloc = renamevars(dbsloc, 'starts', 'onset');

assert(all(dbsloc.native_x(startsWith(dbsloc.name, "dbs_L")) < 0), 'DBS channels and DBS electrodes have different number of rows');
assert(all(dbsloc.native_x(startsWith(dbsloc.name, "dbs_R")) > 0), 'DBS channels and DBS electrodes have different number of rows');


%% Combine electrodes across modalities

concat = {dbsloc}; 
if exist('merloc', 'var'); concat{end+1} = merloc; end
if exist('ecogloc', 'var'); concat{end+1} = ecogloc; end

electrode = bml_annot_rowbind(concat{:}); 
electrode.ends = electrode.starts + electrode.duration; 
electrode = bml_annot_table(electrode); 

%% Transform native coords to mni coords
% taken from Y:\DBS\derivatives\sub-DM-example\preproc\sub_DM_example_ses_intraop_B06_define_electrodes_table_and_MNI_coords.mlx

% DEBUGGING: please ensure that you are using lead_dbs classic
% which ea_setpath
% C:\Program Files\LeadDBS_Classic\leaddbs\ea_setpath.m
assert(contains(which('ea_setpath'), 'LeadDBS_Classic'))

% recalculating mni_nonlinear coords using lead DBS's 
anchor_modality = [PATH_LEADDBS filesep 'anat_t1.nii'];
assert(isfile(anchor_modality),'file %s does not exist',anchor_modality);
warp_file = [PATH_LEADDBS filesep 'glanatInverseComposite.nii.gz'];
assert(isfile(warp_file),'file %s does not exist',warp_file);

leadDBS_mm = [electrode.native_x, electrode.native_y, electrode.native_z]';
leadDBS_mm = [leadDBS_mm; ones(1,size(leadDBS_mm,2))];
[~, leadDBS_vx] = ea_map_coords(leadDBS_mm, anchor_modality);
mni_nonlinear_mm = ea_map_coords(leadDBS_vx, anchor_modality, warp_file,[], 'ANTS',0)';

isCortElec = ismember(electrode.type,{'ECOG'});
electrode(isCortElec,{'mni_x','mni_y','mni_z'}) = ...
    num2cell(mni_nonlinear_mm(isCortElec,1:3));

type= {'ECOG','DBS','SEEG','MICRO'}';
sigma=[1.5,     1,     1,      0.5]';

electrode(:,startsWith(electrode.Properties.VariableNames,'DISTAL'))=[];
cfg=[];
cfg.atlas = 'DISTAL All (Ewert 2017)';
cfg.label_column_basename = 'DISTAL';
cfg.max_assign = 3;
cfg.sigma = bml_map(electrode.type,type,sigma,nan);
cfg.radius = 3;
anat_labels = bml_anat_coord2label(cfg,electrode(:,{'id','mni_x','mni_y','mni_z'}));
electrode = join(electrode,anat_labels);

electrode(:,startsWith(electrode.Properties.VariableNames,'HCPMMP1'))=[];
cfg=[];
cfg.atlas = 'HCPMMP1 (Glasser 2016)';
cfg.label_column_basename = 'HCPMMP1';
cfg.max_assign = 2;
cfg.sigma = bml_map(electrode.type,type,sigma,nan);
cfg.radius = 5;
anat_labels = bml_anat_coord2label(cfg,electrode(:,{'id','mni_x','mni_y','mni_z'}));
electrode = join(electrode, anat_labels);

% electrode = sortrows(electrode, 'name');
electrode = bml_annot_table(electrode);
bml_annot_write_tsv(electrode, [PATH_ANNOT filesep bml_bids_basefname(SUBJECT) '_electrodes.tsv']);


%% Plot mri, pial surface, ecoglocs in NATIVE leaddbs pre-op T1 space 

% lives at Y:\DBS\groupanalyses\ecog-localizations\20240327-trace-native-space-electrodetsv-PLB
% utils/get_leaddbs_opts.m

% contains necessary dependencies
addpath("Y:\DBS\derivatives\sub-DM-example\preproc\utils")

options = get_leaddbs_opts(); 

if exist('atlasname','var')
    options.atlasset=atlasname;
else
    options.atlasset='';
    options.d3.writeatlases=0;
end
options.leadprod='dbs';
options.d3.elrendering=1;
options.d3.exportBB=0;

% renter patient native space space
options.d2.backdrop = 'Patient Pre-OP';
options.native = 1; 
options.patientname = 'leaddbs'; 
options.prefs.patientdir = 'leaddbs'; %change name 


options.root = [PATH_DER_SUB filesep];  %add root, trailing slash is NECESSARY
options.uipatdirs = [PATH_LEADDBS filesep]; %change directory
options.sides = [1 2]; % 1 for right, 2 for left, [1 2] for both

options.sides = [1 2]; % 1 for right, 2 for left, [1 2] for both
if ismember(SUBJECT, {'DM1018'}) % left only
    options.sides = 2; 
elseif ismember(SUBJECT, {''}) % right only
    options.sides = 1;
end

resultfig=ea_elvis(options);
colormap(gray)
hold on
ea_zoomcenter(resultfig.CurrentAxes, [0,0,0], 3);

% % optionally, ROTATE in 3d space and save a couple of views 
% views = []; 
% views.anterior = ea_view(); 
% views.lateral = ea_view(); 
% 
% save('./utils/leaddbs_views', 'views'); 



% Plot pial surface that has been transformed from freesurfer to leaddbs T1
% space
ld = load([PATH_FREESURFER filesep 'cortex_indiv.mat']); 
Vertices = ld.cortex.vert; Faces = ld.cortex.tri;
shg; % delete(Hp); delete(Hc);
Hp = patch('vertices', Vertices,'faces', Faces,...
'facecolor',[.50 .50 .50],'edgecolor','none',...
'facelighting', 'gouraud', 'specularstrength', .50, 'facealpha', 1); Hc = camlight('headlight');
shg

% Plot electrodes 
electrode = bml_annot_read_tsv([PATH_ANNOT filesep bml_bids_basefname(SUBJECT) '_electrodes.tsv']);

idxs_ecog = startsWith(electrode.name, 'ecog'); 
electrode_ecog = electrode(idxs_ecog, :); 
electrode_ecog(:, {'x', 'y', 'z'}) = electrode_ecog(:, {'native_x', 'native_y', 'native_z'});

cfg.h_ax = gca;
cfg.radius = 1;
[h_points] = bml_plot3d_points(cfg, electrode_ecog); 


% plot the hull
ld = load([PATH_FREESURFER filesep 'hull.mat']); 
hull = ld.mask_indices; 

% hHull = plot3(hull(:,1), hull(:,2), hull(:,3),'.');


% Save views
ld = load('./utils/leaddbs_views.mat'); views = ld.views; 

views_str = {'lateral', 'anterior'};

for iview = 1:length(views_str)
    delete((findobj(gca, 'type', 'light')))

    cam_view = views_str{iview};
    ea_view(views.(cam_view));
    camlight('headlight'); camlight('headlight');  % ea_view(views.medial); view(views.medial.az, views.medial.el);
    export_fig([PATH_ECOGLOC filesep bml_bids_basefname(SUBJECT, SESSION)  ...
        '_ecog-leaddbs-native' '_view-' cam_view], '-m3', gcf);

    pause(2)
end

%% Plot ecoglocs in MNI  space 
% Plot MNI cortex
elec = bml_annot_read_tsv([PATH_ANNOT filesep bml_bids_basefname(SUBJECT) '_electrodes.tsv']);

% load standard MNI brain for plotting
ld = load(PATH_MNI_BRAIN, 'Vertices', 'Faces');

cortex = reducepatch(ld.Faces, ld.Vertices, 0.3); 

figure; set(gcf,'Visible','on')

cfg = []; 
cfg.h_ax = gca;
cfg.surface_facealpha = 0.7; 
bml_plot3d_surface(cfg, cortex.vertices, cortex.faces); 
elec.x = elec.mni_x; elec.y = elec.mni_y; elec.z = elec.mni_z;

cfg.h_ax = gca;
cfg.radius = 1;
[h_points] = bml_plot3d_points(cfg, elec); 

% Highlight one electrode
cfg.is_interactive = 1;
elec_subset = elec(contains(elec.name, 'ecog_L114'), :); 
[h_points] = bml_plot3d_points(cfg, elec); 

title({[SUBJECT], ' MNI-space ECoG reconstruction'})

% Save views
ld = load('./utils/leaddbs_views.mat'); views = ld.views; 

views_str = {'lateral', 'anterior'};


for iview = 1:length(views_str)
    delete((findobj(gca, 'type', 'light')))

    cam_view = views_str{iview};
    ea_view(views.(cam_view)); 
    [az, el] = view(); view(az, el+1); shg % have to do this weird reset to get title to show
    camlight('headlight'); % ea_view(views.medial); view(views.medial.az, views.medial.el);
    export_fig([PATH_ECOGLOC filesep bml_bids_basefname(SUBJECT, SESSION)  ...
        '_ecog-leaddbs-mni' '_view-' cam_view], '-m3', gcf);

    pause(2)
end




