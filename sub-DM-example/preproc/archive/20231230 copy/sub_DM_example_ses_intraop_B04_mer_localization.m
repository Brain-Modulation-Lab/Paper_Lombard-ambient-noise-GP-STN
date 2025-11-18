% this script cannot currently be run on apple silicon machines (M1/M2)
%

% add lead and spm12 to the path
% addpath(genpath('C:\MATLAB_external_libs\spm12'));
% addpath(genpath('C:\MATLAB_external_libs\lead_v2.5.3')); % v2.3

% pull this repository from github: https://github.com/Brain-Modulation-Lab/Lead_MER

% then remove lead's version of this directory and add the one from github
% % rmpath('/Users/netstim/Documents/MATLAB/Toolbox/leaddbs/tools/MER/');
% % rmpath('/Users/ly546/Downloads/leaddbs-classic/tools/MER');
% rmpath('Y:\Users\lbullock\MATLAB_external_libs_Turbo20230907\archive\leaddbs-classic-20231125\tools\MER');
% % addpath('/Volumes/Nexus4/Users/lbullock/MATLAB_external_libs_Turbo20230907/Lead_MER-master/');
% cd('/Volumes/Nexus4/DBS/derivatives/');

% WINDOWS/BML Turbo
rmpath('Y:\Users\lbullock\MATLAB_external_libs_Turbo20230907\archive\leaddbs-classic-20231125\tools\MER');
addpath('Y:/Users/lbullock/MATLAB_external_libs_Turbo20230907/Lead_MER-master');
cd('Y:/DBS/derivatives/');

% manually navigate to the DBS directory on the server
subjectdir = [fullfile(pwd) filesep]; % subject directory
subject = 'sub-DM10XX'; 

leaddbs_folder_name = 'Herrington_leaddbs_refined';

%%
which MERState

%% load marker table
markers_input = readtable([subject filesep leaddbs_folder_name filesep 'annot' filesep 'markers.xlsx']);

%% add LFP recordings 3 mm above each MER recording
idx_MER = find(strcmp('MER recording', markers_input.Type));
for i=1:length(idx_MER)
    markers_input(end+1,:) = markers_input(idx_MER(i),:);
    markers_input.id(end) = height(markers_input);
    markers_input.Type{end} = 'LFP recording';
    markers_input.Depth(end) = markers_input.Depth(end)+3;
end

%% load electrode type
load([subject filesep leaddbs_folder_name filesep 'ea_reconstruction.mat']);
if isempty(reco.props(1).elmodel)
    elmodel = reco.props(2).elmodel;
else
    elmodel = reco.props(1).elmodel;
end

%% Setup the options struct.
options.root = subjectdir;
options.patientname = [subject filesep 'Herrington_leaddbs_refined'];

options.patientname = [subject filesep 'Herrington_leaddbs_refined'];
options.uipatdirs = {fullfile(options.root, options.patientname)};
options.native = true;
options.loadrecoforviz = 1;
options.sides = [1 2]; % 1 is right, 2 is left
options.prefs = ea_prefs;
options.elmodel = elmodel;
options = ea_resolve_elspec(options);

%% Initialize the MERState, give it the options, and set it to a default state
temp1 = MERState();
temp1.setOptions(options);
temp1.clearData();
temp1.setDataToDefaults();


%% Determine implanted lead depths
implanted = markers_input(strcmp('lead implanted', markers_input.Notes),:);
% left
try
    temp1.updateDBSDepth('left', implanted.Depth(strcmp('left', implanted.Side)));
    temp1.updateDBSImplantTrack('left', implanted.Tract{strcmp('left', implanted.Side)});
catch
    fprintf('   Couldn''t find implanted depth for left lead!\n')
end
% right
try
    temp1.updateDBSDepth('right', implanted.Depth(strcmp('right', implanted.Side)));
    temp1.updateDBSImplantTrack('right', implanted.Tract{strcmp('right', implanted.Side)});
catch
    fprintf('   Couldn''t find implanted depth for right lead!\n')
end

%% Add markers from the annotation table
MarkerTypes = fieldnames(MERState.MarkerTypes);
for m=1:height(markers_input)
    TypeField = MarkerTypes{cellfun(@(x) strcmp(markers_input.Type(m), MERState.MarkerTypes.(x)), MarkerTypes)};
    temp1.addMarkerAtDepth(markers_input.id(m), markers_input.Side{m}, markers_input.Tract{m}, MERState.MarkerTypes.(TypeField), markers_input.Notes{m}, markers_input.Depth(m));
end

%% Save MER localization object for later.
temp1.save('y');  % Pass 'y' to overwrite file if it exists. Omit to be prompted.

%% export markers in MNI and native space, join with input table and write out new table.
markers = temp1.exportMarkers();
markers.id = (1:height(markers))';
markers = join(markers_input, markers, 'Keys', 'id', 'KeepOneCopy', markers.Properties.VariableNames);
writetable(markers, [subject filesep leaddbs_folder_name filesep 'annot' filesep subject '_MER_coords'], 'delimiter', '\t');

%% Working with data that were already saved (commandline or GUI)
% temp2 = MERState();
% temp2.Config.root = [fullfile(pwd) filesep];
% temp2.Config.patientname = ['DBS3001' filesep 'Anatomy' filesep 'leaddbs_' subject];
% temp2.load();
% markers = temp2.exportMarkers();


%% Plot markers
% optional
marker_types = unique(markers.Type);

colors = {'r', 'g', 'b', 'c', 'y'};

% load atlas
load('C:\MATLAB_external_libs\lead_v2.3/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL Minimal (Ewert 2017)/atlas_index.mat');
%Set target here:
target = 'STN';

switch (target) 
    case 'STN'
        target_id = 1;
    case 'GPi'
        target_id = 5;
    case 'GPe'
        target_id = 15;
end

h = figure; hold on;
% plot STNs
for side=1:2
    Hp = patch('vertices',atlases.fv{target_id,side}.vertices,'faces',atlases.fv{target_id,side}.faces,...
        'facecolor',[.750 .50 .50],'edgecolor','none',...
        'facelighting', 'gouraud', 'specularstrength', .50);
    camlight('headlight','infinite');
    axis on; axis equal
    alpha 0.5
end 
% plot markers
for m=1:length(marker_types)
    coords_mni = cell2mat(table2cell(markers(strcmp(marker_types(m), markers.Type),11:13)));
    plot3(coords_mni(:,1), coords_mni(:,2), coords_mni(:,3), 'color', colors{m}, 'linestyle', 'none', 'marker', '.', 'markersize', 25)
end
legend([strcat('r',target); strcat('l',target); marker_types])

% save figures
saveas(h, [subject filesep 'leaddbs' filesep 'figures', filesep subject, '_MERlocs.fig'], 'fig');
saveas(h, [subject filesep 'leaddbs' filesep 'figures', filesep subject, '_MERlocs.pdf'], 'pdf');

