%% Defining paths
% add brainstorm to path
% addpath(genpath('E:\MATLAB\brainstorm3'));

%% setting paths
SUBJECT='sub-DM10';
SESSION='intraop'; 
PATH_DATA='Y:/DBS';
PATH_CORTEX = [PATH_DATA '/derivatives/' SUBJECT '/freesurfer/freesurfer'];
PATH_LOC = [PATH_DATA '/derivatives/' SUBJECT '/freesurfer/freesurfer/Electrode_Locations'];
PATH_PREPROC = [PATH_DATA '/derivatives/' SUBJECT '/preproc'];
PATH_FREESURFER = [PATH_DATA '/derivatives/' SUBJECT '/freesurfer/freesurfer'];
PATH_MNI_BRAIN = 'C:/MATLAB_external_libs/lead_v2.5.3/templates/space/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexHiRes.mat';
cd(PATH_PREPROC)

% loading cortex surface
load([PATH_CORTEX '/cortex_indiv.mat'])

%% === Step 1: verifying numbering on LEFT ecog strips === %%
LOC_FILE = [PATH_LOC '/CortElecLocL_eq.mat'];

% loading electrode locations on Left side
load(LOC_FILE)
elec = reshape(cell2mat(CortElecLoc),3,length(CortElecLoc))';

%% Plotting
%Displaying cortex
figure;
Hp = patch('Vertices',cortex.vert,'Faces',cortex.tri,...
    'facecolor',[1 1 1],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');
axis off; axis equal
alpha 1

% plot electrodes on cortical surface
hold on
plot3(elec(:,1), elec(:,2), elec(:,3), 'o', 'color', 'r', 'MarkerSize', 15)
text(elec(:,1), elec(:,2), elec(:,3), num2str((1:size(elec,1))')); 

exportgraphics(gcf,['fig/ecog_localization/sub-' SUBJECT '_ses-introap_plot-01-native-cortex-ecog-strips.png'],'Resolution',600); 


%%
% %If Electrodes locations are in incorrect, reorder
% %defining electrode geometry
% Nrow=3;
% NperRow=21;
% Nstrips = 2;
% % 
% %choose one of the following reordering schemes, or create a new one 
%(they only work for 1 electrode at a time)
% 
% % % - the following line 'flips' each row in the proximal-distal direction
% % reorder = (ceil((1:(Nrow*NperRow))/NperRow).*NperRow-(1:(Nrow*NperRow))) ...
% %    + 1 + floor((0:(Nrow*NperRow-1))/NperRow).*NperRow;
% 
% % % - the following line 'rotates' the order 180 deg
% %reorder = linspace(Nrow*NperRow, 1, Nrow*NperRow);
% 
% % - the following lines reorder the electrodes if both electrodes are
% % labeled from the left corner instead of the right corner when using
% % surgeon's view.
% reorder = zeros(1,Nrow*NperRow*Nstrips);
% for i = 1:Nstrips
%     for j = 1:Nrow
%         reorder(1,((i-1)*Nrow*NperRow + (j-1)*NperRow + 1):((i-1)*Nrow*NperRow + (j-1)*NperRow + NperRow)) = ((i-1)*Nrow*NperRow + (Nrow - j)*NperRow + 1):((i-1)*Nrow*NperRow + (Nrow - j)*NperRow + NperRow);
%     end
% end
% 
% % - the following lines switch the order of the electrodes themselves
% reorder = zeros(1,Nrow*NperRow*Nstrips);
% for i = 1:Nstrips
%         reorder(1,((i-1)*Nrow*NperRow + 1):(i*Nrow*NperRow)) = ((Nstrips - i)*Nrow*NperRow + 1):((Nstrips - i + 1)*Nrow*NperRow);
% end
% %
% %reordering loaded electrode table
% elec = elec(reorder,:);
% % verify re-ordering by re-plotting
% 
% %Correcting electrode order in CortEleLocL_eq.mat
% %save backup of localization mat files in folder
% %<subject>/Anatomy/FreeSurfer/preop/Electrode_Locations 
% %as CortElecLocL_eq_<date modified>.mat and CortElecLocR_eq_<date modified>.mat
% %if available. Move to archive folder (create if required).
% 
% % %reloading objects
% clear('CortElecLoc','CortElecLoc0')
% load(LOC_FILE)
% 
% %correcting 
% CortElecLoc=CortElecLoc(reorder);
% %CortElecLoc0=CortElecLoc0(reorder);
% 
% %saving
% %save(LOC_FILE,'CortElecLoc','CortElecLoc0')
% save(LOC_FILE,'CortElecLoc')
% %restart script and comment this section out

%% Saving electrode coordinates
elec_loc = table((1:size(elec,1))',elec(:,1),elec(:,2),elec(:,3),'VariableNames',{'id','nat_x','nat_y','nat_z'});
writetable(elec_loc,[PATH_LOC '/CortElecLocL_NAT.txt'],'Delimiter','\t')

%%
% The following lines transforms native coordinates to mni space using
% leaddbs. keep a copy of each subjects files with each subject also using
% as a template for subsequent files.

ft_defaults
bml_defaults
format long

PATH_LEADDBS = [PATH_DATA '/derivatives/' SUBJECT '/leaddbs'];

%loading electrode table
electrode = bml_annot_read_tsv([PATH_DATA '/derivatives/' SUBJECT filesep 'annot' filesep SUBJECT '_electrodes.tsv']);

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
electrode = join(electrode,anat_labels);

bml_annot_write_tsv(electrode, [PATH_DATA '/derivatives/' SUBJECT filesep 'annot' filesep SUBJECT '_electrodes.tsv']);


%%
% The following lines transforms electrode coordinates and brain surface
% coordinates from individual space to MNI space
% Template script. Keep one for each subject
% Export brainstorm matlab file to matlab as sMRI before running this script

%% Loading required objects

% load standard MNI brain for plotting
load(PATH_MNI_BRAIN, 'Vertices', 'Faces');

% load cortical surface 
load([PATH_FREESURFER filesep 'cortex_indiv.mat'], 'cortex');

% load T1 image, freesurfer function MRIread.m required
f = MRIread([PATH_FREESURFER filesep 'mri\T1.nii']);

% load in sMRI if not already in the workspace
%load([PATH_FREESURFER filesep 'sMRI_' SUBJECT])

% vox2ras transformation matrix
vox2ras = [1,0,0,-128;0,1,0,-128;0,0,1,-128;0,0,0,1];

%% debugging non-appearance of hcpmmp1 atlas labels

figure; % MNI electrode plotted with MNI brain
Hp = patch('vertices', Vertices,'faces', Faces,...
    'facecolor',[.50 .50 .50],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');
axis off; axis equal
alpha 0.5
hold on; plot3(mni_nonlinear_mm(1:126,1), mni_nonlinear_mm(1:126,2), mni_nonlinear_mm(1:126,3), 'b.', 'MarkerSize', 15)

