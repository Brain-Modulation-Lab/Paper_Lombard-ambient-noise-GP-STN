% 20230606 The Chrabaszcz 2019 Jneuro paper shows much stronger STN
% gamma activity than I am currently seeing in the Lombard data. Why? One
% possible reason is that they used Julian Neuman's SPM-based robust
% averaging and smoothing
% (Recoreding setup seems similar)

%%

% ld = load('test_ft2spm.mat'); 
% D = ld.D; 

addpath("C:\Users\LY546\Downloads\wjn_toolbox-master\wjn_toolbox-master"); 

robust = 1; 
D = wjn_average('test_ft2spm.dat', robust); 

% S = []; 
% S.D = 'test_ft2spm.dat'; 
% 
% S.robust = []; 
% S.robust.ks = 3; 
% S.robust.savew = 1; 
% S.robust.removebad = 0; 
% 
% D = spm_eeg_average_TF(S); 

%%
% Plot trials (all conditions)
cfg = [];
cfg.baseline     = [-3 -4];
cfg.channel      = ch_name;
cfg.title        = tmp.title;
cfg.events       = tmp.events;

bml_singleplotTFR(cfg, ft_freqdescriptives([], tmp.TF));

colormap(flipud(cbrewer2('RdBu')));
hC = colorbar;
cl = get(gca, 'clim');
% caxis(min(abs(cl))*[-1 1]);
caxis(tmp.zlim); % or set zlim
hC.Label.String = tmp.zlabel;
% ylim([4 150]); 

title(gca, cfg.title, 'Interpreter','none');
xlabel("Time [s] wrt " + TLOCK);
ylabel('Frequency [Hz]');
xlim(TOI_PLOT);