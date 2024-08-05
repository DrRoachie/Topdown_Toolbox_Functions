
function [tfreq_c] = GetSpectrogram(datadir, Animal, RecDate, Epoch, Condition)
%% Housekeeping
% isSave = 0;
% datadir    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\02_ft_Preprocessed\MrCassius\testToneonset';
% SAVE_DIR    = 'C:\Users\Corey Roach\Desktop\2024_01_12_Top_Down_Analysis'; % save file directory
addpath(genpath(datadir));

%% define data to analyze
% Animal    = 'MrCassius';       % Options: 'MrCassius', 'MrM'; 
% RecDate   = '190418';        
% Epoch     = 'testToneOnset';  % Options: 'testToneOnset', 'preCueOnset', or 'moveOnset'
% Condition = 'OnlyPrior';      % Options: stimulus condition

%% load fieldtrip data formatted by FormatLFP_ft_v2.m
fName = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
load(fullfile(datadir, fName));

%% Set time window 
t_length = 2; % make the trial length 2 sec with zero-padding

if strcmp(Epoch,'testToneOnset')
    t_win = 0.15; % 200 ms 
    t_slidingwin = -0.25:0.01:0.6; % 10-ms sliding window
elseif strcmp(Epoch,'preCueOnset')
    t_length = 2; % make the trial length 2 sec with zero-padding
    t_win = 0.20; % 200 ms
    t_slidingwin = 0.00:0.01:0.80; % 10-ms sliding window
elseif strcmp(Epoch,'moveOnset')
    t_win = 0.20; % 200 ms
    t_slidingwin = -0.40:0.01:0.40; % 10-ms sliding window
end

%% set parameters
params.choice = choice;
params.err = err;
params.pretone = pretone;
params.pretoneLength = pretoneLength;
params.prior = cell2char(prior);
params.SNR = SNR;

%%  select data
iSelect = setStimulusCondition(Condition);

% correct trials
iSelect.err = 'c'; % choose correct trials
data_c = selectData(data,params,iSelect);

% wrong trials
% iSelect.err = 'w'; % choose wrong trials
% data_w = selectData(data,params,iSelect);

%% time frequency analysis
cfg               = [];
cfg.method        = 'mtmconvol';
cfg.output        = 'powandcsd';                          % Options:'powandcsd'; 'fourier'
cfg.foi           = 1:1/(t_length/2):100;
cfg.taper         = 'hanning';                            % Options: 'hanning'; 'dpss'
cfg.t_ftimwin     = ones(1,length(cfg.foi)) .* t_win;     % 250-ms sliding window
cfg.toi           = t_slidingwin; 
cfg.keeptrials    = 'no';                                 % Options: 'no'; 'yes'
cfg.pad           = t_length;

tfreq_c       = ft_freqanalysis(cfg, data_c);
% tfreq_w       = ft_freqanalysis(cfg, data_w);

%% get a baseline for the sakes of zscoring power values 

% if strcmp(Epoch,'testToneOnset')
%     t_win = 0.15;
%     t_slidingwin = -0.45:0.01:-0.35; % 10-ms sliding window (baseline period)
% end
% 
% if strcmp(Epoch,'preCueOnset')
%     t_win = 0.20; % 200 ms
%     t_slidingwin = 0.00:0.01:0.80; % 10-ms sliding window
% end
% 
% if strcmp(Epoch,'moveOnset')
%     t_win = 0.20; % 200 ms
%     t_slidingwin = -0.40:0.01:0.40; % 10-ms sliding window
% end
% 
% % set parameters
% params.choice = choice;
% params.err = err;
% params.pretone = pretone;
% params.pretoneLength = pretoneLength;
% params.prior = cell2char(prior);
% params.SNR = SNR;
% 
% % select data
% iSelect     = setStimulusCondition(Condition);
% %iSelect.err = 'valid'; % choose both correct and wrong trial. 
% iSelect.err = 'c';
% data_baseline = selectData(data,params,iSelect);
% 
% % time frequency analysis
% cfg           = [];
% cfg.method    = 'mtmconvol';
% cfg.output    = 'powandcsd'; %'fourier';
% cfg.foi       = 1:1/(t_length/2):100;
% cfg.taper     = 'hanning'; %'dpss';
% cfg.t_ftimwin = ones(1,length(cfg.foi)) .* t_win; 
% cfg.toi       = t_slidingwin; %-
% cfg.keeptrials = 'yes';
% cfg.pad = t_length;
% 
% tfreq_baseline = ft_freqanalysis(cfg, data_baseline);
% 
% pow = tfreq_baseline.powspctrm; % power spectrum
% crs = tfreq_baseline.crsspctrm; % cross spectrum
% 
% % collapse time
% pow = mean(pow,4);
% crs = mean(crs,4);
% 
% % mean and standard deviation across trials
% mean_pow = squeeze(mean(pow,1));
% std_pow  = squeeze(std(pow,[],1));
% mean_crs = squeeze(mean(crs,1));
% std_crs  = squeeze(std(crs,[],1));
% 
% % reorganize data
% tfreq_baseline.powspctrm_mean = mean_pow;
% tfreq_baseline.powspctrm_std  = std_pow;
% tfreq_baseline.crsspctrm_mean = mean_crs;
% tfreq_baseline.crsspctrm_std  = std_crs;
% tfreq_baseline = rmfield(tfreq_baseline,'powspctrm');
% tfreq_baseline = rmfield(tfreq_baseline,'crsspctrm');
% tfreq_baseline = rmfield(tfreq_baseline,'cumtapcnt');
% tfreq_baseline.time      = tfreq_baseline.time([1 end]);
% 


% plot_spectrogram(tfreq_w,eID_ac);
% plot_spectrogram(tfreq_w,eID_pfc);

% %% save data...
% if isSave==1
%     disp('saving...')
% %     save_file_name = 'coherence_Both';
%     subdir = strcat(Animal,'-',RecDate);
%     save_file_name = strcat('TimeFrequency_',Epoch,'_',Condition,'_pad');
% 
%     if ~isfolder(fullfile(SAVE_DIR,subdir))
%         mkdir(SAVE_DIR,subdir);
%     end
%     save(fullfile(SAVE_DIR,subdir,save_file_name), ...
%         'tfreq_c','tfreq_w',...
%         'tcoh_c','tcoh_w', ...
%         'iSelect');
% 
%     clc;
%     disp('Data saved successfully!');
% end
% 
% 
% %%
% 
% savedir = 'C:\Users\Corey Roach\Desktop';
% save_file_name_1 = 'OnlyPrior_TimeSeries';
% fullPath = fullfile(savedir, save_file_name_1);
% save(fullPath, 'OnlyPrior_TimeSeries');
% 
% %%
% 
% savedir = 'C:\Users\Corey Roach\Desktop';
% save_file_name_2 = 'OnlyPretone_TimeSeries';
% fullPath = fullfile(savedir, save_file_name_2);
% save(fullPath, 'OnlyPretone_TimeSeries');