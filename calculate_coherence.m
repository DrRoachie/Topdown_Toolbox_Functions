function [] = calculate_coherence(Animal, RecDate, Epoch, Condition)

% set parameters required for fieldtrip functions and data selection

params.choice = choice;
params.err = err;
params.pretone = pretone;
params.pretoneLength = pretoneLength;
params.prior = cell2char(prior); 
params.SNR = SNR;


iSelect = setStimulusCondition(Condition);              

iSelect.err = 'c'; % choose correct trials
data_c = selectData(data,params,iSelect);

iSelect.err = 'w'; % choose wrong trials
data_w = selectData(data,params,iSelect);

% non-parametric computation of the cross-spectral density matrix (slow)

cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = 4; %bandwidth of 8 Hz (+-4 Hz) smoothing...
%cfg.pad       = 'nextpow2';
%cfg.foi       = 0:1/t_length:150; % trial length = 1.4 sec
cfg.foilim     = [0 200];

freq_c        = ft_freqanalysis(cfg, data_c);
freq_w        = ft_freqanalysis(cfg, data_w);
fd_c          = ft_freqdescriptives(cfg,freq_c);
fd_w          = ft_freqdescriptives(cfg,freq_w);

nTrial_c = numel(data_c.trial);
nTrial_w = numel(data_w.trial);
nTaper_c = size(freq_c.fourierspctrm,1) / nTrial_c; % number of tapers
nTaper_w = size(freq_w.fourierspctrm,1) / nTrial_w; % number of tapers
nCh = numel(data_c.label); % number of channel
nFreq = numel(freq_c.freq); % number of points in freq

TrialInfo.nTrial.c = nTrial_c;
TrialInfo.nTrial.w = nTrial_w;
TrialInfo.nTaper = nTaper_c;


% compute coherence
cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'*PFC*','*AC*'};
coh_c         = ft_connectivityanalysis(cfg, freq_c);
coh_w         = ft_connectivityanalysis(cfg, freq_w);



% save data...


    subdir = strcat(Animal,'-',RecDate);
    save_file_name = strcat(RecDate,'_','Coherence_',Epoch,'_',Condition);
    
    if ~isfolder(fullfile(savedir, subdir))
        mkdir(savedir,subdir);
    end
    
    save(fullfile(savedir,save_file_name),'coh_c','coh_w','freq_c','freq_w');
    

end



