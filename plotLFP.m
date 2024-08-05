%% 

function [PFC_LFP_avg AC_LFP_avg] = plotLFP(datadir, Animal, Epoch, Condition, Channels)
% This function takes the output file of the Coherence_Analysis script, and back-calculates the raw the LFP using
% the cross spectral matrix. It outputs the  average LFP of Area 1 (usually
% PFC) and Area 2 (usually AC) in the frequency domain and plots. 

% Input
% datadir   = 'data_directory';         % options: locations of coherence spectrum
% Animal    = 'MrCassius';              % options: 'Mr' then animal name
% Epoch     = 'testToneOnset';          % options: 'preCueOnset', 'moveOnset', & 'testToneOnset';
% Condition = 'OnlyPrior';              % options: 'OnlyPrior', 'OnlyPretone';
% Channels  = 'allChan';                % options: 'allChan', 'PFClowmid_ACupmid','PFCdeep_ACupmid', 'PFCdeep_ACdeep'

addpath(genpath(datadir));
fName = strcat(Animal,'_Coherence','_',Epoch,'_', Condition);
%fName_shuffled = strcat(Animal,'_Coherence','_',Epoch,'_', Condition, '_shuffled');
load(fName);
%load(fName_shuffled);

% determine the number of sessions 

    if strcmp(Channels,'allChan') == 1
    PFC_AC_LFP_sessions = length(PFC_AC_Cross_Spectral_Matrix_c);
    end
    
    if strcmp(Channels,'PFClowmid_ACupmid') == 1
    PFC_AC_LFP_sessions = length(PFClowmid_ACupmid_Cross_Spectral_Matrix_c);
    end
    
    if strcmp(Channels,'PFCdeep_ACupmid') == 1
    PFC_AC_LFP_sessions = length(PFCdeep_ACupmid_Cross_Spectral_Matrix_c);
    end
    
    if strcmp(Channels,'PFCdeep_ACdeep') == 1 
    PFC_AC_LFP_sessions = length(PFCdeep_ACdeep_Cross_Spectral_Matrix_c);
    end

% calculate the LFP spectra from the cross spectral matrix 

PFC_AC_LFP_c = [];

for z = 1:PFC_AC_LFP_sessions

cfg            = [];
cfg.method     = 'mtmfft';
cfg.taper      = 'dpss';
cfg.output     = 'fourier';
cfg.tapsmofrq  = 4; %bandwidth of 8 Hz (+-4 Hz) smoothing...
%cfg.pad       = 'nextpow2';
cfg.pad        = 1;
%cfg.foi       = 0:1/t_length:150; % trial length = 1.4 sec
cfg.foilim     = [0 200];
 
    if strcmp(Channels,'allChan') == 1
    LFP_c = ft_freqdescriptives(cfg, PFC_AC_Cross_Spectral_Matrix_c(z));
    end

    if strcmp(Channels,'PFClowmid_ACupmid') == 1
    LFP_c = ft_freqdescriptives(cfg, PFClowmid_ACupmid_Cross_Spectral_Matrix_c(z));
    end
    
    if strcmp(Channels,'PFCdeep_ACupmid') == 1
    LFP_c = ft_freqdescriptives(cfg, PFCdeep_ACupmid_Cross_Spectral_Matrix_c(z));
    end

    if strcmp(Channels,'PFCdeep_ACdeep') == 1 
    LFP_c = ft_freqdescriptives(cfg, PFCdeep_ACdeep_Cross_Spectral_Matrix_c(z));
    end 

PFC_AC_LFP_c = [PFC_AC_LFP_c; LFP_c];

end

% For each area get all the values for each frequency into an array

PFC_LFP_array = [];
AC_LFP_array  = [];


PFC_AC_LFP_c_sessions = length(PFC_AC_LFP_c); 

for zz = 1:PFC_AC_LFP_c_sessions

    if strcmp(Channels,'allChan') == 1
    PFC_LFP_array = [PFC_LFP_array; PFC_AC_LFP_c(zz).powspctrm(1:20, :)];      
    AC_LFP_array  = [AC_LFP_array;  PFC_AC_LFP_c(zz).powspctrm(21:40, :)];    
    end
    
    if strcmp(Channels,'PFClowmid_ACupmid') == 1
    PFC_LFP_array = [PFC_LFP_array; PFC_AC_LFP_c(zz).powspctrm(11:15, :)];    
    AC_LFP_array  = [AC_LFP_array;  PFC_AC_LFP_c(zz).powspctrm(26:30, :)];    
    end
    
    if strcmp(Channels,'PFCdeep_ACdeep') == 1 
    PFC_LFP_array = [PFC_LFP_array; PFC_AC_LFP_c(zz).powspctrm(16:20, :)];      
    AC_LFP_array  = [AC_LFP_array;  PFC_AC_LFP_c(zz).powspctrm(26:30, :)];      
    end 
    
    if strcmp(Channels,'PFCdeep_ACdeep') == 1 
    PFC_LFP_array = [PFC_LFP_array; PFC_AC_LFP_c(zz).powspctrm(16:20, :)];    
    AC_LFP_array  = [AC_LFP_array;  PFC_AC_LFP_c(zz).powspctrm(36:40, :)]; 
    end

end

% Average across the first dimension to collapse across channels 

PFC_LFP_avg = [mean(PFC_LFP_array, 1)];
AC_LFP_avg  = [mean(AC_LFP_array, 1)];


% plot the average LFP 

hold on
plot(PFC_LFP_avg,'LineWidth',3, 'color', 'g');
xlim([4 100])
plot(AC_LFP_avg, 'LineWidth',3, 'color', 'c');
xlim([4 100])
hold off


end