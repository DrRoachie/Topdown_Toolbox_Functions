
function [PFC_AC_Coh_Spectrum,PFC_AC_Coh_fd, TrialInfo] = getDirectoryCoherence_singletrial(session, Condition, Channels, Behavior)

PFC_AC_Coh_Spectrum  = [];
PFC_AC_Coh_fd = [];
TrialInfo = [];

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    for k = 1:length(session)
    
    baseFileName = session(k).name;
    fullFileName = fullfile(session(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    load(baseFileName)
    pat2 = digitsPattern;
    RecDate = extract(baseFileName, pat2);
    
    % set parameters required for fieldtrip functions and data selection
    
    params.choice = choice;
    params.err = err;
    params.pretone = pretone;
    params.pretoneLength = pretoneLength;
    params.prior = cell2char(prior); 
    params.SNR = SNR;
    
    % separate data into right and wrong trials 
    
    iSelect = setStimulusCondition(Condition);              
    
            if strcmp(Behavior,'Correct') == 1
        
            iSelect.err = 'c'; % choose correct trials
            data_c = selectData(data,params,iSelect);
            
            % non-parametric computation of the cross-spectral density matrix (slow)
            
            cfg            = [];
            cfg.method     = 'mtmfft';
            cfg.taper      = 'dpss';
            cfg.output     = 'fourier';
            cfg.tapsmofrq  = 4; %bandwidth of 8 Hz (+-4 Hz) smoothing...
            %cfg.pad       = 'nextpow2';
            cfg.pad        = 1;
            cfg.foilim     = [1 100];
            
            freq_c        = ft_freqanalysis(cfg, data_c);
            fd_c          = ft_freqdescriptives(cfg,freq_c);
            
            nTrial_c = numel(data_c.trial);
            nTaper_c = size(freq_c.fourierspctrm,1) / nTrial_c; % number of tapers
            nCh = numel(data_c.label);                          % number of channel
            nFreq = numel(freq_c.freq);                         % number of points in freq
            
            TI_c.nTrial     = nTrial_c;
            TI_c.nTaper     = nTaper_c;
             
                cfg                   = [];
                cfg.method            = 'coh';
                cfg.channelcmb        = Channels;          
                coh_PFC_AC_c          = ft_connectivityanalysis(cfg, freq_c);
                
                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_c);                     % get the electrode ID of the 1st and 2nd columns        
                f   = coh_PFC_AC_c.freq(1, 1:100);                              % establish freq range, options: coh_c.freq(1, 1:51)
                C_c = reshape_coherence(coh_PFC_AC_c, eID_a, eID_b);            % reshape correct trials
                
                PFC_AC_Coh_Spectrum          = [PFC_AC_Coh_Spectrum, coh_PFC_AC_c];
                PFC_AC_Coh_fd                = [PFC_AC_Coh_fd, fd_c];
                TrialInfo                    = [TrialInfo, TI_c];
    
            end
   
            if strcmp(Behavior,'Wrong') == 1

            iSelect.err = 'w'; % choose correct trials
            data_w = selectData(data,params,iSelect);

            % non-parametric computation of the cross-spectral density matrix (slow)
            
            cfg            = [];
            cfg.method     = 'mtmfft';
            cfg.taper      = 'dpss';
            cfg.output     = 'fourier';
            cfg.tapsmofrq  = 4; %bandwidth of 8 Hz (+-4 Hz) smoothing...
            %cfg.pad       = 'nextpow2';
            cfg.pad        = 1;
            cfg.foilim     = [1 100];
            
            freq_w        = ft_freqanalysis(cfg, data_w);
            fd_w          = ft_freqdescriptives(cfg,freq_w);
            nTrial_w = numel(data_w.trial);
            nTaper_w = size(freq_w.fourierspctrm,1) / nTrial_w; % number of tapers
            nCh = numel(data_w.label);                          % number of channel
            nFreq = numel(freq_w.freq);                         % number of points in freq
            
           
            TI_w.nTrial     = nTrial_w;
            TI_w.nTaper     = nTaper_w;
             
                cfg                   = [];
                cfg.method            = 'coh';
                cfg.channelcmb        = Channels;          
                coh_PFC_AC_w          = ft_connectivityanalysis(cfg, freq_w);
               
                [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_w);                     % get the electrode ID of the 1st and 2nd columns        
                f   = coh_PFC_AC_w.freq(1, 1:100);                              % establish freq range, options: coh_c.freq(1, 1:51)
                C_w = reshape_coherence(coh_PFC_AC_w, eID_a, eID_b);            % reshape wrong trial
              
               
                PFC_AC_Coh_Spectrum         = [PFC_AC_Coh_Spectrum, coh_PFC_AC_w];
                PFC_AC_Coh_fd                = [PFC_AC_Coh_fd, fd_w];
                TrialInfo                    = [TrialInfo, TI_w];
          
            
           
            end
    end


   