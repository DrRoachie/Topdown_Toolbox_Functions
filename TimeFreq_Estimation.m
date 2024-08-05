function [ tfreq, tfreq_scrambled] = TimeFreq_Estimation(Epoch, Behavior, choice, err, pretone, pretoneLength, prior, SNR, Condition, data, data_phase_scrambled)
%UNTITLED Summary of this function goes here
   %   Detailed explanation goes here
   % Set timing and sliding window parameters for the entire epoch 
      

                if strcmp(Epoch,'testToneOnset')
                    t_length     = 2.0; % make the trial length 2 sec with zero-padding
                    t_win        = 0.20; % 200 ms 
                    t_slidingwin = -0.8:0.01:0.60;  % 10-ms sliding window
                
                elseif strcmp(Epoch,'preCueOnset')
                    t_length     = 1.0;             % make the trial length 1 sec with zero-padding
                    t_win        = 0.10; 
                    t_slidingwin = -0.05:0.01:0.85; % 10-ms sliding window
                    %t_slidingwin = 0.01:0.01:0.06; % 10-ms sliding window
                elseif strcmp(Epoch,'moveOnset')
                    t_length     = 1.0;
                    t_win        = 0.10; 
                    t_slidingwin = -0.35:0.01:0.35; % 10-ms sliding window
                
                end
        
                % set parameters
                
                params.choice = choice;
                params.err = err;
                params.pretone = pretone;
                params.pretoneLength = pretoneLength;
                params.prior = cell2char(prior);
                params.SNR = SNR;
            
                %  select data
                
                iSelect = setStimulusCondition(Condition);
                
                % correct trials
                iSelect.err       = 'c'; % choose correct trials
                data_c            = selectData(data,params,iSelect);
                data_scrambled_c  = selectData(data_phase_scrambled,params,iSelect);
                
                % wrong trials
                iSelect.err       = 'w'; % choose wrong trials
                data_w            = selectData(data,params,iSelect);
                data_scrambled_w  = selectData(data_phase_scrambled,params,iSelect);
            
                % This segment calculates the power spectra of every channel and the crossspectral density of all channel combinations for each trial. 
                
                cfg               = [];
                cfg.method        = 'mtmconvol';
                cfg.output        = 'powandcsd';                          % Options:'powandcsd'; 'fourier'
                cfg.foi           = 1:1:100;
                cfg.taper         = 'hanning';                            % Options: 'hanning'; 'dpss'
                cfg.channelcmb    = {'*AC*','*PFC*'};                     % specify channel combination here!
                cfg.t_ftimwin     = ones(1,length(cfg.foi)) .* t_win;     
                cfg.toi           = t_slidingwin; 
                cfg.keeptrials    = 'yes';                                 % Options: 'no'; 'yes'
                cfg.pad           = t_length;
            
                    if strcmp(Behavior,'Correct') == 1 
                        tfreq           = ft_freqanalysis(cfg, data_c);
                        tfreq_scrambled = ft_freqanalysis(cfg, data_scrambled_c);
                    end
            
                    if strcmp(Behavior,'Wrong') == 1 
                        tfreq           = ft_freqanalysis(cfg, data_w);
                        tfreq_scrambled = ft_freqanalysis(cfg, data_scrambled_w);
                    end
        
end