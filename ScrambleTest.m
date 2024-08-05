%function [significant_channels] = ScrambleTest(Condition, Behavior, RecDate, Epoch, data, tfreq, tfreq_scrambled, SAVE_DIR)
function [significant_channels] = ScrambleTest(Epoch, data, tfreq, tfreq_scrambled)
%this function takes the time-freq data series and its matched
%phase-scrambled partner, and runs a two-tailed cluster-corrected test for spectral
%similiarity. It returns a pvalue.
        
        
        % freq_band_list = {'theta'; 'alpha'; 'beta'; 'gamma'; 'highGamma';'fullSpectra'};
        %freq_band_list = {'theta'; 'alpha'; 'beta'};
        freq_band_list = {'theta'};
        nbands = length(freq_band_list);
        significant_channels  = cell(1, nbands); 

        for fb = 1:nbands

            Frequency_Band = freq_band_list{fb};

            % Set timing and sliding window parameters for the entire epoch 
            
                % Do the statistical test between scrambled and unscrambed for selected range 
                % https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/
            
                % prep fieldtrip .cfg 
                
                if   strcmp(Frequency_Band, 'theta')
                     freq_band = [5 8];
                   elseif strcmp(Frequency_Band, 'alpha')
                     freq_band = [8 14];
                   elseif strcmp(Frequency_Band, 'beta')
                     freq_band = [15 30];
                   elseif strcmp(Frequency_Band, 'gamma')
                     freq_band = [31 54];
                   elseif strcmp(Frequency_Band, 'highGamma')
                     freq_band = [66 100];
                   elseif strcmp(Frequency_Band, 'fullSpectra')
                     freq_band = [1 100];
                end
        
                if   strcmp(Epoch,'testToneOnset')
                     TimeWin   = [0 .50];
                elseif strcmp(Epoch,'preCueOnset')
                     TimeWin   = [0 .80]; 
                elseif strcmp(Epoch,'moveOnset')
                     TimeWin   = [0 .3];
                end
            
         % Run the test 
            
                chanNum = length(data.label);
                Chan_SpecEval_temp = cell(chanNum, 1); % Initialize cell array to store results
                
                for ch = 1:chanNum
                
                cfg = [];
                
                % parameters for ft_freqstatistics
                cfg.channel            = data.label{ch};                % Options: Nx1 cell-array with selection of channels (default = 'all')
                cfg.latency            = TimeWin;                       % 
                cfg.frequency          = freq_band;                     % Options: freq band that you care about [begin end] theta
                cfg.avgoverchan        = 'no';                          % Options: 'yes' or 'no' (default = 'no')
                cfg.avgovertime        = 'yes';                         % Options: 'yes' or 'no' (default = 'no')
                cfg.avgoverfreq        = 'no';                          % Options: 'yes' or 'no' (default = 'no')
                cfg.parameter          = 'powspctrm';
                
                ntrials_scrambled                     = size(tfreq_scrambled.powspctrm,1);
                ntrials_data                          = size(tfreq.powspctrm, 1);
                design                                = zeros(1, ntrials_scrambled + ntrials_data);
                design(1, 1:ntrials_data)             = 1;
                design(1, ntrials_data +1:end)        = 2;
                cfg.design                            = design; % design matrix; see https://www.fieldtriptoolbox.org/walkthrough/#non-paired-comparison
                cfg.ivar                              = 1;      % the 1st row codes the independent variable (condition)
                
                % parameters for ft_statistics_montecarlo.m
                cfg.method             = 'montecarlo';                % see function notation ft_freqstatistic 
                cfg.numrandomization   = 300;                        % number of random permutations
                cfg.correctm           = 'cluster';                   % apply multiple-comparison correction, 'no', 'max', cluster', 'tfce', 'bonferroni', 'holm', 'hochberg', 'fdr' (default = 'no')
                cfg.alpha              = 0.05;                        % number, critical value for rejecting the null-hypothesis per tail (default = 0.05)
                cfg.tail               = 0;                           % test the left (-1), right(1) or both(0) tails of the distribution
                cfg.correcttail        = 'prob';                      % see https://www.fieldtriptoolbox.org/faq/why_should_i_use_the_cfg.correcttail_option_when_using_statistics_montecarlo/
                
                % parameters for cluster-based correction
                
                cfg.statistic          = 'indepsamplesT';              % see function notation ft_freqstatistic, between UO (trials) design independent samples
                cfg.clusteralpha       = 0.5;                          % threshold for the sample-specific test, is used for thresholding
                cfg.clusterthreshold   = 'nonparametric_common';       % method for single-sample threshold, 'parametric', 'nonparametric_individual', 'nonparametric_common' (default = 'parametric');
                cfg.clusterstatistic   = 'maxsum';                     % see https://github.com/fieldtrip/fieldtrip/blob/master/ft_statistics_montecarlo.m
                cfg.clustertail        = 0;                            % cfg.tail and cfg.clustertail should be identical
                cfg.neighbors          = [];                           % empty neighbourhood structure, clustering will only be done over frequency and/or time and not over neighbouring channels.
                
                
                [stat] = ft_freqstatistics(cfg, tfreq, tfreq_scrambled);
                Chan_SpecEval_temp{ch} = stat;
                
                end
            
                Chan_SpecEval = struct('output', Chan_SpecEval_temp); % Initialize cell array to store results
    
                % Create a logical array pulling all channels with spectral composition sig different from its phase shuffled counterpart

                neg_prob_array = cell(1, 40);
                pos_prob_array = cell(1, 40);
            
                    for k = 1:numel(Chan_SpecEval)
                        % Check if the variables exist in the output structure
                        if isfield(Chan_SpecEval(k).output, 'negclusters') && ...
                                ~isempty(Chan_SpecEval(k).output.negclusters) && ...
                                isfield(Chan_SpecEval(k).output.negclusters, 'prob')
                            neg_prob_array{k} = Chan_SpecEval(k).output.negclusters.prob;
                        else
                            neg_prob_array{k} = NaN;  % Empty array if field does not exist or is empty
                        end
                    
                        if isfield(Chan_SpecEval(k).output, 'posclusters') && ...
                                ~isempty(Chan_SpecEval(k).output.posclusters) && ...
                                isfield(Chan_SpecEval(k).output.posclusters, 'prob')
                            pos_prob_array{k} = Chan_SpecEval(k).output.posclusters.prob;
                        else
                            pos_prob_array{k} = NaN;  % Empty array if field does not exist or is empty
                        end
                    end
            
                % Combine the neg_prob_array and pos_prob_array into a single matrix
                combined_array = [cell2mat(neg_prob_array); cell2mat(pos_prob_array)];
               
                % Initialize a logical array to store the results
                significant_columns = false(1, 40)';
               
                % Loop through the columns and check if either value is less than 0.05
                for k = 1:size(combined_array, 2)
                    if any(combined_array(:, k) < 0.05)
                        significant_columns(k) = true;
                    end
                end
                
                % Find the indices of significant channels
                significant_indices = find(significant_columns);
                
                % Extract the names of significant channels
                significant_channels{fb} = data.label(significant_indices);

      
        end
                
                % Save the correct output 
                % if isempty(significant_channels)
                %     file_name = fullfile(SAVE_DIR, sprintf('NO_SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band));
                % else
                %     file_name = fullfile(SAVE_DIR, sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band));
                % end
                % 
                % if strcmp(Behavior,'Correct') == 1 || strcmp(Behavior,'Both')==1  
                % save(file_name, 'significant_channels');
                % end 
                % 
                % if strcmp(Behavior,'Wrong') == 1 || strcmp(Behavior,'Both')==1  
                % save(file_name, 'significant_channels');
                % end
        