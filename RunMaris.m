function [data_coh_Maris_Z, data_z_statistic_Maris_threshold, all_cluster_sums_Maris, all_clusters_Maris, data_max_cluster_sum_Maris, Monte_coh_Z_Maris,...
    permutation_zthresholds_Maris,thresholded_Permutation_Maris, monte_all_cluster_sums_Maris, monte_all_clusters_Maris, monte_max_cluster_sum_Maris, pvalue] = RunMaris(PFC_AC_Coh_Spectrum_PriorOnly, PFC_AC_Coh_Spectrum_PretoneOnly,...
                                                                                                                                                                         datadir, sessions, Animal, Behavior, Epoch, Channels,...
                                                                                                                                                                         threshold_value, RandomIteration)
%UNTITLED Summary of this function goes here

% Use the coherence spectra of the prior and pretone conditions to calculate the Maris z-statistic for the data 
        
        data_Maris_zstatistic = []; 
       
            for k = 1:length(sessions)
            Coh_test_statistic_a = atanh(abs(PFC_AC_Coh_Spectrum_PriorOnly(k).cohspctrm)) - (1/(PFC_AC_Coh_Spectrum_PriorOnly(k).dof-2)); % derived from Maris et al 2007 Eq 1 
            Coh_test_statistic_b = atanh(abs(PFC_AC_Coh_Spectrum_PretoneOnly(k).cohspctrm)) - (1/(PFC_AC_Coh_Spectrum_PretoneOnly(k).dof-2));
            Coh_test_statistic_c = sqrt((1/(PFC_AC_Coh_Spectrum_PriorOnly(k).dof-2))+(1/(PFC_AC_Coh_Spectrum_PretoneOnly(k).dof-2)));
            Coh_test_statistic_Zavg = (mean(Coh_test_statistic_a,1)-mean(Coh_test_statistic_b, 1))/Coh_test_statistic_c;
            % Coh_test_statistic_Z = Coh_test_statistic_a-Coh_test_statistic_b/Coh_test_statistic_c;                                            % collapse across channel pairs 
            % Coh_test_statistic_Zavg = mean( Coh_test_statistic_Z , 1);                                                                        % collapse across channel pairs, this line was used till 5/1/24, during the implementation of the ChanWise_SpecEval, we had to collapse earlier in the line prior 
            data_Maris_zstatistic = [data_Maris_zstatistic; Coh_test_statistic_Zavg];
            end    
                                                                                        
% Pool Prior and Pretone Trials 

[RandomPartition_1_Size_c, RandomPartition_2_Size_c, All_Trials_c] = PoolPriorPretone(datadir, sessions);

% Using the pooled data, for each permutation divide the data into two subsets and calculate the z-statistic between the two subsets. 

 % if strcmp(Behavior,'Correct') == 1 
    
    Monte_Carlo_Array_Maris   = [];

    for i = 1:RandomIteration
    
        % generate a random partition equal to the number of trials in the condition with more trials 
    
        totalTrials = size(All_Trials_c.trial, 2);
        randomCols_1 = randperm(totalTrials, RandomPartition_1_Size_c); 
        randomRows_1 = randomCols_1';
    
        for r = 1:RandomPartition_1_Size_c
        Random_Partition_1.time(:, r)        = All_Trials_c.time(:, randomCols_1(r));
        Random_Partition_1.trial(:, r)       = All_Trials_c.trial(:, randomCols_1(r));
        Random_Partition_1.sampleinfo(r, :)  = All_Trials_c.sampleinfo(randomRows_1(r), :);
        end
        
        Random_Partition_1.label      = All_Trials_c(1).label;
        Random_Partition_1.fsample    = All_Trials_c(1).fsample;
        
        % generate a second partition with the remainder of the trials. 
        remainingCols = setdiff(1:totalTrials, randomCols_1);
        randomCols_2  = remainingCols(randperm(length(remainingCols)));
        randomRows_2  = randomCols_2'; 
    
        for rr = 1:RandomPartition_2_Size_c
        Random_Partition_2.time(:, rr)        = All_Trials_c.time(:, randomCols_2(rr));
        Random_Partition_2.trial(:, rr)       = All_Trials_c.trial(:, randomCols_2(rr));
        Random_Partition_2.sampleinfo(rr, :)  = All_Trials_c.sampleinfo(randomRows_2(rr), :);
        end
    
        Random_Partition_2.label      = All_Trials_c(1).label;
        Random_Partition_2.fsample    = All_Trials_c(1).fsample;
    
            % non-parametric computation of the cross-spectral density matrix 

                cfg            = [];
                cfg.method     = 'mtmfft';
                cfg.taper      = 'dpss';
                cfg.output     = 'fourier';
                cfg.tapsmofrq  = 4; %bandwidth of 8 Hz (+-4 Hz) smoothing...
                cfg.pad        = 1;
                cfg.foilim     = [1 100];
                
                freq_1        = ft_freqanalysis(cfg, Random_Partition_1);
                freq_2        = ft_freqanalysis(cfg, Random_Partition_2);
                
                    % calculate the coherence spectra betweeen partition 1 and partition 2 
        
                    cfg            = [];
                    cfg.method     = 'coh';
                    cfg.channelcmb = Channels;
                
                    coh_PFC_AC_1 = ft_connectivityanalysis(cfg, freq_1);       % calculate the coherence spectra of partition 1 
                    f   = coh_PFC_AC_1.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
                    [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1);                % get the electrode ID of the 1st and 2nd columns  
                    clear freq_1
            
                    coh_PFC_AC_2 = ft_connectivityanalysis(cfg, freq_2);       % calculate the coherence spectra of partition 2 
                    clear freq_2
                   
                    C_1 = reshape_coherence(coh_PFC_AC_1, eID_a, eID_b);            % reshape partition 1
                    C_2 = reshape_coherence(coh_PFC_AC_2, eID_a, eID_b);            % reshape partition 2
                    clear coh_PFC_AC_1
                    clear coh_PFC_AC_2
                   
                  
                    % Maris et al 2007 z-statistic 
                   
                    Coh_test_statistic_a = atanh(abs(C_1.cohspctrm)) - (1/(C_1.dof-2));                        % taken from Maris et al 2007 Eq 1 
                    Coh_test_statistic_b = atanh(abs(C_2.cohspctrm)) - (1/(C_2.dof-2));
                    Coh_test_statistic_c = sqrt((1/(C_1.dof-2))+(1/(C_2.dof-2)));
                    Coh_test_statistic_Zavg = (mean(Coh_test_statistic_a,1)-mean(Coh_test_statistic_b, 1))/Coh_test_statistic_c;
                    % Coh_test_statistic_Z = Coh_test_statistic_a-Coh_test_statistic_b/Coh_test_statistic_c;
                    % Coh_test_statistic_Zavg = mean( Coh_test_statistic_Z, 1);                                  % this line collaspes across channel-pair comparisons
                    
                    Monte_Carlo_Array_Maris = [Monte_Carlo_Array_Maris; Coh_test_statistic_Zavg];   
         end  


  %   if strcmp(Behavior,'Wrong') == 1 
  % 
  %   Monte_Carlo_Array_Maris   = [];
  % 
  %   for i = 1:RandomIteration
  % 
  %   % generate a random partition equal to the number of trials in the condition with more trials 
  % 
  %       totalTrials = size(All_Trials_w.trial, 2);
  %       randomCols_1 = randperm(totalTrials, RandomPartition_1_Size_w); 
  %       randomRows_1 = randomCols_1';
  % 
  %       for r = 1:RandomPartition_1_Size_w
  %       Random_Partition_1.time(:, r)        = All_Trials_w.time(:, randomCols_1(r));
  %       Random_Partition_1.trial(:, r)       = All_Trials_w.trial(:, randomCols_1(r));
  %       Random_Partition_1.sampleinfo(r, :)  = All_Trials_w.sampleinfo(randomRows_1(r), :);
  %       end
  % 
  %       Random_Partition_1.label      = All_Trials_w(1).label;
  %       Random_Partition_1.fsample    = All_Trials_w(1).fsample;
  % 
  %       remainingCols = setdiff(1:totalTrials, randomCols_1);
  %       randomCols_2  = remainingCols(randperm(length(remainingCols)));
  %       randomRows_2  = randomCols_2'; 
  % 
  %       for rr = 1:RandomPartition_2_Size_w
  %       Random_Partition_2.time(:, rr)        = All_Trials_w.time(:, randomCols_2(rr));
  %       Random_Partition_2.trial(:, rr)       = All_Trials_w.trial(:, randomCols_2(rr));
  %       Random_Partition_2.sampleinfo(rr, :)  = All_Trials_w.sampleinfo(randomRows_2(rr), :);
  %       end
  % 
  %       Random_Partition_2.label      = All_Trials_w(1).label;
  %       Random_Partition_2.fsample    = All_Trials_w(1).fsample;
  % 
  %           % non-parametric computation of the cross-spectral density matrix (slow)
  % 
  %               cfg            = [];
  %               cfg.method     = 'mtmfft';
  %               cfg.taper      = 'hanning';
  %               cfg.output     = 'fourier';
  %               %cfg.tapsmofrq  = 4; %bandwidth of 8 Hz (+-4 Hz) smoothing...
  %               cfg.pad        = 1;
  %               cfg.foilim     = [1 100];
  % 
  %               freq_1        = ft_freqanalysis(cfg, Random_Partition_1);
  %               freq_2        = ft_freqanalysis(cfg, Random_Partition_2);
  % 
  %                   % calculate the coherence spectra in partition 1 and partition 2 
  % 
  %                   cfg            = [];
  %                   cfg.method     = 'coh';
  %                   cfg.channelcmb = Channels;
  % 
  %                   coh_PFC_AC_1 = ft_connectivityanalysis(cfg, freq_1);       % calculate the coherence spectra of partition 1 
  %                   f   = coh_PFC_AC_1.freq(1, 1:100);                         % establish freq range, options: coh_c.freq(1, 1:51)
  %                   [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1);                % get the electrode ID of the 1st and 2nd columns  
  %                   clear freq_1
  % 
  %                   coh_PFC_AC_2 = ft_connectivityanalysis(cfg, freq_2);       % calculate the coherence spectra of partition 2 
  %                   clear freq_2
  % 
  %                   C_1 = reshape_coherence(coh_PFC_AC_1, eID_a, eID_b);            % reshape partition 1
  %                   C_2 = reshape_coherence(coh_PFC_AC_2, eID_a, eID_b);            % reshape partition 2
  %                   clear coh_PFC_AC_1
  %                   clear coh_PFC_AC_2
  % 
  % 
  %                   % Maris et al 2007 z-statistic 
  % 
  %                   Coh_test_statistic_a = atanh(abs(C_1.cohspctrm)) - (1/(C_1.dof-2));                        % taken from Maris et al 2007 Eq 1 
  %                   Coh_test_statistic_b = atanh(abs(C_2.cohspctrm)) - (1/(C_2.dof-2));
  %                   Coh_test_statistic_c = sqrt((1/(C_1.dof-2))+(1/(C_2.dof-2)));
  %                   Coh_test_statistic_Z = Coh_test_statistic_a-Coh_test_statistic_b/Coh_test_statistic_c;
  %                   Coh_test_statistic_Zavg = mean( Coh_test_statistic_Z, 1);                                  % this line collaspes across channel-pair comparisons
  % 
  %                   Monte_Carlo_Array_Maris = [Monte_Carlo_Array_Maris; Coh_test_statistic_Zavg];   
  %    end  
  % end

 % take the data, threshold, and get the max cluster value (the test statistic)

[data_coh_Maris_Z, data_z_statistic_Maris_threshold, all_cluster_sums_Maris, all_clusters_Maris, data_max_cluster_sum_Maris] = getDataClusterMax(data_Maris_zstatistic,threshold_value);

% take each permutation, threshold, and get the max cluster values (the test statistic). Now you will have the monte carlo distribution comprised of a test-statistic from each permutation

[Monte_coh_Z_Maris,permutation_zthresholds_Maris,thresholded_Permutation_Maris, monte_all_cluster_sums_Maris, monte_all_clusters_Maris, monte_max_cluster_sum_Maris] = getPermClusterMax(Monte_Carlo_Array_Maris,threshold_value);


% calculate the p-value for correct trials 

    pvalue.theta = sum(monte_max_cluster_sum_Maris.theta <= data_max_cluster_sum_Maris.theta) / numel(monte_max_cluster_sum_Maris.theta);
    pvalue.alpha = sum(monte_max_cluster_sum_Maris.alpha <= data_max_cluster_sum_Maris.alpha) / numel(monte_max_cluster_sum_Maris.alpha);
    pvalue.beta = sum(monte_max_cluster_sum_Maris.beta <= data_max_cluster_sum_Maris.beta) / numel(monte_max_cluster_sum_Maris.beta);
    pvalue.gamma = sum(monte_max_cluster_sum_Maris.gamma <= data_max_cluster_sum_Maris.gamma) / numel(monte_max_cluster_sum_Maris.gamma);
    pvalue.highGamma = sum(monte_max_cluster_sum_Maris.highGamma <= data_max_cluster_sum_Maris.highGamma) / numel(monte_max_cluster_sum_Maris.highGamma);
    pvalue.all = sum(monte_max_cluster_sum_Maris.all <= data_max_cluster_sum_Maris.all) / numel(monte_max_cluster_sum_Maris.all);
    
    if pvalue.theta > .5 
    pvalue.theta = 1 - pvalue.theta;
    end 
    
    if pvalue.alpha > .5 
    pvalue.alpha = 1 - pvalue.alpha;
    end 
    
    if pvalue.beta > .5 
    pvalue.beta = 1 - pvalue.beta;
    end 
    
    if pvalue.gamma > .5 
    pvalue.gamma = 1 - pvalue.gamma;
    end 
    
    if pvalue.highGamma > .5 
    pvalue.highGamma = 1 - pvalue.highGamma;
    end
    
    if pvalue.all > .5 
    pvalue.all = 1 - pvalue.all;
    end

 end