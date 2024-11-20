function [] = RunMaris_v3(Animal, RecDate, Epoch, Shared_ChannelPairs, Frequency_Band, Behavior, ...
                          datadir, savedir, sessions, ...
                          RandomIteration, threshold_value)

tic

%UNTITLED RunMaris outfitted to run on multiple channel pairs. Can only be
% run on one session at a time (i.e. length(sessions) = 1). Coherence and
% granger statistical tests running in parallel to save time.

% extract session information
baseFileName = sessions.name;
fullFileName = fullfile(sessions.folder, baseFileName);
fprintf(1, 'Now reading %s/n', fullFileName);
V = load(fullfile(fullFileName));

% set parameters required for fieldtrip functions and data selection
params.choice = V.choice;
params.err = V.err;
params.pretone = V.pretone;
params.pretoneLength = V.pretoneLength;
params.prior = cell2char(V.prior); 
params.SNR = V.SNR;

% filter data for correct PriorOnly
iSelect = setStimulusCondition('OnlyPrior');
iSelect.err = 'c'; % choose correct trials
data_prior = selectData(V.data,params,iSelect);

% filter data for correct PretoneOnly
iSelect = setStimulusCondition('OnlyPretone');
iSelect.err = 'c'; % choose correct trials
data_pretone = selectData(V.data,params,iSelect);

% non-parametric computation of cross-spectral density matrix 
cfg             = [];
cfg.method      = 'mtmfft';
cfg.taper       = 'dpss';
cfg.output      = 'fourier';
cfg.tapsmofrq   = 4;
cfg.pad         =  2;
cfg.pad         = 1;
cfg.foilim      = [0 100];

freq_prior      = ft_freqanalysis(cfg, data_prior);
freq_pretone    = ft_freqanalysis(cfg, data_pretone);

% calculate coherence spectra for data
cfg                   = [];
cfg.method            = 'coh';
cfg.channelcmb        = Shared_ChannelPairs;      

PFC_AC_Coh_Spectrum_PriorOnly       = ft_connectivityanalysis(cfg, freq_prior);
PFC_AC_Coh_Spectrum_PretoneOnly     = ft_connectivityanalysis(cfg, freq_pretone);

% calculate granger spectra for data
cfg             = [];
cfg.method      = 'granger';
cfg.channel     = Shared_ChannelPairs;
cfg.channelcmb  = Shared_ChannelPairs; % for some reason fieldtrip requires both cfg.channel and cfg.channelcmb to fire properly

Granger_PriorOnly       = ft_connectivityanalysis(cfg, freq_prior);
Granger_PretoneOnly     = ft_connectivityanalysis(cfg, freq_pretone);

Granger_PriorOnly.dof   = length(data_prior.trial);
Granger_PretoneOnly.dof = length(data_pretone.trial);

% calculate z statistic for every channel pair

    Coh_data_Maris_zstatistic               = []; 
    Granger_data_Maris_zstatistic_PFC_AC    = [];
    Granger_data_Maris_zstatistic_AC_PFC    = [];

            for k = 1:length(Shared_ChannelPairs(:,1))

                % Use the coherence spectra of the prior and pretone conditions to calculate the Maris z-statistic for the data 
                Coh_test_statistic_a = atanh(abs(PFC_AC_Coh_Spectrum_PriorOnly.cohspctrm(k,:))) - (1/(PFC_AC_Coh_Spectrum_PriorOnly.dof-2)); % derived from Maris et al 2007 Eq 1 
                Coh_test_statistic_b = atanh(abs(PFC_AC_Coh_Spectrum_PretoneOnly.cohspctrm(k,:))) - (1/(PFC_AC_Coh_Spectrum_PretoneOnly.dof-2));
                Coh_test_statistic_c = sqrt((1/(PFC_AC_Coh_Spectrum_PriorOnly.dof-2))+(1/(PFC_AC_Coh_Spectrum_PretoneOnly.dof-2)));
                Coh_test_statistic_Zavg = (Coh_test_statistic_b - Coh_test_statistic_a)/Coh_test_statistic_c;
                Coh_data_Maris_zstatistic  = [Coh_data_Maris_zstatistic; Coh_test_statistic_Zavg];

                % calculate z statistic for granger in both directions
                Send_Chan_lb    = extractAfter(Shared_ChannelPairs{k,1},'*');
                Rec_Chan_lb     = extractAfter(Shared_ChannelPairs{k,2},'*');
                
                Send_Chan_idx   = find(contains(Granger_PriorOnly.label, Send_Chan_lb));    % find the index of the sending/receiving channels to map onto granger spectrum
                Rec_Chan_idx    = find(contains(Granger_PriorOnly.label, Rec_Chan_lb));     % Granger_PriorOnly.label and Granger_PretoneOnly.label are identical
                
                % PFC to AC
                    Granger_Vec_PriorOnly   = squeeze(Granger_PriorOnly.grangerspctrm(Send_Chan_idx, Rec_Chan_idx, :));
                    Granger_Vec_PretoneOnly = squeeze(Granger_PretoneOnly.grangerspctrm(Send_Chan_idx, Rec_Chan_idx, :));
    
                    Granger_test_statistic_a = atanh(abs(Granger_Vec_PriorOnly)) - (1/(Granger_PriorOnly.dof-2)); % derived from Maris et al 2007 Eq 1 
                    Granger_test_statistic_b = atanh(abs(Granger_Vec_PretoneOnly)) - (1/(Granger_PretoneOnly.dof-2));
                    Granger_test_statistic_c = sqrt((1/(Granger_PriorOnly.dof-2))+(1/(Granger_PretoneOnly.dof-2)));
                    Granger_test_statistic_Zavg = (Granger_test_statistic_b- Granger_test_statistic_a)/Granger_test_statistic_c;
                    Granger_data_Maris_zstatistic_PFC_AC = [Granger_data_Maris_zstatistic_PFC_AC; Granger_test_statistic_Zavg']; % squeeze function to get the granger vector changes shape so taking the transpose reshapes it to the proper dimensions

                % AC to PFC
                    Granger_Vec_PriorOnly   = squeeze(Granger_PriorOnly.grangerspctrm(Rec_Chan_idx, Send_Chan_idx, :));
                    Granger_Vec_PretoneOnly = squeeze(Granger_PretoneOnly.grangerspctrm(Rec_Chan_idx, Send_Chan_idx, :));
    
                    Granger_test_statistic_a = atanh(abs(Granger_Vec_PriorOnly)) - (1/(Granger_PriorOnly.dof-2)); % derived from Maris et al 2007 Eq 1 
                    Granger_test_statistic_b = atanh(abs(Granger_Vec_PretoneOnly)) - (1/(Granger_PretoneOnly.dof-2));
                    Granger_test_statistic_c = sqrt((1/(Granger_PriorOnly.dof-2))+(1/(Granger_PretoneOnly.dof-2)));
                    Granger_test_statistic_Zavg = (Granger_test_statistic_b- Granger_test_statistic_a)/Granger_test_statistic_c;
                    Granger_data_Maris_zstatistic_AC_PFC = [Granger_data_Maris_zstatistic_AC_PFC; Granger_test_statistic_Zavg'];
            end    
    
                                                                                        
% Pool Prior and Pretone Trials 

[RandomPartition_1_Size_c, RandomPartition_2_Size_c, All_Trials_c] = PoolPriorPretone(datadir, sessions);

% Using the pooled data, for each permutation divide the data into two subsets and calculate the z-statistic between the two subsets. 

 % if strcmp(Behavior,'Correct') == 1 
    
    Coh_Monte_Carlo_Array_Maris             = [];
    Granger_Monte_Carlo_Array_Maris_PFC_AC  = [];
    Granger_Monte_Carlo_Array_Maris_AC_PFC  = [];

    for i = 1:RandomIteration
    
        fprintf(1, '\nNow runnning iteration number %s\n', num2str(i));

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
                cfg.foilim     = [0 100];
                
                freq_1        = ft_freqanalysis(cfg, Random_Partition_1);
                freq_2        = ft_freqanalysis(cfg, Random_Partition_2);
                
                    % calculate the coherence spectra betweeen partition 1
                    % and partition 2 for each shared channel pair
        
                    cfg            = [];
                    cfg.method     = 'coh';
                    cfg.channelcmb = Shared_ChannelPairs; % COREY START HERE TOMORROW, WHY ARE WE DOING THE ALL CHANNELS EVERY TIME?
                
                    coh_PFC_AC_1 = ft_connectivityanalysis(cfg, freq_1);       % calculate the coherence spectra of partition 1 
                    [eID_a, eID_b]= getAreaLabel(coh_PFC_AC_1);                % get the electrode ID of the 1st and 2nd columns  
            
                    coh_PFC_AC_2 = ft_connectivityanalysis(cfg, freq_2);       % calculate the coherence spectra of partition 2 
                   
                    C_1 = reshape_coherence(coh_PFC_AC_1, eID_a, eID_b);            % reshape partition 1
                    C_2 = reshape_coherence(coh_PFC_AC_2, eID_a, eID_b);            % reshape partition 2

                    clear coh_PFC_AC_1;
                    clear coh_PFC_AC_2;
                    
                    % calculate granger spectra for each channel pair

                    cfg             = [];
                    cfg.method      = 'granger';
                    cfg.channel     = Shared_ChannelPairs;
                    cfg.channelcmb  = Shared_ChannelPairs;
                    
                    G_1             = ft_connectivityanalysis(cfg, freq_1);
                    G_2             = ft_connectivityanalysis(cfg, freq_2);
                    
                    G_1.dof         = length(Random_Partition_1.trial);
                    G_2.dof         = length(Random_Partition_2.trial);
                  
                    clear freq_1;
                    clear freq_2;

                    % calculate z statistic for each channel pair
                   
                    for k = 1:length(Shared_ChannelPairs(:,1))
                        
                        % calculate z statistic for coherence

                        Coh_test_statistic_a = atanh(abs(C_1.cohspctrm(k,:))) - (1/(C_1.dof-2));                        % taken from Maris et al 2007 Eq 1 
                        Coh_test_statistic_b = atanh(abs(C_2.cohspctrm(k,:))) - (1/(C_2.dof-2));
                        Coh_test_statistic_c = sqrt((1/(C_1.dof-2))+(1/(C_2.dof-2)));
                        Coh_test_statistic_Zavg = (mean(Coh_test_statistic_a,1)-mean(Coh_test_statistic_b, 1))/Coh_test_statistic_c;
                        Coh_Monte_Carlo_Array_Maris(i,:,k) = Coh_test_statistic_Zavg;   

                        
                        % calculate z statistic for granger
                       
                        Send_Chan_lb    = extractAfter(Shared_ChannelPairs{k,1},'*');
                        Rec_Chan_lb     = extractAfter(Shared_ChannelPairs{k,2},'*');
                        
                        Send_Chan_idx   = find(contains(G_1.label, Send_Chan_lb));    % find the index of the sending/receiving channels to map onto granger spectrum
                        Rec_Chan_idx    = find(contains(G_1.label, Rec_Chan_lb));     % G_1.label and G_2.label are identical
                        
                        % PFC to AC
                            Granger_Vec_1   = squeeze(G_1.grangerspctrm(Send_Chan_idx, Rec_Chan_idx, :));
                            Granger_Vec_2   = squeeze(G_2.grangerspctrm(Send_Chan_idx, Rec_Chan_idx, :));
            
                            Granger_test_statistic_a = atanh(abs(Granger_Vec_1)) - (1/(G_1.dof-2)); % derived from Maris et al 2007 Eq 1 
                            Granger_test_statistic_b = atanh(abs(Granger_Vec_2)) - (1/(G_2.dof-2));
                            Granger_test_statistic_c = sqrt((1/(Granger_PriorOnly.dof-2))+(1/(Granger_PretoneOnly.dof-2)));
                            Granger_test_statistic_Zavg = (Granger_test_statistic_b- Granger_test_statistic_a)/Granger_test_statistic_c;
                            Granger_Monte_Carlo_Array_Maris_PFC_AC(i,:,k) = Granger_test_statistic_Zavg'; % squeeze function to get the granger vector changes shape so taking the transpose reshapes it to the proper dimensions
        
                        % AC to PFC
                            Granger_Vec_1   = squeeze(G_1.grangerspctrm(Rec_Chan_idx, Send_Chan_idx, :));
                            Granger_Vec_2   = squeeze(G_2.grangerspctrm(Rec_Chan_idx, Send_Chan_idx, :));
            
                            Granger_test_statistic_a = atanh(abs(Granger_Vec_1)) - (1/(G_1.dof-2)); % derived from Maris et al 2007 Eq 1 
                            Granger_test_statistic_b = atanh(abs(Granger_Vec_2)) - (1/(G_2.dof-2));
                            Granger_test_statistic_c = sqrt((1/(Granger_PriorOnly.dof-2))+(1/(Granger_PretoneOnly.dof-2)));
                            Granger_test_statistic_Zavg = (Granger_test_statistic_b- Granger_test_statistic_a)/Granger_test_statistic_c;
                            Granger_Monte_Carlo_Array_Maris_AC_PFC(i,:,k) = Granger_test_statistic_Zavg'; % squeeze function to get the granger vector changes shape so taking the transpose reshapes it to the proper dimensions

                    end

     end  


  
 % take the data, threshold, and get the max cluster value (the test
 % statistic) for each pair

 for cp = 1:length(Shared_ChannelPairs(:,1))

     Current_ChanPair       = Shared_ChannelPairs(cp,:);
     Current_ChanPair_lb    = strcat(Current_ChanPair{1}, '_', Current_ChanPair{2});
     Current_ChanPair_lb    = strrep(Current_ChanPair_lb, '*', '');

% save data for coherence

[data_coh_Maris_Z_c, data_z_statistic_Maris_threshold_c, all_cluster_sums_Maris_c, all_clusters_Maris_c, data_max_cluster_sum_Maris_c] = getDataClusterMax(Coh_data_Maris_zstatistic(cp,:),threshold_value);

% take each permutation, threshold, and get the max cluster values (the test statistic). Now you will have the monte carlo distribution comprised of a test-statistic from each permutation

[Monte_coh_Z_Maris_c,permutation_zthresholds_Maris_c,thresholded_Permutation_Maris_c, monte_all_cluster_sums_Maris_c, monte_all_clusters_Maris_c, monte_max_cluster_sum_Maris_c] = getPermClusterMax(Coh_Monte_Carlo_Array_Maris(:,:,cp),threshold_value);


% calculate the p-value for correct trials 

    pvalue_c.theta = sum(monte_max_cluster_sum_Maris_c.theta <= data_max_cluster_sum_Maris_c.theta) / numel(monte_max_cluster_sum_Maris_c.theta);
    pvalue_c.alpha = sum(monte_max_cluster_sum_Maris_c.alpha <= data_max_cluster_sum_Maris_c.alpha) / numel(monte_max_cluster_sum_Maris_c.alpha);
    pvalue_c.beta = sum(monte_max_cluster_sum_Maris_c.beta <= data_max_cluster_sum_Maris_c.beta) / numel(monte_max_cluster_sum_Maris_c.beta);
    pvalue_c.gamma = sum(monte_max_cluster_sum_Maris_c.gamma <= data_max_cluster_sum_Maris_c.gamma) / numel(monte_max_cluster_sum_Maris_c.gamma);
    pvalue_c.highGamma = sum(monte_max_cluster_sum_Maris_c.highGamma <= data_max_cluster_sum_Maris_c.highGamma) / numel(monte_max_cluster_sum_Maris_c.highGamma);
    pvalue_c.all = sum(monte_max_cluster_sum_Maris_c.all <= data_max_cluster_sum_Maris_c.all) / numel(monte_max_cluster_sum_Maris_c.all);
    
    if pvalue_c.theta > .5 
    pvalue_c.theta = 1 - pvalue_c.theta;
    end 
    
    if pvalue_c.alpha > .5 
    pvalue_c.alpha = 1 - pvalue_c.alpha;
    end 
    
    if pvalue_c.beta > .5 
    pvalue_c.beta = 1 - pvalue_c.beta;
    end 
    
    if pvalue_c.gamma > .5 
    pvalue_c.gamma = 1 - pvalue_c.gamma;
    end 
    
    if pvalue_c.highGamma > .5 
    pvalue_c.highGamma = 1 - pvalue_c.highGamma;
    end
    
    if pvalue_c.all > .5 
    pvalue_c.all = 1 - pvalue_c.all;
    end

    Plot_ThreshedZ(Frequency_Band, data_coh_Maris_Z_c, data_z_statistic_Maris_threshold_c)                  % plot the z-spectra of the two conditions with threshold
    Plot_MonteCarloDist(Frequency_Band,monte_max_cluster_sum_Maris_c,data_max_cluster_sum_Maris_c)          % plot monte-carlo distibution with data line 

    Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_Coh_zthresh.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
    Fig2_name      = sprintf('%s_%s_%s_%s_%s_%s_Coh_montecarlo.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
    save_file_name = sprintf('%s_%s_%s_%s_%s_%s_Coh_data.mat', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);

    saveas(1, fullfile(savedir, Fig1_name));
    saveas(2, fullfile(savedir, Fig2_name));

    save(fullfile(savedir, save_file_name), 'data_coh_Maris_Z_c', 'data_z_statistic_Maris_threshold_c', 'all_cluster_sums_Maris_c', 'all_clusters_Maris_c', 'data_max_cluster_sum_Maris_c',...
                                                        'threshold_value', 'Monte_coh_Z_Maris_c','permutation_zthresholds_Maris_c','thresholded_Permutation_Maris_c',...
                                                        'monte_all_cluster_sums_Maris_c', 'monte_all_clusters_Maris_c', 'monte_max_cluster_sum_Maris_c', 'pvalue_c','-append');

    close all;


    % save data for granger

    % For each channel pair: take the data, threshold, and get the max cluster value (the test statistic) in both directions

    [data_granger_Maris_Z_PFC_AC, data_z_statistic_Maris_threshold_PFC_AC, all_cluster_sums_Maris_PFC_AC, all_clusters_Maris_PFC_AC, data_max_cluster_sum_Maris_PFC_AC] = ...
        getDataClusterMax(Granger_data_Maris_zstatistic_PFC_AC(cp,:),threshold_value);
    
    [data_granger_Maris_Z_AC_PFC, data_z_statistic_Maris_threshold_AC_PFC, all_cluster_sums_Maris_AC_PFC, all_clusters_Maris_AC_PFC, data_max_cluster_sum_Maris_AC_PFC] = ...
        getDataClusterMax(Granger_data_Maris_zstatistic_AC_PFC(cp,:),threshold_value);
    
    % take each permutation, threshold, and get the max cluster values (the test statistic). Now you will have the monte carlo distribution comprised of a test-statistic from each permutation
    
    [Monte_granger_Z_Maris_PFC_AC,permutation_zthresholds_Maris_PFC_AC,thresholded_Permutation_Maris_PFC_AC, monte_all_cluster_sums_Maris_PFC_AC, monte_all_clusters_Maris_PFC_AC, ...
        monte_max_cluster_sum_Maris_PFC_AC] = getPermClusterMax(Granger_Monte_Carlo_Array_Maris_PFC_AC(:,:,cp),threshold_value);
    
    [Monte_granger_Z_Maris_AC_PFC,permutation_zthresholds_Maris_AC_PFC,thresholded_Permutation_Maris_AC_PFC, monte_all_cluster_sums_Maris_AC_PFC, monte_all_clusters_Maris_AC_PFC, ...
        monte_max_cluster_sum_Maris_AC_PFC] = getPermClusterMax(Granger_Monte_Carlo_Array_Maris_AC_PFC(:,:,cp),threshold_value);
    
    % calculate the p-value for correct trials in both directions

    pvalue_PFC_AC.theta = sum(monte_max_cluster_sum_Maris_PFC_AC.theta <= data_max_cluster_sum_Maris_PFC_AC.theta) / numel(monte_max_cluster_sum_Maris_PFC_AC.theta);
    pvalue_PFC_AC.alpha = sum(monte_max_cluster_sum_Maris_PFC_AC.alpha <= data_max_cluster_sum_Maris_PFC_AC.alpha) / numel(monte_max_cluster_sum_Maris_PFC_AC.alpha);
    pvalue_PFC_AC.beta = sum(monte_max_cluster_sum_Maris_PFC_AC.beta <= data_max_cluster_sum_Maris_PFC_AC.beta) / numel(monte_max_cluster_sum_Maris_PFC_AC.beta);
    pvalue_PFC_AC.gamma = sum(monte_max_cluster_sum_Maris_PFC_AC.gamma <= data_max_cluster_sum_Maris_PFC_AC.gamma) / numel(monte_max_cluster_sum_Maris_PFC_AC.gamma);
    pvalue_PFC_AC.highGamma = sum(monte_max_cluster_sum_Maris_PFC_AC.highGamma <= data_max_cluster_sum_Maris_PFC_AC.highGamma) / numel(monte_max_cluster_sum_Maris_PFC_AC.highGamma);
    pvalue_PFC_AC.all = sum(monte_max_cluster_sum_Maris_PFC_AC.all <= data_max_cluster_sum_Maris_PFC_AC.all) / numel(monte_max_cluster_sum_Maris_PFC_AC.all);
    
    if pvalue_PFC_AC.theta > .5 
    pvalue_PFC_AC.theta = 1 - pvalue_PFC_AC.theta;
    end 
    
    if pvalue_PFC_AC.alpha > .5 
    pvalue_PFC_AC.alpha = 1 - pvalue_PFC_AC.alpha;
    end 
    
    if pvalue_PFC_AC.beta > .5 
    pvalue_PFC_AC.beta = 1 - pvalue_PFC_AC.beta;
    end 
    
    if pvalue_PFC_AC.gamma > .5 
    pvalue_PFC_AC.gamma = 1 - pvalue_PFC_AC.gamma;
    end 
    
    if pvalue_PFC_AC.highGamma > .5 
    pvalue_PFC_AC.highGamma = 1 - pvalue_PFC_AC.highGamma;
    end
    
    if pvalue_PFC_AC.all > .5 
    pvalue_PFC_AC.all = 1 - pvalue_PFC_AC.all;
    end

    pvalue_AC_PFC.theta = sum(monte_max_cluster_sum_Maris_AC_PFC.theta <= data_max_cluster_sum_Maris_AC_PFC.theta) / numel(monte_max_cluster_sum_Maris_AC_PFC.theta);
    pvalue_AC_PFC.alpha = sum(monte_max_cluster_sum_Maris_AC_PFC.alpha <= data_max_cluster_sum_Maris_AC_PFC.alpha) / numel(monte_max_cluster_sum_Maris_AC_PFC.alpha);
    pvalue_AC_PFC.beta = sum(monte_max_cluster_sum_Maris_AC_PFC.beta <= data_max_cluster_sum_Maris_AC_PFC.beta) / numel(monte_max_cluster_sum_Maris_AC_PFC.beta);
    pvalue_AC_PFC.gamma = sum(monte_max_cluster_sum_Maris_AC_PFC.gamma <= data_max_cluster_sum_Maris_AC_PFC.gamma) / numel(monte_max_cluster_sum_Maris_AC_PFC.gamma);
    pvalue_AC_PFC.highGamma = sum(monte_max_cluster_sum_Maris_AC_PFC.highGamma <= data_max_cluster_sum_Maris_AC_PFC.highGamma) / numel(monte_max_cluster_sum_Maris_AC_PFC.highGamma);
    pvalue_AC_PFC.all = sum(monte_max_cluster_sum_Maris_AC_PFC.all <= data_max_cluster_sum_Maris_AC_PFC.all) / numel(monte_max_cluster_sum_Maris_AC_PFC.all);
    
    if pvalue_AC_PFC.theta > .5 
    pvalue_AC_PFC.theta = 1 - pvalue_AC_PFC.theta;
    end 
    
    if pvalue_AC_PFC.alpha > .5 
    pvalue_AC_PFC.alpha = 1 - pvalue_AC_PFC.alpha;
    end 
    
    if pvalue_AC_PFC.beta > .5 
    pvalue_AC_PFC.beta = 1 - pvalue_AC_PFC.beta;
    end 
    
    if pvalue_AC_PFC.gamma > .5 
    pvalue_AC_PFC.gamma = 1 - pvalue_AC_PFC.gamma;
    end 
    
    if pvalue_AC_PFC.highGamma > .5 
    pvalue_AC_PFC.highGamma = 1 - pvalue_AC_PFC.highGamma;
    end
    
    if pvalue_AC_PFC.all > .5 
    pvalue_AC_PFC.all = 1 - pvalue_AC_PFC.all;
    end

                % plot and save figures/data for both directions

                Plot_ThreshedZ(Frequency_Band, data_granger_Maris_Z_PFC_AC, data_z_statistic_Maris_threshold_PFC_AC)                                                              % plot the z-spectra of the two conditions with threshold
                Plot_MonteCarloDist(Frequency_Band,monte_max_cluster_sum_Maris_PFC_AC,data_max_cluster_sum_Maris_PFC_AC)                                                          % plot monte-carlo distibution with data line 

                Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_Granger_zthresh_PFC_to_AC.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                Fig2_name      = sprintf('%s_%s_%s_%s_%s_%s_Granger_montecarlo_PFC_to_AC.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                save_file_name = sprintf('%s_%s_%s_%s_%s_%s_Granger_data.mat', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);

                saveas(1, fullfile(savedir, Fig1_name));
                saveas(2, fullfile(savedir, Fig2_name));

                save(fullfile(savedir, save_file_name), 'data_granger_Maris_Z_PFC_AC', 'data_z_statistic_Maris_threshold_PFC_AC', 'all_cluster_sums_Maris_PFC_AC', 'all_clusters_Maris_PFC_AC', ...
                                                        'data_max_cluster_sum_Maris_PFC_AC', 'threshold_value', 'Monte_granger_Z_Maris_PFC_AC','permutation_zthresholds_Maris_PFC_AC', ...
                                                        'thresholded_Permutation_Maris_PFC_AC', 'monte_all_cluster_sums_Maris_PFC_AC', 'monte_all_clusters_Maris_PFC_AC', ...
                                                        'monte_max_cluster_sum_Maris_PFC_AC', 'pvalue_PFC_AC','-append');     % appends to the granger data file that saves in the pairwise connectivity uber

                close all

                Plot_ThreshedZ(Frequency_Band, data_granger_Maris_Z_AC_PFC, data_z_statistic_Maris_threshold_AC_PFC)                                                              % plot the z-spectra of the two conditions with threshold
                Plot_MonteCarloDist(Frequency_Band,monte_max_cluster_sum_Maris_AC_PFC, data_max_cluster_sum_Maris_AC_PFC)                                                          % plot monte-carlo distibution with data line 

                Fig1_name      = sprintf('%s_%s_%s_%s_%s_%s_Granger_zthresh_AC_to_PFC.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);
                Fig2_name      = sprintf('%s_%s_%s_%s_%s_%s_Granger_montecarlo_AC_to_PFC.fig', Animal, RecDate, Epoch, Current_ChanPair_lb, Frequency_Band, Behavior);

                saveas(1, fullfile(savedir, Fig1_name));
                saveas(2, fullfile(savedir, Fig2_name));

                save(fullfile(savedir, save_file_name), 'data_granger_Maris_Z_AC_PFC', 'data_z_statistic_Maris_threshold_AC_PFC', 'all_cluster_sums_Maris_AC_PFC', 'all_clusters_Maris_AC_PFC', ...
                                                        'data_max_cluster_sum_Maris_AC_PFC', 'Monte_granger_Z_Maris_AC_PFC','permutation_zthresholds_Maris_AC_PFC', ...
                                                        'thresholded_Permutation_Maris_AC_PFC', 'monte_all_cluster_sums_Maris_AC_PFC', 'monte_all_clusters_Maris_AC_PFC', ...
                                                        'monte_max_cluster_sum_Maris_AC_PFC', 'pvalue_AC_PFC','-append');     % appends to the granger data file that saves in the pairwise connectivity uber

                close all
 end

 elapsedTime = toc;
 disp(['RunMaris_v3 took ' num2str(elapsedTime) ' seconds.']);

end