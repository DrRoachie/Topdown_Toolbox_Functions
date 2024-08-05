function [data_granger_Maris_Z, data_z_statistic_Maris_threshold, all_cluster_sums_Maris, all_clusters_Maris, data_max_cluster_sum_Maris, Monte_granger_Z_Maris,...
          permutation_zthresholds_Maris,thresholded_Permutation_Maris, monte_all_cluster_sums_Maris, monte_all_clusters_Maris, monte_max_cluster_sum_Maris, ...
          pvalue] = getGrangerMaris(Animal, RecDate, Epoch, Behavior, Frequency_Band, Shared_ChannelPairs, ...
                                    datadir, savedir, sessions, ...
                                    RandomIteration, threshold_value)

% tic;
                                                                                                                                                                         
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

                cfg = [];
                cfg.method = 'mtmfft';
                cfg.taper = 'dpss';
                cfg.output = 'fourier';
                cfg.tapsmofrq = 4;
                cfg.pad = 1;
                cfg.foilim = [0 100];

                freq_1  = ft_freqanalysis(cfg,Random_Partition_1);
                freq_2  = ft_freqanalysis(cfg,Random_Partition_2);

                % Calculate Granger spectra for each channel pair

                    for pr = 1:length(Shared_ChannelPairs)
                        
                        Send_Chan_lb   = Shared_ChannelPairs{pr,1};
                        Rec_Chan_lb    = Shared_ChannelPairs{pr,2};

                        cfg             = [];
                        cfg.method      = 'granger';
                        cfg.channel     = {Send_Chan_lb Rec_Chan_lb};
                        cfg.channelcmb  = {Send_Chan_lb Rec_Chan_lb};

                        G_1             = ft_connectivityanalysis(cfg,freq_1);
                        G_1.dof         = RandomPartition_1_Size_c;

                        G_2             = ft_connectivityanalysis(cfg,freq_2);
                        G_2.dof         = RandomPartition_2_Size_c;

                        % get vector for both directions
                        G_1_Vec_PFC_AC  = squeeze(G_1.grangerspctrm(1,2,:));
                        G_2_Vec_PFC_AC  = squeeze(G_2.grangerspctrm(1,2,:));

                        G_1_Vec_AC_PFC  = squeeze(G_1.grangerspctrm(2,1,:));
                        G_2_Vec_AC_PFC  = squeeze(G_2.grangerspctrm(2,1,:));
                        
                        % Maris et al 2007 z-statistic for both directions
                       
                        Granger_test_statistic_a = atanh(abs(G_1_Vec_PFC_AC)) - (1/(G_1.dof-2));                        % taken from Maris et al 2007 Eq 1 
                        Granger_test_statistic_b = atanh(abs(G_2_Vec_PFC_AC)) - (1/(G_2.dof-2));
                        Granger_test_statistic_c = sqrt((1/(G_1.dof-2))+(1/(G_2.dof-2)));
                        Granger_test_statistic_Zavg = (Granger_test_statistic_b- Granger_test_statistic_a)/Granger_test_statistic_c;
                        
                        Monte_Carlo_Array_Maris_PFC_AC(i,:,pr) = Granger_test_statistic_Zavg';   

                        Granger_test_statistic_a = atanh(abs(G_1_Vec_AC_PFC)) - (1/(G_1.dof-2));                        % taken from Maris et al 2007 Eq 1 
                        Granger_test_statistic_b = atanh(abs(G_2_Vec_AC_PFC)) - (1/(G_2.dof-2));
                        Granger_test_statistic_c = sqrt((1/(G_1.dof-2))+(1/(G_2.dof-2)));
                        Granger_test_statistic_Zavg = (Granger_test_statistic_b- Granger_test_statistic_a)/Granger_test_statistic_c;
                        
                        Monte_Carlo_Array_Maris_AC_PFC(i,:,pr) = Granger_test_statistic_Zavg';   
                    end
    end  

% elapsedTime = toc;
% disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']);

%%
% Define data
Epoch1 = extractBefore(Epoch, 'Onset');
Epoch2 = Epoch;

% Establish file and save directories
datadir = ['/media/arl/SHD/05_Epoc_Cut_2/' Animal '/' Epoch1 '/' RecDate];
addpath(genpath(datadir));
fName = sprintf('%s-%s_bdLFP_%s_ft.mat', Animal, RecDate, Epoch2);

V = load(fName);

% separate prior and pretone data
prior_condition_indices     = find((strcmp(V.prior,'H') | strcmp(V.prior,'L')) & ...      % functionally equivalent to Corey/Taku's iSelect function
                              V.pretone == 'N' & V.choice ~= 'n' & V.err == 'c');
pretone_condition_indices   = find((V.pretone == 'H' | V.pretone =='L') & ...
                            strcmp(V.prior,'N') & V.choice ~= 'n' & V.err == 'c');

prior_selected_trials       = V.trial_id(prior_condition_indices);
pretone_selected_trials     = V.trial_id(pretone_condition_indices);

prior_idx   = find(ismember(V.trial_id,prior_selected_trials)==1);
pretone_idx = find(ismember(V.trial_id,pretone_selected_trials)==1);

prior_data.time         = V.data.time(prior_idx);
prior_data.trial        = V.data.trial(prior_idx);
prior_data.sampleinfo   = V.data.sampleinfo(prior_idx,:);
prior_data.label        = V.data.label;
prior_data.fsample      = V.data.fsample;

pretone_data.time       = V.data.time(pretone_idx);
pretone_data.trial      = V.data.trial(pretone_idx);
pretone_data.sampleinfo = V.data.sampleinfo(pretone_idx,:);
pretone_data.label      = V.data.label;
pretone_data.fsample    = V.data.fsample;

prior_num_trials        = length(prior_data.trial);
pretone_num_trials      = length(pretone_data.trial);

% calculate cross spectral density matrix
cfg             = [];
cfg.method      = 'mtmfft';
cfg.taper       = 'dpss';
cfg.output      = 'fourier';
cfg.tapsmofrq   = 4;
cfg.pad         = 1;
cfg.foilim      = [0 100];

prior_freq      = ft_freqanalysis(cfg,prior_data);
pretone_freq    = ft_freqanalysis(cfg,pretone_data);

% Loop through each shared channel pair
for pr = 1:length(Shared_ChannelPairs)

    Send_Chan_lb   = Shared_ChannelPairs{pr,1};
    Rec_Chan_lb    = Shared_ChannelPairs{pr,2};

    Current_ChanPair_lb = strcat(Send_Chan_lb, '_', Rec_Chan_lb);
    Current_ChanPair_lb = strrep(Current_ChanPair_lb, '*', '');

    cfg             = [];
    cfg.method      = 'granger';
    cfg.channel     = {Send_Chan_lb Rec_Chan_lb};
    cfg.channelcmb  = {Send_Chan_lb Rec_Chan_lb};
    
    Granger_PriorOnly       = ft_connectivityanalysis(cfg,prior_freq);
    Granger_PriorOnly.dof   = prior_num_trials;
    Granger_PretoneOnly     = ft_connectivityanalysis(cfg,pretone_freq);
    Granger_PretoneOnly.dof = pretone_num_trials;

    % get granger vector for both directions
    Granger_Vec_PriorOnly_PFC_AC   = squeeze(Granger_PriorOnly.grangerspctrm(1,2,:));
    Granger_Vec_PretoneOnly_PFC_AC = squeeze(Granger_PretoneOnly.grangerspctrm(1,2,:));

    Granger_Vec_PriorOnly_AC_PFC   = squeeze(Granger_PriorOnly.grangerspctrm(2,1,:));
    Granger_Vec_PretoneOnly_AC_PFC = squeeze(Granger_PretoneOnly.grangerspctrm(2,1,:));
    
    % Use the granger spectra of the prior and pretone conditions to calculate the Maris z-statistic for the data 
        
            Granger_test_statistic_a = atanh(abs(Granger_Vec_PriorOnly_PFC_AC)) - (1/(Granger_PriorOnly.dof-2)); % derived from Maris et al 2007 Eq 1 
            Granger_test_statistic_b = atanh(abs(Granger_Vec_PretoneOnly_PFC_AC)) - (1/(Granger_PretoneOnly.dof-2));
            Granger_test_statistic_c = sqrt((1/(Granger_PriorOnly.dof-2))+(1/(Granger_PretoneOnly.dof-2)));
            Granger_test_statistic_Zavg = (Granger_test_statistic_b- Granger_test_statistic_a)/Granger_test_statistic_c;
            
    data_Maris_zstatistic_PFC_AC = Granger_test_statistic_Zavg';

            Granger_test_statistic_a = atanh(abs(Granger_Vec_PriorOnly_AC_PFC)) - (1/(Granger_PriorOnly.dof-2)); % derived from Maris et al 2007 Eq 1 
            Granger_test_statistic_b = atanh(abs(Granger_Vec_PretoneOnly_AC_PFC)) - (1/(Granger_PretoneOnly.dof-2));
            Granger_test_statistic_c = sqrt((1/(Granger_PriorOnly.dof-2))+(1/(Granger_PretoneOnly.dof-2)));
            Granger_test_statistic_Zavg = (Granger_test_statistic_b- Granger_test_statistic_a)/Granger_test_statistic_c;
            
    data_Maris_zstatistic_AC_PFC = Granger_test_statistic_Zavg';
         

% For each channel pair: take the data, threshold, and get the max cluster value (the test statistic)

[data_granger_Maris_Z_PFC_AC, data_z_statistic_Maris_threshold_PFC_AC, all_cluster_sums_Maris_PFC_AC, all_clusters_Maris_PFC_AC, data_max_cluster_sum_Maris_PFC_AC] = ...
    getDataClusterMax(data_Maris_zstatistic_PFC_AC,threshold_value);

[data_granger_Maris_Z_AC_PFC, data_z_statistic_Maris_threshold_AC_PFC, all_cluster_sums_Maris_AC_PFC, all_clusters_Maris_AC_PFC, data_max_cluster_sum_Maris_AC_PFC] = ...
    getDataClusterMax(data_Maris_zstatistic_AC_PFC,threshold_value);

% take each permutation, threshold, and get the max cluster values (the test statistic). Now you will have the monte carlo distribution comprised of a test-statistic from each permutation

[Monte_granger_Z_Maris_PFC_AC,permutation_zthresholds_Maris_PFC_AC,thresholded_Permutation_Maris_PFC_AC, monte_all_cluster_sums_Maris_PFC_AC, monte_all_clusters_Maris_PFC_AC, ...
    monte_max_cluster_sum_Maris_PFC_AC] = getPermClusterMax(Monte_Carlo_Array_Maris_PFC_AC(:,:,pr),threshold_value);

[Monte_granger_Z_Maris_AC_PFC,permutation_zthresholds_Maris_AC_PFC,thresholded_Permutation_Maris_AC_PFC, monte_all_cluster_sums_Maris_AC_PFC, monte_all_clusters_Maris_AC_PFC, ...
    monte_max_cluster_sum_Maris_AC_PFC] = getPermClusterMax(Monte_Carlo_Array_Maris_AC_PFC(:,:,pr),threshold_value);

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