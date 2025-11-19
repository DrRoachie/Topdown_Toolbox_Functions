%% load data
V                 = load(fullfile('D:\02_Preprocessed\MrCassius\testToneOnset\MrCassius-190330_bdLFP_testToneOnset_ft'));
RecDate           = '190330';
%SNR               = [-11/6, -5/3, 5/3, 11/6, -1.5000, -1.2500, 1.2500, 1.5000]; 
SNR               = [-11/6, -5/3, 5/3, 11/6];

%% Get OnlyPrior Trials 

Condition = 'OnlyPrior'; 

% set parameters required for fieldtrip functions and data selection
             params.choice = V.choice;
             params.err = V.err;
             params.pretone = V.pretone;
             params.pretoneLength = V.pretoneLength;
             params.prior = cell2char(V.prior); 
             params.SNR = V.SNR;
            
             % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 
            
             V = ConvertStim(V, RecDate);
            
             % Adds new field in data structure for OnlyPrior Trials that
             % indicates congruence between target and LED

            if strcmp(Condition,'OnlyPrior') ==1
             
             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.prior));  % Preallocate as a cell array
             
             % Loop through each element in V.prior and V.target

             for i = 1:length(V.prior)
                if strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'H')
                    V.congruency{i} = 'congruent';  % If both are 'H'
                elseif strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'L')
                    V.congruency{i} = 'congruent';  % If both are 'L'
                elseif strcmp(V.prior{i}, 'N')
                    V.congruency{i} = 'neutral';  % If prior is 'N'
                elseif (strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'H'))
                    V.congruency{i} = 'incongruent';  % If prior and target are different
                end
             end
            
             params.congruency = V.congruency;

            end

            % Adds new field in data structure for OnlyPrior Trials that
            % indicates congruence between target and pretone

            if strcmp(Condition,'OnlyPretone') == 1
             
             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.pretone));  % Preallocate as a cell array
             V.pretone = cellstr(V.pretone);
             
             % Loop through each element in V.pretone and V.target

             for i = 1:length(V.pretone)
                if strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'H')
                    V.congruency{i} = 'congruent';  % If both are 'H'
                elseif strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'L')
                    V.congruency{i} = 'congruent';  % If both are 'L'
                elseif strcmp(V.pretone{i}, 'N')
                    V.congruency{i} = 'neutral';  % If prior is 'N'
                elseif (strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'H'))
                    V.congruency{i} = 'incongruent';  % If prior and target are different
                end
             end
            
             params.congruency = V.congruency;

            end

            % Generate Condition 1 
            
            iSelect = setStimulusCondition(Condition);
            iSelect.err = 'c'; % choose correct trials
            iSelect.SNR = SNR;     % Options: Selected up to 4 values: -1.8333, -1.6667, -1.5000, -1.2500, 0, 1.2500, 1.5000, 1.6667, 1.8333 
            iSelect.congruency = 'congruent';                     % Options: 'congruent', 'incongruent', and 'neural'
            OnlyPrior_correct_congruent_data = selectData(V.data,params,iSelect);
            
            % Generate Condition 2 
            
            iSelect = setStimulusCondition(Condition);
            iSelect.err = 'w'; % choose correct trials
            iSelect.SNR = SNR;     % Options: Selected up to 4 values: -1.8333, -1.6667, -1.5000, -1.2500, 0, 1.2500, 1.5000, 1.6667, 1.8333 
            iSelect.congruency = 'incongruent';                     % Options: 'congruent', 'incongruent', and 'neural'
            OnlyPrior_wrong_congruent_data = selectData(V.data,params,iSelect);




%% Choose the data that you want to analyze based on parameters 

Condition = 'OnlyPretone';   

% set parameters required for fieldtrip functions and data selection
             params.choice = V.choice;
             params.err = V.err;
             params.pretone = V.pretone;
             params.pretoneLength = V.pretoneLength;
             params.prior = cell2char(V.prior); 
             params.SNR = V.SNR;
            
             % Recode the frequencies of the target at 'H' or 'L'rather than their frequency values. 
            
             V = ConvertStim(V, RecDate);
            
             % Adds new field in data structure for OnlyPrior Trials that
             % indicates congruence between target and LED

            if strcmp(Condition,'OnlyPrior') ==1
             
             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.prior));  % Preallocate as a cell array
             
             % Loop through each element in V.prior and V.target

             for i = 1:length(V.prior)
                if strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'H')
                    V.congruency{i} = 'congruent';  % If both are 'H'
                elseif strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'L')
                    V.congruency{i} = 'congruent';  % If both are 'L'
                elseif strcmp(V.prior{i}, 'N')
                    V.congruency{i} = 'neutral';  % If prior is 'N'
                elseif (strcmp(V.prior{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.prior{i}, 'L') && strcmp(V.target{i}, 'H'))
                    V.congruency{i} = 'incongruent';  % If prior and target are different
                end
             end
            
             params.congruency = V.congruency;

            end

            % Adds new field in data structure for OnlyPrior Trials that
            % indicates congruence between target and pretone

            if strcmp(Condition,'OnlyPretone') == 1
             
             % Initialize V.congruency as a cell array as a new selection factor for select data function  

             V.congruency = cell(size(V.pretone));  % Preallocate as a cell array
             V.pretone = cellstr(V.pretone);
             
             % Loop through each element in V.pretone and V.target

             for i = 1:length(V.pretone)
                if strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'H')
                    V.congruency{i} = 'congruent';  % If both are 'H'
                elseif strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'L')
                    V.congruency{i} = 'congruent';  % If both are 'L'
                elseif strcmp(V.pretone{i}, 'N')
                    V.congruency{i} = 'neutral';  % If prior is 'N'
                elseif (strcmp(V.pretone{i}, 'H') && strcmp(V.target{i}, 'L')) || (strcmp(V.pretone{i}, 'L') && strcmp(V.target{i}, 'H'))
                    V.congruency{i} = 'incongruent';  % If prior and target are different
                end
             end
            
             params.congruency = V.congruency;

            end

            % Generate Condition 3 
            
            iSelect = setStimulusCondition(Condition);
            iSelect.err = 'c'; % choose correct trials
            iSelect.SNR = SNR;     % Options: Selected up to 4 values: -1.8333, -1.6667, -1.5000, -1.2500, 0, 1.2500, 1.5000, 1.6667, 1.8333 
            iSelect.congruency = 'congruent';                     % Options: 'congruent', 'incongruent', and 'neural'
            OnlyPretone_correct_congruent_data = selectData(V.data,params,iSelect);
            
            % Generate Condition 4 
            
            iSelect = setStimulusCondition(Condition);
            iSelect.err = 'w'; % choose correct trials
            iSelect.SNR = SNR;     % Options: Selected up to 4 values: -1.8333, -1.6667, -1.5000, -1.2500, 0, 1.2500, 1.5000, 1.6667, 1.8333 
            iSelect.congruency = 'incongruent';                     % Options: 'congruent', 'incongruent', and 'neural'
            OnlyPretone_wrong_congruent_data = selectData(V.data,params,iSelect);

%% Average PFC channels and Average AC channels 

% Get the number of trials
OnlyPrior_correct_congruent_n          = numel(OnlyPrior_correct_congruent_data.trial);
OnlyPrior_correct_congruent_dataArray  = cat(3, OnlyPrior_correct_congruent_data.trial{:});
OnlyPrior_wrong_congruent_n            = numel(OnlyPrior_wrong_congruent_data.trial);
OnlyPrior_wrong_congruent_dataArray    = cat(3, OnlyPrior_wrong_congruent_data.trial{:});


% Get the number of trials
OnlyPretone_correct_congruent_n          = numel(OnlyPretone_correct_congruent_data.trial);
OnlyPretone_correct_congruent_dataArray  = cat(3, OnlyPretone_correct_congruent_data.trial{:});
OnlyPretone_wrong_congruent_n            = numel(OnlyPretone_wrong_congruent_data.trial);
OnlyPretone_wrong_congruent_dataArray    = cat(3, OnlyPretone_wrong_congruent_data.trial{:});


%%
% Compute the mean across the third dimension (trials)
OnlyPrior_correct_congruent_TrialCollapsed_PFC   = mean(OnlyPrior_correct_congruent_dataArray(1:20, :, :), 3);
OnlyPrior_correct_congruent_ChannelCollapsed_PFC = mean(OnlyPrior_correct_congruent_TrialCollapsed_PFC, 1);
OnlyPrior_correct_congruent_ChannelCollapsed_PFC_std = std(OnlyPrior_correct_congruent_TrialCollapsed_PFC, 1);

OnlyPrior_correct_congruent_TrialCollapsed_AC   = mean(OnlyPrior_correct_congruent_dataArray(21:40, :, :), 3);
OnlyPrior_correct_congruent_ChannelCollapsed_AC = mean(OnlyPrior_correct_congruent_TrialCollapsed_AC, 1);
OnlyPrior_correct_congruent_ChannelCollapsed_AC_std = std(OnlyPrior_correct_congruent_TrialCollapsed_AC, 1);

%%
% Compute the mean across the third dimension (trials)
OnlyPretone_correct_congruent_TrialCollapsed_PFC   = mean(OnlyPretone_correct_congruent_dataArray(1:20, :, :), 3);
OnlyPretone_correct_congruent_ChannelCollapsed_PFC = mean(OnlyPretone_correct_congruent_TrialCollapsed_PFC, 1);
OnlyPretone_correct_congruent_ChannelCollapsed_PFC_std = std(OnlyPretone_correct_congruent_TrialCollapsed_PFC, 1);

OnlyPretone_correct_congruent_TrialCollapsed_AC   = mean(OnlyPretone_correct_congruent_dataArray(21:40, :, :), 3);
OnlyPretone_correct_congruent_ChannelCollapsed_AC = mean(OnlyPretone_correct_congruent_TrialCollapsed_AC, 1);
OnlyPretone_correct_congruent_ChannelCollapsed_AC_std = std(OnlyPretone_correct_congruent_TrialCollapsed_AC, 1);
%%

OnlyPrior_wrong_congruent_TrialCollapsed_PFC   = mean(OnlyPrior_wrong_congruent_dataArray(1:20, :, :), 3);
OnlyPrior_wrong_congruent_ChannelCollapsed_PFC = mean(OnlyPrior_wrong_congruent_TrialCollapsed_PFC, 1);
OnlyPrior_wrong_congruent_ChannelCollapsed_PFC_std = std(OnlyPrior_wrong_congruent_TrialCollapsed_PFC, 1);

OnlyPrior_wrong_congruent_TrialCollapsed_AC   = mean(OnlyPrior_wrong_congruent_dataArray(21:40, :, :), 3);
OnlyPrior_wrong_congruent_ChannelCollapsed_AC = mean(OnlyPrior_wrong_congruent_TrialCollapsed_AC, 1);
OnlyPrior_wrong_congruent_ChannelCollapsed_AC_std = std(OnlyPrior_wrong_congruent_TrialCollapsed_AC, 1);

%%

OnlyPretone_wrong_congruent_TrialCollapsed_PFC   = mean(OnlyPretone_wrong_congruent_dataArray(1:20, :, :), 3);
OnlyPretone_wrong_congruent_ChannelCollapsed_PFC = mean(OnlyPretone_wrong_congruent_TrialCollapsed_PFC, 1);
OnlyPretone_wrong_congruent_ChannelCollapsed_PFC_std = std(OnlyPretone_wrong_congruent_TrialCollapsed_PFC, 1);

OnlyPretone_wrong_congruent_TrialCollapsed_AC   = mean(OnlyPretone_wrong_congruent_dataArray(21:40, :, :), 3);
OnlyPretone_wrong_congruent_ChannelCollapsed_AC = mean(OnlyPretone_wrong_congruent_TrialCollapsed_AC, 1);
OnlyPretone_wrong_congruent_ChannelCollapsed_AC_std = std(OnlyPretone_wrong_congruent_TrialCollapsed_AC, 1);

%% plot mean LFP correct trials 

figure
hold on
ciplot((OnlyPrior_correct_congruent_ChannelCollapsed_PFC(901:1101)+OnlyPrior_correct_congruent_ChannelCollapsed_PFC_std(901:1101)), (OnlyPrior_correct_congruent_ChannelCollapsed_PFC(901:1101)-OnlyPrior_correct_congruent_ChannelCollapsed_PFC_std(901:1101)), 'green');
plot(OnlyPrior_correct_congruent_ChannelCollapsed_PFC(901:1101), 'color', [0 0 0]);
xlim([1 200])

ciplot((OnlyPrior_correct_congruent_ChannelCollapsed_AC(901:1101)+OnlyPrior_correct_congruent_ChannelCollapsed_AC_std(901:1101)), (OnlyPrior_correct_congruent_ChannelCollapsed_AC(901:1101)-OnlyPrior_correct_congruent_ChannelCollapsed_AC_std(901:1101)), 'blue');
plot(OnlyPrior_correct_congruent_ChannelCollapsed_AC(901:1101), 'color', [0 0 0]);
title('PriorOnly Correct')
xlim([1 200])
ylabel('mean LFP (uV)')
xlabel('time (ms)')
hold off

%% plot mean LFP correct trials 

figure
hold on
ciplot((OnlyPretone_correct_congruent_ChannelCollapsed_PFC(901:1101)+OnlyPretone_correct_congruent_ChannelCollapsed_PFC_std(901:1101)), (OnlyPretone_correct_congruent_ChannelCollapsed_PFC(901:1101)-OnlyPretone_correct_congruent_ChannelCollapsed_PFC_std(901:1101)), 'green');
plot(OnlyPretone_correct_congruent_ChannelCollapsed_PFC(901:1101), 'color', [0 0 0]);
xlim([1 200])

ciplot((OnlyPretone_correct_congruent_ChannelCollapsed_AC(901:1101)+OnlyPretone_correct_congruent_ChannelCollapsed_AC_std(901:1101)), (OnlyPretone_correct_congruent_ChannelCollapsed_AC(901:1101)-OnlyPretone_correct_congruent_ChannelCollapsed_AC_std(901:1101)), 'blue');
plot(OnlyPretone_correct_congruent_ChannelCollapsed_AC(901:1101), 'color', [0 0 0]);
title('PretoneOnly Correct')
xlim([1 200])
ylabel('mean LFP (uV)')
xlabel('time (ms)')
hold off

%% mean LFP wrong trials 

figure()
hold on
ciplot((OnlyPrior_wrong_congruent_ChannelCollapsed_PFC(901:1101)+OnlyPrior_wrong_congruent_ChannelCollapsed_PFC_std(901:1101)), (OnlyPrior_wrong_congruent_ChannelCollapsed_PFC(901:1101)-OnlyPrior_wrong_congruent_ChannelCollapsed_PFC_std(901:1101)), 'green');
plot(OnlyPrior_wrong_congruent_ChannelCollapsed_PFC(901:1101), 'color', [0 0 0]);
xlim([1 200])

ciplot((OnlyPrior_wrong_congruent_ChannelCollapsed_AC(901:1101)+OnlyPrior_wrong_congruent_ChannelCollapsed_AC_std(901:1101)), (OnlyPrior_wrong_congruent_ChannelCollapsed_AC(901:1101)-OnlyPrior_wrong_congruent_ChannelCollapsed_AC_std(901:1101)), 'blue');
plot(OnlyPrior_wrong_congruent_ChannelCollapsed_AC(901:1101), 'color', [0 0 0]);
title('PriorOnly Wrong')
xlim([1 200])
ylabel('mean LFP (uV)')
xlabel('time (ms)')
hold off

%% mean LFP wrong trials 

figure()
hold on
ciplot((OnlyPretone_wrong_congruent_ChannelCollapsed_PFC(901:1101)+OnlyPretone_wrong_congruent_ChannelCollapsed_PFC_std(901:1101)), (OnlyPretone_wrong_congruent_ChannelCollapsed_PFC(901:1101)-OnlyPretone_wrong_congruent_ChannelCollapsed_PFC_std(901:1101)), 'green');
plot(OnlyPretone_wrong_congruent_ChannelCollapsed_PFC(901:1101), 'color', [0 0 0]);
xlim([1 200])

ciplot((OnlyPretone_wrong_congruent_ChannelCollapsed_AC(901:1101)+OnlyPretone_wrong_congruent_ChannelCollapsed_AC_std(901:1101)), (OnlyPretone_wrong_congruent_ChannelCollapsed_AC(901:1101)-OnlyPretone_wrong_congruent_ChannelCollapsed_AC_std(901:1101)), 'blue');
plot(OnlyPretone_wrong_congruent_ChannelCollapsed_AC(901:1101), 'color', [0 0 0]);
title('PretoneOnly Wrong')
xlim([1 200])
ylabel('mean LFP (uV)')
xlabel('time (ms)')
hold off

%% Calculate Coherence

% Define the indices for the 0-200 ms epoch
epoch_start = 901;      % 0 ms corresponds to index 901
epoch_end = 1101;       % 200 ms corresponds to index 1101

% Extract the 200 ms epoch from x and y
OnlyPrior_correct_congruent_ChannelCollapsed_PFC_epoch = OnlyPrior_correct_congruent_ChannelCollapsed_PFC(epoch_start:epoch_end); 
OnlyPrior_correct_congruent_ChannelCollapsed_AC_epoch  = OnlyPrior_correct_congruent_ChannelCollapsed_AC(epoch_start:epoch_end);

window_epoch = hanning(100); % 200 ms window (no need for 1.6 s)
noverlap_epoch = 50; % 50% overlap for 200 ms
nfft_epoch = 1000; % Smaller nfft for higher resolution in the epoch

% Compute coherence on the 200 ms epoch
[Cxy, F] = mscohere(OnlyPrior_correct_congruent_ChannelCollapsed_PFC_epoch, OnlyPrior_correct_congruent_ChannelCollapsed_AC_epoch, window_epoch, noverlap_epoch, nfft_epoch, 1000);

% Extract coherence values for the beta band (13-30 Hz)
freq_band = (F >= 1) & (F <= 30);
OnlyPrior_correct_congruent_ChannelCollapsed_coherence  = Cxy(freq_band);
OnlyPrior_correct_congruent_ChannelCollapsed_freqs  = F(freq_band);


%% Calculate Coherence

% Define the indices for the 0-200 ms epoch
epoch_start = 901;      % 0 ms corresponds to index 901
epoch_end = 1101;       % 200 ms corresponds to index 1101

% Extract the 200 ms epoch from x and y
OnlyPretone_correct_congruent_ChannelCollapsed_PFC_epoch = OnlyPretone_correct_congruent_ChannelCollapsed_PFC(epoch_start:epoch_end); 
OnlyPretone_correct_congruent_ChannelCollapsed_AC_epoch  = OnlyPretone_correct_congruent_ChannelCollapsed_AC(epoch_start:epoch_end);

window_epoch = hanning(100); % 200 ms window (no need for 1.6 s)
noverlap_epoch = 50; % 50% overlap for 200 ms
nfft_epoch = 1000; % Smaller nfft for higher resolution in the epoch

% Compute coherence on the 200 ms epoch
[Cxy, F] = mscohere(OnlyPretone_correct_congruent_ChannelCollapsed_PFC_epoch, OnlyPretone_correct_congruent_ChannelCollapsed_AC_epoch, window_epoch, noverlap_epoch, nfft_epoch, 1000);

% Extract coherence values for the beta band (13-30 Hz)
freq_band = (F >= 1) & (F <= 30);
OnlyPretone_correct_congruent_ChannelCollapsed_coherence  = Cxy(freq_band);
OnlyPretone_correct_congruent_ChannelCollapsed_freqs  = F(freq_band);

%%

% Extract the 200 ms epoch from x and y
OnlyPrior_wrong_congruent_ChannelCollapsed_PFC_epoch = OnlyPrior_wrong_congruent_ChannelCollapsed_PFC(epoch_start:epoch_end); 
OnlyPrior_wrong_congruent_ChannelCollapsed_AC_epoch  = OnlyPrior_wrong_congruent_ChannelCollapsed_AC(epoch_start:epoch_end);

window_epoch = hanning(100); % 200 ms window (no need for 1.6 s)
noverlap_epoch = 50; % 50% overlap for 200 ms
nfft_epoch = 1000; % Smaller nfft for higher resolution in the epoch

% Compute coherence on the 200 ms epoch
[Cxy, F] = mscohere(OnlyPrior_wrong_congruent_ChannelCollapsed_PFC_epoch, OnlyPrior_wrong_congruent_ChannelCollapsed_AC_epoch, window_epoch, noverlap_epoch, nfft_epoch, 1000);

% Extract coherence values for the beta band (13-30 Hz)
freq_band = (F >= 1) & (F <= 30);
OnlyPrior_wrong_congruent_ChannelCollapsed_coherence = Cxy(freq_band);
OnlyPrior_wrong_congruent_ChannelCollapsed_freqs = F(freq_band);

%%

% Extract the 200 ms epoch from x and y
OnlyPretone_wrong_congruent_ChannelCollapsed_PFC_epoch = OnlyPretone_wrong_congruent_ChannelCollapsed_PFC(epoch_start:epoch_end); 
OnlyPretone_wrong_congruent_ChannelCollapsed_AC_epoch  = OnlyPretone_wrong_congruent_ChannelCollapsed_AC(epoch_start:epoch_end);

window_epoch = hanning(100); % 200 ms window (no need for 1.6 s)
noverlap_epoch = 50; % 50% overlap for 200 ms
nfft_epoch = 1000; % Smaller nfft for higher resolution in the epoch

% Compute coherence on the 200 ms epoch
[Cxy, F] = mscohere(OnlyPretone_wrong_congruent_ChannelCollapsed_PFC_epoch, OnlyPretone_wrong_congruent_ChannelCollapsed_AC_epoch, window_epoch, noverlap_epoch, nfft_epoch, 1000);

% Extract coherence values for the beta band (13-30 Hz)
freq_band = (F >= 1) & (F <= 30);
OnlyPretone_wrong_congruent_ChannelCollapsed_coherence = Cxy(freq_band);
OnlyPretone_wrong_congruent_ChannelCollapsed_freqs = F(freq_band);

%%
hold on
plot(OnlyPrior_correct_congruent_ChannelCollapsed_coherence, 'color', 'b', 'LineWidth', 2)
plot(OnlyPrior_wrong_congruent_ChannelCollapsed_coherence, 'color', 'r','LineWidth', 2)
xlabel('Frequency Hz')
ylabel('Coherence Value')
legend({'correct', 'wrong'})
title('PriorOnly Correct versus Wrong Coherence')
hold off

%%
hold on
plot(OnlyPretone_correct_congruent_ChannelCollapsed_coherence, 'color', 'b', 'LineWidth', 2)
plot(OnlyPretone_wrong_congruent_ChannelCollapsed_coherence, 'color', 'r','LineWidth', 2)
xlabel('Frequency Hz')
ylabel('Coherence Value')
legend({'correct', 'wrong'})
title('OnlyPretone Correct versus Wrong Coherence')
hold off



%%
hold on
plot(OnlyPrior_correct_congruent_ChannelCollapsed_coherence, 'color', 'g', 'LineWidth', 2)
plot(OnlyPretone_correct_congruent_ChannelCollapsed_coherence, 'color', 'y','LineWidth', 2)
xlabel('Frequency Hz')
ylabel('Coherence Value')
legend({'OnlyPrior', 'OnlyPretone'})
title('OnlyPrior Correct versus OnlyPretone')
hold off

%% XCorrelation 
 [c_OnlyPrior_correct, lags_OnlyPrior_correct] = xcorr(OnlyPrior_correct_congruent_ChannelCollapsed_PFC_epoch, OnlyPrior_correct_congruent_ChannelCollapsed_AC_epoch, 'coeff');
 [c_OnlyPrior_wrong, lags_OnlyPrior_wrong] = xcorr(OnlyPrior_wrong_congruent_ChannelCollapsed_PFC_epoch, OnlyPrior_wrong_congruent_ChannelCollapsed_AC_epoch, 'coeff');


 %% XCorrelation 
 [c_OnlyPretone_correct, lags_OnlyPretone_correct] = xcorr(OnlyPretone_correct_congruent_ChannelCollapsed_PFC_epoch, OnlyPretone_correct_congruent_ChannelCollapsed_AC_epoch, 'coeff');
 [c_OnlyPretone_wrong, lags_OnlyPretone_wrong] = xcorr(OnlyPretone_wrong_congruent_ChannelCollapsed_PFC_epoch, OnlyPretone_wrong_congruent_ChannelCollapsed_AC_epoch, 'coeff');

 %%

figure()
hold on 
plot(lags_correct, c_correct, 'color', 'b', 'LineWidth', 2);
plot(lags_wrong, c_wrong, 'color', 'r', 'LineWidth', 2);
xlabel('Lags (ms)')
ylabel('Cross Correlation Value')
legend({'correct', 'wrong'})
hold off

%%

figure()
hold on 
plot(lags_OnlyPrior_correct, c_OnlyPrior_correct, 'color', 'g', 'LineWidth', 2);
plot(lags_OnlyPretone_correct, c_OnlyPretone_correct, 'color', 'y', 'LineWidth', 2);
xlabel('Lags (ms)')
ylabel('Cross Correlation Value')
legend({'OnlyPrior', 'OnlyPretone'})
hold off

%%

figure()
hold on 
plot(lags_OnlyPrior_wrong, c_OnlyPrior_wrong, 'color', 'g', 'LineWidth', 2);
plot(lags_OnlyPretone_wrong, c_OnlyPretone_wrong, 'color', 'y', 'LineWidth', 2);
xlabel('Lags (ms)')
ylabel('Cross Correlation Value')
legend({'OnlyPrior', 'OnlyPretone'})
hold off
%% Build Data Array for PSI
xy_OnlyPrior_correct = [OnlyPrior_correct_congruent_ChannelCollapsed_PFC_epoch', OnlyPrior_correct_congruent_ChannelCollapsed_AC_epoch'];
segleng = 500;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [4:8]; % Theta

[psi_OnlyPrior_theta, stdpsi_OnlyPrior_theta, ~, ~] = data2psi(xy_OnlyPrior_correct, segleng, epleng, freqbins);

segleng = 500;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [9:14]; % Theta

[psi_OnlyPrior_alpha, stdpsi_OnlyPrior_alpha, ~, ~] = data2psi(xy_OnlyPrior_correct, segleng, epleng, freqbins);

segleng = 500;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [15:30]; % Theta

[psi_OnlyPrior_beta, stdpsi_OnlyPrior_beta, ~, ~] = data2psi(xy_OnlyPrior_correct, segleng, epleng, freqbins);

%% Build Data Array for PSI
xy_OnlyPretone_correct = [OnlyPretone_correct_congruent_ChannelCollapsed_PFC_epoch', OnlyPretone_correct_congruent_ChannelCollapsed_AC_epoch'];
segleng = 500;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [4:8]; % Theta

[psi_OnlyPretone_theta, stdpsi_OnlyPretone_theta, ~, ~] = data2psi(xy_OnlyPretone_correct, segleng, epleng, freqbins);

segleng = 500;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [9:14]; % Theta

[psi_OnlyPretone_alpha, stdpsi_OnlyPrior_alpha, ~, ~] = data2psi(xy_OnlyPretone_correct, segleng, epleng, freqbins);

segleng = 200;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [15:30]; % Theta

[psi_OnlyPretone_beta, stdpsi_OnlyPretone_beta, ~, ~] = data2psi(xy_OnlyPretone_correct, segleng, epleng, freqbins);


%% Build Data Array for PSI
xy_wrong = [OnlyPrior_wrong_congruent_ChannelCollapsed_PFC_epoch', OnlyPrior_wrong_congruent_ChannelCollapsed_AC_epoch'];
segleng = 500;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [4:8]; % Theta

[psi_wrong_theta, stdpsi_wrong_theta, ~, ~] = data2psi(xy_wrong, segleng, epleng, freqbins);

segleng = 500;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [9:14]; % Theta

[psi_wrong_alpha, stdpsi_wrong_alpha, ~, ~] = data2psi(xy_wrong, segleng, epleng, freqbins);

segleng = 500;  % 1 Hz resolution
epleng = 20;     % Skip standard deviation calculation
freqbins = [15:30]; % Theta

[psi_wrong_beta, stdpsi_wrong_beta, ~, ~] = data2psi(xy_wrong, segleng, epleng, freqbins);