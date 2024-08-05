
function [shuffled_xcorr_result, shuffled_lagResults, xcorr_result, lagResults] = TrialbyTrialXcorr(datadir, Animal, RecDate, Epoch, Condition, Current_ChanPair)
% this function calculates the trial by trial xcorr and lag. 

% When you compute xcorr(x, y), positive lags indicate that x is leading y, and negative lags indicate that y is leading x.

% load fieldtrip data 
addpath(genpath(datadir));
fName = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch,'_ft');
load(fName);

% set parameters for data extractions
params.choice = choice;
params.err = err;
params.pretone = pretone;
params.pretoneLength = pretoneLength;
params.prior = cell2char(prior); %cell2char is an auxillary function that needs to be in your path
params.SNR = SNR;

% Choose the data that you want to analyze based on parameters 

iSelect = setStimulusCondition(Condition);

iSelect.err = 'c'; % choose correct trials
data_c = selectData(data,params,iSelect);
label  = data.label;

numTrial_c = length(data_c.trial);

% Remove the prefix from the channel labels and add an asterisk at the start
modified_labels = cellfun(@(x) ['*' x(4:end)], label, 'UniformOutput', false);

eID_pfc = Current_ChanPair{1, 1};
eID_ac  = Current_ChanPair{1, 2};

% Find the index of the PFC channel
pfc_index = find(strcmp(modified_labels , eID_pfc));

% Find the index of the AC channel
ac_index = find(strcmp(modified_labels , eID_ac ));


%% Loop through each trial and grab the row that corresponds to each channel in the PFC array

% Preallocate the matrix to store the extracted rows
numCells = numel(data_c.trial);
numTimePoints = size(data_c.trial{1}, 2); % Assuming each trial has the same number of time points
PFC_extractedRows = zeros(numCells, numTimePoints);

% Loop through each cell and extract the relevant row
for i = 1:numCells
    % Extract the relevant row based on pfc_index
    PFC_extractedRows(i, :) = data_c.trial{i}(pfc_index, :);
end

%% Loop through each trial and grab the row that corresponds to each channel in the AC array

% Preallocate the matrix to store the extracted rows
numCells = numel(data_c.trial);
numTimePoints = size(data_c.trial{1}, 2); % Assuming each trial has the same number of time points
AC_extractedRows = zeros(numCells, numTimePoints);

% Loop through each cell and extract the relevant row
for ii = 1:numCells
    % Extract the relevant row based on pfc_index
   AC_extractedRows(ii, :) = data_c.trial{ii}(ac_index, :);
end

%% Initialize the cross-correlation trial-by-trial structs

% Assuming PFC_extractedrows and AC_extractedrows are your two arrays
numRows = size(PFC_extractedRows, 1);

% Calculate the length of the cross-correlation result
[c, lags] = xcorr(PFC_extractedRows(1, :), AC_extractedRows(1, :));
numLags = length(c);

% Preallocate arrays to store the cross-correlation results and lags
xcorrResults = zeros(numRows, numLags);
lagsResults = zeros(1, numLags);

% Loop through each row and calculate the cross-correlation
for i = 1:numRows
    % Extract the corresponding rows
    rowPFC = PFC_extractedRows(i, :);
    rowAC = AC_extractedRows(i, :);
    
    % Calculate the cross-correlation
    [c, lags] = xcorr(rowPFC, rowAC, 'coeff');  % The 'coeff' option in xcorr scales the cross-correlation so that the autocorrelation at zero lag 
                                                % (i.e., the cross-correlation of the signal with itself) is equal to 1. 
                                                % This is achieved by dividing the cross-correlation values by the product of 
                                                % the standard deviations of the two signals.
    

    % Store the result in the preallocated arrays
    xcorrResults(i, :) = c;
    
    % Store the lags only once since they are the same for each row
    if i == 1
        lagResults(:) = lags;
    end
end

%% averaging the trial-by-trial cross core

xcorr_result.avg = mean(xcorrResults, 1);
xcorr_result.sem = std(xcorrResults, 0, 1) / sqrt(numTrial_c);

% Calculate 95% confidence intervals
ci_xcorr = 1.96 * xcorr_result.sem;

% Extract lower and upper bounds for xcorr
xcorr_result.lower = xcorr_result.avg - ci_xcorr;
xcorr_result.upper = xcorr_result.avg + ci_xcorr;

%% Trial Shuffle the Extracted Rows to get null 

AC_numCols = size(AC_extractedRows, 2);
AC_shuffledIndices = randperm(AC_numCols);
AC_extractedRows_shuffled = AC_extractedRows(:, AC_shuffledIndices);

PFC_numCols = size(PFC_extractedRows, 2);
PFC_shuffledIndices = randperm(PFC_numCols);
PFC_extractedRows_shuffled = PFC_extractedRows(:, PFC_shuffledIndices);

%% Carry out the cross-correlation on the shuffled data

% Assuming PFC_extractedrows and AC_extractedrows are your two arrays
shuffled_numRows = size(PFC_extractedRows_shuffled, 1);

% Calculate the length of the cross-correlation result
[shuffled_c, shuffled_lags] = xcorr(PFC_extractedRows_shuffled(1, :), AC_extractedRows_shuffled(1, :));
shuffled_numLags = length(shuffled_c);

% Preallocate arrays to store the cross-correlation results and lags
shuffled_xcorrResults = zeros(shuffled_numRows, shuffled_numLags);
shuffled_lagsResults = zeros(1, shuffled_numLags);

% Loop through each row and calculate the cross-correlation
for i = 1:numRows
    % Extract the corresponding rows
    rowPFC_shuffled = PFC_extractedRows_shuffled(i, :);
    rowAC_shuffled = AC_extractedRows_shuffled(i, :);
    
    % Calculate the cross-correlation
    [shuffled_c, shuffled_lags] = xcorr(rowPFC_shuffled, rowAC_shuffled, 'coeff');  % The 'coeff' option in xcorr scales the cross-correlation so that the autocorrelation at zero lag 
                                                % (i.e., the cross-correlation of the signal with itself) is equal to 1. 
                                                % This is achieved by dividing the cross-correlation values by the product of 
                                                % the standard deviations of the two signals.
    

    % Store the result in the preallocated arrays
    shuffled_xcorrResults(i, :) = shuffled_c;
    
    % Store the lags only once since they are the same for each row
    if i == 1
        shuffled_lagResults(:) = shuffled_lags;
    end
end

%% averaging the trial-by-trial cross core

shuffled_xcorr_result.avg = mean(shuffled_xcorrResults, 1);
shuffled_xcorr_result.sem = std(shuffled_xcorrResults, 0, 1) / sqrt(numTrial_c);

% Calculate 95% confidence intervals
shuffled_ci_xcorr = 1.96 * shuffled_xcorr_result.sem;

% Extract lower and upper bounds for xcorr
shuffled_xcorr_result.lower = shuffled_xcorr_result.avg - ci_xcorr;
shuffled_xcorr_result.upper = shuffled_xcorr_result.avg + ci_xcorr;
