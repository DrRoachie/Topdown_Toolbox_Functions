function [selectData] = selectData(data,params,iSelect)
%UNTITLED Summary of this function goes here
%    select data based on given parameter value
%cat

N = numel(data.time); % total number of trial
INDEX = zeros(N,1);


% Check for 'choice' condition (if exists)
if ~isempty(iSelect.choice)
    i = params.choice == iSelect.choice;
else
    i = ones(N,1);
end
INDEX = INDEX + i;

% Check for 'err' condition (valid or specific error)
if ~isempty(iSelect.err)
    if strcmp(iSelect.err,'valid')
        i = ( params.err == 'w' | params.err == 'c' );
    else
        i = params.err == iSelect.err;
    end
else
    i = ones(N,1);
end
INDEX = INDEX + i;

% Check for 'pretone' condition
if ~isempty(iSelect.pretone)
    i = params.pretone == iSelect.pretone;
else
    i = ones(N,1);
end
INDEX = INDEX + i;

% Check for 'pretoneLength' condition
if ~isempty(iSelect.pretoneLength)
    i = params.pretoneLength == iSelect.pretoneLength;
else
    i = ones(N,1);
end
INDEX = INDEX + i;

% Check for 'prior' condition
if ~isempty(iSelect.prior)
    if strcmp(iSelect.prior,'X') % trials in which pretone was presented...
        i = params.prior ~= 'N';
    else
        i = params.prior == iSelect.prior;
    end
else
    i = ones(N,1);
end
INDEX = INDEX + i;

% Check for 'SNR' condition 
if ~isempty(iSelect.SNR)
    if numel(iSelect.SNR) > 1  % If multiple values are provided
        i = ismember(params.SNR, iSelect.SNR);  % Match any of the values in iSelect.SNR
    else
        i = params.SNR == iSelect.SNR;  % Match the single value
    end
else
    i = ones(N, 1);  % Default to all trials if no SNR condition
end
INDEX = INDEX + i;

% Check for 'congruency' condition
if ~isempty(iSelect.congruency)
    if strcmp(iSelect.congruency, 'congruent')
        i = strcmp(params.congruency, 'congruent');
    elseif strcmp(iSelect.congruency, 'incongruent')
        i = strcmp(params.congruency, 'incongruent');
    elseif strcmp(iSelect.congruency, 'neutral')
        i = strcmp(params.congruency, 'neutral');
    else
        i = ones(N, 1);  % Default to all trials if no match
    end
else
    i = ones(N, 1);
end

INDEX = INDEX + i;

% Index to select trials where all conditions are met
iSelectTrial = INDEX == numel(fieldnames(iSelect));  % Matches total number of conditions

% Select trials based on the index
temp_time = data.time(iSelectTrial);
temp_trial = data.trial(iSelectTrial);
temp_sampleinfo = data.sampleinfo(iSelectTrial, :);

% Return selected data
selectData.label = data.label;
selectData.time = temp_time;
selectData.trial = temp_trial;
selectData.fsample = data.fsample;
selectData.sampleinfo = temp_sampleinfo;


end

