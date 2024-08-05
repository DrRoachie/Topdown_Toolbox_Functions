function [selectData] = selectData(data,params,iSelect)
%UNTITLED Summary of this function goes here
%    select data based on given parameter value

N = numel(data.time); % total number of trial
INDEX = zeros(N,1);

if ~isempty(iSelect.choice)
    i = params.choice == iSelect.choice;
else
    i = ones(N,1);
end
INDEX = INDEX + i;

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

if ~isempty(iSelect.pretone)
    i = params.pretone == iSelect.pretone;
else
    i = ones(N,1);
end
INDEX = INDEX + i;

if ~isempty(iSelect.pretoneLength)
    i = params.pretoneLength == iSelect.pretoneLength;
else
    i = ones(N,1);
end
INDEX = INDEX + i;

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

if ~isempty(iSelect.SNR)
    i = params.SNR == iSelect.SNR;
else
    i = ones(N,1);
end
INDEX = INDEX + i;

% index
iSelectTrial = INDEX==6;

% select trials
temp_time = data.time(iSelectTrial);
temp_trial = data.trial(iSelectTrial);
temp_sampleinfo = data.sampleinfo(iSelectTrial,:);

% selected data
selectData.label = data.label;
selectData.time = temp_time;
selectData.trial = temp_trial;
selectData.fsample = data.fsample;
selectData.sampleinfo = temp_sampleinfo;

end

