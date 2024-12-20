%function [RandomPartition_1_Size_c, RandomPartition_2_Size_c, RandomPartition_1_Size_w, RandomPartition_2_Size_w, All_Trials_c, All_Trials_w] = PoolPriorPretone(datadir, sessions)
function [RandomPartition_1_Size_c, RandomPartition_2_Size_c, All_Trials_c] = PoolPriorPretone(datadir, sessions)
% This function pools the  pretone and prior conditions together and
% generates the subset sizes for the Maris MonteCarlo 

% load preprocessed OnlyPrior trials 

addpath(genpath(datadir))

Condition = 'OnlyPrior';

data_OnlyPrior_c = [];

for k = 1:length(sessions) % loops through all the .m files in the directory and pulls the data. if there is only a single .m file then you can process one session at a time

baseFileName = sessions(k).name;
fullFileName = fullfile(sessions(k).folder, baseFileName);
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

iSelect.err = 'c'; % choose correct trials
data_OnlyPrior_c = [data_OnlyPrior_c, selectData(data,params,iSelect)];

% iSelect.err = 'w'; % choose wrong trials
% data_OnlyPrior_w = [data_OnlyPrior_w, selectData(data,params,iSelect)];

end

% load preprocessed OnlyPretone trials 

Condition = 'OnlyPretone';

data_OnlyPretone_c = [];
% data_OnlyPretone_w = [];

for k = 1:length(sessions)

baseFileName = sessions(k).name;
fullFileName = fullfile(sessions(k).folder, baseFileName);
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

iSelect.err = 'c'; % choose correct trials
data_OnlyPretone_c = [data_OnlyPretone_c, selectData(data,params,iSelect)];

% iSelect.err = 'w'; % choose wrong trials
% data_OnlyPretone_w = [data_OnlyPretone_w, selectData(data,params,iSelect)];

end

% concatonate each session for each condition

% concatonate OnlyPrior data across all sessions 
cat_OnlyPrior_c.label       = data_OnlyPrior_c(1).label;
cat_OnlyPrior_c.time        = [data_OnlyPrior_c(:).time];
cat_OnlyPrior_c.trial       = [data_OnlyPrior_c(:).trial];
cat_OnlyPrior_c.fsample     = data_OnlyPrior_c(1).fsample;
cat_OnlyPrior_c.sampleinfo  = vertcat(data_OnlyPrior_c.sampleinfo);

% cat_OnlyPrior_w.label       = data_OnlyPrior_w(1).label;
% cat_OnlyPrior_w.time        = [data_OnlyPrior_w(:).time];
% cat_OnlyPrior_w.trial       = [data_OnlyPrior_w(:).trial];
% cat_OnlyPrior_w.fsample     = data_OnlyPrior_w(1).fsample;
% cat_OnlyPrior_w.sampleinfo  = vertcat(data_OnlyPrior_w.sampleinfo);
% 

% concatonate Pretone data across all sessions 
cat_OnlyPretone_c.label      = data_OnlyPretone_c(1).label;
cat_OnlyPretone_c.time       = [data_OnlyPretone_c(:).time];
cat_OnlyPretone_c.trial      = [data_OnlyPretone_c(:).trial];
cat_OnlyPretone_c.fsample    = data_OnlyPretone_c(1).fsample;
cat_OnlyPretone_c.sampleinfo = vertcat(data_OnlyPretone_c.sampleinfo);

% cat_OnlyPretone_w.label      = data_OnlyPretone_w(1).label;
% cat_OnlyPretone_w.time       = [data_OnlyPretone_w(:).time];
% cat_OnlyPretone_w.trial      = [data_OnlyPretone_w(:).trial];
% cat_OnlyPretone_w.fsample    = data_OnlyPretone_w(1).fsample;
% cat_OnlyPretone_w.sampleinfo = vertcat(data_OnlyPretone_w.sampleinfo);
% 

% Get trial number of pretone and prior conditions and set random partition size

OnlyPrior_c_TrialNum   =  length(cat_OnlyPrior_c.trial);
OnlyPretone_c_TrialNum =  length(cat_OnlyPretone_c.trial);
Total_c_TrialNum = OnlyPrior_c_TrialNum + OnlyPretone_c_TrialNum;
RandomPartition_1_Size_c = max(OnlyPretone_c_TrialNum, OnlyPrior_c_TrialNum);
RandomPartition_2_Size_c = (Total_c_TrialNum - RandomPartition_1_Size_c);

% OnlyPrior_w_TrialNum   =  length(cat_OnlyPrior_w.trial);
% OnlyPretone_w_TrialNum =  length(cat_OnlyPretone_w.trial);
% Total_w_TrialNum = OnlyPrior_w_TrialNum + OnlyPretone_w_TrialNum;
% RandomPartition_1_Size_w = max(OnlyPretone_w_TrialNum, OnlyPrior_w_TrialNum);
% RandomPartition_2_Size_w = (Total_w_TrialNum - RandomPartition_w_1_Size);

% Establish collective set featuring all the data combined into a single pool 

All_Trials_c.label      = data_OnlyPrior_c(1).label;
All_Trials_c.fsample    = data_OnlyPrior_c(1).fsample;
All_Trials_c_struct     = [cat_OnlyPretone_c, cat_OnlyPrior_c];
All_Trials_c.time       = [All_Trials_c_struct(:).time];
All_Trials_c.trial      = [All_Trials_c_struct(:).trial];
All_Trials_c.sampleinfo = vertcat(All_Trials_c_struct.sampleinfo);

% All_Trials_w.label      = data_OnlyPrior_w(1).label;
% All_Trials_w.fsample    = data_OnlyPrior_w(1).fsample;
% All_Trials_w_struct     = [cat_OnlyPretone_w, cat_OnlyPrior_w];
% All_Trials_w.time       = [All_Trials_w_struct(:).time];
% All_Trials_w.trial      = [All_Trials_w_struct(:).trial];
% All_Trials_w.sampleinfo = vertcat(All_Trials_w_struct.sampleinfo);

end