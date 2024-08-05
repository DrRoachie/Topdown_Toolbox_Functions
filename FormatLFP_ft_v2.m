function FormatLFP_ft_v2(Animal,DATA_DIR, SAVE_DIR,RecDate,Epoch)

% modify data structure for fieldtrip toolbox
% tutorial from: www.fieldtriptoolbox.org/tutorial/connectivity
% use fieldtrip functions for preprocessing (resampling, filtering, etc.)
% clear all

isSave      = 1;
isDisplay   = 0;

% set path
% DATA_DIR    = '/Volumes/SHD/MrM_DATA';
% SAVE_DIR    = '/Users/yalecohen/Library/CloudStorage/Dropbox/Takuified';


% define data to analyze
Animal2  = 'MrMiyagi'; %'MrCassius'; 
rmNoise  = 'Y'; % remove line noise 'Y' or 'N'

tic;
disp('loading data...')
fName = strcat(Animal,'-',RecDate,'_LFP_',Epoch); %taku's code for MrCassius
% fName = strcat('/',Animal,'_',RecDate,'/',Animal2,'-',RecDate,'_LFP_',Epoch); % yale's code for Miyagi
load(fullfile(DATA_DIR,fName));

% get parameters to reformat LFP data
lfp    = permute(LFP,[3 2 1]);  % LFP (channel x sample x trial)
lfp    = lfp * 10e6;            % convert unit to uV
t      = timeBin;
fs_old = info.sampFreq;      % original sampling frequency 
fs_new = 1000;               % sampling frequency after downsampling data

% get electrode information
eInfo = countChannels(chanList);

% format bipolar-derived LFP data
% Be sure LFP is (trial) x (sample) x (channel)!!
data = downsample_bdLFP_ft(lfp,t,fs_old,fs_new,rmNoise,eInfo);

stim     = testStim;
trial_id = trialID;

% save data

if isSave==1

sName = strcat(Animal,'-',RecDate,'_bdLFP_',Epoch);
    if strcmp(rmNoise,'Y')
        save_file_name = strcat(sName,'_ft');
    elseif strcmp(rmNoise,'N')
        save_file_name = strcat(sName,'_ft_noFilt');
    end
clc; disp(['saving data... ' save_file_name])

save(fullfile(SAVE_DIR,save_file_name),'data','choice','err','pretone', ...
    'pretoneLength','prior','SNR','stim','trial_id','info','-v7.3');
end

clc; disp('done!')
toc;

% display data (optional)
if isDisplay==1
    cfg = [];
    cfg.viewmode = 'vertical';
    ft_databrowser(cfg,data);
    D = cat(3,data.trial{:});
    mD = mean(D,3);
    t = data.time{1};
    figure;
    subplot(2,1,1); plot(t,mD');
    subplot(2,1,2); plot_spectrum(mD',1000);
end