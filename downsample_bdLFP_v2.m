function [data_bdlfp] = downsample_bdLFP_v2(LFP,t,dsFreq,rmNoise,eInfo)
%downsample_bdLFP Summary of this function goes here
%   LFP     ... LFP data (trial x sample x channel)
%   t       ... time vector corresponding to the LFP
%   dsFreq    ... target frequency usually set 1000 Hz
%   rmNoise ... remove noise ('Y') or not ('N')
%   data_bdlpf ... bipolar derived lfp

lfp = permute(LFP,[3 2 1]); % LFP (channel x sample x trial)
%Fs_original = 24414.0625/8; % sampling frequency of original
% for noise reduction
params.tapers = [3 5]; %[3 5];
params.Fs  = dsFreq;
params.pad = 3;

% lfp = lfp * 10e3; % convert unit to mV
lfp = lfp * 10e6; % convert unit to uV


nTrial = size(lfp,3);

for n=1:nTrial
    clc;
    disp(['processing ' num2str(n) ' out of ' num2str(nTrial) ' trials...'])
    % downsampling
    Y = transpose(lfp(:,:,n)); % n-th trial (sample x channel)
    T = transpose(t/1000); % time in sec
    
    % resample data
    [rs_lfp,tt] = resample(Y,T,dsFreq,'spline');
    
    % baseline correction
    bc_lfp = basecorrectLFP(rs_lfp,tt,[-0.7 -0.65]);
    
    % bipolar derivation (re-referencing)
    [bd_lfp,label] = bipolarLFP(bc_lfp,eInfo);
    
    % noise removal (require chronux)
    if strcmp(rmNoise,'Y')        
        % remove noise (60, 120, 180, and 240 Hz)
        % y = rmlinesc(bd_lfp,params,[],[], 60);
        % y = rmlinesc(y,params,[],[],120);
        % y = rmlinesc(y,params,[],[],180);
        % y = rmlinesc(y,params,[],[],240);

        y = rmlinesc(bd_lfp,params,[],[], [50 70]);
        y = rmlinesc(y,params,[],[],[110 130]);
        y = rmlinesc(y,params,[],[],[170 190]);
        y = rmlinesc(y,params,[],[],[230 250]);
        
        bd_lfp_rmNoise = y;
    elseif strcmp(rmNoise,'N')
        % low-pass filter (cutoff = 150 Hz)
%         y = filtfilt(b_lp,a_lp,bd_lfp);
        bd_lfp_rmNoise = y;
    end
    
    % put the data into cell
    trial{n} = transpose(bd_lfp_rmNoise);
    time{n} = transpose(tt);
end

nSample = size(rs_lfp,1); % samples in resampled data
s_start = 1;
for j=1:nTrial
    s_end = s_start + nSample -1;
    sinfo(j,:) = [s_start s_end];
    % update s_start
    s_start = s_end + 5; % make sure no overlap between sample periods
end

% data_lfp.hdr = hdr;
data_bdlfp.label = label;
data_bdlfp.time = time;
data_bdlfp.trial = trial;
data_bdlfp.fsample = dsFreq;
data_bdlfp.sampleinfo = sinfo;

end

