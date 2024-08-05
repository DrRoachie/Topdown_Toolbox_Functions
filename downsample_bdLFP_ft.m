function [data_bdlfp] = downsample_bdLFP_ft(lfp,t,Fold,Fnew,rmNoise,eInfo)
%downsample_bdLFP Summary of this function goes here
%   lfp     ... LFP data (channel x sample x trial)
%   t       ... time vector corresponding to the LFP
%   Fnew    ... target frequency usually set 1000 Hz
%   rmNoise ... remove noise ('Y') or not ('N')
%   data_bdlpf ... bipolar derived lfp

% for noise reduction
% params.tapers = [3 5]; %[3 5];
% params.Fs  = Fnew;
% params.pad = 3;


% % subtract trial average to get induced activity
% disp('subtracting trial average...')
% mat_mlfp = repmat(mean(lfp,3),1,1,size(lfp,3)); % matrix of trial averaged LFP
% lfp = lfp - mat_mlfp;
% clear mat_mlfp

nTrial = size(lfp,3);
for n=1:nTrial
    clc;
    disp(['processing ' num2str(n) ' out of ' num2str(nTrial) ' trials...'])
    % downsampling
    Y = lfp(:,:,n); % n-th trial (channel x sample)
    T = transpose(t/1000); % time in sec
    
    % low-pass filter before resample (cutoff freq = 250 Hz)
    Y = ft_preproc_lowpassfilter(Y,Fold,250,4);
    Y = Y'; % reshape matrix (sample x channel)
    % resample data
    [rs_lfp,tt] = resample(Y,T,Fnew,'spline');
    rs_lfp = rs_lfp'; tt = tt';
    
    % baseline correction
%     bc_lfp = basecorrectLFP(rs_lfp,tt,[-0.7 -0.65]);
    bc_lfp = ft_preproc_baselinecorrect(rs_lfp,51,100); % [-0.85 -0.80]
    
    % bipolar derivation (re-referencing)
    [bipolar_lfp,label] = bipolarLFP(bc_lfp,eInfo);
    
    % noise removal
    if strcmp(rmNoise,'Y')        
        y = bipolar_lfp;
        % low-pass filter (cutoff = 250 Hz)
%        y = ft_preproc_lowpassfilter(y,Fnew,250,4); % 4th order
        y = ft_preproc_bandpassfilter(y,Fnew,[1 250],4);
        
        % band-stop filter (remove line noise)
        y = ft_preproc_bandstopfilter(y,Fnew,[55 65]); % 4th order filter
        y = ft_preproc_bandstopfilter(y,Fnew,[119 121]);
        y = ft_preproc_bandstopfilter(y,Fnew,[179 181]);
        y = ft_preproc_bandstopfilter(y,Fnew,[239 241]);
%         y = rmlinesc(bipolar_lfp,params,[],[],60);
%         y = rmlinesc(y,params,[],[],120);
%         y = rmlinesc(y,params,[],[],180);
%         y = rmlinesc(y,params,[],[],240);
        
        bipolar_lfp_rmNoise = y;
    elseif strcmp(rmNoise,'N')
        % low-pass filter (cutoff = 250 Hz)
        y = bipolar_lfp;
        y = ft_preproc_lowpassfilter(y,Fnew,250,4); % 4th order
        bipolar_lfp_rmNoise = y;
    end
    
    % put the data into cell
    trial{n} = bipolar_lfp_rmNoise; %transpose(bipolar_lfp_rmNoise);
    time{n} = tt;
end

% insert random number in sample info...
nSample = size(rs_lfp,2); % samples in resampled data
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
data_bdlfp.fsample = Fnew;
data_bdlfp.sampleinfo = sinfo;

end

