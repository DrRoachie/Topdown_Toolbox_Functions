function [zPow] = zscore_spectrogram(tfreq,tfreq_baseline)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pow = tfreq.powspctrm; % power spectrum
mPow = tfreq_baseline.powspctrm_mean; % mean of baseline
sPow = tfreq_baseline.powspctrm_std; % std of baseline 

n_timepoint = size(pow,3);
mat_mPow = repmat(mPow,1,1,n_timepoint);
mat_sPow = repmat(sPow,1,1,n_timepoint);

zPow = ( pow - mat_mPow ) ./ mat_sPow;
end

