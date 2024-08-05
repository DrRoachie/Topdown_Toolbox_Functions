function [bsCorrect_lfp] = basecorrectLFP(lfp,t,basePeriod)
%basecorrectLFP Summary of this function goes here
%   function performs baseline correction of LFP
%   lfp -- lfp data (samples x channel)
%   t   -- time vector
%   basePeriod -- baseline period

bl_lfp = lfp(t>=basePeriod(1)&t<basePeriod(2),:); % cut out baseline activity
mbl_lfp = mean(bl_lfp,1);

nSample = size(lfp,1);
baseline_lfp = ones(nSample,1) * mbl_lfp;

% baseline correction
bsCorrect_lfp = lfp - baseline_lfp;
end

