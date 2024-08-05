function [descriptives_onlyprior_coh,descriptives_onlypretone_coh] = GetCohDescriptives(PFC_AC_Coh_Spectrum_PriorOnly,PFC_AC_Coh_Spectrum_PretoneOnly)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% calculate standard error of OnlyPrior Condition using leave-one-out
% method

num_trials = PFC_AC_Coh_Spectrum_PriorOnly.dof;     % dof field is equal to the number of trials



            % calculate the descriptive of the OnlyPrior Condition for correct trials 
             descriptives_onlyprior_coh.avg           = squeeze(mean(PFC_AC_Coh_Spectrum_PriorOnly.cohspctrm,1));
             descriptives_onlyprior_coh.stdErr        = PFC_AC_Coh_Spectrum_PriorOnly.std';
             descriptives_onlyprior_coh.CIupper       = descriptives_onlyprior_coh.avg  + descriptives_onlyprior_coh.stdErr;      
             descriptives_onlyprior_coh.CIlower       = descriptives_onlyprior_coh.avg  - descriptives_onlyprior_coh.stdErr;
            
             % calculate the descriptives of the OnlyPretone Condition for correct trials 
             descriptives_onlypretone_coh.avg           = squeeze(mean(PFC_AC_Coh_Spectrum_PretoneOnly.cohspctrm,1));
             descriptives_onlypretone_coh.stdErr        = PFC_AC_Coh_Spectrum_PretoneOnly.std';
             descriptives_onlypretone_coh.CIupper       = descriptives_onlypretone_coh.avg  + descriptives_onlypretone_coh.stdErr;      
             descriptives_onlypretone_coh.CIlower       = descriptives_onlypretone_coh.avg  - descriptives_onlypretone_coh.stdErr;
end