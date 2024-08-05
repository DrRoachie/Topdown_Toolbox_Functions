function [iSelect] = setStimulusCondition(Condition)
%UNTITLED Summary of this function goes here
%   set iSelect values depending on Condition

if strcmp(Condition,'Both')
    % condition having both prior and pretone
    iSelect.choice = [];
    iSelect.pretone = [];
    iSelect.pretoneLength = 3;
    iSelect.prior = 'X';
    iSelect.SNR = [];
elseif strcmp(Condition,'OnlyPretone')
    % condition having pretone but no prior
    iSelect.choice = [];
    iSelect.pretone = [];
    iSelect.pretoneLength = 3; 
    iSelect.prior = 'N';
    iSelect.SNR = [];
elseif strcmp(Condition,'OnlyPrior')
    % condition having prior but no pretone
    iSelect.choice = [];
    iSelect.pretone = 'N';
    iSelect.pretoneLength = 0;
    iSelect.prior = 'X';
    iSelect.SNR = [];
end

end

