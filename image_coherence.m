function [] = image_coherence(f,Coherence,string)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n_ch = size(Coherence,1); % number of channel
if n_ch==14 % 16-ch electrode
    list_ch = 2:15;
elseif n_ch==20 % 24-ch electrode
    list_ch = 3:22;
end

imagesc(f(f<150),list_ch,Coherence(:,f<150));
%imagesc(f(f<200),list_ch,Coherence(:,f<200));
xlabel('Frequency [Hz]'); ylabel('Electrode Channel');
colormap jet

s = inputname(2);
if strcmp(s(end),'c')
    s_title = [string ' correct'];
elseif strcmp(s(end),'w')
    s_title = [string ' wrong'];
elseif strcmp(s(end),'d')
    s_title = [string 'difference'];
end
title(s_title);


end

