function [] = plot_coherence(f,Coherence,string,isLabel)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n_ch = size(Coherence,1); % number of channel
c = parula(n_ch); % specify color map to use
% c = jet(n_ch); % specify color map to use

for i=1:n_ch
%     plot(f,Coherence(i,:),'Color',c(i,:)); hold on;
      plot(f,Coherence(i,1:101),'Color',c(i,:)); hold on;
end
xlabel('Frequency [Hz]'); ylabel('Coherence');

if strcmp(isLabel,'y')||strcmp(isLabel,'Y')||isLabel==1
if n_ch==14 % 16-ch electrode
    ch_label = {'ch02','ch03','ch04','ch05','ch06','ch07','ch08','ch09','ch10','ch11','ch12','ch13','ch14','ch15'};
elseif n_ch==20 % 24-ch electrode
    ch_label = {'ch03','ch04','ch05','ch06','ch07','ch08','ch09','ch10','ch11','ch12','ch13','ch14','ch15','ch16','ch17','ch18','ch19','ch20','ch21','ch22'};
end
legend(ch_label);
end

s = inputname(2);
if strcmp(s(end),'c')
    s_title = [string ' correct'];
elseif strcmp(s(end),'w')
    s_title = [string ' wrong'];
end
title(s_title);

box off;

end

