function [] = plot_spectrogram(tfreq,eID)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

zPow  = tfreq.zpowspctrm;
% t     = tfreq.time(81:102);
t     = tfreq.time;
f     = tfreq.freq;
label = tfreq.label;

index = 1:length(label);
idx  = index(contains(label,eID));
n_channel = length(idx); % number of channel in the electrode
if n_channel==0
    error(['cannot find electrode ' eID]);
end
figure('Position',[50 50 1200 700]);
for i=1:n_channel
    if n_channel == 20
        subplot(4,5,i);
    elseif n_channel == 14
        subplot(3,5,i);
    end
    % imagesc(t,f,squeeze(zPow(idx(i),:,81:102)));
    imagesc(t,f,squeeze(zPow(idx(i),:,:)));
    set(gca,'YDir','normal');
    title(label{idx(i)},'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
    caxis([-2 2]);
end
colormap jet


end

