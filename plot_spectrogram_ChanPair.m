function [] = plot_spectrogram_ChanPair(tfreq, Current_ChanPair, Frequency_Band)

powspectrm  = tfreq.powspctrm;
t           = tfreq.time(1:end);

if strcmp(Frequency_Band,'theta')
f           = tfreq.freq(5:8);
label       = tfreq.label;
end

if strcmp(Frequency_Band,'alpha')
f           = tfreq.freq(8:14);
label       = tfreq.label;
end

if strcmp(Frequency_Band,'beta')
f           = tfreq.freq(15:30);
label       = tfreq.label;
end

if strcmp(Frequency_Band,'gamma')
f           = tfreq.freq(31:54);
label       = tfreq.label;
end

if strcmp(Frequency_Band,'highGamma')
f           = tfreq.freq(66:100);
label       = tfreq.label;
end

% Remove the prefix from the channel labels and add an asterisk at the start
modified_labels = cellfun(@(x) ['*' x(4:end)], label, 'UniformOutput', false);

eID_pfc = Current_ChanPair{1, 1};
eID_ac  = Current_ChanPair{1, 2};

% Find the index of the PFC channel
pfc_index = find(strcmp(modified_labels , eID_pfc));

% Find the index of the AC channel
ac_index = find(strcmp(modified_labels , eID_ac ));


if strcmp(Frequency_Band,'theta')
figure 
subplot(1, 2, 1)
imagesc(t,f,squeeze(powspectrm(pfc_index, 5:8,:)));
    set(gca,'YDir','normal');
    title(eID_pfc,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
   
subplot(1, 2, 2) 
    imagesc(t,f,squeeze(powspectrm(ac_index,5:8,:)));
    set(gca,'YDir','normal');
    title(eID_ac,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
 
end

if strcmp(Frequency_Band,'alpha')
figure 
subplot(1, 2, 1)
imagesc(t,f,squeeze(powspectrm(pfc_index, 8:14,:)));
    set(gca,'YDir','normal');
    title(eID_pfc,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
   
subplot(1, 2, 2) 
    imagesc(t,f,squeeze(powspectrm(ac_index, 8:14,:)));
    set(gca,'YDir','normal');
    title(eID_ac,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
 
end


if strcmp(Frequency_Band,'beta')
figure 
subplot(1, 2, 1)
imagesc(t,f,squeeze(powspectrm(pfc_index, 15:30,:)));
    set(gca,'YDir','normal');
    title(eID_pfc,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
   
subplot(1, 2, 2) 
    imagesc(t,f,squeeze(powspectrm(ac_index,15:30,:)));
    set(gca,'YDir','normal');
    title(eID_ac,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
 
end

if strcmp(Frequency_Band,'gamma')
figure 
subplot(1, 2, 1)
imagesc(t,f,squeeze(powspectrm(pfc_index, 31:54,:)));
    set(gca,'YDir','normal');
    title(eID_pfc,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
   
subplot(1, 2, 2) 
    imagesc(t,f,squeeze(powspectrm(ac_index,31:54,:)));
    set(gca,'YDir','normal');
    title(eID_ac,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
 
end

if strcmp(Frequency_Band,'highGamma')
figure 
subplot(1, 2, 1)
imagesc(t,f,squeeze(powspectrm(pfc_index, 66:100,:)));
    set(gca,'YDir','normal');
    title(eID_pfc,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
   
subplot(1, 2, 2) 
    imagesc(t,f,squeeze(powspectrm(ac_index,66:100,:)));
    set(gca,'YDir','normal');
    title(eID_ac,'Interpreter','none');
    xlabel('Time'); ylabel('Frequency');
 
end

colormap jet

end




