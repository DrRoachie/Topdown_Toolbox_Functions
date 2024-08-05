function [] = Plot_CoherenceSpec(Frequency_Band, descriptives_onlyprior_coh, descriptives_onlypretone_coh)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

figure()

if strcmp(Frequency_Band,'theta') == 1 
hold on
ciplot(descriptives_onlyprior_coh.CIupper(1, 5:8) ,descriptives_onlyprior_coh.CIlower(1, 5:8), 'blue');
plot(descriptives_onlyprior_coh.avg(1, 5:8), 'color', [0 0 0]);
ciplot(descriptives_onlypretone_coh.CIupper(1, 5:8) ,descriptives_onlypretone_coh.CIlower(1, 5:8), 'green');
plot(descriptives_onlypretone_coh.avg(1, 5:8), 'color', [0 0 0]);
title ('theta')
ylabel ('coherence')
xlabel ('frequency (Hz)')
xlim([1 4])
xticks(1:4) 
xticklabels({'5';'6';'7'; '8'})
hold off
end

if strcmp(Frequency_Band,'alpha') == 1 
hold on
ciplot(descriptives_onlyprior_coh.CIupper(1, 8:14) ,descriptives_onlyprior_coh.CIlower(1, 8:14), 'blue');
plot(descriptives_onlyprior_coh.avg(1, 8:14), 'color', [0 0 0]);
ciplot(descriptives_onlypretone_coh.CIupper(1, 8:14) ,descriptives_onlypretone_coh.CIlower(1, 8:14), 'green');
plot(descriptives_onlypretone_coh.avg(1, 8:14), 'color', [0 0 0]);
title ('alpha')
ylabel ('coherence')
xlabel ('frequency (Hz)')
xlim([1 7])
xticks(2:2:6);  % Set custom tick positions
xticklabels({'9';'11';'13'})
hold off
end

if strcmp(Frequency_Band,'beta') == 1
hold on
ciplot(descriptives_onlyprior_coh.CIupper(1, 15:30) ,descriptives_onlyprior_coh.CIlower(1, 15:30), 'blue');
plot(descriptives_onlyprior_coh.avg(1, 15:30), 'color', [0 0 0]);
ciplot(descriptives_onlypretone_coh.CIupper(1, 15:30) ,descriptives_onlypretone_coh.CIlower(1, 15:30), 'green');
plot(descriptives_onlypretone_coh.avg(1, 15:30), 'color', [0 0 0]);
title ('beta')
legend('LED Only', '' ,'Pretone Only','')
ylabel ('coherence')
xlabel ('frequency (Hz)')
xlim([1 16])
xticks(2:2:16);  % Set custom tick positions
xticklabels({'16';'18';'20';'22';'24';'26';'28';'30'})
hold off
end

if strcmp(Frequency_Band,'gamma') == 1
hold on
ciplot(descriptives_onlyprior_coh.CIupper(1, 31:54) ,descriptives_onlyprior_coh.CIlower(1, 31:54), 'blue');
plot(descriptives_onlyprior_coh.avg(1, 31:54), 'color', [0 0 0]);
ciplot(descriptives_onlypretone_coh.CIupper(1, 31:54) ,descriptives_onlypretone_coh.CIlower(1, 31:54), 'green');
plot(descriptives_onlypretone_coh.avg(1, 31:54), 'color', [0 0 0]);
title ('gamma')
ylabel ('coherence')
xlabel ('frequency (Hz)')
xlim([1 24])
xticks(2:4:24);  % Set custom tick positions
xticklabels({'32';'36';'40';'44';'48';'52';'54'});
hold off
end

if strcmp(Frequency_Band,'highGamma') == 1
hold on
ciplot(descriptives_onlyprior_coh.CIupper(1, 66:100) ,descriptives_onlyprior_coh.CIlower(1, 66:100), 'blue');
plot(descriptives_onlyprior_coh.avg(1, 66:100), 'color', [0 0 0]);
ciplot(descriptives_onlypretone_coh.CIupper(1, 66:100) ,descriptives_onlypretone_coh.CIlower(1, 66:100), 'green');
plot(descriptives_onlypretone_coh.avg(1, 66:100), 'color', [0 0 0]);
title ('high gamma')
ylabel ('coherence')
xlabel ('frequency (Hz)')
xlim([1 35])
xticks(6:4:32);  % Set custom tick positions
xticklabels({'72';'76';'80';'84';'88';'92';'96'});
hold off
end


if strcmp(Frequency_Band,'fullSpectra') == 1
hold on
ciplot(descriptives_onlyprior_coh.CIupper,descriptives_onlyprior_coh.CIlower, 'blue');
plot(descriptives_onlyprior_coh.avg, 'color', [0 0 0]);
ciplot(descriptives_onlypretone_coh.CIupper,descriptives_onlypretone_coh.CIlower, 'green');
plot(descriptives_onlypretone_coh.avg, 'color', [0 0 0]);
title ('full spectra')
ylabel ('coherence')
xlabel ('frequency (Hz)')
xlim([1 101])
hold off
end

end

