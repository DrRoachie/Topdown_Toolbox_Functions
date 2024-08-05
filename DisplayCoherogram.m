% 6/21/22 rewrite the code to reduce the redundancy of the data
% DATA_DIR = 'G:\LFP\Coherence'; % path for coherence file
DATA_DIR = 'G:\LFP\TimeFrequency'; % path for coherence file

SessionName = 'MrCassius-190421_pad'; %'MrCassius-190421'; %'MrMiyagi-190904';
Epoch       = 'preCueOnset'; %'moveOnset';
Condition   = 'Both';
eID_ac      = 'D2'; % electrode ID in 1st column (usually AC electrode)
eID_pfc     = 'D1'; % electrode ID in 2nd column (usually PFC electrode)

% load data
% fName = strcat('coherence_',Condition);
% load(fullfile(DATA_DIR,fName));
% fName = strcat('Coherence_',Epoch,'_',Condition);
fName = strcat('TimeFrequency_',Epoch,'_',Condition,'_v2');
load(fullfile(DATA_DIR,SessionName,fName));

% % plot power spectrum
% S_w = abs(freq_w.fourierspctrm);
% f = freq_w.freq;
% mS_w = squeeze(mean(S_w,1)); % trial average
% 
% mS_elec_w = mS_w(contains(label,'D4'),:); % choose electrode
% figure;
% plot(f,mS_elec_w);

% choose electrode combination
C_c = reshape_coherogram(tcoh_c,eID_ac,eID_pfc); % correct trial
C_w = reshape_coherogram(tcoh_w,eID_ac,eID_pfc); % wrong trial
f   = tcoh_c.freq;
t   = tcoh_c.time;

% % % average across 1st electrode (eID_ac) % % %
% ch = 13:14; % select PFC channel for display
aveC1_c = squeeze(mean(C_c.cohspctrm_mat,1)); 
aveC1_w = squeeze(mean(C_w.cohspctrm_mat,1));
aveC1_d = aveC1_c - aveC1_w; % difference (correct - wrong)
% string_a = strcat('PFC_D_', eID_pfc(end), ' -- AC_D_', eID_ac(end));
string_a = strcat('PFC_D_', eID_pfc(end));
nCh = size(aveC1_c,1);
% display coherogram
figure('Position',[50 50 1200 700]);
if nCh==20
    for i=1:nCh
        subplot(4,5,i);
        plot_coherogram(t,f,aveC1_c,string_a,i+2);
    end
elseif nCh==14
    for i=1:nCh
        subplot(3,5,i);
        plot_coherogram(t,f,aveC1_c,string_a,i+1);
    end
end


% % % average across 2nd electrode (eID_pfc) % % %
% ch = 8:11; % select AC channel for display

aveC2_c = squeeze(mean(C_c.cohspctrm_mat,2));
aveC2_w = squeeze(mean(C_w.cohspctrm_mat,2));
aveC2_d = aveC2_c - aveC2_w; % difference(correct - wrong)
string_b = strcat('AC_D_', eID_ac(end));
nCh = size(aveC2_c,1);
% display coherogram
figure('Position',[75 75 1200 700]);
if nCh==20
    for i=1:nCh
        subplot(4,5,i);
        plot_coherogram(t,f,aveC2_c,string_b,i+2);
    end
elseif nCh==14
    for i=1:nCh
        subplot(3,5);
        plot_coherogram(t,f,aveC2_c,string_b,i+1);
    end
end

figure('Position',[100 100 1200 700]);
% % average across channels
ch_a = 7:13; % select PFC channel for display
subplot(2,3,1);
plot_coherogram(t,f,aveC1_c,string_a,ch_a);
subplot(2,3,2);
plot_coherogram(t,f,aveC1_w,string_a,ch_a);
subplot(2,3,3);
plot_coherogram(t,f,aveC1_d,string_a,ch_a);

ch_b = 3:6; % select AC channel for display
subplot(2,3,4);
plot_coherogram(t,f,aveC2_c,string_b,ch_b);
subplot(2,3,5);
plot_coherogram(t,f,aveC2_w,string_b,ch_b);
subplot(2,3,6);
plot_coherogram(t,f,aveC2_d,string_b,ch_b);

% % plot coherence between areas...
% nAveCh = 20; % number of averaging channel
% C_c = tcoh_c.cohspctrm;
% C_w = tcoh_w.cohspctrm;
% t   = tcoh_c.time;
% f   = tcoh_c.freq;
% N = size(C_c,1); % total combination of channels
% for i=1:(N/nAveCh)
%     i_start = nAveCh * (i-1) + 1;
%     i_end   = nAveCh * i;
%     aveCa_c(:,:,i) = squeeze(mean(C_c(i_start:i_end,:,:),1));
%     aveCa_w(:,:,i) = squeeze(mean(C_w(i_start:i_end,:,:),1));
%     chID{i} = tcoh_c.labelcmb(i_start,2);
% end
% ACD3_c = aveCa_c(:,:,1:2:end);
% ACD3_w = aveCa_w(:,:,1:2:end);
% ACD4_c = aveCa_c(:,:,2:2:end);
% ACD4_w = aveCa_w(:,:,2:2:end);
% chID_PFC = chID(1:2:end);
% 
% ACD3_PFCD1_c = ACD3_c(:,:,1:14);
% ACD3_PFCD2_c = ACD3_c(:,:,15:28);
% ACD4_PFCD1_c = ACD4_c(:,:,1:14);
% ACD4_PFCD2_c = ACD4_c(:,:,15:28);
% 
% ACD3_PFCD1_w = ACD3_w(:,:,1:14);
% ACD3_PFCD2_w = ACD3_w(:,:,15:28);
% ACD4_PFCD1_w = ACD4_w(:,:,1:14);
% ACD4_PFCD2_w = ACD4_w(:,:,15:28);
% 
% chLabel_PFC = {'ch02','ch03','ch04','ch05','ch06','ch07','ch08','ch09','ch10','ch11','ch12','ch13','ch14','ch15'};
% 
% ch = 12:13; % choose PFC channels ch13 & ch14
% % ch = 1:14;
% figure('Position',[100 100 1250 750]);
% subplot(2,3,1);
% imagesc(t,f,mean(ACD3_PFCD1_c(:,:,ch),3)); set(gca,'YDir','normal');
% % imagesc(t,f,mean(ACD4_PFCD1_c(:,:,ch),3)); set(gca,'YDir','normal');
% xlabel('Time [sec]'); ylabel('Frequency [Hz]');
% title('Correct');
% subplot(2,3,2);
% imagesc(t,f,mean(ACD3_PFCD1_w(:,:,ch),3)); set(gca,'YDir','normal');
% % imagesc(t,f,mean(ACD4_PFCD1_w(:,:,ch),3)); set(gca,'YDir','normal');
% xlabel('Time [sec]'); ylabel('Frequency [Hz]');
% title('Wrong');
% subplot(2,3,3);
% imagesc(t,f,mean(ACD3_PFCD1_c(:,:,ch) - ACD3_PFCD1_w(:,:,ch),3)); set(gca,'YDir','normal');
% % imagesc(t,f,mean(ACD4_PFCD1_c(:,:,ch) - ACD4_PFCD1_w(:,:,ch),3)); set(gca,'YDir','normal');
% xlabel('Time [sec]'); ylabel('Frequency [Hz]');
% title('Difference');
% colormap jet;
% 
% 
% % plot coherence between areas...
% nAveCh = 14; % number of averaging channel
% C2_c = tcoh2_c.cohspctrm;
% C2_w = tcoh2_w.cohspctrm;
% t   = tcoh2_c.time;
% f   = tcoh2_c.freq;
% N = size(C2_c,1); % total combination of channels
% for i=1:(N/nAveCh)
%     i_start = nAveCh * (i-1) + 1;
%     i_end   = nAveCh * i;
%     aveCb_c(:,:,i) = squeeze(mean(C2_c(i_start:i_end,:,:),1));
%     aveCb_w(:,:,i) = squeeze(mean(C2_w(i_start:i_end,:,:),1));
%     chID{i} = coh2_c.labelcmb(i_start,2);
% end
% PFCD1_c = aveCb_c(:,:,1:2:end);
% PFCD1_w = aveCb_w(:,:,1:2:end);
% PFCD2_c = aveCb_c(:,:,2:2:end);
% PFCD2_w = aveCb_w(:,:,2:2:end);
% chID_AC = chID(1:2:end);
% 
% PFCD1_ACD3_c = PFCD1_c(:,:,1:20);
% PFCD1_ACD4_c = PFCD1_c(:,:,21:40);
% PFCD2_ACD3_c = PFCD2_c(:,:,1:20);
% PFCD2_ACD4_c = PFCD2_c(:,:,21:40);
% 
% PFCD1_ACD3_w = PFCD1_w(:,:,1:20);
% PFCD1_ACD4_w = PFCD1_w(:,:,21:40);
% PFCD2_ACD3_w = PFCD2_w(:,:,1:20);
% PFCD2_ACD4_w = PFCD2_w(:,:,21:40);
% 
% chLabel_AC = {'ch03','ch04','ch05','ch06','ch07','ch08','ch09','ch10','ch11','ch12','ch13','ch14','ch15','ch16','ch17','ch18','ch19','ch20','ch21','ch22'};
% 
% ch = 6:9; % choose AC channel ch8 - ch11
% % ch = 3:4; % choose AC channel ch5 - ch6
% % ch = 1:20;
% % figure;
% subplot(2,3,4);
% imagesc(t,f,mean(PFCD1_ACD3_c(:,:,ch),3)); set(gca,'YDir','normal');
% % imagesc(t,f,mean(PFCD1_ACD4_c(:,:,ch),3)); set(gca,'YDir','normal');
% xlabel('Time [sec]'); ylabel('Frequency [Hz]');
% title('Correct');
% subplot(2,3,5);
% imagesc(t,f,mean(PFCD1_ACD3_w(:,:,ch),3)); set(gca,'YDir','normal');
% % imagesc(t,f,mean(PFCD1_ACD4_w(:,:,ch),3)); set(gca,'YDir','normal');
% xlabel('Time [sec]'); ylabel('Frequency [Hz]');
% title('Wrong');
% subplot(2,3,6);
% imagesc(t,f,mean(PFCD1_ACD3_c(:,:,ch) - PFCD1_ACD3_w(:,:,ch),3)); set(gca,'YDir','normal');
% % imagesc(t,f,mean(PFCD1_ACD4_c(:,:,ch) - PFCD1_ACD4_w(:,:,ch),3)); set(gca,'YDir','normal');
% xlabel('Time [sec]'); ylabel('Frequency [Hz]');
% title('Difference');
% colormap jet;