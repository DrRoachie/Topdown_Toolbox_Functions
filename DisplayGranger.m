%%
DATADIR = 'C:\Users\Corey\Documents\00_DATA\Coherence_Analysis_Data';

SessionName = 'MrCassius-190421'; %'MrCassius-190429'; %'MrMiyagi-190904';
Epoch       = 'testToneOnset';
Condition   = 'Both'; %'Pretone';
eID_AC      = 'D2'; % electrode ID in 1st column (usually AC electrode)
eID_PFC     = 'D1'; % electrode ID in 2nd column (usually PFC electrode)


%% load data%
fName = strcat('Coherence_',Epoch,'_',Condition);
%fName = strcat('Frequency_',Epoch,'_',Condition);
load(fullfile(DATADIR,SessionName,fName));

%%
label = granger_c.label;
i_AC  = contains(label,eID_AC);
i_PFC = contains(label,eID_PFC);
index = or(i_AC,i_PFC);
label_AC = label(i_AC);
label_PFC = label(i_PFC);
clear i_AC i_PFC

%% select electrode for the display of Granger causality
G_c = granger_c;
G_w = granger_w;
g_c = G_c.grangerspctrm(index,index,:);
g_w = G_w.grangerspctrm(index,index,:);
label_select = label(index);
G_c.grangerspctrm = g_c;
G_c.label = label_select;
G_c.cfg.channel = label_select;
G_w.grangerspctrm = g_w;
G_w.label = label_select;
G_w.cfg.channel = label_select;

%%

% quick view of Granger causality of all combination...
cfg = [];
cfg.parameter = 'grangerspctrm';
cfg.xlim = [0 150];
% cfg.xlim = [65 120];
cfg.zlim = [0 0.3];
cfg.channel = {'*D1*' '*D2*'}; % choose electrode to show...
ft_connectivityplot(cfg,G_c);

%%
ch_AC  = 1;
ch_PFC = 1;
gg_c = G_c.grangerspctrm;
gg_w = G_w.grangerspctrm;
f = G_c.freq; % frequency
i_AC  = findChannel(label_select,'AC',ch_AC);
i_PFC = findChannel(label_select,'PFC',ch_PFC);
gAC2PFC_c = squeeze(gg_c(i_AC,i_PFC,:));
gAC2PFC_w = squeeze(gg_w(i_AC,i_PFC,:));
gPFC2AC_c = squeeze(gg_c(i_PFC,i_AC,:));
gPFC2AC_w = squeeze(gg_w(i_PFC,i_AC,:));


%% plot
figure;
subplot(2,2,2);
plot(f,gPFC2AC_c); hold on;
plot(f,gPFC2AC_w);
subplot(2,2,3);
plot(f,gAC2PFC_c); hold on;
plot(f,gAC2PFC_w);
