function [] = Plot_Granger(Send_Chan_lb, Rec_Chan_lb, Granger_PriorOnly, Granger_PretoneOnly, STD_Granger_PriorOnly, STD_Granger_PretoneOnly, Null_Thresh_PriorOnly, Null_Thresh_PretoneOnly)

% sender to receiver
figure;
hold on;

x = Granger_PriorOnly.freq;
y = squeeze(Granger_PriorOnly.grangerspctrm(1,2,:))';
error = squeeze(STD_Granger_PriorOnly(1,2,:))';
PriorOnly_PFC_AC = plot(x,y,'r','DisplayName','PriorOnly'); % plot granger
fill([x fliplr(x)],[y+error fliplr(y-error)],'r', 'FaceAlpha', .2,'EdgeColor','none'); % plot confidence interval
plot(x,squeeze(Null_Thresh_PriorOnly(1,2,:)),'r--'); % plot null threshold

x = Granger_PretoneOnly.freq;
y = squeeze(Granger_PretoneOnly.grangerspctrm(1,2,:))';
error = squeeze(STD_Granger_PretoneOnly(1,2,:))';
PretoneOnly_PFC_AC = plot(x,y,'b','DisplayName','PretoneOnly');
fill([x fliplr(x)],[y+error fliplr(y-error)],'b', 'FaceAlpha', .2,'EdgeColor','none');
plot(x,squeeze(Null_Thresh_PretoneOnly(1,2,:)),'b--');

xlabel('Frequency');
ylabel('Granger Value');
title([Send_Chan_lb ' to ' Rec_Chan_lb]);
legend([PriorOnly_PFC_AC PretoneOnly_PFC_AC]);
hold off;

% receiver to sender
figure;
hold on;

x = Granger_PriorOnly.freq;
y = squeeze(Granger_PriorOnly.grangerspctrm(2,1,:))';
error = squeeze(STD_Granger_PriorOnly(2,1,:))';
PriorOnly_AC_PFC = plot(x,y,'r','DisplayName','PriorOnly');
fill([x fliplr(x)],[y+error fliplr(y-error)],'r', 'FaceAlpha', .2,'EdgeColor','none');
plot(x,squeeze(Null_Thresh_PriorOnly(2,1,:)),'r--');

x = Granger_PretoneOnly.freq;
y = squeeze(Granger_PretoneOnly.grangerspctrm(2,1,:))';
error = squeeze(STD_Granger_PretoneOnly(2,1,:))';
PretoneOnly_AC_PFC = plot(x,y,'b','DisplayName','PretoneOnly');
fill([x fliplr(x)],[y+error fliplr(y-error)],'b', 'FaceAlpha', .2,'EdgeColor','none');
plot(x,squeeze(Null_Thresh_PretoneOnly(2,1,:)),'b--');

xlabel('Frequency');
ylabel('Granger Value');
title([Rec_Chan_lb ' to ' Send_Chan_lb]);
legend([PriorOnly_AC_PFC PretoneOnly_AC_PFC]);
hold off;