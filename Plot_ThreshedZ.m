function [] = Plot_ThreshedZ(Frequency_Band, data_coh_Maris_Z, data_z_statistic_Maris_threshold)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
        
        figure()

        if strcmp(Frequency_Band,'theta') == 1 
        plot(abs(data_coh_Maris_Z.theta)) 
        xlim([1 4])
        xticks(1:4) 
        xticklabels({'5';'6';'7'; '8'})
        ylabel ('Maris Z Value')
        xlabel ('frequency (Hz)')
        hold on
        yline(data_z_statistic_Maris_threshold.theta, '--', 'LineWidth', 2)
        hold off
        title('theta')
        end 
        
        if strcmp(Frequency_Band,'alpha') == 1 
        plot(abs(data_coh_Maris_Z.alpha))
        xlim([1 7])
        xticks(2:2:6);  % Set custom tick positions
        xticklabels({'9';'11';'13'})
        ylabel ('Maris Z Value')
        xlabel ('frequency (Hz)')
        hold on
        yline(data_z_statistic_Maris_threshold.alpha, '--', 'LineWidth', 2)
        hold off
        title('alpha')
        end
        
        if strcmp(Frequency_Band,'beta') == 1 
        plot(abs(data_coh_Maris_Z.beta))
        xlim([1 16])
        xticks(2:2:16);  % Set custom tick positions
        xticklabels({'16';'18';'20';'22';'24';'26';'28';'30'})
        ylabel ('Maris Z Value')
        xlabel ('frequency (Hz)')
        hold on
        yline(data_z_statistic_Maris_threshold.beta, '--', 'LineWidth', 2)
        hold off
        title('beta')
        end
        
        if strcmp(Frequency_Band,'gamma') == 1 
        plot(abs(data_coh_Maris_Z.gamma))
        xlim([1 24])
        xticks(2:4:24);  % Set custom tick positions
        xticklabels({'32';'36';'40';'44';'48';'52';'54'});
        ylabel ('Maris Z Value')
        xlabel ('frequency (Hz)')
        hold on
        yline(data_z_statistic_Maris_threshold.gamma, '--', 'LineWidth', 2)
        hold off
        title('gamma')
        end
        
        if strcmp(Frequency_Band,'highGamma') == 1 
        plot(abs(data_coh_Maris_Z.highGamma))
        xlim([1 35])
        xticks(6:4:32);  % Set custom tick positions
        xticklabels({'72';'76';'80';'84';'88';'92';'96'});
        ylabel ('Maris Z Value')
        xlabel ('frequency (Hz)')
        hold on
        yline(data_z_statistic_Maris_threshold.highGamma, '--', 'LineWidth', 2)
        hold off
        title('high Gamma')
        end
        
        if strcmp(Frequency_Band,'fullSpectrum') == 1 
        plot(abs(data_coh_Maris_Z.all))
        ylabel ('Maris Z Value')
        xlabel ('frequency (Hz)')
        hold on
        yline(data_z_statistic_Maris_threshold.all, '--', 'LineWidth', 2)
        hold off
        title('Full Spectrum')
        end


    end
