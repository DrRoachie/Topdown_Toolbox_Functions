function [] = Plot_MonteCarloDist(Frequency_Band,monte_max_cluster_sum_Maris,data_max_cluster_sum_Maris)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

figure()

        if strcmp(Frequency_Band,'theta') == 1 
        histogram(monte_max_cluster_sum_Maris.theta, 100);
        hold on 
        line([data_max_cluster_sum_Maris.theta ,data_max_cluster_sum_Maris.theta ], ylim, 'Color', 'r', 'LineWidth', 2)
        ylabel ('count')
        xlabel ('max cluster value')
        hold off
        title('theta')
        end
        
        if strcmp(Frequency_Band,'alpha') == 1 
        histogram(monte_max_cluster_sum_Maris.alpha, 100);
        hold on 
        line([data_max_cluster_sum_Maris.alpha,data_max_cluster_sum_Maris.alpha ], ylim, 'Color', 'r', 'LineWidth', 2)
        ylabel ('count')
        xlabel ('max cluster value')
        hold off
        title('alpha')
        end
        
        if strcmp(Frequency_Band,'beta') == 1 
        histogram(monte_max_cluster_sum_Maris.beta, 100);
        hold on 
        line([data_max_cluster_sum_Maris.beta ,data_max_cluster_sum_Maris.beta], ylim, 'Color', 'r', 'LineWidth', 2)
        ylabel ('count')
        xlabel ('max cluster value')
        hold off
        title('beta')
        legend('permutation', 'data')
        end
        
        if strcmp(Frequency_Band,'gamma') == 1 
        histogram(monte_max_cluster_sum_Maris.gamma, 100);
        hold on 
        line([data_max_cluster_sum_Maris.gamma ,data_max_cluster_sum_Maris.gamma], ylim, 'Color', 'r', 'LineWidth', 2)
        ylabel ('count')
        xlabel ('max cluster value')
        hold off
        title('gamma')
        end
        
        if strcmp(Frequency_Band,'highGamma') == 1 
        histogram(monte_max_cluster_sum_Maris.highGamma, 100);
        hold on 
        line([data_max_cluster_sum_Maris.highGamma ,data_max_cluster_sum_Maris.highGamma], ylim, 'Color', 'r', 'LineWidth', 2)
        ylabel ('count')
        xlabel ('max cluster value')
        hold off
        title('high gamma')
        end
        
        if strcmp(Frequency_Band,'fullSpectrum') == 1 
        histogram(monte_max_cluster_sum_Maris.all, 100);
        hold on 
        line([data_max_cluster_sum_Maris.all ,data_max_cluster_sum_Maris.all], ylim, 'Color', 'r', 'LineWidth', 2)
        ylabel ('count')
        xlabel ('max cluster value')
        hold off
        title('all')
        end


end


