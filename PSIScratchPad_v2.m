function [stat_1, pvalue] = PSIScratchPad_v2(Frequency_Band, PSI_Condition_1, PSI_Condition_2)

if strcmp(Frequency_Band,'theta') == 1

% Calculate everything you need to do for the data. 

% Extract theta band columns 4-8
PSI_Condition_1_band = PSI_Condition_1.psispctrm(:, 4:8);
PSI_Condition_2_band = PSI_Condition_2.psispctrm(:, 4:8);

PSI_Condition_1_band_col = PSI_Condition_1_band(:);
PSI_Condition_2_band_col = PSI_Condition_2_band(:);

% Average the band
mean_PSI_1 = mean(PSI_Condition_1_band, 2);
mean_PSI_2 = mean(PSI_Condition_2_band, 2);

% Sort out Channels talking in the PFC-to-AC direction 
% Get the channel pairs that are talking in the PFC-to-AC direction for Condition 1,then pull the PSI values those  channel pairs in Condition 2

positive_indices_PSI_1 = find(mean_PSI_1 > 0);

% Pull the values from PSI_1 and PSI_2 at the positive indices
positive_values_PSI_1 = mean_PSI_1(positive_indices_PSI_1);
positive_values_PSI_2 = mean_PSI_2(positive_indices_PSI_1);  % Use the same indices for PSI_2

% Determine the number of positive indices in PSI_1
num_positive_indices = length(positive_indices_PSI_1);

% Sort out Channels talking in the AC-to-PFC direction

negative_indices_PSI_1 = find(mean_PSI_1 < 0);

% Pull the values from cat_1 and cat_2 at the negative indices
negative_values_PSI_1 = mean_PSI_1(negative_indices_PSI_1);
negative_values_PSI_2 = mean_PSI_2(negative_indices_PSI_1);  % Use the same indices for cat_2

% Determine the number of negative indices in PSI_1
num_negative_indices = length(negative_indices_PSI_1);

% Initialize stat_1 with NaN values in case of empty inputs
stat_1.pfc_to_ac_pvalue = NaN;
stat_1.pfc_to_ac_h = NaN;
stat_1.pfc_to_ac_teststatistic = NaN;
stat_1.ac_to_pfc_pvalue = NaN;
stat_1.ac_to_pfc_h = NaN;
stat_1.ac_to_pfc_teststatistic = NaN;
pvalue.pfc_to_ac = NaN;
pvalue.ac_to_pfc = NaN;

% Carry out the Wilcox tests to get zvalues for the comparison of PSI1 and PSI2 for both directions 

% Check if any of the arguments are empty
if ~isempty(positive_values_PSI_1) && ~isempty(positive_values_PSI_2) && ...
   ~isempty(negative_values_PSI_1) && ~isempty(negative_values_PSI_2)
    
    % Perform ranksum test for positive values
    [stat_1.pfc_to_ac_pvalaue, stat_1.pfc_to_ac_h, stat_1.pfc_to_ac_teststatistic] = ...
        ranksum(positive_values_PSI_1, positive_values_PSI_2, 'method', 'approximate');

    % Perform ranksum test for negative values
    [stat_1.ac_to_pfc_pvalue,  stat_1.ac_to_pfc_h, stat_1.ac_to_pfc_teststatistic] = ...
        ranksum(negative_values_PSI_1, negative_values_PSI_2, 'method', 'approximate');
    
    % Get the separate null distributions for the PS

Bootstrap_1 = PSIBootstrap(PSI_Condition_1_band_col, PSI_Condition_2_band_col);
Bootstrap_2 = PSIBootstrap(PSI_Condition_1_band_col, PSI_Condition_2_band_col);

% % Plot the data on the null dist
% 
%         figure
%         subplot(1,2,1)
%         histogram(Bootstrap_1, 100);
%         hold on 
%         line([stat_1.pfc_to_ac_teststatistic.zval ,stat_1.pfc_to_ac_teststatistic.zval], ylim, 'Color', 'r', 'LineWidth', 2)
%         ylabel ('count')
%         xlabel ('rank sum zval')
%         title ('PFC-to-AC')
%         hold off
% 
%         subplot(1,2,2)
%         histogram(Bootstrap_2, 100);
%         hold on 
%         line([stat_1.ac_to_pfc_teststatistic.zval ,stat_1.ac_to_pfc_teststatistic.zval], ylim, 'Color', 'r', 'LineWidth', 2)
%         ylabel ('count')
%         xlabel ('rank sum zvalue')
%         title ('AC-to-PFC')
%         hold off

 pvalue.pfc_to_ac = sum(Bootstrap_1 <= stat_1.pfc_to_ac_teststatistic.zval) / numel(Bootstrap_1);

 if pvalue.pfc_to_ac > .5 
    pvalue.pfc_to_ac = 1 - pvalue.pfc_to_ac;
 end 

 %
 pvalue.ac_to_pfc = sum(Bootstrap_2 <= stat_1.ac_to_pfc_teststatistic.zval) / numel(Bootstrap_2);

 if pvalue.ac_to_pfc > .5 
    pvalue.ac_to_pfc = 1 - pvalue.ac_to_pfc;
 end 

% Get the positive and negative null distribution

% generate random pool of all PSI values in the selected band
C = [PSI_Condition_1_band, PSI_Condition_2_band];
shuffledC_temp = C(randperm(numel(C)));

% Reshape shuffled values back to the size of C
shuffledC = reshape(shuffledC_temp, size(C));

% Randomly select indices
selected_pos_random_values_temp = shuffledC(randperm(length(shuffledC), num_positive_indices));
selected_pos_random_values_1 = selected_pos_random_values_temp(:);
selected_pos_random_values_temp = shuffledC(randperm(length(shuffledC), num_positive_indices));
selected_pos_random_values_2 = selected_pos_random_values_temp(:);

% Randomly select indices
selected_neg_random_values_temp = shuffledC(randperm(length(shuffledC), num_negative_indices));
selected_neg_random_values_1 = selected_neg_random_values_temp(:);
selected_neg_random_values_temp = shuffledC(randperm(length(shuffledC), num_negative_indices));
selected_neg_random_values_2 = selected_neg_random_values_temp(:);

% % Plot the results 
% 
% figure; 
% 
% subplot(2, 2, 1);
% hold on
% scatter(positive_values_PSI_1, positive_values_PSI_2);
% % Draw the line representing a one-to-one relationship (slope = 1)
% %x_limits = [min(positive_values_PSI_1), max(positive_values_PSI_1)];
% %y_limits = [min(positive_values_PSI_2), max(positive_values_PSI_2)];
% 
% xy_min=min([min(positive_values_PSI_1) min(positive_values_PSI_2)]);
% xy_max=max([max(positive_values_PSI_1) max(positive_values_PSI_2)]);
% 
% %y_limits = [min(0), max(0)];
% line([xy_min xy_max], [xy_min xy_max], 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
% xlabel('Condition 1 PSI'); 
% ylabel('Condition 2 PSI');
% title('PFC to AC Communication in Theta')
% hold off 
% 
% subplot(2, 2, 2)
% hold on 
% scatter(negative_values_PSI_1, negative_values_PSI_2);
% xy_min=min([min(negative_values_PSI_1) min(negative_values_PSI_2)]);
% xy_max=max([max(negative_values_PSI_1) max(negative_values_PSI_2)]);
% 
% %y_limits = [min(0), max(0)];
% line([xy_min xy_max], [xy_min xy_max], 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
% xlabel('Condition 1 PSI'); 
% ylabel('Condition 2 PSI'); 
% title('AC to PFC Communication in Theta')
% hold off 
% 
% subplot(2, 2, 3);
% hold on
% scatter(selected_pos_random_values_2, selected_pos_random_values_1);
% xy_min=min([min(selected_pos_random_values_1) min(selected_pos_random_values_2)]);
% xy_max=max([max(selected_pos_random_values_1) max(selected_pos_random_values_2)]);
% 
% %y_limits = [min(0), max(0)];
% line([xy_min xy_max], [xy_min xy_max], 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
% xlabel('NULL Condition 1 PSI'); 
% ylabel('NULL Condition 2 PSI');
% title('NULL PFC to AC Communication in Theta')
% hold off 
% 
% subplot(2, 2, 4);
% hold on 
% scatter(selected_neg_random_values_2, selected_neg_random_values_1);
% xy_min=min([min(selected_neg_random_values_1) min(selected_neg_random_values_2)]);
% xy_max=max([max(selected_neg_random_values_1) max(selected_neg_random_values_2)]);
% %y_limits = [min(0), max(0)];
% line([xy_min xy_max], [xy_min xy_max], 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
% xlabel('NULL Condition 1 PSI'); 
% ylabel('NULL Condition 2 PSI'); 
% title('NULL AC to PFC Communication in Theta')
% hold off 


else
    disp('Skipping ranksum tests due to empty arguments.');

end

end



end

% Get the positive and negative null distribution

% generate random pool of all PSI values in the selected band
% C = [PSI_Condition_1_band, PSI_Condition_2_band];
% shuffledC_temp = C(randperm(numel(C)));
% 
% % Reshape shuffled values back to the size of C
% shuffledC = reshape(shuffledC_temp, size(C));
% 
% % Randomly select indices
% selected_pos_random_values_temp = shuffledC(randperm(length(shuffledC), num_positive_indices));
% selected_pos_random_values_1 = selected_pos_random_values_temp(:);
% selected_pos_random_values_temp = shuffledC(randperm(length(shuffledC), num_positive_indices));
% selected_pos_random_values_2 = selected_pos_random_values_temp(:);
% 
% % Sort out Channels talking in the AC-to-PFC direction
% 
% negative_indices_PSI_1 = find(mean_PSI_1 < 0);
% 
% % Pull the values from cat_1 and cat_2 at the negative indices
% negative_values_PSI_1 = mean_PSI_1(negative_indices_PSI_1);
% negative_values_PSI_2 = mean_PSI_2(negative_indices_PSI_1);  % Use the same indices for cat_2
% 
% % Determine the number of negative indices in PSI_1
% num_negative_indices = length(negative_indices_PSI_1);
% 
% % Randomly select indices
% selected_neg_random_values_temp = shuffledC(randperm(length(shuffledC), num_negative_indices));
% selected_neg_random_values_1 = selected_neg_random_values_temp(:);
% selected_neg_random_values_temp = shuffledC(randperm(length(shuffledC), num_negative_indices));
% selected_neg_random_values_2 = selected_neg_random_values_temp(:);

%% Plot the results 

% figure; 
% 
% subplot(2, 2, 1);
% hold on
% scatter(positive_values_PSI_1, positive_values_PSI_2);
% % Draw the line representing a one-to-one relationship (slope = 1)
% x_limits = [min(positive_values_PSI_1), max(positive_values_PSI_1)];
% y_limits = [min(positive_values_PSI_2), max(positive_values_PSI_2)];
% %y_limits = [min(0), max(0)];
% line(x_limits, y_limits, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
% xlabel('Condition 1 PSI'); 
% ylabel('Condition 2 PSI');
% title('PFC to AC Communication in Theta')
% hold off 
% 
% subplot(2, 2, 2)
% hold on 
% scatter(negative_values_PSI_1, negative_values_PSI_2);
% x_limits = [min(negative_values_PSI_1), max(negative_values_PSI_1)];
% y_limits = [min(negative_values_PSI_2), max(negative_values_PSI_2)];
% %y_limits = [min(0), max(0)];
% line(x_limits, y_limits, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
% xlabel('Condition 1 PSI'); 
% ylabel('Condition 2 PSI'); 
% title('AC to PFC Communication in Theta')
% hold off 
% 
% subplot(2, 2, 3);
% hold on
% scatter(selected_pos_random_values_2, selected_pos_random_values_1);
% x_limits = [min(selected_pos_random_values_2), max(selected_pos_random_values_2)];
% y_limits = [min(selected_pos_random_values_1), max(selected_pos_random_values_1)];
% %y_limits = [min(0), max(0)];
% line(x_limits, y_limits, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
% xlabel('NULL Condition 1 PSI'); 
% ylabel('NULL Condition 2 PSI');
% title('NULL PFC to AC Communication in Theta')
% hold off 
% 
% subplot(2, 2, 4);
% hold on 
% scatter(selected_neg_random_values_2, selected_neg_random_values_1);
% x_limits = [min(selected_neg_random_values_2), max(selected_neg_random_values_2)];
% y_limits = [min(selected_neg_random_values_1), max(selected_neg_random_values_1)];
% %y_limits = [min(0), max(0)];
% line(x_limits, y_limits, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
% xlabel('NULL Condition 1 PSI'); 
% ylabel('NULL Condition 2 PSI'); 
% title('NULL AC to PFC Communication in Theta')
% hold off 


