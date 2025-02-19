function [R_spearman, R_spearmanNull, stat_1] = PSIScratchPad(Frequency_Band, PSI_Condition_1, PSI_Condition_2)

if strcmp(Frequency_Band,'theta') == 1

% Extract theta band columns 4-8
cols_to_plot_1 = PSI_Condition_1.psispctrm(:, 4:8);
cols_to_plot_2 = PSI_Condition_2.psispctrm(:, 4:8);

% generate random pool of all PSI values in the selected band
C = [cols_to_plot_1, cols_to_plot_2];
shuffledC_temp = C(randperm(numel(C)));

% Reshape shuffled values back to the size of C
shuffledC = reshape(shuffledC_temp, size(C));


% average the band
mean_PSI_1 = mean(cols_to_plot_1, 2);
mean_PSI_2 = mean(cols_to_plot_2, 2);

% Sort out Channels talking in the PFC-to-AC direction 
% Get the channel pairs that are talking in the PFC-to-AC direction for Condition 1,then pull the PSI values those  channel pairs in Condition 2

positive_indices_PSI_1 = find(mean_PSI_1 > 0);

% Pull the values from PSI_1 and PSI_2 at the positive indices
positive_values_PSI_1 = mean_PSI_1(positive_indices_PSI_1);
positive_values_PSI_2 = mean_PSI_2(positive_indices_PSI_1);  % Use the same indices for PSI_2

% Determine the number of positive indices in PSI_1
num_positive_indices = length(positive_indices_PSI_1);

% Randomly select indices
selected_pos_random_values_temp = shuffledC(randperm(length(shuffledC), num_positive_indices));
selected_pos_random_values_1 = selected_pos_random_values_temp(:);
selected_pos_random_values_temp = shuffledC(randperm(length(shuffledC), num_positive_indices));
selected_pos_random_values_2 = selected_pos_random_values_temp(:);

% Sort out Channels talking in the AC-to-PFC direction

negative_indices_PSI_1 = find(mean_PSI_1 < 0);

% Pull the values from cat_1 and cat_2 at the negative indices
negative_values_PSI_1 = mean_PSI_1(negative_indices_PSI_1);
negative_values_PSI_2 = mean_PSI_2(negative_indices_PSI_1);  % Use the same indices for cat_2

% Determine the number of negative indices in PSI_1
num_negative_indices = length(negative_indices_PSI_1);

% Randomly select indices
selected_neg_random_values_temp = shuffledC(randperm(length(shuffledC), num_negative_indices));
selected_neg_random_values_1 = selected_neg_random_values_temp(:);
selected_neg_random_values_temp = shuffledC(randperm(length(shuffledC), num_negative_indices));
selected_neg_random_values_2 = selected_neg_random_values_temp(:);

% Sort out Channels talking in the AC-to-PFC direction

% statistics

R_spearman.PFC_to_AC  = corr(positive_values_PSI_1, positive_values_PSI_2, 'Type', 'Spearman');
R_spearman.AC_to_PFC  = corr(negative_values_PSI_1, negative_values_PSI_2, 'Type', 'Spearman');
R_spearmanNull.PFC_to_AC   = corr(selected_pos_random_values_1, selected_pos_random_values_2, 'Type', 'Spearman');
R_spearmanNull.AC_to_PFC    = corr(selected_neg_random_values_1, selected_neg_random_values_2, 'Type', 'Spearman');
[stat_1.pfc_to_ac_pvalaue, stat_1.pfc_to_ac_h, stat_1.pfc_to_ac_teststatistic]  = ranksum(positive_values_PSI_1, positive_values_PSI_2);
[stat_1.ac_to_pfc_pvalue,  stat_1.ac_to_pfc_h, stat_1.ac_to_pfc_teststatistic]  = ranksum(negative_values_PSI_1, negative_values_PSI_2);

% Plot the results 

figure; 

subplot(2, 2, 1);
hold on
scatter(positive_values_PSI_1, positive_values_PSI_2);
% Draw the line representing a one-to-one relationship (slope = 1)
x_limits = [min(positive_values_PSI_1), max(positive_values_PSI_1)];
y_limits = [min(positive_values_PSI_2), max(positive_values_PSI_2)];
%y_limits = [min(0), max(0)];
line(x_limits, y_limits, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
xlabel('Condition 1 PSI'); 
ylabel('Condition 2 PSI');
title('PFC to AC Communication in Theta')
hold off 

subplot(2, 2, 2)
hold on 
scatter(negative_values_PSI_1, negative_values_PSI_2);
x_limits = [min(negative_values_PSI_1), max(negative_values_PSI_1)];
y_limits = [min(negative_values_PSI_2), max(negative_values_PSI_2)];
%y_limits = [min(0), max(0)];
line(x_limits, y_limits, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
xlabel('Condition 1 PSI'); 
ylabel('Condition 2 PSI'); 
title('AC to PFC Communication in Theta')
hold off 

subplot(2, 2, 3);
hold on
scatter(selected_pos_random_values_2, selected_pos_random_values_1);
x_limits = [min(selected_pos_random_values_2), max(selected_pos_random_values_2)];
y_limits = [min(selected_pos_random_values_1), max(selected_pos_random_values_1)];
%y_limits = [min(0), max(0)];
line(x_limits, y_limits, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
xlabel('NULL Condition 1 PSI'); 
ylabel('NULL Condition 2 PSI');
title('NULL PFC to AC Communication in Theta')
hold off 

subplot(2, 2, 4);
hold on 
scatter(selected_neg_random_values_2, selected_neg_random_values_1);
x_limits = [min(selected_neg_random_values_2), max(selected_neg_random_values_2)];
y_limits = [min(selected_neg_random_values_1), max(selected_neg_random_values_1)];
%y_limits = [min(0), max(0)];
line(x_limits, y_limits, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '1:1 Line');
xlabel('NULL Condition 1 PSI'); 
ylabel('NULL Condition 2 PSI'); 
title('NULL AC to PFC Communication in Theta')
hold off 

end

