%function [Summary] = PSIHistogram(averagedGrid_1_pos, averagedGrid_1_neg, averagedGrid_2_pos, averagedGrid_2_neg, nBins) 

    % Ensure all inputs are column vectors
    averagedGrid_1_pos = nonzeros(averagedGrid_1_pos(:));
    averagedGrid_1_neg = nonzeros(averagedGrid_1_neg(:));
    averagedGrid_2_pos = nonzeros(averagedGrid_2_pos(:));
    averagedGrid_2_neg = nonzeros(averagedGrid_2_neg(:));

    % Define shared bin edges based on absolute values of both conditions
    [~, bins] = hist([averagedGrid_1_pos; abs(averagedGrid_1_neg); averagedGrid_2_pos; abs(averagedGrid_2_neg)], nBins);

    % Compute histograms for each dataset
    posCounts_1 = hist(averagedGrid_1_pos, bins);
    negCounts_1 = hist(abs(averagedGrid_1_neg), bins);
    posCounts_2 = hist(averagedGrid_2_pos, bins);
    negCounts_2 = hist(abs(averagedGrid_2_neg), bins);

    % Normalize counts to percentages
    posCounts_1 = posCounts_1 / sum(posCounts_1) * 100;
    negCounts_1 = negCounts_1 / sum(negCounts_1) * 100;
    posCounts_2 = posCounts_2 / sum(posCounts_2) * 100;
    negCounts_2 = negCounts_2 / sum(negCounts_2) * 100;

    % Plot histograms
    figure;
    hold on;

    % Condition 1: Red (Positive) & Soft Red (Negative)
    bar(bins,  posCounts_1, 'FaceColor', [1, 0, 0], 'EdgeColor', 'k', 'BarWidth',    1); % Red with black edges
    bar(bins, -negCounts_1, 'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1); % Soft Red with black edges

    % Condition 2: Blue (Positive) & Soft Blue (Negative)
    bar(-bins,  posCounts_2, 'FaceColor', [0, 0, 1], 'EdgeColor', 'k', 'BarWidth', 1); % Blue with black edges
    bar(-bins, -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

    % Labels and Legend
    xlabel('PSI Value');
    ylabel('Proportion of Sites');
    ylim([-50 50])
    %xlim([-.15 .15])
    legend({'Cong.   top-down', 'Incong. top-down', 'Cong.   bottom-up', 'Incong. bottom-up'}, 'Location', 'SouthEast');

    % Improve aesthetics
    set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Adjust axis font size and line width
    box on;

    hold off;

    %%

PosNegHistogram(PriorOnly_TestTone_HighSNR_Correct_Incongruent_PFC_2_AC, PriorOnly_TestTone_HighSNR_Correct_Incongruent_AC_2_PFC, PretoneOnly_TestTone_HighSNR_Correct_Incongruent_PFC_2_AC, PretoneOnly_TestTone_HighSNR_Correct_Incongruent_AC_2_PFC,20)