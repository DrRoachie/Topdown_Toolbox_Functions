%% Version 1

% function PosNegHistogram(Pos,Neg,nBins)
% 
% %Con1 and Con2 are the two conditions to plot; each is a vector
% 
% %assume always that Pos are positive PSI values
% 
% %assume always that Neg are negative PSI values
% 
% %nbins is the number of bins for the histogram--same for both data sets
% 
% %make sure both are column vectors
% 
% ic =iscolumn(Pos);
% if ~ic
%   Pos=Pos';
% end
% ic =iscolumn(Neg);
% if ~ic
%   Neg=Neg';
% end
% 
% [~ , bins] = hist( [Pos; abs(Neg)] ,nBins);
% posCounts = hist( Pos , bins ); %number of counts per bin for the positive psi
% negCounts = hist( abs(Neg) , bins ); %same but for negative psi
% 
% %normalize
% posCounts=posCounts/sum(posCounts)*100;
% negCounts=negCounts/sum(negCounts)*100;
% 
% bar(bins,posCounts, 'r','BarWidth', 1)
% hold on
% bar(-1*bins,-1*negCounts,'b','BarWidth',1)
% 
%  xlabel('PSI value')
% ylabel('Proportion of sites')

%% Version 2

function [Summary] = PosNegHistogram(averagedGrid_1_pos, averagedGrid_1_neg, averagedGrid_2_pos, averagedGrid_2_neg, nBins) 

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
    bar(bins, posCounts_2, 'FaceColor', [1, 0.6, 0.6], 'FaceAlpha', 0.5,'EdgeColor', 'k', 'BarWidth', 1); % Soft Red with black edges

    % Condition 2: Blue (Positive) & Soft Blue (Negative)
    bar(-bins,  -negCounts_1, 'FaceColor', [0, 0, 1], 'EdgeColor', 'k', 'BarWidth', 1); % Blue with black edges
    bar(-bins, -negCounts_2,  'FaceColor', [0.6, 0.6, 1], 'FaceAlpha', 0.5, 'EdgeColor', 'k', 'BarWidth', 1); % Soft Blue with black edges

    % Labels and Legend
    xlabel('PSI Value');
    ylabel('Proportion of Sites');
    ylim([-100 100])
    legend({'Cong.   top-down', 'Incong. top-down', 'Cong.   bottom-up', 'Incong. bottom-up'}, 'Location', 'SouthEast');

    % Improve aesthetics
    set(gca, 'FontSize', 12, 'LineWidth', 1.5); % Adjust axis font size and line width
    box on;

    hold off;

%% get median, standard deviation, and standard error

Summary.Congruent_PFC_to_AC_median    = median(averagedGrid_1_pos);
Summary.Congruent_PFC_to_AC_std       = std(averagedGrid_1_pos, 0);
Summary.Congruent_PFC_to_AC_stderr    = Summary.Congruent_PFC_to_AC_std / sqrt(length(averagedGrid_1_pos)); 

Summary.Incongruent_PFC_to_AC_median  = median(averagedGrid_2_pos);
Summary.Incongruent_PFC_to_AC_std     = std(averagedGrid_2_pos, 0);
Summary.Incongruent_PFC_to_AC_stderr  = Summary.Incongruent_PFC_to_AC_std  / sqrt(length(averagedGrid_2_pos)); 

Summary.Congruent_AC_to_PFC_median    = median(averagedGrid_1_neg); 
Summary.Congruent_AC_to_PFC_std       = std(averagedGrid_1_neg); 
Summary.Congruent_AC_to_PFC_stderr    = Summary.Congruent_AC_to_PFC_std / sqrt(length(averagedGrid_1_neg)); 

Summary.Incongruent_AC_to_PFC_median  = median(averagedGrid_2_neg);
Summary.Incongruent_AC_to_PFC_std     = std(averagedGrid_2_neg);
Summary.Incongruent_AC_to_PFC_stderr  = Summary.Incongruent_AC_to_PFC_std / sqrt(length(averagedGrid_2_neg));

%% Version 3 

% function PosNegHistogramOutline(averagedGrid_1_pos, averagedGrid_1_neg, averagedGrid_2_pos, averagedGrid_2_neg, nBins)
% 
%     % Ensure all inputs are column vectors
%     averagedGrid_1_pos = averagedGrid_1_pos(:);
%     averagedGrid_1_neg = averagedGrid_1_neg(:);
%     averagedGrid_2_pos = averagedGrid_2_pos(:);
%     averagedGrid_2_neg = averagedGrid_2_neg(:);
% 
%     % Define shared bin edges
%     [counts, bins] = hist([averagedGrid_1_pos; abs(averagedGrid_1_neg); averagedGrid_2_pos; abs(averagedGrid_2_neg)], nBins);
% 
%     % Compute histograms
%     posCounts_1 = hist(averagedGrid_1_pos, bins) / sum(counts) * 100;
%     negCounts_1 = hist(abs(averagedGrid_1_neg), bins) / sum(counts) * 100;
%     posCounts_2 = hist(averagedGrid_2_pos, bins) / sum(counts) * 100;
%     negCounts_2 = hist(abs(averagedGrid_2_neg), bins) / sum(counts) * 100;
% 
%     % Create figure
%     figure;
%     hold on;
% 
%     % Plot histogram outlines using stairs() to create step-like effect
%     stairs(bins, posCounts_1, '-r', 'LineWidth', 2); % Condition 1 Positive (Red)
%     stairs(bins, posCounts_2, '-b', 'LineWidth', 2); % Condition 2 Positive (Blue)
%     stairs(-bins, -negCounts_1, '-', 'Color', [1, 0.6, 0.6], 'LineWidth', 2); % Condition 1 Negative (Soft Red)
%     stairs(-bins, -negCounts_2, '-', 'Color', [0.6, 0.6, 1], 'LineWidth', 2); % Condition 2 Negative (Soft Blue)
% 
%     % Labels and legend
%     xlabel('PSI Value');
%     ylabel('Proportion of Sites');
%     legend({'Cond 1 Pos', 'Cond 2 Pos', 'Cond 1 Neg', 'Cond 2 Neg'}, 'Location', 'NorthEast');
% 
%     % Improve aesthetics
%     set(gca, 'FontSize', 12, 'LineWidth', 1.5);
%     box on;
%     grid on;
% 
%     hold off;
% end
% 
% 
