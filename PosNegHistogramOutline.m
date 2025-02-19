function PosNegHistogramOutline(averagedGrid_1_pos, averagedGrid_1_neg, averagedGrid_2_pos, averagedGrid_2_neg, num_bins)

    % Ensure all inputs are column vectors
    averagedGrid_1_pos = averagedGrid_1_pos(:);
    averagedGrid_1_neg = averagedGrid_1_neg(:);
    averagedGrid_2_pos = averagedGrid_2_pos(:);
    averagedGrid_2_neg = averagedGrid_2_neg(:);

    % Define bin edges based on data range
    bin_edges = linspace(min([averagedGrid_1_pos; averagedGrid_2_pos; abs(averagedGrid_1_neg); abs(averagedGrid_2_neg)]), ...
                         max([averagedGrid_1_pos; averagedGrid_2_pos; abs(averagedGrid_1_neg); abs(averagedGrid_2_neg)]), num_bins);

    % Compute histograms (normalized as probability density function)
    [counts1_pos, edges1] = histcounts(averagedGrid_1_pos, bin_edges, 'Normalization', 'pdf');
    [counts1_neg, edges1_neg] = histcounts(abs(averagedGrid_1_neg), bin_edges, 'Normalization', 'pdf');
    [counts2_pos, edges2] = histcounts(averagedGrid_2_pos, bin_edges, 'Normalization', 'pdf');
    [counts2_neg, edges2_neg] = histcounts(abs(averagedGrid_2_neg), bin_edges, 'Normalization', 'pdf');

    counts1_pos = counts1_pos / sum(counts1_pos) * 100;
    counts1_neg = counts1_neg / sum(counts1_neg) * 100;
    counts2_pos = counts2_pos / sum(counts2_pos) * 100;
    counts2_neg = counts2_neg / sum(counts2_neg) * 100;

    % Compute bin centers
    bin_centers1 = (edges1(1:end-1) + edges1(2:end)) / 2;
    bin_centers2 = (edges2(1:end-1) + edges2(2:end)) / 2;
    
    % Create figure
    figure;
    hold on;

    % Plot positive values for both conditions
    plot(bin_centers1, counts1_pos, '-r', 'LineWidth', 2.5); % Condition 1 Positive (Red)
    plot(bin_centers2, counts2_pos,'Color', [1, 0.6, 0.6], 'LineWidth', 2.5); % Condition 2 Positive (Blue)

    % Plot negative values as mirrored below x-axis
    plot(-bin_centers1, -counts1_neg, '-b', 'LineWidth', 2.5); % Soft Red for Condition 1 Negative
    plot(-bin_centers2, -counts2_neg, 'Color', [0.6, 0.6, 1], 'LineWidth', 2.5); % Soft Blue for Condition 2 Negative

    % Labels and legend
    xlabel('PSI Value');
    ylabel('Proportion of Sites');
    %ylim([-100 100]);
    legend({'Cong.   top-down', 'Incong. top-down', 'Cong.   bottom-up', 'Incong. bottom-up'}, 'Location', 'SouthEast');

    % Improve aesthetics
    set(gca, 'FontSize', 8, 'LineWidth', 1.5);
    box on;
    grid on;

    hold off;
end
