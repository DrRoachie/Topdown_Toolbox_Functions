function [Monte_coh_Z,permutation_zthresholds,thresholded_Permutation, monte_all_cluster_sums, monte_all_clusters, monte_max_cluster_sum] = getPermClusterMax(Monte_Carlo_Array,threshold_value)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Isolate Theta, Alpha, Beta, Gamma, and High-Gamma from the Monte Carlo Array 

Monte_coh_Z.theta     = abs(Monte_Carlo_Array(:, 5:8));    % 4-8   Hz 
Monte_coh_Z.alpha     = abs(Monte_Carlo_Array(:, 8:14));   % 9-14  Hz
Monte_coh_Z.beta      = abs(Monte_Carlo_Array(:, 15:30));  % 15-30 Hz
Monte_coh_Z.gamma     = abs(Monte_Carlo_Array(:, 31:54));  % 31 - 50 Hz
Monte_coh_Z.highGamma = abs(Monte_Carlo_Array(:, 66:100)); % 70 - 100 Hz
Monte_coh_Z.all       = abs(Monte_Carlo_Array);            % 1 - 100 Hz


%% For each band and each iteration calculate the threshold value 

permutation_zthresholds.theta      = quantile(Monte_coh_Z.theta, threshold_value, 2);
permutation_zthresholds.alpha      = quantile(Monte_coh_Z.alpha, threshold_value, 2);
permutation_zthresholds.beta       = quantile(Monte_coh_Z.beta, threshold_value, 2);
permutation_zthresholds.gamma      = quantile(Monte_coh_Z.gamma, threshold_value, 2);
permutation_zthresholds.highGamma  = quantile(Monte_coh_Z.highGamma, threshold_value, 2);
permutation_zthresholds.all        = quantile(Monte_coh_Z.all, threshold_value, 2);

%% threshold the permutation 

for t = 1:size(Monte_coh_Z.theta, 1)
    belowThresholdIndices = Monte_coh_Z.theta(t,:) < permutation_zthresholds.theta(t);
    thresholded_Permutation.theta(t, :) = Monte_coh_Z.theta(t, :);
    thresholded_Permutation.theta(t, belowThresholdIndices) = 0; % sets every value below the 95 quantile to 0;
end

for t = 1:size(Monte_coh_Z.alpha, 1)
    belowThresholdIndices = Monte_coh_Z.alpha(t,:) < permutation_zthresholds.alpha(t);
    thresholded_Permutation.alpha(t, :) = Monte_coh_Z.alpha(t, :);
    thresholded_Permutation.alpha(t, belowThresholdIndices) = 0; % sets every value below the 95 quantile to 0;
end

for t = 1:size(Monte_coh_Z.beta, 1)
    belowThresholdIndices = Monte_coh_Z.beta(t,:) < permutation_zthresholds.beta(t);
    thresholded_Permutation.beta(t, :) = Monte_coh_Z.beta(t, :);
    thresholded_Permutation.beta(t, belowThresholdIndices) = 0; % sets every value below the 95 quantile to 0;
end

for t = 1:size(Monte_coh_Z.gamma, 1)
    belowThresholdIndices = Monte_coh_Z.gamma(t,:) < permutation_zthresholds.gamma(t);
    thresholded_Permutation.gamma(t, :) = Monte_coh_Z.gamma(t, :);
    thresholded_Permutation.gamma(t, belowThresholdIndices) = 0; % sets every value below the 95 quantile to 0;
end

for t = 1:size(Monte_coh_Z.highGamma, 1)
    belowThresholdIndices = Monte_coh_Z.highGamma(t,:) < permutation_zthresholds.highGamma(t);
    thresholded_Permutation.highGamma(t, :) = Monte_coh_Z.highGamma(t, :);
    thresholded_Permutation.highGamma(t, belowThresholdIndices) = 0; % sets every value below the 95 quantile to 0;
end

for t = 1:size(Monte_coh_Z.all, 1)
    belowThresholdIndices = Monte_coh_Z.all(t,:) < permutation_zthresholds.all(t);
    thresholded_Permutation.all(t, :) = Monte_coh_Z.all(t, :);
    thresholded_Permutation.all(t, belowThresholdIndices) = 0; % sets every value below the 95 quantile to 0;
end
%% generate a distribution of z-max values from the thresholded permutation 

[monte_all_cluster_sums.theta , monte_all_clusters.theta, monte_max_cluster_sum.theta] = find_and_sum_clusters(thresholded_Permutation.theta);
[monte_all_cluster_sums.alpha , monte_all_clusters.alpha, monte_max_cluster_sum.alpha] = find_and_sum_clusters(thresholded_Permutation.alpha);
[monte_all_cluster_sums.beta , monte_all_clusters.beta, monte_max_cluster_sum.beta] = find_and_sum_clusters(thresholded_Permutation.beta);
[monte_all_cluster_sums.gamma , monte_all_clusters.gamma, monte_max_cluster_sum.gamma] = find_and_sum_clusters(thresholded_Permutation.gamma);
[monte_all_cluster_sums.highGamma , monte_all_clusters.highGamma, monte_max_cluster_sum.highGamma] = find_and_sum_clusters(thresholded_Permutation.highGamma);
[monte_all_cluster_sums.all , monte_all_clusters.all, monte_max_cluster_sum.all] = find_and_sum_clusters(thresholded_Permutation.all);
end