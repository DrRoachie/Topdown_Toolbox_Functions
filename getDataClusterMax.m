function [data_coh_Z, data_z_statistic_threshold, all_cluster_sums, all_clusters, data_max_cluster_sum] = getDataClusterMax(data_zstatistic,threshold_value)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% isolate the relevant bands of the data-derived zstastics

data_coh_Z.theta     = abs(data_zstatistic(:, 5:8));    % 4-8   Hz 
data_coh_Z.alpha     = abs(data_zstatistic(:, 8:14));   % 8-14  Hz
data_coh_Z.beta      = abs(data_zstatistic(:, 15:30));  % 15-30 Hz
data_coh_Z.gamma     = abs(data_zstatistic(:, 31:54));  % 31 - 50 Hz
data_coh_Z.highGamma = abs(data_zstatistic(:, 66:100)); % 70 - 100 Hz
data_coh_Z.all       = abs(data_zstatistic);            % 1 - 100 Hz

% determine the threshold value for each band
data_z_statistic_threshold.theta     = quantile(data_coh_Z.theta, threshold_value, 2);
data_z_statistic_threshold.alpha     = quantile(data_coh_Z.alpha, threshold_value, 2);
data_z_statistic_threshold.beta      = quantile(data_coh_Z.beta, threshold_value, 2);
data_z_statistic_threshold.gamma     = quantile(data_coh_Z.gamma, threshold_value, 2);
data_z_statistic_threshold.highGamma = quantile(data_coh_Z.highGamma, threshold_value, 2);
data_z_statistic_threshold.all       = quantile(data_coh_Z.all, threshold_value, 2); 

% threshold the data converting values that do not survive threshold to zero then find the zmax cluster for the data

thresholded_data.theta = data_coh_Z.theta;
thresholded_data.theta(thresholded_data.theta < data_z_statistic_threshold.theta) = 0;
[all_cluster_sums.theta , all_clusters.theta, data_max_cluster_sum.theta] = find_and_sum_clusters(thresholded_data.theta);

thresholded_data.alpha = data_coh_Z.alpha;
thresholded_data.alpha(thresholded_data.alpha < data_z_statistic_threshold.alpha) = 0;
[all_cluster_sums.alpha , all_clusters.alpha, data_max_cluster_sum.alpha] = find_and_sum_clusters(thresholded_data.alpha);

thresholded_data.beta = data_coh_Z.beta;
thresholded_data.beta(thresholded_data.beta < data_z_statistic_threshold.beta) = 0;
[all_cluster_sums.beta , all_clusters.beta, data_max_cluster_sum.beta] = find_and_sum_clusters(thresholded_data.beta);

thresholded_data.gamma = data_coh_Z.gamma;
thresholded_data.gamma(thresholded_data.gamma < data_z_statistic_threshold.gamma) = 0;
[all_cluster_sums.gamma , all_clusters.gamma, data_max_cluster_sum.gamma] = find_and_sum_clusters(thresholded_data.gamma);

thresholded_data.highGamma = data_coh_Z.highGamma;
thresholded_data.highGamma(thresholded_data.highGamma < data_z_statistic_threshold.highGamma) = 0;
[all_cluster_sums.highGamma , all_clusters.highGamma, data_max_cluster_sum.highGamma] = find_and_sum_clusters(thresholded_data.highGamma);

thresholded_data.all = data_coh_Z.all;
thresholded_data.all(thresholded_data.all < data_z_statistic_threshold.all) = 0;
[all_cluster_sums.all , all_clusters.all, data_max_cluster_sum.all] = find_and_sum_clusters(thresholded_data.all);


end