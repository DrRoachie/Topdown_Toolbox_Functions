function [all_cluster_sums, all_clusters, max_cluster_sum] = find_and_sum_clusters(matrix)

% input: a matrix of thresholded statistic test values 
% output: all_cluster_sums is the sum of each cluster 
% output: all_clusters is the indicies of clusters of adjacent non-zero values 
% output: max_cluster_sum is largest cluster sum within the band

    [rows, ~] = size(matrix);
    all_clusters = cell(rows, 1);
    all_cluster_sums = cell(rows, 1);
    max_cluster_sum = zeros(rows, 1);

    for i = 1:rows
        row = matrix(i, :);
        non_zero_indices = find(row);
        
        if ~isempty(non_zero_indices)
            clusters = split_clusters(non_zero_indices);
            all_clusters{i} = clusters;
            
            all_cluster_sums{i} = zeros(1, length(clusters));

            for j = 1:length(clusters)
                all_cluster_sums{i}(j) = sum(row(clusters{j}));
            end

            max_cluster_sum(i) = max(all_cluster_sums{i});

            % disp(['Row ', num2str(i), ' Clusters: ', mat2str(clusters)]);
            % disp(['Row ', num2str(i), ' Cluster Sums: ', mat2str(all_cluster_sums{i})]);
            % disp(['Row ', num2str(i), ' Max Cluster Sum: ', num2str(max_cluster_sum(i))]);
        end
    end
end

function clusters = split_clusters(indices)
    clusters = {};
    current_cluster = [indices(1)];

    for i = 2:length(indices)
        if indices(i) == indices(i-1) + 1
            current_cluster = [current_cluster, indices(i)];
        else
            clusters{end+1} = current_cluster;
            current_cluster = [indices(i)];
        end
    end

    clusters{end+1} = current_cluster;
end
