function [labels_a_out, labels_b_out] = SpecBalancedRangeBased(theta_a, theta_b, labels_a, labels_b)
        % inputs:
        %  - theta_a: a vector of theta power values for region A 
        %  - theta_b: a vector of theta power values for region B
        %  - labels_a: corresponding full list of labels for region A
        %  - labels_b: corresponding full list of labels for region B
    
        % outputs:
        % - labels_a_out: selected labels from the region based on overlap or top 25%       
        % - labels_b_out: selected labels from the region based on overlap or bottom 25%       
        %   power range (bottom 25%) or overlapping range
    
        % Define the overlapping range
        min_overlap = max([min(theta_a), min(theta_b)]);
        max_overlap = min([max(theta_a), max(theta_b)]);
    
        % Check if they overlap
        does_overlap = min_overlap <= max_overlap;
    
        % Initialize empty output
        labels_a_out = {};
        labels_b_out = {};
    
        if does_overlap
            % === CASE 1: They OVERLAP → collect channels within the shared range ===
            idx_a = find(theta_a >= min_overlap & theta_a <= max_overlap);
            idx_b = find(theta_b >= min_overlap & theta_b <= max_overlap);
    
            labels_a_out = labels_a(idx_a);  % "small" and "large" are arbitrary here
            labels_b_out = labels_b(idx_b);
    
        else
            % === CASE 2: They DO NOT OVERLAP → do quartile-based selection ===
            range_a = max(theta_a) - min(theta_a);
            range_b = max(theta_b) - min(theta_b);
    
            if range_a < range_b   
                theta_small = theta_a;
                theta_large = theta_b;
                labels_small = labels_a;
                labels_large = labels_b;
            else
                theta_small = theta_b;
                theta_large = theta_a;
                labels_small = labels_b;
                labels_large = labels_a;
            end
    
            % Get top 25% from smaller, bottom 25% from larger
            pct = 0.25;
            n_small = length(theta_small);
            n_large = length(theta_large);
    
            n_top = round(pct * n_small);
            n_bot = round(pct * n_large);
    
            [~, idx_sorted_small] = sort(theta_small, 'ascend');
            [~, idx_sorted_large] = sort(theta_large, 'ascend');
    
            top_idx = idx_sorted_small(end - n_top + 1:end);
            bot_idx = idx_sorted_large(1:n_bot);
    
            labels_a_out = labels_small(top_idx);
            labels_b_out = labels_large(bot_idx);
            
        end
    end
    
