function [padj, order] = holm_bonferroni(p)
    p = p(:); m = numel(p);
    [ps, order] = sort(p, 'ascend');                  % sorted p
    adj = ps .* (m - (1:m)' + 1);                     % Holm scaling
    padj_sorted = min(cummax(adj), 1);                % monotone non-decreasing
    padj = zeros(m,1); padj(order) = padj_sorted;     % back to original order
end
