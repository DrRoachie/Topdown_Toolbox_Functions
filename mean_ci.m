function [m, ci, sdhat] = mean_ci(x, alpha)
    % Mean and 95% CI via classic t-interval (no bootstrap).
    if nargin < 2, alpha = 0.05; end
    x = x(:);
    x = x(~isnan(x));
    n = numel(x);
    if n == 0
        m = NaN; sdhat = NaN; ci = [NaN; NaN]; return;
    end
    m = mean(x);
    sdhat = std(x);
    se = sdhat / sqrt(n);
    tcrit = tinv(1 - alpha/2, max(n-1,1));
    ci = [m - tcrit*se; m + tcrit*se];
end