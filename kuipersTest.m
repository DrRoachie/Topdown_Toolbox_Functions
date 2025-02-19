function [V, pValue] = kuipersTest(angles)
    % Kuiper's Test for Uniformity
    %
    % INPUT:
    %   angles - vector of angles in degrees (0 to 360)
    %
    % OUTPUT:
    %   V      - Kuiper's test statistic
    %   pValue - p-value for the test (approximation for large datasets)

    % Convert angles to radians
    angles = deg2rad(angles);

    % Sort angles
    sortedAngles = sort(angles);

    % Normalize to [0, 1] by dividing by 2*pi
    normAngles = sortedAngles / (2 * pi);

    % Compute the empirical CDF
    n = length(normAngles);
    empiricalCDF = (1:n)' / n;

    % Compute the differences
    D_plus = max(empiricalCDF - normAngles);    % Maximum positive deviation
    D_minus = max(normAngles - [0; empiricalCDF(1:end-1)]); % Maximum negative deviation

    % Compute Kuiper's test statistic
    V = D_plus + D_minus;

    % Compute p-value (approximation for large n)
    lambda = (sqrt(n) + 0.155 + 0.24 / sqrt(n)) * V;
    pValue = 2 * sum((-1).^(1:10) .* exp(-2 * (1:10).^2 * lambda^2));

    % Ensure pValue is in range [0, 1] (numerical stability)
    pValue = max(min(pValue, 1), 0);
end


% 
% fprintf('Kuiper''s V Statistic: %.4f\n', V);
% fprintf('p-Value: %.4f\n', pValue);
% 
% if pValue < 0.05
%     fprintf('Reject the null hypothesis: data is not uniform.\n');
% else
%     fprintf('Fail to reject the null hypothesis: data is uniform.\n');
% end
