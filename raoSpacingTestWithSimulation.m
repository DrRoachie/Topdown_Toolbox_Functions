function [U, pValue] = raoSpacingTestWithSimulation(angles, numSimulations)
    % Rao's Spacing Test with Monte Carlo Simulation for Critical Values
    %
    % INPUT:
    %   angles - vector of angles in degrees (0 to 360)
    %   numSimulations - number of simulations to estimate critical values
    %
    % OUTPUT:
    %   U - Rao's Spacing Test statistic
    %   pValue - estimated p-value

    % Convert angles to radians
    angles = deg2rad(angles);

    % Sort angles
    sortedAngles = sort(angles);

    % Add circular wrap-around
    wrappedAngles = [sortedAngles; sortedAngles(1) + 2 * pi];

    % Compute spacings
    spacings = diff(wrappedAngles);

    % Compute Rao's U statistic
    n = length(angles);
    U = (1 / pi) * sum(abs(spacings - (2 * pi / n)));

    % Monte Carlo simulation to estimate p-value
    simulatedU = zeros(numSimulations, 1);
    for i = 1:numSimulations
        % Generate uniform random angles
        randomAngles = sort(2 * pi * rand(n, 1));
        randomSpacings = diff([randomAngles; randomAngles(1) + 2 * pi]);
        simulatedU(i) = (1 / pi) * sum(abs(randomSpacings - (2 * pi / n)));
    end

    % Compute p-value
    pValue = mean(simulatedU >= U);
end


