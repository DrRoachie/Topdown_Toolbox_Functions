
function [Output_Array] = blurArray(Input_Array, scale_factor)

% decreases the array resolution by the scale factor by averaging the 
% values in a sliding window of size scale_factor x scale_factor 

% scale_factor: scale_factor determines the size of the sliding window, so 
% the output array will be the size of the input array minus (scale_factor - 1). 
% For example, a scale_factor of 2 will decrease a 20x20 array to 19x19.

for row = 1:length(Input_Array(:,1)) - (scale_factor - 1)
    for col = 1:length(Input_Array(:,1)) - (scale_factor - 1)
        % fill in sliding window with values from input array
        window = [];
        for i = 0:scale_factor - 1
            for j = 0:scale_factor - 1
                window(i + 1, j + 1) = Input_Array(row + i, col + j);
            end
        end

        Output_Array(row,col) = mean(window, 'all');
    end
end

% % plot input and output arrays
% figure;
% subplot(1,2,1);
% h1 = heatmap(Input_Array);
% c_lim = get(h1, 'ColorLimits');
% 
% subplot(1,2,2);
% h2 = heatmap(Output_Array);
% set(h2, 'ColorLimits', c_lim);
