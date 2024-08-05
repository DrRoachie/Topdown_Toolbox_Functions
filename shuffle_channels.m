function [shuffled_trial] = shuffle_channels(ordered_array);

ordered_array(randperm(size(ordered_array, 1)), :);

end