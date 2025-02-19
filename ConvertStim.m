
function V = ConvertStim(V, RecDate)


%Initialize the category array
V.target = strings(size(V.stim));

if strcmp(RecDate,'190330') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 395; 396, 1795; 1796, 2279; 2280, 2369]; % Define row ranges
high_values = [1903, 2263, 5382, 9051];    % Corresponding high values
low_values = [141, 141, 951, 2263];       % Corresponding low values
end 

if strcmp(RecDate,'190404') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1522; 1523, 1956; 1957, 2050]; % Define row ranges
high_values = [12800, 6400, 1600];    % Corresponding high values
low_values = [1600, 673, 200];       % Corresponding low values
end 

if strcmp(RecDate,'190413') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1393; 1394, 1827; 1828, 2276]; % Define row ranges
high_values = [2263, 5382, 12800];    % Corresponding high values
low_values = [400, 673, 2263];       % Corresponding low values
end 

if strcmp(RecDate,'190414') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 398; 399, 1529; 1530, 1721; 1722, 2138]; % Define row ranges
high_values = [2263, 3200, 6400, 12800];    % Corresponding high values
low_values = [100, 200, 951, 3200];       % Corresponding low values
end 

if strcmp(RecDate,'190416') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1497; 1498, 1982; 1983, 2292]; % Define row ranges
high_values = [3200, 6400, 12800];    % Corresponding high values
low_values = [141, 1345, 3200];       % Corresponding low values
end 

if strcmp(RecDate,'190418') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1459; 1460, 1876; 1877, 2506]; % Define row ranges
high_values = [3200, 5382, 12800];    % Corresponding high values
low_values = [141, 566, 3200];       % Corresponding low values
end 

if strcmp(RecDate,'190419') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1630; 1631, 2096; 2097, 2658]; % Define row ranges
high_values = [2691, 5383, 12800];    % Corresponding high values
low_values = [336, 673, 2691];       % Corresponding low values
end 

if strcmp(RecDate,'190421') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1452; 1453, 1888; 1889, 2554]; % Define row ranges
high_values = [2263, 7611, 10763];    % Corresponding high values
low_values = [336, 400, 2263];       % Corresponding low values
end 

if strcmp(RecDate,'190423') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1519; 1520, 1975; 1976, 2458]; % Define row ranges
high_values = [2691, 5382, 12800];    % Corresponding high values
low_values = [336, 800, 2691];       % Corresponding low values
end 

if strcmp(RecDate,'190429') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1398]; % Define row ranges
high_values = 1903;    % Corresponding high values
low_values = 238;       % Corresponding low values
end 

if strcmp(RecDate,'190515') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1475; 1476, 1916; 1917, 2851]; % Define row ranges
high_values = [3200, 5382, 12800];    % Corresponding high values
low_values = [119, 673, 3200];       % Corresponding low values
end 

if strcmp(RecDate,'190517') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 586; 587, 1134; 1135, 2045; 2046,2285; 2286, 3206]; % Define row ranges
high_values = [3200, 7611, 12800, 7611, 3200];    % Corresponding high values
low_values = [100, 673, 3200, 673, 100];       % Corresponding low values
end 

if strcmp(RecDate,'190531') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 584; 585, 1113; 1114, 2232; 2233, 2641; 2642, 3225]; % Define row ranges
high_values = [3200, 6400, 12800, 6400, 3200];    % Corresponding high values
low_values = [238, 673, 3200, 673, 238];       % Corresponding low values
end 

if strcmp(RecDate,'190603') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 591; 592, 1128; 1129, 2187; 2188, 2612; 2613, 3531; 3532, 3981; 3982, 4191]; % Define row ranges
high_values = [2263, 6400, 12800, 6400, 2263, 6400, 12800];    % Corresponding high values
low_values = [141, 800, 2263, 800, 141, 800, 2263];       % Corresponding low values
end 

if strcmp(RecDate,'190605') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 819; 820, 1571; 1572, 2767; 2768, 3209; 3210, 4084]; % Define row ranges
high_values = [15222, 5382, 2691, 5382, 15222];    % Corresponding high values
low_values = [2691, 800, 119, 800, 2691];       % Corresponding low values
end 

if strcmp(RecDate,'190703') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 779; 780, 1572; 1573, 2870; 2871, 3199; 3200, 3294]; % Define row ranges
high_values = [10763, 6400, 1600, 6400, 10763];    % Corresponding high values
low_values = [1600, 476, 119, 476, 1600];       % Corresponding low values
end 

if strcmp(RecDate,'190711') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 573; 574, 1111; 1112, 2305; 2306, 2751; 2752, 3256]; % Define row ranges
high_values = [2263, 5382, 12800, 5382, 2263];    % Corresponding high values
low_values = [100, 566, 2263, 566, 100];       % Corresponding low values
end 

if strcmp(RecDate,'190713') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 801; 802, 1486; 1487, 2681; 2682, 2996; 2997, 3462]; % Define row ranges
high_values = [2691, 6400, 12800, 6400, 2691];    % Corresponding high values
low_values = [119, 339, 2691, 336, 119];       % Corresponding low values
end 

if strcmp(RecDate,'190718') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 571; 572, 1102; 1103, 2031; 2032, 2313; 2314, 2285]; % Define row ranges
high_values = [2263, 4525, 12800, 4525, 2263];    % Corresponding high values
low_values = [336, 800, 2263, 800, 336];       % Corresponding low values
end 

if strcmp(RecDate,'190720') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 594; 595, 1239; 1240, 2045; 2046, 2336; 2337, 2789; 2790, 2980; 2981, 3149]; % Define row ranges
high_values = [10763, 6400, 2263, 6400, 10763, 6400, 2263];    % Corresponding high values
low_values = [2263, 1600, 336, 1600, 2263, 1600, 336];       % Corresponding low values
end 

if strcmp(RecDate,'190723') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 597; 598, 1149; 1150, 2009; 2010, 2278; 2279, 2849; 2850, 3050]; % Define row ranges
high_values = [2263, 4525, 12800, 4525, 2263, 4525];    % Corresponding high values
low_values = [200, 400, 2263, 400, 200, 400];       % Corresponding low values
end 

if strcmp(RecDate,'190725') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 567; 568, 1108; 1109, 1993; 1994, 2262; 2263, 2631]; % Define row ranges
high_values = [2263, 4525, 12800, 4525, 2263];    % Corresponding high values
low_values = [400, 800, 2263, 800, 400];       % Corresponding low values
end 

if strcmp(RecDate,'190417') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1604; 1605, 2082; 2083, 3122; 3123, 3577; 3578, 3855]; % Define row ranges
high_values = [3200, 6400, 10763, 6400, 3200];    % Corresponding high values
low_values = [200, 1345, 3200, 1345, 200];       % Corresponding low values
end 

if strcmp(RecDate,'190422') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1009; 1010, 1998; 1999, 3483; 3484, 3933; 3944, 4443]; % Define row ranges
high_values = [3805, 6400, 12800, 6400, 3805];    % Corresponding high values
low_values = [336, 951, 3805, 951, 336];       % Corresponding low values
end 

if strcmp(RecDate,'190425') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 974; 975, 2017; 2018, 3692; 3693, 4108; 4109, 4851]; % Define row ranges
high_values = [2691, 6400, 10763, 6400, 2691];    % Corresponding high values
low_values = [200, 800, 2691, 800, 200];       % Corresponding low values
end 

if strcmp(RecDate,'190425') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 974; 975, 2017; 2018, 3692; 3693, 4108; 4109, 4851]; % Define row ranges
high_values = [2691, 6400, 10763, 6400, 2691];    % Corresponding high values
low_values = [200, 800, 2691, 800, 200];       % Corresponding low values
end 

if strcmp(RecDate,'190427') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1043; 1044, 2073; 2074, 3551; 3552, 4021; 4022, 4820]; % Define row ranges
high_values = [2691, 6400, 12800, 6400, 2691];    % Corresponding high values
low_values = [673, 951, 2691, 951, 673];       % Corresponding low values
end 

if strcmp(RecDate,'190427') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1043; 1044, 2073; 2074, 3551; 3552, 4021; 4022, 4820]; % Define row ranges
high_values = [2691, 6400, 12800, 6400, 2691];    % Corresponding high values
low_values = [673, 951, 2691, 951, 673];       % Corresponding low values
end 

if strcmp(RecDate,'190502') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 999; 1000, 2012; 2013, 3656; 3657, 4135; 4136, 4651]; % Define row ranges
high_values = [2263, 6400, 9051, 6400, 2263];    % Corresponding high values
low_values = [238, 1345, 2263, 1345, 238];       % Corresponding low values
end 

if strcmp(RecDate,'190514') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 40; 41, 1121; 1122, 2010; 2011, 3556; 3557, 4056; 4057, 4522]; % Define row ranges
high_values = [2263, 2263, 7611, 12800, 7611, 2263];    % Corresponding high values
low_values = [336, 119, 673, 2263, 673, 119];       % Corresponding low values
end 

if strcmp(RecDate,'190516') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 932; 933, 1815; 1816, 3267; 3268, 3715; 3716, 4585; 4586, 4999; 5000, 5452]; % Define row ranges
high_values = [2263, 5382, 12800, 5382, 2263, 5382, 12800];    % Corresponding high values
low_values = [336, 673, 2263, 673, 336, 673, 2263];       % Corresponding low values
end 

if strcmp(RecDate,'190525') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 951; 952, 1843; 1844, 3287; 3288, 3664; 3665, 4542; 4543, 4965; 4966, 5024]; % Define row ranges
high_values = [2263, 6400, 12800, 6400, 2263, 6400, 12800];    % Corresponding high values
low_values = [336, 800, 2263, 800, 336, 800, 2263];       % Corresponding low values
end 

if strcmp(RecDate,'190527') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 936; 937, 1915; 1916, 3360; 3361, 3805; 3806, 4706; 4707, 5178; 5179, 5631]; % Define row ranges
high_values = [2263, 5382, 12800, 5382, 2263, 5383, 12800];    % Corresponding high values
low_values = [400, 1345, 2263, 1345, 400, 1345, 2263];       % Corresponding low values
end 

if strcmp(RecDate,'190530') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1017; 1018, 1923; 1924, 3353; 3354, 3751; 3752, 4685; 4686, 5008]; % Define row ranges
high_values = [2263, 6400, 12800, 6400, 2263, 6400];    % Corresponding high values
low_values = [400, 800, 2263, 800, 400, 800];       % Corresponding low values
end 

if strcmp(RecDate,'190601') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1000; 1001, 1904; 1905, 3592; 3593, 3979; 3980, 4893; 4894, 4980]; % Define row ranges
high_values = [2691, 6400, 9051, 6400, 2691, 6400];    % Corresponding high values
low_values = [400, 800, 2691, 800, 400, 800];       % Corresponding low values
end 

if strcmp(RecDate,'190604') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 979; 980, 2010; 2011, 3530; 3531, 4046; 4047, 4981; 4982, 5435]; % Define row ranges
high_values = [12800, 7611, 2263, 7611, 12800, 7611];    % Corresponding high values
low_values = [2263, 1600, 238, 1600, 2263, 1600];       % Corresponding low values
end 

if strcmp(RecDate,'190704') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 975; 976, 1852; 1853, 3305; 3306, 3704; 3705, 4234]; % Define row ranges
high_values = [12800, 6400, 1600, 6400, 12800];    % Corresponding high values
low_values = [1600, 566, 238, 566, 1600];       % Corresponding low values
end 

if strcmp(RecDate,'190709') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1000; 1001, 1894; 1895, 3312; 3513, 3923; 3924, 4420]; % Define row ranges
high_values = [3200, 5382, 7611, 5382, 3200];    % Corresponding high values
low_values = [100, 476, 3200, 476, 100];         % Corresponding low values
end 

if strcmp(RecDate,'190717') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 977; 978, 1827; 1828, 3351; 3352, 3772; 3773, 4247]; % Define row ranges
high_values = [10763, 6400, 2263, 6400, 10763];    % Corresponding high values
low_values = [2263, 400, 100, 400, 2263];         % Corresponding low values
end 

if strcmp(RecDate,'190719') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1080; 1081, 2055; 2056, 3481; 3482, 3899; 3900, 4401]; % Define row ranges
high_values = [10763, 4525, 1903, 4525, 10763];    % Corresponding high values
low_values = [1903, 800, 283, 800, 1903];         % Corresponding low values
end 

if strcmp(RecDate,'190722') == 1
% Define the block ranges and corresponding high/low values
ranges = [1, 1045; 1046, 1881; 1882, 3280; 3281, 3703; 3704, 4154]; % Define row ranges
high_values = [2691, 4525, 12800, 4525, 2691];    % Corresponding high values
low_values = [119, 566, 2691, 566, 199];         % Corresponding low values
end

% Loop through each block range
for block = 1:size(ranges, 1)
    % Extract start and end indices for the block
    start_idx = ranges(block, 1);
    end_idx = ranges(block, 2);
    
    % Get the high and low values for this block
    high_value = high_values(block);
    low_value = low_values(block);
    
    % Assign 'H' or 'L' based on the block definitions
    for i = start_idx:end_idx
        if V.stim(i) == high_value
            V.target(i) = 'H';
        elseif V.stim(i) == low_value
            V.target(i) = 'L';
        else
            V.target(i) = ''; % Leave blank if not matching
        end
    end
end
