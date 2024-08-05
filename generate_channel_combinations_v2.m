
%this function makes a cfg.channelcmb structure fomr ChanWise_SpecEval
%script

    % in AC I have channels 03 - 20 and PFC I have channels 03 -20.  My analysis cares about five channel groupings. 
    % 
    % PFClowmid_ACupmid: PFC (channels 06-10) and AC (channels 11-15)
    % PFCupmid_ACupmid: PFC (channels 06-10) and AC (channels 06-10)
    % PFCdeep_ACupmid: PFC (channels 16-20) and AC (channels 06-10)
    % PFCdeep_ACdeep: PFC (channels 16-20) and AC (channels 16-20)
    % PFCupper_ACupper: PFC (channels 03-07) and AC (channels 03-07)
    % 
    % for each channel grouping, i want to subselect channels listed within the significant_channel_names cell array, and construct cfg.channelcmb from there

% Define channel ranges for each grouping

% Define channel ranges for each grouping
channel_ranges = {
    'PFClowmid_ACupmid', 6:10, 11:15;
    'PFCupmid_ACupmid', 6:10, 6:10;
    'PFCdeep_ACupmid', 16:20, 6:10;
    'PFCdeep_ACdeep', 16:20, 16:20;
    'PFCupper_ACupper', 3:7, 3:7
};

% Initialize cell arrays for each channel grouping
PFClowmid_ACupmid = {};
PFCupmid_ACupmid = {};
PFCdeep_ACupmid = {};
PFCdeep_ACdeep = {};
PFCupper_ACupper = {};

% Iterate through each channel grouping
for i = 1:size(channel_ranges, 1)
    grouping = channel_ranges{i, 1};
    PFC_range = channel_ranges{i, 2};
    AC_range = channel_ranges{i, 3};
    
    % Construct channel names for the current grouping
    PFC_channels = arrayfun(@(x) sprintf('*PFC_ch%02d', x), PFC_range, 'UniformOutput', false);
    AC_channels = arrayfun(@(x) sprintf('*AC_ch%02d', x), AC_range, 'UniformOutput', false);
    
    % Filter channels based on significant_channel_names
    PFC_channels = PFC_channels(ismember(PFC_channels, significant_channels_c));
    AC_channels = AC_channels(ismember(AC_channels, significant_channels_c));
    
    % Create combinations of significant channels for the current grouping
    channelcmb = {};
    for ac = 1:length(AC_channels)
        for pfc = 1:length(PFC_channels)
            channelcmb{end+1, 1} = AC_channels{ac};
            channelcmb{end, 2} = PFC_channels{pfc};
        end
    end
    
    % Assign the combinations to the corresponding variable
    switch grouping
        case 'PFClowmid_ACupmid'
            PFClowmid_ACupmid = channelcmb;
        case 'PFCupmid_ACupmid'
            PFCupmid_ACupmid = channelcmb;
        case 'PFCdeep_ACupmid'
            PFCdeep_ACupmid = channelcmb;
        case 'PFCdeep_ACdeep'
            PFCdeep_ACdeep = channelcmb;
        case 'PFCupper_ACupper'
            PFCupper_ACupper = channelcmb;
    end
end
