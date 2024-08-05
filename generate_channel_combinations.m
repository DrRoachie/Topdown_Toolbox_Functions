function channel_combinations = generate_channel_combinations(ac_start, ac_end, pfc_start, pfc_end)
    %this function makes a cfg.channelcmb structure when specifing a channel range. 


    ac_channels = arrayfun(@(x) sprintf('*AC_ch%02d', x), ac_start:ac_end, 'UniformOutput', false);
    pfc_channels = arrayfun(@(x) sprintf('*PFC_ch%02d', x), pfc_start:pfc_end, 'UniformOutput', false);
    
    channel_combinations = cell(numel(ac_channels) * numel(pfc_channels), 2);
    idx = 1;
    for i = 1:numel(ac_channels)
        for j = 1:numel(pfc_channels)
            channel_combinations{idx, 1} = ac_channels{i};
            channel_combinations{idx, 2} = pfc_channels{j};
            idx = idx + 1;
        end
    end
end