function index = findChannel(label,Area,ch)
%findCannel Summary of this function goes here
%   label -- list of the channel
%   Area  -- brain area (either 'AC' or 'PFC')
%   ch    -- channel number to pick up
if ch<10
    string_ch = ['ch0' num2str(ch)];
else
    string_ch = ['ch' num2str(ch)];
end

search_string = strcat(Area,'_',string_ch);

index = contains(label,search_string);

end

