function [eInfo] = countChannels(chanList)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

LIST_ELECTRODE = {'D1','D2','D3','D4'};

counter = 0; % count number of electrode
for i=1:numel(LIST_ELECTRODE)
    str = LIST_ELECTRODE{i};
    temp = chanList(contains(chanList,str));
    if ~isempty(temp)
        counter = counter + 1;
        for j=1:numel(temp)
            short_name = temp{j};
            str_start = strfind(short_name,str);
            short_name = short_name(str_start:end);
            S{j,1} = short_name;
        end
        chElec{counter} = S;
        clear S
    end
end

nElectrode = counter;
for i=1:nElectrode
    nCh(i) = numel(chElec{i}); 
end

% electrode info
eInfo.list       = chElec;
eInfo.nChannel   = nCh;
eInfo.nElectrode = nElectrode;

end

