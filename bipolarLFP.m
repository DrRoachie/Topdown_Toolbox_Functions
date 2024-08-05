function [bipolar_lfp,label] = bipolarLFP(lfp,eInfo)
%bipolarLFP Summary of this function goes here
%   lfp   ... lfp data (channel x sample)
%   eInfo ... electrode info obtained by countChannels.m
nCh_all = size(lfp,1);
nElectrode = eInfo.nElectrode;

bipolar_lfp = [];
label = {};
for j=1:nElectrode
    nCh_elec = eInfo.nChannel(j);
    temp_list = eInfo.list{j};
    eLFP = lfp(1:nCh_elec,:); % LFP in each electrode
    if nCh_elec==24 % 24ch v-probe
        for i=1:20
            bipolar_lfp_a(i,:) = eLFP(i,:) - eLFP(i+4,:);
        end
        list_a = temp_list(3:22);
    elseif nCh_elec==16 % 16ch v-probe
        for i=1:14
            bipolar_lfp_a(i,:) = eLFP(i,:) - eLFP(i+2,:);
        end
        list_a = temp_list(2:15);
    end
    bipolar_lfp = cat(1,bipolar_lfp,bipolar_lfp_a);
    label = [label; list_a];
    
    lfp(1:nCh_elec,:) = [];
    clear eLFP bipolar_lfp_a list_a
end

if ~isempty(lfp)
    error('something must be wrong...');
end

end

