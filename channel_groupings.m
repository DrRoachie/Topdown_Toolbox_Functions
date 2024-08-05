%% Housekeeping 
% For each session and 

DATA_DIR    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_05_25_Analysis\00_Significant_chans\MrM\LED\190425';
SAVE_DIR    = 'D:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\zz_MetaData\2024_05_25_Analysis\01_Layer_Groupings\MrM\LED\190425';
addpath(genpath(DATA_DIR));

Animal           = 'MrM';                  % Options: 'MrCassius', 'MrM'; 
RecDate          = '190425';        
Epoch            = 'preCueOnset';                % Options: 'testToneOnset', 'preCueOnset', or 'moveOnset'
Behavior         = 'Correct';                    % Options:  'Correct', Wrong, due to memory constraints 

freq_band_list = {'theta'; 'alpha'; 'beta'; 'gamma'; 'highGamma'};
%freq_band_list = {'gamma'};

%% Define the channel groupings that you care about
% these channel groupings are constructed based on
% romanski's belt-PFC anatomy and more the classic cortico-cortico
% connections in V1 

     PFClowmid_ACupmid   =  {'*AC_ch08' '*PFC_ch13'
                             '*AC_ch08' '*PFC_ch14'
                             '*AC_ch08' '*PFC_ch15'
                             '*AC_ch08' '*PFC_ch16'
                             '*AC_ch08' '*PFC_ch17'
                             '*AC_ch09' '*PFC_ch13'
                             '*AC_ch09' '*PFC_ch14'
                             '*AC_ch09' '*PFC_ch15'
                             '*AC_ch09' '*PFC_ch16'
                             '*AC_ch09' '*PFC_ch17'
                             '*AC_ch10' '*PFC_ch13'
                             '*AC_ch10' '*PFC_ch14'
                             '*AC_ch10' '*PFC_ch15'
                             '*AC_ch10' '*PFC_ch16'
                             '*AC_ch10' '*PFC_ch17'
                             '*AC_ch11' '*PFC_ch13'
                             '*AC_ch11' '*PFC_ch14'
                             '*AC_ch11' '*PFC_ch15'
                             '*AC_ch11' '*PFC_ch16'
                             '*AC_ch11' '*PFC_ch17'
                             '*AC_ch12' '*PFC_ch13'
                             '*AC_ch12' '*PFC_ch14'
                             '*AC_ch12' '*PFC_ch15'
                             '*AC_ch12' '*PFC_ch16'
                             '*AC_ch12' '*PFC_ch17'};   

       PFCupmid_ACupmid  =  {'*AC_ch08' '*PFC_ch08'
                             '*AC_ch08' '*PFC_ch09'
                             '*AC_ch08' '*PFC_ch10'
                             '*AC_ch08' '*PFC_ch11'
                             '*AC_ch08' '*PFC_ch12'
                             '*AC_ch09' '*PFC_ch08'
                             '*AC_ch09' '*PFC_ch09'
                             '*AC_ch09' '*PFC_ch10'
                             '*AC_ch09' '*PFC_ch11'
                             '*AC_ch09' '*PFC_ch12'
                             '*AC_ch10' '*PFC_ch08'
                             '*AC_ch10' '*PFC_ch09'
                             '*AC_ch10' '*PFC_ch10'
                             '*AC_ch10' '*PFC_ch11'
                             '*AC_ch10' '*PFC_ch12'
                             '*AC_ch11' '*PFC_ch08'
                             '*AC_ch11' '*PFC_ch09'
                             '*AC_ch11' '*PFC_ch10'
                             '*AC_ch11' '*PFC_ch11'
                             '*AC_ch11' '*PFC_ch12'
                             '*AC_ch12' '*PFC_ch08'
                             '*AC_ch12' '*PFC_ch09'
                             '*AC_ch12' '*PFC_ch10'
                             '*AC_ch12' '*PFC_ch11'
                             '*AC_ch12' '*PFC_ch12'}; 
    
        PFCdeep_ACupmid =  {'*AC_ch08' '*PFC_ch18'
                            '*AC_ch08' '*PFC_ch19'
                            '*AC_ch08' '*PFC_ch20'
                            '*AC_ch08' '*PFC_ch21'
                            '*AC_ch08' '*PFC_ch22'
                            '*AC_ch09' '*PFC_ch18'
                            '*AC_ch09' '*PFC_ch19'
                            '*AC_ch09' '*PFC_ch20'
                            '*AC_ch09' '*PFC_ch21'
                            '*AC_ch09' '*PFC_ch22'
                            '*AC_ch10' '*PFC_ch18'
                            '*AC_ch10' '*PFC_ch19'
                            '*AC_ch10' '*PFC_ch20'
                            '*AC_ch10' '*PFC_ch21'
                            '*AC_ch10' '*PFC_ch22'
                            '*AC_ch11' '*PFC_ch18'
                            '*AC_ch11' '*PFC_ch19'
                            '*AC_ch11' '*PFC_ch20'
                            '*AC_ch11' '*PFC_ch21'
                            '*AC_ch11' '*PFC_ch22'
                            '*AC_ch12' '*PFC_ch18'
                            '*AC_ch12' '*PFC_ch19'
                            '*AC_ch12' '*PFC_ch20'
                            '*AC_ch12' '*PFC_ch21'
                            '*AC_ch12' '*PFC_ch22'};   
        
        PFCdeep_ACdeep =   {'*AC_ch18' '*PFC_ch18'
                            '*AC_ch18' '*PFC_ch19'
                            '*AC_ch18' '*PFC_ch20'
                            '*AC_ch18' '*PFC_ch21'
                            '*AC_ch18' '*PFC_ch22'
                            '*AC_ch19' '*PFC_ch18'
                            '*AC_ch19' '*PFC_ch19'
                            '*AC_ch19' '*PFC_ch20'
                            '*AC_ch19' '*PFC_ch21'
                            '*AC_ch19' '*PFC_ch22'
                            '*AC_ch20' '*PFC_ch18'
                            '*AC_ch20' '*PFC_ch19'
                            '*AC_ch20' '*PFC_ch20'
                            '*AC_ch20' '*PFC_ch21'
                            '*AC_ch20' '*PFC_ch22'
                            '*AC_ch21' '*PFC_ch18'
                            '*AC_ch21' '*PFC_ch19'
                            '*AC_ch21' '*PFC_ch20'
                            '*AC_ch21' '*PFC_ch21'
                            '*AC_ch21' '*PFC_ch22'
                            '*AC_ch22' '*PFC_ch18'
                            '*AC_ch22' '*PFC_ch19'
                            '*AC_ch22' '*PFC_ch20'
                            '*AC_ch22' '*PFC_ch21'
                            '*AC_ch22' '*PFC_ch22'};   
       
        PFCupper_ACupper =  {'*AC_ch03' '*PFC_ch03'
                             '*AC_ch03' '*PFC_ch04'
                             '*AC_ch03' '*PFC_ch05'
                             '*AC_ch03' '*PFC_ch06'
                             '*AC_ch03' '*PFC_ch07'
                             '*AC_ch04' '*PFC_ch03'
                             '*AC_ch04' '*PFC_ch04'
                             '*AC_ch04' '*PFC_ch05'
                             '*AC_ch04' '*PFC_ch06'
                             '*AC_ch04' '*PFC_ch07'
                             '*AC_ch05' '*PFC_ch03'
                             '*AC_ch05' '*PFC_ch04'
                             '*AC_ch05' '*PFC_ch05'
                             '*AC_ch05' '*PFC_ch06'
                             '*AC_ch05' '*PFC_ch07'
                             '*AC_ch06' '*PFC_ch03'
                             '*AC_ch06' '*PFC_ch04'
                             '*AC_ch06' '*PFC_ch05'
                             '*AC_ch06' '*PFC_ch06'
                             '*AC_ch06' '*PFC_ch07'
                             '*AC_ch07' '*PFC_ch03'
                             '*AC_ch07' '*PFC_ch04'
                             '*AC_ch07' '*PFC_ch05'
                             '*AC_ch07' '*PFC_ch06'
                             '*AC_ch07' '*PFC_ch07'};   


%% for each freqency do what is below

for fb = 1:length(freq_band_list)

            Frequency_Band = freq_band_list{fb};

            % Check the ChanWise_SpecEval OnlyPrior output. If the No_SigChans were found generate a null file, else modify the labels.  
            
            Condition = 'OnlyPrior';
            
            NO_PriorOnly_SigChans_fn = sprintf('NO_SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
            nullcheck_1 = dir(fullfile(DATA_DIR, NO_PriorOnly_SigChans_fn));

            if ~isempty(nullcheck_1)
            NULL_fn_1 = fullfile(SAVE_DIR, sprintf('PriorOnly_NO_SigGroups_%s_%s_%s_%s.txt', RecDate, Epoch, Behavior, Frequency_Band));  
            Null_1 = fopen(NULL_fn_1, 'w');
            
            else
    
            PriorOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
            files = dir(fullfile(DATA_DIR, PriorOnly_SigChans_fn));
                
            if ~isempty(files)
            PriorOnly_SigChans = load(fullfile(DATA_DIR, PriorOnly_SigChans_fn), 'significant_channels');
            PriorOnly_modifiedCellArray   = regexprep(PriorOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');
            end
    
            end
            
            % Check the ChanWise_SpecEval OnlyPrior output. If the No_SigChans were found generate a null file, else modify the labels.  

            Condition = 'OnlyPretone';

            NO_PretoneOnly_SigChans_fn = sprintf('NO_SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
            nullcheck_2 = dir(fullfile(DATA_DIR, NO_PretoneOnly_SigChans_fn));

            if ~isempty(nullcheck_2)
            NULL_fn_2 = fullfile(SAVE_DIR, sprintf('PretoneOnly_NO_SigGroups_%s_%s_%s_%s.txt', RecDate, Epoch, Behavior, Frequency_Band));  
            Null_2 = fopen(NULL_fn_2, 'w');
            else

            PretoneOnly_SigChans_fn = sprintf('SigChannels_%s_%s_%s_%s_%s.mat', RecDate, Epoch, Condition, Behavior, Frequency_Band);
            files = dir(fullfile(DATA_DIR, PretoneOnly_SigChans_fn));
            
            if ~isempty(files)
                PretoneOnly_SigChans = load(fullfile(DATA_DIR, PretoneOnly_SigChans_fn), 'significant_channels');
                PretoneOnly_modifiedCellArray = regexprep(PretoneOnly_SigChans.significant_channels , '^(D1_|D2_|D3_|D4_)', '*');
            end

            end

            % Filters out and aves the significant channels out of our preestablished channel groups in the fieldtrip format
            
            % OnlyPrior 
            if exist('PriorOnly_modifiedCellArray', 'var')
                PriorOnly_Chans_Filtered.PFClowmid_ACupmid = filterChannelGrouping(PFClowmid_ACupmid, PriorOnly_modifiedCellArray);
                PriorOnly_Chans_Filtered.PFCupmid_ACupmid  = filterChannelGrouping(PFCupmid_ACupmid, PriorOnly_modifiedCellArray);
                PriorOnly_Chans_Filtered.PFCdeep_ACupmid   = filterChannelGrouping(PFCdeep_ACupmid, PriorOnly_modifiedCellArray);
                PriorOnly_Chans_Filtered.PFCdeep_ACdeep    = filterChannelGrouping(PFCdeep_ACdeep, PriorOnly_modifiedCellArray);
                PriorOnly_Chans_Filtered.PFCupper_ACupper  = filterChannelGrouping(PFCupper_ACupper, PriorOnly_modifiedCellArray);
                Condition = 'OnlyPrior';
                file_name_1 = fullfile(SAVE_DIR, sprintf('PriorOnly_GroupedSigChans_%s_%s_%s_%s.mat', RecDate, Epoch, Behavior, Frequency_Band));     
                save(file_name_1, 'PriorOnly_Chans_Filtered');
            end
            
            % PretoneOnly
            if exist('PretoneOnly_modifiedCellArray', 'var')
                PretoneOnly_Chans_Filtered.PFClowmid_ACupmid = filterChannelGrouping(PFClowmid_ACupmid, PretoneOnly_modifiedCellArray);
                PretoneOnly_Chans_Filtered.PFCupmid_ACupmid  = filterChannelGrouping(PFCupmid_ACupmid, PretoneOnly_modifiedCellArray);
                PretoneOnly_Chans_Filtered.PFCdeep_ACupmid   = filterChannelGrouping(PFCdeep_ACupmid, PretoneOnly_modifiedCellArray);
                PretoneOnly_Chans_Filtered.PFCdeep_ACdeep    = filterChannelGrouping(PFCdeep_ACdeep, PretoneOnly_modifiedCellArray);
                PretoneOnly_Chans_Filtered.PFCupper_ACupper  = filterChannelGrouping(PFCupper_ACupper, PretoneOnly_modifiedCellArray);
                Condition = 'OnlyPretone';
                file_name_2 = fullfile(SAVE_DIR, sprintf('PretoneOnly_GroupedSigChans_%s_%s_%s_%s.mat', RecDate, Epoch, Behavior, Frequency_Band));           
                save(file_name_2,'PretoneOnly_Chans_Filtered');
            end

            clear PretoneOnly_modifiedCellArray
            clear PriorOnly_modifiedCellArray

end

clear all

%% Function to filter the channel groupings
function filteredArray = filterChannelGrouping(channelGrouping, modifiedCellArray)
    filteredArray = cell(0, 2);  % Initialize an empty cell array with two columns
    for i = 1:size(channelGrouping, 1)
        if ismember(channelGrouping{i, 1}, modifiedCellArray) && ismember(channelGrouping{i, 2}, modifiedCellArray)
            filteredArray = [filteredArray; channelGrouping(i, :)];  % Append the row if both elements are in modifiedCellArray
        end
    end
end
