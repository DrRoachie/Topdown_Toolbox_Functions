function [preCueOnset_PSI_grid, testToneOnset_PSI_grid, ...
          preCue_averagedGrid, testTone_averagedGrid, ...
          PFC_2_AC_PreCue, AC_2_PFC_PreCue, ...
          PFC_2_AC_TestTone, AC_2_PFC_TestTone, ...
          group_data] = build_epoch_PSI(Condition, animals, rootdir)
% build_epoch_PSI
% Assemble and summarize PSI values for the Epoch analysis.
%
% INPUTS
%   Condition : char/string (e.g., 'OnlyPrior', 'OnlyPretone', or 'Both')
%               (currently used for bookkeeping; filenames are found by epoch/band)
%   animals   : cellstr of animal names (e.g., {'MrCassius','MrM'})
%   rootdir   : char/string; path to the analysis root directory containing session subfolders
%
% OUTPUTS
%   preCueOnset_PSI_grid   : 20x20 cell, each cell = [n_rows x n_freqs] PSI rows for PreCue
%   testToneOnset_PSI_grid : 20x20 cell, each cell = [n_rows x n_freqs] PSI rows for TestTone
%   preCue_averagedGrid    : 20x20 numeric, per-pair PSI averaged over theta columns (5:8)
%   testTone_averagedGrid  : 20x20 numeric, per-pair PSI averaged over theta columns (5:8)
%   PFC_2_AC_PreCue        : vector of |theta-avg PSI| where PreCue avg > 0  (putative AC←PFC)
%   AC_2_PFC_PreCue        : vector of |theta-avg PSI| where PreCue avg < 0  (putative PFC←AC)
%   PFC_2_AC_TestTone      : vector of |theta-avg PSI| where TestTone avg > 0
%   AC_2_PFC_TestTone      : vector of |theta-avg PSI| where TestTone avg < 0
%   group_data             : 1x4 cell = {PFC_2_AC_PreCue, AC_2_PFC_PreCue, PFC_2_AC_TestTone, AC_2_PFC_TestTone}
%
% NOTES
% - Frequency band fixed to theta (columns 5:8), matching your current pipeline.
% - Empty channel-pair slots remain NaN in the averaged grids.
% - This function doesn’t assume a specific session-folder naming scheme beyond being directories.

% ----------------------------- Book-keeping ------------------------------
% (Condition kept for compatibility/metadata; not used to filter files here)
if nargin < 3
    error('Usage: build_epoch_PSI(Condition, animals, rootdir)');
end
if ischar(animals) || isstring(animals)
    animals = cellstr(animals);
end
Frequency_Band = 'theta';     % fixed per your current analysis
columnRange    = 5:8;         % theta indices

% Channels ch03 .. ch22 (20 each)
channels = arrayfun(@(x) sprintf('ch%02d', x), 3:22, 'UniformOutput', false);

% Initialize grids
preCueOnset_PSI_grid   = cell(20, 20);
testToneOnset_PSI_grid = cell(20, 20);

% Gather session folders (all subdirectories excluding . and ..)
allDir = dir(rootdir);
isValid = [allDir.isdir] & ~ismember({allDir.name},{'.','..'});
sessions = allDir(isValid);

% ---------------- Loop sessions → animals → files → populate grids --------
for i = 1:numel(sessions)
    RecDate = sessions(i).name;

    for j = 1:numel(animals)
        Animal = animals{j};
        session_path = fullfile(rootdir, RecDate, Animal);

        if ~exist(session_path, 'dir')
            continue; % skip missing animal folder for this session
        end

        % find epoch files for theta band
        preCue_file   = dir(fullfile(session_path, ['*', 'preCueOnset_'  Frequency_Band, '*PSI_data.mat']));
        testTone_file = dir(fullfile(session_path, ['*', 'testToneOnset_', Frequency_Band, '*PSI_data.mat']));

        if isempty(preCue_file) || isempty(testTone_file)
            % If either epoch file set is missing, skip this (session,animal)
            continue;
        end

        % We process all found files; map by index without assuming names match
        % If you need strict pairing, you can add matching by shared tokens.
        K = min(numel(preCue_file), numel(testTone_file));
        for k = 1:K
            preCue_filename   = fullfile(session_path, preCue_file(k).name);
            testTone_filename = fullfile(session_path, testTone_file(k).name);

            try
                S1 = load(preCue_filename,   'PSI_preCueOnset');
                S2 = load(testTone_filename, 'PSI_testToneOnset');
            catch
                % Skip if load fails
                continue;
            end
            if ~isfield(S1,'PSI_preCueOnset') || ~isfield(S2,'PSI_testToneOnset')
                continue;
            end

            preCue_Matrix   = S1.PSI_preCueOnset.psispctrm;
            preCue_labels   = S1.PSI_preCueOnset.labelcmb;
            testTone_Matrix = S2.PSI_testToneOnset.psispctrm;
            testTone_labels = S2.PSI_testToneOnset.labelcmb;

            % ----- Fill PreCue grid (append rows per channel-pair) -----
            for rowIdx = 1:size(preCue_Matrix, 1)
                pfcLabel = preCue_labels{rowIdx, 1}; % e.g., 'D1_PFC_ch05'
                acLabel  = preCue_labels{rowIdx, 2}; % e.g., 'D3_AC_ch04'

                pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');
                acChannel  = regexp(acLabel,  'AC_ch\d+',  'match', 'once');

                if isempty(pfcChannel) || isempty(acChannel)
                    continue;
                end

                pfcKey = strrep(pfcChannel, 'PFC_', '');
                acKey  = strrep(acChannel,  'AC_',  '');

                pfcIdx = find(strcmp(channels, pfcKey), 1);
                acIdx  = find(strcmp(channels, acKey),  1);

                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(preCueOnset_PSI_grid{pfcIdx, acIdx})
                        preCueOnset_PSI_grid{pfcIdx, acIdx} = preCue_Matrix(rowIdx, :);
                    else
                        preCueOnset_PSI_grid{pfcIdx, acIdx} = [preCueOnset_PSI_grid{pfcIdx, acIdx}; preCue_Matrix(rowIdx, :)]; %#ok<AGROW>
                    end
                end
            end

            % ----- Fill TestTone grid (append rows per channel-pair) -----
            for rowIdx = 1:size(testTone_Matrix, 1)
                pfcLabel = testTone_labels{rowIdx, 1};
                acLabel  = testTone_labels{rowIdx, 2};

                pfcChannel = regexp(pfcLabel, 'PFC_ch\d+', 'match', 'once');
                acChannel  = regexp(acLabel,  'AC_ch\d+',  'match', 'once');

                if isempty(pfcChannel) || isempty(acChannel)
                    continue;
                end

                pfcKey = strrep(pfcChannel, 'PFC_', '');
                acKey  = strrep(acChannel,  'AC_',  '');

                pfcIdx = find(strcmp(channels, pfcKey), 1);
                acIdx  = find(strcmp(channels, acKey),  1);

                if ~isempty(pfcIdx) && ~isempty(acIdx)
                    if isempty(testToneOnset_PSI_grid{pfcIdx, acIdx})
                        testToneOnset_PSI_grid{pfcIdx, acIdx} = testTone_Matrix(rowIdx, :);
                    else
                        testToneOnset_PSI_grid{pfcIdx, acIdx} = [testToneOnset_PSI_grid{pfcIdx, acIdx}; testTone_Matrix(rowIdx, :)]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

% ---------------------- Average theta columns per pair --------------------
preCue_averagedGrid   = nan(20, 20);
testTone_averagedGrid = nan(20, 20);

for p = 1:20
    for a = 1:20
        if ~isempty(preCueOnset_PSI_grid{p,a})
            preCue_averagedGrid(p,a) = mean(preCueOnset_PSI_grid{p,a}(:, columnRange), 'all', 'omitnan');
        end
        if ~isempty(testToneOnset_PSI_grid{p,a})
            testTone_averagedGrid(p,a) = mean(testToneOnset_PSI_grid{p,a}(:, columnRange), 'all', 'omitnan');
        end
    end
end

% -------------------- Directional splits & grouped data -------------------
abs_preCue   = abs(preCue_averagedGrid);
abs_testTone = abs(testTone_averagedGrid);

mask_nonzero_pre  = preCue_averagedGrid   ~= 0;
mask_nonzero_test = testTone_averagedGrid ~= 0;

PFC_2_AC_PreCue   = abs_preCue(  mask_nonzero_pre  & preCue_averagedGrid   > 0);
AC_2_PFC_PreCue   = abs_preCue(  mask_nonzero_pre  & preCue_averagedGrid   < 0);
PFC_2_AC_TestTone = abs_testTone(mask_nonzero_test & testTone_averagedGrid > 0);
AC_2_PFC_TestTone = abs_testTone(mask_nonzero_test & testTone_averagedGrid < 0);

group_data = {PFC_2_AC_PreCue, AC_2_PFC_PreCue, PFC_2_AC_TestTone, AC_2_PFC_TestTone};

% (Optional) quick summary stats, if needed later:
% means = cellfun(@mean, group_data);
% SEMs  = cellfun(@(x) std(x)/sqrt(max(1,numel(x))), group_data);
% CIs   = 1.96 * SEMs;

end
