function [null_threshold] = getGrangerNullThresh(Send_Chan_lb, Rec_Chan_lb, Bstrp_Iteration, Condition, Frequency_Band, input_data, varargin)

% Computes null threshold for a given session, condition, and frequency band

% Syntax:   null_thresh = getGrangerNullThresh(Send_Chan_lb, Rec_Chan_lb, Bstrp_Iteration, Condition, Frequency_Band, data); 
%           null_thresh = getGrangerNullThresh(Send_Chan_lb, Rec_Chan_lb, Bstrp_Iteration, Condition, Frequency_Band, [], Animal, RecDate, Epoch, datadir);
%           null_thresh = getGrangerNullThresh(Send_Chan_lb, Rec_Chan_lb, Bstrp_Iteration, Condition, Frequency_Band, Animal, RecDate, Epoch, datadir);
% 
% Inputs:
%    data    - (optional) Direct data input. Can be empty if using parameters.
%    Animal  - (optional) Animal identifier
%    RecDate - (optional) Session date
%    Epoch   - (optional) e.g. 'testToneOnset'


    if nargin == 6 && ~isempty(input_data)
        % If data is provided directly, use it.
        V = input_data;
        fprintf('Using provided data.\n');
    elseif nargin == 10 || nargin == 9
        if nargin == 10
            % If parameters are provided, use them to load data.
            Animal  = varargin{1};
            RecDate = varargin{2};
            Epoch   = varargin{3};
            datadir = varargin{4};
        elseif nargin == 9
            % If parameters are provided without data argument.
            Animal  = input_data;
            RecDate = varargin{1};
            Epoch   = varargin{2};
            datadir = varargin{3};
        end
            
            % Establish file and save directories
            addpath(genpath(datadir));
            fName = sprintf('%s-%s_bdLFP_%s_ft.mat', Animal, RecDate, Epoch);
            
            V = load(fullfile(datadir,fName));
    else
        error('Invalid input. Provide either data or parameters (Animal, RecDate, Epoch).');
    end

%% pull frequency band indices
if strcmp(Frequency_Band,'theta')
   Freq_Band_Indices = 5:8;
elseif strcmp(Frequency_Band,'alpha')
   Freq_Band_Indices = 8:14;
elseif strcmp(Frequency_Band,'beta')
   Freq_Band_Indices = 15:30;
elseif strcmp(Frequency_Band,'gamma')
   Freq_Band_Indices = 31:54;
elseif strcmp(Frequency_Band,'highGamma')
   Freq_Band_Indices = 66:100;
end

%% filter out correct prior/pretone only
if strcmp(Condition,'PriorOnly')                                                    % another strategy to filter prior/pretone correct/wrong data
    condition_indices = find((strcmp(V.prior,'H') | strcmp(V.prior,'L')) & ...      % functionally equivalent to Corey/Taku's iSelect chunk
                              V.pretone == 'N' & V.choice ~= 'n' & V.err == 'c');
elseif strcmp(Condition, 'PretoneOnly')
    condition_indices = find((V.pretone == 'H' | V.pretone =='L') & ...
                            strcmp(V.prior,'N') & V.choice ~= 'n' & V.err == 'c');
end
selected_trials = V.trial_id(condition_indices);

idx = find(ismember(V.trial_id,selected_trials)==1);

data.time = V.data.time(idx);
data.trial = V.data.trial(idx);
data.sampleinfo = V.data.sampleinfo(idx,:);
data.label = V.data.label;
data.fsample = V.data.fsample;

%% shuffle each channel and save each shuffled granger
num_trials      = length(data.trial);
length_trial    = length(data.trial{1});

granger_shuff   = zeros(2,2,length(Freq_Band_Indices),Bstrp_Iteration);

for i = 1:Bstrp_Iteration
    
    shuffled_trials = cell(1,num_trials);

    for z = 1:num_trials
        randomcolidx       = randperm(length_trial);
        trial_shuff        = data.trial{z}(:,randomcolidx);
        shuffled_trials{z} = trial_shuff;
    end

    data_shuff = data;
    data_shuff.trial = shuffled_trials;
  
    % cross spectral density matrix computation
    cfg = [];
    cfg.method    = 'mtmfft';
    cfg.taper     = 'dpss';
    cfg.output    = 'fourier';
    cfg.tapsmofrq = 4;
    cfg.pad       = 1;
    cfg.foilim    = [0 100];
    freq_shuff    = ft_freqanalysis(cfg,data_shuff);
    
    % Granger computation
    cfg             = [];
    cfg.method      = 'granger';
    cfg.channel     = {Send_Chan_lb Rec_Chan_lb};
    cfg.channelcmb  = {Send_Chan_lb Rec_Chan_lb};
    
    g                      = ft_connectivityanalysis(cfg,freq_shuff);
    granger_shuff(:,:,:,i) = g.grangerspctrm(:,:,Freq_Band_Indices);
end

% calculate null threshold
null_threshold = mean(granger_shuff, 4);


