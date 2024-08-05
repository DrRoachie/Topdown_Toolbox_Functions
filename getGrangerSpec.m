
function [granger, std_granger] = getGrangerSpec(varargin)

% Computes granger spectrum and confidence intervals for a given session, condition, and frequency band

% Syntax:   granger = getGrangerSpec(freq, data, Send_Chan_lb, Rec_Chan_lb);  *** used this way in Pairwise Connectivity Test Uber
%           granger = getGrangerSpec(Send_Chan_lb, Rec_Chan_lb, Condition, Frequency_Band, [], Animal, RecDate, Epoch, datadir);
%           granger = getGrangerSpec(Send_Chan_lb, Rec_Chan_lb, Condition, Frequency_Band, Animal, RecDate, Epoch, datadir);*** used this way generally
% 
% Inputs:
%    freq           - Cross spectral density matrix used for granger calculation
%    Send_Chan_lb   - Channel label for sender. Used in channel (pair) specification in granger calculation
%    Rec_Chan_lb    - Channel label for receiver. Used in channel (pair) specification in granger calculation
%    Condition      - PriorOnly or PretoneOnly
%    Frequency_Band - 'theta', 'alpha', 'beta', 'gamma', or 'highGamma'
%    data           - Direct data input. Can be empty if using parameters.
%    Animal         - Animal identifier
%    RecDate        - Session date
%    Epoch          - e.g. 'testToneOnset'

if nargin == 5
    % If cross spectral density matrix is provided, use it directly in calculation of granger spectrum
    freq = varargin{1};
    data = varargin{2};
    Send_Chan_lb = varargin{3};
    Rec_Chan_lb  = varargin{4};
    Frequency_Band = varargin{5};

else % Information provided for to calculate cross spectral density

    if nargin == 9
        % If parameters are provided, use them to load data.
        Send_Chan_lb   = varargin{1};
        Rec_Chan_lb    = varargin{2};
        Condition      = varargin{3};
        Frequency_Band = varargin{4};
        Animal  = varargin{6};
        RecDate = varargin{7};
        Epoch   = varargin{8};
        datadir = varargin{9};

    elseif nargin == 8
        % If parameters are provided without data argument.
        Send_Chan_lb   = varargin{1};
        Rec_Chan_lb    = varargin{2};
        Condition      = varargin{3};
        Frequency_Band = varargin{4};
        Animal  = varargin{5};
        RecDate = varargin{6};
        Epoch   = varargin{7};
        datadir = varargin{8};

    end
    % Establish file and save directories
    addpath(genpath(datadir));
    fName = sprintf('%s-%s_bdLFP_%s_ft.mat', Animal, RecDate, Epoch);
    
    V = load(fName);

    
    % filter out correct prior/pretone only
    
    if strcmp(Condition,'PriorOnly')                                                    % another strategy to filter prior/pretone correct/wrong data
        condition_indices = find((strcmp(V.prior,'H') | strcmp(V.prior,'L')) & ...      % functionally equivalent to Corey/Taku's iSelect function
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
    
    % non-parametric computation of cross-spectral density matrix (avged across trials)
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'fourier';
    cfg.tapsmofrq = 4;
    cfg.pad = 1;
    cfg.foilim = [0 100];
    freq = ft_freqanalysis(cfg,data);

end


% Pull frequency bands
    
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

% compute granger

cfg             = [];
cfg.method      = 'granger';
cfg.channel     = {Send_Chan_lb Rec_Chan_lb}; % you need both cfg.channel and cfg.channelcmb, ft wants you to read in all of them and then subselect. We just subselect twice.
cfg.channelcmb  = {Send_Chan_lb Rec_Chan_lb};

granger                 = ft_connectivityanalysis(cfg,freq);

granger.grangerspctrm   = granger.grangerspctrm(:,:,Freq_Band_Indices);     % pull desired frequency band
granger.freq            = granger.freq(Freq_Band_Indices);                  % adjust freq
granger.dof             = length(data.trial);                               % because ft_connectivityanalysis does not output trial number, we define it here. Needed for MonteCarlo Estimation


%% compute variance (5-fold k-1 cross validation method for subselection)

num_trials = length(data.trial);
num_folds    = 5;
kk           = ceil(num_trials/num_folds);
ci_granger   = [];          

for i = 1:kk-1
    rand_indices               = randperm(num_trials, kk);

    data_temp                  = data;
    data_temp.time             = data.time(rand_indices);
    data_temp.trial            = data.trial(rand_indices);
    data_temp.sampleinfo       = data.sampleinfo(rand_indices,:);
  
    cfg             = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'fourier';
    cfg.tapsmofrq   = 4;
    cfg.pad         = 1;
    cfg.foilim      = [0 100];

    freq_temp = ft_freqanalysis(cfg,data_temp);
    
    cfg             = [];
    cfg.method      = 'granger';
    cfg.channel     = {Send_Chan_lb Rec_Chan_lb};
    cfg.channelcmb  = {Send_Chan_lb Rec_Chan_lb};
    
    ci_granger(:,:,:,i) = ft_connectivityanalysis(cfg,freq_temp).grangerspctrm(:,:,Freq_Band_Indices);
end

z = 1.96;
std_granger = z * (std(ci_granger, 0, 4) / sqrt(size(ci_granger, 4)));      % calculate standard deviation for variance



