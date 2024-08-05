function [ coh_new ] = reshape_coherence( coh, iElectrode_AC, iElectrode_PFC )
%reshape_coherence Summary of this function goes here
%   select and reshape coherence calculated by fieldtrip toolbox
%   coh: struct of coherence from fieldtrip
%   iElectrode_AC: AC electrode ID
%   iElectrode_PFC: PFC electrode ID

lb_cmb = coh.labelcmb; % channel combination label
C = coh.cohspctrm; % original matrix of coherence (cannel_cmb x freq)

temp_a = lb_cmb(:,1); % 1st column of lb_cmb (should be AC channels)
temp_b = lb_cmb(:,2); % 2nd column of lb_cmb (should be PFC channels)

% select electrode
i_a = contains(temp_a,iElectrode_AC);
i_b = contains(temp_b,iElectrode_PFC);
C_select = C(and(i_a,i_b),:);
temp_a = temp_a(and(i_a,i_b));
temp_b = temp_b(and(i_a,i_b));
lb_cmb_new = [temp_a temp_b];

lb_a = unique(temp_a); % list of channel (1st column) 
lb_b = unique(temp_b); % list of channel (2nd column)
n_a = numel(lb_a);
n_b = numel(lb_b);
n_fp = size(C_select,2); % number of frequency point 

% reshape matrix
C_new = reshape(C_select,n_a,n_b,n_fp);

% update coh by C_new and lb_cmb_new
coh_new = coh;
coh_new.labelcmb  = lb_cmb_new;
coh_new.label_ac  = lb_a;
coh_new.label_pfc = lb_b; 
coh_new.cohspctrm = C_select; % update
coh_new.cohspctrm_mat = C_new; % in matrix form (AC ch x PFC ch x freq)
% remove updated fields (to avoid confusion)
cfg = coh_new.cfg;
cfg = rmfield(cfg,'channelcmb');
cfg = rmfield(cfg,'channel');
coh_new.cfg = cfg;
end

