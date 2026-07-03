function [LSF_PREF_s, HSF_PREF_s] = ...
    processNeuralDataErrorTrials(dat, slide, bin_size, start_time, end_time, ...
    condition_selector, normalizationTerm_input)
% processNeuralDataErrorTrials Process incorrect trial neural data (simplified)
%
% This function bins spike counts, computes normalized firing rates, and pools
% trials across sessions for LSF and HSF switch trials (PREF direction only).
% Simplified version that only returns what's needed for plotting.
%
% Inputs:
%   dat                      - Pre-processed data table (should be INCORRECT trials)
%   slide                    - Bin slide/step size in msec (e.g., 10)
%   bin_size                 - Bin width in msec (e.g., 100)
%   start_time               - Analysis start time in msec (e.g., -2600)
%   end_time                 - Analysis end time in msec (e.g., 1200)
%   condition_selector       - [numUnits x 1] logical array or cell array from sensitivity_diff_labeled_N
%   normalizationTerm_input  - [numUnits x 1] normalization terms (from correct trials), or [] to calculate
%
% Outputs:
%   LSF_PREF_s  - [num_trials x numBins] LSF switch PREF trials (pooled across sessions)
%   HSF_PREF_s  - [num_trials x numBins] HSF switch PREF trials (pooled across sessions)

%% Set up timing/binning
bins = cat(2, ...
    (start_time:slide:end_time - bin_size)', ...
    (start_time + bin_size:slide:end_time)');
xax = mean(bins,2);
numBins = size(bins, 1);

%% Set up direction selection arrays
sfs = [2 6];
Ladapt = (xax >= -2300 & xax < 0)';
LdotsOn = false(numBins, 2);
ons = {-2300 + 1200 .* (0:sfs(1)), -2300 + 400 .* (0:sfs(2))};

for ff = 1:2
    for ii = 1:2:length(ons{ff})-1
        LdotsOn(xax >= ons{ff}(ii) & xax < ons{ff}(ii+1),ff) = true;
    end
end

Ldots = cat(2, repmat(LdotsOn,1,2), repmat(~LdotsOn,1,2))';

%% Get unique units
uniqueUnitNames = unique(dat.ses_ID);
numUnits = length(uniqueUnitNames);

%% Validate normalization term input
use_custom_norm = false;
if ~isempty(normalizationTerm_input)
    if length(normalizationTerm_input) ~= numUnits
        error('normalizationTerm_input must be empty or have length equal to numUnits (%d)', numUnits);
    end
    use_custom_norm = true;
    normalizationTerm = normalizationTerm_input(:);
else
    normalizationTerm = nan(numUnits, 1);
end

%% Validate condition selector
if ~isempty(condition_selector)
    % Handle if it's a cell array from sensitivity_diff_labeled_N
    if iscell(condition_selector)
        condition_selector = cell2mat(condition_selector(:,2)) >= 0;  % Convert to logical
    end
    
    if length(condition_selector) ~= numUnits
        error('condition_selector must be empty or have length equal to numUnits (%d)', numUnits);
    end
else
    condition_selector = true(numUnits, 1);  % Process all units
end

%% Initialize pooled trial arrays
LSF_PREF_s = [];
HSF_PREF_s = [];

%% Loop through units and pool trials
for uu = 1:numUnits
    
    disp(uniqueUnitNames{uu})
    
    % Filter by condition selector
    if ~condition_selector(uu)
        continue
    end
    
    % Get session data
    Lses = strcmp(dat.ses_ID, uniqueUnitNames{uu}) & dat.coh_final >= 50; 
    if sum(Lses) == 0
        continue
    end    
    
    session_data = dat(Lses, 1:end);
    dirs = [min(session_data.dir_final) max(session_data.dir_final)];
    
    % Set up condition structure
    [SF, ST, SD] = ndgrid(sfs, [0 1], dirs);
    unit = struct( ...
        'dirs',             dirs, ...
        'switch_frequency', num2cell(SF(:)), ...
        'switch_trial',     num2cell(ST(:)), ...
        'start_dir',        num2cell(SD(:)));
    
    dir_rates = [0 0];
    maxRate = nan(1, length(unit));
    
    % Loop through conditions and collect spike data
    for cc = 1:length(unit)
        
        trial_indices = find(session_data.HR == unit(cc).switch_frequency & ...
            session_data.dir_switch == unit(cc).switch_trial & ...
            session_data.dir_prefinal == unit(cc).start_dir);
        
        num_trials = length(trial_indices);
        if num_trials == 0
            continue
        end
        
        unit(cc).spike_rates_std = nan(num_trials, numBins);
        
        % Get binned spike counts
        for tt = 1:num_trials
            
            test_start = session_data.dots_on{trial_indices(tt)}(end);
            spikes = session_data.Unit_1{trial_indices(tt)} - test_start;
            
            trial_end_time = session_data.dots_off(trial_indices(tt)) - test_start;
            spike_rates_raw = nan(1, numBins);
            
            for bb = 1:find(bins(:,2) > trial_end_time,1)-1
                spike_count = sum(spikes >= bins(bb,1) & spikes < bins(bb,2));
                spike_rates_raw(bb) = spike_count * 1000/bin_size;
            end
            
            % Baseline subtract
            baseline_bin = find(bins(:,2) > -2500,1)-1;
            dots_on_bin = find(bins(:,2) > -2400,1)-1;
            baseline_activity = nanmean(spike_rates_raw(baseline_bin:dots_on_bin));
            
            unit(cc).spike_rates_std(tt,:) = spike_rates_raw - baseline_activity;
        end
        
        % Calculate normalization factor from adapting epoch (only if not using custom)
        if ~use_custom_norm
            num_cols = size(unit(cc).spike_rates_std, 2);
            Ladapt_valid = Ladapt(1:min(length(Ladapt), num_cols));
            if sum(Ladapt_valid) > 0
                maxRate(cc) = max(nanmean(unit(cc).spike_rates_std(:, Ladapt_valid)));
            end
        end
    end
    
    % Set normalization term
    if use_custom_norm
        normalizationTerm_use = normalizationTerm(uu);
    else
        normalizationTerm_use = max(maxRate);
    end
    
    % Normalize and determine PREF/NULL
    for cc = 1:length(unit)
        if ~isfield(unit(cc), 'spike_rates_std') || isempty(unit(cc).spike_rates_std)
            continue
        end
        
        unit(cc).spike_rates_std_norm = unit(cc).spike_rates_std / normalizationTerm_use;
        
        % Calculate dir_rates for PREF/NULL determination
        num_cols = size(unit(cc).spike_rates_std, 2);
        Ldots_valid = Ldots(cc, 1:min(size(Ldots, 2), num_cols));
        Ladapt_valid = Ladapt(1:min(length(Ladapt), num_cols));
        
        if sum(Ldots_valid & Ladapt_valid) > 0 && sum(~Ldots_valid & Ladapt_valid) > 0
            dir_rates = dir_rates + [ ...
                nanmean(reshape(unit(cc).spike_rates_std(:, Ldots_valid & Ladapt_valid),[],1)) .* size(unit(cc).spike_rates_std,1)./sum(Lses), ...
                nanmean(reshape(unit(cc).spike_rates_std(:, ~Ldots_valid & Ladapt_valid),[],1)) .* size(unit(cc).spike_rates_std,1)./sum(Lses)];
        end
    end
    
    % Determine PREF/NULL
    if dir_rates(2) > dir_rates(1)
        fprintf('%s = %d PREF\n', uniqueUnitNames{uu}, unit(1).dirs(2))
        [unit.dirs] = deal(unit(1).dirs([2 1]));
    elseif dir_rates(1) > dir_rates(2)
        fprintf('%s = %d PREF\n', uniqueUnitNames{uu}, unit(1).dirs(1))
    end
    
    % Pool switch trials only (PREF direction only)
    % Condition 1: LSF switch PREF
    Lcond_LSF = [unit.switch_trial] == 1 & ...
        [unit.switch_frequency] == sfs(1) & ...
        [unit.start_dir] == unit(1).dirs(1);
    
    if sum(Lcond_LSF) > 0 && isfield(unit(Lcond_LSF), 'spike_rates_std_norm')
        LSF_PREF_s = [LSF_PREF_s; unit(Lcond_LSF).spike_rates_std_norm];
    end
    
    % Condition 3: HSF switch PREF
    Lcond_HSF = [unit.switch_trial] == 1 & ...
        [unit.switch_frequency] == sfs(2) & ...
        [unit.start_dir] == unit(1).dirs(1);
    
    if sum(Lcond_HSF) > 0 && isfield(unit(Lcond_HSF), 'spike_rates_std_norm')
        HSF_PREF_s = [HSF_PREF_s; unit(Lcond_HSF).spike_rates_std_norm];
    end
end

end
