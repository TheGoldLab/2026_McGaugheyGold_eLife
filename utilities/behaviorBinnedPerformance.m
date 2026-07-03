function [LSF_binned_performance, HSF_binned_performance, binned_behavior_diff] = ...
    behaviorBinnedPerformance(dat, trialType)
% behaviorBinnedPerformance Calculate binned behavioral performance
%
% This function bins trials by viewing duration and calculates performance
% (percent correct) for LSF and HSF trials.
%
% Inputs:
%   dat       - Pre-processed data table (should be subset and filtered as needed)
%   trialType - 'switch' (default) or 'nonswitch'
%
% Outputs:
%   LSF_binned_performance - [numSessions x numBins] Performance for LSF trials
%   HSF_binned_performance - [numSessions x numBins] Performance for HSF trials
%   binned_behavior_diff   - [numSessions x numBins] Difference (LSF - HSF)
%
% Binning scheme:
%   Bin 1: 100-225 ms
%   Bin 2: 225-375 ms
%   Bin 3: 375-600 ms
%   Bin 4: 600-1200 ms

if nargin < 2 || isempty(trialType)
    trialType = 'switch';
end

switch trialType
    case 'switch'
        wantSwitch = 1;
    case 'nonswitch'
        wantSwitch = 0;
    otherwise
        error('behaviorBinnedPerformance:badTrialType', ...
            "trialType must be 'switch' or 'nonswitch'");
end

%% Get unique sessions
uniqueSessionNames = unique(dat.ses_ID);
numSessions = length(uniqueSessionNames);
numBins = 4; % Fixed number of bins

%% Initialize output matrices
LSF_binned_performance = nan(numSessions, numBins);
HSF_binned_performance = nan(numSessions, numBins);

%% Loop through sessions
for ss = 1:numSessions
    
    % Get data for this session
    Lses = strcmp(dat.ses_ID, uniqueSessionNames{ss}) & dat.coh_final >= 50;
    
    % Collect data table of binned averages
    % Columns: Switch frequency (HR), Switch(1)/non-switch(0), Viewing duration, RT, Correct(1)/incorrect(0)
    select_dat = table2array(dat(Lses, [2 5 9 11 12]));
    
    %% Assign trials to bins based on viewing duration
    num_trials = size(select_dat, 1);
    bin_idx = nan(num_trials, 1);
    
    for tt = 1:num_trials
        if select_dat(tt, 3) <= 225
            bin_idx(tt, 1) = 1;
        elseif select_dat(tt, 3) > 225 && select_dat(tt, 3) <= 375
            bin_idx(tt, 1) = 2;
        elseif select_dat(tt, 3) > 375 && select_dat(tt, 3) <= 600
            bin_idx(tt, 1) = 3;
        elseif select_dat(tt, 3) > 600
            bin_idx(tt, 1) = 4;
        end
    end
    
    % Add bin index to data
    select_dat_bins = [select_dat, bin_idx];

    %% Create logicals for conditions
    Lsf = select_dat_bins(:, 1) == 2;                % LSF trials
    LtrialType = select_dat_bins(:, 2) == wantSwitch; % Switch or non-switch trials

    % Extract relevant data (viewing duration, RT, correct, bin_idx)
    LSF_trialType = select_dat_bins(Lsf & LtrialType, 3:6);
    HSF_trialType = select_dat_bins(~Lsf & LtrialType, 3:6);

    %% Calculate performance (mean correct) for each bin
    for bb = 1:numBins
        % Column 3 is correct/incorrect (1/0)
        % Column 4 is bin_idx
        LSF_trials_in_bin = LSF_trialType(LSF_trialType(:, 4) == bb, 3);
        HSF_trials_in_bin = HSF_trialType(HSF_trialType(:, 4) == bb, 3);

        LSF_binned_performance(ss, bb) = mean(LSF_trials_in_bin);
        HSF_binned_performance(ss, bb) = mean(HSF_trials_in_bin);
    end
end

%% Calculate difference
binned_behavior_diff = LSF_binned_performance - HSF_binned_performance;

end
