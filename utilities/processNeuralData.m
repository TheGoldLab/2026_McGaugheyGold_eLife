function [spike_rates_avg, spike_rates_SEM, spike_rates_avg_std_norm_avg, ...
    spike_rates_avg_std_norm_SEM, ROC_area, cell_selectivity, ...
    normalizationTerm, uniqueUnitNames, unit_example] = ...
    processNeuralData(dat, slide, bin_size, start_time, end_time, example_unit_idx, normalizationTerm_input)
% processNeuralData Process neural data from DotsRev experiment
%
% This function bins spike counts, computes firing rates (raw and normalized), 
% and calculates ROC areas for directional selectivity analysis.
%
% Inputs:
%   dat              - Pre-processed data table (should be subset and filtered as needed)
%   slide            - Bin slide/step size in msec (e.g., 10)
%   bin_size         - Bin width in msec (e.g., 100)
%   start_time       - Analysis start time in msec (e.g., -2600)
%   end_time         - Analysis end time in msec (e.g., 1200)
%   example_unit_idx - Unit index to save for rasterplot (e.g., 132), or [] for none
%
% Outputs:
%   spike_rates_avg               - [numUnits x numBins x 8] Average firing rates
%   spike_rates_SEM               - [numUnits x numBins x 8] SEMs for firing rates
%   spike_rates_avg_std_norm_avg  - [numUnits x numBins x 8] Normalized firing rates
%   spike_rates_avg_std_norm_SEM  - [numUnits x numBins x 8] SEMs for normalized rates
%   ROC_area                      - [numUnits x numBins x 4] ROC areas
%   cell_selectivity              - [numUnits x 1] Selectivity (PREF - NULL, raw)
%   normalizationTerm             - [numUnits x 1] Max firing rate for normalization
%   uniqueUnitNames             - Cell array of unit identifiers
%   unit_example                  - Unit struct for example neuron (empty if not saved)
%
% Matrix organization (dimension 3 for spike_rates outputs):
%   1. Low switch frequency switch (PREF)
%   2. Low switch frequency switch (NULL)
%   3. High switch frequency switch (PREF)
%   4. High switch frequency switch (NULL)
%   5. Low switch frequency non-switch (NULL)
%   6. Low switch frequency non-switch (PREF)
%   7. High switch frequency non-switch (NULL)
%   8. High switch frequency non-switch (PREF)
%
% PREF/NULL above describes the TEST-epoch direction, matching how every
% caller uses these columns (and matching how the ROC_area loop below
% assigns vals_pref/vals_null). It is NOT simply "trial started on
% unit(1).dirs(dd)": LSF has 1 and HSF has 5 mid-adapting-epoch reversals
% (both odd), so for switch trials (one more reversal at the adapting/test
% boundary) the starting direction and the test-epoch direction end up the
% same, but for non-switch trials (no reversal at that boundary) they end
% up opposite. That's why dd=1 ("starts on dirs(1)", the adapting-epoch
% preferred direction) maps to columns 1/3 (PREF) for switch trials but to
% columns 5/7 (NULL) for non-switch trials -- confirmed empirically against
% Figure 5-Figure Supplement 4 (see Fig5DE_FigSupp2_4.m), whose published,
% non-switch test-epoch PREF time course only matches columns 6/8.
%
% ROC_area matrix organization (dimension 3):
%   1. LSF non-switch trials
%   2. LSF switch trials
%   3. HSF switch trials
%   4. HSF non-switch trials

%% Set up timing/binning (from input parameters)
bins = cat(2, ...
    (start_time:slide:end_time - bin_size)', ...
    (start_time + bin_size:slide:end_time)');
xax = mean(bins,2);
numBins = size(bins, 1);

%% Make selection array to compute PREF/NULL direction
% LdotsOn(:,1) = LSF
% LdotsOn(:,2) = HSF
% Take the adapting epoch + 100 ms latency (-2300:100)
% Fill in dots_on bins for each column

sfs = [2 6];
Ladapt = (xax >= -2300 & xax < 0)';
LdotsOn = false(numBins, 2);
ons = {-2300 + 1200 .* (0:sfs(1)), -2300 + 400 .* (0:sfs(2))};

for ff = 1:2
    for ii = 1:2:length(ons{ff})-1
        LdotsOn(xax >= ons{ff}(ii) & xax < ons{ff}(ii+1),ff) = true;
    end
end

% (2 dirs x 2 switch x 2 hazard)
Ldots = cat(2, repmat(LdotsOn,1,2), repmat(~LdotsOn,1,2))';

%% Get unique session (unit) names
uniqueUnitNames = unique(dat.ses_ID);
numUnits = length(uniqueUnitNames);

%% Validate normalization term input
use_custom_norm = false;
if ~isempty(normalizationTerm_input)
    if length(normalizationTerm_input) ~= numUnits
        error('normalizationTerm_input must be empty or have length equal to numUnits (%d)', numUnits);
    end
    use_custom_norm = true;
    normalizationTerm = normalizationTerm_input(:); % Ensure column vector
else
    normalizationTerm = nans(numUnits,1);
end

%% Initialize output matrices
spike_rates_avg = nans(numUnits, numBins, 8);
spike_rates_SEM = nans(numUnits, numBins, 8);
spike_rates_avg_std_norm_avg = nans(numUnits, numBins, 8);
spike_rates_avg_std_norm_SEM = nans(numUnits, numBins, 8);

ROC_area = nans(numUnits, numBins, 4);

cell_selectivity = nans(numUnits,1);

stim_sz = nans(numUnits,1);

unit_example = [];

%% Loop through units
for uu = 1:numUnits
    
    % Show the name
    disp(uniqueUnitNames{uu})

    % Get the data from this session
    Lses = strcmp(dat.ses_ID, uniqueUnitNames{uu}) & dat.coh_final >= 50; 

    if sum(Lses) == 0
        continue
    end    

    session_data = dat(Lses, 1:end);
    stim_sz(uu,1) = session_data.dot_diam(1);

    dirs = [min(session_data.dir_final) max(session_data.dir_final)];

    % Collect spike bins for 2x2x2 conditions:
    %   Switch frequency (low/high)
    %   Switch trials (yes/no)
    %   Start direction (PREF/NULL)

    [SF, ST, SD] = ndgrid(sfs, [0 1], dirs);
    unit = struct( ...
        'dirs',             dirs, ... % Below will arrange as [PREF, NULL]
        'switch_frequency', num2cell(SF(:)), ...
        'switch_trial',     num2cell(ST(:)), ...
        'start_dir',        num2cell(SD(:)), ...
        'spike_counts',     []);

    dir_rates = [0 0];      % To determine PREF/NULL
    dir_rates_std = [0,0];
    maxRate = nans(1,length(unit));   % To determine normalization factor
    
    % Loop through each condition and collect spike data
    for cc = 1:length(unit)
        
        % Get the appropriate trials
        trial_indices = find(session_data.HR == unit(cc).switch_frequency & ...
            session_data.dir_switch == unit(cc).switch_trial & ...
            session_data.dir_prefinal == unit(cc).start_dir);

        unit(cc).trial_indices = trial_indices;
        num_trials = length(trial_indices);
        
        % Add matrix for binned spike counts and spike rates
        unit(cc).spike_counts = nans(num_trials, numBins);
        unit(cc).spike_rates = nans(num_trials, numBins);
        unit(cc).spike_times = nans(num_trials, 1000);

        % Get the binned spike counts
        for tt = 1:num_trials
            
            % Get spikes wrt dots test epoch start
            test_start = session_data.dots_on{trial_indices(tt)}(end);
            spikes = session_data.Unit_1{trial_indices(tt)} - test_start;
            unit(cc).spike_times(tt,1:size(spikes',2)) = spikes';
            
            % Count in bins
            end_time = session_data.dots_off(trial_indices(tt)) - test_start;

            for bb = 1:find(bins(:,2) > end_time,1)-1
                unit(cc).spike_counts(tt,bb) = sum(spikes >= bins(bb,1) & spikes < bins(bb,2));                
            end

            % Convert to rates (sp/s)
            unit(cc).spike_rates_raw(tt,:) = unit(cc).spike_counts(tt,:) * 1000/bin_size;

            % Standardize firing rates
            baseline_bin = find(bins(:,2) > -2500,1)-1;
            dots_on_bin = find(bins(:,2) > -2400,1)-1;
    
            baseline_activity = nanmean(unit(cc).spike_rates_raw(tt,baseline_bin:dots_on_bin));
            unit(cc).spike_rates_std(tt,:) = unit(cc).spike_rates_raw(tt,:) - baseline_activity;

            % Amend spike rates to be baseline subtracted
            unit(cc).spike_rates(tt,:) =  unit(cc).spike_rates_std(tt,:); 

        end

        % Calculate normalization factor (only if not using custom)
        if ~use_custom_norm
            maxRate(1,cc) = max(nanmean(unit(cc).spike_rates_std(:,Ladapt)));
        end
    end

    % Set normalization term
    if ~use_custom_norm
        normalizationTerm(uu,1) = max(maxRate);
    end

    for cc = 1:length(unit)

        unit(cc).spike_rates_std_norm = unit(cc).spike_rates_std/normalizationTerm(uu,1);

        % Calculate SEMs
        unit(cc).spike_rates_raw_SEM = nanse(unit(cc).spike_rates_raw,1);
        unit(cc).spike_rates_SEM = nanse(unit(cc).spike_rates,1);
        unit(cc).spike_rates_std_norm_SEM = nanse(unit(cc).spike_rates_std_norm,1);
        
        % Update mean dir1/dir2 rates by taking a weighted average of all relevant bins
        dir_rates = dir_rates + [ ...
            nanmean(reshape(unit(cc).spike_rates(:, Ldots(cc,:)&Ladapt),[],1)) .* size(unit(cc).spike_rates,1)./sum(Lses), ...
            nanmean(reshape(unit(cc).spike_rates(:,~Ldots(cc,:)&Ladapt),[],1)) .* size(unit(cc).spike_rates,1)./sum(Lses)];

        % Get difference between motion directions in standard units for across-session comparison
        dir_rates_std = dir_rates_std + [ ...
            nanmean(reshape(unit(cc).spike_rates_std_norm(:, Ldots(cc,:)&Ladapt),[],1)) .* size(unit(cc).spike_rates_std_norm,1)./sum(Lses), ...
            nanmean(reshape(unit(cc).spike_rates_std_norm(:,~Ldots(cc,:)&Ladapt),[],1)) .* size(unit(cc).spike_rates_std_norm,1)./sum(Lses)];
    end

    % Determine PREF/NULL
    if dir_rates(2) > dir_rates(1)
        unit_pref = [uniqueUnitNames{uu}, ' = ' num2str(unit(1).dirs(2)) ' PREF'];
        fprintf('%s\n', unit_pref)   
        [unit.dirs] = deal(unit(1).dirs([2 1]));
    elseif dir_rates(1) > dir_rates(2)
        unit_pref = [uniqueUnitNames{uu}, ' = ' num2str(unit(1).dirs(1)) ' PREF'];
        fprintf('%s\n', unit_pref)  
    end

    %% Save trialwise data for example neuron rasterplot
    if ~isempty(example_unit_idx) && uu == example_unit_idx
        unit_example = unit;
    end

    %% Save trialwise data
    unit_sort = unit;
    pp = 1;
    for ss = 1:2 % switch/no switch
        for ff = 1:2 % switch frequency
            for dd = 1:2
                Lcond = [unit.switch_trial] == 2 - ss & ...
                    [unit.switch_frequency] == sfs(ff) & ...
                    [unit.start_dir] == unit(1).dirs(dd);

                % Save the average firing rate and std. norm. firing rate + SEMs
                spike_rates_avg(uu,:,pp) = nanmean(unit_sort(Lcond).spike_rates);
                spike_rates_SEM(uu,:,pp) = unit_sort(Lcond).spike_rates_SEM;
                spike_rates_avg_std_norm_avg(uu,:,pp) = nanmean(unit_sort(Lcond).spike_rates_std_norm);
                spike_rates_avg_std_norm_SEM(uu,:,pp) = unit_sort(Lcond).spike_rates_std_norm_SEM;

                % Increment index
                pp = pp + 1;
            end
        end
    end

    %% Calculate selectivity
    % Using first 400 ms of PREF stimulus onset (PREF start) + 100 ms latency
    bin1 = find(bins(:,1) == -2300);
    bin2 = find(bins(:,2) == -1900);

    pref_resp = nanmean([spike_rates_avg(uu,bin1:bin2,1), spike_rates_avg(uu,bin1:bin2,3), ...
        spike_rates_avg(uu,bin1:bin2,5), spike_rates_avg(uu,bin1:bin2,7)]);

    null_resp = nanmean([spike_rates_avg(uu,bin1:bin2,2), spike_rates_avg(uu,bin1:bin2,4), ...
        spike_rates_avg(uu,bin1:bin2,6), spike_rates_avg(uu,bin1:bin2,8)]);
        
    cell_selectivity(uu,1) = pref_resp - null_resp;

    %% Calculate ROC areas
    for ss = 1:2 % switch/no switch
        for ff = 1:2 % switch frequency

            for dd = 1:2 % start dir
                Lcond = [unit.switch_trial] == 2 - ss & ...
                    [unit.switch_frequency] == sfs(ff) & ...
                    [unit.start_dir] == unit(1).dirs(dd);

                if ss == 1 && dd == 1
                    vals_pref = unit(Lcond).spike_rates;
                elseif ss == 2 && dd == 1
                    vals_null = unit(Lcond).spike_rates;
                elseif ss == 1 && dd == 2
                    vals_null = unit(Lcond).spike_rates;
                elseif ss == 2 && dd == 2
                    vals_pref = unit(Lcond).spike_rates;
                end
            end

            for b = 1:numBins
                x = vals_pref(:,b);
                y = vals_null(:,b);

                if sum(~isnan(y)) > 0 && sum(~isnan(x)) > 0
                    N = 100;
                    plotflag = 0;
                    [~,~,z] = rocNFullOutput(x,y,N,plotflag); 
                    a(ff,b) = z;
                else
                    a(ff,b) = NaN;
                end
            end

            % Save ROC area matrices:
            %   1. LSF non-switch trials
            %   2. LSF switch trials
            %   3. HSF switch trials
            %   4. HSF non-switch trials
            if ff == 1 && ss == 2
                ROC_area(uu,:,1) = a(ff,:);
            elseif ff == 1 && ss == 1
                ROC_area(uu,:,2) = a(ff,:);
            elseif ff == 2 && ss == 1
                ROC_area(uu,:,3) = a(ff,:);
            elseif ff == 2 && ss == 2
                ROC_area(uu,:,4) = a(ff,:);   
            end
        end
    end
end

end
