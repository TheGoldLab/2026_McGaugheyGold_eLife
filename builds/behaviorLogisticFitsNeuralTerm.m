%% behaviorLogisticFitsNeuralTerm.m
% Fits behavioral data with a time-dependent logistic psychometric function
% that includes a trial-wise MT neural activity term.
%
% Model structure:
%   P(switch) = logistic(β₀ + β₁*dir_bias + β₂*neural×signed_time + β₃*signed_time)
%   where neural term is the average MT firing rate during the test epoch
%   (50-1200 ms after test onset), normalized and averaged across neurons.
%
% Control condition (control = true):
%   Fits standard behavior-only model restricted to PREF direction trials:
%   P(switch) = logistic(β₀ + β₁*dir_bias + β₂*signed_time)
%
% Neural processing:
%   - Bins spike counts in 100ms windows with 10ms slide
%   - Baseline-subtracts using -2500 to -2400ms epoch (100ms before dot-motion onset)
%   - Normalizes by max firing rate during adapting epoch
%   - Determines PREF/NULL directio
%   - Extracts average test-epoch firing rate for each trial
%   - Only PREF direction trials contribute to neural×time interaction
%
% Fits are performed separately for:
%   - Each recording session (neuron)
%   - Each switch-frequency condition (LSF=2, HSF=6)
%
% Outputs:
%   fits - [numUnits x 5 parameters x 2 frequencies] Neural model coefficients
%          β₀: Switch/stay bias
%          β₁: Right/left directional bias
%          β₂: Neural×signed-time interaction 
%          β₃: Sensitivity to viewing duration
%          β₄: Lapse rate
%   fits_control - [numUnits x 4 parameters x 2 frequencies] Control model
%   R_sq - [numUnits x 2 frequencies] Tjur's pseudo-R² for neural model
%   R_sq_control - [numUnits x 2 frequencies] Tjur's pseudo-R² for control
%
% Note: Fitting is computationally intensive. Pre-computed fits are saved in
% data/behaviorFits/ as LogisticFits_NeuralTerm.mat and LogisticFits_NeuralTerm_control.mat
%
% Required functions: 
%   logistFitDotsrevNP.m, logistValDotsrevNP.m, logistErrDotsrev.m (neural model)
%   logistFitDotsrev.m, logistValDotsrev.m, logistErrDotsrev.m  (control model)

%%
% Set paths

cfg = projectDefaults();
cd(cfg.paths.data)

% Load data:

load('mergedTable_proc_neural.mat') % keeps Unit_1, skips pupil traces (see buildMergedTableTiers.m)

% Subset data appropriately

[mergedTableSub] = createDatSubset(mergedTable_proc, 'NP');
dat = mergedTableSub;

clear mergedTable_proc

%% Choose condition to run

control = false; % True will run the behavior-only logistic with no neural term
                 % False will run a logistic that includes a neural term
     
%%
% Get unique session (unit) names:

uniqueUnitNames = unique(dat.ses_ID);
numUnits = length(uniqueUnitNames);

% Set up timing/binning:

slide = 10;     % msec
bin_size = 100; % msec
start_time = -2600;
end_time = 1200;
bins = cat(2, ...
    (start_time:slide:end_time - bin_size)', ...
    (start_time + bin_size:slide:end_time)');
xax = mean(bins,2);
numBins = size(bins, 1);
axlen = length(xax);

% Make selection array to compute PREF/NULL direction
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

%%
% Save the fits

fits = nans(numUnits, 5, 2);
fits_control = nans(numUnits, 4, 2);

% Save the R_sq values

R_sq = nans(numUnits, 1, 2);
R_sq_control = nans(numUnits, 1, 2);

% Plotting controls:

colors = cfg.colors.pair;
line_styles = {'-', ':'};

%%
% Loop through the units:
for uu = 1:numUnits
    
    % Show the name

    disp(uniqueUnitNames{uu})

    % Monkey Ch is excluded from this per-unit neural-term/control fit
    % (also excluded in behaviorLogisticFitsPupilTerm.m and
    % Fig4_FigSupp1B.m's block-split analysis); reason not documented in
    % the original code. Ch is included in the plain behavioral scripts.
    if uniqueUnitNames{uu}(1:2) == "Ch"
        continue
    end

    % Get the data from this session)

    Lses = strcmp(dat.ses_ID, uniqueUnitNames{uu}) & dat.coh_final >= 50; 

    if sum(Lses) == 0
        continue
    end    

    session_data_tbl = dat(Lses, 1:end);

    num_trials = sum(Lses);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% NEURAL %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Collect spike bins for 2x2x2 conditions:
        %   Switch frequency (low/high)
        %   Switch trials (yes/no)
        %   Start direction (PREF/NULL)

    dirs = [min(session_data_tbl.dir_final) max(session_data_tbl.dir_final)];

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

        trial_indices = find(session_data_tbl.HR == unit(cc).switch_frequency & ...
            session_data_tbl.dir_switch == unit(cc).switch_trial & ...
            session_data_tbl.dir_prefinal == unit(cc).start_dir);

        unit(cc).trial_indices = trial_indices;

        num_trials_sub = length(trial_indices);
        
        % Add matrix for binned spike counts and spike rates

        unit(cc).spike_counts = nans(num_trials_sub, numBins);
        unit(cc).spike_rates = nans(num_trials_sub, numBins);
        unit(cc).spike_times = nans(num_trials_sub, 1000);

        % Get the binned spike counts

        for tt = 1:num_trials_sub
            
            % Get spikes wrt dots test epoch start

            test_start = session_data_tbl.dots_on{trial_indices(tt)}(end);
            spikes = session_data_tbl.Unit_1{trial_indices(tt)} - test_start;
            unit(cc).spike_times(tt,1:size(spikes',2)) = spikes';
            
            % Count in bins

            end_time = session_data_tbl.dots_off(trial_indices(tt)) - test_start;

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

            % Store spike rates as baseline subtracted

            unit(cc).spike_rates(tt,:) =  unit(cc).spike_rates_std(tt,:); 

        end

        % Normalize firing rates
        % Normalizing wrt to Adapting Epoch to avoid noise in Testing Epoch which often takes the form of high FR on long trials

        maxRate(1,cc) = max(nanmean(unit(cc).spike_rates_std(:,Ladapt)));

    end

    normalizationTerm(uu,1) = max(maxRate);

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

    % Find testing stimulus average FR (std. norm.) for each trial

    bin1 = 50; % Want to average neurao data starting 50 ms after test-stim onset
    start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);

    FR_testing_stim = nans(num_trials,1);

    trial_idx = 1;

    for ss = 1:size(unit,1)

        FR_stdNorm_sz = size(unit(ss).spike_rates_std_norm,1);

        for ff = 1:FR_stdNorm_sz

            % Store FR and trial index for matching back in session matrix

            FR_testing_stim(trial_idx,1) = nanmean(unit(ss).spike_rates_std_norm(ff,start_bin:end),2);
            FR_testing_stim_idx(trial_idx,1) = unit(ss).trial_indices(ff);

            % Capture PREF trials

            % PREF start/switch trials

            if unit(ss).dirs(1) == unit(ss).start_dir && unit(ss).switch_trial == 1
                FR_testing_stim_dir_idx(trial_idx,1) = 1;

            % NULL start/switch
            
            elseif unit(ss).dirs(2) == unit(ss).start_dir && unit(ss).switch_trial == 0
                FR_testing_stim_dir_idx(trial_idx,1) = 1;

            else
                FR_testing_stim_dir_idx(trial_idx,1) = 0;
            end

            trial_idx = trial_idx+1;
        end
    end

    FR_testing_stim_comb = [FR_testing_stim_idx,FR_testing_stim,FR_testing_stim_dir_idx];

    % Sort rows to match session matrix

    FR_testing_stim_comb = sortrows(FR_testing_stim_comb,1);

    % Append to session_data_tbl

    FR_testing = FR_testing_stim_comb(:,2);
    FR_dir_idx = FR_testing_stim_comb(:,3);

    Ldir = logical(FR_dir_idx);

    session_data_tbl = addvars(session_data_tbl,FR_testing,'After',"Unit_1");
    session_data_tbl = addvars(session_data_tbl,FR_dir_idx,'After',"FR_testing");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% BEHAVIOR %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % Get the data to fit:
    %   1. hazard
    %   2. test coherence
    %   3. prefinal direction
    %   4. final direction
    %   5. test duration
    %   6. choice (1=right, 2=left)
    %   7. Neural term

    % NOTE: indexed by name, not position -- session_data_tbl's column
    % count/order depends on which mergedTable_proc tier was loaded
    % (see buildMergedTableTiers.m), and FR_testing's absolute position
    % shifts accordingly since it's appended via addvars() above.
    session_data_fit = table2array(session_data_tbl(:, ...
        {'HR', 'coh_final', 'dir_prefinal', 'dir_final', 'dur_final', 'choice_final', 'FR_testing'}));
    uniqueHazards = cfg.hazard.codes;

    % Set up data matrix -- see buildBehaviorDataToFit.m /
    % buildInteractionDataToFit.m for column definitions. Both the control
    % and neural-term models here restrict the signed-time/interaction
    % term to PREF-direction trials via FR_dir_idx (unlike
    % behaviorLogisticFitsPupilTerm.m, which uses all trials).
    if control
        [data_to_fit, Lpre] = buildBehaviorDataToFit(session_data_fit, FR_dir_idx);
    else
        [data_to_fit, Lpre] = buildInteractionDataToFit(session_data_fit, session_data_fit(:,7), FR_dir_idx);
    end

    % Fit/plot separately for each hazard

    clf % clear figure

    for hh = 1:2
        
        %Lneural = logical(neural_split(:,hh));

        Lhazard = session_data_fit(:,1) == uniqueHazards(hh);
        
        % Get the fit 

        if control
            [fits_control(uu,:,hh), LL(uu,:,hh)] = logistFitDotsrev(data_to_fit(Lhazard,:));
        else
            [fits(uu,:,hh), LL(uu,:,hh)] = logistFitDotsrevNP(data_to_fit(Lhazard,:));
        end

        % Compute Tjur's D (Tjur's coefficient of disctrimination -- R^2)
             % Compares average fitted probability of the two response outcomes (0 and 1)
             % Calculate mean of predicted probabilities of each outcome and take the difference between the 2 means
             % 0 <= R^2 >= 1 where no disriminating power = 0 and perfect discriminating power = 1

         if control
            fit_y = logistValDotsrev(fits_control(uu,:,hh)', data_to_fit(Lhazard,1:end-1));
            R_sq_control(uu,:,hh) = computeTjursR2(fit_y, data_to_fit(Lhazard,4));
         else
             fit_y = logistValDotsrevNP(fits(uu,:,hh)', data_to_fit(Lhazard,1:end-1));
             R_sq(uu,:,hh) = computeTjursR2(fit_y, data_to_fit(Lhazard,5));
         end

        % Plot separately for predir=L and predir=R

        for xx = 1:2
            figure(2)
            subplot(1,2,xx); hold on; box on;
            plot([0 0], [0 1], 'k--');
            plot(xax([1 end]), [0.5 0.5], 'k--', 'MarkerSize', 15);

            if control
                choice_data = data_to_fit(Lhazard&Lpre(:,xx), 4);
                
                % Show smoothed data
                [sorted_time_axis, I] = sort(data_to_fit(Lhazard&Lpre(:,xx), 3));
            else
                choice_data = data_to_fit(Lhazard&Lpre(:,xx), 5);

                % Show smoothed data
                [sorted_time_axis, I] = sort(data_to_fit(Lhazard&Lpre(:,xx), 4));
            end

            sorted_choice_data = choice_data(I);
            plot(sorted_time_axis, nanrunmean(sorted_choice_data,10), ':', 'Color', colors{hh});

            hold on;
            
            if xx==1
                title(sprintf('Prefinal dir=R'))
                ylabel('Fraction choose Left')
            else
                title(sprintf('Prefinal dir=L'))
                ylabel('Fraction choose Right')

            end
            xlabel('Signed time (-non_switch, +switch)');

            % Show fits

            xax_fits = (-1200:1200)' ;
            show_fit_data = repmat(cat(2, ones(length(xax_fits),2), xax_fits, xax_fits),[1 1 2]);
            show_fit_data(:,2,2) = 0;

            if control
                ys = logistValDotsrev(fits_control(uu,:,hh)', show_fit_data(:,:,xx));

            else
                ys = logistValDotsrevNP(fits(uu,:,hh)', show_fit_data(:,:,xx));
            end

            plot(xax_fits, ys, '-', 'Color', colors{hh}, 'LineWidth', 2);
            axis([xax_fits(1) xax_fits(end) 0 1]);

            drawnow

            hold on
        end     
    end

    if control
        fits_neuralTerm_control(uu,1,:) = fits_control(uu,3,:);
        R_sq_neuralTerm_control(uu,1,:) = R_sq_control(uu,1,:);
    else
        fits_neuralTerm(uu,1,:) = fits(uu,3,:);
        R_sq_neuralTerm(uu,1,:) = R_sq(uu,1,:);
    end

    clear FR_testing_stim_idx FR_testing_stim FR_testing_stim_dir_idx FR_testing FR_dir_idx 
end