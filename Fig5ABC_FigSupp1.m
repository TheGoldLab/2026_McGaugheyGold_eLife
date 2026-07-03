%% Fig5ABC_FigSupp1.m
% Analyzes context-dependent changes in MT discriminability (ROC area)
% and their relationship to behavior across switch-frequency conditions.
%
% ROC area analysis:
%   ROC area quantifies how well single-trial responses discriminate between
%   preferred and null motion directions. Values near 0.5 = no discrimination.
%   Computed separately for LSF and HSF switch trials during test epoch.
%
% Behavioral analysis:
%   Uses behaviorBinnedPerformance.m to calculate performance (fraction correct)
%   in 4 viewing duration bins (100-225, 225-375, 375-600, 600-1200ms), separately
%   for LSF and HSF. Focuses on bin 3 (375-600ms) where differences are maximal.
%
% Main Text Figures:
%
%   Figure 5A: Example neuron (unit 151) showing time course of ROC area for
%     LSF (blue) and HSF (yellow) conditions. X-axis: Time relative to test
%     onset. Y-axis: ROC area (0.5-1.0). Demonstrates individual neuron shows
%     higher selectivity at LSF than HSF.
%
%   Figure 5B: Population comparison of test-epoch ROC area (50-500ms). Scatter
%     plot comparing HSF vs LSF average ROC area. Each point is one neuron.
%     Points above unity line indicate higher discriminability at LSF. Different
%     markers by monkey (circles=An, diamonds=Ch, squares=Mi, N=155 total).
%
%   Figure 5C: Neural-behavioral relationship. Scatter plot showing correlation
%     between neural selectivity difference (LSF-HSF ROC area, 200-400ms) and
%     behavioral performance difference (LSF-HSF accuracy in 375-600ms bin).
%     Each point is one session. Positive correlation indicates that sessions
%     with larger neural selectivity differences show correspondingly larger
%     behavioral performance differences. Includes linear regression fit.
%     Different markers by monkey. 
%
% Figure Supplements:
%
%   Figure 5-Figure Supplement 1A: Same as Figure 5B but separated by individual
%     monkeys in three panels.
%
%   Figure 5-Figure Supplement 1B: Same as Figure 5C but separated by individual
%     monkeys. Three panels showing regression fits and correlation
%     statistics for each monkey separately.
%
% Output organization (from processNeuralData.m):
%   Spike rates: 8 conditions (LSF/HSF × switch/non-switch × PREF/NULL)
%   ROC area: 4 conditions
%     1. LSF non-switch trials
%     2. LSF switch trials  
%     3. HSF switch trials
%     4. HSF non-switch trials
%
% Required data files:
%   - mergedTable_proc.mat
%   - sensitivity_diff_labeled_N.mat
%
% Required functions:
%   - processNeuralData.m (compute ROC area for selectivity)
%   - behaviorBinnedPerformance.m (compute binned performance)
%   - computeCohenDCI.m (calculate effect size)

%% Load data
cfg = projectDefaults();
cd(cfg.paths.data)

load('mergedTable_proc_neural.mat') % keeps Unit_1, skips pupil traces (see buildMergedTableTiers.m)
load('sensitivity_diff_labeled_N.mat') % Neural subset

% Subset data appropriately
% "N" = Neural analysis
[mergedTableSub] = createDatSubset(mergedTable_proc, 'N');
dat = mergedTableSub;

clear mergedTable_proc

% Select correct trials
datCorrect = dat(dat.correct == 1,:);
dat = datCorrect;

% Create monkey indices for subsetting units:
monkeyIdx = getMonkeyIndices(dat);
An = monkeyIdx.An;
Ch = monkeyIdx.Ch;
Mi = monkeyIdx.Mi;

monkeyOrder = {'An', 'Mi', 'Ch'};
monkeyMarkers = struct('An', 'o', 'Mi', 'square', 'Ch', 'diamond');

% Get unique session (unit) names:
uniqueUnitNames = unique(dat.ses_ID);
numUnits = length(uniqueUnitNames);

% Unit plots -- plotting controls:
colors = cfg.colors.pair;

%%
% Process neural data (calling processNeuralData.m)

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

% Output spike matrices organization:
    %   1. Low switch frequency switch (PREF)
    %   2. Low switch frequency switch (NULL)
    %   3. High switch frequency switch (PREF)
    %   4. High switch frequency switch (NULL)
    %   5. Low switch frequency non-switch (NULL) % testing epoch
    %   6. Low switch frequency non-switch (PREF)
    %   7. High switch frequency non-switch (NULL)
    %   8. High switch frequency non-switch (PREF)

% Output ROC area matrix:
    % 1. LSF non-switch trials
    % 2. LSF switch trials
    % 3. HSF switch trials
    % 4. HSF non-switch trials

[spike_rates_avg, spike_rates_SEM, spike_rates_avg_std_norm_avg, ...
    spike_rates_avg_std_norm_SEM, ...
    ROC_area, cell_selectivity, ...
    normalizationTerm, uniqueUnitNames, ...
    unit_example] = processNeuralData(dat, slide, bin_size, start_time, end_time, [], []);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5A %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example unit average ROC area as a function of time relative to test-stimulus onset

uu = 151;

figure(11); clf; hold on; box on;

plot(ROC_area(uu,:,2), 'Color', colors{1}, 'LineWidth', 2)
plot(ROC_area(uu,:,3), 'Color', colors{2}, 'LineWidth', 2)

xlabel('Time relative to test onset (ms)')
ylabel('Average ROC area')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5B %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin1 = 50; % Starting bin for ROC area average
bin2 = 500; % Ending bin for ROC area average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

figure(1); clf; hold on; box on;

x_line = -0.4:0.1:1;
y_line = -0.4:0.1:1;

xlim([0.2 1])
ylim([0.2 1])

plot(x_line, y_line, 'k');

x = mean(ROC_area(Mi,start_bin:end_bin,3), 2, 'omitnan');
y = mean(ROC_area(Mi,start_bin:end_bin,2), 2, 'omitnan');
plot(x,y, 'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

x = mean(ROC_area(An,start_bin:end_bin,3), 2, 'omitnan');
y = mean(ROC_area(An,start_bin:end_bin,2), 2, 'omitnan');
plot(x,y, 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

x = mean(ROC_area(Ch,start_bin:end_bin,3), 2, 'omitnan');
y = mean(ROC_area(Ch,start_bin:end_bin,2), 2, 'omitnan');
plot(x,y, 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

x = mean(ROC_area(:,start_bin:end_bin,3), 2, 'omitnan');
y = mean(ROC_area(:,start_bin:end_bin,2), 2, 'omitnan');

xlabel('ROC area avg. (HSF)')
ylabel('ROC area avg. (LSF)')

[p, h, stats] = signrank(x,y);
d = computeCohenDCI(y, x, 'paired');

title(['Testing Epoch average ROC area ', num2str(bin1), '-', num2str(bin2)], ' ms')
subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5C %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate binned performance difference (LSF-HSF)
% Returns vector of 4 bins:
    % Bin 1: 100-225 ms
    % Bin 2: 225-375 ms
    % Bin 3: 375-600 ms
    % Bin 4: 600-1200 ms

% For performance metric need to change dat to using all trials, not just correct trials
dat = mergedTableSub;

[LSF_switch_binned_performance, HSF_switch_binned_performance, binned_behavior_diff] = ...
    behaviorBinnedPerformance(dat);

% plot controls:
bin1 = 200; % Starting bin for ROC area average
bin2 = 400; % Ending bin for ROC area average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

ROCa_H = mean(ROC_area(:,start_bin:end_bin,3), 2, 'omitnan');
ROCa_L = mean(ROC_area(:,start_bin:end_bin,2), 2, 'omitnan');

ROCa_diff = ROCa_L - ROCa_H; % Difference in average ROC area (LSF-HSF)

figure(18); clf; hold on; box on;

plot(ROCa_diff(Mi), mean(binned_behavior_diff(Mi,3), 2, 'omitnan'), 'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')
plot(ROCa_diff(Ch), mean(binned_behavior_diff(Ch,3), 2, 'omitnan'), 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')
plot(ROCa_diff(An), mean(binned_behavior_diff(An,3), 2, 'omitnan'), 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

mdl = fitlm(ROCa_diff(:), mean(binned_behavior_diff(:,3), 2, 'omitnan'));
p = plot(mdl);
p(1).Visible = 'off';
legend('hide')
title(' ')

xlabel('Adaptation difference (LSF - HSF)')
ylabel('Binned behavioral difference (LSF - HSF)')

[r,p] = corr(ROCa_diff(:), mean(binned_behavior_diff(:,3), 2, 'omitnan'));
subtitle(['p = ', num2str(p), ' & r = ', num2str(r)])

%% Relevant Figure Supplements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5-Figure Supplement 1A %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Average ROC area HSF vs. LSF plotted separately for each animal
% Panel order (An, Mi, Ch) matches the manuscript, not
% monkeyIdx.names' alphabetical (An, Ch, Mi) order.

bin1 = 50; % Starting bin for ROC area average
bin2 = 500; % Ending bin for ROC area average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

x_line = -0.4:0.1:1;
y_line = -0.4:0.1:1;

figure(19); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    isMonkey = monkeyIdx.(monkeyName);

    subplot(1,3,mm); hold on; box on;

    x = mean(ROC_area(isMonkey,start_bin:end_bin,3), 2, 'omitnan');
    y = mean(ROC_area(isMonkey,start_bin:end_bin,2), 2, 'omitnan');

    plot(x_line, y_line, 'k');
    plot(x,y, 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    xlim([0.2 1])
    ylim([0.2 1])

    xlabel('ROC area avg. (HSF)')
    ylabel('ROC area avg. (LSF)')

    [p, h, stats] = signrank(x,y);
    d = computeCohenDCI(y, x, 'paired');

    title(['Testing Epoch average ROC area ', num2str(bin1), '-', num2str(bin2)], ' ms')
    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5-Figure Supplement 1B %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Relationship between context-dependent differences in behavior and context-dependent differences in ROC area
% Plotted separately for each animal

% Get binned performance:
dat = mergedTableSub;

[LSF_switch_binned_performance, HSF_switch_binned_performance, binned_behavior_diff] = ...
    behaviorBinnedPerformance(dat);

% plot controls:
bin1 = 200; % Starting bin for ROC area average
bin2 = 400; % Ending bin for ROC area average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

ROCa_H = mean(ROC_area(:,start_bin:end_bin,3), 2, 'omitnan');
ROCa_L = mean(ROC_area(:,start_bin:end_bin,2), 2, 'omitnan');

ROCa_diff = ROCa_L - ROCa_H; % Difference in average ROC area (LSF-HSF)

% Panel order (An, Mi, Ch) matches the manuscript, not
% monkeyIdx.names' alphabetical (An, Ch, Mi) order.
figure(89); clf;

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    isMonkey = monkeyIdx.(monkeyName);

    subplot(1,3,mm); hold on; box on;

    plot(ROCa_diff(isMonkey), mean(binned_behavior_diff(isMonkey,3), 2, 'omitnan'), 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    mdl = fitlm(ROCa_diff(isMonkey), mean(binned_behavior_diff(isMonkey,3), 2, 'omitnan'));
    p = plot(mdl);
    p(1).Visible = 'off';
    legend('hide')
    title(' ')

    xlabel('Adaptation difference (LSF - HSF)')
    ylabel('Binned behavioral difference (LSF - HSF)')

    [r,p] = corr(ROCa_diff(isMonkey), mean(binned_behavior_diff(isMonkey,3), 2, 'omitnan'));
    subtitle(['p = ', num2str(p), ' & r = ', num2str(r)])
end