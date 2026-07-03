%% Fig4_FigSupp2.m
% Analyzes adaptation in MT neurons by comparing responses to the first
% preferred-direction stimulus vs. the first null-adapted preferred-direction
% stimulus during the adaptation epoch.
%
% Analysis approach:
%   - Processes neural data using processNeuralData.m for correct trials
%   - Classifies neurons as "adapting" or "facilitating" based on change
%     from first stimulus to test-stimulus response:
%       * Adapting: Test response decreases >2.5% from first stimulus
%       * Facilitating: Test response increases >2.5% from first stimulus
%       * Must meet criterion at both LSF and HSF to be classified
%   - Compares first PREF stimulus vs. first null-adapted PREF stimulus
%
% Time epochs analyzed:
%   LSF condition (cond = 1):
%     - First PREF: -2200 to -2000ms
%     - First null-adapted PREF: -1000 to -800ms (after 1200ms null motion)
%
%   HSF condition (cond = 2):
%     - First PREF: -2200 to -2000ms
%     - First null-adapted PREF: -1800 to -1600ms (after 400ms null motion)
%
% Neuron subsets (controlled by 'subset' parameter):
%   - adapt: Neurons classified as adapting at both LSF and HSF
%   - facilitate: Neurons classified as facilitating at both LSF and HSF
%   - all_cells: All recorded neurons (default)
%
% Outputs:
%   Figure 4-Figure Supplement 2: Scatter plot comparing first PREF stimulus
%   response vs. null-adapted PREF stimulus response. Each point is one
%   neuron. Separate markers for each monkey:
%     - Monkey An: circles
%     - Monkey Ch: diamonds
%     - Monkey Mi: squares
%
%   Statistics: Wilcoxon signed-rank test (paired) and Cohen's d effect size
%
% Required data files:
%   - mergedTable_proc.mat
%   - sensitivity_diff_labeled_N.mat (neural subset classification)
%
% Required functions:
%   - processNeuralData.m (extract MT activity)
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

% Create monkey indices for subsetting units.
% Numeric (not logical) indices are needed here because later code
% composes them with a second index, e.g. An(subset(An)).
monkeyIdx = getMonkeyIndices(dat);
An = find(monkeyIdx.An);
Ch = find(monkeyIdx.Ch);
Mi = find(monkeyIdx.Mi);

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

%%
% Process neural data (calling processNeuralData.m)

% Output spike matrices organization:
    %   1. Low switch frequency switch (PREF)
    %   2. Low switch frequency switch (NULL)
    %   3. High switch frequency switch (PREF)
    %   4. High switch frequency switch (NULL)
    %   5. Low switch frequency non-switch (NULL) % testing epoch
    %   6. Low switch frequency non-switch (PREF)
    %   7. High switch frequency non-switch (NULL)
    %   8. High switch frequency non-switch (PREF)

[spike_rates_avg, spike_rates_SEM, spike_rates_avg_std_norm_avg, ...
    spike_rates_avg_std_norm_SEM, ...
    ROC_area, cell_selectivity, ...
    normalizationTerm, uniqueUnitNames, ...
    unit_example] = processNeuralData(dat, slide, bin_size, start_time, end_time, [], []);

% Assign "adapting" vs. "facilitating" groups:

% Bin range for first adapting stimulus
bin1 = -2200;
bin2 = -2000;

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

% Bin range for test-stimulus
bin3 = 200;
bin4 = 400;

test_start_bin = find(bins(:,1) == bin3 - 0.5*bin_size,1);
test_end_bin = find(bins(:,1) == bin4 - 0.5*bin_size,1);

% Calculate "adapting" vs. "facilitating" at HSF
% Thresholds are absolute differences in spike_rates_avg_std_norm_avg,
% which is already normalized to each unit's own max firing rate, so
% 0.025 corresponds to 2.5 percentage points of that per-unit max (not a
% 2.5% relative change in the response itself).
x = mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,3), 2, 'omitnan');
y = mean(spike_rates_avg_std_norm_avg(:,test_start_bin:test_end_bin,3), 2, 'omitnan');

adapt_H = y < x - 0.025;       % Adapting if response drops by > 2.5 percentage points
facilitate_H = y >= x + 0.025; % Facilitating if response rises by > 2.5 percentage points

% Calculate "adapting" vs. "facilitating" at LSF
x = mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,1), 2, 'omitnan');
y = mean(spike_rates_avg_std_norm_avg(:,test_start_bin:test_end_bin,1), 2, 'omitnan');

adapt_L = y < x - 0.025;
facilitate_L = y >= x + 0.025;

% Classify cells that meet "adapting" or "facilitating" criteria at both LSF and HSF
adapt = adapt_H + adapt_L == 2;
facilitate = facilitate_H + facilitate_L == 2;

all_cells = true(numUnits,1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 4-Figure Supplement 2 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Comparison of first preferred-motion stimulus response vs. null-adapted
% preferred-motion response. Published figure has 6 panels: one row per
% neuron subset (All, Adapting, Facilitating) and one column per
% switch-frequency condition (High, Low) -- all computed and plotted here
% in a single pass.

% First preferred-motion stimulus presentation (same for both conditions):
bin_f_stim1 = -2200;
bin_f_stim2 = -2000;

start_bin_f_stim = find(bins(:,1) == bin_f_stim1 - 0.5*bin_size,1);
end_bin_f_stim = find(bins(:,1) == bin_f_stim2 - 0.5*bin_size,1);

subsetGroups = {all_cells, adapt, facilitate};
subsetNames = {'All MT neurons', '"Adapting" neurons', '"Facilitating" neurons'};

condOrder = [2 1]; % Column order: High switch frequency, then Low switch frequency
condNames = {'Low switch frequency', 'High switch frequency'}; % indexed by cond (1=LSF, 2=HSF)

x_line = -0.2:0.1:2;
y_line = -0.2:0.1:2;

figure(6); clf

for rr = 1:numel(subsetGroups)
    subset = subsetGroups{rr};

    for cc = 1:numel(condOrder)
        cond = condOrder(cc);

        if cond == 2

            % High switch frequency
            % Second preferred-motion (first null-adapted) stimulus presentation:
            bin_s_stim1 = -1800;
            bin_s_stim2 = -1600;

            % Switch-trial PREF/NULL pair for HSF (processNeuralData.m dim-3 cols 3,4)
            idx_1 = 3;
            idx_2 = 4;

        elseif cond == 1

            % Low switch frequency
            % Second preferred-motion (first null-adapted) stimulus presentation:
            bin_s_stim1 = -1000;
            bin_s_stim2 = -800;

            % Switch-trial PREF/NULL pair for LSF (processNeuralData.m dim-3 cols 1,2)
            idx_1 = 1;
            idx_2 = 2;
        end

        start_bin_s_stim = find(bins(:,1) == bin_s_stim1 - 0.5*bin_size,1);
        end_bin_s_stim = find(bins(:,1) == bin_s_stim2 - 0.5*bin_size,1);

        subplot(numel(subsetGroups), numel(condOrder), (rr-1)*numel(condOrder)+cc); hold on; box on;

        plot(x_line, y_line, 'k');

        x = mean(spike_rates_avg_std_norm_avg(An(subset(An)),start_bin_f_stim:end_bin_f_stim,idx_1), 2, 'omitnan');
        y = mean(spike_rates_avg_std_norm_avg(An(subset(An)),start_bin_s_stim:end_bin_s_stim,idx_2), 2, 'omitnan');
        plot(x,y, 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

        x = mean(spike_rates_avg_std_norm_avg(Ch(subset(Ch)),start_bin_f_stim:end_bin_f_stim,idx_1), 2, 'omitnan');
        y = mean(spike_rates_avg_std_norm_avg(Ch(subset(Ch)),start_bin_s_stim:end_bin_s_stim,idx_2), 2, 'omitnan');
        plot(x,y, 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

        x = mean(spike_rates_avg_std_norm_avg(Mi(subset(Mi)),start_bin_f_stim:end_bin_f_stim,idx_1), 2, 'omitnan');
        y = mean(spike_rates_avg_std_norm_avg(Mi(subset(Mi)),start_bin_s_stim:end_bin_s_stim,idx_2), 2, 'omitnan');
        plot(x,y,'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

        x = mean(spike_rates_avg_std_norm_avg(subset,start_bin_f_stim:end_bin_f_stim,idx_1), 2, 'omitnan');
        y = mean(spike_rates_avg_std_norm_avg(subset,start_bin_s_stim:end_bin_s_stim,idx_2), 2, 'omitnan');

        [p, h, stats] = signrank(x,y);
        d_val = computeCohenDCI(y, x, 'paired');

        xlabel('First PREF stim. avg. FR')
        ylabel('NULL-adapted PREF stim. avg. FR')

        title(sprintf('%s (n=%d): %s', subsetNames{rr}, sum(subset), condNames{cond}))
        subtitle(['p = ', num2str(p), ' & d = ', num2str(d_val)])

        ylim([-0.2,1])
        xlim([-0.2,1])
    end
end