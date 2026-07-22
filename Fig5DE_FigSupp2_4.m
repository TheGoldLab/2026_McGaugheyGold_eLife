%% Fig5DE_FigSupp2_4.m
% Analyzes relationships between MT neural encoding, discriminability,
% and behavioral performance, separated by sessions where monkeys showed greater
% behavioral sensitivity at LSF vs. HSF.
%
% Session classification:
%   Based on psychometric slope differences (from logistic fits):
%     - LmoresenLSF: Sessions with steeper slopes at LSF (slope_diff > 0)
%     - Lothersen: Sessions with steeper slopes at HSF or equal (slope_diff ≤ 0)
%
% Main Text Figures:
%
%   Figure 5D: Three-panel summary for LmoresenLSF sessions.
%     Panel 1: Population average firing rate time course (0-900ms)
%       - LSF (blue) and HSF (yellow) PREF direction responses
%       - Shows mean±SEM across neurons in these sessions
%       - LSF responses higher throughout test epoch
%
%     Panel 2: ROC area scatter plot (200-400ms average)
%       - HSF vs LSF directional selectivity
%       - Most points above unity line (higher discriminability at LSF)
%       - Different markers by monkey
%
%     Panel 3: Binned behavioral performance
%       - Fraction correct vs viewing duration (4 bins with midpoints: 162.5, 300, 487.5, 900ms)
%       - LSF (blue) and HSF (yellow) with error bars (SEM)
%       - LSF performance higher, especially at intermediate durations
%
%   Figure 5E: Same three-panel format for Lothersen sessions.
%     Tests whether sessions with opposite behavioral pattern show opposite
%     or absent neural patterns.
%
% Figure Supplements:
%
%   Figure 5-Figure Supplement 2: Test-epoch firing rates (200-400ms) separated by
%     monkey AND behavioral sensitivity classification. Six panels (3 monkeys x
%     2 sensitivity groups): row 1 = LmoresenLSF sessions, row 2 = Lothersen
%     sessions. Tests whether neural-behavioral correspondence holds within
%     individual subjects.
%
%   Figure 5-Figure Supplement 4: Exploration of test-stimulus encoding,
%     discriminability, and behavior for non-switch trials. Eight panels:
%     four columns (1: test-stimulus evidence encoding, 2: test-stimulus ROC
%     area 50-200 ms, 3: test-stimulus ROC area 200-400 ms, 4: test-stimulus
%     behavioral performance/fraction correct) by two rows, matching the
%     Figure 5D/E session-group split (row 1 = LmoresenLSF, n=100; row 2 =
%     Lothersen, n=32).
%
% Behavioral bins (from behaviorBinnedPerformance.m):
%   Bin 1: 100-225ms
%   Bin 2: 225-375ms
%   Bin 3: 375-600ms
%   Bin 4: 600-1200ms
%
% Required data files:
%   - mergedTable_proc.mat
%   - sensitivity_diff_labeled_N.mat (psychometric slope classifications)
%
% Required functions:
%   - processNeuralData.m (MT activity and ROC area)
%   - behaviorBinnedPerformance.m (binned accuracy)
%   - computeCohenDCI.m (effect size)

%% Load data
cfg = projectDefaults();

load(fullfile(cfg.paths.data, 'mergedTable_proc_neural.mat')) % keeps Unit_1, skips pupil traces (see buildMergedTableTiers.m)
load(fullfile(cfg.paths.data, 'sensitivity_diff_labeled_N.mat')) % Neural subset

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
% composes them with a second index, e.g. An(cond(An)).
monkeyIdx = getMonkeyIndices(dat);
An = find(monkeyIdx.An);
Ch = find(monkeyIdx.Ch);
Mi = find(monkeyIdx.Mi);

monkeyOrder = {'An', 'Mi', 'Ch'};
monkeyMarkers = struct('An', 'o', 'Mi', 'square', 'Ch', 'diamond');
monkeyNumIdx = struct('An', An, 'Ch', Ch, 'Mi', Mi);

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

%%
% Parse sensitivity (i.e., psychometric slope)

slope_diff = cell2mat(sensitivity_diff_labeled_N(:,2));

LmoresenLSF = slope_diff > 0;
LmoresenHSF = slope_diff < 0;
Lsamesen = slope_diff == 0;

Lother = LmoresenHSF + Lsamesen;

Lothersen = Lother == 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5D-E %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare evidence encoding, discriminability, and behavior, separately for:
%   Figure 5D (figure 99): sessions with greater sensitivity at LSF than HSF (cond = LmoresenLSF)
%   Figure 5E (figure 44): sessions with greater sensitivity at HSF than LSF, or equal (cond = Lothersen)
%
% Panel 2 uses the 200-400ms ROC-area window for both D and E (per this
% file's own header above and the manuscript's Fig 5D/E caption: "(middle)
% MT evidence discriminability (average ROC area, 200-400 ms ...)" with E
% described as "Same as D"). The original script used 50-500ms for D's
% Panel 2, inconsistent with its own header and with E -- fixed here.
% Panels 1/3's viewing-duration XTick still differ between D and E in the
% original script; that difference is preserved, just parameterized per
% panel instead of copy-pasted.

panelCond = {LmoresenLSF, Lothersen};
panelFig = [99 44];
panelROCBin = [200 400; 200 400];            % [bin1 bin2] for Panel 2, per panel
panelTimecourseXTick = [false true];         % whether Panel 1 sets an explicit XTick
panelPerfXTick = {cfg.bins.viewDuration.midpoints, [200,400,600,800,1000]}; % Panel 3 XTick

for ee = 1:2
    cond = panelCond{ee};

    figure(panelFig(ee)); clf

    % Evidence encoding (MT neural activity):
    subplot(1,3,1); hold on; box on;

    testing_dur = 900; % How much of test epoch to display
    pref_only = false;

    end_bin = find(bins(:,1) == testing_dur - 0.5*bin_size,1);

    x_vals_patch = [xax(1:end_bin)' flip(xax(1:end_bin)')]; % X vals -- time relative to test stimulus onset
    x_vals = xax(1:end_bin)';
    x_val_trim = 0; % Testing stim only

    yy = [1,3]; % LSF; HSF

    for cc = 1:2
        dd = yy(cc);
        y_vals = mean(spike_rates_avg_std_norm_avg(cond,1:end_bin,dd), 1, 'omitnan');
        y_vals_SEM = nanse(spike_rates_avg_std_norm_avg(cond,1:end_bin,dd),1);
        y_vals_patch = [y_vals - y_vals_SEM flip(y_vals + y_vals_SEM)];

        p = patch(x_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), y_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), colors{cc});
        set(p,'LineStyle','none')
        alpha(p, 0.2);

        plot(x_vals(~isnan(y_vals) & x_vals > x_val_trim), y_vals(~isnan(y_vals) & x_vals > x_val_trim), 'LineWidth', 0.75, 'Color', colors{cc});
    end

    xlabel('Time relative to test stim onset (ms)')
    ylabel('Avg. firing rate')

    ylim([-0.2 0.8])
    if panelTimecourseXTick(ee)
        set(gca(),'XTick',[200,400,600,800,1000])
    end

    % Evidence discriminability (MT ROC area):
    subplot(1,3,2); hold on; box on;

    bin1 = panelROCBin(ee,1); % Starting bin for ROC area average
    bin2 = panelROCBin(ee,2); % Ending bin for ROC area average

    start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
    end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

    x_line = -0.4:0.1:1;
    y_line = -0.4:0.1:1;

    xlim([0.2 1])
    ylim([0.2 1])

    plot(x_line, y_line, 'k');

    % Monkey Mi:
    x = mean(ROC_area(Mi(cond(Mi)),start_bin:end_bin,3), 2, 'omitnan');
    y = mean(ROC_area(Mi(cond(Mi)),start_bin:end_bin,2), 2, 'omitnan');
    plot(x,y, 'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    % Monkey An:
    x = mean(ROC_area(An(cond(An)),start_bin:end_bin,3), 2, 'omitnan');
    y = mean(ROC_area(An(cond(An)),start_bin:end_bin,2), 2, 'omitnan');
    plot(x,y, 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white', 'MarkerEdgeColor', 'black')

    % Monkey Ch:
    x = mean(ROC_area(Ch(cond(Ch)),start_bin:end_bin,3), 2, 'omitnan');
    y = mean(ROC_area(Ch(cond(Ch)),start_bin:end_bin,2), 2, 'omitnan');
    plot(x,y, 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    x = mean(ROC_area(cond,start_bin:end_bin,3), 2, 'omitnan');
    y = mean(ROC_area(cond,start_bin:end_bin,2), 2, 'omitnan');

    xlabel('ROC area (HSF)')
    ylabel('ROC area (LSF)')

    [p, h, stats] = signrank(x,y);
    d = computeCohenDCI(y, x, 'paired');

    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

    % Performance (avg. binned percent correct)
    subplot(1,3,3); hold on; box on;

    LSF_switch_dat = mean(LSF_switch_binned_performance(cond,:));
    LSF_switch_dat_sem = nanse(LSF_switch_binned_performance(cond,:));

    HSF_switch_dat = mean(HSF_switch_binned_performance(cond,:));
    HSF_switch_dat_sem = nanse(HSF_switch_binned_performance(cond,:));

    errorbar(cfg.bins.viewDuration.midpoints, LSF_switch_dat, LSF_switch_dat_sem, 'Color', colors{1},'LineWidth', 2)
    errorbar(cfg.bins.viewDuration.midpoints, HSF_switch_dat, HSF_switch_dat_sem, 'Color', colors{2}, 'LineWidth', 2)

    ylabel('Fraction correct')
    xlabel('Binned viewing duration (ms)')

    xlim([0 1000])
    ylim([0.4 .9])

    set(gca(),'XTick', panelPerfXTick{ee})
end

%% Relevant Figure Supplements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5-Figure Supplement 2 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test-epoch average activity (HSF vs. LSF), by monkey and by sensitivity
% classification.
% Row 1: sessions more sensitive to evidence at LSF (cond = LmoresenLSF)
% Row 2: sessions more sensitive to evidence at HSF, or equal (cond = Lothersen)
% Panel order (An, Mi, Ch) matches the manuscript, not
% monkeyIdx.names' alphabetical (An, Ch, Mi) order.

bin1 = 200;  % Starting bin for response average
bin2 = 400;  % Ending bin for response average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

sensGroups = {LmoresenLSF, Lothersen};

x_line = -0.2:0.1:2;
y_line = -0.2:0.1:2;

figure(21); clf

for ss = 1:numel(sensGroups)
    cond = sensGroups{ss};

    for mm = 1:numel(monkeyOrder)
        monkeyName = monkeyOrder{mm};
        idx = monkeyNumIdx.(monkeyName);

        subplot(2,3,(ss-1)*3+mm); hold on; box on;

        ylim([-0.2,1])
        xlim([-0.2,1])

        plot(x_line, y_line, 'k');

        x = mean(spike_rates_avg_std_norm_avg(idx(cond(idx)),start_bin:end_bin,3), 2, 'omitnan');
        y = mean(spike_rates_avg_std_norm_avg(idx(cond(idx)),start_bin:end_bin,1), 2, 'omitnan');

        plot(x,y, 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

        xlabel('High switch frequency avg. FR')
        ylabel('Low switch frequency avg. FR')

        [p, h, stats] = signrank(x,y);
        d = computeCohenDCI(y, x, 'paired');

        title(['Testing Epoch average firing rate ', num2str(bin1), '-', num2str(bin2)], ' ms')
        subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5-Figure Supplement 4 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare evidence encoding, discriminability, and behavior for non-switch
% trials. Same row split as Figure 5D-E (row 1 = LmoresenLSF, n=100; row 2 =
% Lothersen, n=32). Columns: (1) test-stimulus evidence encoding (non-switch
% PREF firing rate), (2-3) test-stimulus ROC area in early (50-200 ms) and
% late (200-400 ms) windows, (4) non-switch behavioral performance.
%
% Non-switch PREF/NULL columns of spike_rates_avg_std_norm_avg: 5=LSF
% non-switch NULL, 6=LSF non-switch PREF, 7=HSF non-switch NULL, 8=HSF
% non-switch PREF -- as in the "Output spike matrices organization" block
% above. This is NOT simply "trial started on the adapting-epoch preferred
% direction" (that's what processNeuralData.m's raw dd=1/dd=2 selects):
% because LSF (1 switch) and HSF (5 switches) both have an odd number of
% mid-adapting-epoch reversals, the starting direction only equals the
% test-epoch direction when there's *also* a switch at the adapting/test
% boundary. So dd=1 (starts on the adapting-epoch PREF direction) lands on
% columns 1/3 for switch trials but columns 5/7 for non-switch trials --
% confirmed against this panel's published version (columns 6/8 are the
% ones whose time course actually shows the expected non-switch PREF
% pattern: high baseline, dip at the test-stimulus coherence drop, recovery
% to a sustained plateau). See processNeuralData.m for the full derivation.
% ROC_area's non-switch columns (1=LSF, 4=HSF) don't have this issue -- they
% aren't split by dd, so no start/test-direction mismatch is possible.

% Non-switch behavioral performance, for the Column 4 panels:
dat = mergedTableSub;

[LSF_nonswitch_binned_performance, HSF_nonswitch_binned_performance, ~] = ...
    behaviorBinnedPerformance(dat, 'nonswitch');

rocWindows = [50 200; 200 400]; % [bin1 bin2] for Columns 2 and 3

figure(45); clf
set(gcf, 'Position', [100 100 1500 700]); % 2x4 grid is denser than this file's other panels; widen so titles/subtitles don't collide

for rr = 1:2
    cond = panelCond{rr};

    % Column 1: evidence encoding (non-switch PREF MT activity, LSF vs HSF)
    subplot(2,4,(rr-1)*4+1); hold on; box on;

    testing_dur = 900; % How much of test epoch to display
    end_bin = find(bins(:,1) == testing_dur - 0.5*bin_size,1);

    x_vals_patch = [xax(1:end_bin)' flip(xax(1:end_bin)')];
    x_vals = xax(1:end_bin)';
    x_val_trim = 0; % Testing stim only

    ylim([-0.2 0.8])

    % Shade the two ROC-area windows used in Columns 2-3 (light = early 50-200ms, dark = late 200-400ms)
    patch([50 200 200 50], [-0.2 -0.2 0.8 0.8], [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'LineStyle', 'none');
    patch([200 400 400 200], [-0.2 -0.2 0.8 0.8], [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'LineStyle', 'none');

    yy = [6,8]; % LSF non-switch PREF; HSF non-switch PREF

    for cc = 1:2
        dd = yy(cc);
        y_vals = mean(spike_rates_avg_std_norm_avg(cond,1:end_bin,dd), 1, 'omitnan');
        y_vals_SEM = nanse(spike_rates_avg_std_norm_avg(cond,1:end_bin,dd),1);
        y_vals_patch = [y_vals - y_vals_SEM flip(y_vals + y_vals_SEM)];

        p = patch(x_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), y_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), colors{cc});
        set(p,'LineStyle','none')
        alpha(p, 0.2);

        plot(x_vals(~isnan(y_vals) & x_vals > x_val_trim), y_vals(~isnan(y_vals) & x_vals > x_val_trim), 'LineWidth', 0.75, 'Color', colors{cc});
    end

    xlabel('Time relative to test stim onset (ms)')
    ylabel('Avg. firing rate (non-switch)')

    % Columns 2-3: evidence discriminability (non-switch ROC area), early and late windows
    for ww = 1:2
        subplot(2,4,(rr-1)*4+1+ww); hold on; box on;

        bin1 = rocWindows(ww,1);
        bin2 = rocWindows(ww,2);

        start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
        end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

        x_line = -0.4:0.1:1;
        y_line = -0.4:0.1:1;

        xlim([0.2 1])
        ylim([0.2 1])

        plot(x_line, y_line, 'k');

        % Monkey Mi:
        x = mean(ROC_area(Mi(cond(Mi)),start_bin:end_bin,4), 2, 'omitnan');
        y = mean(ROC_area(Mi(cond(Mi)),start_bin:end_bin,1), 2, 'omitnan');
        plot(x,y, 'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

        % Monkey An:
        x = mean(ROC_area(An(cond(An)),start_bin:end_bin,4), 2, 'omitnan');
        y = mean(ROC_area(An(cond(An)),start_bin:end_bin,1), 2, 'omitnan');
        plot(x,y, 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white', 'MarkerEdgeColor', 'black')

        % Monkey Ch:
        x = mean(ROC_area(Ch(cond(Ch)),start_bin:end_bin,4), 2, 'omitnan');
        y = mean(ROC_area(Ch(cond(Ch)),start_bin:end_bin,1), 2, 'omitnan');
        plot(x,y, 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

        x = mean(ROC_area(cond,start_bin:end_bin,4), 2, 'omitnan');
        y = mean(ROC_area(cond,start_bin:end_bin,1), 2, 'omitnan');

        xlabel('ROC area (HSF, non-switch)')
        ylabel('ROC area (LSF, non-switch)')

        [p, h, stats] = signrank(x,y);
        d = computeCohenDCI(y, x, 'paired');

        title(sprintf('%d-%d ms', bin1, bin2))
        subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])
    end

    % Column 4: non-switch behavioral performance
    subplot(2,4,(rr-1)*4+4); hold on; box on;

    LSF_nonswitch_dat = mean(LSF_nonswitch_binned_performance(cond,:));
    LSF_nonswitch_dat_sem = nanse(LSF_nonswitch_binned_performance(cond,:));

    HSF_nonswitch_dat = mean(HSF_nonswitch_binned_performance(cond,:));
    HSF_nonswitch_dat_sem = nanse(HSF_nonswitch_binned_performance(cond,:));

    errorbar(cfg.bins.viewDuration.midpoints, LSF_nonswitch_dat, LSF_nonswitch_dat_sem, 'Color', colors{1},'LineWidth', 2)
    errorbar(cfg.bins.viewDuration.midpoints, HSF_nonswitch_dat, HSF_nonswitch_dat_sem, 'Color', colors{2}, 'LineWidth', 2)

    ylabel('Fraction correct (non-switch)')
    xlabel('Binned viewing duration (ms)')

    xlim([0 1000])
    ylim([0.4 1])
end
