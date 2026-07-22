%% Fig4_FigSupp1A.m
% Analyzes the emergence and heterogeneity of context-dependent sensory
% adaptation in MT across adaptation and test epochs, comparing "adapting" vs "facilitating"
% neuron subpopulations.
%
% Main Text Figures:
%
%   Figure 4A: Swarm charts comparing population distributions of firing rates
%     for first stimulus vs test stimulus. Two panels:
%     - First stimulus (-2350 to -2000ms): HSF vs LSF (no significant difference)
%     - Test stimulus (50-500ms): HSF vs LSF (significant difference emerges)
%     Yellow dots = HSF, blue dots = LSF. Black/white bars show means.
%
%   Figure 4B: Time course of LSF-HSF difference. Shows mean±SEM for:
%     - Monkey An (blue circles)
%     - Monkey Ch (red diamonds)
%     - Monkey Mi (green squares)
%     - All neurons combined (black triangles)
%     X-axis: First stimulus vs Test stimulus
%     Demonstrates that differential encoding emerges during adaptation period.
%
%   Figure 4C: Example neurons showing opposing dynamics. Two panels:
%     - "Facilitating" neuron (unit 10): Response increases
%     - "Adapting" neuron (unit 77): Response decreases
%     Shows time course with LSF (blue), HSF (yellow)
%
%   Figure 4D: Quantification of adaptation dynamics for HSF condition.
%     Plots Cohen's d comparing first PREF stimulus vs 1st, 2nd, 3rd subsequent
%     PREF stimuli (stimulus numbers 1, 3, 5 in sequence). Two panels:
%     - Facilitating neurons (top): Positive d values, responses increase
%     - Adapting neurons (bottom): Negative d values, responses decrease
%     Error bars show standard error. Demonstrates progressive changes.
%
%   Figures 4E-F: Test-epoch responses (50-500ms) for neuron subpopulations.
%     Scatter plots comparing HSF vs LSF firing rates:
%     - Figure 4E: "Adapting" neurons
%     - Figure 4F: "Facilitating" neurons
%     Different markers by monkey (circles=An, diamonds=Ch, squares=Mi)
%
% Figure Supplements:
%
%   Figure 4-Figure Supplement 1A: Same as Figure 3C (test-epoch HSF vs LSF scatter)
%     but for adapting-stimulus onset epoch (-2350 to -2000ms), separated by
%     monkey. Shows minimal LSF-HSF difference early in trial, demonstrating
%     that context-dependent encoding emerges over time rather than being present
%     from stimulus onset.
%
% Required data files:
%   - mergedTable_proc.mat
%   - sensitivity_diff_labeled_N.mat
%
% Required functions:
%   - processNeuralData.m (extract MT activity and selectivity)
%   - computeCohenDCI.m (calculate effect size)

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

[spike_rates_avg, spike_rates_SEM, spike_rates_avg_std_norm_avg, ...
    spike_rates_avg_std_norm_SEM, ...
    ROC_area, cell_selectivity, ...
    normalizationTerm, uniqueUnitNames, ...
    unit_example] = processNeuralData(dat, slide, bin_size, start_time, end_time, [], []);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 4A %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare population response distribution for first (i.e. adapting stimulus onset) and test stimulus

% First stimulus presentation:
bin1 = -2350; % Starting bin for response average
bin2 = -2000; % Ending bin for response average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

% Cols 3&7 (not 3&8) and 1&5 (not 1&6) are deliberate: this window is
% during the adapting epoch, before test onset, so trials should be
% grouped by adapting-epoch starting direction (dd=1), not by the
% test-epoch PREF/NULL label -- and dd=1 maps to PREF for switch trials
% but NULL for non-switch trials (see processNeuralData.m header note).
HSF_1 = mean([mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,3), 2, 'omitnan'), mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,7), 2, 'omitnan')], 2, 'omitnan');
LSF_1 = mean([mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,1), 2, 'omitnan'), mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,5), 2, 'omitnan')], 2, 'omitnan');
diff_1 = LSF_1 - HSF_1;

% Final stimulus presentation:
bin1 = 50;
bin2 = 500;

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

HSF_3 = mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,3), 2, 'omitnan');
LSF_3 = mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,1), 2, 'omitnan');
diff_3 = LSF_3 - HSF_3;

% Swarm chart for progression of adaptation

% First stimulus:
figure(11); clf;
subplot(1,2,1); hold on; box on;

swarmchart(repmat(1,size(HSF_1)), HSF_1, 60, 'MarkerFaceColor',cfg.colors.HSF, 'MarkerFaceAlpha',0.2,'MarkerEdgeColor',cfg.colors.HSF)
swarmchart(repmat(2,size(LSF_1)), LSF_1, 60, 'MarkerFaceColor',cfg.colors.LSF, 'MarkerFaceAlpha',0.2,'MarkerEdgeColor',cfg.colors.LSF)

x1 = [0.75, 1.25];
y1 = [mean(HSF_1, 'omitnan'), mean(HSF_1, 'omitnan')];

x2 = [1.75, 2.25];
y2 = [mean(LSF_1, 'omitnan'), mean(LSF_1, 'omitnan')];

plot(x1,y1,'k', 'LineWidth',4)
plot(x2,y2,'w', 'LineWidth',4)

[h,p,ci,stats] = ttest(HSF_1(:,1),LSF_1(:,1));

set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'HSF' 'LSF'})
title("First stim.")
subtitle(['p = ', num2str(p)])

% Testing stimulus:
subplot(1,2,2); hold on; box on;

swarmchart(repmat(1,size(HSF_3)), HSF_3, 60, 'MarkerFaceColor',cfg.colors.HSF, 'MarkerFaceAlpha',0.2,'MarkerEdgeColor',cfg.colors.HSF)
swarmchart(repmat(2,size(LSF_3)), LSF_3, 60, 'MarkerFaceColor',cfg.colors.LSF, 'MarkerFaceAlpha',0.2,'MarkerEdgeColor',cfg.colors.LSF)

x1 = [0.75, 1.25];
y1 = [mean(HSF_3, 'omitnan'), mean(HSF_3, 'omitnan')];

x2 = [1.75, 2.25];
y2 = [mean(LSF_3, 'omitnan'), mean(LSF_3, 'omitnan')];

plot(x1,y1,'k', 'LineWidth',4)
plot(x2,y2,'w', 'LineWidth',4)

[h,p,ci,stats] = ttest(HSF_3(:,1),LSF_3(:,1));

set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'HSF' 'LSF'})
title("Testing stim.")
subtitle(['p = ', num2str(p)])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 4B %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare average response for first (i.e. adapting stimulus onset) and test stimulus
% Individual monkey averages plotted separately

figure(57); clf; hold on; box on;

errorbar([1,2], [mean(diff_1(An), 'omitnan'), mean(diff_3(An), 'omitnan')], [nanse(diff_1(An),1), nanse(diff_3(An),1)], 'b-o', 'MarkerSize',12, 'MarkerFaceColor', 'white')
errorbar([1,2], [mean(diff_1(Ch), 'omitnan'), mean(diff_3(Ch), 'omitnan')], [nanse(diff_1(Ch),1), nanse(diff_3(Ch),1)], 'r-diamond', 'MarkerSize',12, 'MarkerFaceColor', 'white')
errorbar([1,2], [mean(diff_1(Mi), 'omitnan'), mean(diff_3(Mi), 'omitnan')], [nanse(diff_1(Mi),1), nanse(diff_3(Mi),1)], 'g-square', 'MarkerSize',12, 'MarkerFaceColor', 'white')

errorbar([1,2], [mean(diff_1, 'omitnan'), mean(diff_3, 'omitnan')], [nanse(diff_1,1), nanse(diff_3,1)], '-^', 'MarkerSize',13, 'MarkerFaceColor', 'black')

set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'First stim.' 'Testing stim.'})

xlim([0.5 2.5])
ylim([-0.025 0.325])

[h,p,ci,stats] = ttest(diff_1,diff_3);
subtitle(['p = ', num2str(p)])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 4C %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example neurons from "facilitating" and "adapting" response groups
% Some plot controls

testing_dur = 900; % How much of test stim to display (ms)
pref_only = false;

start_bin = 1;
end_bin = find(bins(:,1) == testing_dur - 0.5*bin_size,1);

% 8 conditions = LSF/LSF/HSF/HSF/LSF/LSF/HSF/HSF (see processNeuralData.m)
co = cfg.colors.pair([1 1 2 2 1 1 2 2]);

x_vals_patch = [xax(1:end_bin)' flip(xax(1:end_bin)')];
x_vals = xax(1:end_bin)';

% Control conditions plotted from matrix
if pref_only
    cc = [1,3];
else
    cc = [1,2,3,4];
end

figure(24); clf;

% Example "facilitating" neuron:
subplot(2,1,1); hold on; box on;

uu = 10;

for yy = 1:length(cc)
    y_vals = spike_rates_avg_std_norm_avg(uu,1:end_bin,cc(yy));
    y_vals_SEM = spike_rates_avg_std_norm_SEM(uu,1:end_bin,cc(yy));

    x_val_trim = -2500;
    y_vals_patch = [y_vals - y_vals_SEM flip(y_vals + y_vals_SEM)];

    p = patch(x_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), y_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), co{cc(yy)});
    set(p,'LineStyle','none')
    alpha(p, 0.2);

    plot(x_vals(~isnan(y_vals) & x_vals > x_val_trim), y_vals(~isnan(y_vals) & x_vals > x_val_trim), 'LineWidth', 0.75, 'Color', co{cc(yy)});
end

ylim([-0.2 1])

% Example "adapting" neuron:
subplot(2,1,2); hold on; box on;

uu = 77;

for yy = 1:length(cc)
    y_vals = spike_rates_avg_std_norm_avg(uu,1:end_bin,cc(yy));
    y_vals_SEM = spike_rates_avg_std_norm_SEM(uu,1:end_bin,cc(yy));

    x_val_trim = -2500;
    y_vals_patch = [y_vals - y_vals_SEM flip(y_vals + y_vals_SEM)];

    p = patch(x_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), y_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), co{cc(yy)});
    set(p,'LineStyle','none')
    alpha(p, 0.2);

    plot(x_vals(~isnan(y_vals) & x_vals > x_val_trim), y_vals(~isnan(y_vals) & x_vals > x_val_trim), 'LineWidth', 0.75, 'Color', co{cc(yy)});
end

ylim([-0.2 1])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 4D %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Neural response change over repeated HSF pref-motion stimulus presentations
% Averaged and plotted separately for "adapting" and "facilitating" cells

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

%%
% Calculate Cohen's d between first preferred-motion stimulus response and
% subsequent preferred-motion responses at HSF, separately for the
% "adapting" and "facilitating" subsets.

% plot controls:
dirLabel = "PREF";
cond = [2,2,2];  % 2 (HSF); 1 (LSF)
stim1 = 1;       % Relative to first preferred-motion stimulus
stim2 = [1,3,5]; % Comparison stimuli: first, third, fifth (which are first, second, and third pref-motion stimuli)

% Group 1 = "adapting" cells -> d_PREF_a/se_PREF_a, group 2 = "facilitating" cells -> d_PREF_f/se_PREF_f
subsetGroups = {adapt, facilitate};

for gg = 1:numel(subsetGroups)
    subset = subsetGroups{gg};

    % Initialize vector for Cohen's d values
    d = nans(3,1);
    SE = nans(3,1);

    for ii = 1:3

        if stim1 == 1
            bin_f_stim1 = -2200;
            bin_f_stim2 = -2000;

            if cond(ii) == 1 && dirLabel == "PREF"
                idx_1 = 1;
            elseif cond(ii) == 2 && dirLabel == "PREF"
                idx_1 = 3;
            else
                error('Unhandled cond/dirLabel combination for idx_1');
            end
        end

        if stim2(ii) == 1
            bin_s_stim1 = -2200;
            bin_s_stim2 = -2000;

            if cond(ii) == 1 && dirLabel == "PREF"
                idx_2 = 1;
            elseif cond(ii) == 2 && dirLabel == "PREF"
                idx_2 = 3;
            else
                error('Unhandled cond/dirLabel combination for idx_2');
            end

        elseif stim2(ii) == 3
            bin_s_stim1 = -1400;
            bin_s_stim2 = -1200;

            idx_2 = 3;

        elseif stim2(ii) == 5
            bin_s_stim1 = -600;
            bin_s_stim2 = -400;

            idx_2 = 3;
        end

        start_bin_f_stim = find(bins(:,1) == bin_f_stim1 - 0.5*bin_size,1);
        end_bin_f_stim = find(bins(:,1) == bin_f_stim2 - 0.5*bin_size,1);

        start_bin_s_stim = find(bins(:,1) == bin_s_stim1 - 0.5*bin_size,1);
        end_bin_s_stim = find(bins(:,1) == bin_s_stim2 - 0.5*bin_size,1);

        x = mean(spike_rates_avg_std_norm_avg(subset,start_bin_f_stim:end_bin_f_stim,idx_1), 2, 'omitnan');
        y = mean(spike_rates_avg_std_norm_avg(subset,start_bin_s_stim:end_bin_s_stim,idx_2), 2, 'omitnan');

        [p, h, stats] = signrank(x,y);
        [d(ii,1), ~, SE(ii,1)] = computeCohenDCI(y, x, 'paired');

        if isnan(d(ii,1))
            d(ii,1) = 0;
        end
    end

    if gg == 1
        d_PREF_a = d;
        se_PREF_a = SE;
    else
        d_PREF_f = d;
        se_PREF_f = SE;
    end
end

%%
figure(4); clf;

subplot(2,1,1); box on; hold on;

% Facilitating
x = [1,3,5];
y = [d_PREF_f(1), d_PREF_f(2), d_PREF_f(3)];
e = [se_PREF_f(1), se_PREF_f(2), se_PREF_f(3)];

errorbar(x,y,e,'ko','MarkerSize',6,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w','LineStyle', '-');

ylim([0 1.75])
xlim([0.5 5.5])

subplot(2,1,2); box on; hold on;

% Adapting
x = [1,3,5];
y = [d_PREF_a(1), d_PREF_a(2), d_PREF_a(3)];
e = [se_PREF_a(1), se_PREF_a(2), se_PREF_a(3)];

errorbar(x,y,e,'k.','MarkerSize',20,'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k','LineStyle', '-');

ylabel("Cohen's d")
xlabel("PREF motion Stimulus number")

ylim([-1.75 0])
xlim([0.5 5.5])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 4E-F %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test-stimulus average response (HSF vs. LSF)
% Averaged and plotted separately for "adapting" and "facilitating" cells

bin1 = 50;  % Starting bin for test-stimulus response average
bin2 = 500; % Ending bin for test-stimulus response average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

x_line = -0.2:0.1:2;
y_line = -0.2:0.1:2;

% Panel 1 = "adapting" cells, panel 2 = "facilitating" cells
cellGroups = {adapt, facilitate};

figure(33); clf

for gg = 1:numel(cellGroups)
    cond = cellGroups{gg};

    subplot(1,2,gg); hold on; box on;

    plot(x_line, y_line, 'k');

    x = mean(spike_rates_avg_std_norm_avg(Ch(cond(Ch)),start_bin:end_bin,3), 2, 'omitnan');
    y = mean(spike_rates_avg_std_norm_avg(Ch(cond(Ch)),start_bin:end_bin,1), 2, 'omitnan');
    plot(x,y, 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    x = mean(spike_rates_avg_std_norm_avg(Mi(cond(Mi)),start_bin:end_bin,3), 2, 'omitnan');
    y = mean(spike_rates_avg_std_norm_avg(Mi(cond(Mi)),start_bin:end_bin,1), 2, 'omitnan');
    plot(x,y,'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    x = mean(spike_rates_avg_std_norm_avg(An(cond(An)),start_bin:end_bin,3), 2, 'omitnan');
    y = mean(spike_rates_avg_std_norm_avg(An(cond(An)),start_bin:end_bin,1), 2, 'omitnan');
    plot(x,y, 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    xlabel('High switch frequency avg. FR')
    ylabel('Low switch frequency avg. FR')

    x = mean(spike_rates_avg_std_norm_avg(cond,start_bin:end_bin,3), 2, 'omitnan');
    y = mean(spike_rates_avg_std_norm_avg(cond,start_bin:end_bin,1), 2, 'omitnan');

    [p, h, stats] = signrank(x,y);
    d = computeCohenDCI(y, x, 'paired');

    title(['Testing Epoch average firing rate ', num2str(bin1), '-', num2str(bin2)], ' ms')
    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

    ylim([-0.2,1])
    xlim([-0.2,1])
end

%% Relevant Figure Supplements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 4-Figure Supplement 1A %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Average response to adapting-stimulus onset (HSF vs. LSF), by animal
% Panel order (An, Mi, Ch) matches the manuscript, not
% monkeyIdx.names' alphabetical (An, Ch, Mi) order.

bin1 = -2350; % Starting bin for response average (50 ms delay)
bin2 = -2000; % Ending bin for response average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

x_line = -0.2:0.1:2;
y_line = -0.2:0.1:2;

figure(40); clf;

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    isMonkey = monkeyIdx.(monkeyName);

    subplot(1,3,mm); hold on; box on;

    plot(x_line, y_line, 'k');

    x = mean([mean(spike_rates_avg_std_norm_avg(isMonkey,start_bin:end_bin,3), 2, 'omitnan'), mean(spike_rates_avg_std_norm_avg(isMonkey,start_bin:end_bin,7), 2, 'omitnan')], 2, 'omitnan');
    y = mean([mean(spike_rates_avg_std_norm_avg(isMonkey,start_bin:end_bin,1), 2, 'omitnan'), mean(spike_rates_avg_std_norm_avg(isMonkey,start_bin:end_bin,5), 2, 'omitnan')], 2, 'omitnan');

    plot(x,y, 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    [p, h, stats] = signrank(x,y);
    d = computeCohenDCI(y, x, 'paired');

    ylim([-0.2,1])
    xlim([-0.2,1])

    title(['Testing Epoch average firing rate ', num2str(bin1), '-', num2str(bin2)], ' ms')
    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])
end
