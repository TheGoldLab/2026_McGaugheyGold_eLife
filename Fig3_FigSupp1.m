%% Fig3_FigSupp1.m
% Analyzes MT neural responses to preferred-direction motion during test epoch,
% comparing activity between LSF and HSF switch-frequency conditions.
%
% Analysis:
%   - Processes correct trials using processNeuralData.m
%   - Extracts activity during test epoch (50-500ms after test onset)
%   - Compares PREF direction responses between LSF and HSF
%   - Analyzes time course with exponential fits to characterize adaptation dynamics
%
% Main Text Figures:
%
%   Figure 3A: Example single neuron (unit 132). Shows:
%     - Raster plot: Individual spike times for LSF (blue) and HSF (yellow) trials
%       across adaptation and test epochs
%     - PSTH: Peri-stimulus time histogram showing average firing rate (normalized)
%       with SEM shading for LSF and HSF conditions, plus NULL direction responses
%
%   Figure 3B: Population average response across all MT neurons (N=155).
%     Shows mean±SEM normalized firing rate for LSF PREF (blue), HSF PREF (yellow),
%     and NULL directions. Demonstrates population-level differential encoding.
%
%   Figure 3C: Scatter plot comparing average test-epoch (50-500ms) firing rate
%     for HSF vs LSF conditions. Each point is one neuron. Points above unity
%     line indicate higher responses at LSF. Different markers by monkey:
%       - Monkey An: circles (N=55)
%       - Monkey Ch: diamonds (N=13)
%       - Monkey Mi: squares (N=87)
%
%   Figure 3D: Relationship between neuron selectivity (ROC area from adaptation
%     epoch) and magnitude of LSF-HSF firing rate difference during test epoch.
%     Tests whether more selective neurons show larger context-dependent changes.
%     Includes linear regression fit with correlation statistics.
%
% Figure Supplements:
%
%   Figure 3-Figure Supplement 1A: Same as Figure 3C but separated by individual
%     monkeys in three panels. Shows that context-dependent neural responses
%     are present in each subject.
%
%   Figure 3-Figure Supplement 1B: Time course analysis using exponential fits.
%     Fits single exponential (y = a*exp(b*t)) to test-epoch responses (200-600ms).
%     Extracts time constant tau = 1/|b| for LSF and HSF. Compares
%     tau values across conditions (scatter plots by monkey). Tests whether
%     temporal dynamics of neural responses differ between contexts.
%     Tau values capped at 10s
%
% Required data files:
%   - mergedTable_proc.mat (behavioral and neural data)
%   - sensitivity_diff_labeled_N.mat (neural subset classification)
%
% Required functions:
%   - processNeuralData.m (extract MT activity and selectivity)
%   - computeCohenDCI.m (calculate effect size)
%   - fit() (Curve Fitting Toolbox; Figure 3-Figure Supplement 1B only)

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

example_unit_idx = 132; % Unit number to save for example (Figure 3A)

[spike_rates_avg, spike_rates_SEM, spike_rates_avg_std_norm_avg, ...
    spike_rates_avg_std_norm_SEM, ...
    ROC_area, cell_selectivity, ...
    normalizationTerm, uniqueUnitNames, ...
    unit_example] = processNeuralData(dat, slide, bin_size, start_time, end_time, example_unit_idx, []);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 3A %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% RASTER PLOT %%%

figure(8); clf; hold on; box on;

% unit_example(7)/(3): both start-dir groups for LSF switch trials
% (same switch_frequency & switch_trial, differing only in start_dir --
% see the ndgrid(sfs,[0 1],dirs) struct order in processNeuralData.m).
mat_concat_LSF = [unit_example(7).spike_times; unit_example(3).spike_times];

for rr = 1:size(mat_concat_LSF,1)
    plot(mat_concat_LSF(rr,:),rr,'.', 'Color', colors{1}, 'LineStyle', 'none')
end

% unit_example(8)/(4): both start-dir groups for HSF switch trials
mat_concat_HSF = [unit_example(8).spike_times;unit_example(4).spike_times];

for rr = 1:size(mat_concat_HSF,1)
    plot(mat_concat_HSF(rr,:),size(mat_concat_LSF,1) + rr,'.', 'Color',colors{2}, 'LineStyle', 'none')
end

xlim([-2500 800])
ylim([0 size(mat_concat_LSF,1) + rr])

set(gca, 'YTick', [],'XTickLabel',[])

%%
%%% PSTH %%%

% Some plot controls
testing_dur = 900; % How much of test stim to display (ms)
pref_only = false;

end_bin = find(bins(:,1) == testing_dur - 0.5*bin_size,1);

figure(99); clf; hold on; box on;

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

for yy = 1:length(cc)
    y_vals = spike_rates_avg_std_norm_avg(example_unit_idx,1:end_bin,cc(yy));
    y_vals_SEM = spike_rates_avg_std_norm_SEM(example_unit_idx,1:end_bin,cc(yy));

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
%%%%%%%%%%%%%%%%%%%% Figure 3B %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Population average response

% Some plot controls
testing_dur = 900; % How much of test stim to display (ms)
pref_only = false;

end_bin = find(bins(:,1) == testing_dur - 0.5*bin_size,1);

% Set up colors
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

figure(99); clf; hold on; box on;

for yy = 1:length(cc)
    y_vals = nanmean(spike_rates_avg_std_norm_avg(:,1:end_bin,cc(yy)));
    y_vals_SEM = nanse(spike_rates_avg_std_norm_avg(:,1:end_bin,cc(yy)),1);

    x_val_trim = -2500;
    y_vals_patch = [y_vals - y_vals_SEM flip(y_vals + y_vals_SEM)];

    p = patch(x_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), y_vals_patch(~isnan(y_vals_patch) & x_vals_patch > x_val_trim), co{cc(yy)});
    set(p,'LineStyle','none')
    alpha(p, 0.2);

    plot(x_vals(~isnan(y_vals) & x_vals > x_val_trim), y_vals(~isnan(y_vals) & x_vals > x_val_trim), 'LineWidth', 0.75, 'Color', co{cc(yy)});
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 3C %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot controls
bin1 = 50;  % Starting bin for response average
bin2 = 500; % Ending bin for response average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

x_line = -0.2:0.1:2;
y_line = -0.2:0.1:2;

figure(2); clf; hold on; box on;

plot(x_line, y_line, 'k');

x = nanmean(spike_rates_avg_std_norm_avg(Ch,start_bin:end_bin,3),2);
y = nanmean(spike_rates_avg_std_norm_avg(Ch,start_bin:end_bin,1),2);
plot(x,y, 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

x = nanmean(spike_rates_avg_std_norm_avg(Mi,start_bin:end_bin,3),2);
y = nanmean(spike_rates_avg_std_norm_avg(Mi,start_bin:end_bin,1),2);
plot(x,y,'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

x = nanmean(spike_rates_avg_std_norm_avg(An,start_bin:end_bin,3),2);
y = nanmean(spike_rates_avg_std_norm_avg(An,start_bin:end_bin,1),2);
plot(x,y, 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

ylim([-0.2,1])
xlim([-0.2,1])

xlabel('High switch frequency avg. FR')
ylabel('Low switch frequency avg. FR')

x = nanmean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,3),2);
y = nanmean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,1),2);

[p, h, stats] = signrank(x,y);
d = computeCohenDCI(y, x, 'paired');

title(['Testing Epoch average firing rate ', num2str(bin1), '-', num2str(bin2)], ' ms')
subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 3D %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Correlate cell selectivity with avg. difference between LSF and HSF

bin1 = 50;  % Starting bin for response average
bin2 = 500; % Ending bin for response average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

% Vector of difference in response between LSF and HSF to testing stimulus
HSF_resp = nanmean(spike_rates_avg(:,start_bin:end_bin,3),2);
LSF_resp = nanmean(spike_rates_avg(:,start_bin:end_bin,1),2);
Adapt_diff = abs(LSF_resp - HSF_resp);

figure(7); clf; hold on; box on;

x = cell_selectivity(An,:);
y = Adapt_diff(An,:);
plot(x,y, 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

x = cell_selectivity(Ch,:);
y = Adapt_diff(Ch,:);
plot(x,y, 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

x = cell_selectivity(Mi,:);
y = Adapt_diff(Mi,:);
plot(x,y,'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

mdl = fitlm(cell_selectivity,Adapt_diff);
p = plot(mdl);
p(1).Visible = 'off';
legend('hide')
title(' ')

xlabel('Cell selectivity')
ylabel('Avg. firing rate difference (LSF-HSF)')

[r,p] = corr(cell_selectivity, Adapt_diff);
subtitle(['p = ', num2str(p), ' & r = ', num2str(r)])

%% Relevant Figure Supplements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 3-Figure Supplement 1A %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot average test-stimulus activity (HSF vs. LSF) by animal
% Panel order (An, Mi, Ch) matches the manuscript, not
% monkeyIdx.names' alphabetical (An, Ch, Mi) order.

bin1 = 50;  % Starting bin for response average
bin2 = 500; % Ending bin for response average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

x_line = -0.2:0.1:2;
y_line = -0.2:0.1:2;

figure(21); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    isMonkey = monkeyIdx.(monkeyName);

    subplot(1,3,mm); hold on; box on;

    plot(x_line, y_line, 'k');

    x = nanmean(spike_rates_avg_std_norm_avg(isMonkey,start_bin:end_bin,3),2);
    y = nanmean(spike_rates_avg_std_norm_avg(isMonkey,start_bin:end_bin,1),2);
    plot(x,y, 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    [p, h, stats] = signrank(x,y);
    d = computeCohenDCI(y, x, 'paired');

    ylim([-0.2,1])
    xlim([-0.2,1])

    xlabel('High switch frequency avg. FR')
    ylabel('Low switch frequency avg. FR')

    title(['Testing Epoch average firing rate ', num2str(bin1), '-', num2str(bin2)], ' ms')
    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 3-Figure Supplement 1B %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assess timecourse of MT neural response with single exponential fits
% Compare Tau values between HSF and LSF

bin1 = 200;
bin2 = 600;

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);
x_vals = xax(start_bin:end_bin);

cc = [1,3];
plot_it = false;

exp_fits = nan(numUnits,2,2);
exp_fits_tau = nan(numUnits,1,2);

% Loop through units and fit
for uu = 1:numUnits

    if plot_it
        figure(10); clf
    end

    for ff = 1:2

        sf_testing = spike_rates_avg_std_norm_avg(uu,start_bin:end_bin,cc(ff))';
        sf_testing = sf_testing(~isnan(sf_testing),:);

        % Binning:
        b_1 = nanmean(sf_testing(1:14));  % 200-330 ms
        b_2 = nanmean(sf_testing(15:28)); % 340-470 ms
        b_3 = nanmean(sf_testing(29:end));% 480-600 ms

        xData = [1,140,275]'; % x values at bin midpoints
        yData = [b_1,b_2,b_3]';

        f = fit(xData, yData, 'exp1');

        % Save coefficients:
        coefficients = coeffvalues(f);
        exp_fits(uu,1,ff) = coefficients(1);
        exp_fits(uu,2,ff) = coefficients(2);

        % Cap tau values at 10 s
        exp_fits_tau(uu,1,ff) = min(abs(1/coefficients(2)), 10000);

        % Visualize fit:
        if plot_it
            plot(xData,yData, '.-', 'Color', colors{ff}, 'MarkerSize', 20)
            hold on
            plot(f, xData, yData)
            legend('Data', 'Fitted curve')
        end
    end
end

%%
% Use fitted values for plots and stats:
% x = HSF tau (condition cc(2)=3), y = LSF tau (condition cc(1)=1)

x = exp_fits_tau(:,1,2);
y = exp_fits_tau(:,1,1);

x_line = 0:10:10000;
y_line = 0:10:10000;

figure(12); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    isMonkey = monkeyIdx.(monkeyName);

    subplot(1,3,mm); hold on; box on;

    plot(x_line, y_line, 'k')
    plot(x(isMonkey), y(isMonkey),'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',12, 'LineStyle', 'none', 'MarkerFaceColor','white')

    d = computeCohenDCI(y(isMonkey), x(isMonkey), 'paired');
    [p, h, stats] = signrank(x(isMonkey), y(isMonkey));

    title('T values from exponential fit')
    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

    if mm == 1
        xlabel('High switch frequency T')
        ylabel('Low switch frequency T')
    end

    ylim([0 10000])
    xlim([0 10000])
end
