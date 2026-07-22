%% Fig6_FigSupp2AB.m
% NOTE ON FIGURE NUMBERING: this script produces all of Figure 6 (main
% text) as well as Figure 6-Figure Supplement 2, panels A and B (panels C
% and D of that same supplement are produced separately, by
% Fig6_FigSupp2CD.m). Figure 6-Figure Supplement 2 was formerly numbered
% Figure 7-Figure Supplement 1 in the manuscript; it was renumbered because
% these panels are Figure 6B/6D replotted per monkey, so they belong under
% Figure 6's own supplement rather than Figure 7's.
% Analyzes pupil dynamics across switch-frequency conditions and their
% relationship to behavioral sensitivity differences.
%
% Pupil processing:
%   Uses processPupilFull.m to extract evoked pupil responses with:
%   - Low-pass filtering (5 Hz Butterworth)
%   - Outlier removal and interpolation
%   - Z-scoring after subtracting running average
%   - Regression analysis modeling pupil as function of switch-frequency condition and baseline 
%
% Context-stability β coefficient:
%   From regression: Pupil = β₀ + β₁*baseline + β₂*switch-frequency condition + β₃*time
%   β₂ (parameter 3) quantifies differential pupil response to LSF vs HSF
%
% Main Text Figures:
%
%   Figure 6A: Across-session average evoked pupil time courses. Two panels:
%     - Monkey An (left, N=49 sessions)
%     - Monkey Mi (right, N=75 sessions)
%     Shows mean±SEM pupil diameter 
%     LSF (blue) and HSF (yellow) conditions
%
%   Figure 6B: Context-stability β coefficient averaged across monkeys.
%     Shows mean±SEM of β₂ across time.
%     X-axis: Time relative to test-stimulus onset
%     Y-axis: β coefficient magnitude
%     Significance markers (black dots) show bins where β₂ significantly
%     differs from zero (t-test, p<0.05)
%     Peak differences occur 500ms before test onset
%
%   Figure 6C: Time-resolved correlations between pupil β coefficient and
%     behavioral sensitivity difference (LSF-HSF psychometric slope).
%     Three lines:
%     - Monkey Mi (red)
%     - Monkey An (gray)
%     - Combined (black)
%     Significance markers for each (colored dots)
%     Strongest correlations in pre-test period (-500 to 0ms)
%
%   Figure 6D: Summary correlation for pre-test epoch (-500 to 0ms).
%     Scatter plot: X = behavioral slope difference (LSF-HSF sensitivity)
%                   Y = average pupil β coefficient (-500 to 0ms)
%     Each point is one session. Circles = An, Squares = Mi
%     Linear regression fit with correlation statistics
%     Negative correlation: Sessions with larger pupil differences show
%     smaller (more negative) behavioral sensitivity differences
%
% Figure Supplements:
%
%   Figure 6-Figure Supplement 2A: Same as Figure 6B but separated by monkey.
%     Two panels, An (left) and Mi (right).
%
%   Figure 6-Figure Supplement 2B: Same as Figure 6D but separated by monkey.
%     Two panels showing correlation between behavioral sensitivity and
%     pre-test pupil β for An (left) and Mi (right) separately.
%
% Required data files:
%   - mergedTable_proc.mat
%   - sensitivity_diff_labeled_BP.mat (behavioral+pupil subset with slopes)
%
% Required functions:
%   - processPupilFull.m (extract pupil measures and regression coefficients)

%% Load data
cfg = projectDefaults();

load(fullfile(cfg.paths.data, 'mergedTable_proc.mat'))
load(fullfile(cfg.paths.data, 'sensitivity_diff_labeled_BP.mat'))

% Subset data appropriately
% "BP" = Sessions with behavior and pupil data
    % Two sessions from Monkey An had no pupil data recorded
[mergedTableSub] = createDatSubset(mergedTable_proc, 'BP');
dat = mergedTableSub;
clear mergedTable_proc

% Select correct trials
datCorrect = dat(dat.correct == 1,:);
dat = datCorrect;

% Create monkey indices for subsetting sessions.
% Numeric (not logical) indices are needed here because later code
% concatenates them, e.g. fits([An;Mi],...). An and Mi are column vectors
% with different numbers of sessions, so they must be vertically (not
% horizontally) concatenated.
monkeyIdx = getMonkeyIndices(dat);
An = find(monkeyIdx.An);
Mi = find(monkeyIdx.Mi);

% Get unique session names
uniqueSessionNames = unique(dat.ses_ID);

%%
% Process pupil data calling processPupilFull.m
% Binning scheme

slide = 10;
bin_size = 100;
start_time = -200;
end_time = 3600;
bins = cat(2, ...
    (start_time:slide:end_time-bin_size)', ...
    (start_time+bin_size:slide:end_time)');
tax = start_time:end_time;
xax = mean(bins,2);
numBins = size(bins, 1);

save_trial_data = false; % Suppress output from individual sessions saved to path
outPath = ''; % Path for saving data, if needed

[fits, pupil_bin_mean_save_LSF, pupil_bin_mean_save_HSF, ...
    pupil_bin_slope_save_LSF, pupil_bin_slope_save_HSF, ...
    pupil_baseline_save_LSF, pupil_baseline_save_HSF, ...
    ses_pupil_bin_mean_save, ses_pupil_bin_mean_ses_diff] = ...
    processPupilFull(dat, slide, bin_size, start_time, end_time, save_trial_data, outPath);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6A %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot average evoked pupil response for both LSF and HSF conditions
% Plot separately for Monkey An and Monkey Mi

xax_plot = xax - 2400; % Align to testing stimulus onset, not adapting

monkeyOrder = {'An', 'Mi'};
monkeyNumIdx = struct('An', An, 'Mi', Mi);
monkeyMarkers = struct('An', 'o', 'Mi', 'square');
monkeyYlim = struct('An', [-1.05 0], 'Mi', [-2.05 0]); % Different pupil scale per monkey

figure(22); clf;

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    idx = monkeyNumIdx.(monkeyName);

    subplot(1,2,mm); hold on; box on;

    y_vals_HSF = pupil_bin_mean_save_HSF(idx,:);
    y_vals_LSF = pupil_bin_mean_save_LSF(idx,:);

    y_means_HSF = mean(y_vals_HSF,1,'omitnan');
    y_sems_HSF  = nanse(y_vals_HSF,1);

    y_means_LSF = mean(y_vals_LSF,1,'omitnan');
    y_sems_LSF  = nanse(y_vals_LSF,1);

    p=patch([xax_plot; flip(xax_plot)], [y_means_HSF-y_sems_HSF flip(y_means_HSF+y_sems_HSF)]', cfg.colors.HSF );
    set(p,'LineStyle','none')
    alpha(p, 0.2);
    plot(xax_plot, y_means_HSF, 'LineWidth', 2, 'Color', cfg.colors.HSF )

    p=patch([xax_plot; flip(xax_plot)], [y_means_LSF-y_sems_LSF flip(y_means_LSF+y_sems_LSF)]', cfg.colors.LSF );
    set(p,'LineStyle','none')
    alpha(p, 0.2);
    plot(xax_plot, y_means_LSF, 'LineWidth', 2, 'Color', cfg.colors.LSF )

    xlabel('Time since dots on (ms)')
    ylabel('Pupil mean')

    xlim([-2500 700])
    ylim(monkeyYlim.(monkeyName))
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6B %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot context-stability beta coefficient averaged across Monkey An and Monkey Mi
% Beta coefficient averaged within 500 ms window before test-stimulus onset

% Plot all betas

monk = [An;Mi]; % Column vectors of different lengths -- must be vertcat, not horzcat
ymax = .1;

% fits is [numSessions x numBins x 3 params x 2 measures] (processPupilFull.m):
% dim3=3 selects beta2 (hazard-condition coefficient, per header above);
% dim4=1 selects the mean-pupil regression (dim4=2 would be the slope).
for bb = 1:numBins
    [h,p] = ttest(fits(monk,bb,3,1));
    p_val(bb,1) = p;
end

Lpval = p_val < 0.05;

sig_line = repmat(ymax,numBins,1);
sig_line(~Lpval) = nan;

figure(10); cla reset; hold on; box on;

plot([-2400 1000], [0 0], 'k:')
plot([-2400 -2400], ymax.*[-1 1], 'k-');
plot([0 0], ymax.*[-1 1], 'k-');
axis([-2500 700 -ymax ymax]);
ylim([-0.05 0.1])

y_means = mean(fits(monk,:,3,1), 1, 'omitnan');
y_sems  = nanse(fits(monk,:,3,1), 1);

p=patch([xax_plot; flip(xax_plot)], [y_means-y_sems flip(y_means+y_sems)]', cfg.colors.LSF);
set(p,'LineStyle','none')
alpha(p, 0.2);
plot(xax_plot, y_means, 'LineWidth', 2, 'Color', cfg.colors.LSF)

xlabel('Time since dots on (ms)')
ylabel('Context-stability B coef.')

% Add significance

scatter(xax_plot, sig_line, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6C %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot correlation between context-stability beta coefficient as a function of viewing time and behavioral sensitivity (LSF-HSF)
% Correlations run separately for Monkey An and Monkey Mi
% Average correlation across monkeys also included

slope_diff = cell2mat(sensitivity_diff_labeled_BP(:,2)); % Difference in slope terms/sensitivity

% Loop through bins and calculate correlation

HRpupil_behav_corrcoef_Mi = nans(1,numBins);
HRpupil_behav_corrcoef_Mi_p = nans(1,numBins);

HRpupil_behav_corrcoef_An = nans(1,numBins);
HRpupil_behav_corrcoef_An_p = nans(1,numBins);

HRpupil_behav_corrcoef = nans(1,numBins);
HRpupil_behav_corrcoef_p = nans(1,numBins);

for bb = 1:numBins

    [r_Mi,p_Mi] = corr(fits(Mi,bb,3,1), slope_diff(Mi), 'Type', 'Spearman', 'rows', 'complete');

    [r_An,p_An] = corr(fits(An,bb,3,1), slope_diff(An), 'Type', 'Spearman', 'rows', 'complete');
    
    [r,p] = corr(fits([An;Mi],bb,3,1), slope_diff([An;Mi]), 'Type', 'Spearman', 'rows', 'complete');

    HRpupil_behav_corrcoef_Mi(1,bb) = r_Mi;
    HRpupil_behav_corrcoef_Mi_p(1,bb) = p_Mi;

    HRpupil_behav_corrcoef_An(1,bb) = r_An;
    HRpupil_behav_corrcoef_An_p(1,bb) = p_An;

    HRpupil_behav_corrcoef(1,bb) = r;
    HRpupil_behav_corrcoef_p(1,bb) = p;
end

% Prepare significance line:

Lpval_Mi = HRpupil_behav_corrcoef_Mi_p < 0.05;
sig_line_Mi = repmat(0.36,numBins,1);
sig_line_Mi(~Lpval_Mi) = nan;

Lpval_An = HRpupil_behav_corrcoef_An_p < 0.05;
sig_line_An = repmat(0.32,numBins,1);
sig_line_An(~Lpval_An) = nan;

Lpval = HRpupil_behav_corrcoef_p < 0.05;
sig_line = repmat(0.4,numBins,1);
sig_line(~Lpval) = nan;

figure(101); clf; box on; hold on;

plot(xax_plot,HRpupil_behav_corrcoef_Mi, 'LineWidth', 2, 'Color', 'r')
plot(xax_plot,HRpupil_behav_corrcoef_An, 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7])
plot(xax_plot,HRpupil_behav_corrcoef, 'LineWidth', 2, 'Color', 'k')

% Add significance

scatter(xax_plot(1:size(sig_line_Mi,1)), sig_line_Mi, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
scatter(xax_plot(1:size(sig_line_An,1)), sig_line_An, 'MarkerEdgeColor', [0.7, 0.7, 0.7], 'MarkerFaceColor', [0.7, 0.7, 0.7])
scatter(xax_plot(1:size(sig_line,1)), sig_line, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')

ylabel('Correlation coefficient')
xlabel('Time (ms)')

xline(0)
xlim([-2350 700])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6D %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the relationship between behavioral sensitivity (LSF-HSF) and avg. context-stability beta coefficient 
    % Beta coefficient averaged within 500 ms window before test-stimulus onset

% Calculate average of context-stability B coef. (500 ms prior to testing stimuls onset)

bin1 = 1900; % Averaging begins 500 ms prior to test-stimulus onset
bin2 = 2400; % Test-stimulus onset time

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

bin_range = start_bin:end_bin;

figure(11); clf; hold on; box on;

plot(slope_diff(An), mean(fits(An,bin_range,3,1),2),'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white', 'MarkerEdgeColor', 'black')
plot(slope_diff(Mi), mean(fits(Mi,bin_range,3,1),2), 'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

model = fitlm(slope_diff([An;Mi]), mean(fits([An;Mi],bin_range,3,1),2));
p = plot(model);
p(1).Visible='off';
legend('hide')
title('')

[x_2,p_2] = corr(slope_diff([An;Mi]), mean(fits([An;Mi],bin_range,3,1),2), 'Type', 'Spearman');

xlabel('Slope difference (LSF-HSF)')
ylabel('Pupil B term (LSF - HSF)')
subtitle(['p = ', num2str(p_2), ' & r = ', num2str(x_2)])

xlim([-0.05 0.05])

%% Relevant Figure Supplements (Figure 6-Figure Supplement 2, panels A-B; see note at top of file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6-Figure Supplement 2A %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot context-stability beta coefficient
% Separately for Monkey An and Monkey Mi
%
% NOTE: previously hardcoded to `monk = Mi` with a comment saying to
% manually toggle An/Mi and rerun -- so only Mi was ever actually plotted.
% Now loops over both. Also moved off figure(10), which Figure 6B (above)
% also uses -- with the same figure number, this section's `cla reset`
% would silently replace Figure 6B's plot before the script finished
% running, exactly like Fig4_FigSupp1A's original figure(4) clash with its
% own Figure 4D.

ymax = .1;

figure(20); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    monk = monkeyNumIdx.(monkeyName);

    for bb = 1:numBins
        [h,p] = ttest(fits(monk,bb,3,1));
        p_val(bb,1) = p;
    end

    Lpval = p_val < 0.05;

    sig_line = repmat(ymax,numBins,1);
    sig_line(~Lpval) = nan;

    subplot(1,2,mm); hold on; box on;

    plot([-2400 1000], [0 0], 'k:')
    plot([-2400 -2400], ymax.*[-1 1], 'k-');
    plot([0 0], ymax.*[-1 1], 'k-');
    axis([-2500 700 -ymax ymax]);
    ylim([-0.05 0.1])

    y_means = mean(fits(monk,:,3,1), 1, 'omitnan');
    y_sems  = nanse(fits(monk,:,3,1), 1);

    p=patch([xax_plot; flip(xax_plot)], [y_means-y_sems flip(y_means+y_sems)]', cfg.colors.LSF);
    set(p,'LineStyle','none')
    alpha(p, 0.2);
    plot(xax_plot, y_means, 'LineWidth', 2, 'Color', cfg.colors.LSF)

    xlabel('Time since dots on (ms)')
    ylabel('Context-stability B coef.')
    title(sprintf('Monkey %s', monkeyName))

    % Add significance
    scatter(xax_plot, sig_line, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6-Figure Supplement 2B %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the relationship between behavioral sensitivity (LSF-HSF) and avg. context-stability beta coefficient 
    % Beta coefficient averaged within 500 ms window before test-stimulus onset
% Plotting done separately for Monkey An and Monkey Mi

% Calculate average of context-stability B coef. (500 ms prior to testing stimuls onset)

bin1 = 1900; % Averaging begins 500 ms prior to test-stimulus onset
bin2 = 2400; % Test-stimulus onset time

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

bin_range = start_bin:end_bin;

% NOTE: previously used figure(11), the same figure Figure 6D (above) uses
% -- with the same number, this section's `clf` would wipe out Figure 6D's
% plot before the script finished running. Moved to its own figure.
figure(21); clf;

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    idx = monkeyNumIdx.(monkeyName);

    subplot(1,2,mm); hold on; box on;

    plot(slope_diff(idx), mean(fits(idx,bin_range,3,1),2), 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white', 'MarkerEdgeColor', 'black')
    model = fitlm(slope_diff(idx), mean(fits(idx,bin_range,3,1),2));
    p = plot(model);
    p(1).Visible='off';
    legend('hide')
    title('')

    [x_2,p_2] = corr(slope_diff(idx), mean(fits(idx,bin_range,3,1),2), 'Type', 'Spearman');

    xlabel('Slope difference (LSF-HSF)')
    ylabel('Pupil B term (LSF - HSF)')
    subtitle(['p = ', num2str(p_2), ' & r = ', num2str(x_2)])

    xlim([-0.05 0.05])
end