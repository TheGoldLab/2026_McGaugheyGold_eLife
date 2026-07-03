%% Fig5_FigSupp3.m
% Analyzes MT neural activity during incorrect trials, separated by
% switch-frequency condition and behavioral sensitivity.
%
% Analysis approach:
%   - Processes incorrect trials only using processNeuralDataErrorTrials.m
%   - Pools trials across sessions (unlike correct trials which are session-averaged)
%     because individual sessions have too few incorrect trials for reliable estimates
%   - Uses normalization terms from full session analysis for consistent scaling
%   - Separates sessions by behavioral sensitivity: LmoresenLSF (LSF > HSF) vs.
%     Lothersen (HSF >= LSF)
%   - Extracts PREF direction trials only for LSF and HSF conditions
%
% Neural processing:
%   - Bins spike counts in 100ms windows with 10ms slide
%   - Baseline-subtracts and normalizes firing rates
%   - Determines PREF/NULL direction from adaptation epoch
%   - Pools all PREF trials across sessions for each switch-frequency condition
%
% Outputs:
%   Figure 5-Figure Supplement 3: Time course of normalized MT firing rate during
%   test epoch (0-1000ms) for incorrect trials, LSF and HSF switch PREF trials
%   overlaid (blue/yellow with SEM shading). Two panels, one per sensitivity group:
%   - Left: LSF > HSF sessions (LmoresenLSF)
%   - Right: HSF >= LSF sessions (Lothersen)
%
% Required data files:
%   - mergedTable_proc.mat
%   - sensitivity_diff_labeled_N.mat (behavioral sensitivity classification)
%   - normalizationTerm.mat (critical for consistent normalization)
%
% Required functions:
%   - processNeuralDataErrorTrials.m (pooled incorrect trial processing)

%% Load data
cfg = projectDefaults();
cd(cfg.paths.data)

load('mergedTable_proc_neural.mat') % keeps Unit_1, skips pupil traces (see buildMergedTableTiers.m)
load('sensitivity_diff_labeled_N.mat') % Neural subset
load('normalizationTerm.mat') % Normalization term from initial neural processing

% Subset data appropriately
% "N" = Neural analysis
[mergedTableSub] = createDatSubset(mergedTable_proc, 'N');
dat = mergedTableSub;

clear mergedTable_proc

% Select incorrect trials
datIncorrect = dat(dat.correct == 0,:);
dat = datIncorrect;

colors = cfg.colors.pair;

%%
% Process neural data (calling processNeuralDataErrorTrials.m)

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

% Parse sensitivity (i.e., psychometric slope)

slope_diff = cell2mat(sensitivity_diff_labeled_N(:,2));

LmoresenLSF = slope_diff > 0;
LmoresenHSF = slope_diff < 0;
Lsamesen = slope_diff == 0;

Lother = LmoresenHSF + Lsamesen;

Lothersen = Lother == 1;

% Select condition (computed separately for each sensitivity group, one
% panel per group):
    % Sessions more sensitive at LSF: LmoresenLSF (LSF > HSF)
    % Sessions more sensitive at HSF: Lothersen (HSF >= LSF)

condGroups = {LmoresenLSF, Lothersen};
condLabels = {'LSF > HSF', 'HSF \geq LSF'};

% LSF_PREF_s (Low switch frequency preferred-motion trials)
% HSF_PREF_s (High switch frequency preferred-motion trials)
% panelData{gg} = {LSF_PREF_s, HSF_PREF_s} for condGroups{gg}

% condGroups{gg} is ordered/sized by the full ("N") unit list, but dat here
% is incorrect trials only -- processNeuralDataErrorTrials.m re-derives
% numUnits from unique(dat.ses_ID) on that incorrect-only table, so this
% relies on every unit having >=1 incorrect trial (unique() sorts
% alphabetically either way, so order matches as long as the count does).
% If that ever fails, the mismatched-length check inside that function
% raises an explicit error rather than silently misaligning.
panelData = cell(1,2);
for gg = 1:2
    [LSF_PREF_s, HSF_PREF_s] = ...
        processNeuralDataErrorTrials(dat, slide, bin_size, start_time, end_time, ...
        condGroups{gg}, normalizationTerm);
    panelData{gg} = {LSF_PREF_s, HSF_PREF_s};
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 5-Figure Supplement 3 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the output:
% Average MT neural activity for incorrect trials (LSF and HSF switch PREF
% trials overlaid), separated into one panel per sensitivity group:
% Left: LSF > HSF sessions (LmoresenLSF). Right: HSF >= LSF sessions (Lothersen).

testing_dur = 900; % How much of the test epoch to display
end_bin = find(bins(:,1) == testing_dur - 0.5*bin_size, 1);

x_vals = xax(1:end_bin)';
x_vals_patch = [x_vals; flip(x_vals)];
x_val_trim = 0; % Restrict the plot to begin at test-stimulus onset

figure(98); clf

for gg = 1:2
    subplot(1,2,gg); hold on; box on;

    for cc = 1:2
        y_vals = mean(panelData{gg}{cc}(:,1:end_bin), 1, 'omitnan');
        y_vals_SEM = nanse(panelData{gg}{cc}(:,1:end_bin));
        y_vals_patch = [y_vals - y_vals_SEM; flip(y_vals + y_vals_SEM)];

        valid_idx = ~isnan(y_vals_patch) & x_vals_patch > x_val_trim;
        p = patch(x_vals_patch(valid_idx), y_vals_patch(valid_idx), colors{cc});
        set(p, 'LineStyle', 'none', 'FaceAlpha', 0.2);

        valid_idx = ~isnan(y_vals) & x_vals > x_val_trim;
        plot(x_vals(valid_idx), y_vals(valid_idx), 'LineWidth', 2, 'Color', colors{cc});
    end

    xlim([0 1000])
    ylim([-0.2 0.8])
    xlabel('Time from test onset (ms)')
    ylabel('Normalized firing rate')
    title(condLabels{gg})
end