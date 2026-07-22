%% Fig6_FigSupp2CD.m
% Produces Figure 6-Figure Supplement 2, panels C and D (panels A and B of
% that same supplement are produced separately, by Fig6_FigSupp2AB.m; see
% the figure-numbering note at the top of that script). Figure 6-Figure
% Supplement 2 was formerly numbered Figure 7-Figure Supplement 1.
%
% Analyzes relationships between pupil-derived behavioral stability coefficients
% and MT neural activity/selectivity across time.
%
% Analysis approach:
%   - Processes neural data (processNeuralData.m) for correct trials
%   - Extracts MT activity (50-500ms) and ROC area for LSF vs HSF conditions
%   - Processes pupil data (processPupilFull.m) to get regression coefficients
%   - Computes time-resolved Spearman correlations between:
%       A) Pupil switch-frequency β coefficient and MT neural difference (LSF-HSF firing rate)
%       B) Pupil switch-frequency β coefficient and MT selectivity difference (LSF-HSF ROC area)
%   - Analyzes Monkeys An and Mi separately, plus combined
%
% Outputs:
%   Figure 6-Figure Supplement 2C: Time course of correlations between pupil switch-frequency β
%     and MT neural difference (firing rate). Shows:
%     - Monkey Mi (red line)
%     - Monkey An (gray line)
%     - Combined (black line)
%     - Significance markers (p<0.05) for each
%
%   Figure 6-Figure Supplement 2D: Time course of correlations between pupil β
%     and MT selectivity difference (ROC area). Same format as panel C.
%
% Required data files:
%   - mergedTable_proc.mat
%   - sensitivity_diff_labeled_NP.mat (sessions with both neural and pupil)
%
% Required functions:
%   - processNeuralData.m (extract MT activity and selectivity)
%   - processPupilFull.m (extract pupil regression coefficients)

%% Load data
cfg = projectDefaults();

load(fullfile(cfg.paths.data, 'mergedTable_proc.mat'))
load(fullfile(cfg.paths.data, 'sensitivity_diff_labeled_NP.mat')) % Want sessions with neural and pupil data

% Subset data appropriately
% "NP" = Sessions with neural and pupil data
    % Two sessions from Monkey An had no pupil data recorded
[mergedTableSub] = createDatSubset(mergedTable_proc, 'NP');
dat = mergedTableSub;
clear mergedTable_proc

% Select correct trials
datCorrect = dat(dat.correct == 1,:);
dat = datCorrect;

% Create monkey indices for subsetting units.
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
% Process neural data calling processNeuralData.m
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

[spike_rates_avg, spike_rates_SEM, spike_rates_avg_std_norm_avg, ...
    spike_rates_avg_std_norm_SEM, ...
    ROC_area, cell_selectivity, ...
    normalizationTerm, uniqueUnitNames, ...
    unit_example] = processNeuralData(dat, slide, bin_size, start_time, end_time, [], []);

% Calculate LSF and HSF neural activity and ROC area during the testing stimulus:

bin1 = 50;  % Starting bin for response average
bin2 = 500; % Ending bin for response average

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

HSF_neural = mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,3), 2, 'omitnan');
LSF_neural = mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,1), 2, 'omitnan');

HSF_ROCa = mean(ROC_area(:,start_bin:end_bin,3), 2, 'omitnan');
LSF_ROCa = mean(ROC_area(:,start_bin:end_bin,2), 2, 'omitnan');

% Take the difference (to be compared with the difference in context-dependent B coefficient, which is LSF-HSF evoked pupil response)

neural_diff = LSF_neural - HSF_neural;
ROCa_diff = LSF_ROCa - HSF_ROCa;

%%
% Process pupil data calling processPupilFull.m
% Set up timing/binning:

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6-Figure Supplement 2C-D %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot correlation between context-stability beta coefficient as a function
% of viewing time and (C) test-stim neural activity difference, (D)
% test-stim ROC-area difference (both LSF-HSF). Correlations run separately
% for Monkey An and Monkey Mi, plus combined.

xax_plot = xax - 2400; % Align to testing stimulus onset, not adapting

panelDiff = {neural_diff, ROCa_diff};
panelFig = [101 108];

for pp = 1:2
    thisDiff = panelDiff{pp};

    % Loop through bins and calculate correlation
    HRpupil_behav_corrcoef_Mi = nans(1,numBins);
    HRpupil_behav_corrcoef_Mi_p = nans(1,numBins);

    HRpupil_behav_corrcoef_An = nans(1,numBins);
    HRpupil_behav_corrcoef_An_p = nans(1,numBins);

    HRpupil_behav_corrcoef = nans(1,numBins);
    HRpupil_behav_corrcoef_p = nans(1,numBins);

    % fits is [numSessions x numBins x 3 params x 2 measures] (processPupilFull.m):
    % dim3=3 selects the hazard-condition beta coefficient, dim4=1 selects
    % the mean-pupil regression (dim4=2 would be the slope).
    for bb = 1:numBins

        [r_Mi,p_Mi] = corr(fits(Mi,bb,3,1), thisDiff(Mi), 'Type', 'Spearman', 'rows', 'complete');

        [r_An,p_An] = corr(fits(An,bb,3,1), thisDiff(An), 'Type', 'Spearman', 'rows', 'complete');

        [r,p] = corr(fits([An;Mi],bb,3,1), thisDiff([An;Mi]), 'Type', 'Spearman', 'rows', 'complete');

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

    figure(panelFig(pp)); clf; box on; hold on;

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
end
