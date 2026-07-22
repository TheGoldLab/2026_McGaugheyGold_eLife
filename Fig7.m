%% Fig7.m
% Compares the relative contributions of neural activity and pupil-based arousal
% to predicting context-dependent behavioral choices.
%
% Model comparison approach:
%   Three models fit for each session:
%   1. Control (behavior-only): P(switch) = logistic(β₀ + β₁*dir + β₂*time)
%   2. Neural model: P(switch) = logistic(β₀ + β₁*dir + β₂*neural×time + β₃*time)
%   3. Pupil model: P(switch) = logistic(β₀ + β₁*dir + β₂*pupil×time + β₃*time)
%
%   Neural term: Average MT firing rate during test epoch (50-1200ms)
%   Pupil term: Evoked pupil diameter 500ms before test onset (-500 to 0ms)
%
%   Model fits pre-computed by:
%   - behaviorLogisticFitsNeuralTerm.m (neural models)
%   - behaviorLogisticFitsPupilTerm.m (pupil models)
%
% Explanatory power metric:
%   ΔR² = R²_model - R²_control (Tjur's pseudo-R²)
%   Quantifies how much additional variance the neural or pupil term explains
%   beyond the baseline behavior-only model. Positive values indicate the
%   term improves prediction; larger values indicate stronger contribution.
%
% Main Text Figures:
%
%   Figure 7A: Distributions of ΔR² comparing neural vs pupil contributions.
%     Two panels (Monkey An left, Monkey Mi right).
%     
%     Each panel shows two overlaid histograms:
%     - Neural ΔR² (blue histogram with blue dashed line at mean)
%     - Pupil ΔR² (red histogram with red dashed line at mean)
%     
%   Figure 7B: Context-dependent coefficient differences (LSF-HSF).
%     Two panels (Monkey An left, Monkey Mi right).
%     
%     Scatter plots comparing interaction coefficients:
%     X-axis: Neural interaction difference (β₂_LSF - β₂_HSF from neural model)
%     Y-axis: Pupil interaction difference (β₂_LSF - β₂_HSF from pupil model)
%     
%     Each point is one session. Unity line shows where contributions are equal.
%     Points above line indicate larger pupil than neural coefficient differences.
%     
%     Different markers by monkey (circles = An, squares = Mi)
%
% Required data files (pre-computed fits):
%   - LogisticFits_NeuralTerm.mat (neural model coefficients)
%   - LogisticFits_NeuralTerm_control.mat (control for neural sessions)
%   - LogisticFits_NeuralTerm_Rsq.mat (neural model R²)
%   - LogisticFits_NeuralTerm_Rsq_control.mat (control R²)
%   - LogisticFits_PupilTerm.mat (pupil model coefficients)
%   - LogisticFits_PupilTerm_control.mat (control for pupil sessions)
%   - LogisticFits_PupilTerm_Rsq.mat (pupil model R²)
%   - LogisticFits_PupilTerm_Rsq_control.mat (control R²)
%   - sensitivity_diff_labeled_NP.mat (session classifications)
%
% Required functions:
%   - computeCohenDCI.m (calculate effect size)

%% Load data
cfg = projectDefaults();
load(fullfile(cfg.paths.data, 'mergedTable_proc_core.mat')) % behavioral-only: skips Unit_1 and pupil traces (see buildMergedTableTiers.m)
load(fullfile(cfg.paths.data, 'sensitivity_diff_labeled_NP.mat'))

load(fullfile(cfg.paths.fits, 'LogisticFits_NeuralTerm.mat'))
load(fullfile(cfg.paths.fits, 'LogisticFits_NeuralTerm_control.mat'))
load(fullfile(cfg.paths.fits, 'LogisticFits_NeuralTerm_Rsq.mat'))
load(fullfile(cfg.paths.fits, 'LogisticFits_NeuralTerm_Rsq_control.mat'))

load(fullfile(cfg.paths.fits, 'LogisticFits_PupilTerm.mat'))
load(fullfile(cfg.paths.fits, 'LogisticFits_PupilTerm_control.mat'))
load(fullfile(cfg.paths.fits, 'LogisticFits_PupilTerm_Rsq.mat'))
load(fullfile(cfg.paths.fits, 'LogisticFits_PupilTerm_Rsq_control.mat'))

% Subset data appropriately
% "NP" = Sessions with behavior and pupil data
    % Two sessions from Monkey An had no pupil data recorded
[mergedTableSub] = createDatSubset(mergedTable_proc, 'NP');
dat = mergedTableSub;
clear mergedTable_proc

% Create monkey indices for subsetting sessions:
monkeyIdx = getMonkeyIndices(dat);
An = monkeyIdx.An;
Mi = monkeyIdx.Mi;

monkeyOrder = {'An', 'Mi'};
monkeyNumIdx = struct('An', An, 'Mi', Mi);
monkeyMarkers = struct('An', 'o', 'Mi', 'square');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 7A %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate difference in Rsq (difference in explanatory power)
% Neural/pupil Rsq - control Rsq

R_sq_diff_neural = R_sq_neuralTerm - R_sq_neuralTerm_control;
R_sq_diff_pupil = R_sq_pupilTerm - R_sq_pupilTerm_control;

edges = -0.06:0.02:0.12;

% sessionDiff.(monkeyName) holds the per-session neural-pupil difference,
% used below for the An vs. Mi comparison.
sessionDiff = struct();

figure(92); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    idx = monkeyNumIdx.(monkeyName);

    subplot(1,2,mm); hold on; box on;

    x = mean(R_sq_diff_neural(idx,:), 2, 'omitnan');
    y = mean(R_sq_diff_pupil(idx,:), 2, 'omitnan');

    x_avg = mean(x, 'omitnan');
    y_avg = mean(y, 'omitnan');

    histogram(x, edges, 'Normalization', 'probability')
    xline(x_avg, '--b')
    histogram(y, edges, 'Normalization', 'probability')
    xline(y_avg, '--r')

    ylabel('Fraction of sessions')
    subtitle(monkeyName)
    ylim([0 1])

    % NOTE: previously computed via signtest but never displayed anywhere
    % (not printed, not plotted) -- the manuscript figure shows these
    % p-values in-panel, so they're added as text here.
    [p_neural, ~] = signtest(x);
    [p_pupil, ~] = signtest(y);

    text(0.95, 0.9, sprintf('Neural p = %.3g', p_neural), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'Color', 'b')
    text(0.95, 0.8, sprintf('Pupil p = %.3g', p_pupil), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'Color', 'r')

    sessionDiff.(monkeyName) = x - y;
end

[p_AnMi, h_AnMi] = ranksum(sessionDiff.An, sessionDiff.Mi);
fprintf('An vs. Mi neural-pupil difference: p = %.4g (h = %d)\n', p_AnMi, h_AnMi);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 7B %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Context-stability differences in neural/pupil terms
% Plot separately for Monkey An and Monkey Mi
% Different monkeys have different axis ranges below (An: +/-0.08 x, Mi: +/-0.06)

monkeyLim = struct('An', [-0.06 0.08], 'Mi', [-0.06 0.06]);

x_line = -0.2:0.1:2;
y_line = -0.2:0.1:2;

figure(14); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    idx = monkeyNumIdx.(monkeyName);

    subplot(1,2,mm); hold on; box on;

    lims = monkeyLim.(monkeyName);
    ylim(lims)
    xlim(lims)

    plot(x_line, y_line, 'k');

    % fits_neuralTerm/fits_pupilTerm are [numUnits x 1 x 2]; indexing with
    % only 2 subscripts collapses dims 2:3, so column 1 = LSF (hh=1) and
    % column 2 = HSF (hh=2) per cfg.hazard.codes.
    N_L = fits_neuralTerm(idx,1);
    N_H = fits_neuralTerm(idx,2);

    P_L = fits_pupilTerm(idx,1);
    P_H = fits_pupilTerm(idx,2);

    x = N_L-N_H;
    y = P_L-P_H;

    plot(x,y, 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    [p, h, stats] = signrank(x,y);
    d = computeCohenDCI(y, x, 'paired');
    [r,p_corr] = corr(x,y, 'type', 'Spearman');

    xlabel('Neural offset (avg.)')
    ylabel('Pupil offset (avg.)')
    subtitle({['p = ', num2str(p), ' & d = ', num2str(d)], ...
        ['r = ', num2str(r), ' & p = ', num2str(p_corr)]})
end