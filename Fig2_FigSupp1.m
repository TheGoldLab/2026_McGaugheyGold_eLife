%% Fig2_FigSupp1.m
% Compares goodness-of-fit (Tjur's pseudo-R²) for behavioral logistic models
% fitted to real vs. shuffled data, separated by subject.
%
% Analysis:
%   - Loads behavioral fits and shuffle control fits from data/behaviorFits/
%   - Separates sessions by subject (Monkeys An, Mi, Ch)
%   - For each subject:
%       * Plots histogram of R² values for real fits (blue)
%       * Overlays histogram of R² values for shuffled fits (cyan)
%       * Shows mean R² for each distribution (red = real, black = shuffled)
%       * Performs Wilcoxon rank-sum test comparing distributions
%
% Outputs:
%   Figure 2-Figure Supplement 1: Three-panel figure showing R² distributions
%   for each monkey, with statistical comparison of real vs. shuffled fits.
%
% Required data files:
%   - LogisticFits.mat (real behavioral fits)
%   - LogisticFits_Rsq.mat (R² for real fits)
%   - LogisticFits_Shuffle.mat (shuffled control fits)
%   - LogisticFits_Rsq_Shuffle.mat (R² for shuffled fits)

%% Load data
cfg = projectDefaults();

cd(cfg.paths.data)
load('mergedTable_proc_core.mat') % behavioral-only: skips Unit_1 and pupil traces (see buildMergedTableTiers.m)

% Subset data appropriately
% "B" = Behavioral analysis
[mergedTableSub] = createDatSubset(mergedTable_proc, 'B');
dat = mergedTableSub;

cd(cfg.paths.fits)

% Load pseudo R^2 for logistic fit assessment
load('LogisticFits.mat')
load('LogisticFits_Rsq.mat')

load('LogisticFits_Shuffle.mat')
load('LogisticFits_Rsq_Shuffle.mat')

%% Index fits by animal
uniqueSessionNames = unique(dat.ses_ID);
numIterations = 100;

monkeyIdx = getMonkeyIndices(uniqueSessionNames);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 2-Figure Supplement 1 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the distribution of pseudo-R^2 values for real fits vs. shuffled
% fits, one panel per monkey, preserving the An/Mi/Ch panel order used
% throughout the manuscript (not the alphabetical An/Ch/Mi order that
% monkeyIdx.names uses).
monkeyOrder = {'An', 'Mi', 'Ch'};

figure(19); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    isMonkey = monkeyIdx.(monkeyName);

    % Shuffled fits are stacked as numIterations repeats of the
    % per-session block (see behaviorLogisticFitsShuffle.m), so the
    % per-session monkey mask just needs to be tiled to match.
    isMonkeyShuffle = repmat(isMonkey, numIterations, 1);

    subplot(1, 3, mm); hold on; box on;

    % Real fits
    histogram(R_sq(isMonkey,:), 'Normalization', 'probability', 'BinWidth', 0.05)
    xline(mean(R_sq(isMonkey,:), 'all'), 'Color', 'r')

    % Shuffled fits
    histogram(R_sq_shuffle(isMonkeyShuffle,:), 'Normalization', 'probability', 'FaceColor', 'c', 'BinWidth', 0.05)
    xline(mean(R_sq_shuffle(isMonkeyShuffle,:), 'all'), 'Color', 'k')

    ylim([0 1])
    xlim([0 1])
    ylabel('Fraction of sessions')
    xlabel("Tjur's pseudo R^2")
    title(monkeyName)

    % Compare real vs. shuffled distributions, averaging each session's
    % (or shuffle iteration's) R^2 across both hazard conditions
    % (columns 1 and 2).
    R_sq_shuffle_avg = mean(R_sq_shuffle(isMonkeyShuffle,1:2), 2);
    R_sq_avg = mean(R_sq(isMonkey,1:2), 2);

    [p, h, stats] = ranksum(R_sq_shuffle_avg, R_sq_avg);
    fprintf('%s: p = %.3g, h = %d, zval = %.3f\n', monkeyName, p, h, stats.zval);
end
