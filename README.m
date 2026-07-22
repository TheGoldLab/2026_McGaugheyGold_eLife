
% Data and Code for: McGaughey & Gold (2025)
% "Sensory adaptation supports flexible evidence accumulation during perceptial decision making"

% Data and code are available at: https://upenn.box.com/s/uvw64y7xu9zarn3ajd05igmya3a8893l
% Code is available at: https://github.com/TheGoldLab/2026_McGaugheyGold_eLife

% Contact: mckar@pennmedicine.upenn.edu or jigold@pennmedicine.upenn.edu
% Last updated: 07/03/2026

%% Dependencies

% MATLAB toolboxes required:
%   - Statistics and Machine Learning Toolbox (ranksum, signrank, corr,
%     fitlm, ttest, ttest2, and other statistics calls throughout)
%   - Signal Processing Toolbox (butter -- pupil low-pass filtering)
%   - Optimization Toolbox and Global Optimization Toolbox (used inside
%     the logist*Dotsrev* fitting routines called by the
%     behaviorLogisticFits*.m scripts)
%   - Curve Fitting Toolbox (smooth(), fit() -- used in Fig3_FigSupp1.m,
%     processPupilBaseline.m, processPupilFull.m, and
%     behaviorLogisticFitsPupilTerm.m)

% External code (not included in this repo):
%   Two small utility functions come from the Gold lab's shared MATLAB
%   utilities repository, https://github.com/TheGoldLab/Lab_Matlab_Utilities
%   (tested against commit c00a356, 2024-06-13; both live in its general/
%   subfolder):
%       - nanrunmean.m
%       - nans.m
%   Clone that repo and add general/ to the MATLAB path before running
%   anything here. Referencing rather than vendoring these keeps this repo
%   scoped to paper-specific code and picks up upstream fixes automatically;
%   the tradeoff is an extra setup step (clone + addpath) and a small
%   dependency on the external repo keeping those two filenames stable,
%   which is why the commit is pinned above.
%
%   code/utilities/nanse.m is the one exception: it's a modified copy of
%   that same repo's general/nanse.m, vendored (not referenced) because it
%   needed a small compatibility fix for newer MATLAB releases -- see the
%   comment at the top of the file for what changed and why.
%
%   code/utilities/computeCohenDCI.m (Cohen's d, plus a confidence
%   interval and standard error) is also vendored directly rather than
%   referenced, since it wasn't found in Lab_Matlab_Utilities either. It
%   is a strict superset of the older computeCohenD.m that used to live
%   in this repo (identical d calculation, verified against sample data
%   before consolidating) -- every call site now uses this one function,
%   and computeCohenD.m was removed.
%
%   code/utilities/logisticUtilities/logistFitDotsrevUD.m and
%   logistValDotsrevUD.m (a 3-parameter logistic model variant with no
%   time-dependent term) were removed: nothing in this repo called them,
%   they weren't documented anywhere, and they still carried the
%   un-tuned GlobalSearch settings that were fixed everywhere else (see
%   code/builds/README.m "Validation") -- signs they were leftover from
%   an earlier model version. See git history if they're ever needed.

%% Project layout
% code/                   - all MATLAB scripts and functions
%     Fig*.m               - one script per main-text figure (see "Figure scripts" below)
%     builds/              - slow, run-manually scripts that regenerate cached
%                             data (behaviorLogisticFits*.m, buildMergedTableTiers.m);
%                             see code/builds/README.m
%     utilities/           - shared helper functions (see below)
%         logisticUtilities/ - logistFit/logistVal/logistErr functions the
%                             builds/behaviorLogisticFits*.m scripts call
% data/                    - all data assets (see "data/" below)
%     behaviorFits/        - pre-computed logistic fits (see below)

%% code/builds/behaviorLogisticFits.m
    % Includes code for fitting behavioral data with time-dependent logistic function
    % code/utilities/logisticUtilities/ contains functions necessary for fitting
    % As fitting is time intensive, all relevant output is included in data/behaviorFits/

%% code/builds/behaviorLogisticFitsShuffle.m
    % Includes code for fitting behavioral data with shuffled test-stimulus durations and switch/non-switch trial types.
    % Sessions with fit over 100 shuffled iterations (hours of runtime)
    % All relevant output is included in data/behaviorFits/

%% code/builds/behaviorLogisticFitsNeuralTerm.m
    % Includes code for fitting behavioral data with time-dependent logistic function
    % Psychometric function includes a trial-wise neural term (average MT neural activity during the test epoch with 50 ms delay)
    % code/utilities/logisticUtilities/ contains functions necessary for fitting
    % As fitting is time intensive, all relevant output is included in data/behaviorFits/

%% code/builds/behaviorLogisticFitsPupilTerm.m
    % Includes code for fitting behavioral data with time-dependent logistic function
    % Psychometric function includes a trial-wise pupil term (average evoked pupil 500 ms before test-stim onset)
    % code/utilities/logisticUtilities/ contains functions necessary for fitting
    % As fitting is time intensive, all relevant output is included in data/behaviorFits/

%% code/builds/ (folder)
    % Slow, run-manually scripts that regenerate cached data used by the
    % figure scripts (the four behaviorLogisticFits*.m scripts above, plus
    % buildMergedTableTiers.m). Nothing in the figure-generation path calls
    % these directly -- their outputs are pre-computed and checked into
    % data/. See code/builds/README.m for what each script does and when
    % it needs to be re-run.

%% data/behaviorFits/ (folder)
    % Contains all fits, which are loaded in by other scripts to save time/computing
        % Behavior-only (LogisticFits.mat)
        % Behavior-only shuffled (LogisticFits_Shuffle.mat)
        % Behavior + neural term (LogisticFits_NeuralTerm.mat)
        % Behavior-neural term control (LogisticFits_NeuralTerm_control.mat)
        % Behavior + pupil term (LogisticFits_PupilTerm.mat)
        % Behavior-pupil term control (LogisticFits_PupilTerm_control.mat)

    % Contains all Tjur's pseudo-R squared values, which are loaded in by other scipts
        % Behavior-only (LogisticFits_Rsq.mat)
        % Behavior-only shuffled (LogisticFits_Rsq_Shuffle.mat)
        % Behavior + neural term (LogisticFits_NeuralTerm_Rsq.mat)
        % Behavior-neural term control (LogisticFits_NeuralTerm_Rsq_control.mat)
        % Behavior + pupil term (LogisticFits_PupilTerm_Rsq.mat)
        % Behavior-pupil term control (LogisticFits_PupilTerm_Rsq_control.mat)

%% data/
    % mergedTable_proc.mat: trial-wise behavioral, neural, and pupil data
    %     for Monkeys An, Ch, and Mi (80,314 rows x 32 columns). This is
    %     the canonical source table -- 9.5 GB once loaded, because three
    %     columns (pupil_diam, pupil_horiz, pupil_vert) are continuous
    %     per-trial pupil traces and together account for 98.7% of its
    %     size. Only scripts that actually touch pupil traces load this
    %     file directly (Fig6_FigSupp1.m, Fig6_FigSupp2AB.m,
    %     Fig6_FigSupp2CD.m, behaviorLogisticFitsPupilTerm.m, and the
    %     processPupilBaseline.m and processPupilFull.m functions they call).
    %
    % mergedTable_proc_core.mat / mergedTable_proc_neural.mat: lighter
    %     derived caches of the same table, used by every other script so
    %     they don't pay the ~22s cost of loading pupil traces they never
    %     read. Both store their table under the same variable name,
    %     mergedTable_proc, so no other code needs to change based on
    %     which file was loaded.
    %       - core   (~136 MB): drops Unit_1 (spike times) and the pupil
    %                 trace columns. Used by purely behavioral scripts
    %                 (Fig2_FigSupp1.m, Fig2_FigSupp2_3.m, Fig7.m,
    %                 behaviorLogisticFits.m, behaviorLogisticFitsShuffle.m).
    %       - neural (~220 MB): drops only the pupil trace columns, keeps
    %                 Unit_1. Used by scripts that call processNeuralData.m
    %                 or processNeuralDataErrorTrials.m (Fig3_FigSupp1.m,
    %                 Fig4_FigSupp1A/1B/2.m, Fig5ABC_FigSupp1.m,
    %                 Fig5DE_FigSupp2_4.m, Fig5_FigSupp3.m,
    %                 behaviorLogisticFitsNeuralTerm.m).
    %     Both are generated from mergedTable_proc.mat by
    %     code/builds/buildMergedTableTiers.m and are not meant to be
    %     hand-edited; re-run that script if mergedTable_proc.mat itself
    %     is ever regenerated. See the comment at the top of that script
    %     for how "heavy" columns are identified.

    % normalizationTerm.mat: term used for normalization in session-wise neural analysis
        % Needed for Figure 4-Figure Supplement 1B (early vs. late block splits)
        % Maintains consistent normalization when analyzing early vs. late halves of neural data blocks

    % sensitivity_diff_labeled_*.mat: psychometric slope (i.e., behavioral sensitivity)
        % N: Neural data subset (used for neural-only analyses)
        % NP: Neural and pupil data subset (used for neural-pupil analyses)
        % BP: Behavior and pupil data subset (used for behavior-pupil analyses)

%% code/utilities/:
% Core analysis functions used across multiple scripts:

    % projectDefaults.m: Single source of truth for paths, colors, hazard
    %     (switch-frequency) codes, and viewing-duration bins shared across
    %     every script. Every script starts with: cfg = projectDefaults();
    %     Paths are derived from this file's own location (code/utilities/),
    %     so nothing is hardcoded to a particular machine/username.

    % getMonkeyIndices.m: Per-monkey (An/Ch/Mi) row indices/logical masks,
    %     derived from session or unit ID strings (e.g. unique(dat.ses_ID))
    %     rather than hardcoded numeric ranges, since those ranges differ
    %     across data subsets (behavioral vs. neural vs. pupil).
    %     Usage: idx = getMonkeyIndices(dat); fits(idx.An,:,:)
    %     NOTE: Monkey Mi's session/unit IDs use the raw code "MM", not
    %     "Mi" -- see the comment inside this file.

    % createDatSubset: Creates data subsets (B=behavioral, N=neural, P=pupil)
        % Example: createDatSubset(data, 'N') returns neural data subset

    % behaviorBinnedPerformance.m: Bins trials by viewing duration (edges in
    %     cfg.bins.viewDuration) and returns percent-correct per bin, per
    %     session, separately for LSF/HSF (and switch/non-switch via its
    %     trialType argument). Used by Fig2_FigSupp2_3.m, Fig5ABC_FigSupp1.m,
    %     and Fig5DE_FigSupp2_4.m.

    % Behavioral analysis (Logistic regression functions):
        % logistFitDotsrev.m & logistValDotsrev.m: Behavior-only model
        % logistFitDotsrevNP.m & logistValDotsrevNP.m: Behavior + neural/pupil model
        % logistErrDotsrev.m: Error function for both models

    % Neural processing:
        % processNeuralData.m: function that takes in data subset and outputs firing rate, ROC area, etc.
        % processNeuralDataErrorTrials.m: function that processes incorrect trials, averaging across sessions
        % rocNFullOutput: function for ROC area analysis

    % Pupil processing:
        % processPupilBaseline.m: function for baseline pupil analysis (Figure 6-Figure Supplement 1)
        % processPupilFull.m: function for pupil processing and regression analysis (main text Fig. 6)

    % computeCohenDCI.m: Calculate Cohen's d effect size (plus 95% CI and
    %     SE) for paired or independent comparisons. See "Dependencies" above.
        % Usage: d = computeCohenDCI(x, y, 'paired') or computeCohenDCI(x, y, 'independent')
        %        [d, CI, SE] = computeCohenDCI(x, y, 'paired') for the full output

%% Figure scripts
% Script names follow the published manuscript's "Figure X-Figure Supplement Y"
% numbering. When possible, a figure supplement is generated by the same
% script as its related main-text figure (data loading/processing is shared).
%
%   Fig2_FigSupp1.m        - Figure 2-Figure Supplement 1
%   Fig2_FigSupp2_3.m      - Figure 2B-D, Figure 2-Figure Supplements 2 & 3
%   Fig3_FigSupp1.m        - Figure 3A-D, Figure 3-Figure Supplement 1 (A,B)
%   Fig4_FigSupp1A.m       - Figure 4A-F, Figure 4-Figure Supplement 1A
%   Fig4_FigSupp1B.m       - Figure 4-Figure Supplement 1B (run twice: Early_Block, Late_Block)
%   Fig4_FigSupp2.m        - Figure 4-Figure Supplement 2
%   Fig5ABC_FigSupp1.m     - Figure 5A-C, Figure 5-Figure Supplement 1 (A,B)
%   Fig5DE_FigSupp2_4.m    - Figure 5D-E, Figure 5-Figure Supplements 2 & 4
%   Fig5_FigSupp3.m        - Figure 5-Figure Supplement 3
%   Fig6_FigSupp1.m        - Figure 6-Figure Supplement 1 (A,B,C)
%   Fig6_FigSupp2AB.m      - Figure 6A-D, Figure 6-Figure Supplement 2 (A,B)
%   Fig6_FigSupp2CD.m      - Figure 6-Figure Supplement 2 (C,D)
%   Fig7.m                 - Figure 7A-B
%
% Nine near-duplicate scripts left over from an earlier renumbering pass
% (same content, only the Extended-Data figure number in the comments
% differed) were removed; see git history if the old copies are needed.

