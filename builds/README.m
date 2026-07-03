% code/builds/ -- slow, run-manually scripts
%
% Everything in this folder regenerates a cached data file that's already
% checked into data/ (or data/behaviorFits/). None of it runs as part of
% normal figure generation -- the Fig*.m scripts just load the pre-computed
% outputs. Only re-run a script here if its specific upstream source data
% changes; see each entry below for what that trigger is.
%
% Both cfg.paths.functions (code/utilities/) and cfg.paths.logistic
% (code/utilities/logisticUtilities/) are added to the MATLAB path by
% cfg = projectDefaults(), which every script below calls first, so their
% dependencies resolve regardless of your current folder. code/builds/
% itself (cfg.paths.builds) is also added to the path by projectDefaults(),
% so once any script has called it, these can be invoked by name (e.g.
% typing behaviorLogisticFits at the command window) without cd'ing here
% first.

%% buildMergedTableTiers.m
    % Regenerates data/mergedTable_proc_core.mat and
    % data/mergedTable_proc_neural.mat from data/mergedTable_proc.mat by
    % dropping heavy columns (pupil traces, and for the core tier, Unit_1
    % spike times too). See "data/" in the top-level README.m for what
    % each tier is used for.
    %
    % When to run: only if data/mergedTable_proc.mat itself is
    % regenerated. The derived tiers are not meant to be hand-edited.
    % Run this before the behaviorLogisticFits*.m scripts below if
    % mergedTable_proc.mat changed, since three of the four load a tier
    % this script produces.
    % Runtime: dominated by the ~22s load of mergedTable_proc.mat.

%% behaviorLogisticFits.m
    % Fits behavioral data with a time-dependent logistic psychometric
    % function (behavior-only model). Loads mergedTable_proc_core.mat.
    % Outputs: data/behaviorFits/LogisticFits.mat, LogisticFits_Rsq.mat.
    %
    % When to run: only if the underlying behavioral data
    % (mergedTable_proc*.mat) changes.

%% behaviorLogisticFitsShuffle.m
    % Shuffle-control version of behaviorLogisticFits.m: shuffles viewing
    % durations and switch/non-switch trial types and refits, 100
    % iterations per session, to build a null distribution. Loads
    % mergedTable_proc_core.mat.
    % Outputs: data/behaviorFits/LogisticFits_Shuffle.mat,
    % LogisticFits_Rsq_Shuffle.mat.
    %
    % When to run: only if the underlying behavioral data changes.
    % Runtime: hours (100 shuffled refits per session across all sessions).

%% behaviorLogisticFitsNeuralTerm.m
    % Same time-dependent logistic model, plus a trial-wise MT neural
    % activity term (and a PREF-direction-only control fit). Loads
    % mergedTable_proc_neural.mat.
    % Outputs: data/behaviorFits/LogisticFits_NeuralTerm.mat,
    % LogisticFits_NeuralTerm_control.mat, and their _Rsq counterparts.
    %
    % When to run: only if the underlying behavioral or neural data
    % changes.

%% behaviorLogisticFitsPupilTerm.m
    % Same time-dependent logistic model, plus a trial-wise evoked-pupil
    % term (and a control fit). Loads mergedTable_proc.mat directly
    % (not a tier), since it's the only fitting script that touches pupil
    % traces.
    % Outputs: data/behaviorFits/LogisticFits_PupilTerm.mat,
    % LogisticFits_PupilTerm_control.mat, and their _Rsq counterparts.
    %
    % When to run: only if the underlying behavioral or pupil data
    % changes.

%% Dependencies
    % The four behaviorLogisticFits*.m scripts call the logist*Dotsrev*
    % fitting routines in code/utilities/logisticUtilities/, which require
    % the Optimization Toolbox and Global Optimization Toolbox. See the
    % top-level README.m "Dependencies" section for the full list,
    % including the Curve Fitting Toolbox call in
    % behaviorLogisticFitsPupilTerm.m.

%% Shared helper functions (code/utilities/logisticUtilities/)
    % buildBehaviorDataToFit.m, buildInteractionDataToFit.m, and
    % computeTjursR2.m factor out the data_to_fit construction and Tjur's
    % R^2 computation that used to be duplicated (with small
    % session/unit-specific variations) across all four
    % behaviorLogisticFits*.m scripts. Where a script restricts its
    % interaction/signed-time term to PREF-direction trials
    % (behaviorLogisticFitsNeuralTerm.m) vs. using all trials
    % (behaviorLogisticFitsPupilTerm.m), that's now an explicit
    % interactionMask argument at the call site rather than an implicit
    % difference between two copy-pasted blocks -- see "Validation" below
    % for why that distinction matters.

%% Validation (2026-07-03)
    % All four scripts were reviewed for inefficiencies, fixed, run
    % end-to-end against the real project data, and diffed against the
    % pre-computed data/behaviorFits/*.mat files checked into this repo.
    % Full methodology and numbers below; short version: every script's
    % output matches the cache to within the run-to-run variation expected
    % from GlobalSearch's own (unseeded) stochastic search -- except for
    % one real bug, described under behaviorLogisticFitsPupilTerm.m.
    %
    % Shared fix (all four scripts): logistFitDotsrev.m and
    % logistFitDotsrevNP.m called GlobalSearch with NumStageOnePoints=2000,
    % NumTrialPoints=10000 -- 10x MATLAB's own documented defaults, for a
    % smooth 4-5 parameter objective. Tuned to the defaults (200/1000).
    % Measured on one real session/hazard fit: 4.08s -> 0.97s (4.2x),
    % identical optimum (parameter diff 1.2e-4, log-likelihood diff
    % 2.7e-12). This is the dominant lever for all four scripts' runtime.
    %
    % behaviorLogisticFits.m: had a debug-only `for ss = 1 %:numSessions`
    % loop that fit only 1 of 161 sessions; restored to the full loop.
    % Full run: ~200-230s (161 sessions x 2 hazards = 322 fits). vs.
    % cached LogisticFits.mat: r=0.999998-1.000000 on all 4 parameters,
    % r=1.000000 on R^2, NaN/skip pattern matches exactly. One outlier
    % (session 125, hazard 1) on an otherwise-flat bias parameter --
    % expected given the unseeded optimizer, not a new issue.
    %
    % behaviorLogisticFitsShuffle.m: was iteration-outer/session-inner, so
    % the per-session table lookup reran on every one of 100 iterations;
    % reordered to session-outer/iteration-inner (~68s saved, on top of
    % the GlobalSearch fix). Full run (100 iterations x 161 sessions x 2
    % hazards = 32,200 fits): 9,653.8s (2.68 hr) measured, vs. an
    % estimated ~11.3 hr under the original GlobalSearch settings.
    % vs. cached LogisticFits_Shuffle.mat: this is a shuffled null
    % distribution, so rows aren't expected to line up 1:1 (and didn't,
    % by design, once the loop was reordered) -- compared distributions
    % instead. Two-sample KS tests fail to reject "same distribution" for
    % all 4 parameters and R^2 in both hazard conditions (p=0.45-1.00);
    % means and SDs match closely; 0% NaN in both.
    %
    % behaviorLogisticFitsNeuralTerm.m: investigated the per-trial
    % sliding-window spike-count double loop as a suspected bottleneck --
    % measured directly (0.064s for one unit's 191 trials; a vectorized
    % rewrite was actually slower once call overhead is counted) and left
    % unchanged. Full runs: ~295-330s (model), ~190-210s (control). vs.
    % cached data: fits_neuralTerm r=0.997780 (n=280), R_sq r=0.999723;
    % fits_neuralTerm_control r=0.988186 (n=306, 99.7% within 0.01),
    % R_sq_control r=0.998924.
    %
    % behaviorLogisticFitsPupilTerm.m:
    %   1. Vectorized the per-bin baseline regression: was calling
    %      MATLAB's regress() once per time bin (3,801 calls/unit) for
    %      what is, with one predictor and no intercept, a closed-form
    %      OLS-through-the-origin problem. Replaced with one vectorized
    %      computation; matches regress() to machine precision (~9e-16),
    %      measured 14.6x faster on that block.
    %   2. Found and fixed a real bug: the `if control` save block read
    %      fits(uu,3,:)/R_sq(uu,1,:) (only ever populated when
    %      control=false) instead of fits_control/R_sq_control. Every
    %      control run had been silently saving near-empty output --
    %      confirmed the checked-in LogisticFits_PupilTerm_control.mat was
    %      91.5% NaN (26/306 populated). Fixed, regenerated both control
    %      .mat files (now 306/306 populated), and replaced the cached
    %      files (see git history for that commit).
    %   Full runs: ~350-400s (model), ~255-290s (control, bug-fixed). vs.
    %   cached data: fits_pupilTerm R_sq r=0.999183; the fits_pupilTerm
    %   parameter itself is a weaker match overall (r=0.70) but that's
    %   driven by ~1 of 280 units landing at the opposite parameter bound
    %   -- 95.7% of units match within 0.001, 97.5% within 0.01, and a
    %   direct A/B test of original-vs-tuned GlobalSearch settings on a
    %   controlled fit converged identically and reproducibly under both,
    %   ruling out the GlobalSearch tuning as the cause. This is the
    %   known signature of an unseeded global optimizer on a
    %   weakly-identified interaction term for a handful of units,
    %   present in the original pipeline regardless of these changes.
    %   fits_pupilTerm_control (bug-fixed) matches the now-fixed cache
    %   exactly: 0 diff across all 306 entries.
    %
    % Modularization: the data_to_fit construction and R^2 computation
    % above were extracted into buildBehaviorDataToFit.m /
    % buildInteractionDataToFit.m / computeTjursR2.m (see "Shared helper
    % functions" above). Validated two ways: (1) isequaln comparison
    % against the original inline logic across 33 sampled real sessions,
    % covering all four call patterns (masked/unmasked, behavior-only/
    % interaction-term) -- exact match in every case; (2)
    % behaviorLogisticFits.m and both control states of
    % behaviorLogisticFitsNeuralTerm.m / behaviorLogisticFitsPupilTerm.m
    % were re-run end-to-end post-refactor and produced results identical
    % to their pre-refactor runs (same correlations, same diffs, in two
    % cases -- PupilTerm control -- exactly 0 diff). behaviorLogisticFitsShuffle.m
    % was not re-run post-refactor (a full run costs ~2.7 hours); its
    % call pattern is covered by the isequaln test in (1), and its own
    % loop-reorder fix was independently validated via the full run
    % documented above.
