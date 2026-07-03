function cfg = projectDefaults()
%PROJECTDEFAULTS Centralized paths, colors, and task constants shared
% across the McGaughey & Gold (2026) analysis scripts.
%
% cfg = projectDefaults() returns a struct with fields:
%
%   paths.root       - project root (parent of code/ and data/)
%   paths.code       - code/
%   paths.data       - data/                        (mergedTable_proc.mat, etc.)
%   paths.fits       - data/behaviorFits/            (pre-computed logistic fits)
%   paths.functions  - code/utilities/
%   paths.logistic   - code/utilities/logisticUtilities/
%   paths.builds     - code/builds/                  (see code/builds/README.m)
%
%   colors.LSF / colors.HSF   - switch-frequency (context-stability) condition colors
%   colors.pair               - {colors.LSF, colors.HSF}, indexed as co{hh}
%
%   hazard.LSF / hazard.HSF   - HR codes used in mergedTable_proc (2 = LSF, 6 = HSF)
%   hazard.codes              - [hazard.LSF hazard.HSF], indexed as uniqueHazards(hh)
%
%   bins.viewDuration.edges      - [100 225 375 600 1200] ms
%   bins.viewDuration.midpoints  - [162.5 300 487.5 900] ms, for plotting binned performance
%
%   monkeys.names    - {'An','Ch','Mi'}
%   monkeys.markers  - {'o','diamond','square'}, matched to monkeys.names
%
% projectDefaults() also adds paths.functions, paths.logistic, and
% paths.builds to the MATLAB path, so scripts no longer need their own
% addpath() calls.
%
% Paths are derived from this file's own location, so the project can
% live anywhere on disk (no per-machine hardcoded paths).
%
% Usage:
%   cfg = projectDefaults();
%   cd(cfg.paths.data); load('mergedTable_proc.mat')
%   plot(x, y, 'Color', cfg.colors.LSF)

thisFile = mfilename('fullpath');
functionsDir = fileparts(thisFile);

cfg.paths.functions = functionsDir;
cfg.paths.code = fileparts(functionsDir);
cfg.paths.root = fileparts(cfg.paths.code);
cfg.paths.data = fullfile(cfg.paths.root, 'data');
cfg.paths.fits = fullfile(cfg.paths.data, 'behaviorFits');
cfg.paths.logistic = fullfile(functionsDir, 'logisticUtilities');
cfg.paths.builds = fullfile(cfg.paths.code, 'builds');

addpath(cfg.paths.functions);
addpath(cfg.paths.logistic);
addpath(cfg.paths.builds);

% Switch-frequency (context-stability) condition colors:
%   LSF = dark blue/teal, HSF = yellow/orange (matches all published figures)
cfg.colors.LSF = [0.007843137255, 0.1882352941, 0.2784313725];
cfg.colors.HSF = [1, 0.7176470588, 0.01176470588];
cfg.colors.pair = {cfg.colors.LSF, cfg.colors.HSF};

% Hazard-rate (switch-frequency) codes used in mergedTable_proc.HR
cfg.hazard.LSF = 2;
cfg.hazard.HSF = 6;
cfg.hazard.codes = [cfg.hazard.LSF, cfg.hazard.HSF];

% Viewing-duration bins for binned behavioral-performance analyses
% (see behaviorBinnedPerformance.m)
cfg.bins.viewDuration.edges = [100 225 375 600 1200];
cfg.bins.viewDuration.midpoints = [162.5 300 487.5 900];

% Per-monkey plotting markers, used consistently across all figures
cfg.monkeys.names = {'An', 'Ch', 'Mi'};
cfg.monkeys.markers = {'o', 'diamond', 'square'};

end
