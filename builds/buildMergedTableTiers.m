function buildMergedTableTiers()
%buildMergedTableTiers Regenerate the trimmed mergedTable_proc caches.
%
% mergedTable_proc.mat is a single 80,314 x 32 table dominated by three
% pupil-trace columns, which together account for 98.7% of its size (9.4
% of 9.5 GB in memory). Loading the full table costs ~22s regardless of
% whether a script ever touches those columns. Since MATLAB has no way to
% partially load one column of a serialized table object, this script
% pre-splits the table once into two lighter derived caches so scripts
% that don't need pupil traces (or don't need spike times either) can
% skip loading them:
%
%   mergedTable_proc_core.mat   - drops Unit_1 and any column over
%                                  500 MB (behavioral analyses only)
%   mergedTable_proc_neural.mat - drops only columns over 500 MB, keeps
%                                  Unit_1 (neural analyses that don't
%                                  touch pupil)
%   mergedTable_proc.mat        - unchanged; still the only source for
%                                  any script that touches pupil traces
%                                  (processPupilBaseline.m,
%                                  processPupilFull.m,
%                                  behaviorLogisticFitsPupilTerm.m, and
%                                  anything that calls them)
%
% Columns to drop are identified by measured size (>500 MB) rather than
% by hardcoded name strings, because at least one pupil column's real
% name has incidental leading whitespace that doesn't match its expected
% literal ("pupil_vert" vs the actual stored name) -- matching on
% measured size sidesteps that entirely, and Unit_1 is looked up by
% trimmed name match since it's small.
%
% Each derived file stores its table under the same variable name,
% mergedTable_proc, so callers only need to change which filename they
% load -- createDatSubset.m and everything downstream is unaffected.
%
% Run this manually whenever mergedTable_proc.mat itself is regenerated;
% the derived caches are not meant to be hand-edited.

cfg = projectDefaults();
cd(cfg.paths.data)

fprintf('Loading mergedTable_proc.mat (this is the slow, ~22s step)...\n');
t0 = tic;
load('mergedTable_proc.mat', 'mergedTable_proc');
fprintf('  loaded in %.1f s\n', toc(t0));

allNames = mergedTable_proc.Properties.VariableNames;
colBytes = nan(size(allNames));
for k = 1:numel(allNames)
    col = mergedTable_proc.(allNames{k}); %#ok<NASGU> -- read via whos('col') below
    w = whos('col');
    colBytes(k) = w.bytes;
end

heavyColIdx = find(colBytes > 500e6);
unitColIdx = find(strcmp(strtrim(allNames), 'Unit_1'));

assert(numel(heavyColIdx) == 3, ...
    'buildMergedTableTiers:unexpectedHeavyColumns', ...
    'Expected exactly 3 columns over 500MB (the pupil traces), found %d. Aborting rather than guessing.', numel(heavyColIdx));
assert(isscalar(unitColIdx), ...
    'buildMergedTableTiers:unitColumnNotFound', ...
    'Expected exactly 1 column named Unit_1, found %d. Aborting rather than guessing.', numel(unitColIdx));

fprintf('Dropping as "heavy" (>500MB): %s\n', strjoin(allNames(heavyColIdx), ', '));

% --- neural tier: drop only the heavy (pupil) columns ---
neuralTable = mergedTable_proc;
neuralTable(:, heavyColIdx) = [];
mergedTable_proc = neuralTable; %#ok<NASGU> -- saved under this name on purpose
save(fullfile(cfg.paths.data, 'mergedTable_proc_neural.mat'), 'mergedTable_proc', '-v7.3');
fprintf('Wrote mergedTable_proc_neural.mat: %d x %d\n', size(neuralTable,1), size(neuralTable,2));

% --- core tier: also drop Unit_1 ---
coreTable = neuralTable;
coreTable(:, strcmp(strtrim(neuralTable.Properties.VariableNames), 'Unit_1')) = [];
mergedTable_proc = coreTable;
save(fullfile(cfg.paths.data, 'mergedTable_proc_core.mat'), 'mergedTable_proc', '-v7.3');
fprintf('Wrote mergedTable_proc_core.mat: %d x %d\n', size(coreTable,1), size(coreTable,2));

fprintf('Done. mergedTable_proc.mat itself was not modified.\n');

end
