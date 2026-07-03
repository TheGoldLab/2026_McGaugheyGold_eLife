%% Fig4_FigSupp1B.m
% Analyzes adaptation of MT neural selectivity across blocks by comparing
% early vs. late responses within switch-frequency blocks.
%
% Block splitting procedure:
%   - Identifies transitions between LSF (HR=2) and HSF (HR=6) blocks
%   - Splits each block into early (first half) and late (second half) trials
%   - Pools trials across all early blocks and all late blocks separately
%
% Neural analysis:
%   - Processes neural data using processNeuralData.m for early and late blocks
%   - Uses normalization terms from full session analysis (normalizationTerm.mat)
%     to maintain consistent scaling between early/late comparisons
%   - Extracts average firing rate during first pref-motion stimulus presentation (-2350 to -2000ms)
%   - Computes LSF-HSF difference for each block half
%
% Outputs:
%   Figure 4-Figure Supplement 1B: Scatter plot comparing early vs. late block
%   LSF-HSF differences for each neuron (Monkey Mi only).
%
% Statistics:
%   Wilcoxon signed-rank test (paired) comparing early vs. late differences
%   Cohen's d effect size for paired comparison
%
% Note: Analysis restricted to Monkeys An and Mi (Ch excluded).
%
% Required data files:
%   - mergedTable_proc.mat
%   - sensitivity_diff_labeled_N.mat
%   - normalizationTerm.mat (critical for consistent early/late comparison)

%% Load data
cfg = projectDefaults();
cd(cfg.paths.data)

load('mergedTable_proc_neural.mat') % keeps Unit_1, skips pupil traces (see buildMergedTableTiers.m)
load('sensitivity_diff_labeled_N.mat') % Neural subset
load('normalizationTerm.mat') % Normalization term from session-wise neural analysis

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
Ch = find(monkeyIdx.Ch);

normalizationTerm(Ch) = []; % Remove Monkey Ch, who is not used for block analysis
normalizationTerm_input = normalizationTerm;

% Get unique session (unit) names:
uniqueUnitNames = unique(dat.ses_ID);
numUnits = length(uniqueUnitNames);

% Save tables for early and late blocks, separately for LSF and HSF
LSF_Early = [];
LSF_Late = [];
HSF_Early = [];
HSF_Late = [];

%%
% Split blocks

for uu = 1:numUnits

    % Show session name
    disp(uniqueUnitNames{uu})

    if contains(uniqueUnitNames{uu}, "Ch")
        continue
    end

    Lses = strcmp(dat.ses_ID, uniqueUnitNames{uu}) & dat.coh_final >= 50;

    session_data = dat(Lses, 1:end);

    num_trials = size(session_data,1);

    LH_idx = strfind(session_data.HR',[2 6]); % trial indices where an LSF block ends and an HSF block begins
    HL_idx = strfind(session_data.HR',[6 2]); % trial indices where an HSF block ends and an LSF block begins

    num_blocks = size(LH_idx,2) + size(HL_idx,2) + 1;

    % Number of blocks ranges from 2 to 5. The num_blocks branches below
    % assume sessions alternate strictly between LSF/HSF blocks (so exactly
    % one of LH_idx/HL_idx is ever populated where a branch expects only
    % one) -- there's no else/error fallback if that invariant is violated.
    % Find which block is first and create four tables
    % LSF early
    % LSF late
    % HSF early
    % HSF late

    if num_blocks == 2

        if isempty(LH_idx)
            Lfirst = false;
        elseif isempty(HL_idx)
            Lfirst = true;
        end

        % Organize/distribute trials based on which condition came first

        if Lfirst

            L_split = round(LH_idx/2,0);

            L_Early = session_data(1:L_split,:);
            L_Late = session_data(L_split+1:LH_idx,:);

            LSF_Early = [LSF_Early; L_Early];
            LSF_Late = [LSF_Late;L_Late];

            H_split = round((num_trials - LH_idx)/2,0);

            H_Early = session_data(LH_idx+1:LH_idx+H_split,:);
            H_Late = session_data(LH_idx+H_split+1:end,:);

            HSF_Early = [HSF_Early; H_Early];
            HSF_Late = [HSF_Late; H_Late];

        elseif ~Lfirst

            H_split = round(HL_idx/2,0);

            H_Early = session_data(1:H_split,:);
            H_Late = session_data(H_split+1:HL_idx,:);

            HSF_Early = [HSF_Early; H_Early];
            HSF_Late = [HSF_Late; H_Late];

            L_split = round((num_trials - HL_idx)/2,0);

            L_Early = session_data(HL_idx+1:HL_idx+L_split,:);
            L_Late = session_data(HL_idx+L_split+1:end,:);

            LSF_Early = [LSF_Early; L_Early];
            LSF_Late = [LSF_Late;L_Late];
        end

    elseif num_blocks == 3

        if LH_idx < HL_idx
            Lfirst = true;
        elseif HL_idx < LH_idx
            Lfirst = false;
        end

        % Organize/distribute trials based on which condition came first

        if Lfirst

            L_split_1 = round(LH_idx/2,0);
            L_split_2 = round((num_trials - HL_idx)/2,0);

            L_Early_1 = session_data(1:L_split_1,:);
            L_Late_1 = session_data(L_split_1+1:LH_idx,:);

            L_Early_2 = session_data(HL_idx+1:HL_idx+L_split_2,:);
            L_Late_2 = session_data(HL_idx+L_split_2+1:num_trials,:);

            LSF_Early = [LSF_Early; L_Early_1; L_Early_2];
            LSF_Late = [LSF_Late; L_Late_1; L_Late_2];

            H_split = round((HL_idx - LH_idx)/2,0);

            H_Early = session_data(LH_idx+1:LH_idx+H_split,:);
            H_Late = session_data(LH_idx+H_split+1:HL_idx,:);

            HSF_Early = [HSF_Early; H_Early];
            HSF_Late = [HSF_Late; H_Late];

        elseif ~Lfirst

            H_split_1 = round(HL_idx/2,0);
            H_split_2 = round((num_trials - LH_idx)/2,0);

            H_Early_1 = session_data(1:H_split_1,:);
            H_Late_1 = session_data(H_split_1+1:HL_idx,:);

            H_Early_2 = session_data(LH_idx+1:LH_idx+H_split_2,:);
            H_Late_2 = session_data(LH_idx+H_split_2+1:end,:);

            HSF_Early = [HSF_Early; H_Early_1; H_Early_2];
            HSF_Late = [HSF_Late; H_Late_1; H_Late_2];

            L_split = round((LH_idx - HL_idx)/2,0);

            L_Early = session_data(HL_idx+1:HL_idx+L_split,:);
            L_Late = session_data(HL_idx+L_split+1:LH_idx,:);

            LSF_Early = [LSF_Early; L_Early];
            LSF_Late = [LSF_Late; L_Late];
        end

    elseif num_blocks == 4

        if size(LH_idx,2) > size(HL_idx,2)
            Lfirst = true;
        elseif size(HL_idx,2) > size(LH_idx,2)
            Lfirst = false;
        end

        % Organize/distribute trials based on which condition came first

        if Lfirst

            L_split_1 = round(LH_idx(1)/2,0);
            L_split_2 = round((LH_idx(2) - HL_idx)/2,0);

            L_Early_1 = session_data(1:L_split_1,:);
            L_Late_1 = session_data(L_split_1+1:LH_idx(1),:);

            L_Early_2 = session_data(HL_idx+1:HL_idx+L_split_2,:);
            L_Late_2 = session_data(HL_idx+L_split_2+1:LH_idx(2),:);

            LSF_Early = [LSF_Early; L_Early_1; L_Early_2];
            LSF_Late = [LSF_Late; L_Late_1; L_Late_2];

            H_split_1 = round((HL_idx - LH_idx(1))/2,0);
            H_split_2 = round((num_trials - LH_idx(2))/2,0);

            H_Early_1 = session_data(LH_idx(1)+1:LH_idx(1)+H_split_1,:);
            H_Late_1 = session_data(LH_idx(1)+H_split_1+1:HL_idx,:);

            H_Early_2 = session_data(LH_idx(2)+1:LH_idx(2)+H_split_2,:);
            H_Late_2 = session_data(LH_idx(2)+H_split_2+1:num_trials,:);

            HSF_Early = [HSF_Early; H_Early_1; H_Early_2];
            HSF_Late = [HSF_Late; H_Late_1; H_Late_2];

        elseif ~Lfirst

            H_split_1 = round(HL_idx(1)/2,0);
            H_split_2 = round((HL_idx(2) - LH_idx)/2,0);

            H_Early_1 = session_data(1:H_split_1,:);
            H_Late_1 = session_data(H_split_1+1:HL_idx(1),:);

            H_Early_2 = session_data(LH_idx+1:LH_idx+H_split_2,:);
            H_Late_2 = session_data(LH_idx+H_split_2+1:HL_idx(2),:);

            HSF_Early = [HSF_Early; H_Early_1; H_Early_2];
            HSF_Late = [HSF_Late; H_Late_1; H_Late_2];

            L_split_1 = round((LH_idx - HL_idx(1))/2,0);
            L_split_2 = round((num_trials - HL_idx(2))/2,0);

            L_Early_1 = session_data(HL_idx(1)+1:HL_idx(1)+L_split_1,:);
            L_Late_1 = session_data(HL_idx(1)+L_split_1+1:LH_idx,:);

            L_Early_2 = session_data(HL_idx(2)+1:HL_idx(2)+L_split_2,:);
            L_Late_2 = session_data(HL_idx(2)+L_split_2+1:num_trials,:);

            LSF_Early = [LSF_Early; L_Early_1; L_Early_2];
            LSF_Late = [LSF_Late; L_Late_1; L_Late_2];
        end

    elseif num_blocks == 5

        if LH_idx(1) < HL_idx(1)
            Lfirst = true;
        elseif HL_idx(1) < LH_idx(1)
            Lfirst = false;
        end

        % Organize/distribute trials based on which condition came first

        if Lfirst

            L_split_1 = round(LH_idx(1)/2,0);
            L_split_2 = round((LH_idx(2)-HL_idx(1))/2,0);
            L_split_3 = round((num_trials-HL_idx(2))/2,0);

            L_Early_1 = session_data(1:L_split_1,:);
            L_Late_1 = session_data(L_split_1+1:LH_idx(1),:);

            L_Early_2 = session_data(HL_idx(1)+1:HL_idx(1)+L_split_2,:);
            L_Late_2 = session_data(HL_idx(1)+L_split_2+1:LH_idx(2),:);

            L_Early_3 = session_data(HL_idx(2)+1:HL_idx(2)+L_split_3,:);
            L_Late_3 = session_data(HL_idx(2)+L_split_3+1:num_trials,:);

            LSF_Early = [LSF_Early; L_Early_1; L_Early_2; L_Early_3];
            LSF_Late = [LSF_Late; L_Late_1; L_Late_2; L_Late_3];

            H_split_1 = round((HL_idx(1)-LH_idx(1))/2,0);
            H_split_2 = round((HL_idx(2)-LH_idx(2))/2,0);

            H_Early_1 = session_data(LH_idx(1)+1:LH_idx(1)+H_split_1,:);
            H_Late_1 = session_data(LH_idx(1)+1+H_split_1:HL_idx(1),:);

            H_Early_2 = session_data(LH_idx(2)+1:LH_idx(2)+H_split_2,:);
            H_Late_2 = session_data(LH_idx(2)+H_split_2+1:HL_idx(2),:);

            HSF_Early = [HSF_Early; H_Early_1; H_Early_2];
            HSF_Late = [HSF_Late; H_Late_1; H_Late_2];

        elseif ~Lfirst

            H_split_1 = round(HL_idx(1)/2,0);
            H_split_2 = round((HL_idx(2)-LH_idx(1))/2,0);
            H_split_3 = round((num_trials-LH_idx(2))/2,0);

            H_Early_1 = session_data(1:H_split_1,:);
            H_Late_1 = session_data(H_split_1+1:HL_idx(1),:);

            H_Early_2 = session_data(LH_idx(1)+1:LH_idx(1)+H_split_2,:);
            H_Late_2 = session_data(LH_idx(1)+H_split_2+1:HL_idx(2),:);

            H_Early_3 = session_data(LH_idx(2)+1:LH_idx(2)+H_split_3,:);
            H_Late_3 = session_data(LH_idx(2)+H_split_3+1:num_trials,:);

            HSF_Early = [HSF_Early; H_Early_1; H_Early_2; H_Early_3];
            HSF_Late = [HSF_Late; H_Late_1; H_Late_2; H_Late_3];

            L_split_1 = round((LH_idx(1)-HL_idx(1))/2,0);
            L_split_2 = round((LH_idx(2)-HL_idx(2))/2,0);

            L_Early_1 = session_data(HL_idx(1)+1:HL_idx(1)+L_split_1,:);
            L_Late_1 = session_data(HL_idx(1)+1+L_split_1:LH_idx(1),:);

            L_Early_2 = session_data(HL_idx(2)+1:HL_idx(2)+L_split_2,:);
            L_Late_2 = session_data(HL_idx(2)+L_split_2+1:LH_idx(2),:);

            LSF_Early = [LSF_Early; L_Early_1; L_Early_2];
            LSF_Late = [LSF_Late; L_Late_1; L_Late_2];
        end
    end
    clear LH_idx HL_idx
end

Early_Block = [LSF_Early; HSF_Early];
Late_Block = [LSF_Late; HSF_Late];

%%
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

% Process neural data

% Output spike matrices organization:
    %   1. Low switch frequency switch (PREF)
    %   2. Low switch frequency switch (NULL)
    %   3. High switch frequency switch (PREF)
    %   4. High switch frequency switch (NULL)
    %   5. Low switch frequency non-switch (NULL) % testing epoch
    %   6. Low switch frequency non-switch (PREF)
    %   7. High switch frequency non-switch (NULL)
    %   8. High switch frequency non-switch (PREF)

bin1 = -2350;
bin2 = -2000;

start_bin = find(bins(:,1) == bin1 - 0.5*bin_size,1);
end_bin = find(bins(:,1) == bin2 - 0.5*bin_size,1);

% Process each block half (Early, then Late) and save its LSF-HSF
% difference. Both halves are needed by the Figure 4-Figure Supplement 1B
% plot below, so both must run in the same pass.
blockHalves = {Early_Block, Late_Block};

for bb = 1:numel(blockHalves)
    dat = blockHalves{bb};

    % Include normalization term from session-wise neural analysis:
    [spike_rates_avg, spike_rates_SEM, spike_rates_avg_std_norm_avg, ...
        spike_rates_avg_std_norm_SEM, ...
        ROC_area, cell_selectivity, ...
        normalizationTerm, uniqueUnitNames, ...
        unit_example] = processNeuralData(dat, slide, bin_size, start_time, end_time, [], normalizationTerm_input);

    HSF = mean([mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,3), 2, 'omitnan'), mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,7), 2, 'omitnan')], 2, 'omitnan');
    LSF = mean([mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,1), 2, 'omitnan'), mean(spike_rates_avg_std_norm_avg(:,start_bin:end_bin,5), 2, 'omitnan')], 2, 'omitnan');

    if bb == 1
        combined_diff_E = LSF - HSF; % Early block
    else
        combined_diff_L = LSF - HSF; % Late block
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 4-Figure Supplement 1B %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot difference in neural response (LSF-HSF) for Early vs. Late Blocks
% Just want Monkey Mi for this analysis
% (Ch is excluded upstream, and uniqueUnitNames now reflects only the
% An+Mi units returned by processNeuralData for Early_Block/Late_Block,
% so the Monkey Mi mask must be recomputed for this reduced unit list.)

blockMonkeyIdx = getMonkeyIndices(uniqueUnitNames);
Mi = blockMonkeyIdx.Mi;

% plot controls:
figure(1); clf
hold on; box on;

x_line = -0.2:0.1:2;
y_line = -0.2:0.1:2;

ylim([-0.2,0.6])
xlim([-0.2,0.6])

plot(x_line, y_line, 'k');

x = combined_diff_E(Mi);
y = combined_diff_L(Mi);

plot(x,y,'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

xlabel('LSF-HSF Early in block')
ylabel('LSF-HSF Late in block')

[p, h, stats] = signrank(x,y);
d = computeCohenDCI(y, x, 'paired');
subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])