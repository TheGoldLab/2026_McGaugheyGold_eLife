%% Fig6_FigSupp1.m
% Analyzes baseline pupil diameter and its relationship to task engagement
% measured by time to fixation.
%
% Analysis approach:
%   - Extracts baseline pupil using processPupilBaseline.m
%   - Computes trial-wise Spearman correlation between time-to-fixation and
%     baseline pupil for each session 
%   - Compares baseline pupil diameter between LSF and HSF conditions
%   - Analyzes Monkeys An and Mi separately (Ch excluded from pupil analysis)
%
% Outputs:
%   Figure 6-Figure Supplement 1A: Example sessions showing time course of:
%     - Time to fixation (moving average)
%     - Baseline pupil (moving average)
%     - Switch-frequency condition markers (blue=LSF, yellow=HSF)
%     Sessions: Monkey An (session 30), Monkey Mi (session 151)
%
%   Figure 6-Figure Supplement 1B: Distribution of session-wise correlations
%     between time-to-fixation and baseline pupil, separated by monkey.
%     Red dashed lines show mean correlation. One-sample t-tests test
%     whether correlations significantly differ from zero.
%
%   Figure 6-Figure Supplement 1C: Scatter plots comparing LSF vs HSF baseline
%     pupil diameter for each session. Statistics: Wilcoxon signed-rank
%     test (paired) and Cohen's d effect size.
%
% Required data files:
%   - mergedTable_proc.mat
%
% Required functions:
%   - processPupilBaseline.m (extracts baseline pupil measurements)
%   - computeCohenDCI.m (calculates effect size)

%% Load data
cfg = projectDefaults();
cd(cfg.paths.data)

load('mergedTable_proc.mat')

% Subset data appropriately
% "BP" = Sessions with behavior and pupil data
    % Two sessions from Monkey An had no pupil data recorded
[mergedTableSub] = createDatSubset(mergedTable_proc, 'BP');
dat = mergedTableSub;

% Get unique session names
uniqueSessionNames = unique(dat.ses_ID);
numSessions = length(uniqueSessionNames);

% Process pupil data, extracting only baseline pupil diameter
baseline_pupil = processPupilBaseline(dat);

%%
% Loop through sessions, calculating correlation between baseline pupil diameter and time to fixation
% Ignoring Monkey Ch, who was excluded from pupil analysis

for ss = 1:numSessions

    Lses = strcmp(dat.ses_ID, uniqueSessionNames{ss});

    time_to_fix = table2array(dat(Lses, 13));
    HR = table2array(dat(Lses, 2));
    baseline_pupil_ses = baseline_pupil(Lses,:);

    % Save variables and correlations
    TF_baseline(ss,1) = corr(time_to_fix,baseline_pupil_ses, 'type', 'Spearman');

    % Save condition-averaged baseline pupil diameter
    LSF_baseline(ss,1) = mean(baseline_pupil_ses(HR == 2), 'omitnan');
    HSF_baseline(ss,1) = mean(baseline_pupil_ses(HR == 6), 'omitnan');

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6-Figure Supplement 1A %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot example correlations for two sessions, one per monkey: top row =
% time to fixation, bottom row = baseline pupil (matches the manuscript
% layout: An left, Mi right).
    % Monkey An: session 30
    % Monkey Mi: session 151

exampleSessions = [30 151];
exampleNames = {'Monkey An', 'Monkey Mi'};

figure(71); clf

for mm = 1:numel(exampleSessions)
    ss = exampleSessions(mm);

    Lses = strcmp(dat.ses_ID, uniqueSessionNames{ss});

    time_to_fix = table2array(dat(Lses, 13));
    HR = table2array(dat(Lses, 2));
    baseline_pupil_ses = baseline_pupil(Lses,:);
    num_trials = sum(Lses);

    % Create condition vectors
    LSF_on = HR == 2;
    HSF_on = ~LSF_on;

    % Plot time to fixation as a function of trial number
    subplot(2,2,mm); box on; hold on;

    plot(1:num_trials, movmean(time_to_fix,50))

    yl = ylim;
    LSF_line = repmat(yl(2),num_trials,1);
    HSF_line = repmat(yl(2),num_trials,1);

    LSF_line(~LSF_on) = nan;
    HSF_line(~HSF_on) = nan;

    scatter(1:num_trials, LSF_line, 'MarkerEdgeColor', cfg.colors.LSF, 'MarkerFaceColor', cfg.colors.LSF)
    scatter(1:num_trials, HSF_line, 'MarkerEdgeColor', cfg.colors.HSF, 'MarkerFaceColor', cfg.colors.HSF)

    title(sprintf('%s: Time to fixation', exampleNames{mm}))

    % Plot baseline pupil (z-scored) as a function of trial number
    subplot(2,2,2+mm); box on; hold on;

    plot(1:num_trials, movmean(baseline_pupil_ses,50))

    yl = ylim;
    LSF_line = repmat(yl(2),num_trials,1);
    HSF_line = repmat(yl(2),num_trials,1);

    LSF_line(~LSF_on) = nan;
    HSF_line(~HSF_on) = nan;

    scatter(1:num_trials, LSF_line, 'MarkerEdgeColor', cfg.colors.LSF, 'MarkerFaceColor', cfg.colors.LSF)
    scatter(1:num_trials, HSF_line, 'MarkerEdgeColor', cfg.colors.HSF, 'MarkerFaceColor', cfg.colors.HSF)

    title(sprintf('%s: Baseline pupil (z-scored)', exampleNames{mm}))
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6-Figure Supplement 1B %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distribution of trial-wise correlation values (time to fixation and baseline pupil)
% Plotted separately for Monkey An and Monkey Mi

monkeyIdx = getMonkeyIndices(uniqueSessionNames);
An = monkeyIdx.An;
Mi = monkeyIdx.Mi;

monkeyOrder = {'An', 'Mi'};
monkeyMarkers = struct('An', 'o', 'Mi', 'square');
monkeyNumIdx = struct('An', An, 'Mi', Mi);

% Correlations
edges = linspace(0, 0.6, 7);

figure(41); clf; hold on; box on;

% ttestResult.(monkeyName) keeps each monkey's [h,p,ci,stats] separate --
% both calls used to write into the same H/P/CI/STATS variables, so the An
% result was immediately overwritten by the Mi result.
ttestResult = struct();

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    idx = monkeyNumIdx.(monkeyName);

    histogram(TF_baseline(idx), edges, 'Normalization','probability')
    avg = mean(TF_baseline(idx), 'omitnan');
    xline(avg, '--r','LineWidth',3)

    [h, p, ci, stats] = ttest(TF_baseline(idx));
    ttestResult.(monkeyName).h = h;
    ttestResult.(monkeyName).p = p;
    ttestResult.(monkeyName).ci = ci;
    ttestResult.(monkeyName).stats = stats;
end

title('Fixation & baseline pupil')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 6-Figure Supplement 1C %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Comparison of baseline pupil diameter across LSF and HSF conditions
% Plotted separately for Monkey An and Monkey Mi

x_line = -1:.1:2;
y_line = -1:.1:2;

figure(109); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    idx = monkeyNumIdx.(monkeyName);

    subplot(1,2,mm); hold on; box on;

    x = HSF_baseline(idx);
    y = LSF_baseline(idx);

    plot(x_line, y_line, 'k');
    plot(x,y,'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    xlabel('HSF avg. baseline pupil')
    ylabel('LSF avg. baseline pupil')

    [p, h, stats] = signrank(x,y);
    d = computeCohenDCI(y, x, 'paired');

    title('Baseline pupil diameter')
    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

    ylim([0,2])
    xlim([0,2])
end