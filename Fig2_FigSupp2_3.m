%% Fig2_FigSupp2_3.m
% Analyzes behavioral sensitivity as a function of viewing duration across switch-frequency
% conditions using time-dependent logistic psychometric functions.
%
% Analysis:
%   - Loads pre-computed logistic fits from behaviorLogisticFits.m
%   - Extracts β₂ coefficient (time-dependent slope) from fits for LSF vs HSF conditions
%   - Compares slopes across conditions and monkeys
%
% Main Text Figures:
%
%   Figure 2B: Example psychometric functions for one representative session
%     per monkey. Shows smoothed choice data (dotted lines) and fitted
%     logistic curves (solid lines) for LSF (blue) and HSF (yellow).
%     Separate panels for prefinal direction right vs. left. Demonstrates
%     individual session fits with clear time-dependent choice patterns.
%     Example sessions: An (ss=28), Ch (ss=58), Mi (ss=142)
%
%   Figure 2C: Comparison of psychometric slopes (β₂) across conditions.
%     Scatter plot with LSF slope vs HSF slope for each session. Points
%     above unity line indicate greater sensitivity at LSF; below indicates
%     greater sensitivity at HSF. Different markers by monkey:
%       - Monkey An: circles
%       - Monkey Ch: diamonds
%       - Monkey Mi: squares
%     Statistics: Wilcoxon signed-rank test (all sessions) + Cohen's d
%
%   Figure 2D: Behavioral performance on switch trials as a function of
%     test-stimulus duration (4 bins with midpoints: 162.5, 300, 487.5, 900ms).
%     LSF (blue) and HSF (yellow), mean ± SEM across sessions.
%
% Figure Supplements:
%
%   Figure 2-Figure Supplement 2: Same as Figure 2C but separated by individual
%     monkeys in three panels. Shows that context-dependent sensitivity
%     adjustments are present in each subject individually.
%
%   Figure 2-Figure Supplement 3: Behavioral performance on switch trials as a
%     function of viewing duration. Shows that behavioral dynamics are
%     consistent with leaky evidence accumulation.
%
% Required data files:
%   - mergedTable_proc.mat (behavioral data)
%   - LogisticFits.mat (pre-computed regression coefficients)
%
% Required functions:
%   - logistValDotsrev.m (evaluate logistic function)
%   - computeCohenDCI.m (calculate effect size)
%   - behaviorBinnedPerformance.m (compute binned performance)

%% Load data
cfg = projectDefaults();

cd(cfg.paths.data)
load('mergedTable_proc_core.mat') % behavioral-only: skips Unit_1 and pupil traces (see buildMergedTableTiers.m)

cd(cfg.paths.fits)
load('LogisticFits.mat')

% Subset data appropriately
% "B" = Behavioral analysis
[mergedTableSub] = createDatSubset(mergedTable_proc, 'B');
dat = mergedTableSub;

% Get unique session names
uniqueSessionNames = unique(dat.ses_ID);
uniqueHazards = cfg.hazard.codes;

% Assign monkey indices based on session names
monkeyIdx = getMonkeyIndices(uniqueSessionNames);
An = monkeyIdx.An;
Ch = monkeyIdx.Ch;
Mi = monkeyIdx.Mi;

% For plotting
xax = (-1200:1200)';
show_fit_data = repmat(cat(2, ones(length(xax),2), xax),[1 1 2]);
show_fit_data(:,2,2) = 0;
co = cfg.colors.pair;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 2B %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example behavior figures
% One representative session (ss) from each monkey
    % An: ss = 28
    % Ch: ss = 58
    % Mi: ss = 142

% Change session number:
ss = 142;

disp(uniqueSessionNames{ss})

% Get the data to fit:
%   1. hazard
%   2. test coherence
%   3. prefinal direction
%   4. final direction
%   5. test duration
%   6. choice (1=right, 2=left)

Lses = strcmp(dat.ses_ID, uniqueSessionNames{ss}) & ...
    dat.coh_final >= 50;

session_data = table2array(dat(Lses, [2 4 7 6 9 10]));

% For data_to_fit
LdirSwitch = session_data(:,3) ~= session_data(:,4);
LpreDirIsRight = session_data(:,3) < 90 | session_data(:,3) >= 270;
Lpre = [LpreDirIsRight ~LpreDirIsRight];
LchoseRight = session_data(:,6) == 1;

% Set up data matrix
%  1. Switch/stay bias ... column of ones
%  2. Right/left bias ... 0 when prefinal dir=L, 1 when prefinal dir=R
%  3. Signed time ... neg=non-switch, pos=switch
%  4. Choice ... 0/1 = did not/did choose prefinal dir

data_to_fit = ones(sum(Lses), 4);
data_to_fit(~LpreDirIsRight, 2) = 0;
data_to_fit(:,3) = session_data(:,5);
data_to_fit(~LdirSwitch,3) = -data_to_fit(~LdirSwitch,3);
data_to_fit(LpreDirIsRight & LchoseRight,4) = 0; % prefinal dir = R and chose R
data_to_fit(LpreDirIsRight &~LchoseRight,4) = 1; % prefinal dir = R and chose L
data_to_fit(~LpreDirIsRight& LchoseRight,4) = 1; % prefinal dir = L and chose R
data_to_fit(~LpreDirIsRight&~LchoseRight,4) = 0; % prefinal dir = L and chose L

% Fit/plot separately for each hazard
figure(3); clf
for hh = 1:2
    Lhazard = session_data(:,1) == uniqueHazards(hh);

    % Get the fit
    if sum(Lhazard) > 20

         fit_y = logistValDotsrev(fits(ss,:,hh)', data_to_fit(Lhazard,1:end-1));
         R_sq(ss,hh) = abs(mean(fit_y(data_to_fit(Lhazard,4) == 1,:)) -  mean(fit_y(data_to_fit(Lhazard,4) == 0,:)));

        % Plot it separately for predir=L and predir=R
        for xx = 1:2
            subplot(1,2,xx); hold on; box on
            plot([0 0], [0 1], 'k--');
            plot(xax([1 end]), [0.5 0.5], 'k--', 'MarkerSize', 15);

            % Show smoothed data
            [sorted_time_axis, I] = sort(data_to_fit(Lhazard&Lpre(:,xx), 3));
            choice_data = data_to_fit(Lhazard&Lpre(:,xx), 4);
            sorted_choice_data = choice_data(I);
            plot(sorted_time_axis, nanrunmean(sorted_choice_data,5), ':', 'Color', co{hh});

            if xx==1
                title(sprintf('Prefinal dir=R, h=%d', uniqueHazards(hh)))
                ylabel('Fraction choose Left')
            else
                title(sprintf('Prefinal dir=L, h=%d', uniqueHazards(hh)))
                ylabel('Fraction choose Right')
            end

            xlabel('Signed time (-non_switch, +switch)');

            % Show fits
            ys = logistValDotsrev(fits(ss,:,hh)', show_fit_data(:,:,xx));
            plot(xax, ys, '-', 'Color', co{hh}, 'LineWidth', 2);
            axis([xax(1) xax(end) 0 1]);

            % Calculate y-intercept (bias)
            bias(ss,xx,hh) = 0.5 - ys(1201);
        end
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 2C %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare sensitivity/psychometric slope across monkeys

x_line = -0.05:0.01:0.05;
y_line = -0.05:0.01:0.05;

figure(4); clf; box on; hold on;

plot(x_line, y_line, 'k')

plot(fits(An,3,2), fits(An,3,1), 'k', 'Marker', 'o', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')
plot(fits(Ch,3,2), fits(Ch,3,1), 'k', 'Marker', 'diamond', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')
plot(fits(Mi,3,2), fits(Mi,3,1), 'k', 'Marker', 'square', 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

ylabel('Low switch frequency slope')
xlabel('High switch frequency slope')

LSF_slope = fits(:,3,1);
HSF_slope = fits(:,3,2);

[p, h, stats] = signrank(LSF_slope,HSF_slope);
d = computeCohenDCI(LSF_slope, HSF_slope, 'paired');

subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

ylim([0 0.05])
xlim([0 0.05])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 2D %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get binned behavioral accuracy
% Returns vector of 4 bins:
    % Bin 1: 100-225 ms
    % Bin 2: 225-375 ms
    % Bin 3: 375-600 ms
    % Bin 4: 600-1200 ms

[LSF_switch_binned_performance, HSF_switch_binned_performance, binned_behavior_diff] = ...
    behaviorBinnedPerformance(dat);

% Calculate average and SEM

LSF_switch_binned_performance_avg = nanmean(LSF_switch_binned_performance);
HSF_switch_binned_performance_avg = nanmean(HSF_switch_binned_performance);

LSF_switch_binned_performance_avg_sem = nanse(LSF_switch_binned_performance);
HSF_switch_binned_performance_avg_sem = nanse(HSF_switch_binned_performance);

figure(5); clf; box on; hold on;

errorbar(cfg.bins.viewDuration.midpoints, LSF_switch_binned_performance_avg, LSF_switch_binned_performance_avg_sem, 'Color', co{1},'LineWidth', 3)
errorbar(cfg.bins.viewDuration.midpoints, HSF_switch_binned_performance_avg, HSF_switch_binned_performance_avg_sem,  'Color', co{2}, 'LineWidth', 3)

ylabel('Fraction correct')
xlabel('Binned viewing duration (ms)')

%% Relevant Figure Supplements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 2-Figure Supplement 2 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compare sensitivity/psychometric slope for individual animals
% Panel order (An, Mi, Ch) matches the manuscript, not
% monkeyIdx.names' alphabetical (An, Ch, Mi) order.

x_line = -0.05:0.01:0.05;
y_line = -0.05:0.01:0.05;

monkeyOrder = {'An', 'Mi', 'Ch'};
monkeyMarkers = struct('An', 'o', 'Mi', 'square', 'Ch', 'diamond');

figure(14); clf

for mm = 1:numel(monkeyOrder)
    monkeyName = monkeyOrder{mm};
    isMonkey = monkeyIdx.(monkeyName);

    subplot(1,3,mm); box on; hold on;

    plot(x_line, y_line, 'k')

    LSF_slope = fits(isMonkey,3,1);
    HSF_slope = fits(isMonkey,3,2);

    plot(HSF_slope, LSF_slope, 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')

    [p, h, stats] = signrank(LSF_slope, HSF_slope);
    d = computeCohenDCI(LSF_slope, HSF_slope, 'paired');

    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

    ylabel('Low switch frequency slope')
    xlabel('High switch frequency slope')

    ylim([0 0.05])
    xlim([0 0.05])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Figure 2-Figure Supplement 3 %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Behavioral accuracy (LSF vs. HSF) for shorter- (bin 1) and longer-
% duration (bin 3) switch trials, separated by monkey.

x_line = -0.05:0.01:1;
y_line = -0.05:0.01:1;

binsToShow = [1 3];

figure(41); clf

for bb = 1:numel(binsToShow)
    thisBin = binsToShow(bb);

    subplot(1,2,bb); box on; hold on;

    plot(x_line, y_line, 'k')

    for mm = 1:numel(monkeyOrder)
        monkeyName = monkeyOrder{mm};
        isMonkey = monkeyIdx.(monkeyName);

        LSF_binned = LSF_switch_binned_performance(isMonkey, thisBin);
        HSF_binned = HSF_switch_binned_performance(isMonkey, thisBin);

        plot(HSF_binned, LSF_binned, 'k', 'Marker', monkeyMarkers.(monkeyName), 'MarkerSize',10, 'LineStyle', 'none', 'MarkerFaceColor','white')
    end

    LSF_binned = LSF_switch_binned_performance(:,thisBin);
    HSF_binned = HSF_switch_binned_performance(:,thisBin);

    [p, h, stats] = signrank(LSF_binned, HSF_binned);
    d = computeCohenDCI(LSF_binned, HSF_binned, 'paired');

    subtitle(['p = ', num2str(p), ' & d = ', num2str(d)])

    ylabel(sprintf('LSF behavioral accuracy (bin %d)', thisBin))
    xlabel(sprintf('HSF behavioral accuracy (bin %d)', thisBin))

    ylim([0 1])
    xlim([0 1])
end
