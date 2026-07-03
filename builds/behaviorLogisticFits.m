
%% behaviorLogisticFits.m
% Fits behavioral data with a time-dependent logistic psychometric function
% to model choice behavior as a function of viewing duration.
%
% Model structure:
%   P(switch) = logistic(β₀ + β₁*dir_bias + β₂*signed_time)
%   where signed_time is negative for non-switch trials, positive for switch trials
%
% Fits are performed separately for:
%   - Each session
%   - Each switch frequency
%
% Outputs:
%   fits - [numSessions x 4 parameters x 2 switch frequency] Logistic regression coefficients
%          β₀: Switch/stay bias
%          β₁: Right/left directional bias  
%          β₂: Sensitivity to viewing duration (time slope)
%          β₃: Lapse rate
%   R_sq - [numSessions x 2 hazards] Tjur's pseudo-R² (model goodness of fit)
%
% The script also generates and optionally saves psychometric curve plots 
% showing fitted functions overlaid on smoothed behavioral data for each 
% session and hazard rate condition.
%
% Note: Fitting is computationally intensive. Pre-computed fits are saved
% in data/behaviorFits/ as LogisticFits.mat and LogisticFits_Rsq.mat
%
% Required functions: logistFitDotsrev.m, logistValDotsrev.m, logistErrDotsrev.m (in logisticUtilities/ folder)
%%
cfg = projectDefaults();

%%

cd(cfg.paths.data)
load('mergedTable_proc_core.mat') % behavioral-only: skips Unit_1 and pupil traces (see buildMergedTableTiers.m)

%%

% Subset data appropriately
% "B" = Behavioral analysis

[mergedTableSub] = createDatSubset(mergedTable_proc, 'B');
dat = mergedTableSub;

cd(cfg.paths.fits)

%%
% Get unique session names

uniqueSessionNames = unique(dat.ses_ID);
numSessions = length(uniqueSessionNames);
uniqueHazards = cfg.hazard.codes;

% Save the fits
fits = nans(numSessions, 4, 2);
LL = nans(numSessions, 1, 2);

% Save the R_sq values
R_sq = nans(numSessions, 1, 2);

% Save coherence per session
coherences = nans(numSessions, 1);

% For plotting
plot_it = true;
xax = (-1200:1200)';
show_fit_data = repmat(cat(2, ones(length(xax),2), xax),[1 1 2]);
show_fit_data(:,2,2) = 0;
co = cfg.colors.pair;

ses_save = false;
filePath = fullfile(cfg.paths.fits, 'plots'); % only used if ses_save = true
if ses_save && ~isfolder(filePath)
    mkdir(filePath);
end

%%
% Loop through the sessions
for ss = 1:numSessions

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

    % Set up data matrix (bias, right/left bias, signed time, choice) --
    % see buildBehaviorDataToFit.m for column definitions.
    [data_to_fit, Lpre] = buildBehaviorDataToFit(session_data);

    % Save coherence
    coherences(ss) = unique(session_data(:,2));
    
    % Fit/plot separately for each switch-frequency condition
    for hh = 1:2
        Lhazard = session_data(:,1) == uniqueHazards(hh);
        
        % Get the fit 
        if sum(Lhazard) > 20
           
            [fits(ss,:,hh), LL(ss,:,hh)] = logistFitDotsrev(data_to_fit(Lhazard,:));

            % Compute Tjur's D (Tjur's coefficient of disctrimination -- R^2)
                 % Compares average fitted probability of the two response outcomes (0 and 1)
                 % Calculate mean of predicted probabilities of each outcome and take the difference between the 2 means
                 % 0 <= R^2 >= 1 where no disriminating power = 0 and perfect discriminating power = 1

             fit_y = logistValDotsrev(fits(ss,:,hh)', data_to_fit(Lhazard,1:end-1));
             R_sq(ss,hh) = computeTjursR2(fit_y, data_to_fit(Lhazard,4));

            % Plot it separately for pre-final=L and pre=final=R
            if plot_it
                for xx = 1:2
                    figure(1)
                    subplot(1,2,xx); hold on
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
                end
            end
        end

        hold on 
    end

    % Save plots for each session

    if ses_save

    plotName = [char(uniqueSessionNames(ss)), '_LogisticFit'];
    saveas(gcf, fullfile(filePath, plotName), 'svg');

    end

    clf
end

%% Saved output

% fit_mat_name = 'LogisticFits.mat';
% save(fit_mat_name, 'fits')
% 
% R_sq_mat_name = 'LogisticFits_Rsq.mat';
% save(R_sq_mat_name, 'R_sq')
