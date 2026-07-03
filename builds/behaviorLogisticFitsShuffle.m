%% behaviorLogisticFitsShuffle.m
% Shuffle control analysis for behavioral logistic fits. Tests whether the 
% time-dependent effects in the main behavioral model could arise by chance.
%
% Procedure:
%   For each session, randomly shuffles the viewing durations and switch/non-switch trial types.
%   Refits the logistic model to this shuffled data. Process is repeated 
%   100 times per session to generate a null distribution.
%
% Purpose:
%   Establishes a null distribution for comparing against real behavioral fits.
%
% Outputs:
%   fits_shuffle - [numSessions*100 x 4 parameters x 2 switch-frequency conditions] 
%                  Logistic coefficients from shuffled data
%   R_sq_shuffle - [numSessions*100 x 2 switch-frequency conditions] 
%                  Tjur's pseudo-R² for shuffled fits
%
% Note: This analysis is extremely time-intensive (hours of runtime for
% 100 iterations across all sessions). Pre-computed results are saved in
% data/behaviorFits/ as LogisticFits_Shuffle.mat and LogisticFits_Rsq_Shuffle.mat
%
% Required functions: logistFitDotsrev.m, logistValDotsrev.m, logistErrDotsrev.m

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

% Set number of iterations
numIterations = 100;

% Save the fits
fits_shuffle = nans(numSessions*numIterations, 4, 2);
LL_shuffle = nans(numSessions*numIterations, 1, 2);

% Save the R_sq values
R_sq_shuffle = nans(numSessions*numIterations, 1, 2);

% Save coherence per session
coherences = nans(numSessions*numIterations, 1);

% Loop through the sessions
% NOTE: session outer / iteration inner (rather than the reverse) so the
% per-session table subsetting and data_to_fit setup below -- which don't
% depend on the shuffle -- run once per session instead of once per
% (iteration, session) pair, a 100x reduction in that portion of the work.
% This changes which row of fits_shuffle/R_sq_shuffle/coherences a given
% (session, iteration) lands in relative to the original ii-outer/ss-inner
% ordering, but since this is a shuffled null distribution the row order
% carries no meaning -- only the aggregate distribution does.

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
    % see buildBehaviorDataToFit.m for column definitions. Lpre isn't
    % used in this script (no plotting), so it's not captured here.
    data_to_fit_base = buildBehaviorDataToFit(session_data);

    sessionCoherence = unique(session_data(:,2));

    for ii = 1:numIterations

        vv = (ss - 1) * numIterations + ii;

        print_test = ['Running model #', num2str(vv)];
        disp(print_test)

        % Shuffle dur_final column
        data_to_fit = data_to_fit_base;
        idx = randperm(size(data_to_fit,1));
        data_to_fit(:,3) = data_to_fit_base(idx,3);

        % Save coherence
        coherences(vv) = sessionCoherence;

        % Fit separately for each hazard
        for hh = 1:2
            Lhazard = session_data(:,1) == uniqueHazards(hh);

            % Get the fit
            if sum(Lhazard) > 20

                 [fits_shuffle(vv,:,hh), LL_shuffle(vv,:,hh)] = logistFitDotsrev(data_to_fit(Lhazard,:));

                 % Compute Tjur's D (Tjur's coefficient of disctrimination -- R^2)
                    % Compares average fitted probability of the two response outcomes (0 and 1)
                    % Calculate mean of predicted probabilities of each outcome and take the difference between the 2 means
                    % 0 <= R^2 >= 1 where no disriminating power = 0 and perfect discriminating power = 1

                 fit_y = logistValDotsrev(fits_shuffle(vv,:,hh)', data_to_fit(Lhazard,1:end-1));
                 R_sq_shuffle(vv,hh) = computeTjursR2(fit_y, data_to_fit(Lhazard,4));
            end
        end
    end
end

%% Saved output

% fit_mat_name = 'LogisticFits_Shuffle.mat';
% save(fit_mat_name, 'fits_shuffle')
% 
% R_sq_mat_name = 'LogisticFits_Rsq_Shuffle.mat';
% save(R_sq_mat_name, 'R_sq_shuffle')
