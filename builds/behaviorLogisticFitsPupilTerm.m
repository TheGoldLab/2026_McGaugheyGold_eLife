%% behaviorLogisticFitsPupilTerm.m
% Fits behavioral data with a time-dependent logistic psychometric function
% that includes a trial-wise pupil diameter term.
%
% Model structure:
%   P(switch) = logistic(β₀ + β₁*dir_bias + β₂*pupil×time + β₃*signed_time)
%   where pupil term is the evoked pupil response averaged 500ms before 
%   test stimulus onset (baseline-corrected and z-scored).
%
% Control condition (control = true):
%   Fits standard behavior-only model:
%   P(switch) = logistic(β₀ + β₁*dir_bias + β₂*signed_time)
%
% Pupil processing:
%   - Low-pass filters pupil at 5 Hz (Butterworth)
%   - Detects and removes outliers (>2 SD from filtered signal)
%   - Z-scores pupil diameter after subtracting running average
%   - Baseline-corrects using linear regression
%   - Extracts evoked pupil 500ms before test onset (-500 to 0ms)
%   - Clips extreme values to ±1 SD to reduce outlier influence
%
% Fits are performed separately for:
%   - Each recording session
%   - Each switch-frequency condition (LSF=2, HSF=6)
%
% Outputs:
%   fits - [numUnits x 5 parameters x 2 frequencies] Pupil model coefficients
%          β₀: Switch/stay bias
%          β₁: Right/left directional bias
%          β₂: Pupil×time interaction
%          β₃: Sensitivity to viewing duration
%          β₄: Lapse rate
%   fits_control - [numUnits x 4 parameters x 2 frequencies] Control model
%   R_sq - [numUnits x 2 frequencies] Tjur's pseudo-R² for pupil model
%   R_sq_control - [numUnits x 2 frequencies] Tjur's pseudo-R² for control
%
% Note: Fitting is computationally intensive. Pre-computed fits are saved in
% data/behaviorFits/ as LogisticFits_PupilTerm.mat and LogisticFits_PupilTerm_control.mat
%
% Required functions: 
%   logistFitDotsrevNP.m, logistValDotsrevNP.m, logistErrDotsrev.m (pupil model)
%   logistFitDotsrev.m, logistValDotsrev.m, logistErrDotsrev.m  (control model)

%%
% Set paths

cfg = projectDefaults();
cd(cfg.paths.data)

% Load data:

load('mergedTable_proc.mat')

% Subset data appropriately

[mergedTableSub] = createDatSubset(mergedTable_proc, 'NP');
dat = mergedTableSub;

clear mergedTable_proc

%% Choose condition to run

control = false; % True will run the behavior-only logistic with no neural term
                 % False will run a logistic that includes a neural term
%%
% Get unique session (unit) names:

uniqueUnitNames = unique(dat.ses_ID);
numUnits = length(uniqueUnitNames);

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
axlen = length(xax);

tax = start_time:end_time;
taxlen = length(tax);

% Dot switch timing

dot_start_times = {(-2400:-1200:0), (-2400:400:0)};
Ldots_P = false(axlen, 3);
for hh = 1:2
    dir_val = true;
    for tt = 1:length(dot_start_times{hh})-1
        Ldots_P(xax>=dot_start_times{hh}(tt)&xax<dot_start_times{hh}(tt+1), hh) = dir_val;
        dir_val = ~dir_val;
    end
end

Ldots_P(xax >= end_time, 3) = true;

% Filter parameters, pupil smoothing (butterworth low-pass):

srate = 1000;           % EyeLink sample rate
lpass = 5;              % low pass cutoff - 4-10 is a good range
freq = (lpass/srate)*2; % frequency
n_order = 1;            % order of butterworth filter
[butter_pb,butter_pa] = butter(n_order,freq,'low');

% Filter parameters, pupil slope

FILTER_HW = 75;
filter_c  = (3 / (2*FILTER_HW^3+3*FILTER_HW^2+FILTER_HW))*(FILTER_HW:-1:-FILTER_HW);

%%
% Save the fits

fits = nans(numUnits, 5, 2);
fits_control = nans(numUnits, 4, 2);

% Save the R_sq values

R_sq = nans(numUnits, 1, 2);
R_sq_control = nans(numUnits, 1, 2);

% Plotting controls:

colors = cfg.colors.pair;
line_styles = {'-', ':'};

%%
% Loop through the units:
for uu = 1:numUnits
    
    % Show the name

    disp(uniqueUnitNames{uu})

    % Monkey Ch is excluded from this per-unit pupil-term/control fit
    % (also excluded in behaviorLogisticFitsNeuralTerm.m and
    % Fig4_FigSupp1B.m's block-split analysis); reason not documented in
    % the original code. Ch is included in the plain behavioral scripts.
    if uniqueUnitNames{uu}(1:2) == "Ch"
        continue
    end

    % Get the data from this session)

    Lses = strcmp(dat.ses_ID, uniqueUnitNames{uu}) & dat.coh_final >= 50; 

    if sum(Lses) == 0
        continue
    end    

    session_data_tbl = dat(Lses, 1:end);
    session_data = table2array(dat(Lses, [2 4 5 9 12 27 13 19 18 27 28]));

    num_trials = size(session_data,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% PUPIL %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make matrix of pupil/x/y values from 100 ms before dots on to dots off

    pupil_data = nans(num_trials, taxlen);
    eyeX_data = nans(num_trials, taxlen);
    eyeY_data = nans(num_trials, taxlen);
    dots_on_times = dat(Lses, 31).dots_on;
    Lgood = false(num_trials, taxlen); % will pad end with final values
    Lkeep = true(num_trials,1);

    clear tt_array 

    for tt = 1:num_trials
        start_index = max(1, round(dots_on_times{tt}(1) + start_time - session_data(tt,10)));
        end_index = min(start_index + taxlen - 1, round(session_data(tt,9) - session_data(tt,10)));
        from_indices = start_index:end_index;
        to_indices = 1:length(from_indices);
        pupil_data(tt,to_indices) = dat(Lses, 24).pupil_diam{tt}(from_indices);
        Lgood(tt,to_indices) = true;
        eyeX_data(tt,to_indices) = dat(Lses, 25).pupil_horiz{tt}(from_indices);
        eyeY_data(tt,to_indices) = dat(Lses, 26).(' pupil_vert'){tt}(from_indices);

        tt_array(tt,1) = tt;
    end

    % interpolate any missing values (note the double transpose to get
    % interpt1 to work properly)

    pupil_data = pupil_data';
    Lbad = ~isfinite(pupil_data);
    pupil_data(Lbad) = interp1(find(~Lbad), pupil_data(~Lbad), find(Lbad), 'linear', 'extrap');

    % Low-pass filter pupil

    pupil_data = pupil_data';
    pupil_filt = pupil_data;
    pupil_filt(~isfinite(pupil_filt)) = mean(pupil_filt(:), 'omitnan');
    pupil_filt1 = filtfilt(butter_pb,butter_pa,pupil_filt')';

    % Use residuals to find bad values, remove, and interpolate

    pupil_resid = (pupil_data - pupil_filt1)';
    pupil_resid_mean = mean(pupil_resid(:), 'omitnan');
    pupil_resid_std = std(pupil_resid(:), 'omitnan');
    Lbad = (pupil_resid > pupil_resid_mean+(2*pupil_resid_std)) | ...
        (pupil_resid < pupil_resid_mean-(2*pupil_resid_std));
    pupil_clean = pupil_data';
    pupil_clean(Lbad) = nan;
    pupil_clean(Lbad) = interp1(find(~Lbad), pupil_clean(~Lbad), find(Lbad), 'linear', 'extrap');
    pupil_clean = pupil_clean';

    % Re-filter cleaned data

    pupil_filt2 = filtfilt(butter_pb,butter_pa,pupil_clean')';

    % Subtract average within/across trials then z-score

    pupil_data = pupil_filt2;
    pupil_data(~Lgood) = nan;

    % Subtract across trial running-average timecourse

    pupil_data = pupil_data - repmat(smooth(mean(pupil_data, 2, 'omitnan'), 50), 1, size(pupil_data, 2));

    % Subtract average timecourse

    pupil_data = pupil_data - repmat(mean(pupil_data, 'omitnan'), num_trials, 1);

    % Z score:

    pupil_data = (pupil_data - mean(pupil_data(:), 'omitnan'))./ ...
        std(pupil_data(:), 'omitnan');

    % Baseline subtract
    % Take residuals from linear fit with baseline
   
    baselines = mean(pupil_data(:, 100:200), 2, 'omitnan');

    % Vectorized closed-form single-predictor regression through the origin,
    % computed for all bins at once instead of calling regress() once per
    % bin (3801 calls/unit). Since baselines is a single column with no
    % intercept term, regress(y(Lgood), baselines(Lgood)) per bin reduces to
    % b = sum(x.*y)/sum(x.^2), residual = y - x*b -- algebraically identical,
    % validated to match regress() to machine precision (~1e-16).
    Lgood_mat = isfinite(pupil_data);
    b_bin = sum(baselines .* pupil_data .* Lgood_mat, 1, 'omitnan') ./ ...
        sum((baselines.^2) .* Lgood_mat, 1, 'omitnan');
    resid_all = pupil_data - baselines .* b_bin;
    pupil_baseline_resid = nans(size(pupil_data));
    pupil_baseline_resid(Lgood_mat) = resid_all(Lgood_mat);

     pupil_data_evoked = pupil_baseline_resid;

    % Find average pupil diameter during 500 ms pre-testing stimulus

    bin1 = -500;
    bin2 = 0;
    
    Lbin = tax >= bin1 & tax < bin2;

    pupil_evoked_bin = mean(pupil_data_evoked(:, Lbin), 2, 'omitnan');

    % Adjust more extreme values

    L_max = pupil_evoked_bin > 1;
    L_min = pupil_evoked_bin < -1;

    pupil_evoked_bin(L_max) = 1;
    pupil_evoked_bin(L_min) = -1;

    sum_max = sum(L_max);
    sum_min = sum(L_min);
    sum_adj = sum_max + sum_min;

    adjustment_txt = ['Adjusted ', num2str(sum_adj), ' trial values'];
    disp(adjustment_txt)

    % Append to session_data_tbl

    session_data_tbl = addvars(session_data_tbl,pupil_evoked_bin,'After',"Unit_1");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% BEHAVIOR %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     % Get the data to fit:
    %   1. hazard
    %   2. test coherence
    %   3. prefinal direction
    %   4. final direction
    %   5. test duration
    %   6. choice (1=right, 2=left)
    %   7. Pupil term

    session_data_fit = table2array(session_data_tbl(:, [2 4 7 6 9 10 33]));
    uniqueHazards = cfg.hazard.codes;

    % Set up data matrix -- see buildBehaviorDataToFit.m /
    % buildInteractionDataToFit.m for column definitions. Unlike
    % behaviorLogisticFitsNeuralTerm.m, both the control and pupil-term
    % models here use all trials (no PREF-direction restriction).
    if control
        [data_to_fit, Lpre] = buildBehaviorDataToFit(session_data_fit);
    else
        [data_to_fit, Lpre] = buildInteractionDataToFit(session_data_fit, session_data_fit(:,7));
    end

    % Fit/plot separately for each hazard

    clf % clear figure

    for hh = 1:2
        
        %Lneural = logical(neural_split(:,hh));

        Lhazard = session_data_fit(:,1) == uniqueHazards(hh);
        
        % Get the fit 

        if control
            [fits_control(uu,:,hh), LL(uu,:,hh)] = logistFitDotsrev(data_to_fit(Lhazard,:));
        else
            [fits(uu,:,hh), LL(uu,:,hh)] = logistFitDotsrevNP(data_to_fit(Lhazard,:));
        end

        % Compute Tjur's D (Tjur's coefficient of disctrimination -- R^2)
             % Compares average fitted probability of the two response outcomes (0 and 1)
             % Calculate mean of predicted probabilities of each outcome and take the difference between the 2 means
             % 0 <= R^2 >= 1 where no disriminating power = 0 and perfect discriminating power = 1

         if control
            fit_y = logistValDotsrev(fits_control(uu,:,hh)', data_to_fit(Lhazard,1:end-1));
            R_sq_control(uu,:,hh) = computeTjursR2(fit_y, data_to_fit(Lhazard,4));
         else
             fit_y = logistValDotsrevNP(fits(uu,:,hh)', data_to_fit(Lhazard,1:end-1));
             R_sq(uu,:,hh) = computeTjursR2(fit_y, data_to_fit(Lhazard,5));
         end

        % Plot separately for predir=L and predir=R

        for xx = 1:2
            figure(2)
            subplot(1,2,xx); hold on; box on;
            plot([0 0], [0 1], 'k--');
            plot(xax([1 end]), [0.5 0.5], 'k--', 'MarkerSize', 15);

            if control
                choice_data = data_to_fit(Lhazard&Lpre(:,xx), 4);
                
                % Show smoothed data
                [sorted_time_axis, I] = sort(data_to_fit(Lhazard&Lpre(:,xx), 3));
            else
                choice_data = data_to_fit(Lhazard&Lpre(:,xx), 5);

                % Show smoothed data
                [sorted_time_axis, I] = sort(data_to_fit(Lhazard&Lpre(:,xx), 4));
            end

            sorted_choice_data = choice_data(I);
            plot(sorted_time_axis, nanrunmean(sorted_choice_data,10), ':', 'Color', colors{hh});

            hold on;
            
            if xx==1
                title(sprintf('Prefinal dir=R'))
                ylabel('Fraction choose Left')
            else
                title(sprintf('Prefinal dir=L'))
                ylabel('Fraction choose Right')

            end
            xlabel('Signed time (-non_switch, +switch)');

            % Show fits

            xax_fits = (-1200:1200)' ;
            show_fit_data = repmat(cat(2, ones(length(xax_fits),2), xax_fits, xax_fits),[1 1 2]);
            show_fit_data(:,2,2) = 0;

            if control
                ys = logistValDotsrev(fits_control(uu,:,hh)', show_fit_data(:,:,xx));

            else
                ys = logistValDotsrevNP(fits(uu,:,hh)', show_fit_data(:,:,xx));
            end

            plot(xax_fits, ys, '-', 'Color', colors{hh}, 'LineWidth', 2);
            axis([xax_fits(1) xax_fits(end) 0 1]);

            drawnow

            hold on
        end     
    end

    if control
        % NOTE: was reading from fits/R_sq here (a copy-paste leftover from
        % the non-control branch) instead of fits_control/R_sq_control --
        % since fits/R_sq are only ever populated in the else branch below,
        % this silently produced all-NaN fits_pupilTerm_control output.
        fits_pupilTerm_control(uu,1,:) = fits_control(uu,3,:);
        R_sq_pupilTerm_control(uu,1,:) = R_sq_control(uu,1,:);
    else
        fits_pupilTerm(uu,1,:) = fits(uu,3,:);
        R_sq_pupilTerm(uu,1,:) = R_sq(uu,1,:);
    end

    clear pupil_data_evoked pupil_evoked_bin pupil_baseline_resid
end