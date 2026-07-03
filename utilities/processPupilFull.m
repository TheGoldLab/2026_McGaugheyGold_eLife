function [fits, pupil_bin_mean_save_LSF, pupil_bin_mean_save_HSF, ...
    pupil_bin_slope_save_LSF, pupil_bin_slope_save_HSF, ...
    pupil_baseline_save_LSF, pupil_baseline_save_HSF, ...
    ses_pupil_bin_mean_save, ses_pupil_bin_mean_ses_diff] = ...
    processPupilFull(dat, slide, bin_size, start_time, end_time, save_trial_data, outPath)
% processPupilFull Full pupil processing with regression modeling
%
% This function processes pupil data, extracts mean and slope measures,
% and fits a regression model for each session to predict pupil measures
% from hazard rate and baseline.
%
% Inputs:
%   dat              - Pre-processed data table
%   slide            - Bin slide/step size in msec (e.g., 10)
%   bin_size         - Bin width in msec (e.g., 100)
%   start_time       - Analysis start time in msec (e.g., -200)
%   end_time         - Analysis end time in msec (e.g., 3600)
%   save_trial_data  - Boolean, if true saves trial-by-trial data to outPath
%   outPath          - Path for saving trial-by-trial data (required if save_trial_data=true)
%
% Outputs:
%   fits                        - [numSessions x numBins x 3 x 2] Regression coefficients
%   pupil_bin_mean_save_LSF     - [numSessions x numBins] Mean pupil for LSF
%   pupil_bin_mean_save_HSF     - [numSessions x numBins] Mean pupil for HSF
%   pupil_bin_slope_save_LSF    - [numSessions x numBins] Pupil slope for LSF
%   pupil_bin_slope_save_HSF    - [numSessions x numBins] Pupil slope for HSF
%   pupil_baseline_save_LSF     - [numSessions x numBins] Baseline for LSF
%   pupil_baseline_save_HSF     - [numSessions x numBins] Baseline for HSF
%   ses_pupil_bin_mean_save     - [numSessions x numBins x 2] Session averages
%   ses_pupil_bin_mean_ses_diff - [numSessions x numBins x 2] LSF-HSF differences

%% Setup
uniqueSessionNames = unique(dat.ses_ID);
numSessions = length(uniqueSessionNames);
uniqueHazards = [2 6];

%% Binning scheme
bins = cat(2, ...
    (start_time:slide:end_time-bin_size)', ...
    (start_time+bin_size:slide:end_time)');
tax = start_time:end_time;
xax = mean(bins,2);
numBins = size(bins, 1);
axlen = length(tax);

xax_plot = xax - 2400; % Align to testing stimulus onset

%% Filter parameters (butterworth low-pass)
srate = 1000;
lpass = 5;
freq = (lpass/srate)*2;
n_order = 1;
[butter_pb, butter_pa] = butter(n_order, freq, 'low');

%% Filter parameters for pupil slope
FILTER_HW = 75;
filter_c = (3 / (2*FILTER_HW^3 + 3*FILTER_HW^2 + FILTER_HW)) * (FILTER_HW:-1:-FILTER_HW);

%% Initialize outputs
fits = nan(numSessions, numBins, 3, 2);
hazard_order = zeros(numSessions, 1);

pupil_bin_mean_save_LSF = nan(numSessions, numBins);
pupil_bin_mean_save_HSF = nan(numSessions, numBins);
pupil_bin_slope_save_LSF = nan(numSessions, numBins);
pupil_bin_slope_save_HSF = nan(numSessions, numBins);
pupil_baseline_save_LSF = nan(numSessions, numBins);
pupil_baseline_save_HSF = nan(numSessions, numBins);

ses_pupil_bin_mean_save = nan(numSessions, numBins, 2);
ses_pupil_bin_mean_ses_diff = nan(numSessions, numBins, 2);

%% Loop through sessions
for ss = 1:numSessions
    
    Lses = strcmp(dat.ses_ID, uniqueSessionNames{ss}) & dat.coh_final >= 50;
    
    if sum(Lses) == 0
        continue
    end
    
    disp(uniqueSessionNames{ss})
    
    %% Get session data
    session_data = table2array(dat(Lses, [2 4 5 9 12 27 13 19 18 27 28]));
    dirs = table2array(dat(Lses, 'dir_final'));
    unique_dirs = unique(dirs);
    num_trials = size(session_data, 1);
    
    %% Create hazard logicals
    Lhazard = false(num_trials, 2);
    for hh = 1:2
        Lhazard(:,hh) = session_data(:,1) == uniqueHazards(hh);
    end
    
    %% Extract pupil and eye data
    pupil_data = nan(num_trials, axlen);
    eyeX_data = nan(num_trials, axlen);
    eyeY_data = nan(num_trials, axlen);
    dots_on_times = dat(Lses, 31).dots_on;
    Lgood = false(num_trials, axlen);
    
    for tt = 1:num_trials
        start_index = max(1, round(dots_on_times{tt}(1) + start_time - session_data(tt,10)));
        end_index = min(start_index + axlen - 1, round(session_data(tt,9) - session_data(tt,10)));
        from_indices = start_index:end_index;
        to_indices = 1:length(from_indices);
        
        pupil_data(tt, to_indices) = dat(Lses, 24).pupil_diam{tt}(from_indices);
        Lgood(tt, to_indices) = true;
        eyeX_data(tt, to_indices) = dat(Lses, 25).pupil_horiz{tt}(from_indices);
        % Dynamic fieldname with a leading space is deliberate: this column's
        % stored name has incidental whitespace ('pupil_vert' doesn't match
        % it) -- see the comment in buildMergedTableTiers.m.
        eyeY_data(tt, to_indices) = dat(Lses, 26).(' pupil_vert'){tt}(from_indices);
    end
    
    %% Interpolate missing pupil values
    pupil_data = pupil_data';
    Lbad = ~isfinite(pupil_data);
    pupil_data(Lbad) = interp1(find(~Lbad), pupil_data(~Lbad), find(Lbad), 'linear', 'extrap');
    
    %% Low-pass filter pupil
    pupil_data = pupil_data';
    pupil_filt = pupil_data;
    pupil_filt(~isfinite(pupil_filt)) = mean(pupil_filt(:), 'omitnan');
    pupil_filt1 = filtfilt(butter_pb, butter_pa, pupil_filt')';
    
    %% Use residuals to find bad values, remove, and interpolate
    pupil_resid = (pupil_data - pupil_filt1)';
    pupil_resid_mean = mean(pupil_resid(:), 'omitnan');
    pupil_resid_std = std(pupil_resid(:), 'omitnan');
    Lbad = (pupil_resid > pupil_resid_mean + (2*pupil_resid_std)) | ...
           (pupil_resid < pupil_resid_mean - (2*pupil_resid_std));
    
    pupil_clean = pupil_data';
    pupil_clean(Lbad) = nan;
    pupil_clean(Lbad) = interp1(find(~Lbad), pupil_clean(~Lbad), find(Lbad), 'linear', 'extrap');
    pupil_clean = pupil_clean';
    
    %% Re-filter cleaned data
    pupil_filt2 = filtfilt(butter_pb, butter_pa, pupil_clean')';
    
    %% Process eye position data
    eyeX_data = eyeX_data';
    eyeX_data(Lbad) = interp1(find(~Lbad), eyeX_data(~Lbad), find(Lbad));
    eyeX_data = (eyeX_data' - mean(eyeX_data(:), 'omitnan')) ./ std(eyeX_data(:), 'omitnan');
    
    eyeY_data = eyeY_data';
    eyeY_data(Lbad) = interp1(find(~Lbad), eyeY_data(~Lbad), find(Lbad));
    eyeY_data = (eyeY_data' - mean(eyeY_data(:), 'omitnan')) ./ std(eyeY_data(:), 'omitnan');
    
    %% Subtract average and z-score pupil
    pupil_data = pupil_filt2;
    pupil_data(~Lgood) = nan;
    
    % Subtract running average timecourse
    pupil_data = pupil_data - repmat(smooth(mean(pupil_data, 2, 'omitnan'), 50), 1, size(pupil_data, 2));
    
    % Z-score
    pupil_data = (pupil_data - mean(pupil_data(:), 'omitnan')) ./ std(pupil_data(:), 'omitnan');
    
    %% Get running slope of pupil
    pupil_tmp = filter(filter_c, 1, pupil_data')';
    pupil_slope = nan(size(pupil_data));
    pupil_slope(:, (1+FILTER_HW):(end-FILTER_HW-1)) = pupil_tmp(:, length(filter_c)+1:end);
    
    %% Create regression matrix
    % Columns: bias, baseline, hazard rate (LSF=1, HSF=0)
    regression_matrix = cat(2, ...
        ones(num_trials, 2), ...
        double(Lhazard(:,1)));
    
    baselines = cat(2, ...
        mean(pupil_data(:, 100:200), 2, 'omitnan'), ...
        mean(pupil_slope(:, 100:200), 2, 'omitnan'));
    
    %% Loop through time bins and fit regression
    pupil_bin_mean_save = nan(num_trials, numBins, 2);
    pupil_baseline_mean_save = nan(num_trials, numBins, 2);
    
    for bb = 1:numBins
        
        Lbin = tax >= bins(bb,1) & tax < bins(bb,2);
        
        % Fit mean and slope
        for xx = 1:2
            if xx == 1
                pupil_bin_mean = mean(pupil_data(:, Lbin), 2, 'omitnan');
                pupil_bin_mean_save(:,bb,xx) = mean(pupil_data(:, Lbin), 2, 'omitnan');
            else
                pupil_bin_mean = mean(pupil_slope(:, Lbin), 2, 'omitnan');
                pupil_bin_mean_save(:,bb,xx) = mean(pupil_slope(:, Lbin), 2, 'omitnan');
            end
            
            pupil_bin_mean = pupil_bin_mean - baselines(:,xx);
            pupil_bin_mean_save(:,bb,xx) = pupil_bin_mean_save(:,bb,xx) - baselines(:,xx);
            
            Lgood_fit = isfinite(pupil_bin_mean);
            
            regression_matrix(:,2) = baselines(:,xx);
            pupil_baseline_mean_save(:,bb,xx) = baselines(:,xx);
            
            % Run regression
            [fits(ss,bb,:,xx), ~, ~, ~, stats] = regress(pupil_bin_mean(Lgood_fit), ...
                regression_matrix(Lgood_fit, :));
            
            if stats(1) == Inf || stats(1) == -Inf
                stats(1) = nan;
            end
        end
    end
    
    %% Save trial-by-trial data (optional)
    if save_trial_data
        pupil_bin_mean_save_name = [uniqueSessionNames{ss}, '_pupil_bin_mean_trials_correct'];
        save(fullfile(outPath, pupil_bin_mean_save_name), 'pupil_bin_mean_save');
    end
    
    %% Save session-level summaries
    ses_pupil_bin_mean_save(ss,:,:) = mean(pupil_bin_mean_save, 'omitnan');
    ses_pupil_bin_mean_ses_diff(ss,:,:) = mean(pupil_bin_mean_save(Lhazard(:,1),:,:), 'omitnan') - ...
                                           mean(pupil_bin_mean_save(Lhazard(:,2),:,:), 'omitnan');
    
    pupil_bin_mean_save_LSF(ss,:) = mean(pupil_bin_mean_save(Lhazard(:,1),:,1), 'omitnan');
    pupil_bin_mean_save_HSF(ss,:) = mean(pupil_bin_mean_save(Lhazard(:,2),:,1), 'omitnan');
    
    pupil_bin_slope_save_LSF(ss,:) = mean(pupil_bin_mean_save(Lhazard(:,1),:,2), 'omitnan');
    pupil_bin_slope_save_HSF(ss,:) = mean(pupil_bin_mean_save(Lhazard(:,2),:,2), 'omitnan');
    
    pupil_baseline_save_LSF(ss,:) = mean(pupil_baseline_mean_save(Lhazard(:,1),:,1), 'omitnan');
    pupil_baseline_save_HSF(ss,:) = mean(pupil_baseline_mean_save(Lhazard(:,2),:,1), 'omitnan');
    
    %% Track hazard order
    if mean(find(Lhazard(:,2)), 'omitnan') < mean(find(Lhazard(:,1)), 'omitnan')
        hazard_order(ss) = 1;
    end
end

fprintf('Full pupil processing complete for %d sessions\n', numSessions);

end
