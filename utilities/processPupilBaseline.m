function baseline_pupil = processPupilBaseline(dat)
% processPupilBaseline Extract baseline pupil measurements
%
% This function processes pupil data to extract baseline pupil measurements
% for each trial. The baseline is computed as the mean pupil diameter from
% 100-200ms after dots onset (after filtering and z-scoring).
%
% Inputs:
%   dat - Pre-processed data table (should be subset with 'BP' analysis type)
%
% Outputs:
%   baseline_pupil - [num_trials x 1] vector of baseline pupil values

%% Setup
uniqueSessionNames = unique(dat.ses_ID);
numSessions = length(uniqueSessionNames);

% Time parameters
start_time = -200;
end_time = 3600;
tax = start_time:end_time;
axlen = length(tax);

%% Filter parameters (butterworth low-pass)
srate = 1000;
lpass = 5;
freq = (lpass/srate)*2;
n_order = 1;
[butter_pb, butter_pa] = butter(n_order, freq, 'low');

%% Initialize output
baseline_pupil = [];

%% Loop through sessions
for ss = 1:numSessions
    
    Lses = strcmp(dat.ses_ID, uniqueSessionNames{ss});
    
    if sum(Lses) == 0
        continue
    end
    
    disp(uniqueSessionNames{ss})
    
    % Get session data
    session_data = table2array(dat(Lses, [27 18 28])); % trial_begins, dots_off, trial_ends
    num_trials = size(session_data, 1);
    
    % Extract pupil data
    pupil_data = nan(num_trials, axlen);
    dots_on_times = dat(Lses, 31).dots_on;
    Lgood = false(num_trials, axlen);
    
    for tt = 1:num_trials
        start_index = max(1, round(dots_on_times{tt}(1) + start_time - session_data(tt,1)));
        end_index = min(start_index + axlen - 1, round(session_data(tt,2) - session_data(tt,1)));
        from_indices = start_index:end_index;
        to_indices = 1:length(from_indices);
        pupil_data(tt, to_indices) = dat(Lses, 24).pupil_diam{tt}(from_indices);
        Lgood(tt, to_indices) = true;
    end
    
    %% Interpolate missing values
    pupil_data = pupil_data';
    Lbad = ~isfinite(pupil_data);
    pupil_data(Lbad) = interp1(find(~Lbad), pupil_data(~Lbad), find(Lbad), 'linear', 'extrap');
    
    %% Low-pass filter
    pupil_data = pupil_data';
    pupil_filt = pupil_data;
    pupil_filt(~isfinite(pupil_filt)) = mean(pupil_filt(:), 'omitnan');
    pupil_filt1 = filtfilt(butter_pb, butter_pa, pupil_filt')';
    
    %% Use residuals to find and remove outliers
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
    
    %% Subtract average across trials and z-score
    pupil_data = pupil_filt2;
    pupil_data(~Lgood) = nan;
    
    % Subtract running average timecourse
    pupil_data = pupil_data - repmat(smooth(mean(pupil_data, 2, 'omitnan'), 50), 1, size(pupil_data, 2));
    
    % Z-score
    pupil_data = (pupil_data - mean(pupil_data(:), 'omitnan')) ./ std(pupil_data(:), 'omitnan');
    
    %% Extract baseline (100-200ms after dots onset)
    baselines = mean(pupil_data(:, 100:200), 2, 'omitnan');
    
    % Append to output
    baseline_pupil = [baseline_pupil; baselines];
end

fprintf('Baseline pupil extracted for %d trials\n', length(baseline_pupil));

end
