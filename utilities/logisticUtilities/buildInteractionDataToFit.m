function [data_to_fit, Lpre] = buildInteractionDataToFit(session_data, interactionSource, interactionMask)
%BUILDINTERACTIONDATATOFIT Build the 5-column neural/pupil-interaction design matrix.
%
% [data_to_fit, Lpre] = buildInteractionDataToFit(session_data, interactionSource)
% builds the design matrix used by the neural- and pupil-term models:
%   P(switch) = logistic(beta0 + beta1*dir_bias + beta2*interaction*signed_time + beta3*signed_time)
%
% session_data must have columns [hazard, coherence, dir_prefinal,
% dir_final, duration, choice] (only columns 3, 4, 5, 6 are used here) --
% the layout shared by behaviorLogisticFitsNeuralTerm.m and
% behaviorLogisticFitsPupilTerm.m. interactionSource is the per-trial
% neural firing rate (FR_testing) or evoked pupil value (pupil_evoked_bin).
%
% data_to_fit columns:
%   1. Switch/stay bias ... column of ones
%   2. Right/left bias ... 0 when prefinal dir=L, 1 when prefinal dir=R
%   3. Interaction ... interactionSource * signed time
%   4. Signed time ... neg=non-switch, pos=switch
%   5. Choice ... 0/1 = did not/did choose prefinal dir
%
% Lpre is an [n x 2] logical mask, [prefinal-dir-is-right, prefinal-dir-is-left],
% used by the calling script to split plots/data by prefinal direction.
%
% [...] = buildInteractionDataToFit(session_data, interactionSource, interactionMask)
% additionally multiplies the interaction term by a per-trial
% interactionMask, e.g. to restrict behaviorLogisticFitsNeuralTerm.m to
% PREF-direction trials only (interactionMask = FR_dir_idx). Omit or pass
% [] for no restriction (the default), matching
% behaviorLogisticFitsPupilTerm.m.

if nargin < 3 || isempty(interactionMask)
    interactionMask = 1;
end

dirPrefinal = session_data(:,3);
dirFinal = session_data(:,4);
choice = session_data(:,6);

LdirSwitch = dirPrefinal ~= dirFinal;
LpreDirIsRight = dirPrefinal < 90 | dirPrefinal >= 270;
Lpre = [LpreDirIsRight ~LpreDirIsRight];
LchoseRight = choice == 1;

data_to_fit = ones(size(session_data,1), 4);
data_to_fit(~LpreDirIsRight, 2) = 0;
data_to_fit(:,4) = session_data(:,5);
data_to_fit(~LdirSwitch,4) = -data_to_fit(~LdirSwitch,4);
data_to_fit(:,3) = data_to_fit(:,4) .* interactionSource .* interactionMask;
data_to_fit(LpreDirIsRight & LchoseRight,5) = 0; % prefinal dir = R and chose R
data_to_fit(LpreDirIsRight &~LchoseRight,5) = 1; % prefinal dir = R and chose L
data_to_fit(~LpreDirIsRight& LchoseRight,5) = 1; % prefinal dir = L and chose R
data_to_fit(~LpreDirIsRight&~LchoseRight,5) = 0; % prefinal dir = L and chose L

end
