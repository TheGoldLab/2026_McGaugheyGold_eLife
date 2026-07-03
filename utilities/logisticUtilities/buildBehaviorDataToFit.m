function [data_to_fit, Lpre] = buildBehaviorDataToFit(session_data, interactionMask)
%BUILDBEHAVIORDATATOFIT Build the 4-column behavior-only logistic design matrix.
%
% [data_to_fit, Lpre] = buildBehaviorDataToFit(session_data) builds the
% design matrix used by the behavior-only logistic model:
%   P(switch) = logistic(beta0 + beta1*dir_bias + beta2*signed_time)
%
% session_data must have columns [hazard, coherence, dir_prefinal,
% dir_final, duration, choice] (only columns 3, 4, 5, 6 are used here) --
% the layout shared by behaviorLogisticFits.m, behaviorLogisticFitsShuffle.m,
% and the control model in behaviorLogisticFitsNeuralTerm.m /
% behaviorLogisticFitsPupilTerm.m.
%
% data_to_fit columns:
%   1. Switch/stay bias ... column of ones
%   2. Right/left bias ... 0 when prefinal dir=L, 1 when prefinal dir=R
%   3. Signed time ... neg=non-switch, pos=switch
%   4. Choice ... 0/1 = did not/did choose prefinal dir
%
% Lpre is an [n x 2] logical mask, [prefinal-dir-is-right, prefinal-dir-is-left],
% used by the calling script to split plots/data by prefinal direction.
%
% [data_to_fit, Lpre] = buildBehaviorDataToFit(session_data, interactionMask)
% multiplies the signed-time column by a per-trial interactionMask before
% fitting, e.g. to restrict behaviorLogisticFitsNeuralTerm.m's control
% model to PREF-direction trials only (interactionMask = FR_dir_idx).
% Omit or pass [] for no restriction (the default), matching
% behaviorLogisticFitsPupilTerm.m's control model.

if nargin < 2 || isempty(interactionMask)
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
data_to_fit(:,3) = session_data(:,5);
data_to_fit(~LdirSwitch,3) = -data_to_fit(~LdirSwitch,3);
data_to_fit(:,3) = data_to_fit(:,3) .* interactionMask;
data_to_fit(LpreDirIsRight & LchoseRight,4) = 0; % prefinal dir = R and chose R
data_to_fit(LpreDirIsRight &~LchoseRight,4) = 1; % prefinal dir = R and chose L
data_to_fit(~LpreDirIsRight& LchoseRight,4) = 1; % prefinal dir = L and chose R
data_to_fit(~LpreDirIsRight&~LchoseRight,4) = 0; % prefinal dir = L and chose L

end
