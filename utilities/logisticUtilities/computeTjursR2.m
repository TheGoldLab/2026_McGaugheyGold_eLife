function R2 = computeTjursR2(fit_y, choice01)
%COMPUTETJURSR2 Tjur's coefficient of discrimination (pseudo-R^2).
%
% R2 = computeTjursR2(fit_y, choice01) compares the average fitted
% probability of the two observed outcomes (0 and 1): the mean predicted
% probability among actual-1 trials minus the mean predicted probability
% among actual-0 trials. 0 <= R2 <= 1, where 0 = no discriminating power
% and 1 = perfect discrimination.
%
% fit_y: per-trial fitted probabilities (from logistValDotsrev /
%   logistValDotsrevNP).
% choice01: per-trial observed outcome (0/1), the data_to_fit choice
%   column the fit was evaluated against.

R2 = abs(mean(fit_y(choice01 == 1,:)) - mean(fit_y(choice01 == 0,:)));

end
