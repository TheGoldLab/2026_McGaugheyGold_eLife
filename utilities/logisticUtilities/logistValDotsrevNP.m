function ps_ = logistValDotsrevNP(fits, data)
%
% Uses the logistic model to compute p (probability of
% making a particular decision) for the given data (m rows
% of trials x n columns of data categories) and parameters
% (nx1). Assumes logit(p) (that is, ln(p/(1-p)) is a linear
% function of the data.
% fits(1:2) are bias terms, encoded wrt to X (not Y) intercepts
% fits(3) is pupil/neural term
% fits(4) is slope
% fits(5) is upper/lower asymptote (lapse rate)

% Modified to include neural and pupil terms

ps_ = fits(5) + (1 - 2.*fits(5))./(1+exp(-(data(:,end-1:end)-(data(:,1:2)*fits(1:2)))*fits(3:4)));
