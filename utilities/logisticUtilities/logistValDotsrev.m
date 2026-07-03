function ps_ = logistValDotsrev(fits, data)
% function ps_ = logistValDotsrev(fits, data)
%
% Uses the logistic model to compute p (probability of
% making a particular decision) for the given data (m rows
% of trials x n columns of data categories) and parameters
% (nx1). Assumes logit(p) (that is, ln(p/(1-p)) is a linear
% function of the data.
% fits(1:2) are bias terms, encoded wrt to X (not Y) intercepts
% fits(3) is slope
% fits(4) is upper/lower asymptote (lapse rate)


ps_ = fits(4) + (1 - 2.*fits(4))./(1+exp(-fits(3).*(data(:,end)-(data(:,1:2)*fits(1:2)))));

