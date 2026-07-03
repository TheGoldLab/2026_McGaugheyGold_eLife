function [d, CI, SE] = computeCohenDCI(x1, x2, varargin)
%
% call: [d, CI, SE] = computeCohenDCI(x1, x2, varargin)
%
% EFFECT SIZE of the difference between the two
% means of two samples, x1 and x2 (that are vectors),
% computed as "Cohen's d".
%
% If x1 and x2 can be either two independent or paired
% samples, and should be treated accordingly:
%
% [d, CI, SE] = computeCohenDCI(x1, x2, 'independent'); [default]
% [d, CI, SE] = computeCohenDCI(x1, x2, 'paired');
%
% OUTPUTS:
%   d  - Cohen's d effect size
%   CI - 95% confidence interval [lower, upper]
%   SE - Standard error of Cohen's d
%
% Optional: Specify confidence level (default 0.95)
% [d, CI, SE] = computeCohenDCI(x1, x2, 'independent', 0.99);
%
% Note: according to Cohen and Sawilowsky:
%
% d = 0.01 --> very small effect size
% d = 0.20 --> small effect size
% d = 0.50 --> medium effect size
% d = 0.80 --> large effect size
% d = 1.20 --> very large effect size
% d = 2.00 --> huge effect size
%
%
% Ruggero G. Bettinardi (RGB)
% Cellular & System Neurobiology, CRG
% -------------------------------------------------------------------------------------------
%
% Code History:
%
% 25 Jan 2017, RGB: Function is created
% 07 Oct 2025: Added confidence interval and standard error outputs

if nargin < 3
    testType = 'independent';
    confLevel = 0.95;
elseif nargin < 4
    testType = varargin{1};
    confLevel = 0.95;
else
    testType = varargin{1};
    confLevel = varargin{2};
end

% basic quantities:
n1 = numel(x1);
n2 = numel(x2);
mean_x1 = nanmean(x1);
mean_x2 = nanmean(x2);
var_x1 = nanvar(x1);
var_x2 = nanvar(x2);
meanDiff = (mean_x1 - mean_x2);

% select type of test:
isIndependent = strcmp(testType, 'independent');
isPaired = strcmp(testType, 'paired');

% compute 'd' accordingly:
if isIndependent
    sv1 = ((n1-1)*var_x1);
    sv2 = ((n2-1)*var_x2);
    numer = sv1 + sv2;
    denom = (n1 + n2 - 2);
    pooledSD = sqrt(numer / denom); % pooled Standard Deviation
    s = pooledSD; % re-name
    d = meanDiff / s; % Cohen's d (for independent samples)
    
    % Standard error for independent samples
    SE = sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)));
    
    % Degrees of freedom
    df = n1 + n2 - 2;
    
elseif isPaired
    haveNotSameLength = ~isequal(numel(x1), numel(x2));
    if haveNotSameLength
        error('In a paired test, x1 and x2 have to be of same length!')
    end
    
    deltas = x1 - x2; % differences
    sdDeltas = nanstd(deltas); % standard deviation of the differences
    s = sdDeltas; % re-name
    d = meanDiff / s; % Cohen's d (paired version)
    
    % Standard error for paired samples
    n = n1; % same as n2 for paired
    SE = sqrt((1/n) + (d^2 / (2*n)));
    
    % Degrees of freedom
    df = n - 1;
end

% Compute confidence interval
alpha = 1 - confLevel;
t_crit = tinv(1 - alpha/2, df);
CI = [d - t_crit * SE, d + t_crit * SE];

end