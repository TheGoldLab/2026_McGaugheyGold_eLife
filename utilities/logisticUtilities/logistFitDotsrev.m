function [fits_, fval] = logistFitDotsrev(data)

% NEW USING GLOBALSEARCH   
% Define the error function

errFcn = @(fits) logistErrDotsrev(logistValDotsrev(fits', data(:,1:end-1)), data(:,end));

% Guess slope

g_sl = data(:,end-1)\((6 * data(:,end)) - 3);

% Set up the optimization problem

problem = createOptimProblem('fmincon',    ...
   'objective',   errFcn,           ... % Use the objective function
   'x0',          [   0    0 g_sl   0.02], ... % Initial conditions
   'lb',          [-300 -300 0.001  0.01], ... % Parameter lower bounds
   'ub',          [ 300  300 0.05   0.15], ... % Parameter upper bounds
   'options',     optimoptions(@fmincon,    ... % "function minimization with constraints"
   'Algorithm',   'active-set',  ...
   'Display',     'off',         ...
   'MaxIter',     5000,          ...
   'MaxFunEvals', 5000));

% Create a GlobalSearch object
% NumStageOnePoints/NumTrialPoints use MATLAB's own documented defaults
% (200/1000) rather than the 10x-inflated values (2000/10000) this file
% previously carried -- this is a well-behaved, low-dimensional (4-5
% parameter) smooth objective, so the default search breadth is expected
% to be robust; validated against the pre-computed data/behaviorFits/ fits.

gs = GlobalSearch("NumStageOnePoints", 200, 'NumTrialPoints', 1000, 'Display', 'off');
   
% Run it, returning the best-fitting parameter values and the negative-
% log-likelihood returned by the objective function

[fits_, fval] = run(gs,problem);
