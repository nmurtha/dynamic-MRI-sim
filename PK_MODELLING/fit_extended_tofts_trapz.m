function params = fit_extended_tofts_trapz(timesSampled,CtSampled)
%
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           May 12, 2017
%
%  PURPOSE: Fit contrast evolution data to the extended Tofts model:
%             C(t) = Vp*AIF + Ktrans*convolve(exp(-t*kep), AIF)
%           Where kep = Ktrans / Ve. This calls the Tofts model calculator
%           that uses a trapezoidal estimation for convolution integral.
%  
%  INPUTS: 
%   timesSampled  -> Vector of times at which sampled contrast was
%                    obtained.
%
%   CtSampled     -> Vector of contrast evolution data for a single voxel 
%                    of a DCE-MRI dataset.
%
%  OUTPUTS:
%   params        -> Recovered PK parameters in the form [Ktrans, Ve, Vp].
%
%  EXTERNALS:
%   calc_extended_tofts_trapz.m

% Set preferences for the minimizing routine.
options = optimoptions('fmincon',...
    'Algorithm'             , 'interior-point',...
    'MaxIterations'         , 500,...
    'MaxFunctionEvaluations', 400,...
    'OptimalityTolerance'   , 1e-7,...
    'StepTolerance'         , 1e-7,...
    'ConstraintTolerance'   , 1e-6,...
    'Display'               , 'off');

% Set up the minimization function. Want to minimize the sum of the square 
% of the residuals, a nonlinear least-squares fitting problem. 
% NOTE: x represents the fitted parameter vector.
tofts = @calc_extended_tofts_trapz; 
fun   = @(x)sum( (tofts(x,timesSampled) - CtSampled).^2 );

% Set the parameter guesses and upper/lower bounds, as well as linear
% inequality and equality constraints. Notation used here follows the
% documentation for fmincon with Mathworks.
x0    = [1, 0.5, 0.25];            % Initial [Ktrans, Ve, Vp] guess.
A     = [0 1 1];                   % Linear inequality constraint parameter indices.
b     = 1;                         % Inequality constraints, such that A*x <= b.
Aeq   = [];                        % Linear equality constraint parameter indices.
beq   = [];                        % Equality constraints, such that Aeq*x = beq.
lb    = [0 0 0];                   % Lower bound on [Ktrans, Ve, Vp].
ub    = [4 1 1];                   % Upper bound on [Ktrans, Ve, Vp].
nonlcon = [];                      % Nonlinear constraints.

% Minimize the objective function with specified constraints.
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

% Extract fitted PK parameters for output.
params(1) = x(1);       % ktrans
params(2) = x(2);       % Ve
params(3) = x(3);	    % Vp

end