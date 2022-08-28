function Ct = calc_extended_tofts_intg(params, t)
%
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           May 22, 2017          
%
%  PURPOSE: Calculate contrast evolution according to the extended Tofts 
%           model:
%             C(t) = Vp*AIF + Ktrans*convolve(exp(-t*kep), AIF)
%
%             where kep = Ktrans / Ve;
%  
%  INPUTS: 
%   params       -> Vector of PK parameters, [Ktrans, Ve, Vp].
%
%   t            -> Times for which Tofts model is to be calculated.
%
%  OUTPUTS:
%   Ct           -> Vector of contrast evolution data according to extended
%                   Tofts model, units in mmol.
%
%  EXTERNALS:
%   parker_aif.m

% Extract parameters. 
Ktrans = params(1);
Ve     = params(2);
Vp     = params(3);

% Define local function handle for integrand so that it may be passed to 
% numerical integration. Note that tau <= t, since this is called in an
% integral from tau = 0 to t (i.e. the convolution!).
convIntegrand = @(t,tau) parker_aif(tau) .* exp(-(Ktrans/Ve).*(t - tau));

% Get Ct value for each time requested in t. Must use loop because integral
% cannot handle vectors of integration limits...
Ct = zeros(1,numel(t));
for iTime = 1:numel(t)
    % Get numerical estimate to the convolution integral, integrating tau 
    % from 0 to t.
    convIntegral = integral(...               % Numerical integrator...
        @(tau)convIntegrand(t(iTime),tau),... % integrate over tau at t = t(iTime)
        0,...                                 % lower bound of integration on tau
        t(iTime),...                          % upper bound of integration on tau
        'RelTol',1e-6);                       % Tolerance on integration...   
    
    % Assemble the pieces to get the final Tofts model output!
    Ct(iTime) = Vp.*parker_aif(t(iTime)) + Ktrans.*convIntegral;
end
end


