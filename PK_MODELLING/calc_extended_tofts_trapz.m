function Ct = calc_extended_tofts_trapz(params, t)
%
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           May 24, 2017          
%
%  PURPOSE: Calculate contrast evolution according to the extended Tofts 
%           model:
%               C(t) = Vp*AIF + Ktrans*convolve(exp(-t*kep), AIF)
%           where kep = Ktrans / Ve. 
%
%           The numerical integration involved in the convolution is 
%           estimated with a simple sinusoidal rule.
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

% Get AIF.
AIF = parker_aif(t);

% Calculate the contrast evolution by looping through the convolution
% integral at each time, using trapz to estimate the numerical integral at 
% each time point iTime to get C(time = T).
% NOTE: the first point is set manually, since trapz will fuss if tau
% contains only one time point. This is a quick approximation only.
Ct    = zeros(1,numel(t));
Ct(1) = Vp*AIF(1);
for iTime = 2:numel(t)
    
    % Grab times in integrand up to iTime index (0 to t in integral). 
    % Convolution integrates from time tau = 0 to time tau = T to get the 
    % result at time T, so we have to step tau between 0 and T to get C(T).
    tau = t(1:iTime);
    
    % Grab AIF values that correspond to times involved in current step of
    % convolution. 
    % Note: Cp is shorthand for "concentration in plasma", another notation
    % often used for AIF. 
    Cp = AIF(1:iTime);
    
    % Calculate terms in the convolution integrand. Note that tau(end)
    % corresponds to tau = T. 
    integrand = Cp .* exp( -(Ktrans/Ve)*(tau(end)-tau) );
    
    % Approximate the integral over range of times involved.
    currentConv = trapz(tau,integrand);
      
    % Store extended Tofts model Ct value at current time given by iTime.
    Ct(iTime) = Vp*AIF(iTime) + Ktrans.*currentConv;
end

end