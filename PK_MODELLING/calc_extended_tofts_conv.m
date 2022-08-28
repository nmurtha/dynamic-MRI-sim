function Ct = calc_extended_tofts_conv(params, t)
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
%           estimated with conv.
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

% Get time step to scale convolution with (otherwise conv won't approximate
% the actual continuous integral).
if numel(t) > 1
    dT = t(2) - t(1);
else
    dT = 1;
end

% Get Tofts model contrast evolution.
temp = conv( AIF, exp(-(Ktrans/Ve).*t), 'full').*dT;
Ct   = Vp.*AIF + Ktrans.*temp(1:numel(t));

end


