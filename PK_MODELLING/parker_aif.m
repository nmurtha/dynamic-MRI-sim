function [ Cb ] = parker_aif( t )
%function [ Cb ] = parker_aif( t )
%
%Returns a population-averaged AIF. 
% Input
%   t   -> time in seconds
%
% Output
%   Cb  -> blood concentration in mmol
%
% Example usage:
%   t = linspace(0,6,75);
%   Cb = parker_aif(t);
%   plot(t,Cb);
% 
% See "Experimentally-derived functional form for a population-averaged 
%   high-temporal-resolution arterial input function for dynamic contrast-
%   enhanced MRI" Magn Reson Med. 2006 Nov;56(5):993-1000.

% Change dimensions to minutes (as per assumption in Parker et al).
t = t ./ 60;

% parameters
A1      = 0.809;    % mmol*min
A2      = 0.330;    % mmol*min
T1      = 0.17046;  % min
T2      = 0.365;    % min
sig1    = 0.0563;   % min
sig2    = 0.132;    % min
alpha   = 1.050;    % mmol
beta    = 0.1685;   % min^-1
s       = 38.078;   % min^-1
tau     = 0.483;    % min

% Empirically derived function, consisting of two gaussians and an
% exponential decay modulated by a sigmoid. 
Cb = A1*normpdf(t,T1,sig1) +...
     A2*normpdf(t,T2,sig2) +...
     alpha*exp(-beta*t)./(1+exp(-s*(t-tau))); % mmol

end

