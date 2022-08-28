function [errMap] = guardedDividePercentError(reconParam, trueParam)
%
%  WRITTEN BY: Nathan Murtha (nathan.j.murtha@gmail.com)
%              Feb 13, 2017
%
%  PURPOSE: Gets percent error of reconstructed parameter map, avoiding
%           divison by zero. If the true parameter happens to be 0, and the 
%           reconstructed parameter is nonzero, then the absolute error of
%           the parameter is placed in the resulting error map instead (the
%           user MUST be aware of which true parameters were 0 in such a
%           case).
%
%  INPUTS:  
%    reconParam   -> Image that will be divided by another image.
%
%    trueParam    -> Image being divided into another image.
%
%  OUTPUTS:
%    errMap       -> Error map between images. 

% Catch bad input.
assert(isequal(size(reconParam),size(trueParam)),...
    'Size of numerator and denominator images must be the same.');

% Find indices corresponding to nonzero entries in the denominator.
indsNonZeroDen = find(trueParam ~= 0);

% Divide only the nonzero indices. 
errMap = zeros(size(reconParam));

errMap(indsNonZeroDen) = ((reconParam(indsNonZeroDen) - trueParam(indsNonZeroDen))...
    ./ trueParam(indsNonZeroDen)) .* 100;

% Catch cases where numerator was nonzero while denominator WAS zero.
indsZeroDen    = find(trueParam == 0);
indsNonZeroNum = find(reconParam ~= 0);
indsAbsError   = intersect(indsZeroDen,indsNonZeroNum); % contains all indices for which the true param WAS zero, while the recon WAS NOT.

% Put absolute error at indices where true param was 0 while recon was not.
errMap(indsAbsError) = reconParam(indsAbsError) - trueParam(indsAbsError);
