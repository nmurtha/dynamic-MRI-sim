function [ paramMap ] = Get_Tofts_Parameters( reconImages, sliceNumber,...
    effTimes, mask )

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Feb 6, 2017
%
%  PURPOSE: Recovers extended Tofts parameters from reconstructed image
%           timeseries, multiplying the timeseries by a mask to reduce the
%           number of pixels fit (if desired).
%
%  INPUTS:
%    reconImages     -> A 4D stack of reconstructed volume images.
%
%    sliceNumber     -> The 2D slice within reconImages for which
%                       parameters will be fit. I.e., the fitting will
%                       focus on reconImages(:,:,sliceNumber,:).
%
%    effTimes        -> The effective sampling times of each reconstructed
%                       image.
%
%    samplingTimes   -> The vector of sampling times originally used in
%                       generating the contrast curves (i.e. in the
%                       numerical phantom). To minimize inherent error due 
%                       to approximation of numerical integral, the data is
%                       interpolated back onto the fine grid before
%                       fitting. That way the fitter has consistent
%                       numerical error due to integral step size.
%
%    mask            -> A 2D array containing a mask of zeros and ones.
%                       Used to blank out any pixels not in the features of
%                       interest.
%
%  OUTPUTS: 
%    paramMap        -> A map of Ktrans, kep and Vp parameters.                     
%
%   EXTERNALS:
%     parker_aif.m
%     fit_extended_tofts_conv.m
%     fit_extended_tofts_trapz.m

%% ------------------------------------------------------------------------
%  ---------------------   APPLY MASK TO TIME SERIES   --------------------
%  ------------------------------------------------------------------------

% Check input to ensure that the number of time points match the number of
% reconstructed images.
assert( size(reconImages,4) == numel(effTimes),...
    'Number of images must match number of time points.' );

% Get dimensions of reconstructed image array.
numRows   = size(reconImages,1);
numCols   = size(reconImages,2);
numVols   = size(reconImages,4);

% Get the slice time series. 
sliceEvolution = squeeze( reconImages(:,:,sliceNumber,:) );

% Multiply the slice time series by the mask. The effect will be to
% artificially eliminate any non-feature pixels, except if the mask
% contains entirely ones. A tailored mask may help to expedite parameter
% recovery (by avoiding useless pixels), though it certainly represents an
% ideal situation.
for iVol = 1:numVols
    sliceEvolution(:,:,iVol) = mask .* sliceEvolution(:,:,iVol);
end


%% ------------------------------------------------------------------------
%  -----------------------   RECOVER THE PARAMETERS   ---------------------
%  ------------------------------------------------------------------------

% Allocate memory for parameter map. First plane is initial amplitude,
% second plane is decay constant, third plane is the constant offset.
paramMap = zeros(numRows,numCols,3);

% Set threshold for fitter. If pixel timeseries is never higher than this
% threshold, it will not be fitted.
signalRange = max(reconImages(:)) - min(reconImages(:));
thresh = 0.01 * signalRange;

% Fit temporal parameters, pixel-by-pixel over the image.
h = waitbar(0,'Fitting PK Parameters...');
for r = 1:numRows
    waitbar(r/numRows)
    for c = 1:numCols
        
        % Get the "observed" data we are fitting the function to.
        fitData = squeeze( sliceEvolution(r,c,:) )';
        
        % If below threshold, don't fit and move on to next iteration.
        if max( abs(fitData) ) < thresh 
          continue;
        end     
        
        % Perform fit and assign to ParamMap.
        temp = fit_extended_tofts_intg(effTimes,fitData);   
        paramMap(r,c,:) = temp;
    end
end
close(h)


end

