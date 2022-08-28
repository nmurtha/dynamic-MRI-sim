function [ paramMap ] = Get_Exponential_Parameters( reconImages, sliceNumber,...
    effTimes, mask )

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Feb 6, 2017
%
%  PURPOSE: Recovers exponential parameters from reconstructed image
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
%    mask            -> A 2D array containing a mask of zeros and ones.
%
%  OUTPUTS: 
%    paramMap        -> A map of parameters describing the exponential
%                       model S(t) = A*exp(-t/B) + C. 
%                       NOTE: Plane 1 -> A coefficients.
%                             Plane 2 -> B coefficients.
%                             Plane 3 -> C coefficients.                     
%
%   EXTERNALS:
%     [NONE]

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

% Define fit function. Note that params represents the parameters of the
% contrast enhancement model, [amplitude, decay, offset].
fitFun = @(params,effTimes) params(1) .* exp( -1*effTimes./params(2) ) +...
    params(3);

% Set initial guess for parameters, as well as upper and lower bounds.
signalRange = max(reconImages(:)) - min(reconImages(:));
paramsInit     = [signalRange, 100, 0];
paramsLowBound = [0,0,0];
paramsUpBound  = [2*signalRange,Inf,2*signalRange]; % margin for error is allowed in signalRange estimate

% Set threshold for fitter. If pixel timeseries is never higher than this
% threshold, it will not be fitted.
thresh = 0.001 * signalRange;

% Set options for the fitter. 
options.Display           ='none';
options.MaxIterations     = 300;
options.StepTolerance     = 1e-7;
options.FunctionTolerance = 1e-7;

% Fit temporal parameters, pixel-by-pixel over the image.
parfor r = 1:numRows
    for c = 1:numCols
        
        % Get the "observed" data we are fitting the function to.
        fitData = squeeze( sliceEvolution(r,c,:) )';
        
        % If below threshold, don't fit and move on to next iteration.
        if max( abs(fitData) ) < thresh 
          continue;
        end
        
        % Perform fit and assign to ParamMap.
        temp = lsqcurvefit(fitFun,...
            paramsInit,...
            effTimes,...
            fitData,...
            paramsLowBound,...
            paramsUpBound,...
            options);
        
        paramMap(r,c,:) = temp;
    end
end




end

