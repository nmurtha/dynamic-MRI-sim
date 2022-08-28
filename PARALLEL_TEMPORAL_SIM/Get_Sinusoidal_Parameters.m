function [ paramMap ] = Get_Sinusoidal_Parameters( reconImages, sliceNumber,...
    effTimes, mask )

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Feb 6, 2017
%
%  PURPOSE: Recovers linear parameters from reconstructed image
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
%    paramMap        -> A map of parameters describing the sinusoidal model
%                       S(t) = A * sin(2pi/B * (t-C)) + D
%                       NOTE: Plane 1 -> A coefficients.
%                             Plane 2 -> B coefficients.
%                             Plane 3 -> C coefficients.
%                             Plane 4 -> D coefficients.
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
paramMap = zeros(numRows,numCols,4);

% Define fit functions. Will use two methods and take the one that produces 
% the most accurate results. The first form eliminates phase offset using 
% the trig identity sin(w(t-t')) = sin(wt)cos(wt')-sin(wt')cos(wt). The 
% second is of the form sin(w(t-t')). Note that w = 2*pi/period.
fitNoPhase   = @(a,tuniform) a(1).*sin(2*pi/a(2)*tuniform)...
                           - a(3).*cos(2*pi/a(2)*tuniform) + a(4);    
                       
fitWithPhase = @(b,tuniform) b(1).*sin(2*pi/b(2).*(tuniform - b(3)) )+b(4); 

% Set threshold for fitter. If pixel timeseries is never higher than this
% threshold, it will not be fitted.
thresh = 0.001 * max( abs( reconImages(:) ) );

% Set options for the fitter. 
options.Display           ='none';
options.MaxIterations     = 300;
options.StepTolerance     = 1e-7;
options.FunctionTolerance = 1e-7;

% Fit temporal parameters, pixel-by-pixel over the image.
for r = 1:numRows
    for c = 1:numCols
        
        % Get the "observed" data we are fitting the function to.
        fitData = squeeze( sliceEvolution(r,c,:) )';
        
        % If below threshold, don't fit and move on to next iteration.
        if max( abs(fitData) ) < thresh 
          continue;
        end
        
        % Resample data onto uniform, fine grid. This will allow the FFT to
        % be easily taken in the event that the time spacings aren't
        % uniform (as in sampling with CIRCUS patterns).
        desiredFs = (1/(effTimes(2) - effTimes(1))) * 2;      
        [yuniform, tuniform] = resample(fitData,effTimes,...
            desiredFs,'spline');
        
        % Guess initial period by taking the Fourier spectra, and taking
        % the first dominant non-zero frequency (avoids the large zero
        % frequency component introduced by constant offset).
        Y = fft(yuniform); 
        f = ( 0:(length(tuniform)-1) ) * desiredFs/length(tuniform);
        
        freqAxis = f( 1:ceil( length(tuniform)/2 ) ); % Get relevant range of frequencies (only as good as Nyquist)
        fftPos = Y( 1:length(freqAxis) ); % Get relevant portion of frequency spectra (clips out redundant frequencies)
        
        [allMaxima,maximaInds] = findpeaks( abs(fftPos) ); % Get all maxima of spectrum, as well as the indices marking their location.
        ind = find( (allMaxima - max(allMaxima)) == 0 ); % Index in maximaInds marking the index of relevant peak within freqAxis
        if isempty(ind)
            ind = 2;
            initPeriod = abs(1/(freqAxis(ind))); % if frequency resolution is too low, must assume first non-zero frequency is dominant...
        else
            initPeriod = abs(1/(freqAxis(maximaInds(ind)))); % Get period from fft spectra
        end
        
        % Guess initial amplitude based on range of values.
        initAmp = (max(yuniform) - min(yuniform))/2;

        % Guess initial vertical offset based on mean value of data.
        initVertOffset = mean(yuniform);
            
        % Set initial parameter guesses. Assume no phase offset to begin.
        a0 = [initAmp, initPeriod, 0, initVertOffset];
        b0 = [initAmp, initPeriod, 0, initVertOffset];
                
        % Set very generous upper and lower bounds.
        aLowBound = [0, 0, 0, -Inf];
        aUpBound  = [Inf, Inf, Inf, Inf];
        
        bLowBound = [0, 0, -Inf, -Inf];
        bUpBound  = [Inf, Inf, Inf, Inf];
        
        % Perform fit and assign to ParamMap
        a = lsqcurvefit(fitNoPhase, a0,...
                        tuniform, yuniform,...
                        aLowBound, aUpBound,...
                        options);
                    
        b = lsqcurvefit(fitWithPhase, b0,...
                        tuniform, yuniform,...
                        bLowBound, bUpBound,...
                        options);
                    
        % Convert a-parameters, those for the "phase-term-free" fit,
        % back to the form A(1)*sin( A(2)*(t-A(3)) )+A(4).
        A = zeros(1,4);
        A(4) = a(4); % Vertical offset remains unchanged.
        A(2) = a(2); % Period remains unchanged.
        
        if a(1) >= 0 % Adjust arctan based on sign of horizontal component
            A(3) = A(2)/(2*pi) * atan(a(3)/a(1)); % (1/w)*atan(a(3)/a(1))
        else % If horizontal component less than zero, must adjust for left quadrants
            A(3) = (pi + atan(a(3)/a(1))) * A(2)/(2*pi); % atan will return negative angle, so ADD it to pi.
        end
        
        A(1) = sqrt(a(1)^2 + a(3)^2); % magnitude is the root of the quadrature sum of individual magnitudes            
        
        % Choose the reconstruction with the lowest RMSE to timeseries
        % drawn from current voxel.
        recNoPhase   = A(1) .* sin( 2*pi/A(2).*(effTimes - A(3)) ) + A(4);
        recWithPhase = b(1) .* sin( 2*pi/b(2).*(effTimes - b(3)) ) + b(4);

        diffNoPhaseSquare   = (recNoPhase   - fitData).^2;
        diffWithPhaseSquare = (recWithPhase - fitData).^2;
        
        rmseNoPhase   = sqrt( mean( diffNoPhaseSquare(:) ) );
        rmseWithPhase = sqrt( mean( diffWithPhaseSquare(:) ) );       

        if rmseNoPhase > rmseWithPhase % if sin(w(t-t')) better
            paramMap(r,c,:) = b;
        else                           % if sin(wt)+cos(wt) better
            paramMap(r,c,:) = A;
        end
        
    end
end




end

