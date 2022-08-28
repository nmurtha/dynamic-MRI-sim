function [meanTimecourse] = Get_Mean_Timecourse_Parallel(imageSeries, sliceNumber,...
    featureMask)
%
%  WRITTEN BY: Nathan Murtha (nathan.j.murtha@gmail.com)
%              Feb 13, 2017
%
%  PURPOSE: Gets the mean timecourse of all nonzero entries in the
%           appropriate slice of the image series, following application of
%           a mask to isolate the feature of interest. Also returns the std
%           of the timecourse values.
%
%  INPUTS:  
%    imageSeries  -> 4D array of 3D images evolving over time.
%
%    sliceNumber  -> The slice for which the timeseries will be sought,
%                    following application of the featureMask.
%
%    featureMask  -> Mask to isolate the feature of interest.
%
%  OUTPUTS:
%    meanTimecourse -> Error map between images. 
%
%  EXTERNAL FUNCTIONS:
%   [None]

%% ------------------------------------------------------------------------
%  ------------------------   SET UP FOR EXECUTION   ----------------------
%  ------------------------------------------------------------------------

% Pull the appropriate slice evolution from the image series.
sliceEvolution = squeeze( imageSeries(:,:,sliceNumber,:) );

% Initialize matrix for output results.
meanTimecourse = zeros(2,size(sliceEvolution,3));

% Initialize dummy variables for parfor loop (avoids indexing issues).
tempMean = zeros(1,size(sliceEvolution,3));
tempStd  = zeros(1,size(sliceEvolution,3));

%% ------------------------------------------------------------------------
%  -----------------------   GET MEAN TIMECOURSES   -----------------------
%  ------------------------------------------------------------------------

% Loop over the time evolutions and fill with mean and std of images.
parfor iVol = 1:size(sliceEvolution,3)
    % Find nonzero indices of current slice after applying the mask.
    tempTimeCourseData = featureMask .* sliceEvolution(:,:,iVol);
    indsNonZero = find(tempTimeCourseData ~= 0);
    
    % Get mean and std of current feature image.
    tempMean(1,iVol) = mean( tempTimeCourseData(indsNonZero) );
    tempStd(1,iVol)  = std( tempTimeCourseData(indsNonZero) );
end

% Place temporary parfor results into output result.
meanTimecourse(1,:) = tempMean;
meanTimecourse(2,:) = tempStd;

end