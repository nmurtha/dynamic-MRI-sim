function [meanErrors] = Get_Mean_Errors_Parallel(maskedErrs)
%
%  WRITTEN BY: Nathan Murtha (nathan.j.murtha@gmail.com)
%              Feb 17, 2017
%
%  PURPOSE: Gets the mean errors, with std, of the recovered parameters for
%           each feature in a phantom.
%
%  INPUTS:  
%    maskedErrs -> A cell array, one cell per phantom feature, containing 
%                  masked slice percent errors (i.e. all pixels outside of
%                  the feature of interest was set to zero).
%
%  OUTPUTS:
%    meanErrors -> Mean error information for the features. A cell array,
%                  one cell containing a matrix, per feature. The first row
%                  of the contained matrix represents the mean error of
%                  that parameter, the second row represents the std in
%                  that error.
%                  E.g. meanErrors{1}(1,:) -> mean errors for feature 1
%                       meanErrors{1}(2,:) -> std in errors for feature 1
%

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR EXECUTION   -----------------------
%  ------------------------------------------------------------------------

% Initialize the output cell array.
meanErrors = cell(1,length(maskedErrs));

%% ------------------------------------------------------------------------
%  -------------------------   GET MEAN ERRORS  ---------------------------
%  ------------------------------------------------------------------------

% Loop over the BART masked errors, finding the nonzero indices in each
% plane and fiding the mean and std of the corresponding voxels.
parfor iFeat = 1:length(maskedErrs)
    % Find the number of planes in this features' error map.
    numPlanes = size(maskedErrs{iFeat},3);
    
    % For each error plane, find the mean error and std of that error.
    for iPlane = 1:numPlanes
        % Pull the current error planes into temporary variables.
        tempPlane = maskedErrs{iFeat}(:,:,iPlane);
        
        % Find the nonzero indices in the mask.
        indsNonZero = find(tempPlane ~= 0);
        
        % Find the mean error. First row is the mean, second row is the std.     
        meanErrors{iFeat}(1,iPlane) = mean(tempPlane(indsNonZero));       
        meanErrors{iFeat}(2,iPlane) = std(tempPlane(indsNonZero));
    end
    
end