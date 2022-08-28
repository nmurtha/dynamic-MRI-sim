function [score, errMap] = RMSE_ND(imgRef, imgDis)

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           April 26, 2017
%
%  PURPOSE: Calculate the root mean squared error (RMSE) between reference
%           distorted images.         
% 
%  INPUTS:
%   imgRef     -> Reference image (either 2D or 3D array).  
%
%   imgDis     -> Distorted image (either 2D or 3D array).  
% 
%  OUTPUTS:
%   score      -> Score for the RMSE metric of the distorted image.
%
%   errMap     -> Map made from (imgDis - imgRef)
%
%  EXTERNAL FUNCTIONS:
%   [NONE]

%% ------------------------------------------------------------------------
%  ----------------------   SET UP FOR EXECUTION   ------------------------
%  ------------------------------------------------------------------------

% Catch improper inputs on image sizes.
sizeCondition = ( size(imgRef) == size(imgDis) );
assert(min(sizeCondition) == 1,...
    'Reference and distorted image must be same sizes.');

%% ------------------------------------------------------------------------
%  ------------------------   CALCULATE THE RMSE   ------------------------
%  ------------------------------------------------------------------------

% Get error map.
errMap = imgDis - imgRef;

% Get RMSE.
errSquare = errMap.^2;
score     = sqrt(mean(errSquare(:)));

end