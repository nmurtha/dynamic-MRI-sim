function [ pyramidImages ] = getGaussianPyramid( img, numScales, pyrWindow )

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Feb 21, 2017
%
%  PURPOSE: Gets "Gaussian pyramid" breakdown of the image. Not a strict
%           Gaussian pyramid, per se, but follows the lowpass filtering
%           techniques used in the MS-SSIM (see MS-SSIM by Wang).
% 
%  INPUTS:
%   img         -> Image for which pyramid is sought. Can be 2D or 3D.
% 
%   numScales   -> The number of subsampled scales in the pyramid.
%
%   pyrWindow   -> Size of filter to be used in making the pyramid (e.g.
%                  pyrWindow = 2, assuming a square window).
% 
%  OUTPUTS:
%   pyramidImages -> Cell array containing images in the pyramid.
%
%  EXTERNAL FUNCTIONS:
%   [NONE]

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR EXECUTION   -----------------------
%  ------------------------------------------------------------------------

% Get dimensionality of the input image.
dims = numel(size(img));

% Get kernel, choosing between 3D and 2D based on dimensions of the image.
if dims == 2
    kernel = ones(pyrWindow,pyrWindow);
elseif dims == 3
    kernel = ones(pyrWindow,pyrWindow,pyrWindow);
else
    error('Only support 2D and 3D images.')
end
kernel = kernel ./ sum(kernel(:)); % averaging kernel, as per Wang.

% Check that image size is okay for the number of scales requested. 
sizeCondition = min(size(kernel)) <= min(size(img)) * 2^(-(numScales-1));
assert(sizeCondition,'Cant calculate %d scales with image of smallest dimension %d and kernel size %d.',...
    numScales,min(size(img)),min(size(kernel)));

%% ------------------------------------------------------------------------
%  -------------------------   MAKE THE PYRAMID   -------------------------
%  ------------------------------------------------------------------------

% Loop over the higher levels of the lowpass pyramid.
pyramidImages = cell(1,numScales);
pyramidImages{1} = img;
imgCurrent = img;
for iScale = 2:numScales
    % Get lowpass image. Symmetric is used for edge cases here to be
    % consistent with the pyramid scheme in Wang's MS-SSIM code. 
    imLowpass = imfilter(imgCurrent, kernel, 'symmetric');
    
    % Downsample lowpass image by a factor of two.
    if dims == 2
        imgCurrent = imLowpass(1:2:end,1:2:end);
    else 
        imgCurrent = imLowpass(1:2:end,1:2:end,1:2:end);
    end
    
    % Save current lowpass image.
    pyramidImages{iScale} = imgCurrent;
end

end

