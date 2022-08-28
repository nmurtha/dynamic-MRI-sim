function [ score, scaleScores ] = MSSSIM_ND( imgRef, imgDis, numScales,...
    ssimWinDims, pyrWindow, maxGrayScale )

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Feb 21, 2017. Edited April 4, 2017.
%
%  PURPOSE: Extend the MS-SSIM to 3D images. 
%
%           NOTE: See "Multi-Scale Structural Similarity For Image Quality 
%                 Assessment" by Z. Wang et al. 
% 
%  INPUTS:
%   imgRef      -> Reference image (2D or 3D array).  
%   
%   imgDis      -> Distorted image (2D or 3D array).
%
%   numScales   -> The number of scales in the lowpass pyramid that will
%                  be considered in calculating the MS-SSIM score. By
%                  default, is taken as the maximum number of scales that
%                  won't downsample the image to have min dimension less
%                  than the kernel size. 
%                 
%   ssimWinDims  -> A vector, [rows cols slices], that represents the
%                   extent of the gaussian window used to calcualate 
%                   local SSIM properties.
%
%   pyrWindow    -> Size of window to be used in making the
%                   lowpass pyramid. An averaging filter of size pyrWindow 
%                   by pyrWindow is made. By default, pyrWindow = 2.
%
%   maxGrayScale -> Maximum grayscale value in reference image. By
%                   default, taken as max(imgRef(:)). Determines the
%                   numerical stability constants.
% 
%  OUTPUTS:
%   score        -> Score for the MS-SSIM metric of the distorted image.
%
%   scaleScores  -> Score at each scale in the lowpass pyramid.
%
%
%  EXTERNALS:
%   getGaussianPyramid.m
%
%  INTERNALS:
%   getGaussianWeightingFilter

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR EXECUTION   -----------------------
%  ------------------------------------------------------------------------

% -----------------------   SAFETY CHECKS OF INPUTS   ---------------------
% Ensure images are same size.
dimsRef = size(imgRef);
dimsDis = size(imgDis);
assert(isequal(dimsRef, dimsDis),'Images must be the same size.')

% Ensure that no more than 5 pyramid scales are requested.
if exist('numScales','var')
    assert(numScales <= 5,'Only support 5 pyramid scales (empirical MS-SSIM weights run out).')
end

% ------------------------   SET DEFAULTS IF NEEDED   ---------------------
% Check pyramid window side length, set to 2 by default.
if ( ~exist('pyrWindow','var') )
    pyrWindow = 2;
end

% Check number of pyramid scales, set to maximum number possible allowed by
% dimensions of images and size of pyramid kernel by default. 
if (~exist('numScales','var'))
    % Get smallest dimension in images.
    dimMin = min(dimsRef(:));
    
    % Get number of scales that can be accomodated. Note that the
    % requirement leading to this is dimMin*2^-numScales <= pyrWindow. The
    % ceil function is applied to take only an integer number of scales.
    numScales = floor( -log2(pyrWindow/dimMin) ) + 1; % +1 to include base of pyramid, level 0
    numScales = min(5,numScales); % supports no more than 5 scales
end

% Set ssimWinDims to default 11x11(xN) window if not specified, where N is
% equal to either 11 OR the number of slices in the image if the image has
% fewer than 11 slices.
if ~exist('ssimWinDims','var')
    % Get smallest dimension of image at the maximum level in pyramid.
    dimMin = min(dimsRef(:)) * 2^(-(numScales-1)); % -1 to exclude level 0 of pyramid
    
    % Set dimensions of ssim gaussian window
    if ndims(imgRef) == 3 % 3D image
        ssimWinDims = [11 11 min(floor(dimMin),11)];
    elseif ndims(imgRef) == 2 % 2D image
        ssimWinDims = [11 11];
    else
        error('Only support 2D and 3D images.')
    end
end

% Set max grayscale value in image if not specified.
if ~exist('maxGrayScale','var')
    maxGrayScale = max(imgRef(:));
end

% ------------------------   INITIALIZE CONSTANTS   -----------------------
% Make kernel for local ssim calculation. As per SSIM norm, use a gaussian
% window when calculating local luminance, contrast and structure. 
% NOTE: fixed standard deviation of 1.5, as per Wang.
sigma     = 1.5;
gaussFilt = getGaussianWeightingFilter(ssimWinDims,ndims(imgRef),sigma);

% Set default stability parameters from Wang's MS-SSIM paper.
K  = [0.01 0.03];
C1 = (K(1)*maxGrayScale)^2;
C2 = (K(2)*maxGrayScale)^2;
C3 = C2/2;

% Set pyramid scale weights from Wang's MS-SSIM paper. The assumption is
% that luminance, contrast and structure are all equally weighted within a
% single pyramid scale; these weights were empirically found based on that
% assumption.
scaleWeights = [0.0448 0.2856 0.3001 0.2363 0.1333]; 

% Get weights for each scale, stopping at a max of 5 scales. Sum of weights
% is normalized to unity (only the relative weightings are important).
scaleWeights = scaleWeights(1:numScales);
scaleWeights = scaleWeights./sum(scaleWeights);

%% ------------------------------------------------------------------------
%  ------------------------   GET MS-SSIM SCORE   -------------------------
%  ------------------------------------------------------------------------

% Get lowpass (pseudo-Gaussian) Pyramid of each image.
refPyramid = getGaussianPyramid(imgRef, numScales, pyrWindow); 
disPyramid = getGaussianPyramid(imgDis, numScales, pyrWindow);   

% Get structure and contrast maps for each image scale. Luminance is only
% calculated at the highest scale, whereas contrast and structure are
% calculated at each scale.
contrastAndStructMap = cell(1,numScales);
for iScale = 1:numScales
    % Get image at current scale in pyramid.
    refPyrTemp = refPyramid{iScale};
    disPyrTemp = disPyramid{iScale};
    
    % Get mean maps of images. 
    muRef = imfilter(refPyrTemp, gaussFilt, 'replicate');
    muDis = imfilter(disPyrTemp, gaussFilt, 'replicate');
    
    % Get contrast and structure terms for this scale. Note that: 
    %   cov = E[XY] - E[X]E[Y]
    %   var = E[X^2] - E[X]^2
    % where E is the expectation operator. Set any negative values of
    % variance to zero.
    sigmaRefDis  = imfilter(refPyrTemp.*disPyrTemp, gaussFilt, 'replicate') - muRef.*muDis;
    sigmaRef_Sq  = imfilter(refPyrTemp.^2, gaussFilt, 'replicate') - muRef.^2;
    sigmaDis_Sq  = imfilter(disPyrTemp.^2, gaussFilt, 'replicate') - muDis.^2;
    sigmaRef_Sq  = max(0, sigmaRef_Sq);
    sigmaDis_Sq  = max(0, sigmaDis_Sq);
    
    % Get contrast and structure map at current scale.
    contrast  = (2*sqrt(sigmaRef_Sq.*sigmaDis_Sq) + C2) ./ (sigmaRef_Sq + sigmaDis_Sq + C2);
    structure = (sigmaRefDis + C3) ./ (sqrt(sigmaRef_Sq.*sigmaDis_Sq) + C3);
    contrastAndStructMap{iScale} = contrast .* structure;   
    
    % If at the highest scale, get the luminance map.
    if (iScale == numScales)
        luminance = (2*muRef.*muDis + C1)./(muRef.^2 + muDis.^2 + C1);
    end
end

% Pool SSIM maps across scales to get the metric score.
scaleScores = zeros(1,numScales);
for iScale = 1:numScales
    % Get structure and contrast map for this scale.
    tempScaleMap = contrastAndStructMap{iScale};
    
    % If at the final scale, incorporate luminance into score map.
    if (iScale==numScales)
        tempScaleMap = tempScaleMap .* luminance;
    end
    
    % Get score at this scale. 
    scaleScores(iScale) = mean(tempScaleMap(:));
end

% Get final score score, weighting each scale score.
score = prod(scaleScores(1:numScales) .^ scaleWeights(1:numScales));

end

%% ------------------------------------------------------------------------
%  --------------------------   HELPER FUNCTIONS   ------------------------
%  ------------------------------------------------------------------------

function gaussFilt = getGaussianWeightingFilter(filtSize,numDims,sigma)
%
% Returns a Gaussian wighting kernel with size specified by the vector
% filtSize, with standard deviation in weighting elements given by sigma.

% Check that dimensions of image and gaussian kernel are compatible (e.g. 
% not asking for 3D kernel with 2D image).
assert(numel(filtSize) == numDims,...
    'Dimensions of image must match the dimension of gaussian kernel.');

if (numDims == 2) % 2D case
    gaussFilt = fspecial('gaussian',[filtSize(1) filtSize(2)],sigma);
    
elseif (numDims == 3) % 3D case
    % Get size of the 3D Gaussian mask in each dimension. If size is even
    % in a given direction, adjust the kernel accordingly.
    if mod(filtSize(1),2) == 1 % odd size
        ytop  = -floor(filtSize(1)/2);
        ybottom = floor(filtSize(1)/2);
    else % even size
        ytop  = -floor(filtSize(1)/2) + 1;
        ybottom = floor(filtSize(1)/2);
    end
    
    if mod(filtSize(2),2) == 1 % odd size
        xleft  = -floor(filtSize(2)/2);
        xright = floor(filtSize(2)/2);
    else % even size
        xleft  = -floor(filtSize(2)/2) + 1;
        xright = floor(filtSize(2)/2);
    end
    
    if mod(filtSize(3),2) == 1 % odd size
        zrear  = -floor(filtSize(3)/2);
        zfront = floor(filtSize(3)/2);
    else % even size
        zrear  = -floor(filtSize(3)/2) + 1;
        zfront = floor(filtSize(3)/2);
    end
        
    % Get 3D Gaussian mask. If the entries are extremely small, just zero
    % them out.
    [x,y,z]   = ndgrid(ytop:ybottom,xleft:xright,zrear:zfront);
    arg       = -(x.^2 + y.^2 + z.^2)/(2*sigma.^2);
    gaussFilt = exp(arg);
    gaussFilt(gaussFilt<eps*max(gaussFilt(:))) = 0;     
end

% Normalize kernel to unit sum.
sumFilt   = sum(gaussFilt(:));
gaussFilt = gaussFilt/sumFilt;

end
