function [ssimval, ssimmap] = SSIM_ND(imgRef, imgDis,...
    ssimWinDims, exponents, maxGrayScale)
%
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           April 26, 2017
%
%  PURPOSE: Calculates the mean SSIM between a reference and a distorted
%           image.
%   
%           NOTE: This is a modified version of ssim.m (provided by
%                 Mathworks) which allows for assymetric metric window in
%                 3D. See "Image Quality Assessment: From Error Visibility
%                 to Structural Similarity" by Z. Wang et al.
% 
%  INPUTS:
%   imgRef       -> Reference image (either 2D or 3D array).  
%
%   imgDis       -> Distorted image (either 2D or 3D array). 
% 
%   ssimWinDims  -> A vector, [rows cols slices], that represents the
%                   extent of the gaussian window used to calcualate 
%                   local SSIM properties.
%
%   exponents    -> A vector representing relative weights between
%                   luminance, contrast, and structure terms, in that 
%                   order. Default is [1 1 1].
%
%   maxGrayScale -> Maximum grayscale value in reference image. By
%                   default, taken as max(imgRef(:)). Determines the
%                   numerical stability constants.
%
%  OUTPUTS:
%   score        -> Score for the SSIM metric of the distorted image.
%
%   ssimMap      -> SSIM map between the images.
%
%  EXTERNALS:
%   [NONE]
%
%  INTERNALS: 
%    getGaussianWeightingFilter

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR EXECUTION   -----------------------
%  ------------------------------------------------------------------------

% -----------------------   SAFETY CHECKS OF INPUTS   ---------------------
% Check that images are the same size.
dimsRef = size(imgRef);
dimsDis = size(imgDis);
assert(isequal(dimsRef, dimsDis),'Images must be the same size.')

% Check if any exponents are negative.
if exist('exponents','var')
   assert(isequal(exponents>=0,ones(1,3)),'Exponents must be positive');
end

% -----------------------   SET DEFAULTS IF NEEDED   ----------------------
% Set ssimWinDims to default 11x11(xN) window if not specified, where N is
% equal to either 11 OR the number of slices in the image if the image has
% fewer than 11 slices.
if ~exist('ssimWinDims','var')
    if ndims(imgRef) == 3 % 3D image
        ssimWinDims = [11 11 min(dimsRef(3),11)];
    elseif ndims(imgRef) == 2
        ssimWinDims = [11 11];
    else
        error('Only support 2D and 3D images.')
    end
end

% Set exponents to be used in the SSIM (i.e. weights between luminance,
% contrast and structure terms), if not specified. 
if ~exist('exponents','var')
    exponents = [1 1 1];
end

% Set max grayscale of the image if not specified.
if ~exist('maxGrayScale','var')
    maxGrayScale = max(imgRef(:));
end

% ------------------------   INITIALIZE CONSTANTS   -----------------------

% Set stability terms. Wang sets C3 = C2/2 by default. 
K         = [0.01 0.03];
C1        = (K(1)*maxGrayScale)^2;
C2        = (K(2)*maxGrayScale)^2;
C         = [C1 C2 C2/2];

% Make gaussian kernel for local ssim calculations. As per Wang, use a 
% Gaussian window when calculating local luminance, contrast and structure. 
sigma     = 1.5;
gaussFilt = getGaussianWeightingFilter(ssimWinDims,ndims(imgRef),sigma);

%% ------------------------------------------------------------------------
%  ------------------------   GET MEAN SSIM SCORE   -----------------------
%  ------------------------------------------------------------------------

% Get weighted mean maps for the images.
muRef = imfilter(imgRef ,gaussFilt, 'replicate');
muDis = imfilter(imgDis ,gaussFilt, 'replicate');

% Get contrast and structure terms for the images. Note that: 
%   cov = E[XY] - E[X]E[Y]
%   var = E[X^2] - E[X]^2
% where E is the expectation operator. Set any negative values of
% variance to zero.
sigmaRef_Sq = imfilter(imgRef.^2     ,gaussFilt,'replicate') - muRef.^2;
sigmaDis_Sq = imfilter(imgDis.^2     ,gaussFilt,'replicate') - muDis.^2;
sigmaRefDis = imfilter(imgRef.*imgDis,gaussFilt,'replicate') - muRef.*muDis;
sigmaRef_Sq = max(0, sigmaRef_Sq);
sigmaDis_Sq = max(0, sigmaDis_Sq);

% Get luminance (eqn 6), contrast (eqn 9) and structure (eqn 10) terms.
luminance = (2*muRef.*muDis + C(1)) ./ (muRef.^2 + muDis.^2 + C(1));
contrast  = (2*sqrt(sigmaRef_Sq.*sigmaDis_Sq) + C(2)) ./ (sigmaRef_Sq + sigmaDis_Sq + C(2));
structure = (sigmaRefDis + C(3)) ./ (sqrt(sigmaRef_Sq.*sigmaDis_Sq) + C(3));

% Get SSIM map (eqn 12).
ssimmap = (luminance) .^exponents(1) .*... 
          (contrast)  .^exponents(2) .*... 
          (structure) .^exponents(3);

% Pool the SSIM map to get the mean SSIM (eqn 17).
ssimval = mean(ssimmap(:));
    
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

