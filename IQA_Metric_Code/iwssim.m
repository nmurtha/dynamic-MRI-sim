function [IWSSIM, iw_map, IWMSE, IWPSNR]= iwssim(imgRef, imgDis, numScales,...
    dynamicRange, flagIW, K, scaleWeights, windowKernel, windowSizeIW, parentIW, sigmaNsq)

%
% INPUTS:
%  imgRef, imgDis -> images for which the metric is sought.
%  numScales      -> Scales in the pyramid for which metric will be
%                    calculated.
%  dynamicRange   -> Max pixel value in the images.
%  flagIW         -> Flag that turns on or off info weighting.
%  K              -> K parameters for stability constant calculation.
%  scaleWeights   -> weights that will be used at each scale in the
%                    pyramid, for it's SSIM contribution to final score.
%  windowKernel   -> The kernel used for SSIM calculation.
%  windowSizeIW   -> Size of the window used for info mapping.
%  parentIW       -> Flag to include parent coefficient (from scale below)
%                    in the info mapping.
%  sigmaNsq       -> Std of human visual system noise. 
%
% OUTPUTS:
%  IWSSIM         -> the info weighted SSIM.
%  iw_map         -> info maps.
%  IWMSE          -> info weighted mean squared error.
%  IWPSNR         -> info weighted peak signal to noise ratio.

%% ------------------------------------------------------------------------
%  ----------------------   SET UP FOR EXECUTION   ------------------------
%  ------------------------------------------------------------------------

% Set a max value that the iwpsnr can take... 
MAX_PSNR = 1000;

% Must have between 2 and 11 arguments.
if (nargin < 2 || nargin > 11)
    IWSSIM = -Inf;
    IWMSE  = -Inf;
    IWPSNR = -Inf;
    disp('Input argument error');
    return;
end

% Image sizes must match.
if (size(imgRef) ~= size(imgDis))
    IWSSIM = -Inf;
    IWMSE  = -Inf;
    IWPSNR = -Inf;
    disp('Input argument error: images should have the same size');
    return;
end

% Set flag for information weighting to ON if not specified. Otherwise, the
% MS-SSIM will be output.
if (~exist('flagIW','var'))
   flagIW = 1;
end

% Set 5 Laplacian pyramid scales as default.
if (~exist('numScales','var'))
   numScales = 5;
end

% Take default K parameters from SSIM paper.
if (~exist('K','var'))
   K = [0.01 0.03]; % default from [Wang et al, IEEE-TIP 2004]
end

% Take dynamic range to be 255 by default.
if (~exist('dynamicRange','var'))
   dynamicRange = 255;
end

% Take default MS-SSIM scale weights unless specified otherwise.
if (~exist('scaleWeights','var'))
   scaleWeights = [0.0448 0.2856 0.3001 0.2363 0.1333]; % default from [Wang et al, IEEE-Asilomar 2003]
end

% Make Gaussian kernel, unless user specified a kernel (note that user
% specified kernel must be square).
if (exist('windowKernel','var'))
    winSz = size(windowKernel);
    if (winSz(1) ~= winSz(2))
        IWSSIM = -Inf;
        IWMSE  = -Inf;
        IWPSNR = -Inf;
        disp('Window size error: window must be square.');
        return;
    end
    winSize = winSz(1);
else
    winSize = 11; % default from [Wang et al, IEEE-TIP 2004]
    gwin_std = 1.5; % default from [Wang et al, IEEE-TIP 2004]
    windowKernel = fspecial('gaussian', winSize, gwin_std); % default from [Wang et al, IEEE-TIP 2004]
    windowKernel = windowKernel/sum(windowKernel(:));
end

if (~exist('windowSizeIW','var'))
   windowSizeIW = 3; % spatial neighborhood size (%in current scale, 3x3?)
end

if (~exist('parentIW','var'))
   parentIW = 1; % include parent neighbor (% from scale above?)
end

% Set "visual noise" variance.
if (~exist('sigmaNsq','var'))
   sigmaNsq = 0.4; % default from [Sheikh & Bovik, IEEE-TIP 2006]
end

% Set local block size dimensions.
blockSzX = windowSizeIW; 
blockSzY = windowSizeIW;

% Get weights for each scale, stopping at a max of 5 scales. Sum of weights
% is normalized to unity (only the relative weightings are important).
numScales    = min(5,numScales);
scaleWeights = scaleWeights(1:numScales);
scaleWeights = scaleWeights./sum(scaleWeights);

% Ensure that image is large enough such that the kernel can be applied
% even at the smallest scale image in the pyramid.
if (min(size(imgRef)) < winSize*(2^(numScales-1)))
    IWSSIM = -Inf;
    IWMSE  = -Inf;
    IWPSNR = -Inf;
    disp(['Image size too small for ' num2str(numScales) ' scale IW-SSIM evaluation.']);
    return;
end

% Get boundaries for where values in maps will be computed (?)
bound = ceil((winSize-1)/2);
bound1 = bound - floor((blockSzX-1)/2);

%% ------------------------------------------------------------------------
%  -----------------------   CALCULATE THE IWSSIM   -----------------------
%  ------------------------------------------------------------------------

% Get laplacian pyramid of the images.
[pyro,pind]= buildLpyr(double(imgRef),numScales);
[pyrd,pind]= buildLpyr(double(imgDis),numScales);

% Get SSIM maps at each level in the pyramid. 
[cs_map, l_map, se_map] = scale_quality_maps(pyro,pyrd,pind,numScales,K,dynamicRange,windowKernel);

% If info weighting is desired, get info weight maps at each scale.
if (flagIW)
    iw_map = info_content_weight_map(pyro,pyrd,pind,numScales,parentIW,blockSzX,blockSzY,sigmaNsq);
end

for s = 1:numScales
    % Get structure and contrast maps for this scale.
    se = se_map{s};
    cs = cs_map{s};
    
    % If at the top scale, incorporate luminance.
    if (s==numScales)
        cs = cs.*l_map;
    end
    if (flagIW)
        if (s<numScales)
            iw = iw_map{s};
            iw = iw(bound1+1:end-bound1, bound1+1:end-bound1);
        else
            iw = ones(size(cs));
        end
        se = se(bound+1:end-bound, bound+1:end-bound);        
        wmcs(s) = sum(sum(cs.*iw))/sum(sum(iw));
        wmse(s) = sum(sum(se.*iw))/sum(sum(iw));
    else
        wmcs(s) = mean2(cs);
        wmse(s) = mean2(se);
    end
end

% Get metric scores. Note that the last scale in the SSIM was made to have
% luminance contributions, as it should.
IWMSE = prod(wmse(1:numScales).^scaleWeights(1:numScales));
IWPSNR = min(10*log10(255^2/IWMSE), MAX_PSNR);
IWSSIM = prod(wmcs(1:numScales).^scaleWeights(1:numScales));
