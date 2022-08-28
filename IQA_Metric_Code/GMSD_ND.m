function [score, GMS] = GMSD_ND(imgRef, imgDis, downsampleFlag, maxGrayScale)

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Jan 31, 2017
%
%  PURPOSE: Extend the gradient magnitude similarity deviation (GMSD)
%           to 3D images.
%
%           NOTE: See GMSD.m or "Gradient Magnitude Similarity Deviation: 
%                 A Highly Efficient Perceptual Image Quality Index" by W.
%                 Xue et al for GMSD info.
% 
%  INPUTS:
%   imgRef         -> Reference image (either 2D or 3D array).  
%
%   imgDis         -> Distorted image (either 2D or 3D array).  
%   
%   downsampleFlag -> Flag that turns on/off 2x downsampling of image. Off
%                     (0) by default.
%   
%   maxGrayScale   -> Maximum grayscale value in reference image. By
%                     default, taken as max(imgRef(:)). Determines the
%                     numerical stability constant used in GMSD 
%                     calculation.
% 
%  OUTPUTS:
%   score          -> Score for the GMSD metric of the distorted image.
%
%   GMS            -> Gradient Magnitude Similarity (GMS) map.
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

assert(numel(size(imgRef) == 3) || numel(size(imgRef) == 2), ...
    'GMSD calculation only supported for 2D or 3D images.');

% Set downsampling to be *off* by default.
if ~exist('downsampleFlag','var')
    downsampleFlag = 0;
end

% Set max grayscale value in image if not specified.
if ~exist('maxGrayScale','var')
    maxGrayScale = max(imgRef(:));
end

% Set parameter for numerical stability. Note that the authors intend this
% parameter to be 0.26% of the square of the dynamic range of the image.
% This is because they found 0.0026 to be appropriate for images in the
% range of [0,1], so scaling to [0,d] requires C to be scaled as well. In
% this case, C = 0.0026*(d*d), the authors claiming to use d^2 because of
% the denominator in equation 4, m_r^2 + m_d^2. This information was
% obtained through direct correspondence with the authors.
C = 0.0026 * maxGrayScale.^2; 

%% ------------------------------------------------------------------------
%  ---------------------   CASE FOR 3D INPUT IMAGE   ----------------------
%  ------------------------------------------------------------------------

if numel( size(imgRef) ) == 3
    % If user desired downsampling, downsample the image by factor of 2.
    if downsampleFlag == 1
        % Downsample by blurring with 2x2x2 averaging filter and tossing
        % out every second row/col/slice.
        downStep  = 2;
        aveKernel = ones(2,2,2)./8; % Note that each entry has equal weight
        aveTruth  = convn(imgRef, aveKernel, 'same');
        aveRecon  = convn(imgDis, aveKernel, 'same');
        imgRef    = aveTruth(1:downStep:end,1:downStep:end,1:downStep:end);
        imgDis    = aveRecon(1:downStep:end,1:downStep:end,1:downStep:end);
    end
    
    % Define gradient convolution kernel. Use 3D "Prewitt Filters", 
    % extension of those used by Xue et al.
    % NOTE: Author's mention that any gradient calculation method would be
    % acceptable, so could vary this if desired.
    filtX = zeros(3,3,3);
    filtY = zeros(3,3,3);
    filtZ = zeros(3,3,3);   
    for i = 1:3
        filtX(:,:,i) = [1/9  0 -1/9; 1/9 0 -1/9; 1/9   0  -1/9];                   
        filtY(:,:,i) = [1/9 1/9 1/9;  0  0   0 ;-1/9 -1/9 -1/9];                  
        filtZ(:,i,:) = [1/9  0 -1/9; 1/9 0 -1/9; 1/9   0  -1/9];
    end

    % Get gradients of each image by filtering with the gradient filters.
    TruthGx = imfilter(imgRef, filtX, 'replicate'); 
    TruthGy = imfilter(imgRef, filtY, 'replicate');     
    TruthGz = imfilter(imgRef, filtZ, 'replicate');     

    ReconGx = imfilter(imgDis, filtX, 'replicate'); 
    ReconGy = imfilter(imgDis, filtY, 'replicate');     
    ReconGz = imfilter(imgDis, filtZ, 'replicate'); 

    % Get gradient magnitude maps for each image.
    gradientMapTruth = sqrt(TruthGx.^2 + TruthGy.^2 + TruthGz.^2);
    gradientMapRecon = sqrt(ReconGx.^2 + ReconGy.^2 + ReconGz.^2);

    % Get GMS map from gradient magnitude images. In an ideal situation,
    % the images are "maximally similar" and the GMS map is 1 everywhere. 
    GMS = ( 2 * gradientMapTruth .* gradientMapRecon + C ) ./ ...
        ( gradientMapTruth.^2 + gradientMapRecon.^2 + C );

    % Get 3D GMSD score.
    score = std(GMS(:));
    
%% ------------------------------------------------------------------------
%  ---------------------   CASE FOR 2D INPUT IMAGE   ----------------------
%  ------------------------------------------------------------------------

elseif numel( size(imgRef) ) == 2
    if downsampleFlag == 1
        % Downsample by blurring with 2x2 averaging filter and tossing out
        % every second row/col.
        downStep  = 2;
        aveKernel = ones(2,2)./4; % Note that each entry has equal weight
        aveTruth  = convn(imgRef, aveKernel, 'same');
        aveRecon  = convn(imgDis, aveKernel, 'same');
        imgRef    = aveTruth(1:downStep:end,1:downStep:end);
        imgDis    = aveRecon(1:downStep:end,1:downStep:end);
    end
    
    % Define gradient convolution kernel. Use Prewitt Filters.
    % NOTE: Author's mention that any gradient calculation method would be
    % acceptable, so could vary this if desired.
    filtX = [1/3  0 -1/3; 1/3 0 -1/3;  1/3   0  -1/3];
    filtY = [1/3 1/3 1/3;  0  0   0 ; -1/3 -1/3 -1/3];                     
    
    % Get gradients of each image. 
    TruthGx = imfilter(imgRef, filtX, 'replicate'); 
    TruthGy = imfilter(imgRef, filtY, 'replicate');          

    ReconGx = imfilter(imgDis, filtX, 'replicate'); 
    ReconGy = imfilter(imgDis, filtY, 'replicate');     

    % Get gradient magnitude maps for each image.
    gradientMapTruth = sqrt(TruthGx.^2 + TruthGy.^2);
    gradientMapRecon = sqrt(ReconGx.^2 + ReconGy.^2);

    % Get GMS map from gradient magnitude images. In an ideal situation,
    % the images are "maximally similar" and the GMS map is 1 everywhere.
    GMS = ( 2 * gradientMapTruth .* gradientMapRecon + C ) ./ ...
        ( gradientMapTruth.^2 + gradientMapRecon.^2 + C );

    % Get 2D GMSD score.
    score = std(GMS(:));
end

end