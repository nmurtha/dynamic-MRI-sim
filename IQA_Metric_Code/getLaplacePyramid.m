function [ pyramidImages ] = getLaplacePyramid( img, numScales, pyrWindow )

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Feb 21, 2017
%
%  PURPOSE: Gets Laplacian pyramid breakdown of the input image.
% 
%  INPUTS:
%   img         -> Image for which pyramid is sought. Can be 2D or 3D.
% 
%   numScales   -> The number of subsampled scales in the pyramid.
%
%   pyrWindow   -> Size of sides of filter to be used in making the pyramid 
%                  (e.g. pyrWindow = 2, assuming a square window).
% 
%  OUTPUTS:
%   pyramidImages -> Cell array containing images in the pyramid.
%
%  EXTERNALS:
%   mirrorEdgePad.m  -> Pad array with mirror reflections, treating image
%                       edge pixels as mirror surface (written by Nathan 
%                       Murtha).

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR EXECUTION   -----------------------
%  ------------------------------------------------------------------------

% Get 1D binomial kernel, following same method called by Wang through the
% sepspyr package. First generate the 1D kernel as per sepspyr's
% binomialFilter.m. Then get 2D or 3D extensions, assuming seperability.
kernel1D = [0.5 0.5]';
for n=1:pyrWindow-2
    kernel1D = conv([0.5 0.5]', kernel1D);
end
kernel1D = sqrt(2) .* kernel1D;

% Get 2D or 3D kernel, depending on dimensions of the image.
dims = numel(size(img));
if (dims < 3)
    % 2D kernel. 
    kernel = ones(pyrWindow,pyrWindow);
    for i = 1:pyrWindow
        kernel(i,:) = kernel(i,:) * kernel1D(i);
        kernel(:,i) = kernel(:,i) * kernel1D(i);
    end
else 
    % 3D kernel
    kernel = ones(pyrWindow,pyrWindow,pyrWindow);
    for i = 1:pyrWindow
        kernel(i,:,:) = kernel(i,:,:) * kernel1D(i);
        kernel(:,i,:) = kernel(:,i,:) * kernel1D(i);
        kernel(:,:,i) = kernel(:,:,i) * kernel1D(i);
    end
end

% Check that image size is okay for the number of scales requested. 
sizeCondition = min(size(kernel)) <= min(size(img)) * 2^(-(numScales-1));
assert(sizeCondition,'Cant calculate %d scales with image of smallest dimension %d and kernel size %d.',...
    numScales,min(size(img)),min(size(kernel)));

%% ------------------------------------------------------------------------
%  -------------------------   MAKE THE PYRAMID   -------------------------
%  ------------------------------------------------------------------------

% Loop over the higher levels in the pyramid, calculating the bandpass
% filtered entries for all but highest level (which is simply a lowpass
% filtered image).
pyramidImages = cell(1,numScales);
imgCurrent = img;
for iLevel = 1:numScales-1
    % Get lowpass image.
    padImage    = mirrorEdgePad(imgCurrent,kernel); % Get mirror padded array, such that convn later can be used with 'valid' flag and have same padding as used in sepspyr package
    imLowReduce = convn(padImage,kernel,'valid'); % convolve and keep center portion    
    imLowReduce = imLowReduce(1:2:end,1:2:end,1:2:end); % downsample 2x
    
    % Get upsampled version of lowpass image, for bandpass image
    % calculation in the next step.
    imExpand = zeros(size(imgCurrent)); 
    imExpand(1:2:end,1:2:end,1:2:end,1:2:end) = imLowReduce; % Spread pixels of lowpass image into expanded array
    padExpand = mirrorEdgePad(imExpand,kernel); % Pad expanded array with mirror padding
    imExpand = convn(padExpand,kernel,'valid'); % convolve and keep center portion.
    
    % Get bandpass image.
    imBandpass = imgCurrent - imExpand;
    
    % Save current bandpass image into cell array for pyramid.
    pyramidImages{iLevel} = imBandpass;
    
    % Update current image.
    imgCurrent = imLowReduce;
end

% Assign final level entry, a lowpass image, to pyramid.
pyramidImages{numScales} = imgCurrent;

end

