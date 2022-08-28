function padImg = mirrorEdgePad(img,kernel)
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           April 4, 2017.
%
%  PURPOSE: Extend the reflect1 style padding used in the sepspyr package,
%           called by Wang in the IW-SSIM source code. 
%           Pads array "img" with a number of pixels equivalent to the
%           kernel overhang (i.e. number of kernel pixels that are "off
%           image" when convolving at image boundary). E.g. for a kernel
%           spanning 5 pixels, kernel overhang is 2. 
%           
%           Treats the edge pixels of the img array as if it was a mirror
%           surface. See figure below for example of this, where * means
%           padded element and | means original array.
%
%                           *  i  *  h  *  g  *  h  *  i  * 
%                           *  f  *  e  *  d  *  e  *  f  * 
%                           *  c  *  b  |  a  |  b  |  c  |
%                           *  f  *  e  |  d  |  e  |  f  |
%                           *  i  *  h  |  g  |  h  |  i  |  
%
%           In contrast, the "symmetric" padding provided by Mathworks'
%           imfilter treats the virtual edge of the array as a mirror
%           surface, such that the result is behaves like:
% 
%                           *  b  *  a  |  a  |  b  |  c  |                          
%
%
%           NOTE: See rconv2.m in the sepspyr package for reflect1 style
%                 padding used there.
% 
%  INPUTS:
%   img         -> Image to be padded (2D or 3D array).
%
%   kernel      -> The convolution kernel that is to be used with img (2D
%                  or 3D array).
% 
%  OUTPUTS:
%   padImg      -> Mirror-padded version of img.
%
%  EXTERNAL FUNCTIONS:
%    [NONE]

% Get dimensions of the image.
[imgRows, imgCols, imgSlices] = size(img);

% Get portion of kernel that "hangs" off the edge of the image during the
% convolution at edge pixels.
kernelHangRow   = floor(size(kernel,1)/2); % Assuming odd kernel
kernelHangCol   = floor(size(kernel,2)/2); 
if imgSlices > 1
    kernelHangSlice = floor(size(kernel,3)/2); 
end

% Pad array with symmetric padding in each dimension, using padding of size
% kernelHang+1, where the +1 serves as a buffer zone that will be deleted
% to leave a true "mirror" reflection. 
if imgSlices == 1
    padImg = padarray(img,[kernelHangRow+1 kernelHangCol+1],...
        'both','symmetric');
else
    padImg = padarray(img,[kernelHangRow+1 kernelHangCol+1 kernelHangSlice+1],...
        'both','symmetric');
end
    
% Cyclically swap columns that were added adjacent to edge of img towards
% outside of padImg. 
for iCol = kernelHangCol:-1:1
    % Cycle columns on left side of padImg
    padImg(:,[iCol iCol+1],:) = padImg(:,[iCol+1 iCol],:); 
    
    % Cycle columns on right side of padImg
    padImg(:,[imgCols+2*(kernelHangCol+1)-iCol imgCols+2*(kernelHangCol+1)-iCol+1],:) = ...
        padImg(:,[imgCols+2*(kernelHangCol+1)-iCol+1 imgCols+2*(kernelHangCol+1)-iCol],:); 
end

% Cyclically swap rows that were added adjacent to edge of img towards
% outside of padImg. 
for iRow = kernelHangRow:-1:1
    % Cycle rows on top side of padImg
    padImg([iRow iRow+1],:,:) = padImg([iRow+1 iRow],:,:); 
    
    % Cycle rows on bottom side of padImg
    padImg([imgRows+2*(kernelHangRow+1)-iRow imgRows+2*(kernelHangRow+1)-iRow+1],:,:) = ...
        padImg([imgRows+2*(kernelHangRow+1)-iRow+1 imgRows+2*(kernelHangRow+1)-iRow],:,:); 
end

% Cyclically swap slices that were added adjacent to edge of img towards
% outside of padImg. 
if imgSlices > 1
    for iSlice = kernelHangSlice:-1:1
        % Cycle slices on front side of padImg
        padImg(:,:,[iSlice iSlice+1]) = padImg(:,:,[iSlice+1 iSlice]); 

        % Cycle slices on back side of padImg
        padImg(:,:,[imgSlices+2*(kernelHangSlice+1)-iSlice imgSlices+2*(kernelHangSlice+1)-iSlice+1]) = ...
            padImg(:,:,[imgSlices+2*(kernelHangSlice+1)-iSlice+1 imgSlices+2*(kernelHangSlice+1)-iSlice]); 
    end
end

% Delete the surplus outside rows, columns and slices from the outside of 
% padImg, leaving behind a mirror reflection style padded array.
padImg(:,1,:)   = [];
padImg(:,end,:) = [];
padImg(1,:,:)   = [];
padImg(end,:,:) = [];
if imgSlices > 1
    padImg(:,:,1)   = [];
    padImg(:,:,end) = [];
end

end

