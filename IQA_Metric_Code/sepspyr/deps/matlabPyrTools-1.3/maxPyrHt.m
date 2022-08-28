% HEIGHT = maxPyrHt(IMSIZE, FILTSIZE)
%
% Compute maximum pyramid height for given image and filter sizes.
% Specifically: the number of corrDn operations that can be sequentially
% performed when subsampling by a factor of 2.

% Eero Simoncelli, 6/96. Edited by Nathan Murtha to support 3D images, Feb
% 2017.

function height = maxPyrHt(imsz, filtsz)

% Turn vectors into column vectors. Expecting vectors with first non-one
% entry in first index (e.g. [4,1]);
imsz   = imsz(:);
filtsz = filtsz(:);

if any(imsz == 1)         % 1D image, 1D filter
  imsz   = prod(imsz);
  filtsz = prod(filtsz);  % Filter will be reshaped
  
elseif any(filtsz == 1) && numel(imsz) == 2   % 2D image, 1D filter
  filtsz = [filtsz(filtsz ~= 1); filtsz(filtsz ~= 1)];
  
elseif any(filtsz == 1) && numel(imsz) == 3   % 3D image, 1D filter
  filtsz = [filtsz(filtsz ~= 1); filtsz(filtsz ~= 1); filtsz(filtsz ~= 1)];
end

% Recursively call maxPyrHt, counting number of times that the image is
% bigger than the filter after undersampling.
if any(imsz < filtsz)
  height = 0;
else
  height = 1 + maxPyrHt( floor(imsz/2), filtsz ); % recursive call
end

end
