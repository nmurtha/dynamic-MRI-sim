function [overall_mssim, cs_map] = msssim(img1, img2, level, L, K, win, weight, method)

% Multi-scale Structural Similarity Index (MS-SSIM)
% Z. Wang, E. P. Simoncelli and A. C. Bovik, "Multi-scale structural similarity
% for image quality assessment," Invited Paper, IEEE Asilomar Conference on
% Signals, Systems and Computers, Nov. 2003
%
% INPUTS
%  img1, img2 -> images being compared.
%  level      -> Number of scales (levels) in the pyramid breakdown.
%  L          -> max pixel value in image.
%  K          -> K parameters used in stability constant calculation.
%  win        -> Kernel used in convolution for calculating maps.
%  weight     -> weights for each scale (up to 5). 
%  method     -> 'product' or 'sum', for weighted product or sum.
%
% OUTPUTS  
%  overall_mssim -> the final msssim score.
%
% EXTERNALS:
%  ssim_index_new.m

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR EXECUTION   -----------------------
%  ------------------------------------------------------------------------

% Check existance of input variables, and see if they violate any
% restrictions.
if (nargin < 2 || nargin > 8)
   disp('Improper number of inputs')
   overall_mssim = -Inf;
   return;
end

if (~exist('K','var'))
   K = [0.01 0.03];
end

if (~exist('win','var'))
   win = fspecial('gaussian', 11, 1.5);
end

if (~exist('level','var'))
   level = 5;
end

if (~exist('L','var'))
   L = 255;
end

if (~exist('weight','var'))
   weight = [0.0448 0.2856 0.3001 0.2363 0.1333];
end

if (~exist('method','var'))
   method = 'product';
end

if (size(img1) ~= size(img2))
   disp('Images must be the same size')
   overall_mssim = -Inf;
   return;
end

% Ensure image is large enough.
[M, N] = size(img1);
if ((M < 11) || (N < 11))
   overall_mssim = -Inf;
   return
end

if (length(K) ~= 2)
   overall_mssim = -Inf;
   return;
end

if (K(1) < 0 || K(2) < 0)
   overall_mssim = -Inf;
   return;
end

% Ensure that the kernel is appropriately sized and shaped.
[H, W] = size(win);
if ((H*W)<4 || (H>M) || (W>N))
   overall_mssim = -Inf;
   return;
end
   
if (level < 1)
   overall_mssim = -Inf;
   return
end

% Ensure that even at the smallest scale image, the kernel won't be too 
% large.
min_img_width = min(M, N)/(2^(level-1));
max_win_width = max(H, W);
if (min_img_width < max_win_width)
   overall_mssim = -Inf;
   return;
end

% Get weights for each scale, stopping at a max of 5 scales. Sum of weights
% is normalized to unity (only the relative weightings are important).
numScales = min(5,level);
weight = weight(1:numScales);
weight = weight./sum(weight);

if (length(weight) ~= level || sum(weight) == 0)
   overall_mssim = -Inf;
   return;
end

if (method ~= 'wtd_sum' & method ~= 'product')
   overall_mssim = -Inf;
   return;
end

%% ------------------------------------------------------------------------
%  ----------------------   CALCULATE THE MS-SSIM   -----------------------
%  ------------------------------------------------------------------------

% Get image quality maps at each of the pyramid levels.
downsample_filter = ones(2)./4;
im1 = double(img1);
im2 = double(img2);
mssim_array = zeros(1,level);
mcs_array = zeros(1:level);
for l = 1:level
   [mssim_array(l), ~, mcs_array(l), cs_map{l}] = ssim_index_new(im1, im2, K, win, L);

   filtered_im1 = imfilter(im1, downsample_filter, 'symmetric', 'same');
   filtered_im2 = imfilter(im2, downsample_filter, 'symmetric', 'same');
   
   clear im1 im2
   im1 = filtered_im1(1:2:end, 1:2:end);
   im2 = filtered_im2(1:2:end, 1:2:end);
end

% Pool level scores according to method. First scales only consider
% contrast and structure, while last layer includes luminance contribution.
if (method == 'product')
   overall_mssim = prod(mcs_array(1:level-1).^weight(1:level-1)) * (mssim_array(level).^weight(level));
else
   weight = weight./sum(weight);
   overall_mssim = sum(mcs_array(1:level-1).*weight(1:level-1)) + mssim_array(level).*weight(level);
end
