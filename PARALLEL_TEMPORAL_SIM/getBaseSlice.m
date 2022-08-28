function [ base ] = getBaseSlice( nr,nc )
%getBaseSlice Generate base slice for phantom.
%
%   Generates a 2D slice that can be put into a 3D phantom array. Contains
%   a background slice with no features. Background is surrounded by a
%   margin of "air".
%
%    INPUTS
%     nr          -> Number of rows in the output image.
%     nc          -> Number of cols in the output image.
%
%    OUTPUTS
%     base        -> Base slice. Measures (nr x nc). Output is a
%                    mask of 0's and 1's.
%
%    EXTERNALS
%     [none]

% Fill slice with ones initially
base = ones(nr,nc);

% Get margin width for "air" outside phantom
margin_c = round(0.05*nc);
margin_r = round(0.05*nr);

% Set margins equal to 0
base(1:margin_r,1:end) = 0;
base(nr-margin_r:nr,1:end) = 0;
base(1:end,1:margin_c) = 0;
base(1:end,nc-margin_c:nc) = 0;


end

