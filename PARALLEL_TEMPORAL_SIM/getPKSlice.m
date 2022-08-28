function [ outerRings, innerDisks ] = getPKSlice(nr,nc)
%getPKSlice Generate embedded PK features.
%
%   Generates a 2D slice that can be put into a 3D phantom array. Contains
%   embedded disk features, to simulate malignant tissue inside healthy
%   tissue. 
%
%   The number of disk features per row and column of features can be
%   manually changed with the variables ncirc_r and ncirc_c. Additionally,
%   the size of the inner disks relative to the outer rings can be changed
%   with the variables Pmax and Pmin. 
%
%   The outer rings are all combined into one feature, since they're meant
%   to simulate healthy tissue and will use the same PK parameters. The
%   inner disks are groupe by "row" of features, such that one row of
%   disks decreasing in size is one indepedendant feature. This will allow
%   for seperate PK parameters going down the "feature rows" for comparison
%   in a simulation.
% 
%    INPUTS
%     nr          -> Number of rows in the output image.
%     nc          -> Number of cols in the output image.
%
%    OUTPUTS
%     outerRings  -> Outer ring feature. Has holes where inner disks will
%                    sit to avoid additive overlap. Measures (nr x nc).
%                    Output is a mask of 0's and 1's.
%     innerDisks  -> Inner disk features. Fit into holes in outer ring 
%                    feature. Measures (nr x nc x ncirc_r), where the 3rd
%                    dimension indicates an independant feature. Output is
%                    a mask of 0's and 1's.
%
%    EXTERNALS
%     [none]

% Set mask value of features (1 is the appropriate value).
val  = 1; 

% Get coordinates of slice on meshgrid.
[C, R] = meshgrid(1:nc, 1:nr);

% Set number of circles in columns and rows.
ncirc_r = 4; % circles down rows
ncirc_c = 4; % circles down cols

% Set increment stepping between circle centers in each dimension. Note:
% the denominator has a +1 to account for the leftover space when counting
% the equal spaced increments from left to right. Eg. for ncirc_r = 3 we
% get 4 equal gaps in the array plane.
%             |     d     d     d     |
%             |----->----->----->---->|
inc_row = nr / (ncirc_r + 1); % spacing down rows
inc_col = nc / (ncirc_c + 1); % spacing down cols

% Set radius of circle (or ellipse) in each direction. If nr == nc, the 
% shapes will be circles. Note that the increments calculated above are 
% equally spaced between circles, so they can be used to define the radius 
% of each circle. A radius equal to the increments makes the circles just 
% touch, so we must reduce the value slightly to get an appropriate radius 
% for distinct circles.
rad_outer_r = 0.8.* (0.5.*inc_row); % rows spanned
rad_outer_c = 0.8.* (0.5.*inc_col); % columns spanned

% Get radius of inner disks, as a percentage of outer ring radius. A linear
% relationship between index number and percentage of outer radius is
% calculated for the inner disk radius.
Pmax = 0.60;
Pmin = 0.20;
for i = 1:ncirc_r % Get percentages down the rows
   P = (i-1) * ((Pmin-Pmax)/(ncirc_r-1)) + Pmax;
   rad_inner_r(i) = P * rad_outer_r; 
end
for i = 1:ncirc_c % Get percentages down the cols
   P = (i-1) * ((Pmin-Pmax)/(ncirc_c-1)) + Pmax;
   rad_inner_c(i) = P * rad_outer_c; 
end

% Allocate space for outer rings feature and inner rings features. Note
% that there will be as many inner disk feature objects as there are rings
% going down the rows, since each "row" of features will have one common
% set of PK parameters.
outerRings = zeros(nr,nc);
innerDisks = zeros(nr,nc,ncirc_r); % Last index = feature identifier

% Loop down rows and across columns, making outer rings and inner disk
% features.
for r = 1:ncirc_r 
    
    % Allocate array for temp storage of inner disk feature.
    temp = zeros(nr,nc);
    
	for c = 1:ncirc_c
		
        % Get center coordinates for each feature.
        rctr = round(r*inc_row);
		cctr = round(c*inc_col);
        
        % Get indices marking where the features sit.
		inds_outerRing = ((C-cctr)/rad_outer_c).^2 + ((R-rctr)/rad_outer_r).^2 <= 1;
        inds_innerDisk = ((C-cctr)/rad_inner_c(c)).^2 + ((R-rctr)/rad_inner_r(c)).^2 <= 1;
        
        % Use indices to fill feature masks.
		outerRings(inds_outerRing) = val;
        temp(inds_innerDisk) = val;
    end
    
    % Assign one inner disk feature mask. Each "row" of inner disks
    % corresponds to a seperate feature item, for changing the PK
    % parameters.
    innerDisks(:,:,r) = temp;
    
    % Cut hole into outer rings, corresponding to the inner disks.
    outerRings = outerRings - innerDisks(:,:,r);
    
end

end

