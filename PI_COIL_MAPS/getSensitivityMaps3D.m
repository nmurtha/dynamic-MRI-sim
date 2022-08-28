function  [sensMaps, indsPlane] = getSensitivityMaps3D(arrSize,numCoils)

% Get dimensions of the image array, store in vectors for generation of a
% meshgrid.
rows   = 1:arrSize(1);
cols   = 1:arrSize(2);
slices = 1:arrSize(3);

% if user asked for one coil, return uniform sensitivity map, ideal case.
if numCoils == 1
    sensMaps  = ones(arrSize);
    indsPlane = ones(slices(end),cols(end));
    return;
end

% Ensure that the number of coils is even. 
assert(mod(numCoils,2) == 0, 'Number of coils must be even.');

% Generate meshgrid of values based on the dimensions.
[COLS, ROWS, SLICES] = meshgrid(rows,cols,slices);

% Get period of the sinusoid weighting function. Note that with period
% equal to 4 times the extent of the image array, will have a zero-crossing
% of the sinusoid at the edge of image array.
periodRow = 4*rows(end);
periodCol = 4*cols(end);
periodSlc = 4*slices(end);

% Determine coordinates at which "coils" will be placed (assuming point
% source for coils). Coils will be placed as if in some kind of ideal
% cardiac coil, in the superior and inferior planes. Equal number of coils
% in the inferior and superior planes. For simplicity, when the number of
% coils per plane is even, two strips of coils (with equal number of coils)
% are placed; when the number of coils per plane is odd, one coil is placed
% at the center of the plane and the remaining are set in two strips. For
% the cases where there are 1,2,3 coils per plane, the coils are placed
% along the central axis. Some examples (note that ASCII drawings suck):
%      ________________          ________________       ________________
%      |              |          |              |       |              |
%      |              |          |   *      *   |       |   *      *   |
%      |              |          |              |       |       *      |
%      |   *   *   *  |          |              |       |              |
%      |              |          |   *      *   |       |   *      *   |
%      |______________|          |______________|       |______________|
pointsPerPlane = numCoils / 2; % coils place symmetrically above and below object...
if pointsPerPlane <= 3 % If there are 1, 2 or 3 coils per plane.
    % Take middle column as axis for point placement.
    colInds = floor(cols(end)/2) + 1;
    
    % Get the spacing between points. Note that the +1 accounts for the
    % fact that there are 1 more gaps between points than there are actual
    % points:
    %                  |---*---*---*---|   e.g. 3 points but 4 gaps
    sliceStep = slices(end) / (pointsPerPlane + 1);
    
    % Get slice indices by counting the number of steps and incorporating
    % step length between points.
    sliceInds = floor(sliceStep .* [1:pointsPerPlane]);
    
    % Set index of point placements.
    indsPlane = zeros(slices(end),cols(end));
    indsPlane(sliceInds, colInds) = 1;
elseif mod(pointsPerPlane,2) == 1 && pointsPerPlane > 3 % 5, 7, ... coils per plane
    % Mark the middle column.
    midColInd = floor(cols(end)/2)+1;
    
    % Find col ind in the middle of the first "third-section".
    firstColInd = floor(cols(end)/4)+1;
    
    % Find col ind in the middle of the third "third-section".
    thirdColInd = floor(3*cols(end)/4)+1;
    
    % Set the slice direction indices.
    colInds = [firstColInd, midColInd, thirdColInd];
    
    % Get number of points to be placed in the "top third" and "bottom
    % third". We've already decided to place one point in the center column,
    % so have to place (nCoils/2)-1 coils in the remaining area, hence 
    % ((nCoils/2)-1)/2 points per "third".
    numberPointsPerThird = (pointsPerPlane - 1)/2;
    
    % Get "slice step" in each of the "first-section" and "third-section".       
    sliceStep = slices(end) / (numberPointsPerThird+1); 
    
    % Use "slice step" to get slice indices.
    sliceInds = floor(sliceStep .* [1:numberPointsPerThird]);
    
    % Set the indices for each point.
    indsPlane = zeros(slices(end),cols(end));
    indsPlane(sliceInds,colInds(1)) = 1;
    indsPlane(floor(mean([sliceInds(1),sliceInds(end)])) + 1, colInds(2)) = 1;
    indsPlane(sliceInds,colInds(3)) = 1;    
elseif mod(pointsPerPlane,2)==0 && pointsPerPlane > 3
    % Find col ind in the middle of the first "third-section".
    firstColInd = floor(cols(end)/4)+1;
    
    % Find slice ind in the middle of the third "third-section".
    thirdColInd = floor(3*cols(end)/4)+1;
    
    % Get the point place "steps" in each direction. 
    sliceStep = slices(end) / ((pointsPerPlane/2) + 1); 
    
    % Get the indices in each dimensions.
    sliceInds   = floor(sliceStep .* [1:pointsPerPlane/2]);
    colInds = [firstColInd thirdColInd];
    
    % Set indices.
    indsPlane = zeros(slices(end),cols(end));
    indsPlane(sliceInds,colInds) = 1;    
end

% Find the slices and columns for which points are placed. These are
% indicated by entries equal to 1 in indsPlane.
[slicePlacements,colPlacements] = find(indsPlane);

% Get coil sensitivities. Note that the sensitivity profiles of each
% symmetrically placed set of coils in the sup-inf direction is assumed
% equal, so only need to loop over *half* the coils and then use symmetry
% to get the others.
sensMaps = zeros(rows(end),cols(end),slices(end),numCoils);

for iCoil = 1:numCoils/2
    % Get current shift in the column and slice direction.
    colShift   = colPlacements(iCoil);
    sliceShift = slicePlacements(iCoil);
    
    % Get sensitivity profile for current coil in the superior plane. Will
    % use symmetry to get the sensitivity profile in the inferior plane
    % later. Assuming a cosine type sensitivity, relatively flat response
    % near each coil and increasing dropoff with increasing distance.
    tempMap = cos(2*pi/periodRow .* (ROWS)) .*...
        cos(2*pi/periodCol .* (COLS - colShift)) .*...
        cos(2*pi/periodSlc .* (SLICES-sliceShift));
    sensMaps(:,:,:,iCoil) = tempMap;% ./ max(tempMap(:)); % Normalize to have max weight equal to 1.
end

% Pull first numCoils/2 maps out in temporary array, flip the maps upside
% down, and assign the flipped maps to the inferior plane coils.
tempMaps = sensMaps(:,:,:,1:numCoils/2);
sensMaps(:,:,:,numCoils/2+1:end) = flipud(tempMaps);

% Normalize so that sum of squares across each pixel in the maps is 1 (i.e.
% the summed squares of sensitivities across all coils for a SINGLE PIXEL
% in the static image is 1, NOT saying that sum of squared weights in a
% single coil map is 1). 
squaresMaps = sensMaps.^2;
squaresMaps = squaresMaps ./ sum(squaresMaps,4);
sensMaps    = sqrt(squaresMaps);


end
