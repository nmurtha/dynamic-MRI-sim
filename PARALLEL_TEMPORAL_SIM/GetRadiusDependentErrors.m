function [ outerCylErrors, innerCylErrors ] = GetRadiusDependentErrors(...
    errorMap, maskOuter, maskRow1, maskRow2, maskRow3, maskRow4)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Set columns and rows that divide the phantom features into a 4x4 grid.
% THESE PARAMETERS ARE SET FOR A PHANTOM SIZE OF 196x196, MUST BE MANUALLY
% SET.
colLeftSide      = 11;
colFirstQuarter  = 59;
colSecondQuarter = 98;
colThirdQuarter  = 138;
colRightSide     = 185;

rowTop           = 11;
rowFirstQuarter  = 59;
rowSecondQuarter = 98;
rowThirdQuarter  = 138;
rowBottom        = 185;

% Summarize these values in vectors for looping.
rowMarkers = [rowTop rowFirstQuarter rowSecondQuarter rowThirdQuarter rowBottom];
colMarkers = [colLeftSide colFirstQuarter colSecondQuarter colThirdQuarter colRightSide];

% Get composite of all inner row features.
maskInner = maskRow1 + maskRow2 + maskRow3 + maskRow4;

% For each of the sectors between the indices set above, assign the mask
% values to individual masks in an array, one per cylinder.
outerCyls = cell(4,4);
innerCyls = cell(4,4);
for r = 1:numel(rowMarkers) - 1
    for c = 1:numel(colMarkers) - 1

        outerCyls{r,c} = zeros(size(maskOuter));
        outerCyls{r,c}(rowMarkers(r):rowMarkers(r+1),colMarkers(c):colMarkers(c+1)) =...
            maskOuter(rowMarkers(r):rowMarkers(r+1),colMarkers(c):colMarkers(c+1));

        innerCyls{r,c} = zeros(size(maskInner));
        innerCyls{r,c}(rowMarkers(r):rowMarkers(r+1),colMarkers(c):colMarkers(c+1)) =...
            maskInner(rowMarkers(r):rowMarkers(r+1),colMarkers(c):colMarkers(c+1));
    end
end

% Get masked errors by applying individual masks, calculated above, to the
% error maps.
outerCylErrors = cell(4,4);
innerCylErrors = cell(4,4);
for r = 1:numel(rowMarkers) - 1
    for c = 1:numel(colMarkers) - 1

        outerCylErrors{r,c} = outerCyls{r,c} .* errorMap;

        innerCylErrors{r,c} = innerCyls{r,c} .* errorMap;
    end
end


end