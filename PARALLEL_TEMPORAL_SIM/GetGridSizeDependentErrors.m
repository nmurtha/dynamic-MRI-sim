function [ gridErrors ] = GetGridSizeDependentErrors(... 
    errorMap, maskPinGrid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Set columns and rows that divide the phantom grids from each other.
% THESE PARAMETERS ARE SET FOR A PHANTOM SIZE OF 196x196, MUST BE MANUALLY
% SET.
colLeftSide            = 11;
colRightSide           = 185;
colTopRowFirstQuarter  = 46;
colTopRowSecondQuarter = 84;
colTopRowThirdQuarter  = 126;
colMidRowFirstThird    = 64;
colMidRowSecondThird   = 119;
colBotRowFirstHalf     = 94;

rowTop           = 11;
rowBottom        = 185;
rowFirstTier     = 50;
rowSecondTier     = 111;

% Section off portions of the mask for individual pin grids. 
gridMasks = cell(1,9);

gridMasks{1} = zeros(size(maskPinGrid));
gridMasks{1}(rowTop:rowFirstTier,colLeftSide:colTopRowFirstQuarter) = ...
    maskPinGrid(rowTop:rowFirstTier,colLeftSide:colTopRowFirstQuarter);

gridMasks{2} = zeros(size(maskPinGrid));
gridMasks{2}(rowTop:rowFirstTier,colTopRowFirstQuarter:colTopRowSecondQuarter) = ...
    maskPinGrid(rowTop:rowFirstTier,colTopRowFirstQuarter:colTopRowSecondQuarter);

gridMasks{3} = zeros(size(maskPinGrid));
gridMasks{3}(rowTop:rowFirstTier,colTopRowSecondQuarter:colTopRowThirdQuarter) = ...
    maskPinGrid(rowTop:rowFirstTier,colTopRowSecondQuarter:colTopRowThirdQuarter);

gridMasks{4} = zeros(size(maskPinGrid));
gridMasks{4}(rowTop:rowFirstTier,colTopRowThirdQuarter:colRightSide) = ...
    maskPinGrid(rowTop:rowFirstTier,colTopRowThirdQuarter:colRightSide);

gridMasks{5} = zeros(size(maskPinGrid));
gridMasks{5}(rowFirstTier:rowSecondTier,colLeftSide:colMidRowFirstThird) = ...
    maskPinGrid(rowFirstTier:rowSecondTier,colLeftSide:colMidRowFirstThird);

gridMasks{6} = zeros(size(maskPinGrid));
gridMasks{6}(rowFirstTier:rowSecondTier,colMidRowFirstThird:colMidRowSecondThird) = ...
    maskPinGrid(rowFirstTier:rowSecondTier,colMidRowFirstThird:colMidRowSecondThird);

gridMasks{7} = zeros(size(maskPinGrid));
gridMasks{7}(rowFirstTier:rowSecondTier,colMidRowSecondThird:colRightSide) = ...
    maskPinGrid(rowFirstTier:rowSecondTier,colMidRowSecondThird:colRightSide);

gridMasks{8} = zeros(size(maskPinGrid));
gridMasks{8}(rowSecondTier:rowBottom,colLeftSide:colBotRowFirstHalf) = ...
    maskPinGrid(rowSecondTier:rowBottom,colLeftSide:colBotRowFirstHalf);

gridMasks{9} = zeros(size(maskPinGrid));
gridMasks{9}(rowSecondTier:rowBottom,colBotRowFirstHalf:colRightSide) = ...
    maskPinGrid(rowSecondTier:rowBottom,colBotRowFirstHalf:colRightSide);

% Get masked errors by applying individual masks, calculated above, to the
% error map.
gridErrors = cell(1,9);

for i = 1:9
    gridErrors{i} = gridMasks{i} .* errorMap;
end

end