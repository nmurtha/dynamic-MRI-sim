function [ phantBase, phantFeatures ] = Get_CS_Phantom( numRows, numCols,...
                                        samplingTimes, signalRange )
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Jan 22, 2017
%   
%  PURPOSE: Generates a 3D numerical phantom for use in assessing CS 
%           performance in quantitative parameter recovery. Contains 32 
%           slices:
%            Slice 1-2   = empty phantom base
%            Slice 3-5   = dynamic embedded cylinders (exponential)
%            Slice 6-8   = empty phantom base
%            Slice 9-11  = dynamic embedded cylinders (exponential)
%            Slice 12-14 = empty phantom base
%            Slice 15-17 = dynamic embedded cylinders (sinusoidal)
%            Slice 18-20 = empty phantom base
%            Slice 21-23 = static pin grids
%            Slice 24-26 = empty phantom base
%            Slice 27-29 = linearly increasing pin grids
%            Slice 30-32 = empty phantom base
% 
%  INPUTS
%   numRows       -> Number of rows desired in simulated phantom image. A
%                    scalar input.
%
%   numCols       -> Number of columns desired in simulated phantom image.
%                    A scalar input. 
%
%   samplingTimes -> Times at which each frequency encode begins. A vector,
%                    in units of seconds.
%
%   signalRange   -> Maximum value, in arbitrary units, of the magnitude
%                    of feature signals. Assuming that 0 is the min.
%                    Defaults to 255 if value isn't specified by user.
%
%  OUTPUTS
%   phantBase     -> Static base of the phantom. A 3D array with arbitrary 
%                    units.
%
%   phantFeatures -> phantFeatures(i).shape contains the shape of the
%                    i-th phantom feature (3D array, unitless).
%                 -> phantFeatures(i).Ct contains signal intensity 
%                    evolution information for the i-th phantom feature 
%                    (vector, arbitrary units).
%                 -> phantFeatures(i).params contains supplied temporal
%                    parameters for the i-th phantom feature 
%                    (vector, units listed below). 
%                    Exp decays = [initial magnitude (arb units), 
%                                  decay constant (sec),
%                                  constant offset (arb units)].
%                    Sinusoids =  [amplitude (arb units), 
%                                  temporal period (sec), 
%                                  temporal offset (sec),
%                                  constant offset (arb units)].
%                    Linear = [normalized slope (arb unit / sec), 
%                              constant offset (arb units)].
%                    These parameters are manually set in the code below.
%                 
%  EXTERNALS
%   getBaseSlice.m
%   getPKslice.m
%   getPinGridSlice.m

%% ------------------------------------------------------------------------
%  --------------------   SET CONSTANTS FOR PHANTOM   ---------------------
%  ------------------------------------------------------------------------

% Set number of slices. DO NOT CHANGE THIS without manually altering the
% slices in the phantom.
numSlices = 32;

% Set signal intensity range to default if not specified.
if (~exist('signalRange','var'))
    signalRange = 255; % [arb units]
end

%% ------------------------------------------------------------------------
%  -------------------   GATHER PHANTOM SLICE MASKS   ---------------------
%  ------------------------------------------------------------------------
%
%  Phantom slice masks are generated here. Masks are arrays of 0's and
%  1's, where 1's are placed in areas where relevant features will be.
%          maskBaseSlice -> 2D array.
%          maskOuterCyl  -> 2D array.
%          maskInnerCyl  -> 3D array, 4 entries in third dimension.
%
%  These masks will be used to generate the phantFeatures.shape entries, 
%  to which phantFeatures.Ct and phantFeature.params entries will be 
%  assigned. 

% Obtain base slice mask.
maskBaseSlice = getBaseSlice(numRows,numCols);

% Obtain embedded cylinder feature masks.
[ maskOuterCyl, maskInnerCyl ] = getPKSlice(numRows,numCols);

% Obtain pin grid mask.
maskPinSlice = getPinGridSlice(numRows,numCols);


%% ------------------------------------------------------------------------
%  ------------   CREATE BACKGROUND PROFILE FOR EACH SLICE   --------------
%  ------------------------------------------------------------------------
%
%  For each slice with features, the background must have the feature mask
%  itself subtracted to ensure no additive overlap (i.e. we "carve holes" 
%  for the phantom features). 

% Make embedded cylinders slice background profile.
bkrdEmbedCyls = maskBaseSlice - maskOuterCyl;
for iInnerCyl = 1:length( maskInnerCyl(1,1,:) )
    bkrdEmbedCyls = bkrdEmbedCyls - maskInnerCyl(:,:,iInnerCyl);
end

% Make pin grid slice background profile.
bkrdPinFeats = maskBaseSlice - maskPinSlice;


%% ------------------------------------------------------------------------
%  -------------------   ASSEMBLE BASE OF PHANTOM   -----------------------
%  ------------------------------------------------------------------------
%
%  The phantom base is defined by concatenating the above background
%  profiles for each slice type appropriately, making a 3D array consisting
%  of 32 slices, each slice measuring numRows by numCols.
%
%  If the number of slices is to be changed, this is one of the places
%  where the code must manually be altered.

% Allocate array for static phantom base.
phantBase = zeros(numRows,numCols,numSlices); 

% Populate each slice of the phantom.
phantBase(:,:,1)  = maskBaseSlice;
phantBase(:,:,2)  = maskBaseSlice;

phantBase(:,:,3)  = bkrdEmbedCyls;
phantBase(:,:,4)  = bkrdEmbedCyls;
phantBase(:,:,5)  = bkrdEmbedCyls;

phantBase(:,:,6)  = maskBaseSlice;
phantBase(:,:,7)  = maskBaseSlice;
phantBase(:,:,8)  = maskBaseSlice;

phantBase(:,:,9)  = bkrdEmbedCyls;
phantBase(:,:,10) = bkrdEmbedCyls;
phantBase(:,:,11) = bkrdEmbedCyls;

phantBase(:,:,12) = maskBaseSlice;
phantBase(:,:,13) = maskBaseSlice;
phantBase(:,:,14) = maskBaseSlice;

phantBase(:,:,15) = bkrdEmbedCyls;
phantBase(:,:,16) = bkrdEmbedCyls;
phantBase(:,:,17) = bkrdEmbedCyls;

phantBase(:,:,18) = maskBaseSlice;
phantBase(:,:,19) = maskBaseSlice;
phantBase(:,:,20) = maskBaseSlice;

phantBase(:,:,21) = bkrdPinFeats;
phantBase(:,:,22) = bkrdPinFeats;
phantBase(:,:,23) = bkrdPinFeats;

phantBase(:,:,24) = maskBaseSlice;
phantBase(:,:,25) = maskBaseSlice;
phantBase(:,:,26) = maskBaseSlice;

phantBase(:,:,27) = bkrdPinFeats;
phantBase(:,:,28) = bkrdPinFeats;
phantBase(:,:,29) = bkrdPinFeats;

phantBase(:,:,30) = maskBaseSlice;
phantBase(:,:,31) = maskBaseSlice;
phantBase(:,:,32) = maskBaseSlice;


%% ------------------------------------------------------------------------
%  ---------------------   ASSIGN FEATURE SHAPES  -------------------------
%  ------------------------------------------------------------------------
%
%  The feature shape assignments are done here, using the feature maskes
%  created above. The feature shapes are extended to three dimensional 
%  arrays, where all entries are zero outside the feature zones (this is 
%  for convenience in adding the features to the base later, so that array 
%  dimensions agree). The feature shape arrays are then stored in a
%  structure array for output.

% Preallocate structure array for features. This has 3D arrays for the
% shapes, a vector for the contrast evolution, and a vector for the 
% temporal parameters used in making the contrast evolution. Only the first
% field, the 3D shape arrays, is populated in this section.
numFeats = 17;
phantFeatures(1:numFeats) = struct('shape',zeros(numRows,numCols,numSlices),...
    'Ct',zeros(1,length(samplingTimes)),...
    'params',[]);

% Features in slices 3 to 5 will be embedded cylinders.
for iSlice = 3:5
    phantFeatures(1).shape(:,:,iSlice) = maskOuterCyl;
    phantFeatures(2).shape(:,:,iSlice) = maskInnerCyl(:,:,1);
    phantFeatures(3).shape(:,:,iSlice) = maskInnerCyl(:,:,2);
    phantFeatures(4).shape(:,:,iSlice) = maskInnerCyl(:,:,3);
    phantFeatures(5).shape(:,:,iSlice) = maskInnerCyl(:,:,4);
end

% Features in slices 9 to 11 will be embedded cylinders.
for iSlice = 9:11
    phantFeatures(6).shape(:,:,iSlice)  = maskOuterCyl;
    phantFeatures(7).shape(:,:,iSlice)  = maskInnerCyl(:,:,1);
    phantFeatures(8).shape(:,:,iSlice)  = maskInnerCyl(:,:,2);
    phantFeatures(9).shape(:,:,iSlice)  = maskInnerCyl(:,:,3);
    phantFeatures(10).shape(:,:,iSlice) = maskInnerCyl(:,:,4);
end

% Features in slices 15 to 17 will be embedded cylinders.
for iSlice = 15:17
    phantFeatures(11).shape(:,:,iSlice) = maskOuterCyl;
    phantFeatures(12).shape(:,:,iSlice) = maskInnerCyl(:,:,1);
    phantFeatures(13).shape(:,:,iSlice) = maskInnerCyl(:,:,2);
    phantFeatures(14).shape(:,:,iSlice) = maskInnerCyl(:,:,3);
    phantFeatures(15).shape(:,:,iSlice) = maskInnerCyl(:,:,4);
end

% Features in slices 21 to 23 will be pin grids.
for iSlice = 21:23
    phantFeatures(16).shape(:,:,iSlice) = maskPinSlice;
end 

% Features in slice 27 to 29 will be pin grids.
for iSlice = 27:29
    phantFeatures(17).shape(:,:,iSlice) = maskPinSlice;
end


%% ------------------------------------------------------------------------
%  -----------------   ASSIGN FEATURE TEMPORAL DYNAMICS  ------------------
%  ------------------------------------------------------------------------
%
%  The parameters dictating contrast evolution in the features are set
%  here. Slices 3-5 and 9-11 contain exponentially decaying contrast
%  features (the first changing in decay rate and the second changing in
%  initial magnitude) while slices 15-17 contain sinusoidally varying
%  contrast of different frequencies. Slices 21-23 contain static pin
%  grids, while slices 27-29 contain linearly increasing pin grids.
%
%  Recall that the format for temporal parameters goes as:
%       Exp decays -> [initial magnitude (arb units), 
%                      decay constant (sec),
%                      constant offset (arb units)].
%       Sinusoids  -> [amplitude (arb units), 
%                      temporal period (sec), 
%                      temporal offset (sec),
%                      constant offset (arb units)].
%       Linear     -> [normalized slope (arb unit / sec), 
%                      constant offset (arb units)].

% Set base intensity of phantom to be half of the signal range (arbitrary
% choice).
phantBase = phantBase * signalRange/2;

% Features in slices 3 to 5 (exponentially decaying embedded cylinders).
% Decay constants change between rows of cylinder features in this slice.
% Follows the form f(t) = A*exp(-t/B) + C.
phantFeatures(1).params = [floor(0.9*signalRange), 100, floor(0.1*signalRange)]; % outer cylinders all have same dynamics.
phantFeatures(1).Ct = phantFeatures(1).params(1) * ...
    exp(-samplingTimes./phantFeatures(1).params(2)) + ...
    phantFeatures(1).params(3);

phantFeatures(2).params = [floor(0.9*signalRange), 40, floor(0.1*signalRange)]; % top inner cylinders decay "quickly"
phantFeatures(2).Ct = phantFeatures(2).params(1) * ...
    exp(-samplingTimes./phantFeatures(2).params(2)) + ...
    phantFeatures(2).params(3);

phantFeatures(3).params = [floor(0.9*signalRange), 80, floor(0.1*signalRange)]; % second row of inner cylinders
phantFeatures(3).Ct = phantFeatures(3).params(1) * ...
    exp(-samplingTimes./phantFeatures(3).params(2)) + ...
    phantFeatures(3).params(3);

phantFeatures(4).params = [floor(0.9*signalRange), 200, floor(0.1*signalRange)]; % third row of inner cylinders
phantFeatures(4).Ct = phantFeatures(4).params(1) * ...
    exp(-samplingTimes./phantFeatures(4).params(2)) + ...
    phantFeatures(4).params(3);

phantFeatures(5).params = [floor(0.9*signalRange), 300, floor(0.1*signalRange)]; % fourth row of inner cylinders
phantFeatures(5).Ct = phantFeatures(5).params(1) * ...
    exp(-samplingTimes./phantFeatures(5).params(2)) + ...
    phantFeatures(5).params(3);

% Features in slices 9 to 11 (exponentially decaying embedded cylinders).
% Initial signal intensities change between rows of cylinder features in 
% this slice. Follows the form f(t) = A*exp(-t/B) + C.
phantFeatures(6).params = [floor(0.9*signalRange), 100, floor(0.1*signalRange)]; % outer cylinders all have same dynamics.
phantFeatures(6).Ct = phantFeatures(6).params(1) * ...
    exp(-samplingTimes./phantFeatures(6).params(2)) + ...
    phantFeatures(6).params(3);

phantFeatures(7).params = [floor(0.85*signalRange), 100, floor(0.1*signalRange)]; % first row of inner cylinders
phantFeatures(7).Ct = phantFeatures(7).params(1) * ...
    exp(-samplingTimes./phantFeatures(7).params(2)) + ...
    phantFeatures(7).params(3);

phantFeatures(8).params = [floor(0.75*signalRange), 100, floor(0.1*signalRange)]; % second row of inner cylinders
phantFeatures(8).Ct = phantFeatures(8).params(1) * ...
    exp(-samplingTimes./phantFeatures(8).params(2)) + ...
    phantFeatures(8).params(3);

phantFeatures(9).params = [floor(0.65*signalRange), 100, floor(0.1*signalRange)]; % third row of inner cylinders
phantFeatures(9).Ct = phantFeatures(9).params(1) * ...
    exp(-samplingTimes./phantFeatures(9).params(2)) + ...
    phantFeatures(9).params(3);

phantFeatures(10).params = [floor(0.55*signalRange), 100, floor(0.1*signalRange)]; % fourth row of inner cylinders
phantFeatures(10).Ct = phantFeatures(10).params(1) * ...
    exp(-samplingTimes./phantFeatures(10).params(2)) + ...
    phantFeatures(10).params(3);

% Features in slices 15 to 17 (sinusoidally varying embedded cylinders).
% Temporal periods change between rows of cylinder features in this slice.
% Follows the form f(t) = Asin(w(t-t'))+D, where w = 2*pi/period.
phantFeatures(11).params = [0.4*signalRange,60,0,0.5*signalRange]; % outer cylinders all have same dynamics.
phantFeatures(11).Ct = phantFeatures(11).params(1) * ...
    sin(2*pi/phantFeatures(11).params(2)*(samplingTimes-phantFeatures(11).params(3))) + ...
     phantFeatures(11).params(4);

phantFeatures(12).params = [0.4*signalRange,20,0,0.5*signalRange]; % top row of inner cylinders
phantFeatures(12).Ct = phantFeatures(12).params(1) * ...
    sin(2*pi/phantFeatures(12).params(2)*(samplingTimes-phantFeatures(12).params(3))) + ...
     phantFeatures(12).params(4);

phantFeatures(13).params = [0.4*signalRange,40,0,0.5*signalRange]; % second row of inner cylinders
phantFeatures(13).Ct = phantFeatures(13).params(1) *...
    sin(2*pi/phantFeatures(13).params(2)*(samplingTimes-phantFeatures(13).params(3))) + ...
     phantFeatures(13).params(4);

phantFeatures(14).params = [0.4*signalRange,80,0,0.5*signalRange]; % third row of inner cylinders
phantFeatures(14).Ct = phantFeatures(14).params(1) *...
    sin(2*pi/phantFeatures(14).params(2)*(samplingTimes-phantFeatures(14).params(3))) + ...
     phantFeatures(14).params(4);

phantFeatures(15).params = [0.4*signalRange,100,0,0.5*signalRange]; % fourth row of inner cylinders
phantFeatures(15).Ct = phantFeatures(15).params(1) *...
    sin(2*pi/phantFeatures(15).params(2)*(samplingTimes-phantFeatures(15).params(3))) + ...
     phantFeatures(15).params(4);

% Features in slices 21 to 23 (static pin grids). Follows form
% f(t) = (t/t_end) * A + B, such that max value is given by A (due to
% normalizing by t/t_end) and B is a vertical offset.
phantFeatures(16).params = [0,0.75*signalRange]; % all pin grids have same dynamics.
phantFeatures(16).Ct  = (samplingTimes./samplingTimes(end)) .* ...
    phantFeatures(16).params(1) + phantFeatures(16).params(2);
 
% Features in slices 27 to 29 (linearly increasing pin grids). 
% Follows form f(t) = (t/t_end) * A + B, such that max value is given by A 
% (due to normalizing by t/t_end) and B is a vertical offset.
phantFeatures(17).params = [0.9*signalRange,0.1*signalRange]; % all pin grids have same dynamics.
phantFeatures(17).Ct = (samplingTimes./samplingTimes(end)) .* ...
    phantFeatures(17).params(1) + phantFeatures(17).params(2);

end

