function  [sensMaps]=sensitivy_map(arrSize,numCoils)

% SENSITIVITY MAP creates aritficial sensitivity maps
%
%   The map is created such as the Sum of Squares is 1
%
%                      sum(MapW.^2,3)=ones(arrSize)
%
%  INPUTS:
%               arrSize  =  Size of the maps, [rows,cols]  
%               numCoils =  Number of coils,  (even number)
%
%  OUTPUT
%               sensMaps =  Coil sensitivity maps, [rows,cols,numCoils]
%
% EXAMPLE
%              MapW=sensitivy_map([256,256],4);
%
% Santiago Aja-Fernandez
% Parallel Imaging Toolbox
% Feb. 2012
% www.lpi.tel.uva.es/~santi

% Check inputs.
if numel(arrSize)~=2
   error('Input Size must be 2D');
end
if rem(numCoils,2)==1
   error('Number of coils must be even');
end 

% Get size of the image array.
rows = arrSize(1);
cols = arrSize(2);

% Get vector of linearly increasing values, depending on number of rows and
% columns requested. The entries of these vectors will be used to form an
% array of relatively weighted values, where the row/col "closest" to the
% surface has max value and the row/col "furthest" from the surface has min
% value.
rowLinVec = 1:rows;
colLinVec = 1:cols;

% Use linearly increasing vectors specified above to make arrays where the
% row/col closest to one side has max weight, and then the weights decrease
% from there. The values are weighted such that the max weight is pi/2, and
% decreases linearly from there.
weightsRows = repmat(rowLinVec',[1,cols]);
weightsCols = repmat(colLinVec,[rows,1]);
weightsRows = pi/2 .* weightsRows./max(weightsRows(:)); % Scaled such that pi/2 is max value...
weightsCols = pi/2 .* weightsCols./max(weightsCols(:));

% Get the angular placement of each coil, units of radians. 
thetaStep = (2*pi/numCoils);
Theta     = 0:thetaStep:(2*pi-thetaStep); % Don't repeat 2pi, same as zero!
Theta     = Theta(1:end-1);

% Loop over the coils, looping over only half of the number of coils and
% using symmetric to get the opposite coil. 
for ii=1:numCoils/2
  if (Theta(ii)<=pi/2)
      N1=weightsRows.*cos(Theta(ii))+weightsCols.*sin(Theta(ii));
      N1=N1./max(N1(:)).*pi/2; % Normalize such that max value is pi/2.
      sensMaps(:,:,ii)=cos(N1); % Since max range N1 = 0:pi/2, range cos(N1) = 1:0
      sensMaps(:,:,numCoils./2+ii)=sin(N1); % Symmetric coil on opposite quadrant.
  else %>pi/2
      N1=(weightsRows).*abs(cos(Theta(ii)))+(pi/2-weightsCols).*sin(Theta(ii));
      N1=N1./max(N1(:)).*pi/2; 
      sensMaps(:,:,ii)=sin(N1);
      sensMaps(:,:,numCoils./2+ii)=cos(N1);      
  end
end

sensMaps=sensMaps./sqrt(numCoils/2);