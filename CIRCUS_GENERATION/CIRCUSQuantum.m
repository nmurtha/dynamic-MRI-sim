function [phase_table_y,phase_table_z,n_sampled] = CIRCUSQuantum(ny,nz,b,c,q,display)
% Returns the phase table for just one quantum, specified by q argument
% Passing q=-1 uses a random value from 0 to 4*max(ny,nz). 
% NOTE: q in this code corresponds to m in the CIRCUS paper. The choice of
% q (i.e. m) determines which sampling point on the squares we're talking
% about, and hence the "quanta" we want. Consistently taking the first m
% index gives a distinct quanta than consistently taking the third m index,
% for example. 

% Get the largest dimension requested. CIRCUS pattern first generated for
% square NxN grid, then cut down to ny x nz grid later by selection of kz
% lines according to golden angle profile.
N = max(ny,nz);

% Define the Golden Ratio in memory.
r = (1+sqrt(5))/2; 

% Allocate memory for the phase tables. This is an over-allocation, and the
% vectors are trimmed after quanta generation.
phase_table_y = zeros(ny*nz,1);
phase_table_z = zeros(ny*nz,1);

% Create phase table for square NxN grid.
for j = 1:N/2   % The index of each nested square.
    J = 2*j;    % The length of a side on the current square.
    K = 4*J-4;  % The number of points on each square (-4 removes overlaps at corners when adding up side lengths)
   
    % If random quantum is desired, get quanta index. 
    if q<0
      q = randi(4*N,1); % note that q is m in the CIRCUS paper.
    end
    
    % Get indices of sampling points on periphery of current square,
    % applying radial (b) and spiral (c) shifts.
    i = floor( mod( (q+b*J)*(r^-1),1 ) * K ); % eqn 2 in paper. If b = 0, defaults to CIRCUS "base" (eqn 1)
    i = mod(i+ceil(J^c)-1,K);  % eqn 3 in paper, if c = 0 has no effect on i.
    
    % Get the coordinate of the point on the nested square perimeter, and
    % store the coordinates of the point in the phase table vectors.
    [y,z] = findPointOnSquare(i,j,N);
    phase_table_y(j) = y;
    phase_table_z(j) = z;
end

% Trim excess entries in the phase tables, taking only the entries for
% which sampling was performed.
n_sampled = N/2;
phase_table_y = phase_table_y(1:n_sampled);
phase_table_z = phase_table_z(1:n_sampled);

% If a rectangular pattern is desired, grab it.
if nz<ny
  [phase_table_y,phase_table_z,n_sampled] = CIRCUSRect(phase_table_y,phase_table_z,ny,nz);
end

% Get the actual undersampling factor.
R = (ny*nz) / n_sampled;

% Finally, we will use the phase tables to create the matrix which can be
% viewed if the user desires.
if display
    CIRCUS_pattern = zeros(ny,nz);
    for i = 1:n_sampled
      % This way captures duplicates too
      CIRCUS_pattern(phase_table_y(i),phase_table_z(i)) =...
          CIRCUS_pattern(phase_table_y(i),phase_table_z(i))+1;
    end
    imagesc(CIRCUS_pattern)
    title(sprintf('CIRCUS, Quantum %d, R = %.1f, b = %.1f, c = %.1f',q,R,b,c))
end

% Output quanta information to terminal.
%disp(sprintf('quantum %d, size: %d X %d; points sampled: %d; acc: %.2f', q, ny, nz, n_sampled, R))

end

function [y,z] = findPointOnSquare(k,j,N)
% k = sampled point index.
% j = index of current nested square (possible values are 1 to N/2).
% N = side length of largest square

% Get side length of nested square.
J = 2*j;        % J = side length of current nested square

% Get minimum and maximum coordinate values in the y and z direction for
% this square. NOTE: coordinate axis define to have origin at bottom left
% corner of the largest square, such that the nested square is centered at
% (N/2,N/2). The convention here is to have the squares biased toward
% upper-right quadrant (e.g. j = 1 centers squares bottom left corner at
% origin, top right corner in right quadrant), hence the +1 to minCoord. 
minCoord = N/2-j+1;  % Coordinates of bottom-left corner are (minCoord,minCoord). 
maxCoord = N/2+j;    % Coordinates of top-right corner are (maxCoord,maxCoord)

% Traverse the square starting in bottom left, go up to top left, across
% to top right, down to bottom right, back to bottom left. NOTE: the z-axis
% is vertical, the y-axis is horizontal.
if k<(J-1) % Index falls on left side of square, not including top left point
  y = minCoord;
  z = minCoord+k; % Only vertical coordinate increases as index "climbs" left side of square
elseif k<(2*J-1)  % Index falls on the top of the square, not including top right point
  y = minCoord+(k-(J-1)); 
  z = maxCoord;
elseif k<(3*J-2)
  y = maxCoord;
  z = maxCoord-(k-2*(J-1));
else
  y = maxCoord-(k-3*(J-1));
  z = minCoord;
end

end