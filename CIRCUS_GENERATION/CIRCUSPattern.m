function [phase_table_y,phase_table_z,n_sampled,q] =...
    CIRCUSPattern(ny,nz,b,c,mode,N,seed,cal,display)
% Create a [ny,nz] CIRCUS undersampling pattern from multiple quanta.
% Has two modes depending on desired usage.
% Mode 1: Generates a specified number N of quanta
% Mode 2: Generates quanta until a desired undersampling factor N is reached.
% If seed < 0, uses a random quantum each time, otherwise keeps
% increasing the quantum number in sequence starting with seed.

% Initialize output variables.
n_sampled = 0;
phase_table_y = [];
phase_table_z = [];

% Set size of (optionally included) calibration region (i.e. rectangle in
% center of kspace). 
caly = 6;
calz = 6;

% Generate any quanta that are desired, according to the mode specified by
% the user. 
random = (seed < 0);
q = seed;
if N ~= 0
    condition = 1;
else
    condition = 0; % if no quanta desired, skip quanta generation.
end
while (condition)
  if random
    q = q - 1; % Get negative number to return random quanta
  else
    q = q + 1; % Increment quanta counter to next quanta number
  end
  % Get quanta and combine with existing pattern.
  [py_new,pz_new,n_sampled_new] = CIRCUSQuantum(ny,nz,b,c,q,0);
  [phase_table_y,phase_table_z,n_sampled] = CIRCUSCombine(phase_table_y,phase_table_z,n_sampled,py_new,pz_new,n_sampled_new);
  if mode==2
    condition = (n_sampled < ny*nz / N); % Check that appropriate fraction of kspace is sampled
  else
    condition = ((abs(q) - abs(seed)) < N); % Check that q has counted up far enough relative to initial quanta number
  end
end

% Get calibration region and combine if desired.
if cal
  [py_new,pz_new,n_sampled_new] = CIRCUSCalRegion(ny,nz,caly,calz);
  [phase_table_y,phase_table_z,n_sampled] = CIRCUSCombine(phase_table_y,phase_table_z,n_sampled,py_new,pz_new,n_sampled_new);
end

% Sort phase tables from center of kspace, outwards.
[phase_table_y,phase_table_z] = CIRCUSRadialSort(phase_table_y,phase_table_z,ny,nz);

% Actual undersampling factor may differ slightly from desired, calculate.
R = (ny*nz) / n_sampled;

% Display pattern if desired.
if display
    CIRCUS_pattern = zeros(ny,nz);
    for i = 1:n_sampled
      CIRCUS_pattern(phase_table_y(i),phase_table_z(i)) = 1;
    end
    figure; imagesc(CIRCUS_pattern)
    title(sprintf('CIRCUS, R = %.1f, b = %.1f, c = %.1f',R,b,c))
end

%disp(sprintf('size: %d X %d; points sampled: %d; acc: %.2f', ny, nz, n_sampled, R));

end