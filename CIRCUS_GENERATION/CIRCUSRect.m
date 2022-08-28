function [phase_table_y,phase_table_z,n_sampled] = CIRCUSRect(phase_table_y,phase_table_z,ny,nz);
% Convert phase table to a rectangular grid (if square nothing will change)
% Select lines according to the golden ratio profile to get right dimensions.

% If the desired size is a square, then don't do anything.
if (ny==nz)
  return; 
end

% The Golden Ratio
r = (1+sqrt(5))/2; 

% Get the number of sampled points.
n_sampled = size(phase_table_y,1);

% Get index of center points, allocate vector for z-lines that are kept,
% and save the center line. 
n_cntr = ny/2;
z_keep = zeros(nz,1);
z_keep(1) = n_cntr;

% Loop until nz lines have been kept.
n_keep = 1;
while (n_keep < nz)
  % Following eqn 1 in paper, select z-lines using golden ratio profile. 
  z_keep(n_keep+1) = floor(mod(n_cntr*r^-1,1)*ny)+1;
  
  % Increment the index used in eqn for golden ratio profile.
  n_cntr = n_cntr + 1;
  
  % Check for duplicate entries. If none are found, proceed.
  dupe = 0;
  for z_cntr = 0:n_keep-1
    if z_keep(n_keep+1) == z_keep(z_cntr+1)
      dupe=1;
      break;
    end
  end
  if dupe==0
    n_keep = n_keep + 1; 
  end
end

% Set up new Z indices
z_keep_sorted = sort(z_keep);
new_z_index = zeros(ny,1);
for z_cntr = 1:n_keep
  new_z_index(z_keep_sorted(z_cntr)) = z_cntr;
end

%Iterate through phase table and remove anything with a Z coordinate that
%is not in the list of keepers - not efficient but should work
n_sampled_new = 0;
phase_table_ynew = zeros(size(phase_table_y));
phase_table_znew = zeros(size(phase_table_z));
for cntr=1:n_sampled
  keepit = 0;
  for z_cntr = 1:nz
    if phase_table_z(cntr) == z_keep(z_cntr)
      keepit = 1;
      break;
    end
  end
  if keepit
    n_sampled_new = n_sampled_new + 1;
    phase_table_ynew(n_sampled_new) = phase_table_y(cntr);
    phase_table_znew(n_sampled_new) = new_z_index(phase_table_z(cntr));
  end
end
n_sampled = n_sampled_new;
phase_table_y = phase_table_ynew(1:n_sampled);
phase_table_z = phase_table_znew(1:n_sampled);

end
