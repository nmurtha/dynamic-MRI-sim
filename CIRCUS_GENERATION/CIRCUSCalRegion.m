function [phase_table_y,phase_table_z,n_sampled] = CIRCUSCalRegion(ny,nz,caly,calz)
% Returns the phase table for a square calibration region of size caly X calz
% in the center of a matrix of size ny X nz

n_sampled = caly*calz;
phase_table_y = zeros(n_sampled,1);
phase_table_z = zeros(n_sampled,1);

c = 1;
for y=floor(ny/2-caly/2):floor(ny/2+caly/2)-1
  for z=floor(nz/2-calz/2):floor(nz/2+calz/2)-1
    phase_table_y(c) = y;
    phase_table_z(c) = z;
    c = c + 1;
  end
end

end