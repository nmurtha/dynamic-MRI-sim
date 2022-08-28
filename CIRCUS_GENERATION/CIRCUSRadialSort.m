function [phase_table_y,phase_table_z] = CIRCUSRadialSort(phase_table_y,phase_table_z,ny,nz)
% Sort all points in phase table by their distance from the k-space center
% Not as C-like as other CIRCUS routines since we don't need to replicate
% this in EPIC PSD, there's already a function that will do this for us

%dist = sqrt((phase_table_y(:)-ny/2).^2 + (phase_table_z(:)-nz/2).^2);
dist = sqrt((phase_table_y(:)-ny/2).^2 + ((phase_table_z(:)-nz/2)*(ny/nz)).^2); % Nathan added to more equally spread points for narrow nz
[~,indx] = sort(dist);
phase_table_y = phase_table_y(indx);
phase_table_z = phase_table_z(indx);

end
