function [ temp ] = ViewCIRCUS( dim1, dim2, phase1, phase2 )
%ViewCIRCUS Displays the CIRCUS pattern contained in the provided phase
%tables.

temp = zeros(dim1,dim2);
for i = 1:numel(phase1)
    temp(phase1(i),phase2(i)) = 1;
end
%imtool3D(temp)
imagesc(temp)
end

