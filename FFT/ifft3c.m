function res = ifft3c(x)
% res = ifft3c(x)
% orthonormal centered 3D ifft
%

res = sqrt(length(x(:)))*ifftshift(ifftn(fftshift(x)));
