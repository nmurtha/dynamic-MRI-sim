function res = fft3c(x)
% res = ifft3c(x)
% orthonormal centered 3D ifft
%

res = 1/sqrt(length(x(:)))*fftshift(fftn(ifftshift(x)));


