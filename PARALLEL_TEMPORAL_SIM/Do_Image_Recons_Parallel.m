function [ BARTrecons ] = Do_Image_Recons_Parallel( interleavedKspace,...
    BARTregularization )

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Jan 26, 2017
%
%  PURPOSE: Reconstructs images from k-space data using the BART toolkit 
%           for compressed sensing, as well as the conventional IFT 
%           reconstruction.
%
% INPUTS:
%    interleavedKspace  -> 4D array of interleaved kspace data. Fourth
%                          dimension represents the volume. An ndSparse 
%                          array, to save RAM.
%
%    BARTregularization -> Compressed sensing regularization for BART.
%                          E.g. '-R T:7:0:0.03'.
%                          NOTE: See pics.c in BART src directory for
%                                explanation of regularization. For
%                                T:A:B:C, T is the transform while "A is 
%                                transform flags, B is joint threshold 
%                                flags and C is regularization value".
%
% OUTPUTS:
%    BARTrecons        -> 4D array of CS reconstructed images, one entry 
%                         per volume in interleavedKspace.
%
%    IFTrecons         -> 4D array of IFT reconstructed images, one entry 
%                         per volume in interleavedKspace.
%
% EXTERNALS
%   BART toolbox (Lustig et al)
%   ndSparse.m (written by Matt Jacobson)
%   ifft3c.m (written by Lustig, modified by Nathan Murtha)

%% ------------------------------------------------------------------------
%  ---------------------   SET UP FOR EXECUTION   -------------------------
%  ------------------------------------------------------------------------

% Get dimensions of kspace array (i.e. phantom dimensions).
numRows   = size(interleavedKspace,1);
numCols   = size(interleavedKspace,2);
numSlices = size(interleavedKspace,3);
numCoils  = size(interleavedKspace,4);
numVols   = size(interleavedKspace,5);

% Get coil sensitivity maps.
coilMaps = getSensitivityMaps3D([numRows, numCols, numSlices],numCoils);


%% ------------------------------------------------------------------------
%  -------------------   PERFORM IMAGE RECONSTRUCTIONS   ------------------
%  ------------------------------------------------------------------------

% Allocate space for reconstructions.
BARTrecons = zeros(numRows,numCols,numSlices,numVols);
%IFTrecons  = zeros(numRows,numCols,numSlices,numCoils,numVols);

% Loop over interleaved k-space data and reconstruct images.
parfor iVol = 1:numVols  

    % Reconstruct images without compressed sensing.
    %IFTrecons(:,:,:,iVol) = abs(ifft3c( full(interleavedKspace(:,:,:,:,iVol)) ));
    
    % Reconstruct image using BART with user specified regularization.
    % W means wavelet, T means TV, L means local low rank. 
    str = ['pics -S ', BARTregularization];
    
    BARTrecons(:,:,:,iVol) = abs( bart( str,...   % regularization 
        full(interleavedKspace(:,:,:,:,iVol)),... % k-space data
        coilMaps ) );                            % coil sensitivities
end




end

