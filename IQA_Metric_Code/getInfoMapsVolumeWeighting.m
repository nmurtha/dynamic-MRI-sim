function [iw_map]= getInfoMapsVolumeWeighting(pyrRef,pyrDis,...
    numScales,infoWinRows,infoWinCols,infoWinSlices,sigma_nsq)

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           March 9, 2017, edited April 4 2017.
%
%  PURPOSE: Extend Wang's info_content_weight_map.m code to include 3D
%           images. Option to include parent coefficient is excluded for
%           now, coding convenience. 
%
%           NOTE: See info_content_weight_map.m for original code as
%                 provided by Wang.
% 
%  INPUTS:
%   pyrRef        -> Reference image pyramid (cell array, each element 
%                    containing either a 2D or 3D image array). 
%   
%   pyrDis        -> Distorted image pyramid (cell array, each element 
%                    containing either a 2D or 3D image array). 
%
%   numScales     -> The number of subbands in the pyramid that will be 
%                    considered in calculating the score.
%                    Note: info weighting is only calculated for the first
%                          numScales-1 scales. The last scale is simply
%                          filled with ones, such that no effective 
%                          weighting scheme is applied (see IWSSIM paper).
%                          This is because the info weighting applies for a
%                          zero mean Gaussian Scale Mixture model of 
%                          transform coefficients, not valid in the highest
%                          scale of the pyramid.
%                 
%   infoWinRows   -> Size of info weight calculation window, rows.
%
%   infoWinCols   -> Size of info weight calculation window, columns.
%
%   infoWinSlices -> Size of info weight calculation window, slices.
%
%   sigma_nsq     -> Variance of noise in HVS model (see IW-SSIM).
% 
%  OUTPUTS:
%   iw_map        -> Map of info weights for each level of the image 
%                    pyramid (cell array, each cell containing a 2D or 3D 
%                    array).
%
%  EXTERNALS:
%   [NONE]

%% ------------------------------------------------------------------------
%  ----------------------   PERFORM SAFETY CHECKS   -----------------------
%  ------------------------------------------------------------------------

% Check that input images are equal in size. Top level of pyramid is
% sufficient for check, should all be reduced by same factor beyond that.
assert(isequal(size(pyrRef{1}),size(pyrDis{1})),...
    'Images must be same size.')

% Check that the info window is odd in size.
assert(mod(infoWinRows,2)   == 1, 'Info window row span must be odd.');
assert(mod(infoWinCols,2)   == 1, 'Info window column span must be odd.');
assert(mod(infoWinSlices,2) == 1, 'Info window slice span must be odd.');

%% ------------------------------------------------------------------------
% ----------------------   SET DEFAULTS IF NEEDED   -----------------------
% -------------------------------------------------------------------------

% Set number of scales for info weighting to the same as in the pyramid 
% cell array by default. 
if (~exist('numScales','var'))
   numScales = numel(pyrRef);
end

% Set info weighting window to be 3x3x1 for 2D images and 3x3xn for 3D
% images by default, where n is either 3 or 2 (if the image has two
% slices).
if (~exist('infoWinRows','var'))
   infoWinRows = 3;
end
if (~exist('infoWinCols','var'))
   infoWinCols = 3;
end
if (~exist('infoWinSlices','var'))
	infoWinSlices = min( size(pyrRef{1},3), 3 );
end

% Set default human visual system model variance as quoted by Wang.
if (~exist('sigma_nsq','var'))
   sigma_nsq = 0.4; % default from [Sheikh & Bovik, IEEE-TIP 2006]
end 

%% ------------------------------------------------------------------------
%  --------------------   SET CONSTANTS FOR EXECUTION   -------------------
%  ------------------------------------------------------------------------

% Set tolerance parameter for division stability and thresholding.
tol = 1e-15;

%% ------------------------------------------------------------------------
%  ----------------------  GET INFORMATION MAPS   -------------------------
%  ------------------------------------------------------------------------

% Calculate information maps at all but the last scale, as per IW-SSIM
% paper (no info weighting done at highest scale, since info model assumes 
% a zero-mean Gaussian scale mixture and the low-pass portion of the 
% Laplacian pyramid is not zero-mean).
iw_map = cell(1,numScales-1);
for iScale=1:numScales-1
    
    % Get the current pyramid scale transform images.
    imgRef = pyrRef{iScale};
    imgDis = pyrDis{iScale};
    
    % Get size of current pyramid transform image.
    [rows, cols, slices] = size(imgRef);
    dimsPyrImgs          = [rows cols slices];
    
    % Get distortion "gain" estimation window. The gain factor is estimated
    % with a 3x3 window if the image is 2D OR if the user specified
    % in-plane information weighting for a 3D image. The gain factor is
    % estimated with a 3x3x2 or 3x3x3 window otherwise, using 3x3x2 only if
    % the image has two slices. 
    % NOTE: assumption is made in VIF paper (from which IW-SSIM draws from)
    % that distortion types are locally approximated as D = gR + V, hence
    % small window for g estimation.
    if infoWinSlices == 1
        gEstimationWin = ones(3,3,1);  
    else
        gEstimationWin = ones( 3, 3, min(size(imgRef,3), 3) );
    end
    gEstimationWin = gEstimationWin ./ sum(gEstimationWin(:)); % averaging window
    
    % Calculate local mean, variance and covariance maps. Set any negative 
    % values of variance to 0, as per Wang, since variance cannot be
    % negative.
    % NOTE: for 2D images, these arrays will be 2D. For 3D images, these
    % arrays will be 3D. 
    muRefPyr     = imfilter(imgRef, gEstimationWin, 'replicate');
    muDisPyr     = imfilter(imgDis, gEstimationWin, 'replicate');
    covRefDisPyr = imfilter(imgRef.*imgDis, gEstimationWin, 'replicate') - muRefPyr.*muDisPyr;
    varRefPyr    = imfilter(imgRef.^2, gEstimationWin, 'replicate') - muRefPyr.^2;
    varDisPyr    = imfilter(imgDis.^2, gEstimationWin, 'replicate') - muDisPyr.^2;
    varRefPyr(varRefPyr<0) = 0; % Set negative variances to zero.
    varDisPyr(varDisPyr<0) = 0;
    
    % Estimate distortion model noise variance and gain factor. See
    % equation 6, 32 and 33 in the IW-SSIM paper.
    % NOTE: for 2D images, these arrays will be 2D. For 3D images, these
    % arrays will be 3D.  
    distortGain     = covRefDisPyr ./ (varRefPyr+tol); % estimate gain factor, g, noting that cov(X,Y) = expectation(XY^T) 
    distortNoiseVar = (varDisPyr - distortGain.*covRefDisPyr); % estimate noise variance, sigma_V^2  
    
    % Clean up distortion and variance estimates. If values in the maps are
    % below the tolerance, they're likely just zero but caught up in some
    % kind of round-off error (either way, they're insignificant). If the
    % variance in either image is near zero, there's likely no blur effect
    % to model locally, and hence g is zero (see VIF paper). Likewise, if
    % the variance in the distorted image is zero, there's probably no 
    % noise to model.
    distortGain(varRefPyr < tol)     = 0; 
    distortGain(varDisPyr < tol)     = 0; 
    distortNoiseVar(varRefPyr < tol) = varDisPyr(varRefPyr < tol);
    distortNoiseVar(varDisPyr < tol) = 0;   
    
    % Set up for local information calculations. Get the number of pixels
    % that will serve as the center of local neighbourhoods for R and D. In
    % getting these local windows, the sliding window is NOT centered over
    % edge pixels, since that would include pixels in R and D that didn't
    % actually exist in the array (or more precisely, for the code below,
    % would have been wrapped around from the sides and the top). 
    %
    % In the example below, a 5x5 image is partitioned with a 3x3 info
    % window. * denotes members of one local neighbourhood, and ^ denotes
    % members of another local neighbourhood. Each local neighbourhood will
    % constitute a local R and D Gaussian Scale Mixture model for the
    % reference and distorted images. Note that neighbourhoods can share
    % members, since the window is incremented by 1, but neighbourhoods are
    % NOT centered on edge pixels to avoid cases where the neighbourhood
    % would extend past existing values in the image array. To deal with
    % the fact that edge pixels are "trimmed" from the information
    % calculation, the edge values will be replicated at the end of this
    % function so that the size of the information maps matches the size of
    % the input pyramid images at each scale.
    %
    %                       _____________________
    %                       | * | * | * |   |   |
    %                       ---------------------
    %                       | *^| *^| *^|   |   |   Note that there will be
    %                       ---------------------   9 local windows made in
    %                       | *^| *^| *^|   |   |   this case. 
    %                       ---------------------
    %                       | ^ | ^ | ^ |   |   |
    %                       ---------------------
    %                       |   |   |   |   |   |
    %                       ---------------------
    infoWinDims         = [infoWinRows, infoWinCols, infoWinSlices]; 
    numRowsGrouping     = rows-infoWinDims(1)+1;           % avoid rows that cause info window to "hang off edge" 
    numColsGrouping     = cols-infoWinDims(2)+1;           % avoid cols that cause info window to "hang off edge" 
    numSlicesGrouping   = slices-infoWinDims(3)+1;         % avoid slices that cause info window to "hang off edge" 
    numCoefsCenterWin   = numRowsGrouping*numColsGrouping*numSlicesGrouping; % number of coefficients in transform image that are centered in info window for purposes of neighbourhood collections
    numPixInfoWin       = prod(infoWinDims);               % number of coefficients included in each local info calculation
    kernRowsOffCenter   = (infoWinDims(1)-1) / 2;          % overhang of info window in row direction
    kernColsOffCenter   = (infoWinDims(2)-1) / 2;          % overhang of info window in column direction
    kernSlicesOffCenter = (infoWinDims(3)-1) / 2;          % overhand of info window in slice direction
    
    % Rearrange "numCoefsCenterWin" neighborhoods into local window groups, 
    % which constitute the local neighbourhoods for R and D. The resulting 
    % array, localCoefNeighbs, is a 2D array containing sliding window 
    % "groups" along each row (i.e. a single row comprises a group of
    % transform coefficients that form the local R and D models).
    localCoefsNeighbs = zeros(numCoefsCenterWin,numPixInfoWin); % one row for each set of pixels in sliding window progress.
    iNeighb = 0;                                                % index for local neighbourhood
    if numel(dimsPyrImgs) == 2 % 2D image case
        for iRow=-kernRowsOffCenter:kernRowsOffCenter % rather than shift local window, shift the image! Faster.
            for iCol=-kernColsOffCenter:kernColsOffCenter
                iNeighb   = iNeighb + 1;
                imShifted = circshift(imgRef, [iRow iCol]); 
                imShifted = imShifted(...
                    kernRowsOffCenter+1:kernRowsOffCenter+numRowsGrouping,...
                    kernColsOffCenter+1:kernColsOffCenter+numColsGrouping); % Trims edge pixels
                localCoefsNeighbs(:,iNeighb) = imShifted(:);
            end
        end
   else % 3D image case
       for iRow=-kernRowsOffCenter:kernRowsOffCenter % rather than shift kernel, shift the image! Faster.	
           for iCol=-kernColsOffCenter:kernColsOffCenter
               for iSlice=-kernSlicesOffCenter:kernSlicesOffCenter
                   iNeighb   = iNeighb + 1;
                   imShifted = circshift(imgRef, [iRow iCol iSlice]); 
                   imShifted = imShifted(...
                       kernRowsOffCenter+1:kernRowsOffCenter+numRowsGrouping,...
                       kernColsOffCenter+1:kernColsOffCenter+numColsGrouping,...
                       kernSlicesOffCenter+1:kernSlicesOffCenter+numSlicesGrouping); % Trims edge pixels
                   localCoefsNeighbs(:,iNeighb) = imShifted(:);
               end
           end
       end
    end % Could add parent contribution if desired, appropriate resizing required. Omitted for convenience at the moment.

    % Estimate parameters from the GMS assumption of the signal, R = sU.
    %
    % First get the covariance matrix of U, C_u. This is done via eqn 29 in
    % the IW-SSIM paper. Floating point errors may cause some eigenvalues 
    % of C_u to be very small negative numbers, which isn't possible since
    % variance of data is necessarily positive. To correct for this, the
    % negative eigenvalues are set equal to zero, and the remaining
    % eigenvalues are scaled such that the trace of the covariance matrix
    % remains constant (trace(C_u) is equal to the sum of the eigenvalues,
    % and represents a measure of variance of the data, since the 
    % eigenvalues themselves represent variance along the principle axes of
    % the data). See example for a 3x3 eigenvalue matrix below (i.e. L in 
    % the code below):
    %                 _________
    %                | a  0  0 |    suppose that a < 0. Want to preserve: 
    %         C_u =  | 0  b  0 |         trace(C_u) = a + b + c     
    %                | 0  0  c |    set a = 0, and rescale b and c so that: 
    %                -----------         nb + nc = trace(C_u) for int n
    %                                     n = trace(C_u) / (b + c)
    %
    % Then use C_u to estimate s with eqn 30 in the IW-SSIM paper. Because
    % of the cropping of edge pixels during local windowing for information
    % calculation (see loops above), the size of localCoefNeighbs
    % corresponds to an image that is the size of imgRef LESS the pixels we
    % cropped out (e.g. a 256x256 image, windowed into 3x3 local
    % neighbourhoods, becomes 254x254). To account for this, must pad
    % sSquared after estimating it, assuming that edge values don't change
    % significantly from their neighbours. 
    C_u      = (localCoefsNeighbs' * localCoefsNeighbs)/numCoefsCenterWin; 
    [Q,L]    = eig(C_u);                       % eigenvalues in L with orthogonal matrix Q
    L        = diag( diag(L).*(diag(L)>0) )... % this term sets negative eigenvalues to 0...
             * sum( diag(L) ) / (...           % ... and these next three lines are one term that scale the trace of the remaining eigenvalues
               sum(diag(L).*(diag(L)>0) )...
               + ( sum( diag(L).*(diag(L)>0) )==0 )... % NOTE: the sum(*)==0 term exists to prevent division by 0.
               );                              % correct negative eigenvalues, maintaining trace of original C_u
    C_u      = Q*L*Q';                         % reconstruct covariance matrix with adjusted eigenvalues/
    L        = diag(L);                        % store eigenvalues as a column vector
    
    sSquared = (localCoefsNeighbs/C_u).*localCoefsNeighbs/numPixInfoWin;
    sSquared = sum(sSquared,2);
    sSquared = reshape(sSquared,...
        [dimsPyrImgs(1) - 2*kernRowsOffCenter,...
         dimsPyrImgs(2) - 2*kernColsOffCenter,...
         dimsPyrImgs(3) - 2*kernSlicesOffCenter]); % sSquared will be smaller than size of image array due to cropping of edge pixels during info windowing above
    sSquared = padarray(sSquared,...
        [kernRowsOffCenter, kernColsOffCenter, kernSlicesOffCenter],...
        'replicate','both'); % pad edges of sSquared, assuming small padding requirement and slowly changing sSquared values...
     
    % Get info weight map for current scale. See eqn, 28 in the paper. If
    % any of the weights are less than the tolerance, just set them to
    % zero.
    % NOTE: in his original code info_content_weight_map.m, Wang omits the
    % factor of 1/2. It is true that this factor will cancel out when
    % calculating the IW-SSIM, but it is included here to preserve the true
    % value of the mutual information in the event that these maps are used
    % independently of the IW-SSIM.
    infow = zeros(size(imgRef));
    for iEigVal=1:length(L) % sum over each eigenvalue
        infow = infow + ...
            log2(1 + distortNoiseVar./sigma_nsq + (...
            distortNoiseVar./(sigma_nsq^2) + (1+distortGain.^2)/sigma_nsq )...
            .* sSquared.*L(iEigVal) );
    end
    infow = infow ./ 2;
    infow(infow < tol) = 0;
    
    % Assign info map at current scale.
    iw_map{iScale} = infow;
end  

% Populate last scale weight map with all ones, all pixels have equivalent 
% weight. This is done because the last scale of a Laplacian pyramid does
% not fit with the Gaussian Scale Mixture modelling used in the information
% calculation above (i.e. not zero mean with heavy tails).
iw_map{end+1} = ones(size(pyrRef{end}));