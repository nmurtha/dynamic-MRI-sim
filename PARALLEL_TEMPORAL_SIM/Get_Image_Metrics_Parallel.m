function [ RMSEscore, GMSDscore, SSIMscore, MSSSIMscore, IWSSIMscore] = ...
    Get_Image_Metrics_Parallel( phantBase, phantFeatures, ...
        samplingTimes, reconImages, effTimes )

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           May 5, 2017
%
%  PURPOSE: Calculates image quality metrics for the reconstructed images. 
%           It does this by calculating what the phantom looked like at the
%           effective time for which the image was reconstructed, and 
%           comparing the image to that phantom as "ground truth".
%
%           The metrics calculated are:
%              > Root mean squared error (RMSE)
%              > Gradient magnitude similarity deviation (GMSD)
%              > Mean structural similarity index (SSIM)
%              > Multi-scale structural similarity index (MSSSIM)
%              > Information weighted structural similarity index (IWSSIM)
%       
%
%  INPUTS:
%   phantBase      -> Array containing the static background of the dynamic
%                     phantom. A 3D array.
%
%   phantFeatures  -> Structure array containing the information for each
%                     feature in the phantom. Fields are:
%                       > phantFeatures.shape   (3D array, feature shape)
%                       > phantFeatures.Ct      (Temporal evolution vector)
%                       > phantFeatures.params  (Temporal parameters)
%
%   samplingTimes  -> All the times for which sampling takes place. A 
%                     vector in units of seconds.
%  
%   reconImages    -> 4D array of reconstructed image data, where the
%                     fourth dimension corresponds to the volume index.
%                     Assuming magnitude images.
%
%   effTimes       -> The effective time of sampling of images in
%                     reconImages. A vector, in units of seconds.
%
%  OUTPUTS:
%   RMSEscore      -> Vector of RMSE values for each image in reconImages.
%
%   GMSDscore      -> Vector of GMSD values for each image in reconImages.
%
%   SSIMscore      -> Vector of SSIM values for each image in reconImages.
%
%   MSSSIMscore    -> Vector of MSSSIM values for each image in
%                     reconImages.
%
%   IWSSIMscore    -> Vector of IWSSIM values for each image in
%                     reconImages.
%
%
%  EXTERNALS
%   RMSE_ND.m 
%   GMSD_ND.m (extension of GMSD author code into 3D, Nathan Murtha)
%   SSIM_ND.m (modified version of the stock ssim.m, Nathan Murtha)
%   MSSSIM_ND.m (extension of MSSSIM source code to 3D, Nathan Murtha)
%   IWSSIM_ND.m (extension of IWSSIM source code to 3D, Nathan Murtha)

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR EXECUTION   -----------------------
%  ------------------------------------------------------------------------

% Allocate vectors for metric results.
RMSEscore   = zeros(1,size(reconImages,4));
SSIMscore   = zeros(1,size(reconImages,4));
GMSDscore   = zeros(1,size(reconImages,4));
MSSSIMscore = zeros(1,size(reconImages,4));
IWSSIMscore = zeros(1,size(reconImages,4));

% Get the time between k-space phase encoded samples.
dT = samplingTimes(2) - samplingTimes(1);


%% ------------------------------------------------------------------------
%  ----------------------   CALCULATE IMAGE METRICS   ---------------------
%  ------------------------------------------------------------------------

% Loop over each reconstructed image and get the image quality metric.
for iImage = 1:size(reconImages,4)    
    % Must find the time in samplingTimes that is closest to
    % effTimes(iImage), in order to reconstruct the phantom to get "truth".
    % Begin by finding the indices of possible times in samplingTimes that
    % match effTimes(iImage). The effective time will either land between
    % two sampled times, or directly overlap a sampled time. If the
    % effective time equals a sampled time, we'll only have one possible
    % index, and we're done. Otherwise, determine which possible index
    % corresponds to a time that is closer to a sampled time. If the
    % effective time lands *directly* in between two sampled times, then
    % the choice is made (arbitrarily) to bias the time forward to the
    % later sampled time.
    possibleInds = find(abs( effTimes(iImage) - samplingTimes ) <= dT/2); % inds where effTime(iImage) is within dT/2 of samplingTimes(:)    
    if numel(possibleInds) == 1 % if effective time is already a sampled time
        indTrueTime = possibleInds;
    elseif numel(possibleInds) == 2 % if effective time lands between two sampled times
        if abs(samplingTimes(possibleInds(1)) - effTimes(iImage)) < ... % If effective time is closer to left bounding sampled time
                abs(samplingTimes(possibleInds(2)) - effTimes(iImage))
            indTrueTime = possibleInds(1);
        else % if effective time is equally distant from surrounding sampled times OR closer to the right bounding time
            indTrueTime = possibleInds(2);
        end
    end
    
    % Create the true phantom image at the sampled time corresponding to
    % indTrueTime. This is used as "truth" for quality metric calculation.
    phantTrue = phantBase;
    for iFeat = 1:length(phantFeatures)
        phantTrue = phantTrue + ...
            ( phantFeatures(iFeat).shape .* ...
              phantFeatures(iFeat).Ct(indTrueTime) );
    end 
    
    % Calculate the RMSE of the image.
    RMSEscore(1,iImage) = RMSE_ND( phantTrue, reconImages(:,:,:,iImage) );

    % Calculate the GMSD of the image. 
    GMSDscore(1,iImage) = GMSD_ND( phantTrue, reconImages(:,:,:,iImage),...
        0, max(phantTrue(:)));
    
    % Calculate the mean SSIM of the image.
    SSIMscore(1,iImage) = SSIM_ND( phantTrue, reconImages(:,:,:,iImage),...
        [11 11 11],...       % SSIM calculation window dimensions
        [1 1 1],...          % SSIM term exponents 
        max(phantTrue(:)) ); % max image value

    % Calculate the volume based MSSSIM of the image.
    MSSSIMscore(1,iImage) = MSSSIM_ND(phantTrue, reconImages(:,:,:,iImage),...
        3,...         % number of scales
        [11 11 5],... % ssim term calculation window dimensions
        2,...         % side length of square window used in pyramid generation
        max(phantTrue(:)) ); % max image value

    % Calculate the volume based IWSSIM of the image.
    IWSSIMscore(1,iImage) = IWSSIM_ND(phantTrue, reconImages(:,:,:,iImage),...
        3,...         % number of scales
        [11 11 5],... % ssim term calculation window dimensions
        [3 3 3],...   % dimensions of window for information calculation
        5,...         % side length of square window used in pyramid generatio
        max(phantTrue(:)) ); % max image value
end

end
