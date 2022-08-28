function [sampledKspace, timesBegSamp, timesMeanSamp, timesDoneSamp] = Sample_Kspace_Parallel(...
    phantBase, phantFeatures,...
    numKspaceAcq,...
    numCoils,...
    samplingTimes, TR,... 
    CIRCUScalib, CIRCUSquanta,...
    noiseFlag, noisePower)

%% ------------------------------------------------------------------------
%  --------------------   FUNCTION INFORMATION   --------------------------
%  ------------------------------------------------------------------------
% 
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Jan 26, 2017
%
%  PURPOSE: Samples a temporally evolving phantom using alternating
%           calibration region and CIRCUS quanta patterns in pairs, such
%           that the data may be easily interleaved later if the user
%           desires. Sampling follows the following order:
%
%                       C Q C Q C Q .... C Q C Q
%
%           Time stamps marking the beginning, mean, and end time of
%           sampling for each pattern are returned.
%           
%           NOTE: The frequency encode direction is along the length of the
%                 rows of k-space. Two phase encode tables select the row
%                 and slice of k-space that will be sampled, and every 
%                 column along the resulting row is "sampled" instantly 
%                 during one elapsed TR. The assumption here is that the
%                 duration of a frequency encode is very small, and can be
%                 assumed negligible compared to TR. 
%
%  INPUTS:
%   phantBase      -> Array containing the static background of the dynamic
%                     phantom. A 3D array.
%
%   phantFeatures  -> Structure array containing the information of each 
%                     feature in the phantom. Fields are:
%                       > phantFeatures.shape   (3D array with shapes)
%                       > phantFeatures.Ct      (temporal evolution vector)
%                       > phantFeatures.params  (temporal parameters)
%
%   numKspaceAcq   -> Number of acquired kspace calibration/quanta pairs.
%                     NOTE: Collecting equal number of calib and quanta
%                           patterns, such that for numKspaceAcq the user
%                           will actually receive 2*numKspaceAcq sampled
%                           arrays back. This guarantees that sampling
%                           always starts with a calibration region and
%                           ends with a quanta (for convenience). E.g. for
%                           numKspaceAcq = 3:
%                                           C Q C Q C Q
%
%   numCoils       -> Number of coils used in parallel imaging. A scalar
%                     value, either 1 or an even number (i.e. 1,2,4,6,...).
%
%   samplingTimes  -> All the times, in seconds, for which sampling takes 
%                     place.
%
%   TR             -> Repition time, in seconds, for the scan.
%                     NOTE: Only used as an indicator for realistic or
%                           snapshot sampling, since samplingTimes contains
%                           the actual times used in sampling.
%
%   CIRCUScalib    -> Structure array containing the phase tables for the 
%                     calibration region pattern. Contains two vectors:
%                       > CIRCUScalib.phaseTableDim1
%                       > CIRCUScalib.phaseTableDim2
%                     Since all calibration regions are equivalent
%                     trajectories, only one set of phase tables is in
%                     CIRCUScalib.
%
%   CIRCUSquanta    -> Structure array containing the phase tables for the 
%                      effective CIRCUS quanta. Contains two vectors:
%                       > CIRCUSquanta.phaseTableDim1
%                       > CIRCUSquanta.phaseTableDim2
%                      Has as many phase tables as there were pregenerated 
%                      CIRCUS quanta.
%
%   noiseFlag      -> Flag that dictates whether or not phantom image is
%                     corrupted by Rician noise before taking FFT and
%                     sampling. 1 equals "on", 0 equals "off".
%
%   noisePower     -> The approximate "offset" in signal intensity 
%                     introduced by the Rician noise.
%
%  OUTPUTS: 
%   sampledKspace  -> Array of sampled kspace data. A 5D array, with the
%                     first three dimensions corresponding to the sampled
%                     kspace, the fourth dimension corresponding to coil
%                     number, and the fifth corresponding to acquired 
%                     CIRCUS quanta. Array is an ndSparse object, to reduce
%                     RAM usage.
%
%   timesBegSamp   -> Vector of times at which each pattern began being
%                     sampled. Units of seconds.
%
%   timesMeanSamp  -> Vector of times at which each pattern was sampled on
%                     average. Units of seconds.
%
%   timesDoneSamp  -> Vector of times at which each pattern was finished
%                     being sampled. Units of seconds.
%
%  EXTERNAL FUNCTIONS:
%   getSensitivityMaps3D.m
%   ricernd.m (written Ged Ridgway)
%   fft3c.m   (written by Nathan Murtha, trivial extension of fft2c by Lustig)
%   Get_TimeIndex_Identifiers.m 
%   ndSparse.m (written by Matt Jacobson)

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR SAMPLING  -------------------------
%  ------------------------------------------------------------------------

% Grab dimension of the phantom, which correspond to k-space dimensions.
numRows   = size(phantBase,1);
numCols   = size(phantBase,2);
numSlices = size(phantBase,3);

% Get number of pregenerated CIRCUS quanta.
numPregenQuanta = length(CIRCUSquanta);

% Get number of features in the phantom.
numFeatures = length(phantFeatures);

% Get cell array containing vectors of time indices for each sampled CIRCUS
% pattern. This allows sampling with a parfor loop, avoiding indexing
% problems.
timeIdentifiers = Get_TimeIndex_Identifiers( numKspaceAcq, CIRCUScalib,...
    CIRCUSquanta, TR);

% Get coil sensitivity maps.
coilMaps = getSensitivityMaps3D([numRows, numCols, numSlices], numCoils);

%% ------------------------------------------------------------------------
%  -------------------   SAMPLE THE KSPACE (REAL)   -----------------------
%  ------------------------------------------------------------------------

% Realistic sampling is done if TR is nonzero.
if TR ~= 0
    % Set up array to be filled with kspace measurements. Note: the 2* 
    % factor arises because we want to sample pairs of calibration regions 
    % and quanta. Array is ndSparse to conserve memory.
    sampledKspace = ndSparse( ...
        zeros(numRows, numCols, numSlices,... % k-space dimensions
        numCoils,...                          % coil index
        2*numKspaceAcq) );                    % acquired pattern index

    % Allocate vectors marking the beginning, mean, and end time of 
    % sampling for each sampled CIRCUS quanta.
    timesBegSamp  = zeros(1,2*numKspaceAcq);
    timesDoneSamp = zeros(1,2*numKspaceAcq);
    timesMeanSamp = zeros(1,2*numKspaceAcq);

    % Loop over number of acquired quanta (fifth index in output array)
    parfor iAcquired = 1:2*numKspaceAcq       
        % Get the time index identifiers for this pattern (i.e. for 
        % accessing sampling times). Will be a vector of *indices*.
        timeIndices = timeIdentifiers{iAcquired};

        % Determine if the current pattern is a calibration region or a 
        % quanta. Note that calibration regions are always acquired for odd 
        % values of iAcquired, and quanta are always acquired for even 
        % values of iAcquired.
        if mod(iAcquired,2) == 1 % calibration region
            iPattern = -1;
        else % quanta
            iEven = 0.5 * iAcquired; % The index of the current "even" number (e.g. 6 is the 3rd even whole number).
            iPattern = mod(iEven,numPregenQuanta);
            if iPattern == 0 % catch case where mod "resets the clock" to zero.
                iPattern = numPregenQuanta;
            end
        end

        % Assign the current sampling pattern.
        if iPattern > 0
            samplingPattern = CIRCUSquanta(iPattern);
        else
            samplingPattern = CIRCUScalib; 
        end

        % Record time at which current pattern began to be sampled. The
        % same pattern is sampled for all coils "in parallel" (at least, in
        % real life), so all coils begin sampling at the same time.
        timesBegSamp(iAcquired) = samplingTimes( timeIndices(1) );
        
        % Sample k-space for each coil.
        for iCoil = 1:numCoils
            % Assign temporary kspace array to be filled with samples from 
            % the temporally evolving k-space. This array will be saved 
            % into the 5D output array after the sampling loop.
            tempSampKspace = zeros(numRows, numCols, numSlices);
            
            % Sample all the frequency encodes as seen by coil iCoil.
            for iFreqEncode = 1:numel(samplingPattern.phaseTableDim1) % NOTE: phaseTableDim1 and phaseTableDim2 have equal length      
                % Calculate the current phantom from feature temporal 
                % evolutions.
                currentPhant = phantBase;
                for iFeat = 1:numFeatures
                    currentPhant = currentPhant + ...
                        ( phantFeatures(iFeat).shape .* ...
                          phantFeatures(iFeat).Ct(timeIndices(iFreqEncode))...
                        );
                end

                % Add Rician noise to the phantom image if noise flag is  
                % on. Note that the random number generation is seeded by 
                % the time index, so that the noise is easily reproducible 
                % across coils for each volume. 
                if (noiseFlag == 1)
                    rng( timeIndices(iFreqEncode) );
                    currentPhant = ricernd(currentPhant,noisePower);
                end

                % Get k-space of current phantom. NOTE: BART expects
                % orthonormal and centered FFT data, as per the output of 
                % fft3c.m. The k-space is taken from what the current coil
                % "sees", i.e. (sensitivity)*(image).
                kspaceCurrent = fft3c(coilMaps(:,:,:,iCoil).*currentPhant); 

                % Get row and slice location for the current frequency
                % encode, using phase encode coordinate tables.
                y = samplingPattern.phaseTableDim1(iFreqEncode);
                z = samplingPattern.phaseTableDim2(iFreqEncode);

                % Populate current k-space array with frequency encoded  
                % data from the current image.
                tempSampKspace(y,:,z) = kspaceCurrent(y,:,z);
            end
            % Push temporary kspace data to output kspace array.
            % tempSampKspace contains data for the iAcquired-th pattern, as
            % seen by coil iCoil.
            sampledKspace(:,:,:,iCoil,iAcquired) = tempSampKspace;
        end

        % Record time at which iAcquired-th pattern was finished being
        % sampled. Must add 1 to the index because sampling BEGAN at time
        % marked by timeIndices(end). All coils finish sampling at the same
        % time.
        timesDoneSamp(iAcquired) = samplingTimes(timeIndices(end) + 1);

        % Record time at which iAcquired-th pattern was sampled on average.
        % This is done by taking the end time of the current sampling, and  
        % subtracting half the time it took to sample the pattern (rounding
        % if needed to keep the resulting time on the same mesh as 
        % samplingTimes). The coils all sampled the same mean time.
        timesMeanSamp(iAcquired) = timesDoneSamp(iAcquired) - ...
            round( 0.5 * numel(samplingPattern.phaseTableDim1) ) * TR;
    end

%% ------------------------------------------------------------------------
%  ------------------   SAMPLE THE KSPACE (SNAPSHOT)   --------------------
%  ------------------------------------------------------------------------
    
% If TR = 0, sampling the kspace patterns with "snapshot" sampling. This
% means that the entire CIRCUS pattern is sampled instantly, with no
% temporal evolution between frequency encodes, and then a time lag is
% allowed before the next pattern is "snapshotted". 
else
    % Set up array to be filled with kspace measurements. Note: the 2* 
    % factor arises because we want to sample pairs of calibration regions 
    % and quanta. Array is ndSparse to conserve memory.
    sampledKspace = ndSparse(...
        zeros(numRows, numCols, numSlices,... % size of k-space array
        numCoils,...                          % coil index
        2*numKspaceAcq) );                    % acquired pattern index

    % Loop over number of acquired k-space arrays.  
    parfor iAcquired = 1:2*numKspaceAcq
        % Get the time index identifier for this pattern (will be a single 
        % index, due to snapshot sampling at a single instant in time).
        timeIndex = timeIdentifiers{iAcquired};        
        
        % Determine if the current pattern is a calibration region or a 
        % quanta. Note that calibration regions are always acquired for odd
        % values of iAcquired, and quanta are always acquired for even 
        % values of iAcquired.
        if mod(iAcquired,2) == 1 % calibration region
            iPattern = -1;
        else % quanta
            iEven = 0.5 * iAcquired; % The index of the current "even" number (e.g. 6 is the 3rd even whole number).
            iPattern = mod(iEven,numPregenQuanta);
            if iPattern == 0 % catch case where mod "resets the clock" to zero.
                iPattern = numPregenQuanta;
            end
        end

        % Assign the current sampling pattern.
        if iPattern > 0
            samplingPattern = CIRCUSquanta(iPattern);
        else
            samplingPattern = CIRCUScalib;
        end
        
        % Make a 3D mask of the current sampling pattern. First allocate
        % space for the mask, then use the phase tables from
        % samplingPattern to make the first column of the 3D mask, and then
        % extend that 3D mask along the rest of the columns. NOTE: we
        % extend the mask along the columns because that is the direction
        % the frequency encoding runs along.
        sampledMask3D = zeros( size(phantBase) );
        for iPoint = 1:numel(samplingPattern.phaseTableDim1)
            sampledMask3D( samplingPattern.phaseTableDim1(iPoint),...
                1,...
                samplingPattern.phaseTableDim2(iPoint) ) = 1;
        end
        for iCol = 2:size(phantBase,2)
            sampledMask3D(:,iCol,:) = sampledMask3D(:,1,:);
        end

        % Calculate the current phantom from feature temporal evolutions.
        currentPhant = phantBase;
        for iFeat = 1:numFeatures
            currentPhant = currentPhant + ...
                ( phantFeatures(iFeat).shape .* ...
                  phantFeatures(iFeat).Ct( timeIndex ) );
        end

        % Add Rician noise to the phantom image if noise flag is on. Note
        % that the random number generation is seeded by the time index, so
        % that the noise is easily reproducible in each volume. 
        if (noiseFlag == 1)
            rng(timeIndex);
            currentPhant = ricernd(currentPhant,noisePower);
        end
    
        % For each coil, get the k-space data.
        for iCoil = 1:numCoils
            % Get kspace of current image. NOTE: BART expects orthonormal 
            % and centered FFT data, as per the output of fft3c.m. The 
            % k-space is taken from what the current coil "sees", i.e. 
            % (sensitivity)*(image).
            kspaceCurrent = fft3c(coilMaps(:,:,:,iCoil).*currentPhant); 

            % Populate current k-space array with frequency encoded data 
            % from the current image.
            sampledKspace(:,:,:,iCoil,iAcquired) = sampledMask3D .*...
                kspaceCurrent;
        end
    end
    
    % Set the output time vectors. For perfect sampling, the beginning
    % time, mean time, and end sampling time are all the same. All times
    % are the same for each coil, since the coils sample k-space in
    % parallel.
    timesBegSamp  = samplingTimes;
    timesMeanSamp = samplingTimes;
    timesDoneSamp = samplingTimes;
    
end


end