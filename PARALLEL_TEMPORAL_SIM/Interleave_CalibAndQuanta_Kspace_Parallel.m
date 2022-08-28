function [interleavedKspace, aveSampleTimes] = Interleave_CalibAndQuanta_Kspace_Parallel(...
    sampledKspace, patternTimes, numCoils,...
    nCalInter, nQuantaInter, nCalIncrement,timeWeighting)

%% ------------------------------------------------------------------------
%  --------------------   FUNCTION INFORMATION   --------------------------
%  ------------------------------------------------------------------------
% 
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Jan 26, 2017
%
%  PURPOSE: Interleaves acquired CIRCUS calibration regions and quanta,
%           averaging any overlapping points according to their frequency
%           of occurance (e.g. a voxel that was sampled 2 times between 6
%           interleaved patterns is averaged by 2, not 6). Reports
%           effective time of sampling by averaging the sampling times
%           provided for each k-space array.
%
%  INPUTS:
%   sampledKspace  -> Array of sampled k-space data. A 4D array, with the
%                     first three dimensions corresponding to the sampled
%                     k-space and the fourth dimension corresponding to the
%                     time at which it was sampled. Array is an ndSparse
%                     object, to reduce RAM usage.
%
%   patternTimes   -> Vector of times labelling the sampledKspace data. One
%                     time entry, in seconds, per volume in sampledKspace.
%
%   numCoils       -> Number of coils used in parallel imaging. A scalar
%                     value, either 1 or an even number (i.e. 1,2,4,6,...).
%
%   nCalInter      -> Number of calibration regions to interleave into
%                     each pattern.
%
%   nQuantaInter   -> Number of quanta to interleave into each pattern.
%
%   nCalIncrement  -> Incrementation to use between "anchor" calibration
%                     regions while combining data. 
%                     NOTE: See example below.
%                      ncalIncrement = 1 : C  Q  C  Q  C  Q  C  Q  C  Q ... 
%                                          ^     ^     ^     ^     ^
%
%                      nCalIncrement = 2 : C  Q  C  Q  C  Q  C  Q  C  Q ...
%                                          ^           ^           ^
%
%  timeWeighting   -> A string, 'calib' or 'equal', specifying whether the
%                     aveSampleTimes output variable includes times (from
%                     input variable patternTimes) from only calibration
%                     regions or equally weighted between calibration
%                     regions and quanta. 
%
%  OUTPUTS:
%   interleavedKspace -> 5D array of resulting interleaved kspace data. An
%                        ndSparse array, to save RAM. Fourth index
%                        represents coil number, while fifth index
%                        represents interleaved volume number.
%
%   aveSampleTimes    -> Vector containing the average times of the
%                        interleaved kspace data, in seconds.
%
%  EXTERNAL FUNCTIONS:
%   ndSparse.m (written by Matt Jacobson)


%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR EXECUTION   -----------------------
%  ------------------------------------------------------------------------

% Catch case where the user asks to interleave either more calibration
% regions or more quanta than were sampled.
assert(nCalInter <= 0.5*size(sampledKspace,5),...
    'Cannot combine %d calibration regions when only %d are acquired.',...
    nCalInter,0.5*size(sampledKspace,5))

assert(nQuantaInter <= 0.5*size(sampledKspace,5),...
    'Cannot combine %d quanta when only %d are acquired.',...
    nQuantaInter,0.5*size(sampledKspace,5))

% Catch improper timeWeighting specification.
assert(strcmpi(timeWeighting,'calib') || strcmpi(timeWeighting,'equal'),...
    'Must specify either "calib" or "quanta" for timeWeighting (case insensitive).');

if strcmpi(timeWeighting,'calib') == 1
    assert( nCalInter > 0,...
        'Specified time weighting to use calibration regions, but set calibration interleaving amount to 0.');
end

% Allocate space for interleaved kspace arrays and a vector for their 
% average sampled times. 
interleavedKspace = ndSparse(...
    zeros( size(sampledKspace, 1),... % number rows in k-space
           size(sampledKspace, 2),... % number cols in k-space
           size(sampledKspace, 3),... % number slices in k-space
           size(sampledKspace, 4),... % number coils used in parallel imaging
           length(1:nCalIncrement:0.5*size(sampledKspace,5)) ) ); % want to "anchor" interleaves on a C, hence the 0.5* factor.

aveSampleTimes = zeros(1, size(interleavedKspace,5));

% Get vector of which cal regions serve as anchor (in "cal index space").
calAnchors = 1:nCalIncrement:0.5*size(sampledKspace,5);

%% ------------------------------------------------------------------------
%  --------------------   INTERLEAVE THE KSPACE DATA   --------------------
%  ------------------------------------------------------------------------

% Loop over all calibration regions. Note that in the code below, we make
% reference to "iCalib" and "iQuanta", the index of the calibration region
% and quanta respectively. These count the number of calibration regions
% and quanta that we are on relative to the beginning of the data set:
%                           C  Q  C  Q  C  Q  C  Q ...
%                 iTotal:   1  2  3  4  5  6  7  8 ...
%                 iCalib:   1     2     3     4    ...
%                 iQuanta:     1     2     3     4 ...
%
% Note that the following two relations hold, easily derivable by assessing
% the linear relationships between the indices:
%                 iTotal = 2*iCalib - 1   (i.e. the odd numbers)
%                 iTotal = 2*iQuanta      (i.e. the even numbers)
%
% The loop increments the calibration region index according to the
% increment specified by the user, the information for which is contained
% in the vector calAnchors.
parfor iCalib = 1:numel(calAnchors)
    % ----------------  DETERMINE RANGE OF CALIB INDICES  ----------------
    
    % Check if there are sufficient calibration regions "left" of current
    % calibration region index.
    if calAnchors(iCalib) - floor( (nCalInter-1)/2 ) < 1 % floor will bias index to the right here, arbitrary choice
        flagMissingCalLeft = 1;
        missingCalLeft     = 1 - ( calAnchors(iCalib) - floor((nCalInter-1)/2) ); 
    else
        flagMissingCalLeft = 0;
    end
    
    % Check if there are sufficient calibration regions right of the
    % current calibration region index.
    if calAnchors(iCalib) + ceil( (nCalInter-1)/2 ) > 0.5*size(sampledKspace,5) % ceil is used because floor was used above, complimentary notions
        flagMissingCalRight = 1;
        missingCalRight = ( calAnchors(iCalib) + ceil( (nCalInter-1)/2 ) ) - ...
            0.5*size(sampledKspace,5);
    else
        flagMissingCalRight = 0;
    end
    
    % At this point, we know if there are sufficient calibration regions
    % around the current calibration region in each direction, which guides
    % the combination of calibration regions. Must treat three cases to
    % determine the actual bounds on indices of the calibration regions
    % that will be combined. First two cases are if there is a lack
    % of data on the left, or a lack of data on the right (relative to the
    % current calibration region index iCalib). Last case is when there is
    % sufficient data surrounding the current calibration region for a
    % symmetric combination about iCalib.
    if flagMissingCalLeft && ... % When not enough cal regions on left (e.g. for first sampled C)
            (calAnchors(iCalib) + ceil((nCalInter-1)/2) + missingCalLeft) <= 0.5*size(sampledKspace,5)
        iCalLeft  = 1;
        iCalRight = calAnchors(iCalib) + ceil((nCalInter-1)/2) + missingCalLeft; % Append missing data onto the right side of iCalib
    elseif flagMissingCalRight && ... % when not enough cal regions on the right (e.g. for last sampled C)
            (calAnchors(iCalib) - floor((nCalInter-1)/2) - missingCalRight) >= 1
        iCalLeft  = calAnchors(iCalib) - floor((nCalInter-1)/2) - missingCalRight; % Append missing data onto the left side of iCalib
        iCalRight = 0.5*size(sampledKspace,5);
    else % when there are sufficient cal regions on either side of iCalib
        iCalLeft  = calAnchors(iCalib) - floor((nCalInter-1)/2);
        iCalRight = calAnchors(iCalib) + ceil((nCalInter-1)/2);
    end

    
    % --------------  DETERMINE RANGE OF QUANTA INDICES  ----------------
    
    % At this point, we have established which calibration regions will be
    % combined. In a similar matter, must establish which quanta will be
    % combined. Start using iQuanta = iCalib, possible to do so since
    % calibration regions and quanta are always acquired in pairs. Then
    % finding the bounds of quanta indices should proceed the same way as
    % above for bounds on calibration region indices.
    
    % Set current quanta index.
    iQuanta = calAnchors(iCalib);
    
    % Check if there are sufficient quanta left of current quanta index.
    if iQuanta - floor( (nQuantaInter-1)/2 ) < 1  % floor will bias index to the right, arbitrary choice
        flagMissingQuantaLeft = 1;
        missingQuantaLeft     = 1 - ( iQuanta - floor( (nQuantaInter-1)/2 ) );
    else
        flagMissingQuantaLeft = 0;
    end
    
    % Check if there are sufficient quanta right of current quanta index.
    if iQuanta + ceil( (nQuantaInter-1)/2 ) > 0.5*size(sampledKspace,5) % ceil is used because floor was used above, complimentary notions
        flagMissingQuantaRight = 1;
        missingQuantaRight = (iQuanta + ceil( (nQuantaInter-1)/2 )) - ...
            0.5*size(sampledKspace,5);
    else
        flagMissingQuantaRight = 0;
    end
    
    % Similar to the case above with calibration region combinations, the
    % three cases for quanta combination.  
    if flagMissingQuantaLeft && ... % quanta missing on left side (e.g. for first sampled quanta)
            (iQuanta + ceil((nQuantaInter-1)/2) + missingQuantaLeft) <= 0.5*size(sampledKspace,5)
        iQuanleft  = 1;
        iQuanright = iQuanta + ceil((nQuantaInter-1)/2) + missingQuantaLeft; % Append missing data to right side
    elseif flagMissingQuantaRight && ... % quanta missing on right side (e.g. for last sampled quanta)
            (iQuanta - floor((nQuantaInter-1)/2) - missingQuantaRight) >= 1
        iQuanleft  = iQuanta - floor((nQuantaInter-1)/2) - missingQuantaRight; % Append missing data to left side
        iQuanright = 0.5*size(sampledKspace,5);
    else % when there are sufficient cal regions on either side of iQuanta
        iQuanleft  = iQuanta - floor((nQuantaInter-1)/2);
        iQuanright = iQuanta + ceil((nQuantaInter-1)/2);
    end   
    
    % -------------------  SUM DATA AND FIND OVERLAPS  -------------------
    
    % At this point, we've determined the indices of the calibration
    % regions and quanta that will be combined to make the current 
    % interleaved array. Have to now pull those patterns out from the 
    % sampledKspace array, add them up, and average each voxel by its 
    % frequency of occurance. To do this, we convert each "calibration 
    % index" and "quanta index" into the "true index" of the sampled
    % k-space array.
    %
    % First, declare a 3D array to hold the sum of kspace data (that will
    % later be averaged), as well as a 2D matrix whos rows hold information
    % regarding whether or not a linear index in the currently examined 3D
    % array of data corresponds to an element with nonzero data. Also
    % declare a variable for holding sum of the sampled times for each
    % pattern, in order to calculate the average time of the resulting
    % interleaved pattern.
    
    for iCoil = 1:numCoils
        % Allocate 3D array to hold sum of all interleaved k-space arrays.
        sumDataKspace = zeros( size( sampledKspace(:,:,:,1,1) ) );

        % Allocate 2D array whos rows will indicate whether or not a linear
        % index in each interleaved pattern contains nonzero data (i.e. along
        % each row, the columns represent linear indices of points in
        % consituent k-space arrays that are being combined, and take values of
        % 0 or 1 depending on whether or not that linear index is nonzero).
        linIndsWithData = zeros(...
            (iCalRight-iCalLeft+1) + (iQuanright-iQuanleft+1),... % As many rows as there are combined patterns
            numel(sampledKspace(:,:,:,1,1)) ); % As many columns as elements in k-space array

        % Initialize variable to hold the sum of sampling times for the
        % constituent k-space data, such that the times may be averaged later
        % to make an effective time for the interleaved result.
        sumTimes = 0;
        flagAddedCalib = 0;

        % Add all calibration region data for current coil to the summed data 
        % array, and record which linear index in each pattern held nonzero
        % data (allows us to average each voxel by frequency of data occurance
        % later).
        for iC = iCalLeft:iCalRight
            % Get index of volume in sampledKspace corresponding to the current
            % calibration region index.
            iTotal = 2*iC - 1; % odd indices in sampledKspace

            % Pull calibration region data from sampledKspace.
            calData = sampledKspace(:,:,:,iCoil,iTotal);   

            % Get effective time of calData.
            flagAddedCalib = 1;
            sumTimes = sumTimes + patternTimes(iTotal);

            % Add calibration region data to summed array of data.
            sumDataKspace = sumDataKspace + calData;

            % Find all of the linear indices of calData with nonzero entries, 
            % and mark those linear indices in the 2D matrix of linear index
            % information, linIndsWithData.
            linIndsWithData(iC-iCalLeft+1, calData~=0) = 1;        
        end

        % Repeat above, but for quanta data now. Note that sumDataKspace now
        % contains the sum of all calibration data that is being interleaved.
        for iQ = iQuanleft:iQuanright
            % Get index of volume in sampledKspace corresponding to the current
            % calibration region index.
            iTotal = 2*iQ; % even indices in sampledKspace

            % Pull calibration region data from sampledKspace.
            quantaData = sampledKspace(:,:,:,iCoil,iTotal);

            % Get effective time of quantaData.
            if flagAddedCalib == 0 || strcmpi(timeWeighting,'equal')
                sumTimes = sumTimes + patternTimes(iTotal);
            end

            % Add quanta region data to summed array of data.
            sumDataKspace = sumDataKspace + quantaData;

            % Find all of the linear indices of quantaData with nonzero 
            % entries, and mark those linear indices in the 2D matrix of linear 
            % index information. Note that we must increment row index in 
            % linear index matrix to avoid overwriting data from the 
            % calibration regions, hence the (iCalRight-iCalLeft+1) term.
            linIndsWithData((iCalRight-iCalLeft+1)+iQ-iQuanleft+1,...
                quantaData~=0) = 1;        
        end   

        % --------------  INTERLEAVE DATA AND AVERAGE OVERLAPS  --------------

        % At this stage, have an array full of the sum of all k-space data that
        % will be interleaved, as well as a matrix whos rows carry information
        % about whether or not a linear index for an constituent array of data
        % contained data or not. Will average the data in sumDataKspace by
        % taking linIndsWithData, summing along columns to determine frequency
        % that each linear index held data, replacing zeros in the resulting
        % summed vector by ones to avoid divisions by zero, and dividing each
        % element of sumDataKspace by the appropriate index frequency from the
        % resulting summed vector, indFreqOccurData.
        %
        % NOTE: If you don't understand why we're summing linIndsWithData,
        % consider the picture below. Each row represents the linear index of a
        % consituent array that was summed into sumDataKspace. Entries of 1 in
        % a row mean that the linear index of that array held nonzero data,
        % hence a sampled kspace entry. By stacking all this linear index info
        % as rows in a 2D matrix, we can look along the columns to determine if
        % a given linear index was sampled many times. Example for combining
        % one calibration region and two quanta may be:
        % 
        %               C  | 0  1  1  0  0  1  0  1  1  1 ... 
        %               Q  | 0  1  0  0  1  1  1  1  0  0 ...
        %               Q  | 1  0  0  0  0  1  1  1  1  0 ...
        %                  ------------------------------ ...
        % indFreqOccurData   1  2  1  0  1  3  2  3  2  1 ...   

        % Determine frequency of occurance of data at each linear index.
        indFreqOccurData = sum(linIndsWithData,1);

        % To avoid division by 0, find and replace all zeros in
        % indFreqOccurData with ones. This won't affect the average of
        % overlapping k-space data, because division by one has no effect.
        indFreqOccurData(indFreqOccurData == 0) = 1;

        % Average each voxel in sumDataKspace by the number of times which it
        % was sampled.
        sumDataKspace(:) = sumDataKspace(:) ./ indFreqOccurData(:);

        % Store the interleaved, averaged k-space data.
        interleavedKspace(:,:,:,iCoil,iCalib) = sumDataKspace;
    end % end iCoil loop
    
    % Get average time of the interleaved data.
    if flagAddedCalib == 1 && strcmpi(timeWeighting,'calib')
        aveSampleTimes(iCalib) = sumTimes / ...
            ((iCalRight-iCalLeft+1));
    else
        aveSampleTimes(iCalib) = sumTimes / ...
            ((iCalRight-iCalLeft+1) + (iQuanright-iQuanleft+1));
    end
    
end % end of interleaving loop! 

%% ------------------------------------------------------------------------
%  ----------------------   TRIM EQUIVALENT ARRAYS   ----------------------
%  ------------------------------------------------------------------------
%
%  For some combinations of nCalInter and nQuantaInter, due to the wrapping 
%  pf interleaved arrays near the periphery of sampledKspace,
%  interleavedKspace may have equivalent volumes near the beginning and end
%  of the volume data. That is to say, for example, that interleaved volume
%  1 and 2 may be identical, volumes end and end-1 may be identical, and so 
%  on. 
%
%  To avoid returning redundant data, this section of code "trims the fat"
%  and returns only unique combinations of interleaved kspace data. It does
%  this by looping over interleavedKspace volumes, beginning to end, and
%  marking indices of volumes that should be deleted. If two arrays are
%  equivalent, the *later* array is marked for deletion (hence the loop
%  over 1:size(interleavedKspace,4)-1, -1 to account for the fact that loop
%  extends its look to size(interleavedKspace,4)).
%
%  For example, if 6 pairs of calibration regions and quanta are acquired,
%  with nCalIncrement = 1, nCal = 3, nQuanta = 4, we should get the first
%  two and last two volumes of interleavedKspace as identical, and so
%  delIndices = [2,6]. See figure below:
%
%                       C  Q  C  Q  C  Q  C  Q  C  Q  C  Q
%  interleaved vol:     1  1  1  1  1  1     1 
%                       2  2  2  2  2  2     2 (equals vol 1, delete)
%                          3  3  3  3  3  3  3
%                                4  4  4  4  4  4  4     
%                                      5  5  5  5  5  5  5
%                                      6  6  6  6  6  6  6 (equals vol 5, delete)

% Set counter for index of deleted indices. Also initialize a vector for
% indices of volumes to be deleted. It's empty at first, but will be
% updated if any redundant volumes exist.
iDelete    = 0;
delIndices = [];

% Loop over interleavedKspace volumes, deleting redundant data. Should
% suffice to check only the first coil for equality in interleaved array,
% since all others followed same interleaving path.
for iVol = 1:size(interleavedKspace,5)-1
    currentVol = interleavedKspace(:,:,:,1,iVol);  
    nextVol    = interleavedKspace(:,:,:,1,iVol+1);
    if isequal(currentVol, nextVol) == 1
       iDelete = iDelete + 1;
       delIndices(iDelete) = iVol + 1; % Want to delete the later occurance of equivalent data
    end
end

% Delete the redundant data, if any exists.
if numel(delIndices) > 0 % check if delIndices ~= []
    fprintf('Deleting interleaved volume %d, redundant data.\n',delIndices)
    interleavedKspace(:,:,:,:,delIndices) = [];
    aveSampleTimes(delIndices) = [];
end



end