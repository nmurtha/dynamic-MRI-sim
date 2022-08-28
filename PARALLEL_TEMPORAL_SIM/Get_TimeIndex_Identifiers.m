function [ timeIdentifiers ] = Get_TimeIndex_Identifiers( numKspaceAcq,...
    CIRCUScalib, CIRCUSquanta, TR)

%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Feb 15, 2017
% 
%  PURPOSE: Finds what time indices will be sampled by each pattern, so
%           that a set of identifying time indices can be used in a parfor
%           sampling loop (must pass specific values to each worker in a
%           parfor loop, can't have an "overall time index").
%
%           Works by counting the number of frequency encodes to be sampled
%           in each case, and simply populating a vector of linearly
%           increasing integers for each pattern. In the special case of
%           snapshot sampling, time indices are a single number,
%           corresponding to the number of the current sampled array.
%
%  INPUTS: 
%   numKspaceAcq   -> Number of acquired kspace calibration/quanta pairs.
%                     NOTE: Collecting equal number of calib and quanta
%                           patterns, such that for numKspaceAcq the user
%                           will actually receive 2*numKspaceAcq sampled
%                           arrays back. This guarantees that sampling
%                           always starts with a calibration region and
%                           ends with a quanta, for convenience. E.g. for
%                           numKspaceAcq = 3:
%                                           C Q C Q C Q
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
%  OUTPUTS:
%   timeIdentifiers -> Cell array containing vectors of time indices for
%                      each sampled kspace array. 
%                      NOTE: This means, for example, that the first C will
%                            consist of time indices 1 to n, and the next
%                            sampled Q will consist of indices (n+1) to m,
%                            and so on. For the case of snapshot sampling,
%                            each "vector" will actually be a single
%                            number.
%
%  EXTERNALS:
%   [none]

%% ------------------------------------------------------------------------
%  -----------------------   SET UP FOR THE RUN   -------------------------
%  ------------------------------------------------------------------------

% Set up a variable to hold the running sum of time indices. This is used
% in the for loop as a place holder for the time indices.
runningSumTimeInds = 0;

% Get number of pregenerated quanta.
numPregenQuanta = length(CIRCUSquanta);

% Allocate cell array for results.
timeIdentifiers = cell(1,2*numKspaceAcq);

%% ------------------------------------------------------------------------
%  ----------   GET TIME INDEX IDENTIFIERS FOR "REAL" SAMPLING   ----------
%  ------------------------------------------------------------------------

% Real sampling is performed if TR > 0
if TR ~= 0
    % Loop over all patterns to be acquired, counting up TR's.
    for iAcquired = 1:2*numKspaceAcq

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

        % Count the number of points in the current pattern. This will be
        % the number of entries in its phase table.
        numPoints = numel(samplingPattern.phaseTableDim1);

        % Populate vector of time indices for this particular pattern.
        timeIdentifiers{iAcquired} = runningSumTimeInds + (1:numPoints);

        % Update time index tracker.
        runningSumTimeInds = runningSumTimeInds + numPoints;

    end % end of for loop

%% ------------------------------------------------------------------------
%  --------   GET TIME INDEX IDENTIFIERS FOR "SNAPSHOT" SAMPLING   --------
%  ------------------------------------------------------------------------
    
else
    % Loop over all patterns to be acquired, assigning sampling index as
    % the single time index for snapshot sampling (because all sampling
    % occurs in a single instant).
    for iAcquired = 1:2*numKspaceAcq
        % Populate vector of time indices for this particular sampled pattern.
        timeIdentifiers{iAcquired} = iAcquired;
    end
    
end

end

