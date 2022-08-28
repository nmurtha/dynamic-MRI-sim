function [ samplingTimes ] = Get_Sampled_Times( numKspaceAcq, CIRCUScalib,...
                CIRCUSquanta, TR, dTimePerfect )
            
%% ------------------------------------------------------------------------
%  --------------------   FUNCTION INFORMATION   --------------------------
%  ------------------------------------------------------------------------
% 
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Jan 26, 2017
%
%  PURPOSE: Takes in the phase tables for sampled calibration region and
%           CIRCUS quanta and calculates each sampled time, storing them in
%           a vector. The acquisition order is assumed to be equal numbers
%           of calibration regions (C) and quanta (Q):
%                           C Q C Q C Q ... C Q C Q   
%
%           NOTE: If TR equals 0, the function assumes that snapshot
%                 sampling of kspace will be used. In this case, only one
%                 sampling time per CIRCUS pattern will be generated,
%                 starting from 0 seconds and proceeding with a lag of
%                 dTimePerfect seconds.
%
%  INPUTS:
%   numKspaceAcq   -> Number of acquired CIRCUS calibration/quanta pairs. A
%                     scalar input.
%                     NOTE: Collecting equal number of calib and quanta
%                           patterns, such that for numKspaceAcq the user
%                           will actually receive 2*numKspaceAcq sampled
%                           arrays back. This guarantees that sampling
%                           always starts with a calibration region and
%                           ends with a quanta (for convenience). E.g. for
%                           numKspaceAcq = 3:
%                                           C Q C Q C Q
%
%   CIRCUScalib    -> Structure array containing the phase tables for the
%                     calibration region pattern. Contains two vectors:
%                       > CIRCUScalib.phaseTableDim1
%                       > CIRCUScalib.phaseTableDim2
%
%   CIRCUSquanta   -> Structure array containing the phase tables for the
%                     CIRCUS quanta patterns. Contains two vectors:
%                       > CIRCUSquanta.phaseTableDim1
%                       > CIRCUSquanta.phaseTableDim2
%
%   TR             -> Repition time, in seconds, between frequency encodes.
%                     A scalar input.
%
%   dTimePerfect   -> Time lag, in seconds, between snapshot sampled
%                     kspace patterns. A scalar input.
%                     NOTE: This parameter is only ever used if TR = 0.
%
%  OUTPUTS: 
%   samplingTimes  -> All the times, in seconds, for which sampling takes 
%                     place. A vector output.
%
%  EXTERNALS:
%   [NONE]

%% ------------------------------------------------------------------------
%  -------------------   CHECK INPUTS FOR BAD CASES   ---------------------
%  ------------------------------------------------------------------------

% Ensure that TR and dTimePerfect are positive, can go backwards in time!
assert(TR >= 0, 'TR must be either 0 or a positive number.');
assert(dTimePerfect > 0, 'dTimePerfect must be exclusively greater than 0.');

%% ------------------------------------------------------------------------
%  -----------------   SAMPLED TIMES WITH NONZERO TR   --------------------
%  ------------------------------------------------------------------------

% If TR is nonzero, count TR's for each CIRCUS pattern.
if TR ~= 0
    % Set initial simulated time value to 0 seconds. Will update in loop.
    simulatedTimeLength = 0; % [sec]

    % Get number of pregenerated quanta. Assuming here that all calibration
    % regions are equivalent k-space trajectories.
    numPregenQuanta = length(CIRCUSquanta);

    % Loop over all patterns to be acquired, counting up TR's. Note: the 2*
    % factor arises because calibration regions and quanta are collected as
    % pairs.
    for iAcquired = 1:2*numKspaceAcq    
        % Determine if the current pattern is a calibration region or a 
        % quanta. Note that calibration regions are always acquired for odd 
        % values of iAcquired, and quanta are always acquired for even 
        % values of iAcquired.
        if mod(iAcquired,2) == 1 % calibration region
            iPattern = -1;
        else % quanta
            iEven = 0.5 * iAcquired; % The index of the current "even" number (e.g. 6 is the 3rd even whole number).
            iPattern = mod(iEven,numPregenQuanta); % get which quanta pattern to use, re-using quanta if sampling exceeds pregenerated numbers.
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

       % Count next TR's from the current sampling pattern.
       simulatedTimeLength = simulatedTimeLength +...
           TR * numel(samplingPattern.phaseTableDim1); % one TR per phase encode.
    end % end of iAcquired loop.

    % Get vector of sampled times. The sampling time starts at 0 seconds
    % (introduction of bolus) and runs until all points have been sampled, 
    % with sampling occuring every TR seconds. 
    samplingTimes = 0:TR:simulatedTimeLength; 

    % Ensure that samplingTimes(end) actually equals simulatedTimeLength.
    % Floating point errors may cause miscount of one TR in bad cases 
    % (rare, but happens just often enough crash a run, so we handle it).
    while (abs(samplingTimes(end) - simulatedTimeLength) > 0.01*TR)
        if (samplingTimes(end) < simulatedTimeLength) % if undercounted, add time to fix
            samplingTimes(end+1) = samplingTimes(end) + TR;
        elseif(samplingTimes(end) > simulatedTimeLength) % if overcounted, remove excess time to fix
            samplingTimes(end) = [];
        end
    end
    samplingTimes(end) = simulatedTimeLength; % ensure exact equivalence

%% ------------------------------------------------------------------------
%  -------------------   SAMPLED TIMES WITH ZERO TR   ---------------------
%  ------------------------------------------------------------------------    
    
% If TR = 0 then we generate times for "snapshot sampled" patterns.
% This means that at each time step, the entire CIRCUS pattern will be
% sampled in an instant and then the time will "lag" by dTimePerfect before
% another perfect sampling takes place. 
else
    samplingTimes = 0:dTimePerfect:dTimePerfect*(2*numKspaceAcq-1); % subtract 1 since we start from time 0 seconds.
end % end of if statement set above.


end