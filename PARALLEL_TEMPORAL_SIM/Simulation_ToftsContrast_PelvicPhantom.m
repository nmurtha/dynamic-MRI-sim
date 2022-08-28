%  ------------------------------------------------------------------------
%  ---------------------   ABOUT THIS SIMULATION   ------------------------
%  ------------------------------------------------------------------------
%
% WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%
% PURPOSE: Numerical phantom simulation to broadly assess the effects of
%          acquisition and reconstruction methodologies on CS
%          reconstructions and accuracy of recovered parameters modelling
%          temporal dynamics.
%
%          This version includes the use of the ndSparse class, in an
%          effort to avoid saving data to matfiles during execution. This
%          should greatly increase execution time. Most of the large arrays
%          are ndSparse, and must have full(*) appended to them in order to
%          use them for certain functions (e.g. with BART). 
%
%          New to this version is the use of parallel processing via parfor 
%          loops. 
%
%
% EXTERNALS:
%   -> Generate_CIRCUS_PhaseTable.m
%   -> Get_Sampled_Times.m
%   -> Get_Pelvic_Phantom.m
%   -> Sample_Kspace_Parallel.m
%   -> Interleave_CalibAndQuanta_Parallel.m
%   -> Do_Image_Recons_Parallel.m
%   -> Get_Image_Metrics_Parallel.m
%   -> Get_Tofts_Parameters.m
%   -> guardedDividePercentError.m
%   -> Get_Mean_Errors_Parallel.m
%   -> Get_Mean_Timecourse_Parallel.m
%   -> ndSparse.m (written by Matt Jacobson)

%% ------------------------------------------------------------------------
%  ----------------------   SET UP FOR EXECUTION   ------------------------
%  ------------------------------------------------------------------------
%
%  This section sets up for the execution of the simulation by clearing the
%  memory, starting a parellel pool for parallel processing (user can
%  adjust the number of CPU's requested), directing the MATLAB environment
%  to the BART toolbox, and creating a timestamped folder for the data to
%  be saved within.

% Clear memory, close windows, clear junk text.
clear;
close all;
clc;

% Set BART toolbox path (i.e. top level directory of BART toolkit).
setenv('TOOLBOX_PATH','/Users/nmurtha/bart-0.3.01');
%setenv('TOOLBOX_PATH','/bmrl/matlab/toolbox/bart-0.3.00')

% Update user with simulation progress.
disp(strcat('SIMULATION: Setting up for execution. (', ...
    datestr(datetime('now')), ')' ))

% Start timer for overall runtime.
tic

% Initialize a local parallel pool with N_CPU CPU's.
N_CPU = 2;
POOL = parpool('local',N_CPU); 

% Make timestamped directory to save simulation results in.
timeStamp = datetime;
timeStr   = datestr(timeStamp,'mmm_dd_yyyy_HH-MM-SS'); % format datetime to string
fileStr   = strcat('Sim_',timeStr); % create directory name for results
mkdir(fileStr); % create the directory

% Clear unnecessary variables.
clear N_CPU timeStamp timeStr

%% ------------------------------------------------------------------------
%  ---------------------   SET SIMULATION PARAMETERS   --------------------
%  ------------------------------------------------------------------------
% 
%  User defined simulation parameters are set here. 
%
%  Note: For flags, 1 = on and 0 = off. 

% Update user with simulation progress.
disp(strcat('SIMULATION: Setting simulation parameters. (', ...
    datestr(datetime('now')), ')' ))

% Set Rician noise properties. Flag = 0 means noiseless, flag = 1 will
% distribute the phantom IMAGE according to Rician distribution.
flagRicianNoise = 0;
noisePower      = 1000; 

% Set the acquisition parameters. NOTE: If TR is set to 0 seconds,
% simulation will enter "snapshot mode" (i.e. instant sampling of CIRCUS
% patterns in k-space, each pattern seperated by some time dTimePerfect set 
% below).
TR              = 0.000;  % [sec], one TR between frequency encodes (i.e. points in CIRCUS phase tables).
dTimePerfect    = 1;      % [sec], time lag between snapshot sampled patterns (only used if TR = 0).
numKspaceAcq    = 30;      % Pairs of calib/quanta sampling patterns to acquire (e.g. numKspaceAcq = 2 samples 2 pairs of CQ)

numCoils        = 1;      % Number of coils to use in parallel imaging.
 
numCIRCUSPregen = 1;    % Number of CIRCUS patterns to pregenerate (will be looped through if needed).
CIRCUSmode      = 2;      % =1 specifies # quanta to return, =2 specifies undersampling factor.
CIRCUSreturn    = 1;     % quanta per pattern for mode=1, undersampling factor for mode=2.
CIRCUSb         = 40;    % radial CIRCUS parameter ("shear" in quanta).
CIRCUSc         = 1.5;   % spiral CIRCUS parameter ("twist" in quanta).  
CIRCUSseed      = 0;     % seed < 0 will produce random quanta.

% Set CIRCUS interleaving parameters.
numCalInterleave    = 0;      % number of calibration regions to interleave.
numQuantaInterleave = 1;      % number of quanta to interleave.
numCalIncrement     = 1;      % jump between calibration regions when interleaving.
patternTimePriority = 'equal';% 'calib' to set effective times based on calibration regions, 'equal' to include quanta with calibration regions 

% Set BART CS parameters. E.g. '-R W:7:0:0.02 -R L:7:7:0.01'. Syntax goes
% as '-R A:B:C:D', where:
%     A -> Sparsity model for regularization (i.e. T, L or W)
%     B -> Dimension to perform transform along, bitmasked (e.g. 7 =
%          dimensions 0, 1 and 2, the spatial dimensions). 
%     C -> Joint-thresholding dimensions (I think it places dependencies on
%          the sparse model BETWEEN the dimensions of the transform...?).
%          For T and W, set to 0. For L, set to 7.
%     D -> Regularization weight (e.g. 0.01).
BARTregularization = '-R W:7:0:0.00'; 

%% ------------------------------------------------------------------------
%  --------------------  PREGENERATE CIRCUS PATTERNS  --------------------- 
%  ------------------------------------------------------------------------
%  
%  CIRCUS patterns are pregenerated here. Phase table information for the 
%  calibration region and the quanta are stored in structure arrays.
%  The structure arrays produced here, and their fields, are as follows:
%   
%   CIRCUScalib   -> CIRCUScalib.phaseTableDim1 (vector of coordinates)
%                 -> CIRCUScalib.phaseTableDim2 (vector of coordinates)
%
%   CIRCUSquanta  -> CIRCUSquanta(i).phaseTableDim1 contains the phase 
%                    table coordinates in the first specified dimension for 
%                    the i-th quanta (vector of coordinates).
%                 -> CIRCUSquanta(i).phaseTableDim2 contains the phase 
%                    table coordinates in the second specified dimension  
%                    for the i-th quanta (vector of coordinates).

% Update user with simulation progress.
disp(strcat('SIMULATION: Pregenerating CIRCUS patterns. (',...
    datestr(datetime('now')), ')' ))

% Get pregenerated CIRCUS pattern phase tables. Returns one calibration
% region, and numQuantaPregen "effective" quanta containing as many single
% CIRCUS quanta as specified by the variables CIRCUSmode and CIRCUSreturn
% (e.g. each entry of CIRCUSquanta may have been composed by adding
% together 10 single CIRCUS quanta, making a larger "effective" quanta).
[ CIRCUScalib, CIRCUSquanta ] = Generate_CIRCUS_PhaseTables(...
    256, 30,... % frequency encodes will run "in-slice", so phase encodes depend on number of rows and slices
    CIRCUSb, CIRCUSc,...             % the CIRCUS radial and spiral parameters.
    CIRCUSmode, CIRCUSreturn,...     % parameters dictating the undersampling performed by CIRCUS.
    CIRCUSseed, numCIRCUSPregen);    % Quanta seed, and number of patterns to generate in total.

%% ------------------------------------------------------------------------
%  ----------------   POPULATE VECTOR OF SAMPLING TIMES   -----------------
%  ------------------------------------------------------------------------
%
%  A vector of sampling times is populated here, by counting up the number
%  of repitition times that will have elapsed when sampling all the
%  patterns. The function called here acts to clean up any round-off error
%  introduced when counting up the number of TR's that elapse. The vector 
%  produced here is called samplingTimes. Units are in seconds. 
%
%  The last entry in samplingTimes is the time at which the last frequency
%  encode was *finished* being sampled, and not a time at which sampling
%  actually started. So, samplingTimes(end) is not a "time of interest" in
%  terms of sampling a frequency encode, but is included because it 
%  represents the end of the elapsed simulation time. See the figure below. 
%
%  Freq Encode samples:      K         K          K          K         ...
%                time:   0 -----> TR -----> 2TR -----> 3TR -----> 4TR  ...
%
%  To sample four k-space points, the clock will end at 4 TR even though no
%  point *began* to be sampled at time 4 TR.
%
%  If TR = 0, the simulation assumes "snapshot mode". This involves 
%  sampling each CIRCUS pattern instantly at times seperated by a lag of 
%  dTimePerfect (i.e. a TR of 0 ms, with a delay allowed between 
%  acquisitions for contrast evolution). In this case, one sampling time is 
%  generated per CIRCUS pattern in the following manner:
%  
%    CIRCUS Patterns:    C           Q          C          Q      ...
%               time:    0           dT         2dT        3dT    ...

% Update with simulation progress.
disp(strcat('SIMULATION: Generating sampling times vector. (',...
    datestr(datetime('now')), ')' ))

% Get sampling times.
samplingTimes = Get_Sampled_Times( numKspaceAcq,... % Number of patterns that will be sampled.
    CIRCUScalib, CIRCUSquanta,...                   % CIRCUS pattern phase tables.
    TR, dTimePerfect );                             % Time for each sampling.
            
%% ------------------------------------------------------------------------
%  --------------------   GENERATE NUMERICAL PHANTOM   --------------------
%  ------------------------------------------------------------------------
%
%  A numerical phantom is generated here. The phantom is composed of a
%  static base, upon which dynamic features are placed. 
%
%  Specifically, the phantom base is a static 3D array of a constant
%  background value, in which holes have been "carved out" at locations
%  where the dynamic features will be "inserted". The phantom features
%  constitute one field of the phantFeatures structure array (see below),
%  and are themselves 3D arrays that are nonzero only in the region
%  occupied by the feature; in this fashion, the 3D array representing the
%  feature can be added to the 3D array representing the phantom base. To
%  simulate contrast evolution in the feature, the phantFeatures structure
%  array contains a Ct field, which is a vector representing contrast
%  values over time that can be multiplied (one entry at a time) into the
%  feature shape array.
%
%  The variables returned are:
%     phantBase     -> static base of the phantom, in a 3D array.
%
%     phantFeatures -> phantFeatures.shape contains the shape of the
%                      phantom feature in a 3D array.
%                   -> phantFeatures.Ct contains signal intensity 
%                      evolution information for the phantom feature, 
%                      in arbitary units, in a vector.
%                   -> phantFeatures.params contains the true temporal
%                      parameters for the feature in a vector.                      

% Update user with simulation progress.
disp(strcat('SIMULATION: Getting numerical phantom. (',...
    datestr(datetime('now')), ')' ))

% Get the phantom.
[ phantBase, phantFeatures ] = Get_Pelvic_Phantom(...
    samplingTimes,...
    [0.5 0.2 0.1]); % specified PK parameters, [Ktrans Ve Vp].

%% ------------------------------------------------------------------------
%  -----------------   SAMPLE KSPACE OF DYNAMIC PHANTOM   -----------------
%  ------------------------------------------------------------------------
% 
%  The k-space of the dynamically evolving phantom is sampled, with either
%  realistic sampling (for TR > 0) or "snapshot sampling" (for TR = 0).
%  For realistic sampling, the CIRCUS acquisitions proceed in sequential
%  pairs of calibration region (C) and quanta (Q):
%                        C  Q  C  Q ... C  Q  C  Q
%
%  When interleaving the k-space data later on, a time must be associated
%  with each k-space array sampled here, such that an "effective" time for
%  the interleaved volume may be assigned. The user can choose to have the
%  sampled time for the data here correspond to the beginning, end, or mean
%  sampled time of each k-space array.
%
%  NOTE: Snapshot sampling still follows the CQCQCQ... sampling pattern
%  with CIRCUS, which means that odd volume will still only be a
%  calibration region. If you try to view a volume that you expect to have 
%  R = 1, and it looks like garbage, make sure that you haven't tried to an
%  image which came from a reconstructed calibration region! 
%
%  The outputs are as follows:
%   sampledKspace   -> Array of sampled kspace data. A 5D array, with the
%                      first three dimensions corresponding to the sampled
%                      kspace, the fourth dimension corresponding to coil
%                      number, and the fifth corresponding to acquired 
%                      CIRCUS quanta. Array is an ndSparse object, to 
%                      reduce RAM usage.
%    timesBegSamp   -> Vector of times, in seconds, in which each pattern
%                      *began* to be sampled.
%    timesMeanSamp  -> Vector of times, in seconds, in which each pattern
%                      was sampled *on average*.
%    timesEndSamp   -> Vector of times, in seconds, in which each pattern
%                      was *finished* being sampled.

% Update user with simulation progress.
disp(strcat('SIMULATION: Sampling kspace. (',...
    datestr(datetime('now')), ')' ))

% Sample the k-space of the dynamically evolving phantom.
[sampledKspace, timesBegSamp, timesMeanSamp, timesDoneSamp] = Sample_Kspace_Parallel(...
    phantBase, phantFeatures,...  % phantom base and feature information
    numKspaceAcq,...              % number of C Q pairs to acquire
    numCoils,...                  % number of coils used in parallel imaging
    samplingTimes, TR,...         % sampling time info. TR used only to indicate realistic or snapshot sampling.
    CIRCUScalib, CIRCUSquanta,... % CIRCUS phase tables.
    flagRicianNoise, noisePower); % Add noise, if desired.

%% ------------------------------------------------------------------------
%  --------------------   INTERLEAVE KSPACE DATA   ------------------------
%  ------------------------------------------------------------------------
%
%  The k-space data, composed of sampled calibration regions and quanta,
%  are interleaved. The algorithm anchors itself on the current calibration
%  region (which are looped over), and finds the quanta nearest to it in
%  accordance with the number of quanta that the user wants combined; it
%  likewise does this for calibration regions if the user requests that
%  calibration regions be combined. The step between anchor calibration
%  regions is given by numCalIncrement. 
%
%  For example, if 2 calibration region were to be combined with 3 quanta,
%  the algorithm would be the following in the sequence (note that *
%  represents the anchor calibration region).
%                              v   *   v   v   v
%                  C   Q   C   Q   C   Q   C   Q   C   Q   C   Q
%
%  Any points that overlap in the interleaved arrays are averaged by their
%  frequency of overlap (i.e. if 5 arrays were combined into the final
%  array, and one point saw entries in 2 of the 5 interleaved array, the
%  value stored in the final array at that point would be the sum of the 
%  two overlapping points divided by 2, not by 5).
%
%  NOTE: if you want to perfectly sampled k-space with R = 1, this
%        interleaving code can be used to omit calibration regions and just
%        "combine" 1 quanta that was snapshot sampled with an R = 1
%        undersampling. Simply set CIRCUSmode = 2, CIRCUSreturn = 1,
%        numCalInterleave = 0, numQuantaInterleave = 1, numCalIncrement =
%        1.
%
%  The output variables are: 
%     interleavedKspace -> 5D array of interleaved kspace data, averaging
%                          any overlapping points by their frequency of
%                          occurance. Fourth index represents coil, fifth
%                          represents the interleaved volume number.
%     effSampleTimes    -> Vector of effective sampling times for the
%                          interleaved kspace arrays, obtained by averaging
%                          the times of the consituent patterns specified
%                          in the second argument of the Interleave...
%                          function.

% Update user
disp(strcat('SIMULATION: Interleaving kspace. (',...
    datestr(datetime('now')), ')' ))

[interleavedKspace, effSampleTimes] =  Interleave_CalibAndQuanta_Kspace_Parallel(...
    sampledKspace, timesBegSamp,...           % sampled k-space and its associated time stamp
    numCoils,...                              % number of coils used in parallel imaging
    numCalInterleave, numQuantaInterleave,... % number of calibration regions and quanta to interleave
    numCalIncrement,...                       % number of calibration regions to skip between "anchors"
    patternTimePriority);                     % scheme for generating effSampleTimes

% Save the sampled k-space to file and clear it from RAM, won't be needed
% from now on.
save([fileStr '/SampledKspaceData.mat'],...
    'sampledKspace','timesBegSamp','timesMeanSamp','timesDoneSamp',...
    '-v7.3')
clear sampledKspace timesBegSamp timesMeanSamp timesDoneSamp

%% ------------------------------------------------------------------------
%  ---------------------  PERFORM CS RECONSTRUCTION  ----------------------
%  ------------------------------------------------------------------------
%
%  BART is used to acquire CS reconstructed images on the interleaved
%  k-space data here. Conventional image reconstructions are also returned.
%      BARTrecons      -> 4D array of CS reconstructed images, one entry  
%                         per volume in interleavedKspace.
%
%      IFTrecons       -> 4D array of IFT reconstructed images, one entry  
%                         per volume in interleavedKspace.
% Update user
disp(strcat('SIMULATION: Reconstructing images with BART. (',...
    datestr(datetime('now')), ')' ))

BARTrecons = Do_Image_Recons_Parallel( interleavedKspace,...
    BARTregularization );

% Save the interleaved k-space to file and clear it from RAM, won't be 
% needed from now on.
save([fileStr '/InterleavedKspaceData.mat'],...
    'interleavedKspace',...
    '-v7.3')
clear interleavedKspace

%% ------------------------------------------------------------------------
%  -------------------------   GET IMAGE METRICS   ------------------------
%  ------------------------------------------------------------------------
%
%  Metrics are calculated for each image here, by comparing the
%  reconstructed image to what the numerical should have looked like at the
%  effective sampling time of the interleaved k-space array from which an
%  image reconstruction was obtained. The metrics are:
%     -> Root mean squared error (RMSE).
%     -> Gradient Magnitude Similary Deviation (GMSD)
%     -> Mean Structural Similarity Index (MSSIM).
%     -> Multi-Scale Structural Similarity Index (MS-SSIM)
%     -> Information Weighted Structural Similarity Index (IW-SSIM)

% Update user
disp(strcat('SIMULATION: Calculating image metrics. (',...
    datestr(datetime('now')), ')' ))

% Get image metrics from BART CS reconstructions. 
[ RMSE_BART, GMSD_BART, SSIM_BART, MSSSIM_BART, IWSSIM_BART] = ...
     Get_Image_Metrics_Parallel( phantBase, phantFeatures, ...
         samplingTimes, BARTrecons, effSampleTimes );

% Get image metrics from IFT reconstructions. 
% [ RMSE_IFT, GMSD_IFT, SSIM_IFT, MSSSIM_IFT, IWSSIM_IFT] = ...
%     Get_Image_Metrics_Parallel( phantBase, phantFeatures, ...
%         samplingTimes, IFTrecons, effSampleTimes );
    
    
%% ------------------------------------------------------------------------
%  -------------------   RECOVER TEMPORAL PARAMETERS   --------------------
%  ------------------------------------------------------------------------
%
%  Fitting to the contrast enhancement curves is performed here, to recover
%  the parameters describing the signal evolution.
%
%  Parameters for each style of feature are:
%    Exponential -> [Amplitude, Decay Time, Constant Offset]
%    Sinusoidal  -> [Amplitude, Period, Temporal Offset, Constant Offset]
%    Linear      -> [Max Increase Value, Constant Offset]

% Update user
disp(strcat('SIMULATION: Recovering temporal parameters. (',...
    datestr(datetime('now')), ')' ))

% Get recovered parameters from slice 4 features (embedded, exponential).
disp(strcat('SIMULATION: Recovering temporal parameters from slice 4. (',...
    datestr(datetime('now')), ')' ))

paramSlice4_BART = Get_Tofts_Parameters( BARTrecons, 15,...
    effSampleTimes, phantFeatures.shape(:,:,15) );
% paramSlice4_IFT = Get_Exponential_Parameters( IFTrecons, 4,...
%     effSampleTimes, maskEmbedCyl );

% Get recovered parameters from slice 10 features (embedded, exponential).
disp(strcat('SIMULATION: Recovering temporal parameters from slice 10. (',...
    datestr(datetime('now')), ')' ))

paramSlice10_BART = Get_Tofts_Parameters( BARTrecons, 10,...
    effSampleTimes, maskEmbedCyl );
% paramSlice10_IFT = Get_Exponential_Parameters( IFTrecons, 10,...
%     effSampleTimes, maskEmbedCyl );
 
% Get recovered parameters from slice 16 features (embedded, sinusoidal).
disp(strcat('SIMULATION: Recovering temporal parameters from slice 16. (',...
    datestr(datetime('now')), ')' ))

paramSlice16_BART = Get_Tofts_Parameters( BARTrecons, 16,...
    effSampleTimes, maskEmbedCyl );
% paramSlice16_IFT = Get_Sinusoidal_Parameters( IFTrecons, 16,...
%     effSampleTimes, maskEmbedCyl );

% Get recovered parameters from slice 22 features (pin grids, linear).
disp(strcat('SIMULATION: Recovering temporal parameters from slice 22. (',...
    datestr(datetime('now')), ')' ))

paramSlice22_BART = Get_Linear_Parameters( BARTrecons, 22,...
    effSampleTimes, maskPingrids );
% paramSlice22_IFT = Get_Linear_Parameters( IFTrecons, 22,...
%     effSampleTimes, maskPingrids );

% Get recovered parameters from slice 28 features (pin grids, linear).
disp(strcat('SIMULATION: Recovering temporal parameters from slice 28. (',...
    datestr(datetime('now')), ')' ))

paramSlice28_BART = Get_Linear_Parameters( BARTrecons, 28,...
    effSampleTimes, maskPingrids );
% paramSlice28_IFT = Get_Linear_Parameters( IFTrecons, 28,...
%     effSampleTimes, maskPingrids );

% Clear unnecessary variables.
clear iFeat maskEmbedCyl

%% ------------------------------------------------------------------------
%  ---------------------   GET TRUE PARAMETER MAPS   ----------------------
%  ------------------------------------------------------------------------

% Update user
disp(strcat('SIMULATION: Constructing true parameter maps for features. (',...
    datestr(datetime('now')), ')' ))

% Get true parameters from slice 4 features (embedded, exponential).
paramSlice4_True = zeros(phantNumRows,phantNumCols,3);
for i = 1:3
    for j = 1:5
        paramSlice4_True(:,:,i) = paramSlice4_True(:,:,i) + ...
            phantFeatures(j).shape(:,:,4) .* phantFeatures(j).params(i);
    end
end

% Get true parameters from slice 10 features (embedded, exponential).
paramSlice10_True = zeros(phantNumRows,phantNumCols,3);
for i = 1:3
    for j = 6:10
        paramSlice10_True(:,:,i) = paramSlice10_True(:,:,i) + ...
            phantFeatures(j).shape(:,:,10) .* phantFeatures(j).params(i);
    end
end

% Get true parameters from slice 16 features (embedded, sinusoid).
paramSlice16_True = zeros(phantNumRows,phantNumCols,3);
for i = 1:3
    for j = 11:15
        paramSlice16_True(:,:,i) = paramSlice16_True(:,:,i) + ...
            phantFeatures(j).shape(:,:,16) .* phantFeatures(j).params(i);
    end
end

% Get true parameters from slice 22 features (pin grids, linear).
paramSlice22_True = zeros(phantNumRows,phantNumCols,2);
for i = 1:2
    paramSlice22_True(:,:,i) = paramSlice22_True(:,:,i) + ...
        phantFeatures(16).shape(:,:,22) .* phantFeatures(16).params(i);
end

% Get true parameters from slice 28 features (pin grids, linear).
paramSlice28_True = zeros(phantNumRows,phantNumCols,2);
for i = 1:2
    paramSlice28_True(:,:,i) = paramSlice28_True(:,:,i) + ...
        phantFeatures(17).shape(:,:,28) .* phantFeatures(17).params(i);
end

% Clear unnecessary variables.
clear i j

%% ------------------------------------------------------------------------
%  ----------------------   GET PERCENT ERROR MAPS   ----------------------
%  ------------------------------------------------------------------------

% Update user
disp(strcat('SIMULATION: Getting percent error maps for feature parameters. (',...
    datestr(datetime('now')), ')' ))

% Get percent error maps for slice 4.
percentErrSlice4_BART = guardedDividePercentError(...
    paramSlice4_BART, paramSlice4_True );
% percentErrSlice4_IFT = guardedDividePercentError(...
%     paramSlice4_IFT, paramSlice4_True );

% Get percent error maps for slice 10.
percentErrSlice10_BART = guardedDividePercentError(...
    paramSlice10_BART, paramSlice10_True );
% percentErrSlice10_IFT = guardedDividePercentError(...
%     paramSlice10_IFT, paramSlice10_True );

% Get percent error maps for slice 16.
percentErrSlice16_BART = guardedDividePercentError(...
    paramSlice16_BART, paramSlice16_True );
% percentErrSlice16_IFT = guardedDividePercentError(...
%     paramSlice16_IFT, paramSlice16_True );

% Get percent error maps for slice 22.
percentErrSlice22_BART = guardedDividePercentError(...
    paramSlice22_BART, paramSlice22_True );
% percentErrSlice22_IFT = guardedDividePercentError(...
%     paramSlice22_IFT, paramSlice22_True );

% Get percent error maps for slice 28.
percentErrSlice28_BART = guardedDividePercentError(...
    paramSlice28_BART, paramSlice28_True );
% percentErrSlice28_IFT = guardedDividePercentError(...
%     paramSlice28_IFT, paramSlice28_True );

%% ------------------------------------------------------------------------
%  -----------------   GET AVERAGE ERROR OF EACH FEATURE   ----------------
%  ------------------------------------------------------------------------

% Update user
disp(strcat('SIMULATION: Getting average error of each feature parameter. (',...
    datestr(datetime('now')), ')' ))

% Get masks for the embedded cylinder features.
maskOuterCyl  = phantFeatures(1).shape(:,:,4);
maskInnerRow1 = phantFeatures(2).shape(:,:,4);
maskInnerRow2 = phantFeatures(3).shape(:,:,4);
maskInnerRow3 = phantFeatures(4).shape(:,:,4);
maskInnerRow4 = phantFeatures(5).shape(:,:,4);

% Apply mask to each feature error map for the BART results.
maskedErrFeats_BART = cell(1,length(phantFeatures));

maskedErrFeats_BART{1} = percentErrSlice4_BART .* maskOuterCyl;  
maskedErrFeats_BART{2} = percentErrSlice4_BART .* maskInnerRow1;
maskedErrFeats_BART{3} = percentErrSlice4_BART .* maskInnerRow2;
maskedErrFeats_BART{4} = percentErrSlice4_BART .* maskInnerRow3;
maskedErrFeats_BART{5} = percentErrSlice4_BART .* maskInnerRow4;

maskedErrFeats_BART{6}  = percentErrSlice10_BART .* maskOuterCyl;  
maskedErrFeats_BART{7}  = percentErrSlice10_BART .* maskInnerRow1;
maskedErrFeats_BART{8}  = percentErrSlice10_BART .* maskInnerRow2;
maskedErrFeats_BART{9}  = percentErrSlice10_BART .* maskInnerRow3;
maskedErrFeats_BART{10} = percentErrSlice10_BART .* maskInnerRow4;

maskedErrFeats_BART{11} = percentErrSlice16_BART .* maskOuterCyl;  
maskedErrFeats_BART{12} = percentErrSlice16_BART .* maskInnerRow1;
maskedErrFeats_BART{13} = percentErrSlice16_BART .* maskInnerRow2;
maskedErrFeats_BART{14} = percentErrSlice16_BART .* maskInnerRow3;
maskedErrFeats_BART{15} = percentErrSlice16_BART .* maskInnerRow4;

maskedErrFeats_BART{16} = percentErrSlice22_BART .* maskPingrids;
maskedErrFeats_BART{17} = percentErrSlice28_BART .* maskPingrids;

% Apply mask to each feature error map for the IFT results.
% maskedErrFeats_IFT = cell(1,length(phantFeatures));
% 
% maskedErrFeats_IFT{1} = percentErrSlice4_IFT .* maskOuterCyl;  
% maskedErrFeats_IFT{2} = percentErrSlice4_IFT .* maskInnerRow1;
% maskedErrFeats_IFT{3} = percentErrSlice4_IFT .* maskInnerRow2;
% maskedErrFeats_IFT{4} = percentErrSlice4_IFT .* maskInnerRow3;
% maskedErrFeats_IFT{5} = percentErrSlice4_IFT .* maskInnerRow4;
% 
% maskedErrFeats_IFT{6}  = percentErrSlice10_IFT .* maskOuterCyl;  
% maskedErrFeats_IFT{7}  = percentErrSlice10_IFT .* maskInnerRow1;
% maskedErrFeats_IFT{8}  = percentErrSlice10_IFT .* maskInnerRow2;
% maskedErrFeats_IFT{9}  = percentErrSlice10_IFT .* maskInnerRow3;
% maskedErrFeats_IFT{10} = percentErrSlice10_IFT .* maskInnerRow4;
% 
% maskedErrFeats_IFT{11} = percentErrSlice16_IFT .* maskOuterCyl;  
% maskedErrFeats_IFT{12} = percentErrSlice16_IFT .* maskInnerRow1;
% maskedErrFeats_IFT{13} = percentErrSlice16_IFT .* maskInnerRow2;
% maskedErrFeats_IFT{14} = percentErrSlice16_IFT .* maskInnerRow3;
% maskedErrFeats_IFT{15} = percentErrSlice16_IFT .* maskInnerRow4;
% 
% maskedErrFeats_IFT{16} = percentErrSlice22_IFT .* maskPingrids;
% maskedErrFeats_IFT{17} = percentErrSlice28_IFT .* maskPingrids;

% Get the mean error information.
meanErrsBART = Get_Mean_Errors_Parallel(maskedErrFeats_BART);
%meanErrsIFT  = Get_Mean_Errors_Parallel(maskedErrFeats_IFT);    

% Get shape-specific error information.
[ outerCylErrorsSlice4_BART, innerCylErrorsSlice4_BART ] = GetRadiusDependentErrors(...
maskedErrFeats_BART{1}+maskedErrFeats_BART{2}+maskedErrFeats_BART{3}+maskedErrFeats_BART{4}+maskedErrFeats_BART{5},...
    maskOuterCyl, maskInnerRow1, maskInnerRow2, maskInnerRow3, maskInnerRow4);
% [ outerCylErrorsSlice4_IFT, innerCylErrorsSlice4_IFT ] = GetRadiusDependentErrors(...
% maskedErrFeats_IFT{1}+maskedErrFeats_IFT{2}+maskedErrFeats_IFT{3}+maskedErrFeats_IFT{4}+maskedErrFeats_IFT{5},...
%     maskOuterCyl, maskInnerRow1, maskInnerRow2, maskInnerRow3, maskInnerRow4);

[ outerCylErrorsSlice10_BART, innerCylErrorsSlice10_BART ] = GetRadiusDependentErrors(...
maskedErrFeats_BART{6}+maskedErrFeats_BART{7}+maskedErrFeats_BART{8}+maskedErrFeats_BART{9}+maskedErrFeats_BART{10},...
    maskOuterCyl, maskInnerRow1, maskInnerRow2, maskInnerRow3, maskInnerRow4);
% [ outerCylErrorsSlice10_IFT, innerCylErrorsSlice10_IFT ] = GetRadiusDependentErrors(...
% maskedErrFeats_IFT{6}+maskedErrFeats_IFT{7}+maskedErrFeats_IFT{8}+maskedErrFeats_IFT{9}+maskedErrFeats_IFT{10},...
%     maskOuterCyl, maskInnerRow1, maskInnerRow2, maskInnerRow3, maskInnerRow4);

[ outerCylErrorsSlice16_BART, innerCylErrorsSlice16_BART ] = GetRadiusDependentErrors(...
maskedErrFeats_BART{11}+maskedErrFeats_BART{12}+maskedErrFeats_BART{13}+maskedErrFeats_BART{14}+maskedErrFeats_BART{15},...
    maskOuterCyl, maskInnerRow1, maskInnerRow2, maskInnerRow3, maskInnerRow4);
% [ outerCylErrorsSlice16_IFT, innerCylErrorsSlice16_IFT ] = GetRadiusDependentErrors(...
% maskedErrFeats_IFT{11}+maskedErrFeats_IFT{12}+maskedErrFeats_IFT{13}+maskedErrFeats_IFT{14}+maskedErrFeats_IFT{15},...
%     maskOuterCyl, maskInnerRow1, maskInnerRow2, maskInnerRow3, maskInnerRow4);

gridErrorsSlice22_BART = GetGridSizeDependentErrors(maskedErrFeats_BART{16},...
    maskPingrids);
%gridErrorsSlice22_IFT = GetGridSizeDependentErrors(maskedErrFeats_IFT{16},...
%    maskPingrids);

gridErrorsSlice28_BART = GetGridSizeDependentErrors(maskedErrFeats_BART{17},...
    maskPingrids);
%gridErrorsSlice28_IFT = GetGridSizeDependentErrors(maskedErrFeats_IFT{17},...
%    maskPingrids);

% Clean up unneccesary variables
clear maskedErrFeats_IFT maskedErrFeats_BART

%% ------------------------------------------------------------------------
%  -----------------   GET AVERAGE FEATURE TIMECOURSES   ------------------
%  ------------------------------------------------------------------------

% Update user
disp(strcat('SIMULATION: Getting average feature timecourses. (',...
    datestr(datetime('now')), ')' ))

% Get average timecourses of feature 1 through 5.
aveTimeCourseFeat1BART = Get_Mean_Timecourse_Parallel(BARTrecons, 4,...
    maskOuterCyl);
aveTimeCourseFeat2BART = Get_Mean_Timecourse_Parallel(BARTrecons, 4,...
    maskInnerRow1);
aveTimeCourseFeat3BART = Get_Mean_Timecourse_Parallel(BARTrecons, 4,...
    maskInnerRow2);
aveTimeCourseFeat4BART = Get_Mean_Timecourse_Parallel(BARTrecons, 4,...
    maskInnerRow3);
aveTimeCourseFeat5BART = Get_Mean_Timecourse_Parallel(BARTrecons, 4,...
    maskInnerRow4);

% aveTimeCourseFeat1IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 4,...
%     maskOuterCyl);
% aveTimeCourseFeat2IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 4,...
%     maskInnerRow1);
% aveTimeCourseFeat3IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 4,...
%     maskInnerRow2);
% aveTimeCourseFeat4IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 4,...
%     maskInnerRow3);
% aveTimeCourseFeat5IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 4,...
%     maskInnerRow4);

% Get average timecourses of feature 6 through 10.
aveTimeCourseFeat6BART = Get_Mean_Timecourse_Parallel(BARTrecons, 10,...
    maskOuterCyl);
aveTimeCourseFeat7BART = Get_Mean_Timecourse_Parallel(BARTrecons, 10,...
    maskInnerRow1);
aveTimeCourseFeat8BART = Get_Mean_Timecourse_Parallel(BARTrecons, 10,...
    maskInnerRow2);
aveTimeCourseFeat9BART = Get_Mean_Timecourse_Parallel(BARTrecons, 10,...
    maskInnerRow3);
aveTimeCourseFeat10BART = Get_Mean_Timecourse_Parallel(BARTrecons, 10,...
    maskInnerRow4);

% aveTimeCourseFeat6IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 10,...
%     maskOuterCyl);
% aveTimeCourseFeat7IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 10,...
%     maskInnerRow1);
% aveTimeCourseFeat8IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 10,...
%     maskInnerRow2);
% aveTimeCourseFeat9IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 10,...
%     maskInnerRow3);
% aveTimeCourseFeat10IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 10,...
%     maskInnerRow4);

% Get average timecourses of feature 10 through 15.
aveTimeCourseFeat11BART = Get_Mean_Timecourse_Parallel(BARTrecons, 16,...
    maskOuterCyl);
aveTimeCourseFeat12BART = Get_Mean_Timecourse_Parallel(BARTrecons, 16,...
    maskInnerRow1);
aveTimeCourseFeat13BART = Get_Mean_Timecourse_Parallel(BARTrecons, 16,...
    maskInnerRow2);
aveTimeCourseFeat14BART = Get_Mean_Timecourse_Parallel(BARTrecons, 16,...
    maskInnerRow3);
aveTimeCourseFeat15BART = Get_Mean_Timecourse_Parallel(BARTrecons, 16,...
    maskInnerRow4);

% aveTimeCourseFeat11IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 16,...
%     maskOuterCyl);
% aveTimeCourseFeat12IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 16,...
%     maskInnerRow1);
% aveTimeCourseFeat13IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 16,...
%     maskInnerRow2);
% aveTimeCourseFeat14IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 16,...
%     maskInnerRow3);
% aveTimeCourseFeat15IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 16,...
%     maskInnerRow4);

% Get average timecourses of feature 16 and 17.
aveTimeCourseFeat16BART = Get_Mean_Timecourse_Parallel(BARTrecons, 22,...
    maskPingrids);
aveTimeCourseFeat17BART = Get_Mean_Timecourse_Parallel(BARTrecons, 28,...
    maskPingrids);

% aveTimeCourseFeat16IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 22,...
%     maskPingrids);
% aveTimeCourseFeat17IFT = Get_Mean_Timecourse_Parallel(IFTrecons, 28,...
%     maskPingrids);

% Clear unnecessary variables.
clear maskOuterCyl maskInnerRow1 maskInnerRow2 maskInnerRow3 maskInnerRow4 ...
      maskPingrids 

%% ------------------------------------------------------------------------
%  --------------------   CLEAN UP FOR END OF EXECUTION   -----------------
%  ------------------------------------------------------------------------

% Update user
disp(strcat('SIMULATION: Ending execution, saving data. (',...
    datestr(datetime('now')), ')' ))

% Save workspace results.
save([fileStr '/IntermediateResults.mat'],'-v7.3')

%delete(POOL);
toc

