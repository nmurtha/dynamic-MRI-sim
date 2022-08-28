function [CIRCUScalib, CIRCUSquanta] = Generate_CIRCUS_PhaseTables(...
    phaseEncodeDim1, phaseEncodeDim2,...
    CIRCUSb, CIRCUSc,...
    CIRCUSmode, CIRCUSreturn,...
    CIRCUSseed, numEffQuanta)

%% ------------------------------------------------------------------------
%  --------------------   FUNCTION INFORMATION   --------------------------
%  ------------------------------------------------------------------------
% 
%  WRITTEN: Nathan Murtha (nathan.j.murtha@gmail.com)
%           Jan 25, 2017
%
%  PURPOSE: Generates phase tables for CIRCUS k-space sampling patterns.
%           Intent is to seperately generate calibration region and
%           "effective" quanta patterns composed of a user defined
%           number of true single quanta. This allows flexible acquisition
%           and interleaving of k-space data later. 
%
%           See "Accelerated MRI with CIRcular Cartesian UnderSampling
%           (CIRCUS): a variable density Cartesian sampling Strategy for
%           compressed sensing and parallel imaging" by Jing Liu and David
%           Saloner.
%
%  INPUTS:
%   phaseEncodeDim1 -> Dimension of first phase encode direction (e.g. the
%                      number of rows in the image array). A scalar input.
%
%   phaseEncodeDim2 -> Dimension of second phase encode direction (e.g. the
%                      number of slices in the image array). A scalar 
%                      input.
%                      NOTE: Required that phaseEncodeDim2 <
%                            phaseEncodeDim1.
%
%   CIRCUSb         -> The CIRCUS radial parameter, b. Intuitively, it is 
%                      associated with the "shear" in CIRCUS quanta. A
%                      scalar input.
%
%   CIRCUSc         -> The CIRCUS spiral parameter, c. Intuitively, it is 
%                      associated with the "twist" in CIRCUS quanta. A
%                      scalar input.
%
%   CIRCUSmode      -> 1: Return specified number of quanta.
%                      2: Generate quanta until specified
%                      undersampling factor is reached. A scalar input.
%
%   CIRCUSreturn    -> If mode 1 is used, represents the number of quanta
%                      that will be returned in the effective patterns.
%                      If mode 2 is used, represents the goal undersampling
%                      factor that is aimed for in the effective patterns.
%                      A scalar input.
%
%   CIRCUSseed      -> The seed for quanta generation. Seeds < 0 will
%                      cause the return of random quanta. Seeds >= 0 will
%                      return sequential quanta. A scalar input.
%                      NOTE: Due to the variable handling in
%                            CIRCUSPattern.m, a seed of "n" will cause the 
%                            return of the "n+1st" quanta. E.g. seed 0 
%                            gives quanta 1.
%
%   numEffQuanta    -> Number of "effective" quanta to be generated. A 
%                      scalar input. 
%                      NOTE: Effective quanta refers to the fact that these
%                           patterns may be a composite of several single 
%                           quanta. Each effective quanta has as many 
%                           single quanta as the variables CIRCUSmode and 
%                           CIRCUSreturn dictate. 
% 
%  OUTPUTS: 
%   CIRCUScalib    -> Structure array output containing the phase tables 
%                     for the calibration region pattern. Contains two
%                     vectors:
%                       > CIRCUScalib.phaseTableDim1
%                       > CIRCUScalib.phaseTableDim2
%                     Since all calibration regions are equivalent
%                     trajectories, only one set of phase tables is made in
%                     CIRCUScalib.
%
%   CIRCUSquanta    -> Structure array output containing the phase tables 
%                      for the effective CIRCUS quanta. Contains two 
%                      vectors:
%                       > CIRCUSquanta.phaseTableDim1
%                       > CIRCUSquanta.phaseTableDim2
%                      Has as many phase tables as numEffQuanta.
%
%  EXTERNALS:
%   CIRCUSPattern.m (written James Rioux, edited Nathan Murtha)

%% ------------------------------------------------------------------------
%  --------------   GENERATE CALIBRATION REGION PHASE TABLES   ------------
%  ------------------------------------------------------------------------

% Allocate structure array for phase tables.
CIRCUScalib = struct('phaseTableDim1',[],...
                     'phaseTableDim2',[]);

% Acquire the phase tables.
[CIRCUScalib.phaseTableDim1,CIRCUScalib.phaseTableDim2,~,~] = ...
    CIRCUSPattern(phaseEncodeDim1, phaseEncodeDim2, ...
                  CIRCUSb, CIRCUSc, ...
                  1, 0, ...
                  0, 1, 0);
                 
%% ------------------------------------------------------------------------
%  ---------------   GENERATE EFFECTIVE QUANTA PHASE TABLES   -------------
%  ------------------------------------------------------------------------                 
              
% Allocate structure array for phase tables.
CIRCUSquanta = struct('phaseTableDim1',[],...
                      'phaseTableDim2',[]);
                  
% Generate the user requested number of effective quanta. Note that the
% seed for quanta generation is updated with each call to CIRCUSPattern.
for iPattern = 1:numEffQuanta
    [CIRCUSquanta(iPattern).phaseTableDim1,...
     CIRCUSquanta(iPattern).phaseTableDim2,~,CIRCUSseed] =... % seed updated in every iteration to advance quanta generation
        CIRCUSPattern(phaseEncodeDim1, phaseEncodeDim2, ...
                      CIRCUSb, CIRCUSc, ...
                      CIRCUSmode, CIRCUSreturn,...
                      CIRCUSseed,0,0);
end
              
end





