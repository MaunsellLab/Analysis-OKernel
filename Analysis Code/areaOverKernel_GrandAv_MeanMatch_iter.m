function [lumP, gabP] = ...
    areaOverKernel_GrandAv_MeanMatch_iter(bootSamps)
% Iteratively Compute Area Over Kernel for Group Data
% Over the entire -400 to 400 ms window shown in figures for stimulus aligned 
% Bootstrap to do the test for each bin
% Spit out significant or not for each bin
% Only Consider Animals For Which Data was obtained for both stimuli

% Inputs:
% invert = 0,1 whether to invert the AOK from negative to positive
% bootSamps = How many bootstraps to create CIs
% analysisStartMS = relative to stimOnset (time = 0)
% analysisDurMS = How many bins (ms) to compute the deflection in the
% kernel after the start of analysis
% dp_cut, true false, whether to sub-select sessions based on delta dprime
% (effect of LED).
% nSigma = number of sigma on delta d' distribution to drop. 


%% Load Table and Directory for Stim Profiles

% Location of the Final Stim Profiles and Tables
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';

% Starts with V1
lumFolder = strcat(analysisDir, 'V1',' Lum MM/'); % Mean Matched
gabFolder = strcat(analysisDir, 'V1',' Gabor MM/');

% Load Tables for each stimulus type
load([lumFolder 'masterTable.mat'], 'V1Lum_DDP');
lumT = V1Lum_DDP;
load([gabFolder 'masterTable.mat'], 'V1Gab_DDP');
gabT = V1Gab_DDP;
clear V1Gab_DDP V1Lum_DDP;
%% Grab all the stim Profiles one mouse then create kernel
% Use Limits to Subselect Sessions From Each Mouse
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;

% Init
analysisDurMS    = 25; 
analysisStartBin = 26;
analysisEndBin   = 775; 

%% Both Gabor and Luminance

load('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Lum MM/stimProfiles/lumProfiles_DDP.mat');

lumProfiles.hitProfiles = [lumProfiles_DDP.hitProfiles];
lumProfiles.missProfiles = [lumProfiles_DDP.missProfiles];
lumProfiles.earlyProfiles = [lumProfiles_DDP.earlyProfiles];
lumProfiles.RTProfiles = [lumProfiles_DDP.RTProfiles];
lumProfiles.stimRTProfiles = [lumProfiles_DDP.stimRTProfiles];
clear lumProfiles_DDP;

% Luminance Kernels
lumKernel = [lumProfiles.hitProfiles; -lumProfiles.missProfiles];

load('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Gabor MM/stimProfiles/gabProfiles_DDP.mat');
gabProfiles.hitProfiles = [gabProfiles_DDP.hitProfiles];
gabProfiles.missProfiles = [gabProfiles_DDP.missProfiles];
gabProfiles.earlyProfiles = [gabProfiles_DDP.earlyProfiles];
gabProfiles.RTProfiles = [gabProfiles_DDP.RTProfiles];
gabProfiles.stimRTProfiles = [gabProfiles_DDP.stimRTProfiles];
clear gabProfiles_DDP;

% Gabor Kernels
gabKernel = [gabProfiles.hitProfiles; -gabProfiles.missProfiles];

%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
lum_nHit = size(lumProfiles.hitProfiles,1);
gab_nHit = size(gabProfiles.hitProfiles,1);
lum_nTotal = size(lumKernel,1);
gab_nTotal = size(gabKernel,1);

% Init Archives For BootStrap Samples
lumBootsAOK_V1 = zeros(bootSamps, 800);
gabBootsAOK_V1 = zeros(bootSamps, 800);

% Bin to Start Computing AOK
lookBack = floor(analysisDurMS/2);
startBin = analysisStartBin-lookBack;

for binNum = analysisStartBin:analysisEndBin
    for bootNum = 1:bootSamps
        % Samps For This Round of BootStrapping
        lumSamps = [randsample(lum_nHit,lum_nHit,true)'...
            randsample([lum_nHit+1:lum_nTotal],lum_nTotal-lum_nHit,true)]';
        gabSamps = [randsample(gab_nHit,gab_nHit,true)'...
            randsample([gab_nHit+1:gab_nTotal],gab_nTotal - gab_nHit,true)]';
        % Take Samples w/ replacement
        lumBoot_V1 = lumKernel(lumSamps,:);
        gabBoot_V1 = gabKernel(gabSamps,:);

        lumBootsAOK_V1(bootNum,binNum) = sum(mean(-1*lumBoot_V1(:,startBin:startBin+analysisDurMS-1)));
        gabBootsAOK_V1(bootNum,binNum) = sum(mean(-1*gabBoot_V1(:,startBin:startBin+analysisDurMS-1)));
    end
    startBin = startBin+1;
end

%% %% Find Bins with Significant AOK

lumP = zeros(1, 800);
gabP = zeros(1, 800);

for binNum = analysisStartBin:analysisEndBin
    lumP(1,binNum) = (size(lumBootsAOK_V1,1) - sum(lumBootsAOK_V1(:,binNum)>0))/size(lumBootsAOK_V1,1);
    gabP(1,binNum) = (size(gabBootsAOK_V1,1) - sum(gabBootsAOK_V1(:,binNum)>0))/size(gabBootsAOK_V1,1);
end
