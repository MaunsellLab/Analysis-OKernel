function [v1p_lum, scp_lum, v1p_gab, scp_gab] = ...
    areaOverKernel_GrandAv_iter(bootSamps, dp_cut, nSigma)
% Iteratively Compute Area Over Kernel for Group Data
% Over the entire -400 to 400 ms window shown in figures. 
% Bootstrap to do the test for each bin
% Spit out significant or not for each bin

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
lumFolder = strcat(analysisDir, 'V1',' Lum/');
gabFolder = strcat(analysisDir, 'V1',' Gabor/');

% Load Tables for each stimulus type
load([lumFolder 'masterTable.mat'], 'T');
lumT = T;
load([gabFolder 'masterTable.mat'], 'T');
gabT = T;
clear T;
%% Basic SetUp

limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;

analysisDurMS    = 25; % Was 100
analysisStartBin = 26;
analysisEndBin   = 775; 
%% V1 Luminance

% First For Luminance
limits.animal = {'1960', '2015', '2083', '2126', '2220', '2221'};
U = selectUsingLimits(lumT, limits);
condition = ['V1' ' ' 'Lum'];

% Fitted Delta D' Distribution for V1 Luminance
if dp_cut == 1
    mu =  -0.011;
    sigma =  0.2326;
    deltaD = [U.stimDPrime]-[U.noStimDPrime];
    bins = -1:0.1:1;
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
    nDropped = height(U) - sum(keepIdx);
else
     keepIdx = ones(height(U),1);
end

% Clean Up U
U = U(keepIdx,:);

% Extract optoProfiles
lumProfiles.hitProfiles = [];
lumProfiles.missProfiles = [];
lumProfiles.earlyProfiles = [];
lumProfiles.RTProfiles = [];
lumProfiles.stimRTProfiles = [];

% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
    lumProfiles.hitProfiles = [lumProfiles.hitProfiles; stimProfiles.hitProfiles];
    lumProfiles.missProfiles = [lumProfiles.missProfiles; stimProfiles.missProfiles];
    lumProfiles.earlyProfiles = [lumProfiles.earlyProfiles; stimProfiles.earlyProfiles];
    lumProfiles.RTProfiles = [lumProfiles.RTProfiles; stimProfiles.RTProfiles];
    lumProfiles.stimRTProfiles = [lumProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end

%% Repeat For V1 Gabor
limits.animal = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
U = selectUsingLimits(gabT, limits);
condition = ['V1' ' ' 'Gabor'];

% Fitted Delta D' Distribution for V1 Gabor
if dp_cut == 1
    mu =  -0.0034;
    sigma = 0.277;
    deltaD = [U.stimDPrime]-[U.noStimDPrime];
    bins = -1:0.1:1;
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
    nDropped = height(U) - sum(keepIdx);
else
     keepIdx = ones(height(U),1);
end

% Clean Up Table
U = U(keepIdx,:);

% Extract optoProfiles
gabProfiles.hitProfiles = [];
gabProfiles.missProfiles = [];
gabProfiles.earlyProfiles = [];
gabProfiles.RTProfiles = [];
gabProfiles.stimRTProfiles = [];

% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
    gabProfiles.hitProfiles = [gabProfiles.hitProfiles; stimProfiles.hitProfiles];
    gabProfiles.missProfiles = [gabProfiles.missProfiles; stimProfiles.missProfiles];
    gabProfiles.earlyProfiles = [gabProfiles.earlyProfiles; stimProfiles.earlyProfiles];
    gabProfiles.RTProfiles = [gabProfiles.RTProfiles; stimProfiles.RTProfiles];
    gabProfiles.stimRTProfiles = [gabProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end

%% Collate Traces and compute n's

% Collate All Traces 
gabKernel = [gabProfiles.hitProfiles; -gabProfiles.missProfiles];
lumKernel = [lumProfiles.hitProfiles; -lumProfiles.missProfiles];

% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
lum_nHit = size(lumProfiles.hitProfiles,1);
gab_nHit = size(gabProfiles.hitProfiles,1);
lum_nTotal = size(lumKernel,1);
gab_nTotal = size(gabKernel,1);


%% Bootstrap over each bin

% Init Archives For BootStrap Samples
lumBootsAOK_V1 = zeros(bootSamps, 800);
gabBootsAOK_V1 = zeros(bootSamps, 800);

% Bin to Start Computing AOK
startBin = analysisStartBin;

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
        % Advance Start Bin
    end
    startBin = startBin+1;
end

%% Find Bins with Significant AOK

v1p_lum = zeros(1, 800);
v1p_gab = zeros(1, 800);

for binNum = analysisStartBin:analysisEndBin
    v1p_lum(1,binNum) = (size(lumBootsAOK_V1,1) - sum(lumBootsAOK_V1(:,binNum)>0))/size(lumBootsAOK_V1,1);
    v1p_gab(1,binNum) = (size(gabBootsAOK_V1,1) - sum(gabBootsAOK_V1(:,binNum)>0))/size(gabBootsAOK_V1,1);
end

%% Repeat for SC
lumFolder = strcat(analysisDir, 'SC',' Lum/');
gabFolder = strcat(analysisDir, 'SC',' Gabor/');
clear lumT gabT;

% Load Tables for each stimulus type
load([lumFolder 'masterTable.mat'], 'T');
lumT = T;
load([gabFolder 'masterTable.mat'], 'T');
gabT = T;
clear T;
%% Both Gabor and Luminance
% First For Luminance
limits.animal = {'1458', '1548', '1674', '1675', '1902', '1905', '2057', '2058', '2063', '2169', '2236'};
U = selectUsingLimits(lumT, limits);
condition = ['SC' ' ' 'Lum'];

% Fitted Delta D' Distribution for SC Luminance
if dp_cut == 1
    mu = -0.03914;
    sigma = 0.2715;
    deltaD = [U.stimDPrime]-[U.noStimDPrime];
    bins = -1:0.1:1;
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
    nDropped = height(U) - sum(keepIdx);
else
    keepIdx = ones(height(U),1);
end

% Clean Up Table
U = U(keepIdx,:);

% Extract optoProfiles
lumProfiles.hitProfiles = [];
lumProfiles.missProfiles = [];
lumProfiles.earlyProfiles = [];
lumProfiles.RTProfiles = [];
lumProfiles.stimRTProfiles = [];

% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
    lumProfiles.hitProfiles = [lumProfiles.hitProfiles; stimProfiles.hitProfiles];
    lumProfiles.missProfiles = [lumProfiles.missProfiles; stimProfiles.missProfiles];
    lumProfiles.earlyProfiles = [lumProfiles.earlyProfiles; stimProfiles.earlyProfiles];
    lumProfiles.RTProfiles = [lumProfiles.RTProfiles; stimProfiles.RTProfiles];
    lumProfiles.stimRTProfiles = [lumProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end

% Repeat For Gabor
limits.animal = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
U = selectUsingLimits(gabT, limits);
condition = ['SC' ' ' 'Gabor'];

% Fitted Delta D' Distribution for SC Gabor
if dp_cut == 1
    mu = 0.0056;
    sigma = 0.287;
    deltaD = [U.stimDPrime]-[U.noStimDPrime];
    bins = -1:0.1:1;
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
    nDropped = height(U) - sum(keepIdx);
else
    keepIdx = ones(height(U),1);
end

% Clean Up Table
U = U(keepIdx,:);

% Extract optoProfiles
gabProfiles.hitProfiles = [];
gabProfiles.missProfiles = [];
gabProfiles.earlyProfiles = [];
gabProfiles.RTProfiles = [];
gabProfiles.stimRTProfiles = [];

% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
    gabProfiles.hitProfiles = [gabProfiles.hitProfiles; stimProfiles.hitProfiles];
    gabProfiles.missProfiles = [gabProfiles.missProfiles; stimProfiles.missProfiles];
    gabProfiles.earlyProfiles = [gabProfiles.earlyProfiles; stimProfiles.earlyProfiles];
    gabProfiles.RTProfiles = [gabProfiles.RTProfiles; stimProfiles.RTProfiles];
    gabProfiles.stimRTProfiles = [gabProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end

%% Collate Traces and n's

% Collate All Traces
gabKernel = [gabProfiles.hitProfiles; -gabProfiles.missProfiles];
lumKernel = [lumProfiles.hitProfiles; -lumProfiles.missProfiles];

% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
lum_nHit = size(lumProfiles.hitProfiles,1);
gab_nHit = size(gabProfiles.hitProfiles,1);
lum_nTotal = size(lumKernel,1);
gab_nTotal = size(gabKernel,1);

%% Bootstrap over each bin


% Init Archives For BootStrap Samples
lumBootsAOK_SC = zeros(bootSamps, 800);
gabBootsAOK_SC = zeros(bootSamps, 800);

startBin = analysisStartBin;

for binNum = analysisStartBin:analysisEndBin
    for bootNum = 1:bootSamps
        % Samps For This Round of BootStrapping
        lumSamps = [randsample(lum_nHit,lum_nHit,true)'...
            randsample([lum_nHit+1:lum_nTotal],lum_nTotal-lum_nHit,true)]';
        gabSamps = [randsample(gab_nHit,gab_nHit,true)'...
            randsample([gab_nHit+1:gab_nTotal],gab_nTotal - gab_nHit,true)]';
        % Take Samples w/ replacement
        lumBoot_SC = lumKernel(lumSamps,:);
        gabBoot_SC = gabKernel(gabSamps,:);
        
        % Compute AOK at this Bin
        lumBootsAOK_SC(bootNum,binNum) = sum(mean(-1*lumBoot_SC(:,startBin:startBin+analysisDurMS-1)));
        gabBootsAOK_SC(bootNum,binNum) = sum(mean(-1*gabBoot_SC(:,startBin:startBin+analysisDurMS-1)));
        % Advance Start Bin
    end
    startBin = startBin+1; % Advance Analysis Window
end

%% Evaluate Results over each Bin of the Bootstrap
scp_lum = zeros(1, 800);
scp_gab = zeros(1, 800);

for binNum = analysisStartBin:analysisEndBin
    scp_lum(1,binNum) = (size(lumBootsAOK_SC,1) - sum(lumBootsAOK_SC(:,binNum)>0))/size(lumBootsAOK_SC,1);
    scp_gab(1,binNum) = (size(gabBootsAOK_SC,1) - sum(gabBootsAOK_SC(:,binNum)>0))/size(gabBootsAOK_SC,1);
end


