function [v1p, scp] = ...
    areaOverKernel_GrandAv_MOTOR_iter(bootSamps, dp_cut, nSigma, kernType)
% Iteratively Compute Area Over Kernel for Group Data
% Over the entire -600 to 200 ms window shown in figures for motor aligned 
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

% kernType
% Can take either 'FA' or 'RT' to generate data for 

% Analysis window is 100 ms before lever release.
%% Load Table and Directory for Stim Profiles
% Location of the Final Stim Profiles and Tables
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';
% Starts with V1
v1Folder = strcat(analysisDir, 'V1',' FA/');
scFolder = strcat(analysisDir, 'SC',' FA/');
% Load Tables for each stimulus type
load([v1Folder 'masterTable.mat'], 'T');
v1T = T;
load([scFolder 'masterTable.mat'], 'T');
scT = T;
clear T;
%% Grab all the stim Profiles one mouse then create kernel
% Use Limits to Subselect Sessions From Each Mouse
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;

% Init
analysisDurMS    = 25; 
analysisStartBin = 26;
analysisEndBin   = 775; 
%% Get V1 Traces, Subselect by deltad' if dp_cut == 1
% First For Luminance
limits.animal = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
U = selectUsingLimits(v1T, limits);
condition = ['V1' ' ' 'FA'];

% Fitted Delta D' Distribution for V1 Data Combined
if dp_cut == 1
    mu =  -0.006;
    sigma =  0.2584;
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
v1Profiles.hitProfiles = [];
v1Profiles.missProfiles = [];
v1Profiles.earlyProfiles = [];
v1Profiles.RTProfiles = [];
v1Profiles.stimRTProfiles = [];

% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
    v1Profiles.hitProfiles = [v1Profiles.hitProfiles; stimProfiles.hitProfiles];
    v1Profiles.missProfiles = [v1Profiles.missProfiles; stimProfiles.missProfiles];
    v1Profiles.earlyProfiles = [v1Profiles.earlyProfiles; stimProfiles.earlyProfiles];
    v1Profiles.RTProfiles = [v1Profiles.RTProfiles; stimProfiles.RTProfiles];
    v1Profiles.stimRTProfiles = [v1Profiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end


%% Repeat For SC
limits.animal = {'1458', '1548', '1674', '1675', '1902', '1905', '2057', '2058', '2063', '2169', '2236'};
U = selectUsingLimits(scT, limits);
condition = ['SC' ' ' 'FA'];

% Fitted Delta D' Distribution for SC Data Combined
if dp_cut == 1
    mu =  -0.023;
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
scProfiles.hitProfiles = [];
scProfiles.missProfiles = [];
scProfiles.earlyProfiles = [];
scProfiles.RTProfiles = [];
scProfiles.stimRTProfiles = [];

% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
    scProfiles.hitProfiles = [scProfiles.hitProfiles; stimProfiles.hitProfiles];
    scProfiles.missProfiles = [scProfiles.missProfiles; stimProfiles.missProfiles];
    scProfiles.earlyProfiles = [scProfiles.earlyProfiles; stimProfiles.earlyProfiles];
    scProfiles.RTProfiles = [scProfiles.RTProfiles; stimProfiles.RTProfiles];
    scProfiles.stimRTProfiles = [scProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end

%% Collate All Traces

if strcmp(kernType, 'FA')
    scKernel = [scProfiles.earlyProfiles];
    v1Kernel = [v1Profiles.earlyProfiles];
elseif strcmp(kernType, 'RT')
    scKernel = [scProfiles.RTProfiles];
    v1Kernel = [v1Profiles.RTProfiles];
end

% Number of Outcomes
v1_nTotal = size(v1Kernel,1);
sc_nTotal = size(scKernel,1);
%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)

% Init Archives For BootStrap Samples
BootsAOK_V1 = zeros(bootSamps, 800);
BootsAOK_SC = zeros(bootSamps, 800);

% Bin to Start Computing AOK
lookBack = floor(analysisDurMS/2);
startBin = analysisStartBin-lookBack;

for binNum = analysisStartBin:analysisEndBin
    for bootNum = 1:bootSamps
        % Samps For This Round of BootStrapping
        v1Samps = randsample(v1_nTotal,v1_nTotal,true)';
        scSamps = randsample(sc_nTotal,sc_nTotal,true)';

        % Take Samples w/ replacement
        Boot_V1 = v1Kernel(v1Samps,:);
        Boot_SC = scKernel(scSamps,:);

        BootsAOK_V1(bootNum,binNum) = sum(mean(-1*Boot_V1(:,startBin:startBin+analysisDurMS-1)));
        BootsAOK_SC(bootNum,binNum) = sum(mean(-1*Boot_SC(:,startBin:startBin+analysisDurMS-1)));
    end
    startBin = startBin + 1;
end
%% Evaluate Results
scp = zeros(1, 800);
v1p = zeros(1, 800);

for binNum = analysisStartBin:analysisEndBin
    scp(1,binNum) = (size(BootsAOK_SC,1) - sum(BootsAOK_SC(:,binNum)>0))/size(BootsAOK_SC,1);
    v1p(1,binNum) = (size(BootsAOK_V1,1) - sum(BootsAOK_V1(:,binNum)>0))/size(BootsAOK_V1,1);
end

end
