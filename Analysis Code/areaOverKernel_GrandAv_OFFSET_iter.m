function [v1p_off, scp_off] = ...
    areaOverKernel_GrandAv_OFFSET_iter(bootSamps)
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


%% Load Table and Directory for Stim Profiles

% Location of the Final Stim Profiles and Tables
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';

% Starts with V1
V1Folder = strcat(analysisDir, 'V1',' Offset/');
SCFolder = strcat(analysisDir, 'SC',' Offset/');

% Load Tables for each stimulus type
load([V1Folder 'masterTable.mat'], 'T');
v1T = T;
load([SCFolder 'masterTable.mat'], 'T');
scT = T;
clear T;

%% Grab all the stim Profiles one mouse then create kernel
% Use Limits to Subselect Sessions From Each Mouse
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;

analysisDurMS    = 25; % Was 100
analysisStartBin = 26;
analysisEndBin   = 775; 
%% V1 Offset
limits.animal = {'2016', '2083', '2126', '2220', '2221'};
U = selectUsingLimits(v1T, limits);
condition = ['V1' ' ' 'Offset'];

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

% total Kernel
v1Kernel = [v1Profiles.hitProfiles; -v1Profiles.missProfiles];

% Compute Mean AOK during the analysis Window

%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
nHit = size(v1Profiles.hitProfiles,1);
v1_nTotal = size(v1Kernel,1);

% Init Archives For BootStrap Samples
BootsAOK_V1 = zeros(bootSamps, 800);
lookBack = floor(analysisDurMS/2);
startBin = analysisStartBin-lookBack;

for binNum = analysisStartBin:analysisEndBin;
    for bootNum = 1:bootSamps
        % Samps For This Round of BootStrapping
        v1Samps = [randsample(nHit,nHit,true)'...
            randsample([nHit+1:v1_nTotal],v1_nTotal-nHit,true)]';
        % Take Samples w/ replacement
        Boot_V1 = v1Kernel(v1Samps,:);


        BootsAOK_V1(bootNum,binNum) = sum(mean(-1*Boot_V1(:,startBin:startBin+analysisDurMS-1)));
    end
    startBin = startBin+1;
end



%% Repeat for SC
limits.animal = {'1674', '1675', '1902', '2057', '2063', '2236'};
SCFolder = strcat(analysisDir, 'SC',' Offset/');
load([SCFolder 'masterTable.mat'], 'T');
scT = T;
clear T;

%% Both Gabor and Luminance
% First For Luminance
U = selectUsingLimits(scT, limits);
condition = ['SC' ' ' 'Offset'];

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

% Collate All Traces From This Mouse
scKernel = [scProfiles.hitProfiles; -scProfiles.missProfiles];

%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
nHit = size(scProfiles.hitProfiles,1);
sc_nTotal = size(scKernel,1);

% Init Archives For BootStrap Samples
BootsAOK_SC = zeros(bootSamps, 800);
lookBack = floor(analysisDurMS/2);
startBin = analysisStartBin-lookBack;

for binNum = analysisStartBin:analysisEndBin
    for bootNum = 1:bootSamps
        % Samps For This Round of BootStrapping
        scSamps = [randsample(nHit,nHit,true)'...
            randsample([nHit+1:sc_nTotal],sc_nTotal-nHit,true)]';
        % Take Samples w/ replacement
        Boot_SC = scKernel(scSamps,:);

        BootsAOK_SC(bootNum,binNum) = sum(mean(-1*Boot_SC(:,startBin:startBin+analysisDurMS-1)));
    end
    startBin = startBin+1; % Advance Analysis Window
end
%% Evaluate Results over each Bin of the Bootstrap
scp_off = zeros(1, 800);
v1p_off = zeros(1, 800);

for binNum = analysisStartBin:analysisEndBin
    v1p_off(1,binNum) = (size(BootsAOK_V1,1) - sum(BootsAOK_V1(:,binNum)>0))/size(BootsAOK_V1,1);
    scp_off(1,binNum) = (size(BootsAOK_SC,1) - sum(BootsAOK_SC(:,binNum)>0))/size(BootsAOK_SC,1);
end


end

