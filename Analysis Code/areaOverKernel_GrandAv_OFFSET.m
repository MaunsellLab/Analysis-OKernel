function [BootsAOK_V1, BootsAOK_SC] = ...
    areaOverKernel_GrandAv_OFFSET(bootSamps, analysisStartMS, analysisDurMS)
% Compute Area Over Kernel for Group Data
% Bootstrap 95%CI per mouse
% Plot Scatter of results
% Only Consider Animals For Which Data was obtained for both stimuli

% Inputs:
% invert = 0,1 whether to invert the AOK from negative to positive
% bootSamps = How many bootstraps to create CIs
% analysisStartMS = relative to stimOnset (time = 0)
% analysisDurMS = How many bins (ms) to compute the deflection in the
% kernel after the start of analysis

analysisStartBin = 401+analysisStartMS; % Stim on is at bin 401
analysisEndBin = analysisStartBin + analysisDurMS;
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

% Init
AOK_V1 = zeros(1,1);
AOK_V1_95CI  = zeros(1,2);
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

AOK_V1(1,1) = sum(mean(-1*v1Kernel(:,analysisStartBin:analysisEndBin)));

%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
nHit = size(v1Profiles.hitProfiles,1);
v1_nTotal = size(v1Kernel,1);

% Init Archives For BootStrap Samples
BootsAOK_V1 = zeros(bootSamps, 1);

for bootNum = 1:bootSamps
    % Samps For This Round of BootStrapping
    v1Samps = [randsample(nHit,nHit,true)'...
        randsample([nHit+1:v1_nTotal],v1_nTotal-nHit,true)]';
    % Take Samples w/ replacement
    Boot_V1 = v1Kernel(v1Samps,:);


    BootsAOK_V1(bootNum,1) = sum(mean(-1*Boot_V1(:,analysisStartBin:analysisEndBin)));
end

AOK_V1_95CI(1,:) = quantile(BootsAOK_V1,[0.025 0.975]);

%% Repeat for SC
limits.animal = {'1674', '1675', '1902', '2057', '2063', '2236'};
SCFolder = strcat(analysisDir, 'SC',' Offset/');
load([SCFolder 'masterTable.mat'], 'T');
scT = T;
clear T;

% Init
AOK_SC = zeros(1,1);
AOK_SC_95CI  = zeros(1,2);

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


AOK_SC(1,1) = sum(mean(-1*scKernel(:,analysisStartBin:analysisEndBin)));

%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
nHit = size(scProfiles.hitProfiles,1);
sc_nTotal = size(scKernel,1);

% Init Archives For BootStrap Samples
BootsAOK_SC = zeros(bootSamps, 1);

for bootNum = 1:bootSamps
    % Samps For This Round of BootStrapping
    scSamps = [randsample(nHit,nHit,true)'...
        randsample([nHit+1:sc_nTotal],sc_nTotal-nHit,true)]';
    % Take Samples w/ replacement
    Boot_SC = scKernel(scSamps,:);
    
    BootsAOK_SC(bootNum,1) = sum(mean(-1*Boot_SC(:,analysisStartBin:analysisEndBin)));


end

AOK_SC_95CI(1,:) = quantile(BootsAOK_SC,[0.025 0.975]);

%% Make Scatter Plot of Results

V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];

%% Histograms
% Colors
V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];
bins = -3:0.2:3;
% Figure
figure('Position', [10 10 1000 500]);
axis square;
hold on;
histogram(BootsAOK_V1,bins, 'FaceColor', V1Color, 'FaceAlpha', 0.5, 'Normalization','probability');
histogram(BootsAOK_SC,bins, 'FaceColor', SCColor, 'FaceAlpha', 0.5, 'Normalization','probability');
plot([0 0], [0 max(ylim)], 'Color', 'k', 'LineStyle', '--', 'LineWidth',1);
plot([median(BootsAOK_V1) median(BootsAOK_V1)], [0 max(ylim)], 'Color', V1Color, 'LineStyle', '--', 'LineWidth',1);
plot([median(BootsAOK_SC) median(BootsAOK_SC)], [0 max(ylim)], 'Color', SCColor, 'LineStyle', '--', 'LineWidth',1);
title('Gabor: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
ylabel('Probability');
xlabel('Area Over The Kernel');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
ylim([0 0.35]);
xlim([-2.5 2.5]);
legend('V1', 'SC', 'Location','northwest');
hold off;

% Bootstraps Different from 0?
v1p = (length(BootsAOK_V1) - sum(BootsAOK_V1>0))/length(BootsAOK_V1)
scp = (length(BootsAOK_SC) - sum(BootsAOK_SC>0))/length(BootsAOK_SC)


end

