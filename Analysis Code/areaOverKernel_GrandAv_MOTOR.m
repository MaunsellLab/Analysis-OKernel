function [BootsAOK_V1, BootsAOK_SC] = ...
    areaOverKernel_GrandAv_MOTOR(bootSamps, analysisStartMS, analysisDurMS, dp_cut, nSigma, kernType)
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
% dp_cut, true false, whether to sub-select sessions based on delta dprime
% (effect of LED).
% nSigma = number of sigma on delta d' distribution to drop. 

% kernType
% Can take either 'FA' or 'RT' to generate data for 

analysisStartBin = 501+analysisStartMS; % Lever Release is at bin 601
analysisEndBin = analysisStartBin + analysisDurMS;
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
AOK_V1 = zeros(1,1);
AOK_SC = zeros(1, 1);
AOK_V1_95CI  = zeros(1,2);
AOK_SC_95CI  = zeros(1,2);
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

AOK_V1(1,1) = sum(mean(-1*scKernel(:,analysisStartBin:analysisEndBin)));
AOK_SC(1,1) = sum(mean(-1*v1Kernel(:,analysisStartBin:analysisEndBin)));
%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
v1_nTotal = size(v1Kernel,1);
sc_nTotal = size(scKernel,1);
% Init Archives For BootStrap Samples
BootsAOK_V1 = zeros(bootSamps, 1);
BootsAOK_SC = zeros(bootSamps, 1);

for bootNum = 1:bootSamps
    % Samps For This Round of BootStrapping
    v1Samps = randsample(v1_nTotal,v1_nTotal,true)';
    scSamps = randsample(sc_nTotal,sc_nTotal,true)';
   
    % Take Samples w/ replacement
    Boot_V1 = v1Kernel(v1Samps,:);
    Boot_SC = scKernel(scSamps,:);


    BootsAOK_V1(bootNum,1) = sum(mean(-1*Boot_V1(:,analysisStartBin:analysisEndBin)));
    BootsAOK_SC(bootNum,1) = sum(mean(-1*Boot_SC(:,analysisStartBin:analysisEndBin)));


end
AOK_V1_95CI(1,:) = quantile(BootsAOK_V1,[0.025 0.975]);
AOK_SC_95CI(1,:) = quantile(BootsAOK_SC,[0.025 0.975]);

%% PLOT HISTOGRAMS

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
title('FA: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
ylabel('Probability');
xlabel('Area Over The Kernel');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
ylim([0 0.4]);
xlim([-2.5 2.5]);
legend('V1', 'SC', 'Location','northwest');
hold off;

% Different from 0?
v1p = signrank(BootsAOK_V1);
scp = signrank(BootsAOK_SC);

end
