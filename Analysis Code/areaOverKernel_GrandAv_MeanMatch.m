function [lumBootsAOK_V1, gabBootsAOK_V1] = ...
    areaOverKernel_GrandAv_MeanMatch(bootSamps, analysisStartMS, analysisDurMS)
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


analysisStartBin = 401+analysisStartMS; % Stim on is at bin 401
analysisEndBin = analysisStartBin + analysisDurMS;
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
gabAOK_V1 = zeros(1,1);
lumAOK_V1 = zeros(1, 1);
lumAOK_V1_95CI  = zeros(1,2);
gabAOK_V1_95CI  = zeros(1,2);

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


% Invert 
gabAOK_V1(1,1) = sum(mean(-1*gabKernel(:,analysisStartBin:analysisEndBin)));
lumAOK_V1(1,1) = sum(mean(-1*lumKernel(:,analysisStartBin:analysisEndBin)));


%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
lum_nHit = size(lumProfiles.hitProfiles,1);
gab_nHit = size(gabProfiles.hitProfiles,1);
lum_nTotal = size(lumKernel,1);
gab_nTotal = size(gabKernel,1);

% Init Archives For BootStrap Samples
lumBootsAOK_V1 = zeros(bootSamps, 1);
gabBootsAOK_V1 = zeros(bootSamps, 1);

for bootNum = 1:bootSamps
    % Samps For This Round of BootStrapping
    lumSamps = [randsample(lum_nHit,lum_nHit,true)'...
        randsample([lum_nHit+1:lum_nTotal],lum_nTotal-lum_nHit,true)]';
    gabSamps = [randsample(gab_nHit,gab_nHit,true)'...
        randsample([gab_nHit+1:gab_nTotal],gab_nTotal - gab_nHit,true)]';
    % Take Samples w/ replacement
    lumBoot_V1 = lumKernel(lumSamps,:);
    gabBoot_V1 = gabKernel(gabSamps,:);


    lumBootsAOK_V1(bootNum,1) = sum(mean(-1*lumBoot_V1(:,analysisStartBin:analysisEndBin)));
    gabBootsAOK_V1(bootNum,1) = sum(mean(-1*gabBoot_V1(:,analysisStartBin:analysisEndBin)));


end
lumAOK_V1_95CI(1,:) = quantile(lumBootsAOK_V1,[0.025 0.975]);
gabAOK_V1_95CI(1,:) = quantile(gabBootsAOK_V1,[0.025 0.975]);


%% Histograms

% Colors
V1Color = [0.8500 0.3250 0.0980];
bins = -3:0.2:3;

% Figure
figure('Position', [10 10 1000 500]);
axis square;
hold on;
histogram(gabBootsAOK_V1,bins, 'FaceColor', V1Color, 'FaceAlpha', 0.5, 'Normalization','probability');
histogram(lumBootsAOK_V1,bins, 'FaceColor', V1Color*0.2, 'FaceAlpha', 0.5, 'Normalization','probability');
title('V1: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
ylabel('Probability');
xlabel('Area Over The Kernel');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
ylim([0 0.3]);
xlim([-2.6 2.6]);
plot([0 0], [0 max(ylim)], 'Color', 'k', 'LineStyle', '--', 'LineWidth',1);
plot([median(gabBootsAOK_V1) median(gabBootsAOK_V1)], [0 max(ylim)], 'Color', V1Color, 'LineStyle', '--', 'LineWidth',1);
plot([median(lumBootsAOK_V1) median(lumBootsAOK_V1)], [0 max(ylim)], 'Color', V1Color*0.2, 'LineStyle', '--', 'LineWidth',1);
legend('Gabor', 'Luminance', '','','','Location','northwest');
hold off;

% Different from 0?
signrank(lumBootsAOK_V1)
signrank(gabBootsAOK_V1)

