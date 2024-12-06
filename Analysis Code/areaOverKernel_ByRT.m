function [lumBootsAOK_V1, lumBootsAOK_SC, gabBootsAOK_V1, gabBootsAOK_SC] = ...
    areaOverKernel_ByRT(bootSamps, analysisStartMS, analysisDurMS, dp_cut, nSigma)
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
lumFolder = strcat(analysisDir, 'V1',' Lum/');
gabFolder = strcat(analysisDir, 'V1',' Gabor/');

% Load Tables for each stimulus type
load([lumFolder 'masterTable.mat'], 'T');
lumT = T;
load([gabFolder 'masterTable.mat'], 'T');
gabT = T;
clear T;

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

% Get RTs
lumRTs_V1 =  cell2mat([U.stimCorrectRTs]');

%% Repeat For Gabor
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

gabRTs_V1 =  cell2mat([U.stimCorrectRTs]');


% Collate All Traces From This Mouse
gabKernel = [gabProfiles.hitProfiles];
lumKernel = [lumProfiles.hitProfiles];

medGabRT_V1 = median(gabRTs_V1)
medLumRT_V1 = median(lumRTs_V1)

% Compute Mean AOK during the analysis Window
gabAOK_V1(1,1) = sum(mean(-1*gabKernel(:,analysisStartBin:analysisEndBin)));
lumAOK_V1(1,1) = sum(mean(-1*lumKernel(:,analysisStartBin:analysisEndBin)));


%% Labels for Median Splits

medLogLumV1 = (lumRTs_V1 <= medLumRT_V1);
medLogGabV1 = (gabRTs_V1 <= medGabRT_V1);
belowMedGab = gabKernel(medLogGabV1,:); 
aboveMedGab = gabKernel(~medLogGabV1,:); 

%% Plot


V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];
bins = -3:0.2:3;

RTbins = 0:25:600;

% Figure
figure('Position', [10 10 1000 500]);
subplot(1,2,2); 
axis square;
hold on;
histogram(gabRTs_V1(medLogGabV1),RTbins, 'FaceColor', V1Color, 'FaceAlpha', 0.5, 'Normalization','probability');
histogram(gabRTs_V1(~medLogGabV1),RTbins, 'FaceColor', SCColor, 'FaceAlpha', 0.5, 'Normalization','probability');
ylabel('Probability');
xlabel('Reaction Time (ms)');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
box off;
title('Reaction Times');
legend('< Median', '> Median', 'Location','northeast');


% Set Limits
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
limits.yAxis = 0.5;
limits.yAxis = 0.0;
limits.numBoot = 250;

% Stuff For Plots
[plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
yLabel = 'Normalized Power';






ax = gca;
ax.YColor = [0 0 0];
ax.TickDir = 'out';
ax.XGrid = 'on';
ax.FontSize = 14;
xlim(ax, [0, bins]);
ylabel(yLabel);

set(gca,'XTick', [0, -plotStartMS, bins]);
set(gca, 'XTickLabel', {sprintf('%d', plotStartMS), '0', sprintf('%d', plotEndMS)});
set(gca, 'LineWidth', 1);
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '--';
ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
xlabel('Time Relative to Stimulus');









box off;
set(gca, 'TickDir', 'out');
ylabel('Probability');
xlabel('Area Over The Kernel');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
ylim([0 0.4]);
plot([0 0], [0 max(ylim)], 'Color', 'k', 'LineStyle', '--', 'LineWidth',1);
plot([median(gabBootsAOK_V1) median(gabBootsAOK_V1)], [0 max(ylim)], 'Color', V1Color, 'LineStyle', '--', 'LineWidth',1);
plot([median(gabBootsAOK_SC) median(gabBootsAOK_SC)], [0 max(ylim)], 'Color', SCColor, 'LineStyle', '--', 'LineWidth',1);
xlim([-2.5 2.5]);
legend('V1', 'SC', 'Location','northwest');
hold off;






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
    lumSamps = randsample(lum_nHit,lum_nHit,true);
    gabSamps = randsample(gab_nHit,gab_nHit,true);
    % Take Samples w/ replacement
    lumBoot_V1 = lumKernel(lumSamps,:);
    gabBoot_V1 = gabKernel(gabSamps,:);


    lumBootsAOK_V1(bootNum,1) = sum(mean(-1*lumBoot_V1(:,analysisStartBin:analysisEndBin)));
    gabBootsAOK_V1(bootNum,1) = sum(mean(-1*gabBoot_V1(:,analysisStartBin:analysisEndBin)));


end

lumAOK_V1_95CI(1,:) = quantile(lumBootsAOK_V1,[0.025 0.975]);
gabAOK_V1_95CI(1,:) = quantile(gabBootsAOK_V1,[0.025 0.975]);

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

% Init
gabAOK_SC = zeros(1,1);
lumAOK_SC = zeros(1, 1);
lumAOK_SC_95CI  = zeros(1,2);
gabAOK_SC_95CI  = zeros(1,2);

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

lumRTs_SC =  cell2mat([U.stimCorrectRTs]');


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


gabRTs_SC =  cell2mat([U.stimCorrectRTs]');


% Collate All Traces From This Mouse
gabKernel = [gabProfiles.hitProfiles];
lumKernel = [lumProfiles.hitProfiles];

medGabRT_SC = median(gabRTs_SC)
medLumRT_SC = median(lumRTs_SC)

gabAOK_SC(1,1) = sum(mean(-1*gabKernel(:,analysisStartBin:analysisEndBin)));
lumAOK_SC(1,1) = sum(mean(-1*lumKernel(:,analysisStartBin:analysisEndBin)));



%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
lum_nHit = size(lumProfiles.hitProfiles,1);
gab_nHit = size(gabProfiles.hitProfiles,1);
lum_nTotal = size(lumKernel,1);
gab_nTotal = size(gabKernel,1);

% Init Archives For BootStrap Samples
lumBootsAOK_SC = zeros(bootSamps, 1);
gabBootsAOK_SC = zeros(bootSamps, 1);

for bootNum = 1:bootSamps
    % Samps For This Round of BootStrapping
    lumSamps = randsample(lum_nHit,lum_nHit,true);
    gabSamps = randsample(gab_nHit,gab_nHit,true);
    % Take Samples w/ replacement
    lumBoot_SC = lumKernel(lumSamps,:);
    gabBoot_SC = gabKernel(gabSamps,:);
    

    lumBootsAOK_SC(bootNum,1) = sum(mean(-1*lumBoot_SC(:,analysisStartBin:analysisEndBin)));
    gabBootsAOK_SC(bootNum,1) = sum(mean(-1*gabBoot_SC(:,analysisStartBin:analysisEndBin)));
end

lumAOK_SC_95CI(1,:) = quantile(lumBootsAOK_SC,[0.025 0.975]);
gabAOK_SC_95CI(1,:) = quantile(gabBootsAOK_SC,[0.025 0.975]);


%% Histograms
% Colors
V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];
bins = -3:0.2:3;
% Figure
figure('Position', [10 10 1000 500]);
subplot(1,2,1); %Gabors First
axis square;
hold on;
histogram(gabBootsAOK_V1,bins, 'FaceColor', V1Color, 'FaceAlpha', 0.5, 'Normalization','probability');
histogram(gabBootsAOK_SC,bins, 'FaceColor', SCColor, 'FaceAlpha', 0.5, 'Normalization','probability');
title('Gabor: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
ylabel('Probability');
xlabel('Area Over The Kernel');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
ylim([0 0.4]);
plot([0 0], [0 max(ylim)], 'Color', 'k', 'LineStyle', '--', 'LineWidth',1);
plot([median(gabBootsAOK_V1) median(gabBootsAOK_V1)], [0 max(ylim)], 'Color', V1Color, 'LineStyle', '--', 'LineWidth',1);
plot([median(gabBootsAOK_SC) median(gabBootsAOK_SC)], [0 max(ylim)], 'Color', SCColor, 'LineStyle', '--', 'LineWidth',1);
xlim([-2.5 2.5]);
legend('V1', 'SC', 'Location','northwest');
hold off;

subplot(1,2,2); % Luminance
axis square;
hold on;
histogram(lumBootsAOK_V1,bins, 'FaceColor', V1Color, 'FaceAlpha', 0.5, 'Normalization','probability');
histogram(lumBootsAOK_SC,bins, 'FaceColor', SCColor, 'FaceAlpha', 0.5,'Normalization','probability');
title('Luminance: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
ylabel('Probability');
xlabel('Area Over The Kernel');
ylim([0 0.4]);
plot([0 0], [0 max(ylim)], 'Color', 'k', 'LineStyle', '--', 'LineWidth',1);
plot([median(lumBootsAOK_V1) median(lumBootsAOK_V1)], [0 max(ylim)], 'Color', V1Color, 'LineStyle', '--', 'LineWidth',1);
plot([median(lumBootsAOK_SC) median(lumBootsAOK_SC)], [0 max(ylim)], 'Color', SCColor, 'LineStyle', '--', 'LineWidth',1);
xlim([-2.5 2.5]);
legend('V1', 'SC', 'Location','northwest');
hold off;


% Different from 0?
v1p_lum = (length(lumBootsAOK_V1) - sum(lumBootsAOK_V1>0))/length(lumBootsAOK_V1)
scp_lum = (length(lumBootsAOK_SC) - sum(lumBootsAOK_SC>0))/length(lumBootsAOK_SC)

v1p_gab = (length(gabBootsAOK_V1) - sum(gabBootsAOK_V1>0))/length(gabBootsAOK_V1)
scp_gab = (length(gabBootsAOK_SC) - sum(gabBootsAOK_SC>0))/length(gabBootsAOK_SC)




