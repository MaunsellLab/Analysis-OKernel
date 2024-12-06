% Run meanMatchedKernels.m to do the mean matching and save the profiles. 
% This script just loads the profiles and makes the kernel plots with RTs

%% Get Luminance Profiles

% Location of Profiles and Table
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';

% Stuff For Plots
V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];
[plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
yLabel = 'Normalized Power';
% Set Limits
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
limits.yAxis = 0.5;
limits.yAxis = 0.0;
limits.numBoot = 250;

% Load Luminance Tables and Profiles
s = ' ';
cond = 'Lum MM';
lumFolder = [analysisDir,'V1',s,cond,'/'];
load([lumFolder 'masterTable.mat'], 'V1Lum_DDP');
lumT = V1Lum_DDP; clear V1Lum_DDP;

U = selectUsingLimits(lumT, limits);
condition = ['V1' ' ' cond];

% Get all RTs
lumHitRTs = [];
for sessionNum = 1:height(lumT)
    lumHitRTs = [lumHitRTs, cell2mat(lumT.stimCorrectRTs(sessionNum))];
end

% Extract optoProfiles
lumProfiles.hitProfiles = [];
lumProfiles.missProfiles = [];
lumProfiles.earlyProfiles = [];
lumProfiles.RTProfiles = [];
lumProfiles.stimRTProfiles = [];

% Load Pre-created Table Matched for Delta D-Prime
load('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Lum MM/stimProfiles/lumProfiles_DDP.mat');
lumProfiles.hitProfiles = [lumProfiles_DDP.hitProfiles];
lumProfiles.missProfiles = [lumProfiles_DDP.missProfiles];
lumProfiles.earlyProfiles = [lumProfiles_DDP.earlyProfiles];
lumProfiles.RTProfiles = [lumProfiles_DDP.RTProfiles];
lumProfiles.stimRTProfiles = [lumProfiles_DDP.stimRTProfiles];


%% Get Gabor Profiles

% Load Gabor Tables and Profiles
cond = 'Gabor MM';
gabFolder = [analysisDir,'V1',s,cond,'/'];
load([gabFolder 'masterTable.mat'], 'V1Gab_DDP');
gabT = V1Gab_DDP; clear V1Gab_DDP;

U = selectUsingLimits(gabT, limits);
condition = ['V1' ' ' cond];

% Get all RTs
gabHitRTs = [];
for sessionNum = 1:height(gabT)
    gabHitRTs = [gabHitRTs, cell2mat(gabT.stimCorrectRTs(sessionNum))];
end

% Extract optoProfiles
gabProfiles.hitProfiles = [];
gabProfiles.missProfiles = [];
gabProfiles.earlyProfiles = [];
gabProfiles.RTProfiles = [];
gabProfiles.stimRTProfiles = [];

load('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Gabor MM/stimProfiles/gabProfiles_DDP.mat');
gabProfiles.hitProfiles = [gabProfiles_DDP.hitProfiles];
gabProfiles.missProfiles = [gabProfiles_DDP.missProfiles];
gabProfiles.earlyProfiles = [gabProfiles_DDP.earlyProfiles];
gabProfiles.RTProfiles = [gabProfiles_DDP.RTProfiles];
gabProfiles.stimRTProfiles = [gabProfiles_DDP.stimRTProfiles];

%% Stats of Performance and Collect Traces
nHitsLum  = size(lumProfiles.hitProfiles, 1);
nMissLum  = size(lumProfiles.missProfiles, 1);
nEarlyLum = size(lumProfiles.earlyProfiles, 1);

nHitsGab  = size(gabProfiles.hitProfiles, 1);
nMissGab  = size(gabProfiles.missProfiles, 1);
nEarlyGab = size(gabProfiles.earlyProfiles, 1);

% Luminance Kernel
lumKernel = [lumProfiles.hitProfiles; -lumProfiles.missProfiles];
lumKernel_hit = lumProfiles.hitProfiles / 2 + 0.5;
lumKernel_miss = lumProfiles.missProfiles / 2 + 0.5;
lumKernel_RT = lumProfiles.RTProfiles / 2 + 0.5;
lumKernel_early = lumProfiles.earlyProfiles / 2 + 0.5;

% Gabor Kernels
gabKernel = [gabProfiles.hitProfiles; -gabProfiles.missProfiles];
gabKernel_hit = gabProfiles.hitProfiles / 2 + 0.5;
gabKernel_miss = gabProfiles.missProfiles / 2 + 0.5;
gabKernel_RT = gabProfiles.RTProfiles / 2 + 0.5;
gabKernel_early = gabProfiles.earlyProfiles / 2 + 0.5;
%% Prep Plots and BootStrap

% Design Filter For Kernels
sampleFreqHz = 1000;
filterLP = designfilt('lowpassfir', 'PassbandFrequency', 90 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 90 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod','equiripple');

% Bootstrap Gabor
gabBootKernel = bootstrp(250, @mean, gabKernel);
gabPCs = prctile(gabBootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
gabPCMeans = mean(gabPCs, 2);
CIs = zeros(3, size(gabBootKernel, 2));
for c = 1:3
    gabCIs(c, :) = filtfilt(filterLP, gabPCs(c, :) - gabPCMeans(c)) + gabPCMeans(c);
end
gabx = 1:size(gabCIs, 2);

% Bootstrap Luminance
lumBootKernel = bootstrp(250, @mean, lumKernel);
lumPCs = prctile(lumBootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
lumPCMeans = mean(lumPCs, 2);
CIs = zeros(3, size(lumBootKernel, 2));
for c = 1:3
    lumCIs(c, :) = filtfilt(filterLP, lumPCs(c, :) - lumPCMeans(c)) + lumPCMeans(c);
end
lumx = 1:size(lumCIs, 2);

%% Plot

bins = size(gabBootKernel, 2);
figure('Units', 'inches', 'Position', [3, 1, 8, 10]);
% Gabor
subplot(2,1,1);
yyaxis left
plot([0, bins], [limits.yAxis, limits.yAxis], 'Color', 'k', 'LineStyle','--');
hold on;
plot(gabx, gabCIs(2, :), 'Color', V1Color, 'LineStyle', '-', 'LineWidth',2);
gabx2 = [gabx, fliplr(gabx)];
gabfillCI = [gabCIs(1, :), fliplr(gabCIs(3, :))];
fill(gabx2, gabfillCI, V1Color, 'lineStyle', '-', 'edgeColor', V1Color, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
ax = gca;
ax.YColor = [0 0 0];
ax.TickDir = 'out';
ax.XGrid = 'on';
ax.FontSize = 14;
bins = size(gabKernel, 2);
xlim(ax, [0, bins]);
ylabel(yLabel);
set(gca,'XTick', [0, -plotStartMS, bins]);
set(gca, 'XTickLabel', {sprintf('%d', plotStartMS), '0', sprintf('%d', plotEndMS)});
set(gca, 'LineWidth', 1);
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '--';
ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
xlabel('Time Relative to Stimulus');

% Plot RTs on right y-axes
edges = [0:25:800];
yyaxis right
histogram(gabHitRTs+400,edges, 'Normalization', 'probability', 'FaceAlpha', 0.1, 'FaceColor', V1Color, 'edgeAlpha', 0.5); % Have to shift by +400 to align with zero on the x-axis
ylabel('probability of behavioral response');
plotTitle = sprintf('Gab Kernel (n=%d)', nHitsGab+nMissGab);
title(plotTitle);
hold off;

% Lum
subplot(2,1,2);
bins = size(lumBootKernel, 2);
yyaxis left
ax.YColor = [0 0 0];
plot([0, bins], [limits.yAxis, limits.yAxis],  'Color', 'k', 'LineStyle','--');
hold on;
plot(lumx, lumCIs(2, :), 'Color', V1Color, 'LineStyle', '-', 'LineWidth', 2);
ax.YColor = [0 0 0];
lumx2 = [lumx, fliplr(lumx)];
lumfillCI = [lumCIs(1, :), fliplr(lumCIs(3, :))];
fill(lumx2, lumfillCI, V1Color, 'lineStyle', '-', 'edgeColor', V1Color, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
ax = gca;
ax.YColor = 'k';
ax.TickDir = 'out';
ax.FontSize = 14;
bins = size(lumKernel, 2);
xlim(ax, [0, bins]);
ax.XGrid = 'on';
ylabel(yLabel);
title(plotTitle);

set(gca,'XTick', [0, -plotStartMS, bins]);
set(gca, 'XTickLabel', {sprintf('%d', plotStartMS), '0', sprintf('%d', plotEndMS)});
set(gca, 'LineWidth', 1);
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '--';
ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
xlabel('Time Relative to Stimulus');
% Plot RTs on right y-axes
edges = [0:25:800];
yyaxis right
histogram(lumHitRTs+400,edges, 'Normalization', 'probability', 'FaceAlpha', 0.1, 'FaceColor', V1Color, 'edgeAlpha', 0.5); % Have to shift by +400 to align with zero on the x-axis
ylabel('probability of behavioral response');
plotTitle = sprintf('Lum Kernel (n=%d)', nHitsLum+nMissLum);
title(plotTitle);
hold off;

% Put plots on same axes
%lims = [min(min(llim1, llim2)) max(max(llim1, llim2))];
lims = [-0.04 0.02];
subplot(2,1,1)
yyaxis left
ylim(lims);
yyaxis right
ax.YColor = 'k';
ylim([0 1]);
set(gca,'YTick', [0, 0.25]);

subplot(2,1,2);
yyaxis left
ax.YColor = 'k';
ylim(lims);
yyaxis right
ax.YColor = 'k';
ylim([0 1]);
set(gca,'YTick', [0, 0.25]);
