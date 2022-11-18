% plot Kernel and RTs

% Directory With all the final Tables and Profiles
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';

% Set Limits
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
limits.yAxis = 0.5;
limits.yAxis = 0.0;
limits.numBoot = 250;


% Load Tables for each stimulus type
lumFolder = strcat(analysisDir, brainArea,' Lum/');
load([lumFolder 'masterTable.mat'], 'T');
lumT = T;

gabFolder = strcat(analysisDir, brainArea,' Gabor/');
load([gabFolder 'masterTable.mat'], 'T');
gabT = T;
clear T;



% total kernel trials weighted across all trials. We need to multiple the weighted sum by 2 because it is effectively
% a mean of the hit and miss kernels, not a difference. By taking the mean, we lose the doubling that we should get
% from the opposing effects.  This has been validated in simulations.
h = figure;
set(h, 'Units', 'inches', 'Position', [25, 1.25, 8.5, 11]);
clf;
ylabel = 'Normalized Power';
plotTitle = sprintf('Total Kernel (n=%d)', numHits + numMisses);
hitMissBoot = [stimProfiles.hitProfiles; -stimProfiles.missProfiles];
doOneBootPlot(hitMissBoot, limits, 'stim', plotStartMS, plotEndMS, plotTitle, '');