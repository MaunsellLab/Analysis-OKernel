% AOK by criterion
%% Set Up
analysisStartMS = 0;
analysisDurMS = 100;
analysisStartBin = 401+analysisStartMS; % Stim on is at bin 401
analysisEndBin = analysisStartBin + analysisDurMS;
%% Location of the Final Stim Profiles and Tables For V1
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
%% Get StimProfiles and Compute AOK For V1 Luminance
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
limits.animal = {'1960', '2015', '2083', '2126', '2220', '2221'};
U = selectUsingLimits(lumT, limits);
condition = ['V1' ' ' 'Lum'];
V1LumD = ([U.stimDPrime] + [U.noStimDPrime])/2;
V1LumD = [U.noStimDPrime];
V1LumC = ([U.stimC] + [U.noStimC])/2;
V1LumC = [U.noStimC];

% Init AOK Storage
V1LumAOKs = zeros(height(U),1);
% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    lumKernel = []; % Init
    lumKernel = [stimProfiles.hitProfiles; -stimProfiles.missProfiles];
    % Compute AOK for this sessions
    % Computes average Kernel for corresponding session, inverts (-1) and
    % sums across the analysis window
    V1LumAOKs(i,1) = -1*sum(mean(lumKernel(:,analysisStartBin:analysisEndBin))); % invert
end
%% Get StimProfiles and Compute AOK For V1 Gabor
limits.animal = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
U = selectUsingLimits(gabT, limits);
condition = ['V1' ' ' 'Gabor'];
V1GabD = ([U.stimDPrime] + [U.noStimDPrime])/2;
V1GabD = [U.noStimDPrime];
V1GabC = ([U.stimC] + [U.noStimC])/2;
V1GabC = [U.noStimC];

% Init AOK Storage
V1GabAOKs = zeros(height(U),1);
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    gabKernel = []; % Init
    gabKernel = [stimProfiles.hitProfiles; -stimProfiles.missProfiles];
    % Compute AOK for this sessions
    % Computes average Kernel for corresponding session, inverts (-1) and
    % sums across the analysis window
    V1GabAOKs(i,1) = -1*sum(mean(gabKernel(:,analysisStartBin:analysisEndBin)));
end
%% Get SC Data
lumFolder = strcat(analysisDir, 'SC',' Lum/');
gabFolder = strcat(analysisDir, 'SC',' Gabor/');
clear lumT gabT;
% Load Tables for each stimulus type
load([lumFolder 'masterTable.mat'], 'T');
lumT = T;
load([gabFolder 'masterTable.mat'], 'T');
gabT = T;
clear T;
%% SC Luminance
limits.animal = {'1458', '1548', '1674', '1675', '1902', '1905', '2057', '2058', '2063', '2169', '2236'};
U = selectUsingLimits(lumT, limits);
condition = ['SC' ' ' 'Lum'];
SCLumD = ([U.stimDPrime] + [U.noStimDPrime])/2;
SCLumD = [U.noStimDPrime];
SCLumC = ([U.stimC] + [U.noStimC])/2;
SCLumC = [U.noStimC];
% Init
SCLumAOKs = zeros(height(U),1);
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    lumKernel = []; % Init
    lumKernel = [stimProfiles.hitProfiles; -stimProfiles.missProfiles];
    % Compute AOK for this sessions
    % Computes average Kernel for corresponding session, inverts (-1) and
    % sums across the analysis window
    SCLumAOKs(i,1) = -1*sum(mean(lumKernel(:,analysisStartBin:analysisEndBin))); % invert
end
%% SC Gabor
limits.animal = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
U = selectUsingLimits(gabT, limits);
condition = ['SC' ' ' 'Gabor'];
SCGabD = ([U.stimDPrime] + [U.noStimDPrime])/2;
SCGabD = [U.noStimDPrime];
SCGabC = ([U.stimC] + [U.noStimC])/2;
SCGabC = [U.noStimC];
% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    gabKernel = []; % Init
    gabKernel = [stimProfiles.hitProfiles; -stimProfiles.missProfiles];
    % Compute AOK for this sessions
    % Computes average Kernel for corresponding session, inverts (-1) and
    % sums across the analysis window
    SCGabAOKs(i,1) = -1*sum(mean(gabKernel(:,analysisStartBin:analysisEndBin)));
end
%% Correlation

% drop = find(isnan(V1GabD));
% V1GabAOKs(drop) = [];
% V1GabD(drop) = [];
% drop = find(isnan(SCGabD));
% SCGabAOKs(drop) = [];
% SCGabD(drop) = [];

[rhoD1,pD1] = corr(V1GabAOKs, V1GabD);
[rhoD2,pD2] = corr(SCGabAOKs, SCGabD);
[rhoD3,pD3] = corr(V1LumAOKs, V1LumD);
[rhoD4,pD4] = corr(SCLumAOKs, SCLumD);

% Fit Linear Function
f1 = fitlm(V1GabAOKs, V1GabD);
f2 = fitlm(SCGabAOKs, SCGabD);
f3 = fitlm(V1LumAOKs, V1LumD);
f4 = fitlm(SCLumAOKs, SCLumD);

% AOK By D
xs = -15:0.01:15;

figure('Position',[10 10 800 1000]);
subplot(2,2,1);
hold on;
scatter(V1GabAOKs, V1GabD, 30, 'black', 'filled');
plot(xs, f1.Coefficients.Estimate(1) + xs.*f1.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Gabor: AOK by d''');
xlabel('Area Over the Kernel');
ylabel('d''');
axis square;
xlim([-15 15]);
ylim([0.5 4.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [1 2 3 4]);
box off;
hold off;

subplot(2,2,2);
hold on;
scatter(V1LumAOKs, V1LumD, 30, 'black', 'filled');
plot(xs, f3.Coefficients.Estimate(1) + xs.*f3.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Luminance: AOK by d''');
xlabel('Area Over the Kernel');
ylabel('d''');
axis square;
xlim([-15 15]);
ylim([0.5 4.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [1 2 3 4]);
box off;
hold off;


subplot(2,2,3);
hold on;
scatter(SCGabAOKs, SCGabD, 30, 'black', 'filled');
plot(xs, f2.Coefficients.Estimate(1) + xs.*f2.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Gabor: AOK by d''');
xlabel('Area Over the Kernel');
ylabel('d''');
axis square;
xlim([-15 15]);
ylim([0.5 4.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [1 2 3 4]);
box off;
hold off;

subplot(2,2,4);
hold on;
scatter(SCLumAOKs, SCLumD, 30, 'black', 'filled');
plot(xs, f4.Coefficients.Estimate(1) + xs.*f4.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Luminance: AOK by d''');
xlabel('Area Over the Kernel');
ylabel('d''');
axis square;
xlim([-15 15]);
ylim([0.5 4.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [1 2 3 4]);
box off;
hold off;
 
%% AOK By C

% eliminate NaNs
drop = find(isnan(V1GabC));
V1GabAOKs(drop) = [];
V1GabC(drop) = [];
drop = find(isnan(SCGabC));
SCGabAOKs(drop) = [];
SCGabC(drop) = [];


[rhoC1,pC1] = corr(V1GabAOKs, V1GabC);
[rhoC2,pC2] = corr(SCGabAOKs, SCGabC);
[rhoC3,pC3] = corr(V1LumAOKs, V1LumC);
[rhoC4,pC4] = corr(SCLumAOKs, SCLumC);

% Fit Linear Function
f1 = fitlm(V1GabAOKs, V1GabC);
f2 = fitlm(SCGabAOKs, SCGabC);
f3 = fitlm(V1LumAOKs, V1LumC);
f4 = fitlm(SCLumAOKs, SCLumC);


figure('Position',[10 10 800 1000]);
subplot(2,2,1);
hold on;
scatter(V1GabAOKs, V1GabC, 30, 'black', 'filled');
plot(xs, f1.Coefficients.Estimate(1) + xs.*f1.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Gabor: AOK by criterion');
xlabel('Area Over the Kernel');
ylabel('criterion');
axis square;
xlim([-15 15]);
ylim([-0.5 1.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [-0.5 0 0.5 1.0 1.5]);
box off;
hold off;

subplot(2,2,2);
hold on;
scatter(V1LumAOKs, V1LumC, 30, 'black', 'filled');
plot(xs, f3.Coefficients.Estimate(1) + xs.*f3.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Luminance: AOK by criterion');
xlabel('Area Over the Kernel');
ylabel('criterion');
axis square;
xlim([-15 15]);
ylim([-0.5 1.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [-0.5 0 0.5 1.0 1.5]);
box off;
hold off;


subplot(2,2,3);
hold on;
scatter(SCGabAOKs, SCGabC, 30, 'black', 'filled');
plot(xs, f2.Coefficients.Estimate(1) + xs.*f2.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Gabor: AOK by criterion');
xlabel('Area Over the Kernel');
ylabel('criterion');
axis square;
xlim([-15 15]);
ylim([-0.5 1.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [-0.5 0 0.5 1.0 1.5]);
box off;
hold off;

subplot(2,2,4);
hold on;
scatter(SCLumAOKs, SCLumC, 30, 'black', 'filled');
plot(xs, f4.Coefficients.Estimate(1) + xs.*f4.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Luminance: AOK by criterion');
xlabel('Area Over the Kernel');
ylabel('criterion');
axis square;
xlim([-15 15]);
ylim([-0.5 1.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [-0.5 0 0.5 1.0 1.5]);
box off;
hold off;

