% AOK by d'
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
V1LumDelta = [U.stimDPrime] - [U.noStimDPrime];

% Idx for which sessions were excluded.
mu =  -0.011;
sigma =  0.2326;
deltaD = [U.stimDPrime]-[U.noStimDPrime];
cutoff = mu + 1*sigma;
keepIdxV1Lum = V1LumDelta <= cutoff;

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
V1GabDelta = [U.stimDPrime] - [U.noStimDPrime];

% Idx for exluded sessions
mu =  -0.0034;
sigma = 0.277;
cutoff = mu + 1*sigma;
keepIdxV1Gab = V1GabDelta <= cutoff;

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
SCLumDelta = [U.stimDPrime] - [U.noStimDPrime];

% Idx For Excluded Sessions
mu = -0.03914;
sigma = 0.2715;
cutoff = mu + 1*sigma;
keepIdxSCLum = SCLumDelta <= cutoff;

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
SCGabDelta = [U.stimDPrime] - [U.noStimDPrime];

% Idx of excluded sessions
mu = 0.0056;
sigma = 0.287;
bins = -1:0.1:1;
cutoff = mu + 1*sigma;
keepIdxSCGab = SCGabDelta <= cutoff;

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
[rho,p] = corr(V1GabAOKs, V1GabDelta);
[rho,p] = corr(SCGabAOKs, SCGabDelta);
[rho,p] = corr(V1LumAOKs, V1LumDelta);
[rho,p] = corr(SCLumAOKs, SCLumDelta);

% Fit Linear Function
f1 = fitlm(V1GabAOKs, V1GabDelta);
f2 = fitlm(SCGabAOKs, SCGabDelta);
f3 = fitlm(V1LumAOKs, V1LumDelta);
f4 = fitlm(SCLumAOKs, SCLumDelta);
%% Scatters
xs = -15:0.01:15;

figure('Position',[10 10 800 1000]);
subplot(2,2,1);
hold on;
scatter(V1GabAOKs(keepIdxV1Gab), V1GabDelta(keepIdxV1Gab), 30, 'black', 'filled');
scatter(V1GabAOKs(~keepIdxV1Gab), V1GabDelta(~keepIdxV1Gab), 30, 'black');
plot(xs, f1.Coefficients.Estimate(1) + xs.*f1.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Gabor: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-15 15]);
ylim([-0.8 0.8]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [-0.8 -0.4 0 0.4 0.8]);
box off;
hold off;

subplot(2,2,2);
hold on;
scatter(V1LumAOKs(keepIdxV1Lum), V1LumDelta(keepIdxV1Lum), 30, 'black', 'filled');
scatter(V1LumAOKs(~keepIdxV1Lum), V1LumDelta(~keepIdxV1Lum), 30, 'black');
plot(xs, f3.Coefficients.Estimate(1) + xs.*f3.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Luminance: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-15 15]);
ylim([-0.8 0.8]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [-0.8 -0.4 0 0.4 0.8]);
box off;
hold off;


subplot(2,2,3);
hold on;
scatter(SCGabAOKs(keepIdxSCGab), SCGabDelta(keepIdxSCGab), 30, 'black', 'filled');
scatter(SCGabAOKs(~keepIdxSCGab), SCGabDelta(~keepIdxSCGab), 30, 'black');
plot(xs, f2.Coefficients.Estimate(1) + xs.*f2.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Gabor: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-15 15]);
ylim([-0.8 0.8]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [-0.8 -0.4 0 0.4 0.8]);
box off;
hold off;

subplot(2,2,4);
hold on;
scatter(SCLumAOKs(keepIdxSCLum), SCLumDelta(keepIdxSCLum), 30, 'black', 'filled');
scatter(SCLumAOKs(~keepIdxSCLum), SCLumDelta(~keepIdxSCLum), 30, 'black');
plot(xs, f4.Coefficients.Estimate(1) + xs.*f4.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Luminance: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-15 15]);
ylim([-0.8 0.8]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-15 0 15]);
set(gca, 'YTick', [-0.8 -0.4 0 0.4 0.8]);
box off;
hold off;
 


