% AOK by average C on tested stimulus
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
animals = {'1960', '2015', '2083', '2126', '2220', '2221'};
condition = ['V1' ' ' 'Lum'];

% Init AOK and d' Storage
V1LumAOKs = zeros(length(animals),1);
V1LumCs = zeros(length(animals),1);

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(lumT, limits);
    V1LumCs(a,1) = nanmean(([U.stimC]+[U.noStimC])/2);
    lumKernel = []; % Init
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        lumKernel = [lumKernel; [stimProfiles.hitProfiles; -stimProfiles.missProfiles]];
    end
    % Compute AOK for this sessions
    % Computes average Kernel for corresponding session, inverts (-1) and
    % sums across the analysis window
    V1LumAOKs(a,1) = -1*sum(mean(lumKernel(:,analysisStartBin:analysisEndBin))); % invert
end
%% Get StimProfiles and Compute AOK For V1 Gabor
animals = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
condition = ['V1' ' ' 'Gabor'];

% Init AOK Storage
V1GabAOKs = zeros(length(animals),1);
V1GabCs = zeros(length(animals),1);

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(gabT, limits);
    V1GabCs(a,1) = nanmean(([U.stimC]+[U.noStimC])/2);
    gabKernel = []; % Init
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        gabKernel = [gabKernel; [stimProfiles.hitProfiles; -stimProfiles.missProfiles]];
    end
    V1GabAOKs(a,1) = -1*sum(mean(gabKernel(:,analysisStartBin:analysisEndBin))); % invert
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
animals = {'1458', '1548', '1674', '1675', '1902', '1905', '2057', '2058', '2063', '2169', '2236'};
condition = ['SC' ' ' 'Lum'];

% Init AOK and d' Storage
SCLumAOKs = zeros(length(animals),1);
SCLumCs = zeros(length(animals),1);

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(lumT, limits);
    SCLumCs(a,1) = nanmean(([U.stimC]+[U.noStimC])/2);
    lumKernel = []; % Init
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        lumKernel = [lumKernel; [stimProfiles.hitProfiles; -stimProfiles.missProfiles]];
    end
    SCLumAOKs(a,1) = -1*sum(mean(lumKernel(:,analysisStartBin:analysisEndBin))); % invert
end
%% SC Gabor
animals = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
condition = ['SC' ' ' 'Gabor'];
% Init AOK and d' Storage
SCGabAOKs = zeros(length(animals),1);
SCGabCs = zeros(length(animals),1);

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(gabT, limits);
    SCGabCs(a,1) = nanmean(([U.stimC]+[U.noStimC])/2);
    gabKernel = []; % Init
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        gabKernel = [gabKernel; [stimProfiles.hitProfiles; -stimProfiles.missProfiles]];
    end
    SCGabAOKs(a,1) = -1*sum(mean(gabKernel(:,analysisStartBin:analysisEndBin))); % invert
end
%% Scatters and Correlation
[rho1,p1] = corr(V1GabAOKs, V1GabCs);
[rho2,p2] = corr(SCGabAOKs, SCGabCs);
[rho3,p3] = corr(V1LumAOKs, V1LumCs);
[rho4,p4] = corr(SCLumAOKs, SCLumCs);

% Fit Linear Function
f1 = fitlm(V1GabAOKs, V1GabCs);
f2 = fitlm(SCGabAOKs, SCGabCs);
f3 = fitlm(V1LumAOKs, V1LumCs);
f4 = fitlm(SCLumAOKs, SCLumCs);

% plot
xs = -5:0.01:5;

figure('Position',[10 10 800 1000]);
subplot(2,2,1);
hold on;
scatter(V1GabAOKs, V1GabCs, 50, 'black', 'filled');
plot(xs, f1.Coefficients.Estimate(1) + xs.*f1.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Gabor: AOK by C');
xlabel('Area Over the Kernel');
ylabel('Criterion');
axis square;
xlim([-5 5]);
ylim([0 1]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-5 0 5]);
set(gca, 'YTick', [0 1]);
box off;
hold off;

subplot(2,2,2);
hold on;
scatter(V1LumAOKs, V1LumCs, 50, 'black', 'filled');
plot(xs, f3.Coefficients.Estimate(1) + xs.*f3.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Luminance: AOK by C');
xlabel('Area Over the Kernel');
ylabel('Criterion');
axis square;
xlim([-5 5]);
ylim([0 1]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-5 0 5]);
set(gca, 'YTick', [0 1]);
box off;
hold off;


subplot(2,2,3);
hold on;
scatter(SCGabAOKs, SCGabCs, 50, 'black', 'filled');
plot(xs, f2.Coefficients.Estimate(1) + xs.*f2.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Gabor: AOK by C');
xlabel('Area Over the Kernel');
ylabel('Criterion');
axis square;
xlim([-5 5]);
ylim([0 1]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-5 0 5]);
set(gca, 'YTick', [0 1]);
box off;
hold off;

subplot(2,2,4);
hold on;
scatter(SCLumAOKs, SCLumCs, 50, 'black', 'filled');
plot(xs, f4.Coefficients.Estimate(1) + xs.*f4.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Luminance: AOK by C');
xlabel('Area Over the Kernel');
ylabel('Criterion');
axis square;
xlim([-5 5]);
ylim([0 1]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-5 0 5]);
set(gca, 'YTick', [0 1]);
box off;
hold off;
 
