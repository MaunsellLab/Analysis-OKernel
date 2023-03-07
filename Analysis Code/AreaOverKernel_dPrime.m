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
%% Scatters
figure('Position',[10 10 800 1000]);
subplot(2,2,1);
hold on;
scatter(V1LumAOKs, V1LumDelta, 30, 'black', 'filled');
title('V1 Luminance: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-15 15]);
ylim([-0.8 0.8]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
box off;
hold off;

subplot(2,2,2);
hold on;
scatter(V1GabAOKs, V1GabDelta, 30, 'black', 'filled');
title('V1 Gabor: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-15 15]);
ylim([-0.8 0.8]);
set(gca, 'FontSize', 16);
set(gca, 'TickDir', 'out');
box off;
hold off; 

subplot(2,2,3);
hold on;
scatter(SCLumAOKs, SCLumDelta, 30, 'black', 'filled');
title('SC Luminance: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-15 15]);
ylim([-0.8 0.8]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
box off;
hold off;

subplot(2,2,4);
hold on;
scatter(SCGabAOKs, SCGabDelta, 30, 'black', 'filled');
title('SC Gabor: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-15 15]);
ylim([-0.8 0.8]);
set(gca, 'FontSize', 16);
set(gca, 'TickDir', 'out');
box off;
hold off; 
