function areaOverKernel(brainArea, bootSamps, analysisStartMS, analysisDurMS)
% Compute Area Over Kernel for Individual Mice
% Bootstrap 95%CI per mouse
% Plot Scatter of results

% Inputs:
% brainArea = 'V1' or 'SC'
% bootSamps = How many bootstraps to create CIs
% analysisStartMS = relative to stimOnset (time = 0)
% analysisDurMS = How many bins (ms) to compute the deflection in the
% kernel after the start of analysis

analysisStartBin = 401+analysisStartMS; % Stim on is at bin 401
analysisEndBin = analysisStartBin + analysisDurMS;
%% Load Table and Directory for Stim Profiles
% Location of the Final Stim Profiles and Tables
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';
lumFolder = strcat(analysisDir, brainArea,' Lum/');
gabFolder = strcat(analysisDir, brainArea,' Gabor/');

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

% Only Consider Animals For Which Data was obtained for both stimuli
if strcmp(brainArea,'SC')
    animals = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
elseif strcmp(brainArea,'V1')
    animals = {'1960', '2015', '2083', '2126', '2220', '2221'};
end

% Init
gabAOK = zeros(1,length(animals));
lumAOK = zeros(1, length(animals));
lumAOK_95CI  = zeros(length(animals),2);
gabAOK_95CI  = zeros(length(animals),2);

%% Loop Through this for Both Gabor and Luminance
% First For Luminance
for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};

    U = selectUsingLimits(lumT, limits);
    condition = [brainArea ' ' 'Lum'];

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

    % Repeat For Gabor
    U = selectUsingLimits(gabT, limits);
    condition = [brainArea ' ' 'Gabor'];
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

    % Collate All Traces From This Mouse
    gabKernel = [gabProfiles.hitProfiles; -gabProfiles.missProfiles];
    lumKernel = [lumProfiles.hitProfiles; -lumProfiles.missProfiles];

    % Compute Mean AOK during the analysis Window
    gabAOK(1,a) = sum(mean(gabKernel(:,analysisStartBin:analysisEndBin)));
    lumAOK(1,a) = sum(mean(lumKernel(:,analysisStartBin:analysisEndBin)));

    %% Bootstrap CIs from Stim Profiles (bootNum = 10000)
    % How Many of Each Outcome Type to control bootstrapping to match
    % experimental proportions
    lum_nHit = size(lumProfiles.hitProfiles,1);
    gab_nHit = size(gabProfiles.hitProfiles,1);
    lum_nTotal = size(lumKernel,1);
    gab_nTotal = size(gabKernel,1);

    % Init Archives For BootStrap Samples
    lumBootsAOK = zeros(bootSamps, 1);
    gabBootsAOK = zeros(bootSamps, 1);


    for bootNum = 1:bootSamps
        % Samps For This Round of BootStrapping
        lumSamps = [randsample(lum_nHit,lum_nHit,true)'...
            randsample([lum_nHit+1:lum_nTotal],lum_nTotal-lum_nHit,true)]';
        gabSamps = [randsample(gab_nHit,gab_nHit,true)'...
            randsample([gab_nHit+1:gab_nTotal],gab_nTotal - gab_nHit,true)]';
        % Take Samples w/ replacement
        lumBoot = lumKernel(lumSamps,:);
        gabBoot = gabKernel(gabSamps,:);
        % Compute 95%CI during the analysis Window From BootStraps
        lumBootsAOK(bootNum,1) = sum(mean(lumBoot(:,analysisStartBin:analysisEndBin)));
        gabBootsAOK(bootNum,1) = sum(mean(gabBoot(:,analysisStartBin:analysisEndBin)));
    end
    lumAOK_95CI(a,:) = quantile(lumBootsAOK,[0.025 0.975]);
    gabAOK_95CI(a,:) = quantile(gabBootsAOK,[0.025 0.975]);
end

%% Make Scatter Plot of Results
figure;
hold on;
axis square;
scatter(lumAOK, gabAOK, 60, 'black','filled');
% Plot 95% CIs
for ani = 1:length(animals)
%     plot([lumAOK(ani)-lumAOK_95CI(ani,1)  lumAOK(ani)+lumAOK_95CI(ani,2)], [gabAOK(ani) gabAOK(ani)], 'LineStyle','-', 'LineWidth',1, 'Color', 'k');
%     plot([lumAOK(ani) lumAOK(ani)], [gabAOK(ani)-gabAOK_95CI(ani,1) gabAOK(ani)+gabAOK_95CI(ani,2)], 'LineStyle','-', 'LineWidth',1, 'Color', 'k');

    plot([lumAOK_95CI(ani,1)  lumAOK_95CI(ani,2)], [gabAOK(ani) gabAOK(ani)], 'LineStyle','-', 'LineWidth',1, 'Color', 'k');
    plot([lumAOK(ani) lumAOK(ani)], [gabAOK_95CI(ani,1) gabAOK_95CI(ani,2)], 'LineStyle','-', 'LineWidth',1, 'Color', 'k');
end
% Lines That Divide Quadrants into Effect on Performance
plot([-4 4], [0 0], 'LineStyle','--', 'LineWidth',0.5, 'Color', 'k');
plot([0 0], [-4 4], 'LineStyle','--', 'LineWidth',0.5, 'Color', 'k');

xlabel('Luminance AOK');
ylabel('Gabor AOK');
title(strcat(brainArea,' Area Over Kernel'));
box off;
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
xlim([-4 4]);
ylim([-4 4]);
hold off;
end




