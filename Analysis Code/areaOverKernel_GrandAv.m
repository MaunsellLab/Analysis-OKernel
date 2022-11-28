function areaOverKernel_GrandAv(invert, bootSamps, analysisStartMS, analysisDurMS)
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
condition = ['V1' ' ' 'Gabor'];

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
if invert == 0
    gabAOK_V1(1,1) = sum(mean(gabKernel(:,analysisStartBin:analysisEndBin)));
    lumAOK_V1(1,1) = sum(mean(lumKernel(:,analysisStartBin:analysisEndBin)));
else
    gabAOK_V1(1,1) = sum(mean(-1*gabKernel(:,analysisStartBin:analysisEndBin)));
    lumAOK_V1(1,1) = sum(mean(-1*lumKernel(:,analysisStartBin:analysisEndBin)));
end

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

    if invert == 0
        % Compute 95%CI during the analysis Window From BootStraps
        lumBootsAOK_V1(bootNum,1) = sum(mean(lumBoot_V1(:,analysisStartBin:analysisEndBin)));
        gabBootsAOK_V1(bootNum,1) = sum(mean(gabBoot_V1(:,analysisStartBin:analysisEndBin)));
    else 
        lumBootsAOK_V1(bootNum,1) = sum(mean(-1*lumBoot_V1(:,analysisStartBin:analysisEndBin)));
        gabBootsAOK_V1(bootNum,1) = sum(mean(-1*gabBoot_V1(:,analysisStartBin:analysisEndBin)));
    end

end
lumAOK_V1_95CI(1,:) = quantile(lumBootsAOK_V1,[0.025 0.975]);
gabAOK_V1_95CI(1,:) = quantile(gabBootsAOK_V1,[0.025 0.975]);

%% Repeat for SC
limits.animal = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
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
U = selectUsingLimits(lumT, limits);
condition = ['SC' ' ' 'Lum'];

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
condition = ['SC' ' ' 'Gabor'];
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

if invert == 0
    % Compute Mean AOK during the analysis Window
    gabAOK_SC(1,1) = sum(mean(gabKernel(:,analysisStartBin:analysisEndBin)));
    lumAOK_SC(1,1) = sum(mean(lumKernel(:,analysisStartBin:analysisEndBin)));
else 
    gabAOK_SC(1,1) = sum(mean(-1*gabKernel(:,analysisStartBin:analysisEndBin)));
    lumAOK_SC(1,1) = sum(mean(-1*lumKernel(:,analysisStartBin:analysisEndBin)));
end
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
    lumSamps = [randsample(lum_nHit,lum_nHit,true)'...
        randsample([lum_nHit+1:lum_nTotal],lum_nTotal-lum_nHit,true)]';
    gabSamps = [randsample(gab_nHit,gab_nHit,true)'...
        randsample([gab_nHit+1:gab_nTotal],gab_nTotal - gab_nHit,true)]';
    % Take Samples w/ replacement
    lumBoot_SC = lumKernel(lumSamps,:);
    gabBoot_SC = gabKernel(gabSamps,:);
    
    if invert == 0
        % Compute 95%CI during the analysis Window From BootStraps
        lumBootsAOK_SC(bootNum,1) = sum(mean(lumBoot_SC(:,analysisStartBin:analysisEndBin)));
        gabBootsAOK_SC(bootNum,1) = sum(mean(gabBoot_SC(:,analysisStartBin:analysisEndBin)));
    else
        lumBootsAOK_SC(bootNum,1) = sum(mean(-1*lumBoot_SC(:,analysisStartBin:analysisEndBin)));
        gabBootsAOK_SC(bootNum,1) = sum(mean(-1*gabBoot_SC(:,analysisStartBin:analysisEndBin)));
    end

end

lumAOK_SC_95CI(1,:) = quantile(lumBootsAOK_SC,[0.025 0.975]);
gabAOK_SC_95CI(1,:) = quantile(gabBootsAOK_SC,[0.025 0.975]);

%% Make Scatter Plot of Results

V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];

figure;
hold on;
axis square;
if invert == 0 % Scatter Plot of AOK
    % Plot V1
    scatter(lumAOK_V1, gabAOK_V1, 100, V1Color, 'filled', 'DisplayName', 'V1');
    % Plot SC
    scatter(lumAOK_SC, gabAOK_SC, 100, SCColor, 'filled', 'DisplayName', 'SC');
    % Plot 95% CIs V1
    plot([lumAOK_V1_95CI(1,1)  lumAOK_V1_95CI(1,2)], [gabAOK_V1 gabAOK_V1], 'LineStyle','-', 'LineWidth',2, 'Color', V1Color);
    plot([lumAOK_V1 lumAOK_V1], [gabAOK_V1_95CI(1,1) gabAOK_V1_95CI(1,2)], 'LineStyle','-', 'LineWidth',2, 'Color', V1Color);
    % Plot 95% CIs SC
    plot([lumAOK_SC_95CI(1,1)  lumAOK_SC_95CI(1,2)], [gabAOK_SC gabAOK_SC], 'LineStyle','-', 'LineWidth',2, 'Color', SCColor);
    plot([lumAOK_SC lumAOK_SC], [gabAOK_SC_95CI(1,1) gabAOK_SC_95CI(1,2)], 'LineStyle','-', 'LineWidth',2, 'Color', SCColor);
    % Lines That Divide Quadrants into Effect on Performance
    plot([-4 4], [0 0], 'LineStyle','--', 'LineWidth',0.5, 'Color', 'k');
    plot([0 0], [-4 4], 'LineStyle','--', 'LineWidth',0.5, 'Color', 'k');
    xlabel('Luminance AOK');
    ylabel('Gabor AOK');
    title('Area Over Kernel');
    box off;
    set(gca, 'TickDir', 'out');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1);
    xlim([-3 3]);
    ylim([-3 3]);
    legend('V1', 'SC')
    hold off;
else % Bar plot
    b = bar([1,2,3,4], [lumAOK_V1, gabAOK_V1, lumAOK_SC, gabAOK_SC], 0.8);
    hold on;
    b.FaceColor = 'flat';
    % Color Bars
    b.CData(1,:) = V1Color;
    b.CData(2,:) = V1Color;
    b.CData(3,:) = SCColor;
    b.CData(4,:) = SCColor;
    % Plot CIs
    plot([1 1], [lumAOK_V1_95CI(1) lumAOK_V1_95CI(2)], 'LineWidth', 2, 'Color', 'k');
    plot([2 2], [gabAOK_V1_95CI(1) gabAOK_V1_95CI(2)], 'LineWidth', 2, 'Color', 'k');
    plot([3 3], [lumAOK_SC_95CI(1) lumAOK_SC_95CI(2)], 'LineWidth', 2, 'Color', 'k');
    plot([4 4], [gabAOK_SC_95CI(1) gabAOK_SC_95CI(2)], 'LineWidth', 2, 'Color', 'k');
    % Customize
    xlim([0.5 4.5]);
    ylim([-0.75 2.5]);
    %xlabel('Condition');
    ylabel('Area Over Kernel');
    title('AOK');
    box off;
    set(gca, 'TickDir', 'out');
    set(gca, 'FontSize', 14);
    set(gca, 'LineWidth', 1);
    set(gca, 'XTick', [1 2 3 4]);
    NL = "\newline";
    labelArray = {'Lum'+NL+' V1', 'Gabor'+NL+'  V1', 'Lum'+NL+' SC', 'Gabor'+NL+'  SC'};
    set(gca, 'XTickLabel', labelArray);
    set(gca,'YTick', [-0.5 0 2.5]);
    hold off;
end




