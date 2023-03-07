function [BootsAOK_V1, lumBootsAOK_SC, gabBootsAOK_V1, gabBootsAOK_SC] = ...
    areaOverKernel_GrandAv_OFFSET(invert, bootSamps, analysisStartMS, analysisDurMS)
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
V1Folder = strcat(analysisDir, 'V1',' Offset/');
SCFolder = strcat(analysisDir, 'SC',' Offset/');

% Load Tables for each stimulus type
load([V1Folder 'masterTable.mat'], 'T');
v1T = T;
load([SCFolder 'masterTable.mat'], 'T');
scT = T;
clear T;

%% Grab all the stim Profiles one mouse then create kernel
% Use Limits to Subselect Sessions From Each Mouse
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;

% Init
AOK_V1 = zeros(1,1);
AOK_V1_95CI  = zeros(1,2);
%% V1 Offset
limits.animal = {'2016', '2083', '2126', '2220', '2221'};
U = selectUsingLimits(v1T, limits);
condition = ['V1' ' ' 'Offset'];

% Extract optoProfiles
v1Profiles.hitProfiles = [];
v1Profiles.missProfiles = [];
v1Profiles.earlyProfiles = [];
v1Profiles.RTProfiles = [];
v1Profiles.stimRTProfiles = [];

% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
    v1Profiles.hitProfiles = [v1Profiles.hitProfiles; stimProfiles.hitProfiles];
    v1Profiles.missProfiles = [v1Profiles.missProfiles; stimProfiles.missProfiles];
    v1Profiles.earlyProfiles = [v1Profiles.earlyProfiles; stimProfiles.earlyProfiles];
    v1Profiles.RTProfiles = [v1Profiles.RTProfiles; stimProfiles.RTProfiles];
    v1Profiles.stimRTProfiles = [v1Profiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end

% total Kernel
v1Kernel = [v1Profiles.hitProfiles; -v1Profiles.missProfiles];

% Compute Mean AOK during the analysis Window
if invert == 0
    AOK_V1(1,1) = sum(mean(v1Kernel(:,analysisStartBin:analysisEndBin)));
else
    AOK_V1(1,1) = sum(mean(-1*v1Kernel(:,analysisStartBin:analysisEndBin)));
end
%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
nHit = size(v1Profiles.hitProfiles,1);
v1_nTotal = size(v1Kernel,1);

% Init Archives For BootStrap Samples
BootsAOK_V1 = zeros(bootSamps, 1);

for bootNum = 1:bootSamps
    % Samps For This Round of BootStrapping
    v1Samps = [randsample(nHit,nHit,true)'...
        randsample([nHit+1:v1_nTotal],v1_nTotal-nHit,true)]';
    % Take Samples w/ replacement
    Boot_V1 = v1Kernel(v1Samps,:);

    if invert == 0
        % Compute 95%CI during the analysis Window From BootStraps
        BootsAOK_V1(bootNum,1) = sum(mean(Boot_V1(:,analysisStartBin:analysisEndBin)));
    else 
        BootsAOK_V1(bootNum,1) = sum(mean(-1*Boot_V1(:,analysisStartBin:analysisEndBin)));
    end

end
AOK_V1_95CI(1,:) = quantile(BootsAOK_V1,[0.025 0.975]);

%% Repeat for SC
limits.animal = {'1674', '1675', '1902', '2057', '2063', '2236'};
SCFolder = strcat(analysisDir, 'SC',' Offset/');
load([SCFolder 'masterTable.mat'], 'T');
scT = T;
clear T;

% Init
AOK_SC = zeros(1,1);
AOK_SC_95CI  = zeros(1,2);

%% Both Gabor and Luminance
% First For Luminance
U = selectUsingLimits(scT, limits);
condition = ['SC' ' ' 'Offset'];

% Extract optoProfiles
scProfiles.hitProfiles = [];
scProfiles.missProfiles = [];
scProfiles.earlyProfiles = [];
scProfiles.RTProfiles = [];
scProfiles.stimRTProfiles = [];

% get all Profiles from this mouse
for i = 1:height(U)
    load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
        condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
    scProfiles.hitProfiles = [scProfiles.hitProfiles; stimProfiles.hitProfiles];
    scProfiles.missProfiles = [scProfiles.missProfiles; stimProfiles.missProfiles];
    scProfiles.earlyProfiles = [scProfiles.earlyProfiles; stimProfiles.earlyProfiles];
    scProfiles.RTProfiles = [scProfiles.RTProfiles; stimProfiles.RTProfiles];
    scProfiles.stimRTProfiles = [scProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end

% Collate All Traces From This Mouse
scKernel = [scProfiles.hitProfiles; -scProfiles.missProfiles];

if invert == 0
    % Compute Mean AOK during the analysis Window
    AOK_SC(1,1) = sum(mean(scKernel(:,analysisStartBin:analysisEndBin)));
else 
    AOK_SC(1,1) = sum(mean(-1*scKernel(:,analysisStartBin:analysisEndBin)));
end
%% Bootstrap CIs from Stim Profiles (bootNum = num bootstraps)
% How Many of Each Outcome Type to control bootstrapping to match
% experimental proportions
nHit = size(scProfiles.hitProfiles,1);
sc_nTotal = size(scKernel,1);

% Init Archives For BootStrap Samples
BootsAOK_SC = zeros(bootSamps, 1);

for bootNum = 1:bootSamps
    % Samps For This Round of BootStrapping
    scSamps = [randsample(nHit,nHit,true)'...
        randsample([nHit+1:sc_nTotal],sc_nTotal-nHit,true)]';
    % Take Samples w/ replacement
    Boot_SC = scKernel(scSamps,:);
    
    if invert == 0
        % Compute 95% CI during the analysis Window From BootStraps
        BootsAOK_SC(bootNum,1) = sum(mean(Boot_SC(:,analysisStartBin:analysisEndBin)));
    else
        BootsAOK_SC(bootNum,1) = sum(mean(-1*Boot_SC(:,analysisStartBin:analysisEndBin)));
    end

end

AOK_SC_95CI(1,:) = quantile(BootsAOK_SC,[0.025 0.975]);

%% Make Scatter Plot of Results

V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];

figure;
hold on;
axis square;
if invert == 0 % Scatter Plot of AOK
    % Plot V1
    scatter(lumAOK_V1, AOK_V1, 100, V1Color, 'filled', 'DisplayName', 'V1');
    % Plot SC
    scatter(lumAOK_SC, AOK_SC, 100, SCColor, 'filled', 'DisplayName', 'SC');
    % Plot 95% CIs V1
    plot([lumAOK_V1_95CI(1,1)  lumAOK_V1_95CI(1,2)], [AOK_V1 AOK_V1], 'LineStyle','-', 'LineWidth',2, 'Color', V1Color);
    plot([lumAOK_V1 lumAOK_V1], [AOK_V1_95CI(1,1) AOK_V1_95CI(1,2)], 'LineStyle','-', 'LineWidth',2, 'Color', V1Color);
    % Plot 95% CIs SC
    plot([lumAOK_SC_95CI(1,1)  lumAOK_SC_95CI(1,2)], [AOK_SC AOK_SC], 'LineStyle','-', 'LineWidth',2, 'Color', SCColor);
    plot([lumAOK_SC lumAOK_SC], [AOK_SC_95CI(1,1) AOK_SC_95CI(1,2)], 'LineStyle','-', 'LineWidth',2, 'Color', SCColor);
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
    b = bar([1,2,3,4], [AOK_V1, lumAOK_V1, AOK_SC, lumAOK_SC], 0.8);
    hold on;
    b.FaceColor = 'flat';
    % Color Bars
    b.CData(1,:) = V1Color;
    b.CData(2,:) = V1Color;
    b.CData(3,:) = SCColor;
    b.CData(4,:) = SCColor;
    % Plot CIs
    plot([1 1], [AOK_V1_95CI(1) AOK_V1_95CI(2)], 'LineWidth', 2, 'Color', 'k');
    plot([2 2], [lumAOK_V1_95CI(1) lumAOK_V1_95CI(2)], 'LineWidth', 2, 'Color', 'k');
    plot([3 3], [AOK_SC_95CI(1) AOK_SC_95CI(2)], 'LineWidth', 2, 'Color', 'k');
    plot([4 4], [lumAOK_SC_95CI(1) lumAOK_SC_95CI(2)], 'LineWidth', 2, 'Color', 'k');
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

%% Histograms
% Colors
V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];
bins = -3:0.2:3;
% Figure
figure('Position', [10 10 1000 500]);
axis square;
hold on;
histogram(BootsAOK_V1,bins, 'FaceColor', V1Color, 'FaceAlpha', 0.5, 'Normalization','probability');
histogram(BootsAOK_SC,bins, 'FaceColor', SCColor, 'FaceAlpha', 0.5, 'Normalization','probability');
title('Gabor: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
ylabel('Probability');
xlabel('Area Over The Kernel');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
ylim([0 0.35]);
xlim([-2.5 2.5]);
legend('V1', 'SC', 'Location','northwest');
hold off;


% Different from 0?
signrank(BootsAOK_SC)
signrank(BootsAOK_V1)

median(BootsAOK_V1)
