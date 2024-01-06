% Computes the Hit Rate and Kernel for Traces Where the opto stimulus is
% either on/off during the windowSpanMS
% Sweeps this computation over the entire -400 to 400 ms window surrounding
% the stimulus onset.

%% Set Up Analysis Span
windowSpanMS     = 5; % Tunes the Length over which to look for power on/off
                      % 25 ms was the size of a single opto frame
stimulusWindowMS = 800; % Window over which the stim kernels were computed
compBins         = stimulusWindowMS/windowSpanMS; % Number of Bins
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
%% Loop Over the Analysis Window, computing delta hit rate for each Bin

limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
animals = {'1960', '2015', '2083', '2126', '2220', '2221'};
condition = ['V1' ' ' 'Lum'];

% Column 1 = Power Off During Window
% Column 2 = Power On  During Window
V1Lum_hitRate = zeros(compBins,2,length(animals));

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(lumT, limits);

    % get all hit and miss Profiles from this mouse
    V1lumKernel_H = []; % Init
    V1lumKernel_M = []; % Init

    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        V1lumKernel_H = [V1lumKernel_H; stimProfiles.hitProfiles];
        V1lumKernel_M = [V1lumKernel_M; stimProfiles.missProfiles];
    end

    for binNum = 1:compBins
        startBin     = 1+(windowSpanMS*(binNum-1));
        endBin       = (windowSpanMS*binNum);

        % List of trials for hits/misses with power on during critical window
        hitListON   = any(V1lumKernel_H(:,startBin:endBin) == 1,2);
        hitListOFF  = any(V1lumKernel_H(:,startBin:endBin) ~= 1,2);
        missListON  = any(V1lumKernel_M(:,startBin:endBin) == 1,2);
        missListOFF = any(V1lumKernel_M(:,startBin:endBin) ~= 1,2);

        % Power Off During Window
        V1Lum_hitRate(binNum,1,a) = 100*(sum(hitListOFF)/(sum(hitListOFF)+sum(missListOFF)));
        % POwer ON During Window
        V1Lum_hitRate(binNum,2,a) = 100*(sum(hitListON)/(sum(hitListON)+sum(missListON)));
    end
end


%% Get StimProfiles and Compute AOK For V1 Gabor
animals = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
condition = ['V1' ' ' 'Gabor'];

% Column 1 = Power Off During Window
% Column 2 = Power On  During Window
V1Gab_hitRate = zeros(compBins,2,length(animals));

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(gabT, limits);

    V1gabKernel_H = []; % Init
    V1gabKernel_M = []; % Init

    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        V1gabKernel_H = [V1gabKernel_H; stimProfiles.hitProfiles];
        V1gabKernel_M = [V1gabKernel_M; stimProfiles.missProfiles];
    end

    for binNum = 1:compBins
        startBin     = 1+(windowSpanMS*(binNum-1));
        endBin       = (windowSpanMS*binNum);

        % List of trials for hits/misses with power on during critical window
        hitListON   = any(V1gabKernel_H(:,startBin:endBin) == 1,2);
        hitListOFF  = any(V1gabKernel_H(:,startBin:endBin) ~= 1,2);
        missListON  = any(V1gabKernel_M(:,startBin:endBin) == 1,2);
        missListOFF = any(V1gabKernel_M(:,startBin:endBin) ~= 1,2);

        % Power Off During Window
        V1Gab_hitRate(binNum,1,a) = 100*(sum(hitListOFF)/(sum(hitListOFF)+sum(missListOFF)));
        % POwer ON During Window
        V1Gab_hitRate(binNum,2,a) = 100*(sum(hitListON)/(sum(hitListON)+sum(missListON)));
    end

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

% Column 1 = Power Off During Window
% Column 2 = Power On  During Window
SCLum_hitRate = zeros(compBins,2,length(animals));

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(lumT, limits);

    SClumKernel_H = []; % Init
    SClumKernel_M = []; % Init

    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        SClumKernel_H = [SClumKernel_H; stimProfiles.hitProfiles];
        SClumKernel_M = [SClumKernel_M; stimProfiles.missProfiles];
    end

    for binNum = 1:compBins
        startBin     = 1+(windowSpanMS*(binNum-1));
        endBin       = (windowSpanMS*binNum);

        % List of trials for hits/misses with power on during critical window
        hitListON   = any(SClumKernel_H(:,startBin:endBin) == 1,2);
        hitListOFF  = any(SClumKernel_H(:,startBin:endBin) ~= 1,2);
        missListON  = any(SClumKernel_M(:,startBin:endBin) == 1,2);
        missListOFF = any(SClumKernel_M(:,startBin:endBin) ~= 1,2);

        % Power Off During Window
        SCLum_hitRate(binNum,1,a) = 100*(sum(hitListOFF)/(sum(hitListOFF)+sum(missListOFF)));
        % POwer ON During Window
        SCLum_hitRate(binNum,2,a) = 100*(sum(hitListON)/(sum(hitListON)+sum(missListON)));

    end
end


%% SC Gabor
animals = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
condition = ['SC' ' ' 'Gabor'];

% Column 1 = Power Off During Window
% Column 2 = Power On  During Window
SCGab_hitRate = zeros(compBins,2,length(animals));

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(gabT, limits);

    SCgabKernel_H = []; % Init
    SCgabKernel_M = []; % Init    
    
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        SCgabKernel_H = [SCgabKernel_H; stimProfiles.hitProfiles];
        SCgabKernel_M = [SCgabKernel_M; stimProfiles.missProfiles];
    end

    for binNum = 1:compBins
        startBin     = 1+(windowSpanMS*(binNum-1));
        endBin       = (windowSpanMS*binNum);

        % List of trials for hits/misses with power on during critical window
        hitListON   = any(SCgabKernel_H(:,startBin:endBin) == 1,2);
        hitListOFF  = any(SCgabKernel_H(:,startBin:endBin) ~= 1,2);
        missListON  = any(SCgabKernel_M(:,startBin:endBin) == 1,2);
        missListOFF = any(SCgabKernel_M(:,startBin:endBin) ~= 1,2);

        % Power Off During Window
        SCGab_hitRate(binNum,1,a) = 100*(sum(hitListOFF)/(sum(hitListOFF)+sum(missListOFF)));
        % POwer ON During Window
        SCGab_hitRate(binNum,2,a) = 100*(sum(hitListON)/(sum(hitListON)+sum(missListON)));
    end
end

%% Line Plots of Hit Rates Over Time
% Hit Rate Normalized To No-Power Condition

% Colors
V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];

% Common y-scale for all plots
yBottom = 0.90;
yTop    = 1.10;
plotBins = 1:windowSpanMS:800;
fillX = [plotBins, fliplr(plotBins)];

% Mean 
gabMean = mean(V1Gab_hitRate(:,2,:)./V1Gab_hitRate(:,1,:),3);
lumMean = mean(V1Lum_hitRate(:,2,:)./V1Lum_hitRate(:,1,:),3);
% SEM
gabSEM  = (std(squeeze(V1Gab_hitRate(:,2,:)./V1Gab_hitRate(:,1,:)),0,2)/sqrt(size(V1Gab_hitRate,3)));
lumSEM  = (std(squeeze(V1Lum_hitRate(:,2,:)./V1Lum_hitRate(:,1,:)),0,2)/sqrt(size(V1Lum_hitRate,3)));

% rolling t-test, difference from 1
v1l = squeeze(V1Lum_hitRate(:,2,:)./V1Lum_hitRate(:,1,:));
v1g = squeeze(V1Gab_hitRate(:,2,:)./V1Gab_hitRate(:,1,:));
v1l_P = zeros(length(plotBins),1);
v1g_P = zeros(length(plotBins),1);
for i = 1:length(plotBins)
    [v1l_P(i,1),~] = ttest(v1l(i,:),1,"Alpha",0.05);
    [v1g_P(i,1),~] = ttest(v1g(i,:),1,"Alpha",0.05);
end

% Set Non-Significant Bins to 0
v1g_P(v1g_P==0) = nan;
v1l_P(v1l_P==0) = nan;


figure('Position',[10 10 800 1000]);
subplot(2,2,1);
hold on;
axis square;
% Horizontal Line at y = 1
plot([1 800],[1 1], 'LineWidth',0.5, 'Color','k', 'LineStyle','--')
% LinePlot
plot(plotBins, gabMean, 'LineWidth',1,'Color',V1Color);
% SEM Fill
gabfillCI = [gabMean-gabSEM; flipud(gabMean+gabSEM)]';
fill(fillX, gabfillCI, V1Color, 'lineStyle', '-', 'edgeColor', V1Color, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
% Plot Significant Bins
scatter(plotBins, 1.08*v1g_P, 60,'black','filled');

% Customize
ylim([yBottom yTop]);
title('V1 Contrast: Hit Rate Relative to LED OFF');
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xticks([0 400 800]);
xticklabels({'-400', '0', '400'});
ax = gca;
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '--';
ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
xlabel('Time (ms)');
ylabel('Hit Rate Relative to LED OFF');

% V1 Luminance
subplot(2,2,2);
hold on;
axis square;
% Horizontal Line at y = 1
plot([1 800],[1 1], 'LineWidth',0.5, 'Color','k', 'LineStyle','--');
% LinePlot
plot(1:5:800, lumMean, 'LineWidth',1,'Color',V1Color);
% SEM Fill
fillX = [plotBins, fliplr(plotBins)];
lumfillCI = [lumMean-lumSEM; flipud(lumMean+lumSEM)]';
fill(fillX, lumfillCI, V1Color, 'lineStyle', '-', 'edgeColor', V1Color, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
% Plot Significant Bins
scatter(plotBins, 1.08*v1l_P, 60,'black','filled');

% Customize
ylim([yBottom yTop]);
title('V1 Luminance: Hit Rate Relative to LED OFF');
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xticks([0 400 800]);
xticklabels({'-400', '0', '400'});
ax = gca;
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '--';
ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
xlabel('Time (ms)');
ylabel('Hit Rate Relative to LED OFF');

% Now SC Data

% Mean 
gabMean = mean(SCGab_hitRate(:,2,:)./SCGab_hitRate(:,1,:),3);
lumMean = mean(SCLum_hitRate(:,2,:)./SCLum_hitRate(:,1,:),3);
% SEM
gabSEM  = (std(squeeze(SCGab_hitRate(:,2,:)./SCGab_hitRate(:,1,:)),0,2)/sqrt(size(SCGab_hitRate,3)));
lumSEM  = (std(squeeze(SCLum_hitRate(:,2,:)./SCLum_hitRate(:,1,:)),0,2)/sqrt(size(SCLum_hitRate,3)));

% rolling t-test, difference from 1
scl = squeeze(SCLum_hitRate(:,2,:)./SCLum_hitRate(:,1,:));
scg = squeeze(SCGab_hitRate(:,2,:)./SCGab_hitRate(:,1,:));
scl_P = zeros(length(plotBins),1);
scg_P = zeros(length(plotBins),1);
for i = 1:length(plotBins)
    [scl_P(i,1),~] = ttest(scl(i,:),1,"Alpha",0.05);
    [scg_P(i,1),~] = ttest(scg(i,:),1,"Alpha",0.05);
end

% Set Non-Significant Bins to 0
scg_P(scg_P==0) = nan;
scl_P(scl_P==0) = nan;

% Repeated Measures ANOVA
scGab = squeeze(SCGab_hitRate(:,2,:)./SCGab_hitRate(:,1,:));
scLum = squeeze(SCLum_hitRate(:,2,:)./SCLum_hitRate(:,1,:));




subplot(2,2,3);
hold on;
axis square;
% Horizontal Line at y = 1
plot([1 800],[1 1], 'LineWidth',0.5, 'Color','k', 'LineStyle','--');
% LinePlot
plot(1:5:800, gabMean, 'LineWidth',1,'Color',SCColor);
% SEM Fill
fillX = [plotBins, fliplr(plotBins)];
gabfillCI = [gabMean-gabSEM; flipud(gabMean+gabSEM)]';
fill(fillX, gabfillCI, SCColor, 'lineStyle', '-', 'edgeColor', SCColor, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
% Plot Significant Bins
scatter(plotBins, 1.08*scg_P, 60,'black','filled');

% Customize
ylim([yBottom yTop]);
title('SC Contrast: Hit Rate Relative to LED OFF');
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xticks([0 400 800]);
xticklabels({'-400', '0', '400'});
ax = gca;
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '--';
ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
xlabel('Time (ms)');
ylabel('Hit Rate Relative to LED OFF');

% SC Luminance
subplot(2,2,4);
hold on;
axis square;
% Horizontal Line at y = 1
plot([1 800],[1 1], 'LineWidth',0.5, 'Color','k', 'LineStyle','--');
% LinePlot
plot(1:5:800, lumMean, 'LineWidth',2,'Color',SCColor);
% SEM Fill
fillX = [plotBins, fliplr(plotBins)]';
lumfillCI = [lumMean-lumSEM; flipud(lumMean+lumSEM)]';
fill(fillX, lumfillCI, SCColor, 'lineStyle', '-', 'edgeColor', SCColor, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
% Plot Significant Bins
scatter(plotBins, 1.08*scl_P, 60,'black','filled');
% Customize
ylim([yBottom yTop]);
title('SC Luminance: Hit Rate Relative to LED OFF');
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xticks([0 400 800]);
xticklabels({'-400', '0', '400'});
ax = gca;
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '--';
ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
xlabel('Time (ms)');
ylabel('Hit Rate Relative to LED OFF');


