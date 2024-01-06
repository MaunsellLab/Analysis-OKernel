%;; Computes the Hit Rate and Kernel for Traces Where the opto stimulus is
% either on/off during the windowSpanMS
windowSpanMS = 5; % Tunes the Length over which to look for power on/off
        % 25 ms was the size of a single opto frame
stimOnBin    = 401; % Bin Of the Kernels where stim turns on;
centerBin    = 0; % Center of Test Window (average peak time of NBK == 51 ms across conditions).
leftSpan     = floor(windowSpanMS/2); % 
rightSpan    = ceil(windowSpanMS/2);
startBin     = stimOnBin + centerBin - leftSpan; % Start of Analysis Window after stim onset
endBin       = stimOnBin + centerBin + rightSpan; % End of Analysis Windown after stim onset
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

% Column 1 = Power Off During Window
% Column 2 = Power On  During Window
V1Lum_hitRate = zeros(length(animals),2);

% Master List Across All Animals
V1lumKernel_ON = [];
V1lumKernel_OFF = [];

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(lumT, limits);
    
    V1lumKernel_H = []; % Init
    V1lumKernel_M = []; % Init
   
    % get all hit and miss Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        V1lumKernel_H = [V1lumKernel_H; stimProfiles.hitProfiles];
        V1lumKernel_M = [V1lumKernel_M; stimProfiles.missProfiles];
    end
    
    % List of trials for hits/misses with power on during critical window
    hitListON   = any(V1lumKernel_H(:,startBin:endBin) == 1,2);
    hitListOFF  = any(V1lumKernel_H(:,startBin:endBin) ~= 1,2);
    missListON  = any(V1lumKernel_M(:,startBin:endBin) == 1,2);
    missListOFF = any(V1lumKernel_M(:,startBin:endBin) ~= 1,2);

    % Power Off During Window
    V1Lum_hitRate(a,1) = 100*(sum(hitListOFF)/(sum(hitListOFF)+sum(missListOFF)));
    % POwer ON During Window
    V1Lum_hitRate(a,2) = 100*(sum(hitListON)/(sum(hitListON)+sum(missListON)));

    % Add Kernels to Master List
    V1lumKernel_ON  = [V1lumKernel_ON; [V1lumKernel_H(hitListON,:); -V1lumKernel_M(missListON,:)]];
    V1lumKernel_OFF = [V1lumKernel_OFF; [V1lumKernel_H(hitListOFF,:); -V1lumKernel_M(missListOFF,:)]];
end


%% Get StimProfiles and Compute AOK For V1 Gabor
animals = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
condition = ['V1' ' ' 'Gabor'];

% Column 1 = Power Off During Window
% Column 2 = Power On  During Window
V1Gab_hitRate = zeros(length(animals),2);

% Master List Across All Animals
V1gabKernel_ON = [];
V1gabKernel_OFF = [];

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
   
    % List of trials for hits/misses with power on during critical window
    hitListON   = any(V1gabKernel_H(:,startBin:endBin) == 1,2);
    hitListOFF  = any(V1gabKernel_H(:,startBin:endBin) ~= 1,2);
    missListON  = any(V1gabKernel_M(:,startBin:endBin) == 1,2);
    missListOFF = any(V1gabKernel_M(:,startBin:endBin) ~= 1,2);

    % Power Off During Window
    V1Gab_hitRate(a,1) = 100*(sum(hitListOFF)/(sum(hitListOFF)+sum(missListOFF)));
    % POwer ON During Window
    V1Gab_hitRate(a,2) = 100*(sum(hitListON)/(sum(hitListON)+sum(missListON)));

    % Add Kernels to Master List
    V1gabKernel_ON  = [V1gabKernel_ON; [V1gabKernel_H(hitListON,:); -V1gabKernel_M(missListON,:)]];
    V1gabKernel_OFF = [V1gabKernel_OFF; [V1gabKernel_H(hitListOFF,:); -V1gabKernel_M(missListOFF,:)]];
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
SCLum_hitRate = zeros(length(animals),2);

% Master List Across All Animals
SClumKernel_ON = [];
SClumKernel_OFF = [];

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
    
    % List of trials for hits/misses with power on during critical window
    hitListON   = any(SClumKernel_H(:,startBin:endBin) == 1,2);
    hitListOFF  = any(SClumKernel_H(:,startBin:endBin) ~= 1,2);
    missListON  = any(SClumKernel_M(:,startBin:endBin) == 1,2);
    missListOFF = any(SClumKernel_M(:,startBin:endBin) ~= 1,2);

    % Power Off During Window
    SCLum_hitRate(a,1) = 100*(sum(hitListOFF)/(sum(hitListOFF)+sum(missListOFF)));
    % POwer ON During Window
    SCLum_hitRate(a,2) = 100*(sum(hitListON)/(sum(hitListON)+sum(missListON)));

    % Add Kernels to Master List
    SClumKernel_ON  = [SClumKernel_ON; [SClumKernel_H(hitListON,:); -SClumKernel_M(missListON,:)]];
    SClumKernel_OFF = [SClumKernel_OFF; [SClumKernel_H(hitListOFF,:); -SClumKernel_M(missListOFF,:)]];

end

%% SC Gabor
animals = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
condition = ['SC' ' ' 'Gabor'];

% Column 1 = Power Off During Window
% Column 2 = Power On  During Window
SCGab_hitRate = zeros(length(animals),2);

% Master List Across All Animals
SCgabKernel_ON = [];
SCgabKernel_OFF = [];

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
    
    % List of trials for hits/misses with power on during critical window
    hitListON   = any(SCgabKernel_H(:,startBin:endBin) == 1,2);
    hitListOFF  = any(SCgabKernel_H(:,startBin:endBin) ~= 1,2);
    missListON  = any(SCgabKernel_M(:,startBin:endBin) == 1,2);
    missListOFF = any(SCgabKernel_M(:,startBin:endBin) ~= 1,2);

    % Power Off During Window
    SCGab_hitRate(a,1) = 100*(sum(hitListOFF)/(sum(hitListOFF)+sum(missListOFF)));
    % POwer ON During Window
    SCGab_hitRate(a,2) = 100*(sum(hitListON)/(sum(hitListON)+sum(missListON)));

    % Add Kernels to Master List
    SCgabKernel_ON  = [SClumKernel_ON; [SCgabKernel_H(hitListON,:); -SCgabKernel_M(missListON,:)]];
    SCgabKernel_OFF = [SClumKernel_OFF; [SCgabKernel_H(hitListOFF,:); -SCgabKernel_M(missListOFF,:)]];

end

%% Bar Plots with scatter overlaid

% Hit Rate Normalized To No-Power Condition


% Common y-scale for all plots
yBottom = 0.85;
yTop    = 1.05;

figure('Position',[10 10 800 1000]);
% V1 Contrast
subplot(2,2,1);
hold on;
% Bar Means
bar([1 2], mean(V1Gab_hitRate./V1Gab_hitRate(:,1),1), 'FaceColor', 'w', 'LineWidth',1, 'BarWidth', 0.4);
title('V1 Contrast: Off vs. On');
ylabel('Hit Rate Relative to LED OFF');
ylim([yBottom yTop]);
xticks([1 2]);
xticklabels({'LED OFF', 'LED ON'});
% Connecting Lines, one per mouse
for i = 1:size(V1Gab_hitRate,1)
    line([1 2], [V1Gab_hitRate(i,1)/V1Gab_hitRate(i,1), V1Gab_hitRate(i,2)/V1Gab_hitRate(i,1)], 'LineWidth', 0.5, 'Color', '[0.5 0.5 0.5]');
end
% Individual Mice values
scatter(ones(size(V1Gab_hitRate,1),1),[V1Gab_hitRate(:,1)./V1Gab_hitRate(:,1)],70, 'k');
scatter(1+ones(size(V1Gab_hitRate,1),1),[V1Gab_hitRate(:,2)./V1Gab_hitRate(:,1)],70, 'k');
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
yticks([0.9 1.0 1.1]);
xlim([0.5 2.5]);
box off;
hold off;

% SC Contrast
subplot(2,2,3);
hold on;
% Bar Means
bar([1 2], mean(SCGab_hitRate./SCGab_hitRate(:,1),1), 'FaceColor', 'w', 'LineWidth',1, 'BarWidth', 0.4);
title('SC Contrast: Off vs. On');
ylabel('Hit Rate Relative to LED OFF');
ylim([yBottom yTop]);
xticks([1 2]);
xticklabels({'LED OFF', 'LED ON'});
% Connecting Lines, one per mouse
for i = 1:size(SCGab_hitRate,1)
    line([1 2], [SCGab_hitRate(i,1)/SCGab_hitRate(i,1), SCGab_hitRate(i,2)/SCGab_hitRate(i,1)], 'LineWidth', 0.5, 'Color', '[0.5 0.5 0.5]');
end
% Individual Mice values
scatter(ones(size(SCGab_hitRate,1),1),[SCGab_hitRate(:,1)./SCGab_hitRate(:,1)],70, 'k');
scatter(1+ones(size(SCGab_hitRate,1),1),[SCGab_hitRate(:,2)./SCGab_hitRate(:,1)],70, 'k');
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
yticks([0.9 1.0 1.1]);
xlim([0.5 2.5]);
box off;
hold off;


% V1 Luminance
subplot(2,2,2);
hold on;
% Bar Means
bar([1 2], mean(V1Lum_hitRate./V1Lum_hitRate(:,1),1), 'FaceColor', 'w', 'LineWidth',1, 'BarWidth', 0.4);
title('V1 Luminance: Off vs. On');
ylabel('Hit Rate Relative to LED OFF');
ylim([yBottom yTop]);
xticks([1 2]);
xticklabels({'LED OFF', 'LED ON'});
% Connecting Lines, one per mouse
for i = 1:size(V1Lum_hitRate,1)
    line([1 2], [V1Lum_hitRate(i,1)/V1Lum_hitRate(i,1), V1Lum_hitRate(i,2)/V1Lum_hitRate(i,1)], 'LineWidth', 0.5, 'Color', '[0.5 0.5 0.5]');
end
% Individual Mice values
scatter(ones(size(V1Lum_hitRate,1),1),[V1Lum_hitRate(:,1)./V1Lum_hitRate(:,1)],70, 'k');
scatter(1+ones(size(V1Lum_hitRate,1),1),[V1Lum_hitRate(:,2)./V1Lum_hitRate(:,1)],70, 'k');
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
yticks([0.9 1.0 1.1]);
xlim([0.5 2.5]);
box off;
hold off;

% SC Luminance
subplot(2,2,4);
hold on;
% Bar Means
bar([1 2], mean(SCLum_hitRate./SCLum_hitRate(:,1),1), 'FaceColor', 'w', 'LineWidth',1, 'BarWidth', 0.4);
title('SC Luminance: Off vs. On');
ylabel('Hit Rate Relative to LED OFF');
ylim([yBottom yTop]);
xticks([1 2]);
xticklabels({'LED OFF', 'LED ON'});
% Connecting Lines, one per mouse
for i = 1:size(SCLum_hitRate,1)
    line([1 2], [SCLum_hitRate(i,1)/SCLum_hitRate(i,1), SCLum_hitRate(i,2)/SCLum_hitRate(i,1)], 'LineWidth', 0.5, 'Color', '[0.5 0.5 0.5]');
end
% Individual Mice values
scatter(ones(size(SCLum_hitRate,1),1),[SCLum_hitRate(:,1)./SCLum_hitRate(:,1)],70, 'k');
scatter(1+ones(size(SCLum_hitRate,1),1),[SCLum_hitRate(:,2)./SCLum_hitRate(:,1)],70, 'k');
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
yticks([0.9 1.0 1.1]);
xlim([0.5 2.5]);
box off;
hold off;

% Stats On Difference From 1
[~,p1] = ttest(V1Gab_hitRate(:,2)./V1Gab_hitRate(:,1),1,'Alpha',0.05);
[~,p2] = ttest(SCGab_hitRate(:,2)./SCGab_hitRate(:,1),1,'Alpha',0.05);
[~,p3] = ttest(V1Lum_hitRate(:,2)./V1Lum_hitRate(:,1),1,'Alpha',0.05);
[~,p4] = ttest(SCLum_hitRate(:,2)./SCLum_hitRate(:,1),1,'Alpha',0.05);

%% Bootstrapped Kernel Profile of Selected Trials

% Filter For Kernels
sampleFreqHz = 1000;
filterLP = designfilt('lowpassfir', 'PassbandFrequency', 90 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 90 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod','equiripple');

% BootStrapped ON
BootKernel_ON  = bootstrp(250, @mean, SClumKernel_ON);
ONPCs = prctile(BootKernel_ON, [15.9, 50, 84.1]);              % +/- 1 SEM
ONPCMeans = mean(ONPCs, 2);
% Bootstrapped OFF
BootKernel_OFF = bootstrp(250, @mean, SClumKernel_OFF);
OFFPCs = prctile(BootKernel_OFF, [15.9, 50, 84.1]);              % +/- 1 SEM
OFFPCMeans = mean(OFFPCs, 2);

% CIs
ONCIs = zeros(3, size(BootKernel_ON, 2));
OFFCIs = zeros(3, size(BootKernel_OFF, 2));
for c = 1:3
    ONCIs(c, :) = filtfilt(filterLP, ONPCs(c, :) - ONPCMeans(c)) + ONPCMeans(c);
    OFFCIs(c, :) = filtfilt(filterLP, OFFPCs(c, :) - OFFPCMeans(c)) + OFFPCMeans(c);
end
%% Figure Showing Kernel Profile of Selected Trials

% Figure Set Up
bins = 800;
ONx  = 1:size(ONCIs, 2);
OFFx = 1:size(OFFCIs,2); 
[plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
yLabel = 'Normalized Power';

figure('Units', 'inches', 'Position', [3, 1, 8, 8]);
axis square
plot([0, bins], [0, 0], 'Color', 'k', 'LineStyle','--');
hold on;
% Kernel for OFF Trials
plot(OFFx, OFFCIs(2, :), 'Color', '[0.5 0.5 0.5]', 'LineStyle', '-', 'LineWidth',2);
OFFx2 = [OFFx, fliplr(OFFx)];
OFFfillCI = [OFFCIs(1, :), fliplr(OFFCIs(3, :))];
fill(OFFx2, OFFfillCI, 'k', 'lineStyle', '-', 'edgeColor', '[0.5, 0.5, 0.5]', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
% Kernel for ON Trials
plot(ONx, ONCIs(2, :), 'Color', 'k', 'LineStyle', '-', 'LineWidth',2);
ONx2 = [ONx, fliplr(ONx)];
ONfillCI = [ONCIs(1, :), fliplr(ONCIs(3, :))];
fill(ONx2, ONfillCI, 'k', 'lineStyle', '-', 'edgeColor', 'k', 'edgeAlpha', 0.5, 'faceAlpha', 0.10);

ax = gca;
ax.YColor = [0 0 0];
ax.TickDir = 'out';
ax.XGrid = 'on';
ax.FontSize = 14;
xlim(ax, [0, bins]);
ylabel(yLabel);
set(gca,'XTick', [0, -plotStartMS, bins]);
set(gca, 'XTickLabel', {sprintf('%d', plotStartMS), '0', sprintf('%d', plotEndMS)});
set(gca, 'LineWidth', 1);
set(gca,'YTick', [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]);
ax.XMinorGrid = 'on';
ax.MinorGridLineStyle = '--';
ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
xlabel('Time Relative to Stimulus');
ylim([-0.3 0.3]);
legend('','LED OFF', '', 'LED ON');
hold off;
