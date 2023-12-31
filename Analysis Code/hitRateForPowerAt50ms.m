% Computes the Hit Rate and Kernel for Traces Where the opto stimulus is
% either on/off during the windowSpanMS
StimOnBin = 401; % Bin Of the Kernels where stim turns on;
windowSpanMS = 10; % Tunes the Length over which to look for power on/off
startBin = StimOnBin + 45; % Start of Analysis Window after stim onset
endBin   = startBin + windowSpanMS; 
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



%% Bar Plots



% plot
xs = -5:0.01:5;

figure('Position',[10 10 800 1000]);
subplot(2,2,1);
hold on;
scatter(V1GabAOKs, V1GabDeltaDs, 50, 'black', 'filled');
plot(xs, f1.Coefficients.Estimate(1) + xs.*f1.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Gabor: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-5 5]);
ylim([-0.5 0.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-5 0 5]);
set(gca, 'YTick', [-0.8 -0.4 0 0.4 0.8]);
box off;
hold off;

subplot(2,2,2);
hold on;
scatter(V1LumAOKs, V1LumDeltaDs, 50, 'black', 'filled');
plot(xs, f3.Coefficients.Estimate(1) + xs.*f3.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('V1 Luminance: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-5 5]);
ylim([-0.5 0.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-5 0 5]);
set(gca, 'YTick', [-0.8 -0.4 0 0.4 0.8]);
box off;
hold off;


subplot(2,2,3);
hold on;
scatter(SCGabAOKs, SCGabDeltaDs, 50, 'black', 'filled');
plot(xs, f2.Coefficients.Estimate(1) + xs.*f2.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Gabor: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-5 5]);
ylim([-0.5 0.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-5 0 5]);
set(gca, 'YTick', [-0.8 -0.4 0 0.4 0.8]);
box off;
hold off;

subplot(2,2,4);
hold on;
scatter(SCLumAOKs, SCLumDeltaDs, 50, 'black', 'filled');
plot(xs, f4.Coefficients.Estimate(1) + xs.*f4.Coefficients.Estimate(2), 'Color', 'r', 'LineStyle', '--')
title('SC Luminance: AOK by Delta d''');
xlabel('Area Over the Kernel');
ylabel('delta d''');
axis square;
xlim([-5 5]);
ylim([-0.5 0.5]);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
set(gca, 'XTick', [-5 0 5]);
set(gca, 'YTick', [-0.8 -0.4 0 0.4 0.8]);
box off;
hold off;
 

