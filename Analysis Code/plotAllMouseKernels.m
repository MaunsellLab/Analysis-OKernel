% Plots individual animal kernels vertically
%% Set Up
dp_cut = 1; % paper cutoff for delta d'
nSigma = 1;
timeBins = 1:1:800;
limitY = [-0.1 0.1];
plotH = 5000;
plotW = 1000;
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

% Subselect based on delta d'
A = selectUsingLimits(lumT, limits);
% Fitted Delta D' Distribution for V1 Luminance
if dp_cut == 1
    mu =  -0.011;
    sigma =  0.2326;
    deltaD = [A.stimDPrime]-[A.noStimDPrime];
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
else
     keepIdx = ones(height(A),1);
end
% Drop Sessions that don't meet criteria
lumT = A(keepIdx,:);

% Number of mice to set spacing
nMice = length(animals);
V1Lum = cell(nMice,2);

for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    % subselect for this mouse
    U = selectUsingLimits(lumT, limits);
    D = nanmean([U.stimDPrime]-[U.noStimDPrime]);
    V1Lum{a,2} = D;
    
    lumKernel = []; % Init
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        lumKernel = [lumKernel; [stimProfiles.hitProfiles; -stimProfiles.missProfiles]];
    end
    V1Lum{a,1} = lumKernel; % Write to cell array.
end

% Rank By Delta d'
[values, ranks] = sort([V1Lum{:,2}],'ascend');

figure('Position',[10 10 plotW plotH]);
for a=1:nMice
    % compute SEM for this mouse = sqrt(pq/n)
    SEM = 2*sqrt(.5*(1-.5)/size(V1Lum{ranks(a),1},1));
    subplot(6,2,a);
    hold on;
    ylim(limitY); % The range for each kernel is [-0.1 - 0.1]
    plot(timeBins, mean(V1Lum{ranks(a),1},1), 'LineWidth',2, 'Color', 'k');
    % Plot SEM as shaded region
    X = [1 800 800 1]; Y = [-SEM -SEM SEM SEM];
    fill(X,Y,[0.9 0.9 0.9],'EdgeColor','k', 'FaceAlpha', 0.25);
    % Line For StimOnset
    plot([400 400], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k');
    % Plot Minor Tick Dashed Lines
    plot([100 100], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([200 200], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([300 300], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([500 500], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([600 600], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([700 700], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    % plot Horizontal Line at 0.0
    plot([1 800], [0.0 0.0], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k');
    titleString = [animals{ranks(a)}, '; ', 'delta d:''',  ' ' , num2str(round(V1Lum{ranks(a),2},3)), '; ', num2str(size(V1Lum{ranks(a),1},1)), ' trials'];
    title(titleString);
    set(gca, 'XTick', [1 401 800]);
    set(gca, 'XTickLabel', {'-400', '0', '400'});
    set(gca, 'TickDir', 'out');
end
hold off;
%% Get StimProfiles and Compute AOK For V1 Gabor
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
animals = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
condition = ['V1' ' ' 'Gabor'];

% Subselect based on delta d'
A = selectUsingLimits(gabT, limits);

% Fitted Delta D' Distribution for V1 Gabor
if dp_cut == 1
    mu =  -0.0034;
    sigma = 0.277;
    deltaD = [A.stimDPrime]-[A.noStimDPrime];
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
else
     keepIdx = ones(height(A),1);
end
% Drop Sessions that don't meet criteria
gabT = A(keepIdx,:);

% Number of mice to set spacing
nMice = length(animals);
V1Gab = cell(nMice,2);



for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(gabT, limits);
    D = nanmean([U.stimDPrime]-[U.noStimDPrime]);
    V1Gab{a,2} = D;
    gabKernel = []; % Init
    
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        gabKernel = [gabKernel; [stimProfiles.hitProfiles; -stimProfiles.missProfiles]];
    end
    V1Gab{a,1} = gabKernel;
end

% Rank By Delta d'
[values, ranks] = sort([V1Gab{:,2}],'ascend');


figure('Position',[10 10 plotW plotH]);
for a = 1:nMice
    % compute SEM for this mouse = sqrt(pq/n)
    SEM = 2*sqrt(.5*(1-.5)/size(V1Gab{ranks(a),1},1));
    subplot(6,2,a);
    hold on;
    ylim(limitY); % The range for each kernel is [-0.1 - 0.1] 
    plot(timeBins, mean(V1Gab{ranks(a),1},1), 'LineWidth',2, 'Color', 'k');
    % Plot SEM as shaded region
    X = [1 800 800 1]; Y = [-SEM -SEM SEM SEM];
    fill(X,Y,[0.9 0.9 0.9],'EdgeColor','k', 'FaceAlpha', 0.25);
    % Line For StimOnset
    plot([400 400], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k');
    % Plot Minor Tick Dashed Lines
    plot([100 100], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([200 200], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([300 300], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([500 500], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([600 600], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([700 700], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    % plot Horizontal Line at 0.0
    plot([1 800], [0.0 0.0], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k');
    titleString = [animals{ranks(a)}, '; ', 'delta d:''',  ' ' , num2str(round(V1Gab{ranks(a),2},3)), '; ', num2str(size(V1Gab{ranks(a),1},1)), ' trials'];
    title(titleString);
    set(gca, 'XTick', [1 401 800]);
    set(gca, 'XTickLabel', {'-400', '0', '400'});
    set(gca, 'TickDir', 'out');
end
hold off;
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
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
animals = {'1458', '1548', '1674', '1675', '1902', '1905', '2057', '2058', '2063', '2169', '2236'};
condition = ['SC' ' ' 'Lum'];

A = selectUsingLimits(lumT, limits);

% Fitted Delta D' Distribution for SC Luminance
if dp_cut == 1
    mu = -0.03914;
    sigma = 0.2715;
    deltaD = [A.stimDPrime]-[A.noStimDPrime];
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
else
    keepIdx = ones(height(A),1);
end

% Drop Sessions that don't meet criteria
lumT = A(keepIdx,:);

% Number of mice to set spacing
nMice = length(animals);
SCLum = cell(nMice,2);



for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(lumT, limits);
    D = nanmean([U.stimDPrime]-[U.noStimDPrime]);
    SCLum{a,2} = D;
    lumKernel = []; % Init
    
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        lumKernel = [lumKernel; [stimProfiles.hitProfiles; -stimProfiles.missProfiles]];
    end
    SCLum{a,1} = lumKernel;
end

% Rank By Delta d'
[values, ranks] = sort([SCLum{:,2}],'ascend');

figure('Position',[10 10 plotW plotH]);
for a = 1:nMice
    % compute SEM for this mouse = sqrt(pq/n)
    SEM = 2*sqrt(.5*(1-.5)/size(SCLum{ranks(a),1},1));
    subplot(6,2,a);
    hold on;
    ylim(limitY); % The range for each kernel is [-0.1 - 0.1] 
    plot(timeBins, mean(SCLum{ranks(a),1},1), 'LineWidth',2, 'Color', 'k');
    % Plot SEM as shaded region
    X = [1 800 800 1]; Y = [-SEM -SEM SEM SEM];
    fill(X,Y,[0.9 0.9 0.9],'EdgeColor','k', 'FaceAlpha', 0.25);
    % Line For StimOnset
    plot([400 400], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k');
    % Plot Minor Tick Dashed Lines
    plot([100 100], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([200 200], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([300 300], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([500 500], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([600 600], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([700 700], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    % plot Horizontal Line at 0.0
    plot([1 800], [0.0 0.0], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k');
    titleString = [animals{ranks(a)}, '; ', 'delta d:''',  ' ' , num2str(round(SCLum{ranks(a),2},3)), '; ', num2str(size(SCLum{ranks(a),1},1)), ' trials'];
    title(titleString);
    set(gca, 'XTick', [1 401 800]);
    set(gca, 'XTickLabel', {'-400', '0', '400'});
    set(gca, 'TickDir', 'out');
end
hold off;
%% SC Gabor
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
animals = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
condition = ['SC' ' ' 'Gabor'];

A = selectUsingLimits(gabT, limits);

% Fitted Delta D' Distribution for SC Gabor
if dp_cut == 1
    mu = 0.0056;
    sigma = 0.287;
    deltaD = [A.stimDPrime]-[A.noStimDPrime];
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
else
    keepIdx = ones(height(A),1);
end

% Clean Up Table
gabT = A(keepIdx,:);

% Number of mice to set spacing
nMice = length(animals);
SCGab = cell(nMice,2);

figure('Position',[10 10 plotW plotH]);
title('SC Gabor');
for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(gabT, limits);
    D = nanmean([U.stimDPrime]-[U.noStimDPrime]);
    SCGab{a,2} = D;
    gabKernel = []; % Init
    
    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
        gabKernel = [gabKernel; [stimProfiles.hitProfiles; -stimProfiles.missProfiles]];
    end
    SCGab{a,1} = gabKernel;
end

% Rank By Delta d'
[values, ranks] = sort([SCGab{:,2}],'ascend');

for a = 1:nMice
    % compute SEM for this mouse = sqrt(pq/n)
    SEM = 2*sqrt(.5*(1-.5)/size(SCGab{ranks(a),1},1));
    subplot(6,2,a);
    hold on;
    ylim(limitY); % The range for each kernel is [-0.1 - 0.1] 
    plot(timeBins, mean(SCGab{ranks(a),1},1), 'LineWidth',2, 'Color', 'k');
    % Plot SEM as shaded region
    X = [1 800 800 1]; Y = [-SEM -SEM SEM SEM];
    fill(X,Y,[0.9 0.9 0.9],'EdgeColor','k', 'FaceAlpha', 0.25);
    % Line For StimOnset
    plot([400 400], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k');
    % Plot Minor Tick Dashed Lines
    plot([100 100], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([200 200], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([300 300], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([500 500], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([600 600], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    plot([700 700], [-0.1 0.1], 'LineWidth', 1, 'LineStyle', '--', 'Color', '[0.5 0.5 0.5]');
    % plot Horizontal Line at 0.0
    plot([1 800], [0.0 0.0], 'LineWidth', 1, 'LineStyle', '--', 'Color', 'k');
    titleString = [animals{ranks(a)}, '; ', 'delta d:''',  ' ' , num2str(round(SCGab{ranks(a),2},3)), '; ', num2str(size(SCGab{ranks(a),1},1)), ' trials'];
    title(titleString);
    set(gca, 'XTick', [1 401 800]);
    set(gca, 'XTickLabel', {'-400', '0', '400'});
    set(gca, 'TickDir', 'out');
end
hold off;