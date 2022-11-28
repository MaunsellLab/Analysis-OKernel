% plotOptoEffect

% Loads the tables of session by session data plots scatter of opto vs.
% no-optp d'
V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';
s = ' ';
%% Load Tables and Extract Data

% SC Luminance
scFolder = [analysisDir,'SC',s,'Lum','/'];
load([scFolder 'masterTable.mat'], 'T');
scLum = T;
clear T scFolder;

% Extract mice numbers, d's, c's 
lumAnimals = [scLum.animal];
lumUnique = unique(lumAnimals);
% Extract d' and c
lum_d_TU = [scLum.topUpDPrime];
lum_d_ST = [scLum.stimDPrime];
lum_d_NS = [scLum.noStimDPrime];
lum_c_TU = [scLum.topUpC];
lum_c_ST = [scLum.stimC];
lum_c_NS = [scLum.noStimC];

% SC Gabor
scFolder = [analysisDir,'SC',s,'Gabor','/'];
load([scFolder 'masterTable.mat'], 'T');
scGab = T;
clear T;

% Extract mice numbers, d's, c's 
gabAnimals = [scGab.animal];
gabUnique = unique(gabAnimals);
% Extract d' and c
gab_d_TU = [scGab.topUpDPrime];
gab_d_ST = [scGab.stimDPrime];
gab_d_NS = [scGab.noStimDPrime];
gab_c_TU = [scGab.topUpC];
gab_c_ST = [scGab.stimC];
gab_c_NS = [scGab.noStimC];

%% D-Primes
% Summary Stats of Top Up d'
scLumD_TU = nanmedian(lum_d_TU);
scGabD_TU = nanmedian(gab_d_TU);
% Quantiles
quantile(lum_d_TU, [.25 .75]);
quantile(gab_d_TU, [.25 .75]);

% Summary Stats of No Stim d'
scLumD_NS = nanmedian(lum_d_NS);
scGabD_NS = nanmedian(gab_d_NS);
% Quantiles
quantile(lum_d_NS, [.25 .75]);
quantile(gab_d_NS, [.25 .75]);

% Summary Stats of Stim d'
scLumD_ST = nanmedian(lum_d_ST);
scGabD_ST = nanmedian(gab_d_ST);
% Quantiles
quantile(lum_d_ST, [.25 .75]);
quantile(gab_d_ST, [.25 .75]);

% Delta d'
scLum_DeltaD = nanmedian(lum_d_ST - lum_d_NS);
scGab_DeltaD = nanmedian(gab_d_ST - gab_d_NS);
% Quantiles
quantile(lum_d_ST - lum_d_NS, [.25 .75]);
quantile(gab_d_ST - gab_d_NS, [.25 .75]);
%% Criterions
% Summary Stats of Top Up d'
scLumC_TU = nanmedian(lum_c_TU);
scGabC_TU = nanmedian(gab_c_TU);
% Quantiles
quantile(lum_c_TU, [.25 .75]);
quantile(gab_c_TU, [.25 .75]);

% Summary Stats of No Stim d'
scLumC_NS = nanmedian(lum_c_NS);
scGabC_NS = nanmedian(gab_c_NS);
% Quantiles
quantile(lum_c_NS, [.25 .75]);
quantile(gab_c_NS, [.25 .75]);

% Summary Stats of Stim d'
scLumC_ST = nanmedian(lum_c_ST);
scGabC_ST = nanmedian(gab_c_ST);
% Quantiles
quantile(lum_c_ST, [.25 .75]);
quantile(gab_c_ST, [.25 .75]);

% Delta d'
scLum_DeltaC = nanmedian(lum_c_ST - lum_c_NS);
scGab_DeltaC = nanmedian(gab_c_ST - gab_c_NS);
% Quantiles
quantile(lum_c_ST - lum_c_NS, [.25 .75]);
quantile(gab_c_ST - gab_c_NS, [.25 .75]);


% Scatter Plots of Each 
%% SC Luminance
figure('Position', [100 100 1000 1000]);
subplot(2,2,1);
hold on;
axis square;
scatter(lum_d_NS, lum_d_ST, 20, 'filled', 'MarkerFaceColor', SCColor);
plot([0 4], [0 4], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 4]);
ylim([0 4]);
set(gca,'XTick', [0 2 4]);
set(gca,'YTick', [0 2 4]);
title('SC Luminance d-prime');
xlabel('Unstimulated');
ylabel('Stimulated');
hold off;

% Delta d' Histogram
subplot(2,2,3);
hold on;
histogram(lum_d_ST - lum_d_NS, 'FaceColor',SCColor, 'FaceAlpha', 0.5);
Ys = ylim;
plot([scLum_DeltaD scLum_DeltaD], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
%plot([0 0], [Ys(1) Ys(2)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-1 1]);
ylabel('Counts');
xlabel('Delta d-prime');
title('SC Luminance: Delta d');
hold off;


subplot(2,2,2);
hold on;
axis square;
scatter(lum_c_NS, lum_c_ST, 20, 'filled', 'MarkerFaceColor', SCColor);
plot([0 1.5], [0 1.5], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 1.5]);
ylim([0 1.5]);
set(gca,'XTick', [0 1.5]);
set(gca,'YTick', [0 1.5]);
title('SC Luminance Criterion');
xlabel('Unstimulated');
ylabel('Stimulated');
hold off;

% Delta C Histogram
subplot(2,2,4);
hold on;
histogram(lum_c_ST - lum_c_NS, 'FaceColor',SCColor, 'FaceAlpha', 0.5);
Ys = ylim;
plot([scLum_DeltaC scLum_DeltaC], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
%plot([0 0], [Ys(1) Ys(2)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-0.5 0.5]);
ylabel('Counts');
xlabel('Delta C');
title('SC Luminance: Delta C');
hold off;

%% SC Gabor

figure('Position', [100 100 1000 1000]);
subplot(2,2,1);
hold on;
axis square;
scatter(gab_d_NS, gab_d_ST, 20, 'filled', 'MarkerFaceColor', SCColor);
plot([0 4], [0 4], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 4]);
ylim([0 4]);
set(gca,'XTick', [0 2 4]);
set(gca,'YTick', [0 2 4]);
title('SC Gabor d-prime');
xlabel('Unstimulated');
ylabel('Stimulated');
hold off;

% Delta d' Histogram
subplot(2,2,3);
hold on;
histogram(gab_d_ST - gab_d_NS, 'FaceColor',SCColor, 'FaceAlpha', 0.5);
Ys = ylim;
plot([scGab_DeltaD scGab_DeltaD], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
% plot([0 0], [Ys(1) Ys(2)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-1 1]);
ylabel('Counts');
xlabel('Delta d-prime');
title('SC Gabor: Delta d');
hold off;


subplot(2,2,2);
hold on;
axis square;
scatter(gab_c_NS, gab_c_ST, 20, 'filled', 'MarkerFaceColor', SCColor);
plot([0 1.5], [0 1.5], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 1.5]);
ylim([0 1.5]);
set(gca,'XTick', [0 1.5]);
set(gca,'YTick', [0 1.5]);
title('SC Gabor Criterion');
xlabel('Unstimulated');
ylabel('Stimulated');
hold off;

% Delta C Histogram
subplot(2,2,4);
hold on;
histogram(gab_c_ST - gab_c_NS, 'FaceColor',SCColor, 'FaceAlpha', 0.5);
Ys = ylim;
plot([scGab_DeltaC scGab_DeltaC], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
% plot([0 0], [Ys(1) Ys(2)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-0.5 0.5]);
ylabel('Counts');
xlabel('Delta C');
title('SC Gabor: Delta C');
hold off;


%% V1 Data

v1Folder = [analysisDir,'V1',s,'Lum','/'];
load([v1Folder 'masterTable.mat'], 'T');
v1Lum = T;
clear T v1Folder;

% Extract mice numbers, d's, c's 
lumAnimals = [v1Lum.animal];
lumUnique = unique(lumAnimals);
% Extract d' and c
lum_d_TU = [v1Lum.topUpDPrime];
lum_d_ST = [v1Lum.stimDPrime];
lum_d_NS = [v1Lum.noStimDPrime];
lum_c_TU = [v1Lum.topUpC];
lum_c_ST = [v1Lum.stimC];
lum_c_NS = [v1Lum.noStimC];

% SC Gabor
v1Folder = [analysisDir,'V1',s,'Gabor','/'];
load([v1Folder 'masterTable.mat'], 'T');
v1Gab = T;
clear T;

% Extract mice numbers, d's, c's 
gabAnimals = [v1Gab.animal];
gabUnique = unique(lumAnimals);
% Extract d' and c
gab_d_TU = [v1Gab.topUpDPrime];
gab_d_ST = [v1Gab.stimDPrime];
gab_d_NS = [v1Gab.noStimDPrime];
gab_c_TU = [v1Gab.topUpC];
gab_c_ST = [v1Gab.stimC];
gab_c_NS = [v1Gab.noStimC];


%% D-Primes
% Summary Stats of Top Up d'
v1LumD_TU = nanmedian(lum_d_TU);
v1GabD_TU = nanmedian(gab_d_TU);
% Quantiles
quantile(lum_d_TU, [.25 .75]);
quantile(gab_d_TU, [.25 .75]);

% Summary Stats of No Stim d'
v1LumD_NS = nanmedian(lum_d_NS);
v1GabD_NS = nanmedian(gab_d_NS);
% Quantile
quantile(lum_d_NS, [.25 .75]);
quantile(gab_d_NS, [.25 .75]);

% Summary Stats of Stim d'
v1LumD_ST = nanmedian(lum_d_ST);
v1GabD_ST = nanmedian(gab_d_ST);
% Quantiles
quantile(lum_d_ST, [.25 .75]);
quantile(gab_d_ST, [.25 .75]);

% Delta d'
v1Lum_DeltaD = nanmedian(lum_d_ST - lum_d_NS);
v1Gab_DeltaD = nanmedian(gab_d_ST - gab_d_NS);
% Quantiles
quantile(lum_d_ST - lum_d_NS, [.25 .75]);
quantile(gab_d_ST - gab_d_NS, [.25 .75]);
%% Criterions
% Summary Stats of Top Up d'
v1LumC_TU = nanmedian(lum_c_TU);
v1GabC_TU = nanmedian(gab_c_TU);
% Quantiles
quantile(lum_c_TU, [.25 .75]);
quantile(gab_c_TU, [.25 .75]);

% Summary Stats of No Stim d'
v1LumC_NS = nanmedian(lum_c_NS);
v1GabC_NS = nanmedian(gab_c_NS);
% Quantiles
quantile(lum_c_NS, [.25 .75]);
quantile(gab_c_NS, [.25 .75]);

% Summary Stats of Stim d'
v1LumC_ST = nanmedian(lum_c_ST);
v1GabC_ST = nanmedian(gab_c_ST);
% Quantiles
quantile(lum_c_ST, [.25 .75]);
quantile(gab_c_ST, [.25 .75]);

% Delta d'
v1Lum_DeltaC = nanmedian(lum_c_ST - lum_c_NS);
v1Gab_DeltaC = nanmedian(gab_c_ST - gab_c_NS);
% Quantiles
quantile(lum_c_ST - lum_c_NS, [.25 .75]);
quantile(gab_c_ST - gab_c_NS, [.25 .75]);

%% Plots of V1

% Luminance
figure('Position', [100 100 1000 1000]);
subplot(2,2,1);
hold on;
axis square;
scatter(lum_d_NS, lum_d_ST, 20, 'filled', 'MarkerFaceColor', V1Color);
plot([0 5], [0 5], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 4.75]);
ylim([0 4.75]);
set(gca,'XTick', [0 2 4]);
set(gca,'YTick', [0 2 4]);
title('V1 Luminance d-prime');
xlabel('Unstimulated');
ylabel('Stimulated');
hold off;

% Delta d' Histogram
subplot(2,2,3);
hold on;
histogram(lum_d_ST - lum_d_NS, 'FaceColor',V1Color, 'FaceAlpha', 0.5);
Ys = ylim;
plot([v1Lum_DeltaD v1Lum_DeltaD], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
%plot([0 0], [Ys(1) Ys(2)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-1 1]);
ylabel('Counts');
xlabel('Delta d-prime');
title('V1 Luminance: Delta d');
hold off;


subplot(2,2,2);
hold on;
axis square;
scatter(lum_c_NS, lum_c_ST, 20, 'filled', 'MarkerFaceColor', V1Color);
plot([0 2], [0 2], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 2]);
ylim([0 2]);
set(gca,'XTick', [0 2]);
set(gca,'YTick', [0 2]);
title('V1 Luminance Criterion');
xlabel('Unstimulated');
ylabel('Stimulated');
hold off;

% Delta C Histogram
subplot(2,2,4);
hold on;
histogram(lum_c_ST - lum_c_NS, 'FaceColor',V1Color, 'FaceAlpha', 0.5);
Ys = ylim;
plot([v1Lum_DeltaC v1Lum_DeltaC], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
%plot([0 0], [Ys(1) Ys(2)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-0.5 0.5]);
ylabel('Counts');
xlabel('Delta C');
title('V1 Luminance: Delta C');
hold off;

%% V1 Gabor

figure('Position', [100 100 1000 1000]);
subplot(2,2,1);
hold on;
axis square;
scatter(gab_d_NS, gab_d_ST, 20, 'filled', 'MarkerFaceColor', V1Color);
plot([0 4], [0 4], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 4]);
ylim([0 4]);
set(gca,'XTick', [0 2 4]);
set(gca,'YTick', [0 2 4]);
title('V1 Gabor d-prime');
xlabel('Unstimulated');
ylabel('Stimulated');
hold off;

% Delta d' Histogram
subplot(2,2,3);
hold on;
histogram(gab_d_ST - gab_d_NS, 'FaceColor',V1Color, 'FaceAlpha', 0.5);
Ys = ylim;
plot([v1Gab_DeltaD v1Gab_DeltaD], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
% plot([0 0], [Ys(1) Ys(2)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-1 1]);
ylabel('Counts');
xlabel('Delta d-prime');
title('V1 Gabor: Delta d');
hold off;

subplot(2,2,2);
hold on;
axis square;
scatter(gab_c_NS, gab_c_ST, 20, 'filled', 'MarkerFaceColor', V1Color);
plot([0 1.5], [0 1.5], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 1.5]);
ylim([0 1.5]);
set(gca,'XTick', [0 1.5]);
set(gca,'YTick', [0 1.5]);
title('V1 Gabor Criterion');
xlabel('Unstimulated');
ylabel('Stimulated');
hold off;

% Delta C Histogram
subplot(2,2,4);
hold on;
histogram(gab_c_ST - gab_c_NS, 'FaceColor',V1Color, 'FaceAlpha', 0.5);
Ys = ylim;
plot([v1Gab_DeltaC v1Gab_DeltaC], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
% plot([0 0], [Ys(1) Ys(2)], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '-');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-0.5 0.5]);
ylabel('Counts');
xlabel('Delta C');
title('V1 Gabor: Delta C');
hold off;

%% Repeat For False Alarm Rates

% Combine Across Stimuli
scFA_Stim = [scLum.stimPFA; scGab.stimPFA];
scFA_NoStim = [scLum.noStimPFA; scGab.noStimPFA];
v1FA_Stim = [v1Lum.stimPFA; v1Gab.stimPFA];
v1FA_NoStim = [v1Lum.noStimPFA; v1Gab.noStimPFA];

% SC delta FA
nanmedian(scFA_Stim - scFA_NoStim);
quantile(scFA_Stim - scFA_NoStim, [.25 .75])
% SC no Stim FA
nanmedian(scFA_NoStim)
quantile(scFA_NoStim, [.25 .75])
% SC Stim FA
nanmedian(scFA_Stim)
quantile(scFA_Stim, [.25 .75])


figure('Position', [100 100 1000 1000]);
subplot(2,2,1);
hold on;
axis square;
scatter(scFA_NoStim(1:length(scLum.noStimPFA)), scFA_Stim(1:length(scLum.stimPFA)), 25, 'x', 'MarkerEdgeColor', SCColor);
scatter(scFA_NoStim(length(scLum.noStimPFA)+1:end), scFA_Stim(length(scLum.noStimPFA)+1:end), 40, 'o', 'MarkerEdgeColor', SCColor);
plot([0 0.2], [0 0.2], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 0.2]);
ylim([0 0.2]);
set(gca,'XTick', [0 0.2]);
set(gca,'YTick', [0 0.2]);
title('SC False Alarm');
xlabel('Unstimulated');
ylabel('Stimulated');
legend('Luminance', 'Gabor', 'Location', 'northwest');
hold off;

% Delta FA
subplot(2,2,3)
hold on;
histogram(scFA_Stim - scFA_NoStim, 'FaceColor', SCColor, 'FaceAlpha', 0.5);
Ys = ylim;
plot([nanmedian(scFA_Stim - scFA_NoStim) nanmedian(scFA_Stim - scFA_NoStim)], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-0.1 0.1]);
ylabel('Counts');
xlabel('Delta C');
title('SC: delta FA');
hold off;

% V1 delta FA
nanmedian(v1FA_Stim - v1FA_NoStim)
quantile(v1FA_Stim - v1FA_NoStim, [.25 .75])
% SC no Stim FA
nanmedian(v1FA_NoStim)
quantile(v1FA_NoStim, [.25 .75])
% SC Stim FA
nanmedian(v1FA_Stim)
quantile(v1FA_Stim, [.25 .75])

subplot(2,2,2);
hold on;
axis square;
scatter(v1FA_NoStim(1:length(v1Lum.noStimPFA)), v1FA_Stim(1:length(v1Lum.stimPFA)), 25, 'x', 'MarkerEdgeColor', V1Color);
scatter(v1FA_NoStim(length(v1Lum.noStimPFA)+1:end), v1FA_Stim(length(v1Lum.noStimPFA)+1:end), 40, 'o', 'MarkerEdgeColor', V1Color);
plot([0 0.2], [0 0.2], 'LineStyle','--', 'LineWidth', 1, 'Color', 'k');
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([0 0.15]);
ylim([0 0.15]);
set(gca,'XTick', [0 0.15]);
set(gca,'YTick', [0 0.15]);
title('V1 False Alarm');
xlabel('Unstimulated');
ylabel('Stimulated');
legend('Luminance', 'Gabor', 'Location', 'northwest');
hold off;

% Delta FA
subplot(2,2,4)
hold on;
histogram(v1FA_Stim - v1FA_NoStim, 'FaceColor', V1Color, 'FaceAlpha', 0.5);
Ys = ylim;
plot([nanmedian(v1FA_Stim - v1FA_NoStim) nanmedian(v1FA_Stim - v1FA_NoStim)], [Ys(1) Ys(2)], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-'); % Mark Median
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
xlim([-0.1 0.1]);
ylabel('Counts');
xlabel('Delta C');
title('V1: delta FA');
hold off;



