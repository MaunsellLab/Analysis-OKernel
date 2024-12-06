% Use After running areaOverKernel_GrandAv_iter which bootstraps across
% each time bin to compute the AOK distribution. 
% This script checks which bins were significant to plot 
% line of significant AOKs

% In the original version of this I didn't test 100 ms at the beginning 
% and end so we can set those to nans 

% First 100 untested bins
% V1 Data
v1p_lum(1,1:25) = nan;
v1p_gab(1,1:25) = nan;
% SC Data
scp_lum(1,1:25) = nan;
scp_gab(1,1:25) = nan;

% last 100 bins
% V1 Data
v1p_lum(1,776:end) = nan;
v1p_gab(1,776:end) = nan;
% SC Data
scp_lum(1,776:end) = nan;
scp_gab(1,776:end) = nan;



% Convert Vectors to T/F
% Convert to Logical
v1Lum_sig = v1p_lum < 0.001;
v1Gab_sig = v1p_gab < 0.001;

% Convert to Logical
scLum_sig = scp_lum < 0.001;
scGab_sig = scp_gab < 0.001;

%% Histograms

% Time Series
bins = 1:1:800;

% Figure
figure('Units', 'inches', 'Position', [3, 1, 8, 10]);

subplot(2,2,1);
hold on;
plot([400 400], [0.5 1.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth',0.5);
scatter(bins,v1Gab_sig, 20, "red", "filled");
title('V1 Gabor: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
xlim([1 800]);
ylim([0.5 1.5]);
hold off;

subplot(2,2,2);
hold on;
plot([400 400], [0.5 1.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth',0.5);
scatter(bins,v1Lum_sig, 20, "red", "filled");
title('V1 Luminance: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
xlim([1 800]);
ylim([0.5 1.5]);
hold off;

subplot(2,2,3);
hold on;
plot([400 400], [0.5 1.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth',0.5);
scatter(bins,scGab_sig, 20, "red", "filled");
title('SC Gabor: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
xlim([1 800]);
ylim([0.5 1.5]);
hold off;

subplot(2,2,4);
hold on;
plot([400 400], [0.5 1.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth',0.5);
scatter(bins,scLum_sig, 20, "red", "filled");
title('SC Gabor: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
xlim([1 800]);
ylim([0.5 1.5]);
hold off;










