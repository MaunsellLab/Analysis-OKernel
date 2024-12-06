% Use After running areaOverKernel_GrandAv_OFFSET_iter which bootstraps across
% each time bin to compute the AOK distribution. 
% This script checks which bins were significant to plot 
% line of significant AOKs

% Because it's a rolling 25 ms window, I didn't test 25 ms at the beginning 
% and end so we can set those to nans 

% First 25 untested bins
% V1 Data
v1p_off(1,1:25) = nan;
% SC Data
scp_off(1,1:25) = nan;

% last 25 bins
% V1 Data
v1p_off(1,776:800) = nan;
% SC Data
scp_off(1,776:800) = nan;

%% Convert Vectors to T/F
% Convert to Logical
v1_sig = v1p_off < 0.001;
sc_sig = scp_off < 0.001;

%% Histograms

% Time Series
bins = 1:1:800;

% Figure
figure('Units', 'inches', 'Position', [3, 1, 8, 5]);

subplot(1,2,1);
hold on;
plot([400 400], [0.5 1.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth',0.5);
scatter(bins,v1_sig, 20, "red", "filled");
title('V1 Gabor: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
xlim([1 800]);
ylim([0.5 1.5]);
hold off;

subplot(1,2,2);
hold on;
plot([400 400], [0.5 1.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth',0.5);
scatter(bins,sc_sig, 20, "red", "filled");
title('SC Gabor: AOK Bootstraps');
box off;
set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);
xlim([1 800]);
ylim([0.5 1.5]);
hold off;






