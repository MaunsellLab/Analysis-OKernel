function plotKernelandRT(stimulus, kernelType, dp_cut, nSigma, dp_cut_at0)
%stimulus = 'Lum';
%stimulus = 'Gabor';
%stimulus = 'Offset';
%stimulus = 'FA';

% kernelType = 'stimulus'
% kernelType = 'response'

V1Color = [0.8500 0.3250 0.0980];
SCColor = [0.3010 0.7450 0.9330];

% Function takes a stimulus condition and plots the kernels for comparisons
% across V1 and SC
% 'Lum'
% 'Gabor'
% 'Offset'
% 'FA'

% dp_cut
% if true, it drops sessions that have delta d-prime > mu +nSigma*sigma

% dp_cut_at0
% if true, drops any session with a delta d' > 0.
% Overrides dp_cut
if dp_cut_at0 == 1
    dp_cut = 0;
end
%% Set Up
% Directory With all the final Tables and Profiles
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';
% Set Limits
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;
limits.yAxis = 0.5;
limits.yAxis = 0.0;
limits.numBoot = 250;

% Stuff For Plots
[plotStartMS, plotEndMS, plotRTStartMS] = plotLimits();
yLabel = 'Normalized Power';
%% Load Tables and Profiles for each Brain Area
%% SC
s = ' ';
scFolder = [analysisDir,'SC',s,stimulus,'/'];
load([scFolder 'masterTable.mat'], 'T');
scT = T;

U = selectUsingLimits(scT, limits);
condition = ['SC' ' ' stimulus];

% Sub-Select Sessions if you want to do this based on delta dPrime.
if dp_cut == 1
 % Previously made fits of histcounts to find mean and sigma
    if strcmp(stimulus, 'Lum')
        mu = -0.03914;
        sigma = 0.2715;
    elseif strcmp(stimulus, 'Gabor')
        mu = 0.0056;
        sigma = 0.287;
    elseif strcmp(stimulus, 'Offset') % no d' on filter controls (offset)
        mu = 0;
        sigma = 1000;
    elseif strcmp(stimulus, 'FA')
        mu = -0.023;
        sigma = 0.277;
    end
    x =  -1:0.01:1;
    % Compute delta d' cutoff;
    deltaD = [U.stimDPrime]-[U.noStimDPrime];
    bins = -1:0.1:1;
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
    nDropped = height(U) - sum(keepIdx);

elseif dp_cut_at0 == 1 % Drop sessions with delta d' > 0
    deltaD = [U.stimDPrime]-[U.noStimDPrime];
    keepIdx = deltaD < 0;
    nDropped = height(U) - sum(keepIdx);

else % Keep all sessions
    keepIdx = ones(height(U),1);
end

% Clean Up Table
U = U(keepIdx,:);

% Get all RTs
scHitRTs = [];
for sessionNum = 1:height(scT)
    if keepIdx(sessionNum) == 1
        scHitRTs = [scHitRTs, cell2mat(scT.stimCorrectRTs(sessionNum))];
    else 
        continue
    end
end

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

%% Repeat For V1
v1Folder = [analysisDir,'V1',s,stimulus,'/'];
load([v1Folder 'masterTable.mat'], 'T');
v1T = T;
clear T;

U = selectUsingLimits(v1T, limits);
condition = ['V1' ' ' stimulus];

% Sub-Select Sessions if you want to do this based on delta dPrime.
if dp_cut == 1
    % Previously made fits of histcounts to find mean and sigma
    if strcmp(stimulus, 'Lum')
        mu =  -0.011;
        sigma =  0.2326;
    elseif strcmp(stimulus, 'Gabor')
        mu =  -0.0034;
        sigma = 0.277;
    elseif strcmp(stimulus, 'Offset') % no d' on filter controls (offset)
        mu = 0;
        sigma = 1000;
    elseif strcmp(stimulus, 'FA')
        mu = -0.06;
        sigma = 0.258;
    end
    x =  -1:0.01:1;
    % Compute delta d' cutoff;
    deltaD = [U.stimDPrime]-[U.noStimDPrime];
    bins = -1:0.1:1;
    cutoff = mu + nSigma*sigma;
    keepIdx = deltaD <= cutoff;
    nDropped = height(U) - sum(keepIdx);

elseif dp_cut_at0 == 1 % Drop sessions with delta d' > 0
    deltaD = [U.stimDPrime]-[U.noStimDPrime];
    keepIdx = deltaD < 0;
    nDropped = height(U) - sum(keepIdx);

else % Keep all sessions
    keepIdx = ones(height(U),1);
end

% Clean Up Table
U = U(keepIdx,:);

% Get all RTs
v1HitRTs = [];
for sessionNum = 1:height(v1T)
    if keepIdx(sessionNum) == 1;
        v1HitRTs = [v1HitRTs, cell2mat(v1T.stimCorrectRTs(sessionNum))];
    else
        continue
    end
end


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


%% Stats of Performance
nHitsSC  = size(scProfiles.hitProfiles, 1);
nMissSC  = size(scProfiles.missProfiles, 1);
nEarlySC = size(scProfiles.earlyProfiles, 1);

nHitsV1  = size(v1Profiles.hitProfiles, 1);
nMissV1  = size(v1Profiles.missProfiles, 1);
nEarlyV1 = size(v1Profiles.earlyProfiles, 1);
%% Collate All Traces

% SC Kernel
scKernel = [scProfiles.hitProfiles; -scProfiles.missProfiles];
scKernel_hit = scProfiles.hitProfiles / 2 + 0.5;
scKernel_miss = scProfiles.missProfiles / 2 + 0.5;
scKernel_RT = scProfiles.RTProfiles / 2 + 0.5;
scKernel_early = scProfiles.earlyProfiles / 2 + 0.5;

% V1 Kernels
v1Kernel = [v1Profiles.hitProfiles; -v1Profiles.missProfiles];
v1Kernel_hit = v1Profiles.hitProfiles / 2 + 0.5;
v1Kernel_miss = v1Profiles.missProfiles / 2 + 0.5;
v1Kernel_RT = v1Profiles.RTProfiles / 2 + 0.5;
v1Kernel_early = v1Profiles.earlyProfiles / 2 + 0.5;

%% Bootstrap and Filter

% Design Filter For Kernels
sampleFreqHz = 1000;
filterLP = designfilt('lowpassfir', 'PassbandFrequency', 90 / sampleFreqHz, ...
    'StopbandFrequency', 2 * 90 / sampleFreqHz, 'PassbandRipple', 1, 'StopbandAttenuation', 60, ...
    'DesignMethod','equiripple');

% Different Plots depending on whether you are looking at
% stimulus/motor kernels

if strcmp(kernelType, 'stimulus') % Stimulus Aligned Combined Kernel
    % Bootstrap V1
    v1BootKernel = bootstrp(250, @mean, v1Kernel);
    v1PCs = prctile(v1BootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
    v1PCMeans = mean(v1PCs, 2);
    CIs = zeros(3, size(v1BootKernel, 2));
    for c = 1:3
        v1CIs(c, :) = filtfilt(filterLP, v1PCs(c, :) - v1PCMeans(c)) + v1PCMeans(c);
    end
    v1x = 1:size(v1CIs, 2);

    % Bootstrap SC
    scBootKernel = bootstrp(250, @mean, scKernel);
    scPCs = prctile(scBootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
    scPCMeans = mean(scPCs, 2);
    CIs = zeros(3, size(scBootKernel, 2));
    for c = 1:3
        scCIs(c, :) = filtfilt(filterLP, scPCs(c, :) - scPCMeans(c)) + scPCMeans(c);
    end
    scx = 1:size(scCIs, 2);

    %% Plot
    bins = size(v1BootKernel, 2);
    figure('Units', 'inches', 'Position', [3, 1, 8, 10]);
    % V1
    subplot(2,1,1);
    yyaxis left
    plot([0, bins], [limits.yAxis, limits.yAxis], 'Color', 'k', 'LineStyle','--');
    hold on;
    plot(v1x, v1CIs(2, :), 'Color', V1Color, 'LineStyle', '-', 'LineWidth',2);
    v1x2 = [v1x, fliplr(v1x)];
    v1fillCI = [v1CIs(1, :), fliplr(v1CIs(3, :))];
    fill(v1x2, v1fillCI, V1Color, 'lineStyle', '-', 'edgeColor', V1Color, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
    ax = gca;
    ax.YColor = [0 0 0];
    ax.TickDir = 'out';
    ax.XGrid = 'on';
    ax.FontSize = 14;
    bins = size(v1Kernel, 2);
    xlim(ax, [0, bins]);
    ylabel(yLabel);
    %llim1 = ylim;

    set(gca,'XTick', [0, -plotStartMS, bins]);
    set(gca, 'XTickLabel', {sprintf('%d', plotStartMS), '0', sprintf('%d', plotEndMS)});
    set(gca, 'LineWidth', 1);
    ax.XMinorGrid = 'on';
    ax.MinorGridLineStyle = '--';
    ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
    xlabel('Time Relative to Stimulus');
    
    
    % Plot RTs on right y-axes
    edges = [0:25:800];
    yyaxis right
    histogram(v1HitRTs+400,edges, 'Normalization', 'probability', 'FaceAlpha', 0.1, 'FaceColor', V1Color, 'edgeAlpha', 0.5); % Have to shift by +400 to align with zero on the x-axis
    ylabel('probability of behavioral response');
    plotTitle = sprintf('V1 Kernel (n=%d)', nHitsV1+nMissV1);
    title(plotTitle);
    hold off;

    % SC
    subplot(2,1,2);
    bins = size(scBootKernel, 2);
    plotTitle = sprintf('SC Hit Kernel (n=%d)', nHitsSC);
    yyaxis left
    ax.YColor = [0 0 0];
    plot([0, bins], [limits.yAxis, limits.yAxis],  'Color', 'k', 'LineStyle','--');
    hold on;
    plot(scx, scCIs(2, :), 'Color', SCColor, 'LineStyle', '-', 'LineWidth', 2);
    ax.YColor = [0 0 0];
    scx2 = [scx, fliplr(scx)];
    scfillCI = [scCIs(1, :), fliplr(scCIs(3, :))];
    fill(scx2, scfillCI, SCColor, 'lineStyle', '-', 'edgeColor', SCColor, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
    ax = gca;
    ax.YColor = 'k';
    ax.TickDir = 'out';
    ax.FontSize = 14;
    bins = size(scKernel, 2);
    xlim(ax, [0, bins]);
    %llim2 = ylim;
    ax.XGrid = 'on';
    ylabel(yLabel);
    title(plotTitle);


    set(gca,'XTick', [0, -plotStartMS, bins]);
    set(gca, 'XTickLabel', {sprintf('%d', plotStartMS), '0', sprintf('%d', plotEndMS)});
    set(gca, 'LineWidth', 1);
    ax.XMinorGrid = 'on';
    ax.MinorGridLineStyle = '--';
    ax.XAxis.MinorTickValues = [0 100 200 300 500 600 700 800];
    xlabel('Time Relative to Stimulus');
    % Plot RTs on right y-axes
    edges = [0:25:800];
    yyaxis right
    histogram(scHitRTs+400,edges, 'Normalization', 'probability', 'FaceAlpha', 0.1, 'FaceColor', SCColor, 'edgeAlpha', 0.5); % Have to shift by +400 to align with zero on the x-axis
    ylabel('probability of behavioral response');
    plotTitle = sprintf('SC Kernel (n=%d)', nHitsSC+nMissSC);
    title(plotTitle);
    hold off;

    % Put plots on same axes
    %lims = [min(min(llim1, llim2)) max(max(llim1, llim2))];
    lims = [-0.04 0.02];
    subplot(2,1,1)
    yyaxis left
    ylim(lims);
    yyaxis right
    ax.YColor = 'k';
    ylim([0 1]);
    set(gca,'YTick', [0, 0.25]);

    subplot(2,1,2);
    yyaxis left
    ax.YColor = 'k';
    ylim(lims);
    yyaxis right
    ax.YColor = 'k';
    ylim([0 1]);
    set(gca,'YTick', [0, 0.25]);

elseif strcmp(kernelType, 'response') % FA and RT Kernels
   
    % False Alarms
    v1BootKernel = bootstrp(250, @mean, v1Kernel_early);
    v1PCs = prctile(v1BootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
    v1PCMeans = mean(v1PCs, 2);
    CIs = zeros(3, size(v1BootKernel, 2));
    for c = 1:3
        v1CIs(c, :) = filtfilt(filterLP, v1PCs(c, :) - v1PCMeans(c)) + v1PCMeans(c);
    end
    v1x = 1:size(v1CIs, 2);

    % Bootstrap SC
    scBootKernel = bootstrp(250, @mean, scKernel_early);
    scPCs = prctile(scBootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
    scPCMeans = mean(scPCs, 2);
    CIs = zeros(3, size(scBootKernel, 2));
    for c = 1:3
        scCIs(c, :) = filtfilt(filterLP, scPCs(c, :) - scPCMeans(c)) + scPCMeans(c);
    end
    scx = 1:size(scCIs, 2);

    %% Plot
    figure('Units', 'inches', 'Position', [3, 1, 8, 10]);
    % V1
    bins = size(v1BootKernel, 2);
    subplot(2,1,1);
    plot([0, bins], [0.5, 0.5], 'Color', 'k', 'LineStyle','--');
    hold on;
    plot(v1x, v1CIs(2, :), 'Color', V1Color, 'LineStyle', '-', 'LineWidth', 2);
    v1x2 = [v1x, fliplr(v1x)];
    v1fillCI = [v1CIs(1, :), fliplr(v1CIs(3, :))];
    fill(v1x2, v1fillCI, V1Color, 'lineStyle', '-', 'edgeColor', V1Color, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
    ax = gca;
    ax.YColor = [0 0 0];
    ax.TickDir = 'out';
    ax.XGrid = 'on';
    ax.FontSize = 14;
    bins = size(v1Kernel_early, 2);
    xlim(ax, [0, bins]);
    ylabel(yLabel);
    %llim1 = ylim;
    set(gca,'XTick', [0, -plotStartMS+200, bins]);
    set(gca, 'XTickLabel', {sprintf('%d', plotRTStartMS), '0', sprintf('%d', plotRTStartMS + plotEndMS - plotStartMS)});
    set(gca, 'LineWidth', 1);
    ax.XMinorGrid = 'on';
    ax.MinorGridLineStyle = '--';
    ax.XAxis.MinorTickValues = [0 100 200 300 400 500 700 800];
    xlabel('Time Relative to Lever Release');
    plotTitle = sprintf('V1 FA Kernel (n=%d)', nEarlyV1);
    title(plotTitle);
    hold off;

    % SC
    subplot(2,1,2);
    bins = size(scBootKernel, 2);
    ax.YColor = [0 0 0];
    plot([0, bins], [0.5, 0.5],  'Color', 'k', 'LineStyle','--');
    hold on;
    plot(scx, scCIs(2, :), 'Color', SCColor, 'LineStyle', '-', 'LineWidth', 2);
    ax.YColor = [0 0 0];
    scx2 = [scx, fliplr(scx)];
    scfillCI = [scCIs(1, :), fliplr(scCIs(3, :))];
    fill(scx2, scfillCI, SCColor, 'lineStyle', '-', 'edgeColor', SCColor, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
    ax = gca;
    ax.YColor = [0 0 0];
    ax.TickDir = 'out';
    ax.FontSize = 14;
    bins = size(scKernel, 2);
    xlim(ax, [0, bins]);
    %llim2 = ylim;
    ax.XGrid = 'on';
    ylabel(yLabel);
    plotTitle = sprintf('SC FA Kernel (n=%d)', nEarlySC);
    title(plotTitle);
    set(gca,'XTick', [0, -plotStartMS+200, bins]);
    set(gca, 'XTickLabel', {sprintf('%d', plotRTStartMS), '0', sprintf('%d', plotRTStartMS + plotEndMS - plotStartMS)});
    set(gca, 'LineWidth', 1);
    ax.XMinorGrid = 'on';
    ax.MinorGridLineStyle = '--';
    ax.XAxis.MinorTickValues = [0 100 200 300 400 500 700 800];
    xlabel('Time Relative to Lever Release');
    hold off;
    % Put plots on same axes
    %lims = [min(min(llim1, llim2)) max(max(llim1, llim2))];
    lims = [0.47 0.53];
    subplot(2,1,1);
    ylim(lims);
    subplot(2,1,2);
    ylim(lims);  


    %% Repeat For Reaction Times
    v1BootKernel = bootstrp(250, @mean, v1Kernel_RT);
    v1PCs = prctile(v1BootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
    v1PCMeans = mean(v1PCs, 2);
    CIs = zeros(3, size(v1BootKernel, 2));
    for c = 1:3
        v1CIs(c, :) = filtfilt(filterLP, v1PCs(c, :) - v1PCMeans(c)) + v1PCMeans(c);
    end
    v1x = 1:size(v1CIs, 2);

    % Bootstrap SC
    scBootKernel = bootstrp(250, @mean, scKernel_RT);
    scPCs = prctile(scBootKernel, [15.9, 50, 84.1]);              % +/- 1 SEM
    scPCMeans = mean(scPCs, 2);
    CIs = zeros(3, size(scBootKernel, 2));
    for c = 1:3
        scCIs(c, :) = filtfilt(filterLP, scPCs(c, :) - scPCMeans(c)) + scPCMeans(c);
    end
    scx = 1:size(scCIs, 2);

    % Plot
    figure('Units', 'inches', 'Position', [3, 1, 8, 10]);
    % V1
    bins = size(v1BootKernel, 2);
    subplot(2,1,1);
    plot([0, bins], [0.5, 0.5], 'Color', 'k', 'LineStyle','--');
    hold on;
    plot(v1x, v1CIs(2, :), 'Color', V1Color, 'LineStyle', '-', 'LineWidth', 2);
    v1x2 = [v1x, fliplr(v1x)];
    v1fillCI = [v1CIs(1, :), fliplr(v1CIs(3, :))];
    fill(v1x2, v1fillCI, V1Color, 'lineStyle', '-', 'edgeColor', V1Color, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
    ax = gca;
    ax.YColor = [0 0 0];
    ax.TickDir = 'out';
    ax.XGrid = 'on';
    ax.FontSize = 14;
    bins = size(v1Kernel_RT, 2);
    xlim(ax, [0, bins]);
    ylabel(yLabel);
    %llim1 = ylim;
    set(gca,'XTick', [0, -plotStartMS+200, bins]);
    set(gca, 'XTickLabel', {sprintf('%d', plotRTStartMS), '0', sprintf('%d', plotRTStartMS + plotEndMS - plotStartMS)});
    set(gca, 'LineWidth', 1);
    ax.XMinorGrid = 'on';
    ax.MinorGridLineStyle = '--';
    ax.XAxis.MinorTickValues = [0 100 200 300 400 500 700 800];
    xlabel('Time Relative to Lever Release');
    plotTitle = sprintf('V1 RT Kernel (n=%d)', nHitsV1);
    title(plotTitle);
    hold off;

    % SC
    subplot(2,1,2);
    bins = size(scBootKernel, 2);
    ax.YColor = [0 0 0];
    plot([0, bins], [0.5, 0.5],  'Color', 'k', 'LineStyle','--');
    hold on;
    plot(scx, scCIs(2, :), 'Color', SCColor, 'LineStyle', '-', 'LineWidth', 2);
    ax.YColor = [0 0 0];
    scx2 = [scx, fliplr(scx)];
    scfillCI = [scCIs(1, :), fliplr(scCIs(3, :))];
    fill(scx2, scfillCI, SCColor, 'lineStyle', '-', 'edgeColor', SCColor, 'edgeAlpha', 0.5, 'faceAlpha', 0.10);
    ax = gca;
    ax.YColor = 'k';
    ax.TickDir = 'out';
    ax.FontSize = 14;
    bins = size(scKernel, 2);
    xlim(ax, [0, bins]);
    %llim2 = ylim;
    ax.XGrid = 'on';
    ylabel(yLabel);
    plotTitle = sprintf('SC RT Kernel (n=%d)', nHitsSC);
    title(plotTitle);
    set(gca,'XTick', [0, -plotStartMS+200, bins]);
    set(gca, 'XTickLabel', {sprintf('%d', plotRTStartMS), '0', sprintf('%d', plotRTStartMS + plotEndMS - plotStartMS)});
    set(gca, 'LineWidth', 1);
    ax.XMinorGrid = 'on';
    ax.MinorGridLineStyle = '--';
    ax.XAxis.MinorTickValues = [0 100 200 300 400 500 700 800];
    xlabel('Time Relative to Lever Release');
    hold off;
    % Put plots on same axes
    %lims = [min(min(llim1, llim2)) max(max(llim1, llim2))];
    lims = [0.47 0.53];
    subplot(2,1,1);
    ylim(lims);
    subplot(2,1,2);
    ylim(lims);  

end



