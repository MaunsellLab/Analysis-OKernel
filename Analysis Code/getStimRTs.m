function [noStimRTs, stimRTs] = getStimRTs(brainArea)

[~, analysisDirName, tableName] = whichData();
load([analysisDirName, tableName], 'T');
limits = setLimits('All');
limits.rampMS = 0;

if strcmp(brainArea,'SC')
    % Full Kernel For Aligned Luminance Patch
    limits.animal = {'1458', '1548', '1674', '1675', '1902', '1905'};
    % VGATs
    %limits.animal = {'2057', '2058', '2063'};
    
    % Control Kernel For Offset Luminance Patch
    % limits.animal = {'1674', '1675', '1902'};
elseif strcmp(brainArea,'V1')
    limits.animal = {'1951', '1956', '1960', '1995', '2015', '2016', '2018', '2080'};
end

U = selectUsingLimits(T, limits);

if size(U, 1) == 0
    return;
end

stim = [U.stimCorrectRTs];
noStim = [U.noStimCorrectRTs];
faStim = [U.stimEarlyRTs];
faNoStim = [U.noStimEarlyRTs];

stimRTs = [];
noStimRTs = [];

sessions = length(stim);

for sessionNum = 1:sessions
    stimRTs = [stimRTs, stim{sessionNum,1}, faStim{sessionNum,1}];
    noStimRTs = [noStimRTs, stim{sessionNum,1}, faNoStim{sessionNum,1}];
end

edges = -1000:20:1000;
% Make Histogram of all reaction times
figure;
ax = gca;
axis square;
hold on;
ax.YAxis.Visible = 'off';
yyaxis right;
histogram(stimRTs, edges);
xlabel('time (ms)');
ylabel('counts');
xlim([-400 400]);
ax.TickDir = 'out';
ax.YLim = [0 5000];
ax.FontSize = 14;
ax.YColor = 'k';
ax.LineWidth = 1;
ax.XTick = [-400 -200 0 200 400];
yl = ylim;
ax.YTick = [yl];
title('Hit Reaction Times');
box off;
hold off;

% Compare CDFs
figure;
ax = gca;
axis square;
hold on;
cdfplot(noStimRTs);
cdfplot(stimRTs);
xlabel('time (ms)');
ylabel('Cumulative Probability of Behavioral Response');
xlim([-400 400]);
ax.TickDir = 'out';
ax.FontSize = 14;
ax.YColor = 'k';
ax.LineWidth = 1;
%ax.XTick = [-400 -200 0 200 400];
ax.YLim = [0 1];
%ax.YTick = [yl];
title('RT CDFs');
legend('noStim', 'Stim');
box off;
hold off;


end





