function plotDPrimes(brainArea)
% Plot all d Primes Over Sessions For Each Animal

% brainArea input is a string: 
% 'SC' for superior colliculus  
% 'V1' for visual cortex  


rampMS = 0;
[~, analysisDirName, tableName] = whichData();
load([analysisDirName, tableName], 'T');
limits = setLimits('All');
limits.rampMS = 0;

if strcmp(brainArea,'SC')
    % Full Kernel For Aligned Luminance Patch
    limits.animal = {'1458', '1548', '1634',  '1674', '1675', '1902', '1905'};
    % animals = {'2054','2057', '2058',  '2060', '2061', '2063'};
    % Control Kernel For Offset Luminance Patch
    % animals = {'1674', '1675', '1902'};
    
    % Full Kernel For Aligned Gabor Patch
    % animals = {'1458', '1548', '1674'};
    % animals = {'2054', '2058', '2060', '2063'};
    
elseif strcmp(brainArea,'V1')
    limits.animal = {'1951', '1952', '1956', '1957', '1958', '1960', '1995', '1996', '1998', '2018'};
    % animals = {'1958', '1960', '1995', '1998'};
end

% Select Subset of mice
U = selectUsingLimits(T, limits);
nMice = length(limits.animal);

% init vars
dPrimes=cell(nMice,2);
deltaDPrime = cell(nMice,1);
mice = [U.animal];
sessionCounts = zeros(nMice,1);



for mouseNum = 1:nMice
    mouseIdx=strcmp(limits.animal{1,mouseNum},mice);
    dPrimes{mouseNum,1} = [U.noStimDPrime(mouseIdx)];
    dPrimes{mouseNum,2} = [U.stimDPrime(mouseIdx)];
    deltaDPrime{mouseNum,1} = [U.stimDPrime(mouseIdx)] - [U.noStimDPrime(mouseIdx)];
    sessionCounts(mouseNum) = length(deltaDPrime{mouseNum,1});
end


%%%%%%%%% Plot scatter plot of d-primes from each session (stim v unstim)
% Set Axes Limits
minAxes = 0.5;
maxAxes = 3.2;

% Paired t-test for each mouse d-prime
sigDPrime = zeros(1,nMice);

figure(1);
axis square;
hold on;
% Scatter 
for mouseNum = 1:nMice
    scatter(dPrimes{mouseNum,1},dPrimes{mouseNum,2}, 25, 'filled')
    [~, sigDPrime(1,mouseNum)] = ttest(dPrimes{mouseNum,1},dPrimes{mouseNum,2});  
end
xlim([minAxes maxAxes]);
xlabel('Unstimulated d-prime');
ylim([minAxes maxAxes]);
xlabel('Stimulated d-prime');
plot([minAxes maxAxes],[minAxes maxAxes], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);
title('d-prime');
set(gca, 'FontSize',14);
set(gca, 'TickDir', 'out');
legend(unique(mice), 'Location', 'northwest');
box off
hold off;
    
%%%%%%%%%%%%%% Plot delta d-prime by session to look for changes
minSessions = min(sessionCounts);
dPrimeMat = zeros(nMice, minSessions);

for mouseNum = 1:nMice
    dPrimeMat(mouseNum,:) = deltaDPrime{mouseNum,1}(1:minSessions,1);
end

figure(2);
axis square;
hold on;
plot(1:minSessions, mean(dPrimeMat,1), 'LineWidth', 1.5, 'Color', 'k');
% No Effect Line
plot([0 minSessions+1],[0 0], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.5);

% Error For each point
for nPoints = 1:minSessions
    plot([nPoints nPoints], [mean(dPrimeMat(:,nPoints))+std(dPrimeMat(:,nPoints))/sqrt(nMice)...
        mean(dPrimeMat(:,nPoints))-std(dPrimeMat(:,nPoints))/sqrt(nMice)],...
        'Color', 'k', 'LineWidth', 1);
end

title('d-prime over sessions');
set(gca, 'FontSize',14);
set(gca, 'TickDir', 'out');
xlabel('Session Number');
ylabel('delta d-prime');
ylim([-0.3 0.1]);
xlim([0 minSessions+1]);
box off
hold off;






end