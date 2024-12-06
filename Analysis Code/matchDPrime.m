function [V1GabIdx_DDP, V1LumIdx_DDP, V1GabIdx_DP, V1LumIdx_DP, V1GabIdx_TU, V1LumIdx_TU] = matchDPrime(V1Gab, V1Lum);

% Outputs: logical indexes for mean matched d prime sessions
% Does so for matching based on
% d-prime on TopUp (TU)
% d-prime on unstimulated test Stimulus (DP)
% delta d-prime stim-unstim (DDP)

% Inputs: Data MasterTables for Gabor (Gab) and Luminance (Lum)
% Create MasterTables using whichdata() and preProcessAll

% Use Outputs to subselect using getOptoProfiles.

% Compare Kernels for Mean Matched d' Primes
% First You need to create master tables for each stimulus type
% Use preProcessAll and whichData to specify
% Table created in the Analysis-OKernel/Mat Files/masterTable.m
% Transfer and save them to a new folder and rename according to stim type
% Here they are labeled V1Gab and V1Lum

%% 
nGabSessions = size(V1Gab,1);
nLumSessions = size(V1Lum,1);

% Extract d-primes

% Top Up
gabTopUpDP = V1Gab.topUpDPrime;
lumTopUpDP = V1Lum.topUpDPrime;
% Unstimulated
gabUnstimDP = V1Gab.noStimDPrime;
lumUnstimDP = V1Lum.noStimDPrime;
% Delta d'
gabDeltaDP = V1Gab.stimDPrime - V1Gab.noStimDPrime; 
lumDeltaDP = V1Lum.stimDPrime - V1Lum.noStimDPrime; 

%% Examine Distribution of d-primes
figure('Position', [100 100 2000 1000]);

% Top Up
subplot(1,3,1);
axis square;
hold on;
edges = floor(min([gabTopUpDP; lumTopUpDP])):0.25:ceil(max([gabTopUpDP; lumTopUpDP]));
histogram(gabTopUpDP, edges, 'Normalization', 'probability');
histogram(lumTopUpDP, edges, 'Normalization', 'probability');
ylabel('probability');
xlabel('d-prime');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
box off;
legend('V1 gabor', 'V1 luminance');
title('d-prime: Top Up');
hold off;

% Unstimulated
subplot(1,3,2);
axis square;
hold on;
edges = floor(min([gabUnstimDP; lumUnstimDP])):0.25:ceil(max([gabUnstimDP; lumUnstimDP]));
histogram(gabUnstimDP, edges, 'Normalization', 'probability');
histogram(lumUnstimDP, edges, 'Normalization', 'probability');
ylabel('probability');
xlabel('d-prime');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
box off;
legend('gabor', 'luminance');
legend('V1 gabor', 'V1 luminance');
hold off;

% Delta d'
subplot(1,3,3);
axis square;
hold on;
edges = floor(min([gabDeltaDP; lumDeltaDP])):0.1:ceil(max([gabDeltaDP; lumDeltaDP]));
histogram(gabDeltaDP, edges, 'Normalization', 'probability');
histogram(lumDeltaDP, edges, 'Normalization', 'probability');
ylabel('probability');
xlabel('d-prime');
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'TickDir', 'out');
box off;
legend('V1 gabor', 'V1 luminance');
title('delta d-prime');
hold off;

%% Mean Matched For Top Up

% Top Up
edges = floor(min([gabTopUpDP; lumTopUpDP])):0.25:ceil(max([gabTopUpDP; lumTopUpDP]));
[nLumTopUp,~] = histcounts(lumTopUpDP, edges);
[nGabTopUp,~] = histcounts(gabTopUpDP, edges);
nTopUps = min([nLumTopUp; nGabTopUp],[],1);

% Init
gabLocs  = [];
lumLocs  = [];

% Draw Sessions From Tables That Meet These Criteria
for edge = 1:length(nTopUps);
    if nTopUps(edge) == 0
        continue
    else
        % How Many Sessions to Draw in this bin?
        numVals = nTopUps(edge);
        % Create Logical For the Sessions in This Range
        V1GabIdx = V1Gab.topUpDPrime >= edges(edge) & V1Gab.topUpDPrime < edges(edge+1)';
        V1LumIdx = V1Lum.topUpDPrime >= edges(edge) & V1Lum.topUpDPrime < edges(edge+1)';
        % Locations where true, subsampled for mean matching
        gabLocs = [gabLocs; randsample(find(V1GabIdx==1),numVals)];
        lumLocs = [lumLocs; randsample(find(V1LumIdx==1), numVals)];
    end
end

% Init
V1GabIdx_TU = logical(zeros(size(V1Gab,1),1));
V1LumIdx_TU = logical(zeros(size(V1Lum,1),1));
% Change for sampled values
V1GabIdx_TU(gabLocs)=1;
V1LumIdx_TU(lumLocs)=1;   

%% Mean Matched For Unstimulated Test Stimulus

% Unstim d-prime
edges = floor(min([gabUnstimDP; lumUnstimDP])):0.25:ceil(max([gabUnstimDP; lumUnstimDP]));
[nLumUnstimDP,~] = histcounts(lumUnstimDP, edges);
[nGabUnstimDP,~] = histcounts(gabUnstimDP, edges);
nUnstim = min([nLumUnstimDP; nGabUnstimDP],[],1);

% Init
gabLocs  = [];
lumLocs  = [];

% Draw Sessions From Tables That Meet These Criteria
for edge = 1:length(nUnstim);
 
    if nUnstim(edge) == 0
        continue
    else
        % How Many Sessions to Draw in this bin?
        numVals = nUnstim(edge);
        % Create Logical For the Sessions in This Range
        V1GabIdx = V1Gab.noStimDPrime >= edges(edge) & V1Gab.noStimDPrime < edges(edge+1)';
        V1LumIdx = V1Lum.noStimDPrime >= edges(edge) & V1Lum.noStimDPrime < edges(edge+1)';
        % Locations where true, subsampled for mean matching
        gabLocs = [gabLocs; randsample(find(V1GabIdx==1),numVals)];
        lumLocs = [lumLocs; randsample(find(V1LumIdx==1), numVals)];
    end
end

% Init
V1GabIdx_DP = logical(zeros(size(V1Gab,1),1));
V1LumIdx_DP = logical(zeros(size(V1Lum,1),1));
% Change for sampled values
V1GabIdx_DP(gabLocs)=1;
V1LumIdx_DP(lumLocs)=1; 

%% Delta d-prime
edges = floor(min([gabDeltaDP; lumDeltaDP])):0.1:ceil(max([gabDeltaDP; lumDeltaDP]));
[nLumDeltaDP,~] = histcounts(lumDeltaDP, edges);
[nGabDeltaDP,~] = histcounts(gabDeltaDP, edges);
ndeltaDP = min([nLumDeltaDP; nGabDeltaDP],[],1);

% Init
gabLocs  = [];
lumLocs  = [];

% Draw Sessions From Tables That Meet These Criteria
for edge = 1:length(ndeltaDP);
 
    if ndeltaDP(edge) == 0
        continue
    else
        % How Many Sessions to Draw in this bin?
        numVals = ndeltaDP(edge);
        % Create Logical For the Sessions in This Range
        V1GabIdx = (V1Gab.stimDPrime - V1Gab.noStimDPrime) >= edges(edge) & (V1Gab.stimDPrime - V1Gab.noStimDPrime) < edges(edge+1)';
        V1LumIdx = (V1Lum.stimDPrime - V1Lum.noStimDPrime) >= edges(edge) & (V1Lum.stimDPrime - V1Lum.noStimDPrime) < edges(edge+1)';
        % Locations where true, subsampled for mean matching
        gabLocs = [gabLocs; randsample(find(V1GabIdx==1),numVals)];
        lumLocs = [lumLocs; randsample(find(V1LumIdx==1), numVals)];
    end
end

% Init
V1GabIdx_DDP = logical(zeros(size(V1Gab,1),1));
V1LumIdx_DDP = logical(zeros(size(V1Lum,1),1));
% Change for sampled values
V1GabIdx_DDP(gabLocs)=1;
V1LumIdx_DDP(lumLocs)=1; 


end


