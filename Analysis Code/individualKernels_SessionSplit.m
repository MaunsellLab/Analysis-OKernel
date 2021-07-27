function individualKernels_SessionSplit
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

rampMS = 0;
animals = {'1462', '1463'};

[~, analysisDataDir, tableName] = whichData();
load([analysisDataDir tableName], 'T');

limits = setLimits('All');
limits.rampMS = rampMS;
limits.half{1} = '1st Half';
limits.half{2} = '2nd Half';
for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(T, limits);
    if height(U) == 0
        continue;
    end
    stimProfiles = getOptoProfiles(U);
    
    % Num Trials to Split On
    nHits = floor(size(stimProfiles.hitProfiles,1)/2);
    nMiss = floor(size(stimProfiles.missProfiles,1)/2);
    nFA = floor(size(stimProfiles.earlyProfiles,1)/2);
    
    % Split stimProfiles by half
    for i = 1:2
        if i == 1
            profiles.hitProfiles = stimProfiles.hitProfiles(1:nHits,:);
            profiles.missProfiles = stimProfiles.missProfiles(1:nMiss,:);
            profiles.earlyProfiles = stimProfiles.earlyProfiles(1:nFA,:);
            profiles.RTProfiles = stimProfiles.RTProfiles(1:nHits,:);
            profiles.stimRTProfiles = stimProfiles.stimRTProfiles(1:nHits,:);
            
            plotKernelPage(U, limits, profiles);
            saveas(gcf, sprintf('%sFigures/Kernels/%d %s %s.pdf', analysisDataDir, limits.animal, limits.half{i}));
        else
            profiles.hitProfiles = stimProfiles.hitProfiles(nHits+1:end,:);
            profiles.missProfiles = stimProfiles.missProfiles(nMiss+1:end,:);
            profiles.earlyProfiles = stimProfiles.earlyProfiles(nFA+1:end,:);
            profiles.RTProfiles = stimProfiles.RTProfiles(nHits+1:end,:);
            profiles.stimRTProfiles = stimProfiles.stimRTProfiles(nHits+1:end,:);
            
            plotKernelPage(U, limits, profiles);
            saveas(gcf, sprintf('%sFigures/Kernels/%d %s %s.pdf', analysisDataDir, limits.animal, limits.half{i}));

        end
        
    end
end
end
