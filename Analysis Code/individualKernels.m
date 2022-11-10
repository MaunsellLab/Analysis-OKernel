function individualKernels(brainArea, stimType)
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

% brainArea input is a string:
% 'SC' for superior colliculus
% 'V1' for visual cortex

% stimType input is also a string:
% 'Lum' for luminance patch
% 'Gabor' for Gabor patch
% 'Off' for Offset Control (misaligned)
% 'FA' for combined Gabor/Lum to look at impact on motor

% if not, throw errors
assert(ischar(brainArea),'brainArea must be text: SC or V1');
assert(ischar(stimType),'stimType must be text: Lum, or Gabor or Offset or FA');


rampMS = 0;

condition = [brainArea ' ' stimType];
[~, analysisDataDir, tableName] = whichData(condition);

if contains(computerName(), 'maunsell')
    load([analysisDirName, tableName], 'T');
else
    load([analysisDirName, condition, tableName], 'T');
end


limits = setLimits('All');
limits.rampMS = rampMS;

if strcmp(brainArea,'SC')

    if strcmp(stimType,'Lum')
        % Full Kernel For Aligned Luminance Patch
        animals = {'1458', '1548', '1674', '1675', '1902', '2057', '2058', '2063', '2169', '2236'};
        % animals = {'2063'};
    elseif strcmp(stimType,'Gabor')
        % Full Kernel For Aligned Gabor Patch
        animals = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
    elseif strcmp(stimType, 'Offset')
        % Control Kernel For Offset Luminance Patch
        animals = {'1674', '1675', '1902', '2057', '2063', '2236'};
        % animals = {'2057', '2063', '2236'};
    elseif strcmp(stimType, 'FA')
        % Control Kernel For Motor Impact (combined across visual stimuli)
        animals = {'1458', '1548', '1674', '1675', '1902', '2057', '2058', '2063', '2169', '2236'};
    end


elseif strcmp(brainArea,'V1')

    if strcmp(stimType,'Lum')
        % V1 Luminance
        animals = {'1960', '2015', '2083', '2126', '2220', '2221'};
    elseif strcmp(stimType,'Gabor')
        % V1 Gabor
        animals = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
    elseif strcmp(stimType,'Offset')
        % V1 Offset
        animals = {'2016', '2083', '2126', '2220', '2221'};
    elseif strcmp(stimType, 'FA')
        % Control Kernel For Motor Impact (Combined Across Stimuli)
        animals = {'1960', '2015', '2016', '2083', '2126', '2206', '2207', '2220', '2221'};
    end
end


for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(T, limits);
    if height(U) == 0
        continue
    end
    stimProfiles = getOptoProfiles(U);
    plotKernelPage(U, limits, stimProfiles);
    saveas(gcf, sprintf('%sFigures/Kernels/%d %s.pdf', analysisDataDir, limits.rampMS, limits.animal));
end
% This isn't plotting anything right now, eventually when we have many
% mice we will want to plot all the kernels on the same axes.
%   figure(2);
%   sameAxisScaling('both', 1, 3, 1:length(animals));
%   saveas(gcf, sprintf('%s Analysis/Figures/Kernels/Ramp %d Individuals.pdf', dataDirName, limits.rampMS));
end
