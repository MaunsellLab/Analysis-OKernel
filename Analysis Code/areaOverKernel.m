function [~] = areaOverKernel(brainArea)
% Compute Area Over Kernel for Individual Mice 
% Bootstrap 95%CI per mouse 
% Plot Scatter of results

% Location of the Final Stim Profiles and Tables
analysisDir = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';
lumFolder = strcat(analysisDir, brainArea,' Lum/');
gabFolder = strcat(analysisDir, brainArea,' Gabor/');

% Load Tables for each stimulus type
load([lumFolder 'masterTable.mat'], 'T');
lumT = T;
load([gabFolder 'masterTable.mat'], 'T');
gabT = T;
clear T;

% Use Limits to Subselect Sessions From Each Mouse
limits = setLimits('All');
rampMS = 0;
limits.rampMS = rampMS;

% Only Consider Animals For Which Data was obtained for both stimuli
if strcmp(brainArea,'SC')
    animals = {'1548', '1674', '2057', '2058', '2063', '2169', '2236'};
elseif strcmp(brainArea,'V1')
    animals = {'1960', '2015', '2083', '2126', '2220', '2221'};
end

% Loop Through this for Both Stimulus Types

% First For Luminance
for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(lumT, limits);
    condition = [brainArea ' ' 'Lum'];
    
    % Extract optoProfiles
    lumProfiles.hitProfiles = [];
    lumProfiles.missProfiles = [];
    lumProfiles.earlyProfiles = [];
    lumProfiles.RTProfiles = [];
    lumProfiles.stimRTProfiles = [];

    % get all Profiles from this mouse
    for i = 1:height(U)
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', U.animal(i), '/', U.date(i), '.mat'), 'stimProfiles');
    % Append Profiles
        lumProfiles.hitProfiles = [lumProfiles.hitProfiles; stimProfiles.hitProfiles];
        lumProfiles.missProfiles = [lumProfiles.missProfiles; stimProfiles.missProfiles];
        lumProfiles.earlyProfiles = [lumProfiles.earlyProfiles; stimProfiles.earlyProfiles];
        lumProfiles.RTProfiles = [lumProfiles.RTProfiles; stimProfiles.RTProfiles];
        lumProfiles.stimRTProfiles = [lumProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
    end








end





% Load Table and Directory for Stim Profiles
% First Grab all the stim Profiles one mouse then create kernel
% Sum difference from 0.5 across first 100 ms (maybe tune this number)
% Bootstrap CIs from Stim Profiles (bootNum = 10000) 






end


