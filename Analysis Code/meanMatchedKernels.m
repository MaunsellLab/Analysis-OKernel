function meanMatchedKernels()

%% Start with Gabor Data
% Specify Data and Make Master Table for Gabor
preProcessAll(30, 'V1 Gabor');
% Load Master Table and Rename
analysisDirName = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';
load([analysisDirName 'masterTable.mat'], 'T');
V1Gab = T;
clear T;
% Save Gabor Master Table to Dropbox
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/V1Gabor/masterTable/V1Gab.mat','V1Gab');
%% Repeat for Luminance Data
preProcessAll(30, 'V1 Lum');
% Load Master Table and Rename
load([analysisDirName 'masterTable.mat'], 'T');
V1Lum = T;
clear T;
% Save Luminance Master Table to Dropbox
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/V1Lum/masterTable/V1Lum.mat','V1Lum');
%% Use Only Mice that have done sessions both stimuli
gabAnimals = unique([V1Gab.animal]);
lumAnimals = unique([V1Lum.animal]);
mice = intersect(lumAnimals,gabAnimals);

gabSessLog = zeros(1, height(V1Gab));
lumSessLog = zeros(1, height(V1Lum));

for i = 1:height(V1Gab)
    if ismember(V1Gab.animal(i), mice)
        gabSessLog(1,i) = 1;
    else 
        continue
    end
end

for i = 1:height(V1Lum)
    if ismember(V1Lum.animal(i), mice)
        lumSessLog(1,i) = 1;
    else
        continue
    end
end

% Convert To Logical
gabSessLog = logical(gabSessLog)';
lumSessLog = logical(lumSessLog)';

% Sub-Select
V1Gab = V1Gab(gabSessLog,:);
V1Lum = V1Lum(lumSessLog,:);

%% Run matchDPrime to get indicies of sessions for Mean Matching
[V1GabIdx_DDP, V1LumIdx_DDP, V1GabIdx_DP,...
    V1LumIdx_DP, V1GabIdx_TU, V1LumIdx_TU] = matchDPrime(V1Gab, V1Lum);\
%% SubSelect Sessions Based on Matching
% This is dorky and inefficient but otherwise would have to make a bunch of changes to the
% way the stim profile creation is handled for all other analyses. 

% Do for V1 First
% preProcessAll(30, 'V1 Gabor');
V1Gab_TU = V1Gab(V1GabIdx_TU,:);
V1Gab_DP = V1Gab(V1GabIdx_DP,:);
V1Gab_DDP = V1Gab(V1GabIdx_DDP,:);

% Stim Profiles Subset According to Mean Matching Conditions
gabProfiles_TU  = getOptoProfiles(V1Gab_TU, 'V1 Gabor MM');
gabProfiles_DP  = getOptoProfiles(V1Gab_DP, 'V1 Gabor MM');
gabProfiles_DDP = getOptoProfiles(V1Gab_DDP, 'V1 Gabor MM');

% Save a copy for later
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Gabor MM/stimProfiles/gabProfiles_TU.mat','gabProfiles_TU');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Gabor MM/stimProfiles/gabProfiles_DP.mat','gabProfiles_DP');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Gabor MM/stimProfiles/gabProfiles_DDP.mat','gabProfiles_DDP');
% Save Table
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Gabor MM/masterTable.mat','V1Gab_DDP');

% Repeat for Luminance
% preProcessAll(30, 'V1 Lum');
V1Lum_TU = V1Lum(V1LumIdx_TU,:);
V1Lum_DP = V1Lum(V1LumIdx_DP,:);
V1Lum_DDP = V1Lum(V1LumIdx_DDP,:);

% Stim Profiles Subset According to Mean Matching Conditions
lumProfiles_TU  = getOptoProfiles(V1Lum_TU, 'V1 Lum MM');
lumProfiles_DP  = getOptoProfiles(V1Lum_DP, 'V1 Lum MM');
lumProfiles_DDP = getOptoProfiles(V1Lum_DDP, 'V1 Lum MM');

save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Lum MM/stimProfiles/lumProfiles_TU.mat','lumProfiles_TU');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Lum MM/stimProfiles/lumProfiles_DP.mat','lumProfiles_DP');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Lum MM/stimProfiles/lumProfiles_DDP.mat','lumProfiles_DDP');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/V1 Lum MM/masterTable.mat','V1Lum_DDP');
end




