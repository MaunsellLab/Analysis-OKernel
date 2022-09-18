function meanMatchedKernels()

%% Start with Gabor Data
% Specify Data and Make Master Table for Gabor
preProcessAll(30, 'V1 Gabor');
% Load Master Table and Rename
analysisDirName = '/Users/Shared/Analysis/Analysis-OKernel/Mat Files/';
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
%% Run matchDPrime to get indicies of sessions for Mean Matching
[V1GabIdx_DDP, V1LumIdx_DDP, V1GabIdx_DP,...
    V1LumIdx_DP, V1GabIdx_TU, V1LumIdx_TU] = matchDPrime(V1Gab, V1Lum);
%% SubSelect Sessions Based on Matching
% This is dorky and inefficient but otherwise would have to make a bunch of changes to the
% way the stim profile creation is handled for all other analyses. 

% Do for V1 First
preProcessAll(30, 'V1 Gabor');
V1Gab_TU = V1Gab(V1GabIdx_TU,:);
V1Gab_DP = V1Gab(V1GabIdx_DP,:);
V1Gab_DDP = V1Gab(V1GabIdx_DDP,:);

% Stim Profiles Subset According to Mean Matching Conditions
gabProfiles_TU  = getOptoProfiles(V1Gab_TU);
gabProfiles_DP  = getOptoProfiles(V1Gab_DP);
gabProfiles_DDP = getOptoProfiles(V1Gab_DDP);

% Save a copy for later
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/V1Gabor/stimProfiles/gabProfiles_TU.mat','gabProfiles_TU');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/V1Gabor/stimProfiles/gabProfiles_DP.mat','gabProfiles_DP');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/V1Gabor/stimProfiles/gabProfiles_DDP.mat','gabProfiles_DDP');

% Repeat for Luminance
preProcessAll(30, 'V1 Lum');
V1Lum_TU = V1Lum(V1LumIdx_TU,:);
V1Lum_DP = V1Lum(V1LumIdx_DP,:);
V1Lum_DDP = V1Lum(V1LumIdx_DDP,:);

% Stim Profiles Subset According to Mean Matching Conditions
lumProfiles_TU  = getOptoProfiles(V1Lum_TU);
lumProfiles_DP  = getOptoProfiles(V1Lum_DP);
lumProfiles_DDP = getOptoProfiles(V1Lum_DDP);

save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/V1Lum/stimProfiles/lumProfiles_TU.mat','lumProfiles_TU');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/V1Lum/stimProfiles/lumProfiles_DP.mat','lumProfiles_DP');
save('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/V1Lum/stimProfiles/lumProfiles_DDP.mat','lumProfiles_DDP');
%% Make Kernel Plots

limits = setLimits('All');
limits.rampMS = 0;

% Gabor Matched for top up
plotKernelPage(V1Gab_TU, limits, gabProfiles_TU);
saveas(gcf, '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/Figures/V1_Gabor_TU.pdf');
% D-prime
plotKernelPage(V1Gab_DP, limits, gabProfiles_DP);
saveas(gcf, '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/Figures/V1_Gabor_DP.pdf');
% delta d-prime
plotKernelPage(V1Gab_DDP, limits, gabProfiles_DDP);
saveas(gcf, '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/Figures/V1_Gabor_DDP.pdf');


% Now For Luminance
% top up
plotKernelPage(V1Lum_TU, limits, lumProfiles_TU);
saveas(gcf, '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/Figures/V1_Lum_TU.pdf');
% D-prime
plotKernelPage(V1Lum_DP, limits, lumProfiles_DP);
saveas(gcf, '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/Figures/V1_Lum_DP.pdf');
% delta d-prime
plotKernelPage(V1Lum_DDP, limits, lumProfiles_DDP);
saveas(gcf, '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/meanMatched/Figures/V1_Lum_DDP.pdf');


end
