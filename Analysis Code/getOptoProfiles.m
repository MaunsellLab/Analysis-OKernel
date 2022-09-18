function optoProfiles = getOptoProfiles(T)
            
optoProfiles.hitProfiles = [];
optoProfiles.missProfiles = [];
optoProfiles.earlyProfiles = [];
optoProfiles.RTProfiles = [];
optoProfiles.stimRTProfiles = [];
% [~, analysisDirName] = whichData();

for i = 1:height(T)
  load(strcat('/Users/Shared/Analysis/Analysis-OKernel', '/Mat Files/Stim Profiles/', T.animal(i), '/', T.date(i), '.mat'), 'stimProfiles'); 
  optoProfiles.hitProfiles = [optoProfiles.hitProfiles; stimProfiles.hitProfiles];
  optoProfiles.missProfiles = [optoProfiles.missProfiles; stimProfiles.missProfiles];
  optoProfiles.earlyProfiles = [optoProfiles.earlyProfiles; stimProfiles.earlyProfiles];
  optoProfiles.RTProfiles = [optoProfiles.RTProfiles; stimProfiles.RTProfiles];
  optoProfiles.stimRTProfiles = [optoProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end