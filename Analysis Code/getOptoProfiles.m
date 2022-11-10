function optoProfiles = getOptoProfiles(T)
            
optoProfiles.hitProfiles = [];
optoProfiles.missProfiles = [];
optoProfiles.earlyProfiles = [];
optoProfiles.RTProfiles = [];
optoProfiles.stimRTProfiles = [];

if contains(computerName(), 'maunsell')
    [~, analysisDirName] = whichData();
else 
    [~, analysisDirName] = whichData(condition);
end

for i = 1:height(T)
    if contains(computerName(), 'maunsell')
        load(strcat('/Users/Shared/Analysis/Analysis-OKernel', '/Mat Files/Stim Profiles/', T.animal(i), '/', T.date(i), '.mat'), 'stimProfiles');
    else
        load(strcat('/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/',...
            condition, '/Stim Profiles/', T.animal(i), '/', T.date(i), '.mat'), 'stimProfiles');
    end

    optoProfiles.hitProfiles = [optoProfiles.hitProfiles; stimProfiles.hitProfiles];
    optoProfiles.missProfiles = [optoProfiles.missProfiles; stimProfiles.missProfiles];
    optoProfiles.earlyProfiles = [optoProfiles.earlyProfiles; stimProfiles.earlyProfiles];
    optoProfiles.RTProfiles = [optoProfiles.RTProfiles; stimProfiles.RTProfiles];
    optoProfiles.stimRTProfiles = [optoProfiles.stimRTProfiles; stimProfiles.stimRTProfiles];
end