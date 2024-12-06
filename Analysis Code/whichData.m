function [dataDirName, analysisDirName, tableName] = whichData(condition)

% 'Condition specifies what data you are going to batch analyze 
% possible inputs:
% SC Gabor'
% 'SC Lum'
% 'SC Offset'
% 'V1 Gabor'
% 'V1 Lum'
% 'V1 Offset'
% 'V1 Lum MM' (Mean Matched)
% 'V1 Gabor MM' (Mean Matched) 

if contains(computerName(), 'maunsell')
	dataDirName = '../../../Data/OKernel/';
	analysisDirName = '../../Analysis-OKernel/';
    tableName = sprintf('Mat Files/masterTable.mat');
else    
	dataDirName = ['/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/BehavData/' condition '/'];
	%analysisDirName = '../../Analysis-OKernel/';
    analysisDirName = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/MasterFiles/';
    tableName = sprintf('/masterTable.mat');
end

