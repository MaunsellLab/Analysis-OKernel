function [dataDirName, analysisDirName, tableName] = whichData(condition)

% 'Condition specifies what data you are going to batch analyze 
% possible inputs:
% SC Gabor'
% 'SC Lum'
% 'SC Offset'
% 'V1 Gabor'
% 'V1 Lum'
% 'V1 Offset'

tableName = sprintf('Mat Files/masterTable.mat');
if contains(computerName(), 'maunsell')
	dataDirName = '../../../Data/OKernel/';
	analysisDirName = '../../Analysis-OKernel/';
else    
	dataDirName = ['/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/BehavData/' condition '/'];
	analysisDirName = '../../Analysis-OKernel/';
end

