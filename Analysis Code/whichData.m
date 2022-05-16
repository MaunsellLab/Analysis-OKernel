function [dataDirName, analysisDirName, tableName] = whichData()

% tableName = sprintf('Analysis/Mat Files/masterTable.mat');
% dirName = '/Users/Shared/Data/SCernel/';

tableName = sprintf('Mat Files/masterTable.mat');
if contains(computerName(), 'maunsell')
	dataDirName = '../../../Data/OKernel/';
	analysisDirName = '../../Analysis-OKernel/';
else
	dataDirName = '/Users/jacksoncone/Dropbox/PostDoctoral Projects/!Experiments/Colliculus/BehavData/SC Lum/';
	analysisDirName = '../../Analysis-OKernel/';
end

