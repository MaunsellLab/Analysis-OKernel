function grandAverageKernels(brainArea)
% Plot all the kernels for the grand average across selected sessions

% brainArea input is a string: 
% 'SC' for superior colliculus  
% 'V1' for visual cortex

  [~, analysisDirName, tableName] = whichData();
  load([analysisDirName, tableName], 'T');
  limits = setLimits('All');
  limits.rampMS = 0;
 
  if strcmp(brainArea,'SC')
      % Full Kernel For Aligned Luminance Patch
      limits.animal = {'1458', '1548', '1674', '1675', '1902', '1905'};
      % VGATs
      %limits.animal = {'2057', '2058', '2063'};
      
      % Control Kernel For Offset Luminance Patch
      % limits.animal = {'1674', '1675', '1902'};
  elseif strcmp(brainArea,'V1')
      limits.animal = {'1951', '1956', '1960', '1995', '2015', '2016', '2018', '2080'};
  end
  
  
  U = selectUsingLimits(T, limits);
  if size(U, 1) == 0
    return;
  end
  stimProfiles = getOptoProfiles(U);
  plotKernelPage(U, limits, stimProfiles);
  saveas(gcf, sprintf('%sFigures/Kernels/%s.pdf', analysisDirName, limits.animal{1}));
end
