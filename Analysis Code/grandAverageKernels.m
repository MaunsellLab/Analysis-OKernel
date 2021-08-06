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
      limits.animal = {'1458', '1520', '1548'};
  elseif strcmp(brainArea,'V1')
      limits.animal = {'1462', '1463'};
  end
  
  
  U = selectUsingLimits(T, limits);
  if size(U, 1) == 0
    return;
  end
  stimProfiles = getOptoProfiles(U);
  plotKernelPage(U, limits, stimProfiles);
  saveas(gcf, sprintf('%sFigures/Kernels/%s.pdf', analysisDirName, limits.animal{1}));
end
