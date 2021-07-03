function grandAverageKernels
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

  [~, analysisDirName, tableName] = whichData();
  load([analysisDirName, tableName], 'T');
  limits = setLimits('All');
  limits.rampMS = 0;
  U = selectUsingLimits(T, limits);
  if size(U, 1) == 0
    return;
  end
  stimProfiles = getOptoProfiles(U);
  plotKernelPage(U, limits, stimProfiles);
  saveas(gcf, sprintf('%sFigures/Kernels/%s.pdf', analysisDirName, limits.animal{1}));
end
