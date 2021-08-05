function individualKernels
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

  rampMS = 0;
  %animals = {'1458', '1462', '1463', '1548'};
  animals = {'1462', '1463'};
  
  [~, analysisDataDir, tableName] = whichData();
  load([analysisDataDir tableName], 'T');
  
  limits = setLimits('All');
  limits.rampMS = rampMS;
  
  for a = 1:length(animals)
    limits.aniNum = a;
    limits.animal = animals{a};
    U = selectUsingLimits(T, limits);
    if height(U) == 0
      continue
    end
    stimProfiles = getOptoProfiles(U);
    plotKernelPage(U, limits, stimProfiles);
    saveas(gcf, sprintf('%sFigures/Kernels/%d %s.pdf', analysisDataDir, limits.rampMS, limits.animal));
  end
  % This isn't plotting anything right now, eventually when we have many
  % mice we will want to plot all the kernels on the same axes.
%   figure(2); 
%   sameAxisScaling('both', 1, 3, 1:length(animals));
%   saveas(gcf, sprintf('%s Analysis/Figures/Kernels/Ramp %d Individuals.pdf', dataDirName, limits.rampMS));
end
