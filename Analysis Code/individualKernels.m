function individualKernels(brainArea)
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

  % brainArea input is a string:
  % 'SC' for superior colliculus
  % 'V1' for visual cortex

  rampMS = 0;
  
  [~, analysisDataDir, tableName] = whichData();
  load([analysisDataDir tableName], 'T');
  
  limits = setLimits('All');
  limits.rampMS = rampMS;
  
  if strcmp(brainArea,'SC')
      % Full Kernel For Aligned Luminance Patch
      % animals = {'1458', '1548', '1634',  '1674', '1675', '1902', '1905'};
      animals = {'2057', '2058', '2063'};
      % Control Kernel For Offset Luminance Patch
      % animals = {'1674', '1675', '1902'};
      % animals = {'1634'};
      
      % Full Kernel For Aligned Gabor Patch
      % animals = {'1458', '1548', '1674'};
      % animals = {'2054', '2058', '2060', '2063'};
      
  elseif strcmp(brainArea,'V1')
      animals = {'1956', '1960', '1995', '1996', '1998', '2015', '2016', '2017', '2018', '2078', '2080', '2083', '2086', '2126', '2128'};
      % animals = {'1958', '1960', '1995', '1998'};
  end
  
 
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
