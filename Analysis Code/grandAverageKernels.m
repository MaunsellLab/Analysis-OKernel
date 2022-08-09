function grandAverageKernels(brainArea, stimType)
% Plot all the kernels for the grand average across selected sessions

% brainArea input is a string:
% 'SC' for superior colliculus
% 'V1' for visual cortex

% stimType input is also a string:
% 'Lum' for luminance patch
% 'Gab' for Gabor patch
% 'Off' for Offset Control (misaligned)

% if not, throw errors
assert(ischar(brainArea),'brainArea must be text: SC or V1');
assert(ischar(stimType),'stimType must be text: Lum, or Gab or Off');
 
 
  [~, analysisDirName, tableName] = whichData();
  load([analysisDirName, tableName], 'T');
  limits = setLimits('All');
  limits.rampMS = 0;
 
  if strcmp(brainArea,'SC')
      
      if strcmp(stimType,'Lum')
          % Full Kernel For Aligned Luminance Patch
          limits.animal = {'1458', '1548', '1674', '1675', '1902', '2057', '2058', '2063', '2169', '2236'};
      elseif strcmp(stimType,'Gab')
          % Kernel for Aligned Gabor Patch
          limits.animal = {'1548', '1674', '2057', '2058', '2063', '2169'};
      elseif strcmp(stimType,'Off')
          % Control Kernel For Offset Luminance Patch
          limits.animal = {'1674', '1675', '1902', '2057', '2063'};
      end
      
      
  elseif strcmp(brainArea,'V1')
      if strcmp(stimType,'Lum')
          % V1 Luminance
          limits.animal = {'1960','2015','2083','2126'};
      elseif strcmp(stimType,'Gab')
          % V1 Gabor
          limits.animal = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
      elseif strcmp(stimType,'Off')
          % V1 Offset
          limits.animal = {'2016','2083', '2126', '2220'};
      end
      
  end
  
  
  U = selectUsingLimits(T, limits);
  if size(U, 1) == 0
    return;
  end
  stimProfiles = getOptoProfiles(U);
  plotKernelPage(U, limits, stimProfiles);
  saveas(gcf, sprintf('%sFigures/Kernels/%s.pdf', analysisDirName, limits.animal{1}));
end
