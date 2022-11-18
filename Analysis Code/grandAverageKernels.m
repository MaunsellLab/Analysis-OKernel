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
assert(ischar(stimType),'stimType must be text: Lum, or Gabor Offset, or FA');

condition = [brainArea ' ' stimType];
[~, analysisDirName, tableName] = whichData(condition);

if contains(computerName(), 'maunsell')
    load([analysisDirName, tableName], 'T');
else
    load([analysisDirName, condition, tableName], 'T');
end

limits = setLimits('All');
limits.rampMS = 0;
 
  if strcmp(brainArea,'SC')
      if strcmp(stimType,'Lum')
          % Full Kernel For Aligned Luminance Patch
          limits.animal = {'1458', '1548', '1674', '1675', '1902', '2057', '2058', '2063', '2169', '2236'};
      elseif strcmp(stimType,'Gabor')
          % Kernel for Aligned Gabor Patch
          limits.animal = {'1548', '1674', '2057', '2058', '2063', '2169'};
      elseif strcmp(stimType,'Offset')
          % Control Kernel For Offset Luminance Patch
          limits.animal = {'1674', '1675', '1902', '2057', '2063', '2036'};
      elseif strcmp(stimType, 'FA')
          % Control Kernel for Locomotor Effects (Combines Lum/Gabor)
          limits.animal = {'1458', '1548', '1674', '1675', '1902', '2057', '2058', '2063', '2169', '2236'};     
      end
      
      
  elseif strcmp(brainArea,'V1')
      if strcmp(stimType,'Lum')
          % V1 Luminance
          limits.animal = {'1960', '2015', '2083', '2126', '2220', '2221'};
      elseif strcmp(stimType,'Gabor')
          % V1 Gabor
          limits.animal = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
          % V1 Gabor Who Also Have Luminance Data
          % limits.animal = {'1960', '2015', '2083', '2126', '2220', '2221'};
      elseif strcmp(stimType,'Offset')
          % V1 Offset
          limits.animal = {'2016','2083', '2126', '2220', '2221'};
      elseif strcmp(stimType, 'FA')
          % Control Kernel for Motor Impact (Combines Luminance and Gabor)
          limits.animal = {'1960', '2015', '2016', '2083', '2126', '2207', '2220', '2221'};
      end
      
  end
  
  
  U = selectUsingLimits(T, limits);
  if size(U, 1) == 0
    return;
  end
  stimProfiles = getOptoProfiles(U, condition);
  plotKernelPage(U, limits, stimProfiles);
  saveas(gcf, sprintf('%sFigures/Kernels/%s.pdf', analysisDirName, limits.animal{1}));
end
