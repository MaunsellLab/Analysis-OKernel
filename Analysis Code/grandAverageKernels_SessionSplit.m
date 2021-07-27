function grandAverageKernels_SessionSplit
% Plot all the kernels for the grand average across selected sessions, once for steps and once for ramps

  [~, analysisDirName, tableName] = whichData();
  load([analysisDirName, tableName], 'T');
  limits = setLimits('All');
  limits.rampMS = 0;
  limits.half{1} = '1st Half';
  limits.half{2} = '2nd Half';
  U = selectUsingLimits(T, limits);
  if size(U, 1) == 0
    return;
  end
  stimProfiles = getOptoProfiles(U);
  
  % Grab the outcomes from each session, will use to split sessions in half
  nSessions = size(U,1);
  nHits = U.stimCorrects;
  nMiss = U.stimFails;
  nFA = U.stimEarlies;

  % Kernel Profiles from 1st half of sessions
  profiles1.hitProfiles = [];
  profiles1.missProfiles = [];
  profiles1.earlyProfiles = [];
  profiles1.RTProfiles = [];
  profiles1.stimRTProfiles = [];
  
  % Kernel Profiles from 2nd half of sessions
  profiles2.hitProfiles = [];
  profiles2.missProfiles = [];
  profiles2.earlyProfiles = [];
  profiles2.RTProfiles = [];
  profiles2.stimRTProfiles = [];
  
  hitCounter = 1;
  missCounter = 1;
  faCounter = 1;
  
  for i = 1:nSessions

      % Hits
      profiles1.hitProfiles = [profiles1.hitProfiles;...
          stimProfiles.hitProfiles(hitCounter:hitCounter+floor(nHits(i)/2),:)];
      
      profiles2.hitProfiles = [profiles2.hitProfiles;...
          stimProfiles.hitProfiles(hitCounter+floor(nHits(i)/2)+1:hitCounter+floor(nHits(i))-1,:)];
      
      % Aligned to RT
      profiles1.RTProfiles = [profiles1.RTProfiles;...
          stimProfiles.RTProfiles(hitCounter:hitCounter+floor(nHits(i)/2),:)];
      
      profiles2.RTProfiles = [profiles2.RTProfiles;...
          stimProfiles.RTProfiles(hitCounter+floor(nHits(i)/2)+1:hitCounter+floor(nHits(i))-1,:)];
      
      % Aligned to Stim-RT
      profiles1.stimRTProfiles = [profiles1.stimRTProfiles;...
          stimProfiles.stimRTProfiles(hitCounter:hitCounter+floor(nHits(i)/2),:)];
      
      profiles2.stimRTProfiles = [profiles2.stimRTProfiles;...
          stimProfiles.stimRTProfiles(hitCounter+floor(nHits(i)/2)+1:hitCounter+floor(nHits(i))-1,:)];
      
      hitCounter = hitCounter + nHits(i);
      
      % Misses
      profiles1.missProfiles = [profiles1.missProfiles;...
          stimProfiles.missProfiles(missCounter:missCounter+floor(nMiss(i)/2),:)];
          
      profiles2.missProfiles = [profiles2.missProfiles;...
          stimProfiles.missProfiles(missCounter+floor(nMiss(i)/2)+1:missCounter+floor(nMiss(i))-1,:)];

      missCounter = missCounter + nMiss(i);
  
      % FA
      profiles1.earlyProfiles = [profiles1.earlyProfiles;...
          stimProfiles.earlyProfiles(faCounter:faCounter+floor(nFA(i)/2),:)];
      
      profiles2.earlyProfiles = [profiles2.earlyProfiles;...
          stimProfiles.earlyProfiles(faCounter+floor(nFA(i)/2)+1:faCounter+floor(nFA(i))-1,:)];
      
      faCounter = faCounter + nFA(i);
  end
  
 
  for i = 1:2
      if i == 1
          plotKernelPage(U, limits, profiles1);
          saveas(gcf, sprintf('%sFigures/Kernels/%s.pdf', analysisDirName, limits.half{i}));
      else
          plotKernelPage(U, limits, profiles2);
          saveas(gcf, sprintf('%sFigures/Kernels/%s.pdf', analysisDirName, limits.half{i})); 
          
      end
      

  end
  

end
