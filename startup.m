disp 'lol matlab'

%try 
    %options = optimset('Display','off');
    %a = lsqcurvefit(@(a,x) a*sin(x),1,[1,2,3],[1,2,3],0,2,options);
    %clear options;
    %clear a;
%catch
    %disp(['Opt. Toolbox license checkout failed.']);
    %clear exc;
%end

set(0,'DefaultAxesFontSize',16, ...
      'DefaultTextFontSize',16, ...
      'DefaultLineLineWidth',1.5, ...
      'DefaultAxesBox', 'on');
