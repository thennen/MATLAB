% Keep trying to get the optimization license

n = 1;
options = optimset('Display','off');
fail = 1;
while fail
    fail = 0;
    try b = lsqcurvefit(@(a,x) a*sin(x),1,[1,2,3],[1,2,3],0,2,options);
    catch exc
        fail = 1;
        getReport(exc)
        disp(['Opt. Toolbox license checkout failed.  Attempt ' num2str(n)]) 
        n = n + 1;
        pause(10);
        continue
    end
end
disp('Optimization Toolbox License Checkout Succeeded')
clear('options','exc','fail','n','b');