function sharrockplot(identifier)
% Plot Sharrock extrapolation Hc vs log(1/sweeprate)

basevars = evalin('base','who');
varnames = findbestmatch(basevars,[identifier '_PKDynamic']);
varname = varnames{1};
SharrockExtrap = evalin('base',[varname '.SharrockExtrap']);
Analysis = evalin('base',[varname '.Analysis']);

exponent = Analysis.n;
exponent = num2cell(exponent);
parse = @(x) ['n = ' num2str(x)];
exponent = cellfun(parse,exponent,'UniformOutput',false);
%f0 = Analysis.f0;
%f0 = num2cell(f0);
%f0 = cellfun(@num2str, f0, 'UniformOutput', false);


figure;
plot(log(SharrockExtrap.InverseRate), SharrockExtrap.Hc, '.', 'MarkerSize', 15, 'DisplayName', 'PKdata');
xlabel('log(1/SweepRate)');
ylabel('Hc');
title([varname ' Sharrock Fit']);
hold all;
% for now, assuming only exponent can change
plot(log(SharrockExtrap.Fit_InverseRate),real(SharrockExtrap.Fit_Hc), 'DisplayName', exponent)


end

