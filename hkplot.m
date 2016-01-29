function hkplot(IPHys)
%HKPLOT plot the line fits for Hk calculation, given analyzed IPHys struct.
    %fig = figure;
    hold all;
    try
        mask1 = IPHys.Field > 0;
        mask2 = IPHys.Field < 0;
        plot(IPHys.Field(mask1), IPHys.Analysis.HkLine1(mask1), ':', 'color', 'red')
        plot(IPHys.Field(mask2), IPHys.Analysis.HkLine2(mask2), ':', 'color', 'red')
        plot(IPHys.Field, IPHys.Analysis.RawMsLine1, ':', 'color', 'red')
        plot(IPHys.Field, IPHys.Analysis.RawMsLine2, ':', 'color', 'red')
        plot(IPHys.Field, IPHys.RawCorrX, '-', 'color', 'blue')
    catch
        disp('hkplot failed')
    end

end