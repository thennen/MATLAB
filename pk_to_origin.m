function originout = pk_to_origin(pkstruct, x, y)
    len = length(pkstruct.(x));
    leni = cellfun(@length, pkstruct.(x));
    
    originout = zeros(max(leni), 2*len);
    originout(originout == 0) = NaN;
    
    for i=1:len
        originout(1:leni(i),2*i-1) = pkstruct.(x){i};
       	originout(1:leni(i),2*i) = pkstruct.(y){i};
    end


    %pareval('PKRecoil','().RemnantM(().RemnantM==0)=NaN;')

    %pareval('PKRecoil','().eSFDfit = polyfit(().RemnantM(40:end),().TagawaDeltaH_ext(40:end),1)')
    %pareval('PKRecoil','().eSFDfitvals = polyval(().eSFDfit,().StartingM)')
    %pareval('PKRecoil','().NormDeltaH_ext = (().TagawaDeltaH_ext-(-1-().eSFDfit(2))/().eSFDfit(1))') %remove any offset from 0 at startingM = -1
    %pareval('PKRecoil','().NormDeltaH_ext = ().NormDeltaH_ext./().eSFDfit(1)')  %normalize based on fit

    %RemM_eSFD = extractexcel('PKRecoil','RemnantM','TagawaDeltaH_ext');

    %RemM_eSFD_Norm = extractexcel('PKRecoil','RemnantM','NormDeltaH_ext');
end