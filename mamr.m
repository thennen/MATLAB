function DataOut = mamr(DataIn)

Field = DataIn.Field;
V = DataIn.V_Hall;
Length = length(Field);

% Assume field sweep starts increasing from zero


if abs(Field(1) - Field(end)) < 2
    % Correct drift by linear subtraction
    CorrV = V + (V(1) - V(end))/(Length-1)*[0:Length-1].';
    
    
    [~,Imax] = max(Field);
    [~,Imin] = min(Field);
    if CorrV(Imax) < CorrV(Imin)
        CorrV = -CorrV; % reorient so positive field gives positive M
    end
   
    
    % Zero field normalization method
    norm1 = CorrV(1);
    norm2 = lininterp1(flipud(Field(Imax:Imin)),flipud(CorrV(Imax:Imin)),0);
    %CorrV = (CorrV - (norm1 + norm2)/2)./abs((norm1-norm2)/2);
    
    % Saturation fit method
    % Need to rearrange hys array to standard sweep so that normalizehys function works
    % properly
    order = [Imax:Length 1:Imax-1];
    iorder = 1:Length;
    iorder(order) = iorder;
    orderField = Field(order);
    orderV = CorrV(order);
    
    FitRange = FindFitRange(DataIn);
    orderCorrV = normalizehys(orderField,orderV,FitRange);
    
    CorrV = orderCorrV(iorder);
    
    
    % Find Hc 
    Hc1 = lininterp1(CorrV(1:Imax),Field(1:Imax),0);
    Hc2 = lininterp1(flipud(CorrV(Imax:Imin)),flipud(Field(Imax:Imin)),0);
    Hc = (Hc1-Hc2)/2;
    
    % Find Hn  (only positive values for now)
    %HnDef = 0.9*(norm1+norm2)/2; % 0.9 of remnant value
    HnDef = 0.8;
   
    [~,Izero] = min(abs(Field(Imax:Imin))); % Find index of array where field has returned to its nearest to zero value
    Izero = Izero + Imax -1;  % Izero was determined with Imax offset
    Hn1 = lininterp1(CorrV(1:Imax),Field(1:Imax),-HnDef);
    Hn2 = lininterp1(flipud(CorrV(Izero:Imin)),flipud(Field(Izero:Imin)),HnDef);
    Hn = (Hn1-Hn2)/2;
    
    % Find Hs
    HsDef = 0.8;
    Hs1 = lininterp1(CorrV(1:Imax),Field(1:Imax),HsDef);
    Hs2 = lininterp1(flipud(CorrV(Izero:Imin)),flipud(Field(Izero:Imin)),-HsDef);
    Hs = (Hs1-Hs2)/2;
    
    
    
    
    DataOut = DataIn;
    DataOut.CorrV = CorrV;
    DataOut.Hc = Hc;
    DataOut.Hc1 = Hc1;
    DataOut.Hc2 = Hc2;
    DataOut.Hn1 = Hn1;
    DataOut.Hn2 = Hn2;
    DataOut.Hn = Hn;
    DataOut.Hs = Hs;
else
    warning('Starting and ending fields differ by %.1f',Field(1)-Field(end))
    CorrV = V;
    
    DataOut = DataIn;
end



end

function FitRange = FindFitRange(DataIn)
% Searches DataIn and base workspace for a variable named FitRange =
% [FitRange(1) FitRange(2)], and if found returns that variable's value,
% otherwise returns a default value.
    MaxField = max(DataIn.Field);
    FitRangeSearch = {'FitRange','fitrange','Fitrange'};
    if sum(isfield(DataIn,FitRangeSearch))
        FitRangeName = FitRangeSearch(isfield(DataIn,FitRangeSearch) == 1);
        eval(['FitRange = DataIn.' FitRangeName{1} ';']);
    else
        string = 'who (''-regexp'' ,''\w*(F|f)(it|IT)(R|r)(ANGE|ange)\w*'')';
        FitRangeVar = evalin('base', string);                                       % Looks for variable in base workspace which specifies the fit range to be used
        if ~isempty(FitRangeVar)                                                    % Should have form FitRange = [FitRange(1) FitRange(2)]
            MaybeFitRange = evalin('base', FitRangeVar{1});
            if length(MaybeFitRange) == 2 && MaybeFitRange(1) >= -5000 && MaybeFitRange(1) < 40000 && MaybeFitRange(2) > 0 && MaybeFitRange(2) < 40000 && MaybeFitRange(1) < MaybeFitRange(2)
                FitRange = MaybeFitRange;
            else
                disp('There is a fit range variable with an unexpected format. Using default values.')
                FitRange = [0.75*MaxField MaxField];
            end
        else
            FitRange = [0.75*MaxField MaxField];                                     % If FitRange is not found, use this by default
        end
    end
end

