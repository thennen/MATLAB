function DataOut = NormalizeRemDecay(DataIn)

% MaxField = max(DataIn.Field{1}(:,1));
FitRange = [0 6000];

FitMask1 = DataIn.Field{1} > FitRange(1) & DataIn.Field{1} < FitRange(2);
FitMask2 = DataIn.Field{1} < -FitRange(1) & DataIn.Field{1} > -FitRange(2);

FitLine1 = polyfit(DataIn.Field{1}(FitMask1),DataIn.Kerr_x{1}(FitMask1),0);
FitLine2 = polyfit(DataIn.Field{1}(FitMask2),DataIn.Kerr_x{1}(FitMask2),0);

% Define and initialize DataOut vector
DataOut = [0 0];

% Calculate amplitude
DataOut(1) = abs((FitLine1(1) - FitLine2(1))/2);

% Calculate offset
DataOut(2) = (FitLine1(1) + FitLine2(1))/2;

end