%Trying to correct SUL signal in minor loops

SlopeRange = [2000,6000];

BranchAvgField = (Field{2} - Field{1})/2;
BranchAvgKerr = (Kerr{2}-Kerr{1})/2;

%FitMask1 = Field{1} >= SlopeRange(1) & Field{1} <- SlopeRange(2);
FitMask = BranchAvgField >= -SlopeRange(2) & BranchAvgField <= -SlopeRange(1);

%FitSlope1 = polyfit(Field{1}(FitMask1),Kerr{1}(FitMask1),1);
FitSlope = polyfit(BranchAvgField(FitMask),BranchAvgKerr(FitMask),1);

% extrapolate sul signal for field <  sloperange(1)
region1 = BranchAvgField < 0 & BranchAvgField > -SlopeRange(1);   %  extrap region
region2 = BranchAvgField <= -SlopeRange(1);                      % rest of neg region

SULsignal = [BranchAvgKerr(region2), polyval(FitSlope,BranchAvgField(region1))];
%remove offset
SULsignal = SULsignal - FitSlope(2);
%mirror
SULsignalField = [BranchAvgField(region1 | region2), -fliplr(BranchAvgField(region1 | region2))];
SULsignal = [SULsignal, -fliplr(SULsignal)];

%subtract SUL signal from all loops

for i=1:length(Field)
    Kerr{i} = Kerr{i} - interp1(SULsignalField,SULsignal,Field{i},'linear','extrap');
end

%something else

for i=1:6
    H8_PKDynamic.KerrTest(:,i) = H8_PKDynamic.Kerr(:,i) - interp1(H8_PKRecoil.SULsignalField,H8_PKRecoil.SULsignal,H8_PKDynamic.Field(:,i),'linear','extrap');
end
    