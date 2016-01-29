function [NormHys, Amplitude, Offset, FitLine1, FitLine2] = branchnormalizehys(Field,HysVals,FitRange)
% Centers and normalizes hysteresis loop
% No slope correction
% FitRange = [+minimum field, +maximum field]

MinField = min(Field);
MaxField = max(Field);

BranchMask = true(length(Field),1).';
BranchMask(floor(length(Field)/2):end) = false;
if Field(1) < 0
   BranchMask = ~BranchMask;                                         	% In case loop starts from negative fields
end

FitMask1 = (Field > FitRange(1)) & (Field < FitRange(2)) & BranchMask;                  %
FitMask2 = Field < -FitRange(1) & Field > -FitRange(2) & ~BranchMask;                   %

if any(FitMask1)
    FitLine1 = polyfit(Field(FitMask1),HysVals(FitMask1),0);            % Fit line for red loop positive saturation
else
    FitLine1 = 0;
end
if any(FitMask2)
    FitLine2 = polyfit(Field(FitMask2),HysVals(FitMask2),0);            % Fit line for red loop negative saturation
else
    FitLine2 = 0;
end

Amplitude = abs((FitLine1(1) - FitLine2(1))/2);                         %
Offset = (FitLine1(1) + FitLine2(1))/2;                                 %
NormHys = (HysVals - Offset) ./ Amplitude;                              % Remove offset and normalize loop


end