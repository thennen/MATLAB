function [NormHys, Amplitude, Offset, SatLine1, SatLine2, Slope] = normalizehys2(Field,HysVals,SlopeRange,SatRange)
% Centers and normalizes hysteresis loop
% Slope correction

MinField = min(Field);
MaxField = max(Field);


%% Subtract Slope

if length(SlopeRange) == 2
    SlopeFitMask1 = Field > SlopeRange(1) & Field < SlopeRange(2);                  % Mask everything outside of the specified fit field ranges
    SlopeFitMask2 = Field < -SlopeRange(1) & Field > -SlopeRange(2);                %
elseif length(SlopeRange) == 4
    SlopeFitMask1 = Field > SlopeRange(3) & Field < SlopeRange(4);                  % Mask everything outside of the specified fit field ranges
    SlopeFitMask2 = Field > SlopeRange(1) & Field < SlopeRange(2);                  %
end

%BranchMask = true(length(Field),1).';
%BranchMask(floor(length(Field)/2):end) = false;
%if Field(1) < 0
%    BranchMask = ~BranchMask;                                         	% In case loop starts from negative fields
%end

%FitMask1 = (Field > FitRange(1)) & (Field < FitRange(2)) & BranchMask; %
%FitMask2 = Field < -FitRange(1) & Field > -FitRange(2) & ~BranchMask;  %


if any(SlopeFitMask1)
    SlopeLine1 = polyfit(Field(SlopeFitMask1),HysVals(SlopeFitMask1),1);            % Fit line for red loop positive saturation
else
    SlopeLine1 = 0;
end
if any(SlopeFitMask2)
    SlopeLine2 = polyfit(Field(SlopeFitMask2),HysVals(SlopeFitMask2),1);            % Fit line for red loop negative saturation
else
    SlopeLine2 = SlopeLine1;
end

Slope = (SlopeLine2(1) + SlopeLine1(1))/2;
NormHys = HysVals - Slope*Field;

%% Normalize after slope correction

if length(SlopeRange) == 2
    SatFitMask1 = Field > SatRange(1) & Field < SatRange(2);                  % Mask everything outside of the specified fit field ranges
    SatFitMask2 = Field < -SatRange(1) & Field > -SatRange(2);                %
elseif length(SlopeRange) == 4
    SatFitMask1 = Field > SatRange(3) & Field < SatRange(4);                  % Mask everything outside of the specified fit field ranges
    SatFitMask2 = Field > SatRange(1) & Field < SatRange(2);                %
end

SatLine1 = polyfit(Field(SlopeFitMask1),NormHys(SlopeFitMask1),0);
if any(SatFitMask2)
    SatLine2 = polyfit(Field(SlopeFitMask2),NormHys(SlopeFitMask2),0);
else
    SatLine2 = 0;
end

Amplitude = abs((SatLine1(1) - SatLine2(1))/2);                         %
Offset = (SatLine1(1) + SatLine2(1))/2;                                 %
NormHys = (NormHys - Offset) ./ Amplitude;                              % Remove offset and normalize loop


end
