function DataOut = DMSDCD(DataIn,varargin)
% DMSDCD  Analyze dms DCD data.
% Determines saturation levels of DCD curve by fitting parabolas to the first 
% and last fitpercent of the curve.  By default, the curve offset is determined 
% by these fits, but can also be specified directly as an input argument.  
% Specifying the offset directly is useful when you want to ensure the same 
% value is used for multiple DCD curves, as in angular dependence measurements.
%
% Sample1_DCD = DMSDCD(Sample1_DCD) Analyzes DCD data for single data structure
% Sample1_DCD = DMSDCD(Sample1_DCD, 0.31) Offset specified by input argument
%
% Analysis produces the following fields:
%
% Ms         Amplitude of DCD curve
% Offset     Offset used in correction
% CorrPerp   Corrected DCD curve, in perpendicular direction
% Fangle     Angle of applied field wrt perp axis
% Fitlines   Structure containing fit details
% Hn
% Hc

if nargin > 2
    error('Too many arguments.');
end

fitpercent = 0.15;
field = DataIn.Field;
maxfield = max(field);
% Calculate signal in perpendicular direction, from X and Y 
perp = DataIn.X .* sin(DataIn.Angle*pi/180) + DataIn.Y.*cos(DataIn.Angle*pi/180);

%% Fit parabolas to ends of curve
fitmask1 = field <= maxfield*fitpercent;
fitmask2 = field <= maxfield & field >= maxfield*(1-fitpercent);
% Temporary override
%fitmask1 = field<10;
%fitmask2 = field>400;
fit1 = polyfit(field(fitmask1),perp(fitmask1),2);
fit2 = polyfit(field(fitmask2),perp(fitmask2),0);
M1 = fit1(3);
M2 = fit2;
Ms = (M2-M1)/2;

if nargin == 1
    offset = (M2+M1)/2;
else
    % This is for inputting offset directly
    offset = varargin{1};
end
corrperp = perp - offset;

fitmask3 = corrperp <= Ms*.25 & corrperp >= -Ms*.25;

% If there are fewer than two data points in hc fit range
if sum(fitmask3) < 2
% Find points nearest 0
fitmask3(find(corrperp>0,2,'first')) = true;
fitmask3(find(corrperp<0,2,'last')) = true;
end

    fit3 = polyfit(field(fitmask3),corrperp(fitmask3),1);
    Hc = -fit3(2)/fit3(1);
    % Just use the nearest points surrounding M=0
    % !! overrides previous Hc calculation !!
    % Hc = lininterp1(corrperp,field,0)
    % Hn = interp1(corrperp,field,-Ms*0.90);
    
DataOut = DataIn;
DataOut.Ms = Ms;
DataOut.Offset = offset;
DataOut.Hn = Hn;
DataOut.CorrPerp = corrperp;
DataOut.Hc = Hc;
DataOut.Fitlines.RawMsLine1(1:length(field)) = Ms;
DataOut.Fitlines.RawMsLine2 = -DataOut.Fitlines.RawMsLine1;
DataOut.Fitlines.RawHcLine = fit3(2) + fit3(1)*field;

DataOut.FitLines(1:length(field),1) = Ms;
DataOut.FitLines(1:length(field),2) = -Ms;
DataOut.FitLines(1:length(field),3) = fit3(2) + fit3(1)*field;
DataOut.FitLines(1:length(field),4) = 0;

DataOut.Fangle = 90 - round(DataIn.Angle(1));
end

