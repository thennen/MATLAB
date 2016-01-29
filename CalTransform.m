function DMSStructOut = CalTransform(DMSStructIn,NewCal)
% CALTRANSFORM Calibration transformation for DMS data
% DMSStructIn is the structure containing the DMS data, which must include
% Angle, X, Y, and MeasurementParameters
% MeasurementParameters and NewCal are structures containing fields
% EmuPerVolt, YCC, SSC, Alpha, and Beta, corresponding to the calibration

% Transform back to coil voltages using calibration data
if(~isfield(DMSStructIn,'MeasurementParameters'))
    error('calibration data not present or input is not a DMS struct')
end

DMSStructOut = DMSStructIn;

theta = DMSStructIn.Angle * pi / 180;

for i = 1:length(theta)
    %define rotation matrices
    R{i} = [cos(theta(i)) sin(theta(i)); -sin(theta(i)) cos(theta(i))];
end

EmuPerVolt = DMSStructIn.MeasurementParameters.EmuPerVolt;
YCC = DMSStructIn.MeasurementParameters.YCC;
SSC = DMSStructIn.MeasurementParameters.SSC;
Alpha = DMSStructIn.MeasurementParameters.Alpha * pi/180;
Beta = DMSStructIn.MeasurementParameters.Beta * pi/180;

T = EmuPerVolt .* [cos(Alpha) cos(Beta)*YCC; sin(Alpha)*SSC sin(Beta)*YCC*SSC];

NewEmuPerVolt = NewCal.EmuPerVolt;
NewYCC = NewCal.YCC;
NewSSC = NewCal.SSC;
NewAlpha = NewCal.Alpha * pi/180;
NewBeta = NewCal.Beta * pi/180;

NewT = NewEmuPerVolt .* [cos(NewAlpha) cos(NewBeta)*NewYCC; sin(NewAlpha)*NewSSC sin(NewBeta)*NewYCC*NewSSC];

%rotate to Par and Perp frame
for i = 1:length(theta)
    ParPerp(i,1:2) =  R{i} \ [DMSStructIn.X(i); DMSStructIn.Y(i)];
end
%Par = ParPerp(1,:).';
%Perp = ParPerp(2,:).';

%undo previous calibration, apply new calibration
NewParPerp = (NewT * (T \ ParPerp.')).';

%rotate back to X and Y frame
for i = 1:length(theta)
    NewXY(:,i) =  R{i} * NewParPerp(i,1:2).';
end

NewX = NewXY(1,:).';
NewY = NewXY(2,:).';

DMSStructOut.X = NewX;
DMSStructOut.Y = NewY;
DMSStructOut.MeasurementParameters = NewCal;

end

