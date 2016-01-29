function TorqueOut = DMSTorque(TorqueIn,GlassTorque,varargin)
%UNTITLED Summary of this function goes here

if nargin > 3
    error('too many arguments.');
end
TorqueOut = TorqueIn;

TorqueOut.Analysis.GlassTemp = GlassTorque.Temp;


Angle = round(TorqueIn.Angle);
TField = TorqueIn.Field(1);
X = TorqueIn.X;
Y = TorqueIn.Y;
GlassX = GlassTorque.X; % glass better have the same field angles and TField as sample
GlassY = GlassTorque.Y;
if nargin == 2
    %CorrX = X - (GlassX - GlassX(Angle == 90) - 7.2622e-08*TField - 8.9258e-06);    % If there is no perphys, use nominal glass background
    %CorrX = X - GlassX;
    CorrX = X - interp1(GlassTorque.Angle, GlassX, Angle,'linear','extrap');
else
    PerpHys = varargin{1};
    if max(PerpHys.Field) < TField                                          % If perp hys didn't go to fields as high as torque field, fit line and extrapolate, otherwise interp1
        BGLine = polyfit(PerpHys.Field,PerpHys.Background,1);
        CorrX = X - (GlassX - GlassX(Angle == 90) + BGLine(1)*TField + BGLine(2));
    else
        FieldLength = length(PerpHys.Field);                                                    %
        PerpHys.Field(1:floor(FieldLength/2)) = PerpHys.Field(1:floor(FieldLength/2)) + 0.001; % silly way to prevent interp1 from complaining about distinct x values.
        GlassX = interp1(GlassTorque.Angle,GlassX,Angle,'linear','extrap');
        CorrX = X - (GlassX - GlassX(Angle == 90) + interp1(PerpHys.Field,PerpHys.Background,TField,'linear','extrap'));
    end
    
    if isfield(PerpHys,'Ms')
        TorqueOut.Ms = PerpHys.Ms;
        TorqueOut.Thickness = PerpHys.Thickness;
        TorqueOut.Area = PerpHys.Area;
        
        Kshape = PerpHys.Ms^2*2*pi/1000000;
    end
end
%CorrY = Y - GlassY;
CorrY = Y - interp1(GlassTorque.Angle, GlassY, Angle,'linear','extrap');
CorrAngle = Angle + atan(CorrY./CorrX).*180/pi;
%FFTangle = [0:5:355].';
FFTangle = [0:5:175].';
FFTY = interp1(CorrAngle,CorrY,FFTangle,'linear','extrap');
FFT = fft(FFTY);

N2amp = abs(FFT(2))/18;
PFN2amp = -imag(FFT(2))/18; % phase = 0
N2angle = angle(FFT(2)) * 180 / pi;
N2component = FFT(1)/18 + FFT(2)/18*exp(2*pi*1i*2*FFTangle/360);
N2component = interp1(FFTangle,N2component,CorrAngle,'linear','extrap');

%Japan method
%JapanRange = (Angle < 150 & Angle > 30) | (Angle < 330 & Angle > 210);
%JapanN2amp = lsqcurvefit(@sin2,[.0001],CorrAngle(JapanRange),CorrY(JapanRange));
%JapanN2amp = lsqcurvefit(@sin2,[.0001 0],CorrAngle,CorrY);
%JapanN2component = sin2(JapanN2amp,CorrAngle);
%TorqueOut.JapanN2amp = JapanN2amp(1);
%TorqueOut.JapanN2component = JapanN2component;


TorqueOut.N2amp = N2amp;
TorqueOut.PFN2amp = PFN2amp;
TorqueOut.CorrX = CorrX;
TorqueOut.CorrY = CorrY;
TorqueOut.CorrAngle = CorrAngle;
TorqueOut.N2component = N2component;
TorqueOut.N2angle = N2angle;


if isfield(TorqueOut,'Area') && isfield(TorqueOut,'Thickness') && isfield(TorqueOut,'Ms') % may cause problems at some point
    TorqueOut.K1p = N2amp*TField/1000/TorqueOut.Area/TorqueOut.Thickness;
    TorqueOut.K1 = Kshape + N2amp*TField/1000/TorqueOut.Area/TorqueOut.Thickness;
    TorqueOut.Hk = TorqueOut.K1*2000/TorqueOut.Ms;
end


%for n=0:16
%FFTPlot = FFTPlot + FFT(n+1)/36*exp(2*pi*i*n*FFTangle/360)
%end

end

