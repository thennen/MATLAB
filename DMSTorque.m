function TorqueOut = DMSTorque(TorqueIn,GlassTorque,varargin)
%DMSTORQUE Analyses VSM Torque Data
% DMSTorque(Sample_Torque, Glass_Torque, Sample_PerpHys)
% DMSTorque(Sample_Torque, Glass_Torque, Sample_PerpHys, 'demag',0.4)
% 
% PerpHys used for Ms and for subtracting the glass background, since the 
% glass used for Glass_Torque will in general have a different thickness than 
% the glass substrate of the sample
%
% Sample thickness, area, and Ms can be specified as fields in TorqueIn OR in hysteresis 
% data.  For example, if the following fields are set, Sample1_Torque will 
% find and use them:
% Sample1_PerpHys.Area = 50.1
% Sample1_PerpHys.Thickness = 8
% DMSHys('Sample1_PerpHys')       (to create Sample1_PerpHys.Ms)
%
% The following fields will be added to the data structure:
% ---------------------------------------------------------
% K1amp         Sin(2theta) amplitude of the component of M, perp to H
% CorrX         BG corrected X
% CorrY         BG corrected Y
% MAngle        Angle of the M vector
% K1component   N=2 fourier component, or fit result on CorrY
% K1angle       Phase of N=2 fourier component
%
% If Thickness, Area, and Ms are given, the following fields are also added:
% ----------------------------------------------------------
% K1eff         Effective anisotropy = K1amp*TorqueField/Volume
% K1            Demag corrected anisotropy = K1eff - 2piNMs^2
% Hk            Anisotropy Field = 2*K1/Ms


%% Default options

% Demag factor
options.demag = 1;
% Plot torque curve and fit at the end
options.plot = 0;
% Force Sin2x fit instead of FFT
options.fit = 0;

% Parse input arguments to options
[options, paramopts] = THargparse(varargin, options);

% Set analysis parameters from key fields in structure
% only if they weren't set by input argument
fnames = fieldnames(options);
for i=1:length(fnames)
    % If there is a field with an option as a name AND that option was not
    % set by input parameters
    if isfield(TorqueIn, fnames(i)) && ~any(strcmpi(fnames(i), paramopts))
        options.(fnames{i}) = TorqueIn.(fnames{i});
    end
end

varname = inputname(1);
TorqueOut = TorqueIn;
TorqueOut.Analysis.GlassTemp = GlassTorque.Temp;

Angle = round(TorqueIn.Angle);
TField = TorqueIn.Field(1);
X = TorqueIn.X;
Y = TorqueIn.Y;
% Glass better have the same TField as sample
GlassX = GlassTorque.X;
GlassY = GlassTorque.Y;


%% Do background correction
if nargin > 2 && isstruct(varargin{1})
    % Assume varargin{1} is PerpHys
    PerpHys = varargin{1};
    % Look for Ms, Thickness, and Area in PerpHys data, copy to TorqueOut
    if isfield(PerpHys, 'Area') && isfield(PerpHys, 'Thickness') && isfield(PerpHys, 'Ms')
        TorqueOut.Ms = PerpHys.Ms;
        TorqueOut.Thickness = PerpHys.Thickness;
        TorqueOut.Area = PerpHys.Area;
    end
    
    % If perp hys didn't go to fields as high as torque field, fit line and extrapolate, otherwise interp1
    if max(PerpHys.Field) < TField
        BGLine = polyfit(PerpHys.Field,PerpHys.Background,1);
        BGX = (GlassX - GlassX(Angle == 90) + BGLine(1)*TField + BGLine(2));
    else
        FieldLength = length(PerpHys.Field);
        % Silly way to prevent interp1 from complaining about distinct x values.
        PerpHys.Field(1:floor(FieldLength/2)) = PerpHys.Field(1:floor(FieldLength/2)) + 0.0001;
        GlassX = interp1(GlassTorque.Angle,GlassX,Angle,'linear','extrap');
        % To calculate X background, start with GlassX, adding diamagnetic constant
        % determined by PerpHys
        BGX = (GlassX - GlassX(Angle == 90) + interp1(PerpHys.Field,PerpHys.Background,TField,'linear','extrap'));
    end
else
    % If there is no perphys, use glass torque directly as glass background of sample
    % This will result in some error when the glass piece measured in GlassTorque has
    % a different thickness than the glass substrate on the sample.
    BGX = interp1(GlassTorque.Angle, GlassX, Angle,'linear','extrap');
end
CorrX = X - BGX;
% Diamagnetic component doesn't matter for Y, just subtract GlassY
BGY = interp1(GlassTorque.Angle, GlassY, Angle,'linear','extrap');
CorrY = Y - BGY;


%% Do FFT to determine amplitude of the sin(2theta) term
% Need at least one cycle to do an FFT, meaning a field angle range of at least 180 deg.
% Use 360 degrees of data if it's there, else use 180 degrees of data if it's there
% else try to fit a sin(2theta) function to the data, which requires license checkout.
MAngle = Angle + atan(CorrY./CorrX).*180/pi;
minAngle = min(MAngle);
maxAngle = max(MAngle);
diffAngle = max(MAngle) - min(MAngle);
if ~options.fit && diffAngle >= 360
    % Prepare arrays for FFT 
    FFTangle = linspace(min(Angle),min(Angle)+360, 256);
    FFTY = interp1(MAngle,CorrY,FFTangle,'linear','extrap');
    FFT = fft(FFTY);
    K1amp = abs(FFT(3))/128;
    K1angle = angle(FFT(3)) * 180 / pi;
    K1component = FFT(1)/256 + FFT(3)/128*exp(2*pi*1i*2*MAngle/360);
    K1component = real(K1component);
elseif ~options.fit && diffAngle >=180
    FFTangle = linspace(min(Angle),min(Angle)+180, 128);
    FFTY = interp1(MAngle,CorrY,FFTangle,'linear','extrap');
    FFT = fft(FFTY);
    K1amp = abs(FFT(2))/64;
    K1angle = angle(FFT(2)) * 180 / pi;
    K1component = FFT(1)/128 + FFT(2)/64*exp(2*pi*1i*MAngle/180);
    K1component = real(K1component);
else
    % Try to fit sin2 function
    % phase currently locked to 0 degrees
    K1amp = lsqcurvefit(@sin2, [50], MAngle, CorrY);
    K1component = sin2(K1amp, MAngle);
    K1angle = 0;
end

%Japan method
%JapanRange = (Angle < 150 & Angle > 30) | (Angle < 330 & Angle > 210);
%JapanK1amp = lsqcurvefit(@sin2,[.0001],MAngle(JapanRange),CorrY(JapanRange));
%JapanK1component = sin2(JapanK1amp,MAngle);
%TorqueOut.JapanK1amp = JapanK1amp(1);
%TorqueOut.JapanK1component = JapanK1component;

TorqueOut.K1amp = K1amp;
TorqueOut.CorrX = CorrX;
TorqueOut.CorrY = CorrY;
TorqueOut.MAngle = MAngle;
TorqueOut.K1component = K1component;
TorqueOut.K1angle = K1angle;
TorqueOut.Analysis.options = options;
TorqueOut.Analysis.Angle = Angle;
TorqueOut.Analysis.BGX = BGX;
TorqueOut.Analysis.BGY = BGY;


% If Area, thickness, and Ms are defined, calculate anisotropy energy/fields
if isfield(TorqueOut,'Area') && isfield(TorqueOut,'Thickness') && isfield(TorqueOut,'Ms')
    Thickness = TorqueOut.Thickness;
    Area = TorqueOut.Thickness;
    Ms = TorqueOut.Ms;
    Kshape = options.demag * Ms^2 * 2 * pi/1000000;
    TorqueOut.K1eff = K1amp*TField/1000/Area/Thickness;
    TorqueOut.K1 = Kshape + K1amp*TField/1000/Area/Thickness;
    TorqueOut.Hk = TorqueOut.K1*2000/Ms;
end

% Make a simple plot if option is set
if options.plot
    figure
    plot(TorqueOut.MAngle, TorqueOut.CorrY, TorqueOut.MAngle, TorqueOut.K1component)
    xlabel('M Angle (degrees)')
    ylabel('Moment Perp to H (uemu)')
    title(strrep(varname, '_', ' '))
end

end

function F = sin2(P,xdata)
%function for fitting Asin(2x+B) to some data.  xdata should be given in
%degrees.
A = P(1);
B = 0;
F = A*sin(2*(xdata+B)*pi/180);
end

