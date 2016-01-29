function Out = findperppar( In )
%FINDPERPPAR calculate perpendicular and parallel signal from X and Y in
% DMS data structure.

Out = In;
Out.Perp = In.X .* sin(In.Angle*pi/180) + In.Y.*cos(In.Angle*pi/180);
Out.Par = In.X .* cos(In.Angle*pi/180) - In.Y.*sin(In.Angle*pi/180);
Out.Magnitude = (In.X.^2 + In.Y.^2).^(1/2);

if isfield(In,'X2')
    Out.Perp2 = In.X2 .* sin(In.Angle*pi/180) + In.Y2.*cos(In.Angle*pi/180);
    Out.Par2 = In.X2 .* cos(In.Angle*pi/180) - In.Y2.*sin(In.Angle*pi/180);
    Out.Magnitude2 = (In.X2.^2 + In.Y2.^2).^(1/2);
end

if isfield(In,'RawCorrX')
    Out.RawCorrPerp = In.RawCorrX .* sin(In.Angle*pi/180) + In.RawCorrY.*cos(In.Angle*pi/180);
    Out.RawCorrPar = In.RawCorrX .* cos(In.Angle*pi/180) - In.RawCorrY.*sin(In.Angle*pi/180);
    Out.RawCorrMagnitude = (In.RawCorrX.^2 + In.RawCorrY.^2).^(1/2);
end

if isfield(In,'CorrY')
    Out.CorrPerp = In.CorrX .* sin(In.Angle*pi/180) + In.CorrY.*cos(In.Angle*pi/180);
    Out.CorrPar = In.CorrX .* cos(In.Angle*pi/180) - In.CorrY.*sin(In.Angle*pi/180);
    Out.CorrMagnitude = (In.CorrX.^2 + In.CorrY.^2).^(1/2);
end

if isfield(In,'Xinfield')
    Out.Perpinfield = In.Xinfield .* sin(In.Angle*pi/180) + In.Yinfield.*cos(In.Angle*pi/180);
    Out.Parinfield = In.Xinfield .* cos(In.Angle*pi/180) - In.Yinfield.*sin(In.Angle*pi/180);
end

end

