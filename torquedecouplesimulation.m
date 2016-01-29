Ms = [600 400];
Ku = [6 4];
Thickness = [4.5 3.25];

%Ms = [600 250];
%Ku = [6 2];
%Thickness = [8 8];
MsRaw = Ms.*Thickness*50/1000;
MAngle = [0:5:360].';

MsWeightAvg = sum(Ms.*Thickness)/sum(Thickness);
KuWeightAvg = sum(Ku.*Thickness)/sum(Thickness);

TorqueWeightAvg = (KuWeightAvg - 2*pi*MsWeightAvg^2/1000000)*50*sum(Thickness)/20*sin(2*pi*MAngle/180);

YSum = MAngle - MAngle;
XSum = YSum;

YSumavgdemag = MAngle - MAngle;
XSumavgdemag = YSum;

for i = 1:length(Ku)
   Y(:,i) = (Ku(i)-2*pi*Ms(i)^2/1000000)*50*Thickness(i)/20*sin(2*pi*MAngle/180);
   X(:,i) = (MsRaw(i)^2 - Y(:,i).^2).^(1/2);
   
   Yavgdemag(:,i) = (Ku(i)-2*pi*MsWeightAvg^2/1000000)*50*Thickness(i)/20*sin(2*pi*MAngle/180);
   Xavgdemag(:,i) = (MsRaw(i)^2 - Yavgdemag(:,i).^2).^(1/2);
   
   DAngle(:,i) = 180/pi*asin(Y(:,i)/MsRaw(i));
   DAngleavgdemag(:,i) = 180/pi*asin(Yavgdemag(:,i)/MsRaw(i));
   
   
   FAngle(:,i) = MAngle-180/pi*asin(Y(:,i)/MsRaw(i));
   FAngleavgdemag(:,i) = MAngle-180/pi*asin(Yavgdemag(:,i)/MsRaw(i));
   
   InterpY(:,i) = interp1(FAngle(:,i),Y(:,i),0:5:360);
   InterpX(:,i) = interp1(FAngle(:,i),X(:,i),0:5:360);
   
   InterpYavgdemag(:,i) = interp1(FAngleavgdemag(:,i),Yavgdemag(:,i),0:5:360);
   InterpXavgdemag(:,i) = interp1(FAngleavgdemag(:,i),Xavgdemag(:,i),0:5:360);
   
   %TorqueSum = TorqueSum + InterpTorque(:,i);
   YSum = YSum + InterpY(:,i);
   XSum = XSum + InterpX(:,i);
   
   YSumavgdemag = YSumavgdemag + InterpYavgdemag(:,i);
   XSumavgdemag = XSumavgdemag + InterpXavgdemag(:,i);
end

MagSum = (XSum.^2 + YSum.^2).^(1/2);
SumMAngle = 180/pi*atan(YSum./XSum) + MAngle;

MagSumavgdemag = (XSumavgdemag.^2 + YSumavgdemag.^2).^(1/2);
SumMAngleavgdemag = 180/pi*atan(YSumavgdemag./XSumavgdemag) + MAngle;

FFTAngle = [0:5:355].';
FFTY = interp1(SumMAngle,YSum,FFTAngle);
FFT = fft(FFTY);
N2amp = abs(FFT(3))/36;
PFN2amp = -imag(FFT(3))/36; % phase = 0
N2Angle = angle(FFT(3)) * 180 / pi;
N2component = FFT(1)/36 + FFT(3)/36*exp(2*pi*1i*2*FFTAngle/360);
N2component = interp1(FFTAngle,N2component,SumMAngle,'linear','extrap');

KuDecoupled = N2amp*20/50/sum(Thickness) + 2*pi*(sum(Ms.*Thickness)/sum(Thickness))^2/1000000;

FFTYavgdemag = interp1(SumMAngleavgdemag,YSumavgdemag,FFTAngle);
FFTavgdemag = fft(FFTYavgdemag);
N2ampavgdemag = abs(FFTavgdemag(3))/36;
PFN2ampavgdemag = -imag(FFTavgdemag(3))/36; % phase = 0
N2Angleavgdemag = angle(FFTavgdemag(3)) * 180 / pi;
N2componentavgdemag = FFTavgdemag(1)/36 + FFTavgdemag(3)/36*exp(2*pi*1i*2*FFTAngle/360);
N2componentavgdemag = interp1(FFTAngle,N2componentavgdemag,SumMAngleavgdemag,'linear','extrap');

KuDecoupledavgdemag = N2ampavgdemag*20/50/sum(Thickness) + 2*pi*(sum(Ms.*Thickness)/sum(Thickness))^2/1000000;

%plot(MAngle,Y)
%hold all
%plot(MAngle,InterpTorque)
%plot(MAngleOfTorqueSum,TorqueSum,MAngleOfTorqueSum,N2component,MAngle,TorqueWeightAvg)
plot(SumMAngle,YSum,MAngle,TorqueWeightAvg,SumMAngleavgdemag,YSumavgdemag)
