function  DataOut = RemDecayAnalysis_v2(MajorLoops,RemDecay,~,Info)

%{
***********************************************************************
*     This function plots:                                            *
*       1) remanent decay (PK at H = 0 versus time)                   *
*       2) viscosity S versus M0                                      *
*       3) viscosity S versus applied field Happ                      *
*     PK data are centered around the x-axis and normalized using     *
*     values from the first major loop. S and M0 are calculated       *
*     from the fit of normalized PK vs log(t).                        *
*                                                                     *
***********************************************************************
%}


% MAJOR LOOPS PLOT and IMAGE SAVING
PlotMajorLoops(MajorLoops,Info);



% NORMALIZATION
% (3 possible ways of doing it: uncomment the desired one)

% Define and initialize vector to store parameters for normalization
NormFactors = [0 0];

% (1) Using first half major loop, calculates amplitude and offset values needed for
% normalization BY FITTING A HORIZONTAL LINE IN SATURATION REGIME (HIGH FIELD)
%NormFactors = NormalizeRemDecay(MajorLoops);

% (2) Using first half major loop, calculates amplitude and offset values needed for
% normalization BY USING MAX AND MIN PK VALUES
%NormFactors(1) = (max(MajorLoops.Kerr_x{1}) - min(MajorLoops.Kerr_x{2}))/2;     % Amplitude (half peak-to-peak)
%NormFactors(2) = (max(MajorLoops.Kerr_x{1}) + min(MajorLoops.Kerr_x{2}))/2;     % Offset

% (3) Using first half major loop, calculates amplitude and offset values needed for
% normalization BY CONSIDERING VALUES AT H = 0
NormFactors(1) = (MajorLoops.Kerr_x{1}(185) - MajorLoops.Kerr_x{2}(185))/2;      % Amplitude (half peak-to-peak)
NormFactors(2) = (MajorLoops.Kerr_x{1}(185) + MajorLoops.Kerr_x{2}(185))/2;      % Offset



% Store number of columns (i.e., number of cells) of cell array RemDecay.Kerr_x
RemDecay_Norm.size = size(RemDecay.Kerr_x,2);

% Set starting time close to zero for each Remanent Decay measurement
% and normalize PK values using parameters estimated from major loop
for i = 1:RemDecay_Norm.size
    RemDecay_Norm.Time{i} = RemDecay.Time{i} - RemDecay.Time{i}(1) + 0.2;
    RemDecay_Norm.Kerr_x{i} = (RemDecay.Kerr_x{i} - NormFactors(2)) ./ NormFactors(1);
end    

% Store normalized data for Remanent Decay in workspace
DataOut.RemDecay_Norm = RemDecay_Norm;




% Plot remanent magnetization as a function of time
figure;
for i = 1:RemDecay_Norm.size
  plot(log(RemDecay_Norm.Time{i}),RemDecay_Norm.Kerr_x{i},'Linewidth',1);
  hold all;
end
xlabel('log(t)');
ylabel('PK');



% Define and initialize cell array to store output from fits
FitLine_masked = cell(1,RemDecay_Norm.size);

% Linear fit of Remanent Decay data vs Log(time) data in the range 10 to 300s
for i = 1:RemDecay_Norm.size
    FitMask = RemDecay_Norm.Time{i} > 10 & RemDecay_Norm.Time{i} < 300;
    FitLine_masked{i} = polyfit(log(RemDecay_Norm.Time{i}(FitMask)),RemDecay_Norm.Kerr_x{i}(FitMask),1);
    clear FitMask;
end

% Create and initialize M0 and S variables
M0 = zeros(RemDecay_Norm.size,1);
S = zeros(RemDecay_Norm.size,1);

% From masked fit (i.e., only in a specific range of times), assign M0 (y-axis intercept)
% and Viscosity (slope) values to the respective one-column vectors
% Also apply viscosity scaling correction (230x)
for i = 1:RemDecay_Norm.size
  M0(i,1) = FitLine_masked{i}(1,2);
  S(i,1) = 230*FitLine_masked{i}(1,1);
end

% Define line vector of applied fields and makes it a column vector
% It would be better to get the vector from the input file, as it is it needs
% to be changed when recipe changes
Happ = [0 -250 -500 -1000:-1000:-9000]';



% Plot and save Viscosity vs Applied Field graph
h3 = figure;
plot(Happ,S,'o','MarkerSize',4,'Linewidth',2);
xlabel('H_{app}');
ylabel('S (% remanent decay / dec)');
print(h3,'-djpeg','-r300',[Info.DiskID '_RemDec_vs_Happ']);

Xinterp = -1 : 0.01 : 1;
Yinterp = interp1(M0,S,Xinterp,'linear','extrap');


% Plot and save Viscosity vs M0 graph
h4 = figure;
plot(M0,S,'o','MarkerSize',4,'Linewidth',2);
hold all;
plot(Xinterp,Yinterp);
xlim([-1 1]);
xlabel('M_0');
ylabel('S (% remanent decay / dec)');
text('Position',[0.9,Yinterp(191)],'HorizontalAlignment','Right','String','S at Ms = 0.9 \rightarrow');
print(h4,'-djpeg','-r300',[Info.DiskID '_RemDec_vs_M0']);

% Yinterp(191) is the PK value at Ms = 0.9
DataOut.Ms_09 = Yinterp(191);

end