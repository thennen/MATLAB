function DataOut = pkrecoil(DataIn , varargin)
% This function analyses PK recoil data.  DataIn should be a struct with
% cell array fields 'Field', 'Time','Kerr_x', and 'Kerr_y', with each cell of
% the cell arrays corresponding to a field sweep measurement starting from a different reverse field.
% DataOut is a struct which contains all the fields in DataIn, with
% additional fields resulting from the analysis routine.
% 
% Version 2/12/12
%
% Field Name        Description
% ----------        -----------
% Kerr
% Hn
% Hc
% Hs
% SFD
% RemnantM
% StartingM
% SFD_ext
% SFD_int
% Nd
% Cluster
% H_Interp 
% M_Interp
% DeltaH_int
% DeltaH_ext
% TagawaDeltaH_ext
% TagawaDeltaH_int
% iField

% TODO: add option to correct for amplitude drift by comparing last minor
% loop with major loop

DataOut = DataIn;

RedOrBlue = 1; % Choose red or blue laser data for analysis. 0 for red, 1 for blue
kink = false;
SlopeCorr = false;
fast = false;
makeplot = false;
drift = false;

if nargin > 1
    for i=1:(nargin-1)
        argi = varargin{i};
        if strcmpi(argi,'Red')
            RedOrBlue = 0;
        elseif strcmpi(argi,'Blue')
            RedOrBlue = 1;
        elseif strcmpi(argi,'fast')
            fast = true;
        elseif strcmpi(argi,'kink')
            %kink correction (SUL signal subtraction)
            kink = true;
        elseif strcmpi(argi,'Slope')
            SlopeCorr = true;
        elseif strcmpi(argi,'drift')
            drift = true;
        elseif strcmp(argi,'plot')
            makeplot = true;
        end
    end
end

% Determines fit range for fitting saturation levels
FitRange = FindFitRange(DataIn);

if ~isstruct(DataIn)
    error('Input is not a struct.')
end

Field = DataIn.Field;
nloops = length(Field);

if RedOrBlue == 0
    Kerr = DataIn.Kerr_y;                             
elseif RedOrBlue == 1
    Kerr = DataIn.Kerr_x;
end

DataOut.Analysis.FitRange = FitRange;
DataOut.Analysis.varargin = varargin;

%% Center and normalize loops (based on major loop)
[a ,Ampl1, Offs1, UOffs1] = normalizehys(Field{1},Kerr{1},FitRange);
[a ,Ampl2, Offs2, UOffs2] = normalizehys(Field{2},Kerr{2},FitRange);
Ampl = (Ampl1 + Ampl2)/2;
Offs = (Offs1 + Offs2)/2;
UOffs = (UOffs1 + UOffs2)/2;
Ampldiscrep = abs(Ampl1-Ampl2)*100/Ampl;
%Offsdiscrep = abs(Offs1-Offs2)*100/Ampl;
Kerr{1} = (Kerr{1} - Offs) ./ Ampl;
Kerr{2} = (Kerr{2} - Offs) ./ Ampl;

if  Ampldiscrep > 1
    warning('Amplitude of the two major loop sweeps differ by %.2f%%',Ampldiscrep)
    disp('Using second sweep for loop normalization')
    Ampl = Ampl2;
    Offs = Offs2;
    UOffs = UOffs2;
    % Could also try to correct the first sweep, by scaling or by subtracting a
    % linear drift
end
 
for i = 3:nloops
    % Use the Ampl determined from major loop to scale each minor loop, but
    % determine offset independently
    [a,a,a,UOffsi] = normalizehys(Field{i},Kerr{i},FitRange);
    Kerr{i} = (Kerr{i} - UOffsi + UOffs - Offs) ./ Ampl;
end

%% Optional Slope Correction
if SlopeCorr == 1;
    [a ,Ampl1, Offs1, UOffs1, a,Slope1] = normalizehys2(Field{1},Kerr{1},FitRange,FitRange);
    [a ,Ampl2, Offs2, UOffs2, a,Slope2] = normalizehys2(Field{2},Kerr{2},FitRange,FitRange);
    Ampl = (Ampl1 + Ampl2)/2;
    Offs = (Offs1 + Offs2)/2;
    Slope = (Slope1 + Slope2)/2;
    UOffs = (UOffs1 + UOffs2)/2;
    Kerr{1} = (Kerr{1}-Slope*Field{1} - Offs) ./ Ampl; %
    Kerr{2} = (Kerr{2}-Slope*Field{2} - Offs) ./ Ampl; % do the same thing to both major loop branches
    
    for i = 3:nloops
        Kerr{i} = (Kerr{i}-Slope*Field{i})./Ampl; % use major loop normalization
        
        [a,a,a,UOffsi,a] = normalizehys(Field{i},Kerr{i},FitRange); % Align saturation regions independently
        Kerr{i} = Kerr{i} - UOffsi + 1;
        %Kerr{i} = Kerr{i} - Slopei*Field{i};
        %Kerr{i} = (Kerr{i} - Slopei*Field{i} - UOffsi + UOffs);
        %Kerr{i} = Kerr{i} ./Ampl;
        %Kerr{i} = (Kerr{i}-Slopei*Field{i} - UOffsi + UOffs - Offs);
        %Kerr{i} = (Kerr{i}-Slopei*Field{i} - UOffsi + UOffs - Offs) ./ Ampl;
    end
end

%% Kink correction (SUL signal subtraction)
if kink == 1;
    SlopeRange = [2000,5000];

    BranchAvgField = (Field{2} - Field{1})/2;
    BranchAvgKerr = (Kerr{2}-Kerr{1})/2;
    FitMask = BranchAvgField >= -SlopeRange(2) & BranchAvgField <= -SlopeRange(1);
    FitSlope = polyfit(BranchAvgField(FitMask),BranchAvgKerr(FitMask),1);

    % extrapolate sul signal for field <  sloperange(1)
    region1 = BranchAvgField < 0 & BranchAvgField > -SlopeRange(1);   %  extrap region
    region2 = BranchAvgField <= -SlopeRange(1);                      % rest of neg region

    SULsignal = [BranchAvgKerr(region2), polyval(FitSlope,BranchAvgField(region1))];
    
    %remove offset
    SULsignal = SULsignal - FitSlope(2);
    %mirror
    SULfield = [BranchAvgField(region1 | region2), -fliplr(BranchAvgField(region1 | region2))];
    SULsignal = [SULsignal, -fliplr(SULsignal)];
    
    DataOut.Analysis.SULfield = SULfield;
    DataOut.Analysis.SULsignal = SULsignal;
    
    %subtract SUL signal from all loops
    for i=1:nloops
        Kerr{i} = Kerr{i} - interp1(SULfield,SULsignal,Field{i},'linear','extrap');
    end
    
    % renormalize (there's probably a better way to do this)
    [a ,Ampl1, Offs1, UOffs1] = normalizehys(Field{1},Kerr{1},FitRange);
    [a ,Ampl2, Offs2, UOffs2] = normalizehys(Field{2},Kerr{2},FitRange);
    Ampl = (Ampl1 + Ampl2)/2;
    Offs = (Offs1 + Offs2)/2;
    UOffs = (UOffs1 + UOffs2)/2;
    %Ampldiscrep = abs(Ampl1-Ampl2)*100/Ampl;
    %Offsdiscrep = abs(Offs1-Offs2)*100/Ampl;
    
    for i=1:nloops
        Kerr{i} = (Kerr{i} - Offs) ./ Ampl;
    end
    
end

%% Subtract linear drift
if drift
    % Currently no flexibility for different order minor loops
    FitMask = Field{2} >= -FitRange(2) & Field{2} <= -FitRange(1);
    DriftFit1 = polyfit(Field{2}(FitMask),Kerr{2}(FitMask),0);
    FitMask = Field{end} >= -FitRange(2) & Field{end} <= -FitRange(1);
    DriftFit2 = polyfit(Field{end}(FitMask),Kerr{end}(FitMask),0);
    
    DriftTotal = DriftFit2(1) - DriftFit1(1);
    Driftinc = DriftTotal / (nloops-2);
    
    % Subtract drift, constant value per minor loop
    for i = 3:nloops
        Kerr{i} = Kerr{i} - (i-2)*Driftinc;
    end
end


%% Find loop parameters (Hn, Hc, Hs) 
HnDef = 0.9;
HsDef = 0.9;

% Put major loop arrays in the right order for interpolation
[MajField1, order1] = sort(Field{1});
MajLoop1 = Kerr{1}(order1);
[MajField2, order2] = sort(Field{2});
MajLoop2 = Kerr{2}(order2);

Hc1 = lininterp1(MajLoop1,MajField1,0);
Hc2 = lininterp1(MajLoop2,MajField2,0);
Hc = abs(-Hc1 + Hc2)/2;

% Determine if MajLoop1 or MajLoop2 is the rightloop
% Free to measure loops in either order
if Hc1 > Hc2
    RLoopField = MajField1;
    RLoop = MajLoop1;
    LLoopField = MajField2;
    LLoop = MajLoop2;
else
    RLoopField = MajField2;
    RLoop = MajLoop2;
    LLoopField = MajField1;
    LLoop = MajLoop1;
end

% Threshold method for Hn, Hs
Hn1 = lininterp1(LLoop, LLoopField, HnDef);
Hn2 = lininterp1(RLoop, RLoopField, -HnDef);
Hn = (-Hn1 + Hn2) / 2;

Hs1 = lininterp1(LLoop, LLoopField, -HsDef);
Hs2 = lininterp1(RLoop, RLoopField, HsDef);
Hs = (-Hs1 + Hs2) / 2;

SFD1 = lininterp1(LLoop, LLoopField,0.5) - lininterp1(LLoop, LLoopField,-0.5);
SFD2 = lininterp1(RLoop, RLoopField,0.5) - lininterp1(RLoop, RLoopField,-0.5);
SFD = (SFD1 + SFD2)/2;

%% Extract remnant values, initial values
Fieldmask = cellfun(@(a)a<1000 & a>-1000,Field,'UniformOutput',0);             %create logical masks for fitting remnant curves near zero field
CrossesZero = cellfun(@(b) sum(b) > 2,Fieldmask);                            %logical array indicating whether or not each loop has enough field values near zero for a parabola fit

for i=1:nloops
        StartingM(i) = Kerr{i}(1);
    if CrossesZero(i) == 1    
        Fieldi = Field{i};
        Loopi = Kerr{i};
        Fieldmaski = Fieldmask{i};
        ParaFit = polyfit(Fieldi(Fieldmaski),Loopi(Fieldmaski),2);                 %fit parabola to fieldmask region
        RemnantM(i) = ParaFit(3);
    else
        %RemnantM(i) = NaN;
        RemnantM(i) = 0;
    end
end

%% Find H(M=1/2) for each minor loop  (Mr < Ms/2)

if fast
    M_Interp = linspace(-.98,.98,99).'; %M levels to extract H(M) for each curve
else
    M_Interp = linspace(-.99,.99,397).';  %M levels to extract H(M) for each curve
end

MTagawa = 0.5;                        %M level for Tagawa analysis
if ~any(M_Interp == MTagawa)
   error('M_Interp must contain the Tagawa level') 
end
if ~all(M_Interp + flipud(M_Interp) <= eps)
    error('M_Interp must be symmetric around zero')
end

for k = 1:length(M_Interp)
    H_Interp(k,1) = lininterp1(MajLoop1,MajField1,M_Interp(k));
    H_Interp(k,2) = lininterp1(MajLoop2,MajField2,M_Interp(k));
end 

% Identify right loop for H interp
if Hc1 > Hc2
    H_InterpRight = H_Interp(:,1);
else
    H_InterpRight = H_Interp(:,2);
end

for i = 3:nloops
    for k = 1:length(M_Interp)
        if ~fast && StartingM(i) < M_Interp(k)
            H_Interp(k,i) = lininterp1(Kerr{i},Field{i},M_Interp(k));
        elseif fast && abs(StartingM(i)) < 0.1 && StartingM(i) < M_Interp(k)
            % fast mode only does interp for curves near M = 0
            H_Interp(k,i) = lininterp1(Kerr{i},Field{i},M_Interp(k));
        else
        	H_Interp(k,i) = NaN;
        end
    end
    
    HalfReturnVal(i) = (1 + StartingM(i))/2;
    %HalfReturnVal(i) = (1 + RemnantM(i))/2;
    
    H_HalfReturned(i) = lininterp1(Kerr{i},Field{i},HalfReturnVal(i));
    H_MajorHalfReturned(i) = lininterp1(RLoop, RLoopField,-HalfReturnVal(i));
    
    DeltaH_int(:,i) = H_InterpRight - H_Interp(:,i);
    DeltaH_ext(:,i) = H_Interp(:,i) - flipud(H_InterpRight); % is this even a thing?
    
    
end

TagawaDeltaH_ext = H_HalfReturned - H_MajorHalfReturned;
TagawaDeltaH_int = DeltaH_int(find(M_Interp == MTagawa,1),:);
TagawaDeltaH_ext(1:2) = NaN;
TagawaDeltaH_int(1:2) = NaN;

%% Interpolate for entire Mr = 0 loop
%MinorBelow = find(RemnantM(3:end)<0,1,'first')+ 2;
%MinorAbove = find(RemnantM(3:end)>0,1,'last') + 2;
%WeightBelow = -RemnantM(MinorBelow);
%WeightAbove = RemnantM(MinorAbove);
%HcrMinorLoopH = (H_Interp(:,MinorBelow).*WeightAbove + H_Interp(:,MinorAbove).*WeightBelow)./(WeightAbove+WeightBelow);  % Not sure if this is working properly

%% iSFD and Hd for Mr = 0 loop

%create arrays for interpolation, sorted for interpolation function
[StartM_interp, order] = sort(StartingM(3:end));
TagawaDeltaH_ext_interp = TagawaDeltaH_ext(3:end);
TagawaDeltaH_ext_interp = TagawaDeltaH_ext_interp(order);
TagawaDeltaH_int_interp = TagawaDeltaH_int(3:end);
TagawaDeltaH_int_interp = TagawaDeltaH_int_interp(order);

% simple linear interp from nearest neighbors
SFD_ext = lininterp1(StartM_interp,TagawaDeltaH_ext_interp,0);
SFD_int = lininterp1(StartM_interp,TagawaDeltaH_int_interp,0);

% Parabola interp (for noisy data, but shouldn't hurt clean data)
%polyfitmask = StartM_interp < .25 & StartM_interp > -0.25;
%SFD_ext_poly = polyfit(StartM_interp(polyfitmask),TagawaDeltaH_ext_interp(polyfitmask),2);
%SFD_ext = polyval(SFD_ext_poly,0);
%SFD_int_poly = polyfit(StartM_interp(polyfitmask),TagawaDeltaH_int_interp(polyfitmask),2);
%SFD_int = polyval(SFD_int_poly,0);

% calculate N and cluster if Ms and thickness are given

if isfield(DataIn,'Ms') && isfield(DataIn,'Thickness')
    Ms = DataIn.Ms;
    Thickness = DataIn.Thickness;
    Nd = SFD_ext/(4*pi*Ms);
    Cluster = Thickness * sqrt(1-Nd^2)/Nd;
   
    for i=1:nloops
        iField{i} = Field{i} - Nd*4*pi*Ms*Kerr{i};
    end
end

%% Assign output values to DataOut 

OutputVars = {
              'Kerr'
              'Hn'
              'Hc'
              'Hc1'
              'Hc2'
              'Hs'
              'SFD'
              'RemnantM'
              'StartingM'
              'SFD_ext'
              'SFD_int'
              'Nd'
              'Cluster'
              'H_Interp' 
              'M_Interp'
              'HcrMinorLoopH'
              'DeltaH_int'
              'DeltaH_ext'
              'TagawaDeltaH_ext'
              'TagawaDeltaH_int'
              'iField'
              };

for i=1:length(OutputVars)
    if exist(OutputVars{i},'var')
        eval(['DataOut.' OutputVars{i} '=' OutputVars{i} ';']);
    end
end

if makeplot
    % preliminary plotting routine
    figure;
    splot(DataOut,'Field','Kerr',jet)
end

end
