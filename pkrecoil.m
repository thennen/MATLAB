function DataOut = pkrecoil(DataIn , varargin)
% This function analyses PK recoil data.  DataIn should be a struct with cell
% array fields 'Field', 'Time','Kerr_x', and 'Kerr_y', with each cell of the
% cell arrays corresponding to a field sweep measurement starting from a
% different reverse field.  DataOut is a struct which contains all the fields
% in DataIn, with additional fields resulting from the analysis routine.
% 
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

% TODO: Recoil slope

% Get the name of the input variable
varname = inputname(1);

% Error if input is not a struct. 
if ~isstruct(DataIn)                                    
    error(['Input ' varname ' is not a struct.'])
end

% Field should be a cell array containing the field array for each sweep
Field = DataIn.Field;
nloops = length(Field);
MaxField = max(Field{1});

% Set default options
options.red = false;    % blue is default
options.kink = false;   % true, false, [a, b]
options.slope = false;  % true (avg), 2 (asym), 3 (left), 4 (right)
options.fast = true;    % true, false
options.plot = false;   % true, false
options.drift = true;   % true, false
% Fitrange priority:  Command line, Field, Base Var, default
options.fitrange = FindFitRange(DataIn);

% Parse input arguments TH style
options = THargparse(varargin, options);

FitRange = options.fitrange;

if ~isstruct(DataIn)
    error('Input is not a struct.')
end

if options.red == true
    Kerr = DataIn.Kerr_y;
else
    Kerr = DataIn.Kerr_x;
end

% convert to 4 value fitrange
if length(FitRange) == 2
    FitRange = [-FitRange(2) -FitRange(1) FitRange];
elseif length(FitRange) ~= 4
    error('Something wrong with FitRange')
end

% Range used to correct offset drift in minor loops
%MinorFitRange = FitRange;
MinorFitRange = [MaxField*0.85, MaxField];

%% Basic center and normalize loops (based on major loop)
% Fit each major hys branch, with possible asymmetric fit range
[~, Ampl1, Offs1, UOffs1] = normalizehys(Field{1},Kerr{1},-fliplr(FitRange));
[~, Ampl2, Offs2, UOffs2] = normalizehys(Field{2},Kerr{2},FitRange);
Ampl = (Ampl1 + Ampl2)/2;
Offs = (Offs1 + Offs2)/2;
UOffs = (UOffs1 + UOffs2)/2;
Ampldiscrep = abs(Ampl1-Ampl2)*100/Ampl;
%Offsdiscrep = abs(Offs1-Offs2)*100/Ampl;
if  Ampldiscrep > 1
    disp('Amplitude of the two major loop sweeps differ by %.2f%%',Ampldiscrep)
    disp('Using second sweep for loop normalization')
    Ampl = Ampl2;
    Offs = Offs2;
    UOffs = UOffs2;
    % Could also try to correct the first sweep, by scaling or by subtracting a
    % linear drift
end
for i = 1:nloops
    Kerr{i} = (Kerr{i} - Offs) ./ Ampl;
end
% Correct offset drift
[~, ~, ~, MUOffs] = normalizehys(Field{1}, Kerr{1},MinorFitRange);
for i = 3:nloops
    % correct for offset drift by matching avg M levels in MinorFitRange
    [~, ~, ~, UOffsi] = normalizehys(Field{i},Kerr{i},MinorFitRange);
    Kerr{i} = Kerr{i} - UOffsi + MUOffs;
end

%% Optional Slope Correction
if options.slope == 1;
    %Normal single slope
    %[~, Ampl1, Offs1, UOffs1, ~,Slope1] = normalizehys(Field{1},Kerr{1},[-20000,-16500,10000,20000],'slope',true);
    %[~, Ampl2, Offs2, UOffs2, ~,Slope2] = normalizehys(Field{2},Kerr{2},[-20000,-10000,16500,20000],'slope',true);
    [~, Ampl1, Offs1, UOffs1, ~,Slope1] = normalizehys(Field{1},Kerr{1}, FitRange,'slope',true);
    [~, Ampl2, Offs2, UOffs2, ~,Slope2] = normalizehys(Field{2},Kerr{2}, FitRange,'slope',true);
    Ampl = (Ampl1 + Ampl2)/2;
    Offs = (Offs1 + Offs2)/2;
    Slope = (Slope1 + Slope2)/2;
    UOffs = (UOffs1 + UOffs2)/2;
    Kerr{1} = (Kerr{1}-Slope*Field{1} - Offs) ./ Ampl; %
    Kerr{2} = (Kerr{2}-Slope*Field{2} - Offs) ./ Ampl; % do the same thing to both major loop branches
    for i = 3:nloops
       % Use major loop normalization
       Kerr{i} = (Kerr{i}-Slope*Field{i})./Ampl;
       % Align saturation regions independently
       [~, ~, ~,UOffsi, ~] = normalizehys(Field{i},Kerr{i},FitRange); 
       Kerr{i} = Kerr{i} - UOffsi + 1;
    end
elseif options.slope == 2
    % Asymmetric slope
    MajField = [Field{1} Field{2}];
    MajKerr = [Kerr{1} Kerr{2}];
    [~, Ampl, Offs, ~, ~, asymslope] = normalizehys(MajField,MajKerr,[6000,15000],'slope','asym','branch',true);
    for i = 1:nloops
        SlopeBG = (Field{i} >= 0).*Field{i}.*asymslope(1) + (Field{i} < 0).*Field{i}.*asymslope(2);
        Kerr{i} = ((Kerr{i} - SlopeBG) - Offs) ./ Ampl;
    end
elseif options.slope == 3
    % to do
elseif options.slope == 4
    % to do
end

%% Kink correction (SUL signal subtraction)
if options.kink;
    if length(options.kink) == 2
        SlopeRange = options.kink;
    else
        SlopeRange = [2000,5000];
    end

    BranchAvgField = (Field{2} - Field{1})/2;
    BranchAvgKerr = (Kerr{2}-Kerr{1})/2;
    FitMask = BranchAvgField >= -SlopeRange(2) & BranchAvgField <= -SlopeRange(1);
    FitSlope = polyfit(BranchAvgField(FitMask),BranchAvgKerr(FitMask),1);

    % extrapolate sul signal for field <  sloperange(1)
    %  extrap region
    region1 = BranchAvgField < 0 & BranchAvgField > -SlopeRange(1);
    % rest of neg region
    region2 = BranchAvgField <= -SlopeRange(1);

    SULsignal = [BranchAvgKerr(region2), polyval(FitSlope,BranchAvgField(region1))];
    
    %remove offset
    SULsignal = SULsignal - FitSlope(2);
    %mirror
    SULfield = [BranchAvgField(region1 | region2), -fliplr(BranchAvgField(region1 | region2))];
    SULsignal = [SULsignal, -fliplr(SULsignal)];
    
    Analysis.SULfield = SULfield;
    Analysis.SULsignal = SULsignal;
    
    %subtract SUL signal from all loops
    for i=1:nloops
        Kerr{i} = Kerr{i} - interp1(SULfield,SULsignal,Field{i},'linear','extrap');
    end
    
    % renormalize (there's probably a better way to do this)
    [~, Ampl1, Offs1, UOffs1] = normalizehys(Field{1},Kerr{1},-fliplr(FitRange));
    [~, Ampl2, Offs2, UOffs2] = normalizehys(Field{2},Kerr{2},FitRange);
    Ampl = (Ampl1 + Ampl2)/2;
    Offs = (Offs1 + Offs2)/2;
    UOffs = (UOffs1 + UOffs2)/2;
    %Ampldiscrep = abs(Ampl1-Ampl2)*100/Ampl;
    %Offsdiscrep = abs(Offs1-Offs2)*100/Ampl;
    
    for i=1:nloops
        Kerr{i} = (Kerr{i} - Offs) ./ Ampl;
    end
    
end

%% Subtract linear amplitude drift
if options.drift
    % Currently no flexibility for different order minor loops
    % Try to determine if last loop reached saturation
    % estimate Hs, check if reverse field is significantly above that
    FitMask1 = Field{2} >= FitRange(1) & Field{2} <= FitRange(2);
    DriftFit1 = polyfit(Field{2}(FitMask1),Kerr{2}(FitMask1),1);
    FitMask2 = abs(Kerr{1}) < 0.25;
    DriftFit2 = polyfit(Field{1}(FitMask2),Kerr{1}(FitMask2),1);
    estim_Hs = (DriftFit2(2) - DriftFit1(2)) / (DriftFit1(1) - DriftFit2(1));
    DriftTotal = mean(Kerr{end}(FitMask1)) - mean(Kerr{2}(FitMask1));
    
    if DataIn.ReverseFields(end) > estim_Hs - 1000
        disp('Last minor loop starts nowhere near saturation, no drift correction applied')
    else
        Driftinc = DriftTotal / (nloops-2);
        for i = 3:nloops
            Kerr{i} = (Kerr{i}-1) * (1 + (i-2)*Driftinc/2) + 1;
        end
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
HcDiff = abs(Hc1 + Hc2);
if HcDiff > 150
    disp([varname ' has large difference in left/right Hc: ' num2str(HcDiff) ' Oe'])
end

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
% create logical masks for fitting remnant curves near zero field
Fieldmask = cellfun(@(a)a<1000 & a>-1000,Field,'UniformOutput',0);
% logical array indicating whether or not each loop has enough field values
% near zero for a parabola fit
CrossesZero = cellfun(@(b) sum(b) > 2,Fieldmask);

for i=1:nloops
        StartingM(i) = Kerr{i}(1);
    if CrossesZero(i) == 1
        Fieldi = Field{i};
        Loopi = Kerr{i};
        Fieldmaski = Fieldmask{i};
        %fit parabola to fieldmask region
        ParaFit = polyfit(Fieldi(Fieldmaski),Loopi(Fieldmaski),2);
        RemnantM(i) = ParaFit(3);
    else
        %RemnantM(i) = NaN;
        RemnantM(i) = 0;
    end
end

%% Find H(M=1/2) for each minor loop  (Mr < Ms/2)

if options.fast
    %M levels to extract H(M) for each curve
    M_Interp = linspace(-.98,.98,99).';
else
    %M levels to extract H(M) for each curve
    M_Interp = linspace(-.99,.99,397).';
end

%M level for Tagawa analysis
MTagawa = 0.5;
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

H_HalfReturned = NaN(nloops,1);
H_MajorHalfReturned = NaN(nloops,1);
DeltaH_int = NaN(length(M_Interp),nloops);

for i = 3:nloops
    for k = 1:length(M_Interp)
        if ~options.fast && StartingM(i) < M_Interp(k)
            H_Interp(k,i) = lininterp1(Kerr{i},Field{i},M_Interp(k));
        elseif options.fast && abs(StartingM(i)) < 0.1 && StartingM(i) < M_Interp(k)
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
    % is this even a thing?
    DeltaH_ext(:,i) = H_Interp(:,i) - flipud(H_InterpRight);
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

% Create arrays for interpolation, sorted for interpolation function
[StartM_interp, order] = sort(StartingM(3:end));
TagawaDeltaH_ext_interp = TagawaDeltaH_ext(3:end);
TagawaDeltaH_ext_interp = TagawaDeltaH_ext_interp(order);
TagawaDeltaH_int_interp = TagawaDeltaH_int(3:end);
TagawaDeltaH_int_interp = TagawaDeltaH_int_interp(order);

% Simple linear interp from nearest neighbors
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
DataOut = DataIn;

Analysis.FitRange = FitRange;
Analysis.varargin = varargin;
Analysis.options = options;

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
              'Analysis'
              };

for i=1:length(OutputVars)
    if exist(OutputVars{i},'var')
        eval(['DataOut.' OutputVars{i} '=' OutputVars{i} ';']);
    end
end

if options.plot
    % preliminary plotting routine
    plotvarname = strrep(varname, '_', ' ');
    figure;
    splot(DataOut,'Field','Kerr',jet);
    title(plotvarname);
end

end
