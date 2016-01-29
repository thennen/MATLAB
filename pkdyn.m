function DataOut = pkdyn(DataIn,varargin)
%PKDYN analyze PK Dynamic data.
% DataIn should be a struct with matrix fields 'Field', 'Time','Kerr_x', and
% 'Kerr_y', with each column of the matrices corresponding to a field sweep
% measurement done at different rates.
%
% Analysis options can be set by command line arguments, or by defining a field
% inside the data structure.  e.g.
% >>sample1_PKDyn = pkdyn(sample1_PKDyn, 'red', 'f0', 10^10, 'n', 1.5)
% is equivalent to 
% >>sample1_PKDyn.red = 1
% >>sample1_PKDyn.f0 = 10^10
% >>sample1_PKDyn.n = 1.5
% >>sample1_PKDyn = pkdyn(sample1_PKdyn)
%
% One of the Sharrock parameters f0 and n can take multiple values, and the 
% Sharrock fit will be done using all of the values.  e.g.
% >>sample1_PKDyn = pkdyn(sample1_PKDyn, 'n', [1:.25:2])
% >>sample1_PKDyn.KVkT_Hn
%   68.4088   77.9067   86.7319   94.9345  102.5615
% 
% 
% DataOut is a struct which contains all the fields in DataIn, with
% additional fields resulting from the analysis routine.
% 
%
% Field Name        Description   *columns represent varying Sharrock parameters
% ----------        -----------
% Kerr              Centered and normalized red or blue kerr loops in each column
% Hc                1d array of Hc values for each kerr loop
% SweepRate         1d array of Sweep rates for each loop
% KVkT              KVkT value extracted for each SharrockM level, *
% KVkT_Hc           KVkT interpolated at Hc, *
% KVkt_Hn           KVkT interpolated at Hn, *
% H0                Sharrock short time extrapolated H vs M, *
% H0_Hn             Sharrock short time extrapolated H at Hn
% H0_Hc             Sharrock short time extrapolated H at Hc
% SFD_0             Sharrock short time SFD, *
% Analysis.FitRange Field range used for fitting the saturation region for loop normalization
% Analysis.n        n value(s) used for Sharrock fits
% Analysis.f0       f0 value(s) used for Sharrock fits
% Analysis.varargin Arguments passed to analysis function
%
% If DataIn contains field 'SFD_ext' (automatically added by pk function if recoil data 
% is present) the following fields will also be created:
%
% iField            2d array, columns contain demag corrected field for each loop
% iH0               Demag corrected sharrock short time extrapolated H vs M, *
% iSFD_0            Demag corrected Sharrock short time SFD, *           
% iH0_Hn            Demag corrected Sharrock short time extrapolated H at Hn
% iH0_Hc            Demag corrected Sharrock short time extrapolated H at Hc
% iKVkT_Hn          Demag corrected KVkT(Hn)
% iKVkT_Hc          Demag corrected KVkT(Hc)
% iKVkT             Demag corrected KVkT for all SharrockM
% Eb_Hn             Eb(Hn, Hd = 4piN(Hndef))
% ...

% Get the name of the input variable
varname = inputname(1);

% Error if input is not a struct. 
if ~isstruct(DataIn)                                    
    error(['Input ' varname ' is not a struct.'])
end

%% Default analysis parameters
options.red = false;      % True for red, false for blue
options.kink = false;     % Remove "kink" from sul signal
options.slope = false;    % Not currently used
options.fast = false;     % Speed up analysis 
options.plot = false;     % Make a plot at the end of analysis
options.fitrange = FindFitRange(DataIn);
options.f0 = 10^9;
options.n = 2;

% Normalized magnetization level used to define the nucleation point
HnDef = 0.9;
% Definition for Sharrock fit function at bottom

% Parse input arguments to options
[options, paramopts] = THargparse(varargin, options);

% Set analysis parameters from key fields in structure
% only if they weren't set by input argument
fnames = fieldnames(options);
for i=1:length(fnames)
    % if there is a field with an option as a name and that option was not
    % set by input parameters
    if isfield(DataIn, fnames(i)) && ~any(strcmpi(fnames(i), paramopts))
        options.(fnames{i}) = torquein.(fnames{i});
    end
end

% Save some typing later...
f0 = options.f0;
n = options.n;
FitRange = options.fitrange;
Field = DataIn.Field;
Time = DataIn.Time;                                   

% Determine if Sharrock extrap. is to be repeated for a series of f0 or n
% values, as determined by the length of the f0 or n variables
f0length = length(f0);
nlength = length(n);
if f0length > 1 && nlength > 1
    error('Cannot vary two Sharrock parameters at a time')
elseif f0length > 1
    NParameters = f0length;
    SVaryingParameter = 'f0';
elseif nlength > 1
    NParameters = nlength;
    SVaryingParameter = 'n';
else
    NParameters = 1;
    SVaryingParameter = 'none';
end

% Values of M to do Sharrock extrapolation
% Should contain +-0.5 and HnDef
if options.fast
    SharrockM = [-HnDef -0.5 0 0.5 HnDef].';
else
    SharrockM = [-0.98:0.02:0.98].';     
end
MNum = length(SharrockM);
                                                           
% Analyze blue or red data
if options.red == 0
    Kerr = DataIn.Kerr_x;                             
elseif options.red == 1
    Kerr = DataIn.Kerr_y;                             
end

%% Calculate sweep rates using timestamp
% Use sweep rate calculated as average point-by-point rate, skipping the last 
% 20 points, where the field may stop changing
dH = Field(2:end,:) - Field(1:end-1,:);               
dt = Time(2:end,:) - Time(1:end-1,:);                 
RateMatrix = dH ./ dt;                                
SweepRate = abs(mean(RateMatrix(1:end-20,:)));

%% Loop through kerr loops, contained in columns of 'Kerr', normalize and find Hc, Hn
for j = 1:size(Field,2)
    Fieldj = Field(:,j);
    Loopj = Kerr(:,j);
    
    % Center and normalize loop j
    [Loopj, ampj, offsj] = normalizehys(Fieldj,Loopj,FitRange);
    Kerr(:,j) = Loopj;  % Add to output matrix
    amp(:,j) = ampj;    % amplitude of each loop if curious
    offs(:,j) = offsj;  % offs of each loop
    
    if options.kink && isfield(DataIn,'Analysis')
    % Look for SUL signal, which may have been copied to the Analysis field by the pk 
    % function if it found a matching recoil variable
        if isfield(DataIn.Analysis,'SULsignal')
            % Subtract SUL signal
            Loopj = Loopj - interp1(DataIn.Analysis.SULfield,DataIn.Analysis.SULsignal,Fieldj,'linear','extrap');
            Kerr(:,j) = normalizehys(Fieldj,Loopj,FitRange);
        end
    end
    
    % Interpolate loops for M values for use in sharrock extrap
    [UniqueLoopj,m1,n1] = unique(Loopj);
    UniqueFieldj = Fieldj(m1);
    SharrockH(:,j) = interp1(UniqueLoopj,UniqueFieldj,SharrockM);
                                                                              
    % Find Hc by line fit interpolation
    if length(FitRange) == 2
        % Mask loops outside of  -.25 < M < .25
        HcFitMask = abs(Fieldj) < FitRange(2) & Loopj < .25 & Loopj > -.25;
    elseif length(FitRange) == 4
        HcFitMask = abs(Fieldj) < FitRange(4) & Loopj < .25 & Loopj > -.25;
    end
    % Fit lines for use in Hc calculation
    HcFitLine = polyfit(Fieldj(HcFitMask),Loopj(HcFitMask),1);
    Hc(j) = -HcFitLine(2)/HcFitLine(1);
    
    % Find Hn (fit a parabola to Hn region
    HnFitMask = abs(-Loopj - HnDef) < 0.05;
    HnFitPoly = polyfit(Fieldj(HnFitMask),Loopj(HnFitMask),2);
    Hn(j) = (-HnFitPoly(2)+sqrt(HnFitPoly(2)^2 - 4*HnFitPoly(1)*(HnFitPoly(3)+HnDef)))/(2*HnFitPoly(1));
    
end

%% Sharrock extrapolation
% Initialize progress bar
progressbari(0, MNum*NParameters);

a0 = [50, 8000];   % Define start values for fit: a(1) = KV/kT,  a(2) = H0
ub = [500, 30000]; % Upper bound
lb = [0, 0];       % Lower bound
ia0 = [50, 8000];  % Start values for demag corrected extrapolation

% Suppresses the display of certain fit notifications
fitoptions = optimset('Display', 'off');
%fitoptions = optimset(); % uncomment if you want to display them

if isfield(DataIn, 'SFD_ext')
    % Check to see if SFD_ext from corresponding recoil data was added
    % (usually by the pk function).  If so, demag correction is done.
    SFD_ext = DataIn.SFD_ext;
    SharrockHd = SFD_ext*SharrockM;
    iField = Field - SFD_ext*Kerr;
    DemagGiven = true;
else
     DemagGiven = false;
end

% Some variables for plotting the sharrock fit
SharrockExtrap.Fit_Rate = logspace(1,11).';
SharrockExtrap.Fit_InverseRate = 1./SharrockExtrap.Fit_Rate;
SharrockExtrap.Hc = Hc;
SharrockExtrap.Hn = Hn;
SharrockExtrap.Rate = SweepRate;
SharrockExtrap.InverseRate = 1./SweepRate;

% Preallocate arrays for speed
KVkT = zeros(MNum,NParameters);
H0 = zeros(MNum,NParameters);
KVkT_Hn = zeros(1,NParameters);
KVkT_Hc = zeros(1,NParameters);
SFD_0 = zeros(1,NParameters);
H0_Hn = zeros(1,NParameters);
H0_Hc = zeros(1,NParameters);
if DemagGiven
    iKVkT = zeros(MNum,NParameters);
    iH0 = zeros(MNum,NParameters);
    iSFD_0 = zeros(1,NParameters);
    iKVkT_Hn = zeros(1,NParameters);
    iKVkT_Hc = zeros(1,NParameters);
    iH0_Hn = zeros(1,NParameters);
    iH0_Hc = zeros(1,NParameters);
    Eb_Hn = zeros(1,NParameters);
end

for p = 1:NParameters
    % Sharrock function for lsqcurvefit needs fit parameters in an array (a)
    switch SVaryingParameter
        case 'f0'
            FitSharrock = @(a, R) Sharrock(a(1), a(2), R, f0(p), n);
        case 'n'
            FitSharrock = @(a, R) Sharrock(a(1), a(2), R, f0, n(p));
        case 'none'
            FitSharrock = @(a, R) Sharrock(a(1), a(2), R, f0, n);
    end
    
    % For every M level specified to do the Sharrock extrapolation
    for m = 1:MNum
        % Sharrock fit at SharrockM(m)
        fitparams = lsqcurvefit(FitSharrock,a0,abs(SweepRate),SharrockH(m,:),lb,ub,fitoptions);
        KVkT(m,p) = fitparams(1);
        H0(m,p) = fitparams(2);
        % Start next fit at the parameters found in the previous fit (might save time)
        a0 = fitparams;
        
        if SharrockM(m) == 0
            % Save the fit functions as determined at Hc, for plotting
            SharrockExtrap.Fit_Hc(:,p) = FitSharrock(a0,SharrockExtrap.Fit_Rate);
        end
        
        if SharrockM(m) == -HnDef
            % Save the fit functions as determined at Hn, for plotting
            SharrockExtrap.Fit_Hn(:,p) = FitSharrock(a0,SharrockExtrap.Fit_Rate);
        end
        
        % If Hd is given, also do Sharrock extrapolation with demag correction
        if DemagGiven
            % Sharrock fit at SharrockM(m)
            ifitparams = lsqcurvefit(FitSharrock, ia0, abs(SweepRate), (SharrockH(m, :) - SharrockHd(m)), lb, ub, fitoptions);
            iKVkT(m, p) = ifitparams(1);
            iH0(m, p) = ifitparams(2);
            ia0 = ifitparams;
            
            if SharrockM(m) == 0
                % Save the fit functions as determined at Hc, for plotting
                SharrockExtrap.Fit_iHc(:, p) = FitSharrock(ia0, SharrockExtrap.Fit_Rate);
            end
            
            if SharrockM(m) == -HnDef
                % Save the fit functions as determined at Hn, for plotting
                SharrockExtrap.Fit_iHn(:, p) = FitSharrock(ia0, SharrockExtrap.Fit_Rate);
            end
        end
        progressbari(m+(p-1)*MNum, MNum*NParameters);
    end
    
    KVkT_Hn(p) = KVkT(SharrockM == -HnDef, p);
    KVkT_Hc(p) = KVkT(SharrockM == 0, p);
    SFD_0(p) = H0(SharrockM == 0.5, p) - H0(SharrockM == -0.5, p);
    H0_Hn(p) = H0(SharrockM == -HnDef, p);
    H0_Hc(p) = H0(SharrockM == 0, p);
    
    if DemagGiven
        iSFD_0(p) = iH0(SharrockM == 0.5, p) - iH0(SharrockM == -0.5, p);
        iKVkT_Hn(p) = iKVkT(SharrockM == -HnDef, p);
        iKVkT_Hc(p) = iKVkT(SharrockM == 0, p);
        iH0_Hn(p) = iH0(SharrockM == -HnDef, p);
        iH0_Hc(p) = iH0(SharrockM == 0, p);
        SharrockExtrap.iHn = Hn + HnDef*SFD_ext;
        
        % Calculate energy barrier in field
        if strcmp(SVaryingParameter,'n')
            Eb_Hn(p) = iKVkT_Hn(p) * (1 - HnDef*SFD_ext/iH0_Hn(p))^n(p);
        else
            Eb_Hn(p) = iKVkT_Hn(p) * (1 - HnDef*SFD_ext/iH0_Hn(p))^n;
        end
    end
end
     

%% Assign output values to DataOut

DataOut = DataIn;

OutputVars = {
              'iField'
              'Hc'
              'Hn'
              'Kerr'
              'SweepRate'
              'SharrockM'
              'KVkT'
              'KVkT_Hc'
              'KVkT_Hn'
              'H0'
              'H0_Hn'
              'H0_Hc'
              'SFD_0'
              'iH0'
              'iH0_Hn'
              'iH0_Hc'
              'iSFD_0'
              'iKVkT_Hn'
              'iKVkT_Hc'
              'iKVkT'
              'Eb_Hn'
              'SharrockExtrap'
              'amp'
              'offs'
              };

% Add OutputVars to structure if they exist
for i=1:length(OutputVars)
    if exist(OutputVars{i},'var')
        eval(['DataOut.' OutputVars{i} '=' OutputVars{i} ';']);
    end
end

DataOut.Analysis.FitRange = FitRange;
DataOut.Analysis.n = n;
DataOut.Analysis.f0 = f0;
DataOut.Analysis.varargin = varargin;
DataOut.Analysis.red = options.red;
DataOut.Analysis.options = options;

if options.plot
    % preliminary plotting routine
    plotvarname = strrep(varname, '_', ' ');
    figure;
    subplot(121);
    splot(DataOut,'Field','Kerr');
    title(plotvarname);
    subplot(122);
    sharrockplot(DataOut);
    title(plotvarname);
end

end

function SharrockVal = Sharrock(KVkT,H0,R,f0,n)
    SharrockVal = H0 * (1 - ( (1/KVkT)*log(f0*H0/2 / KVkT ./ R)).^(1/n));  
end

function progressbari(CurrentVal,MaxVal)
if CurrentVal == 0;
    fprintf(1,'[          ]');
else
    progress = floor(CurrentVal*10 / MaxVal); % 0-10
    
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b'); % Erase bar
    for j = 1:progress
        fprintf(1,'-'); % Print 1-10 '-'
    end
    for k = 1:(10-progress)
        fprintf(1,' ');
    end
    fprintf(1,']');
    
    if CurrentVal == MaxVal
        fprintf('\n')
    end
end

end


