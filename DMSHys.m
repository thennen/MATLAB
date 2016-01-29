function DMSHys(Input,varargin)
% DMSHYS  Analyze DMS hysteresis data
% Subtracts background and makes various corrections to hysteresis data, and calculates relevant parameters
% 
% Analysis options can be set by passing arguments.  Arguments are parsed in the style of THargparse
% e.g. DMSHys('Sample1_PerpHys', 'drift', 'fitrange', [10000, 15000])
%
% options may also be set on a per-sample basis by defining fields corresponding to option names
% e.g. if you have a variable Sample1_PerpHys
% Sample1_PerpHys.drift = true
% Sample1_PerpHys.fitrange = [12000, 14000]
%

%% Set default analysis parameters

% Use only the branches of the hysteresis loop coming from saturation
options.fitbranch = true;
% Minor loop, no amplitude normalization, aligns saturation with previous
options.minor = false;
% Subtract background
options.subtract = true;
% Use the glass data to determine background
options.useglass = true;
% Use the same glass as the previously analyzed data (slope subtraction may differ)
options.sameglass = false;
% Use the exact same background as previously analyzed data
options.samebg = false;
% Sample has SUL
options.sul = false;
% Enter visualsubtraction routine for custom subtraction
options.visualsubtraction = false;
% Try to correct drift
options.drift = false;
% Plot corrected loop at the end
options.plot = false;
% Field range to use for fits
options.fitrange = []; % Default determined by maxfield later
% Field range to use for Hk calculation
options.hkfitrange = [3000, 10000];

%% Check for parameters in arguments
[options, paramopts] = THargparse(varargin, options);

%% Create cell array of all the variable names to be processed
switch class(Input)
    case 'char'
        % Do a search of the workspace for variables containing the input 
        % character string
        DMSvars = findbasevars(Input, 'struct');
        if isempty(DMSvars)
            error('No DMS variables match your search string')
        end
    case 'cell'
        % Cell array should be filled with sample names
        numvalidvars = 0;
        for i = 1:size(Input,1)
            VarMatches = findbasevars(Input{i}, 'struct');
            if ~isempty(VarMatches)
                numvalidvars = numvalidvars+1;
                % Only first match is used
                DMSvars{numvalidvars} = VarMatches{1}; 
            end
        end
        % if cell array has more than one column, assume that the other 
        % columns are parameters to be added to sample variables, and that the 
        % parameter names are contained in the first row
        if size(Input,2) > 1
            setfields(Input);
        end
    otherwise
        error('Invalid input variable type')
end

%% Loop through the array of determined variables for hysteresis routine
for z=1:length(DMSvars)
    DataIn = evalin('base',DMSvars{z});
    DataOut = DataIn;
    MaxField = max(DataIn.Field);
    Fangle = round(DataIn.Angle(1));
    % Just for convenience
    Field = DataIn.Field;
    X = DataIn.X;
    Y = DataIn.Y;
    
    % Take analysis parameters from key fields in structure
    fnames = fieldnames(options);
    for i=1:length(fnames)
        % If there is a field with an option as a name AND that option was not
        % set by input parameters
        if isfield(DataIn, fnames(i)) && ~any(strcmpi(fnames(i), paramopts))
            options.(fnames{i}) = DataIn.(fnames{i});
        end
    end
    
    % Determine fit range
    % If fitrange was not determined by arguments, or given as a field, set default
    if isempty(options.fitrange)
        options.fitrange = [0.5*MaxField MaxField];
    else
        % TODO: Check if fitrange is valid
    end

    %% subtract linear drift (or proportional to temperature)
    if options.drift
        Xdrift = X(end) - X(1);

        % linear drift
        %range = [0:1:length(X)-1].';
        %X = X - Xdrift*range/length(X);
        
        % temperature drift
        Tdrift = DataIn.Temp(end) - DataIn.Temp(1);
        X = X - Xdrift/Tdrift * (DataIn.Temp - DataIn.Temp(1));
    end

    if ~options.samebg
            % Determine which data to use in the tail fits, and do a line fit for each tail

            BranchMask = logical(ones(length(Field),1));
            BranchMask(floor(length(Field)/2):end) = false;
            if Field(1) < 0
                % In case loop starts from negative fields
                BranchMask = ~BranchMask;
            end

            if options.minor
                FitMask1 = Field > options.fitrange(1) & Field < options.fitrange(2);
                FitMask2 = FitMask1;
            elseif options.fitbranch
                % Fit only the hysteresis branches coming from saturation
                FitMask1 = Field > options.fitrange(1) & Field < options.fitrange(2) & BranchMask;
                FitMask2 = Field < -options.fitrange(1) & Field > -options.fitrange(2) & ~BranchMask;
            else
                % Fit everything in the specified field ranges
                FitMask1 = Field > options.fitrange(1) & Field < options.fitrange(2);
                FitMask2 = Field < -options.fitrange(1) & Field > -options.fitrange(2);
        end

        FitLine1 = polyfit(Field(FitMask1),X(FitMask1),1);
        FitLine2 = polyfit(Field(FitMask2),X(FitMask2),1);
        BGSlope = (FitLine1(1) + FitLine2(1))/2;
        BGSlopeDiff = (FitLine1(1) - FitLine2(1))/2;
        BGIntercept = (FitLine1(2) + FitLine2(2))/2;                            

        %% Look for glass data, calculate background

        % Regular expression to search for glass hysteresis loop
        % in the base workspace
        string = 'who (''-regexp'' ,''\w*(G|g)(lass|LASS)\w*(H|h)(YS|ys)\w*'')';
        GlassFile = evalin('base', string);

        FoundGlass = 0;
        if ~isempty(GlassFile) && options.useglass && options.subtract
            % Search for glass variable with the same field angle as DataIn
            for i = 1:length(GlassFile)
                if evalin('base',['round(' GlassFile{i} '.Angle(1))']) == Fangle
                    GlassFile = GlassFile{i};
                    FoundGlass = 1;
                    break;
                end    
            end
            if FoundGlass == 1
                % Copy glass data to variable
                GlassHys = evalin('base', GlassFile);

                if isfield(GlassHys,'PFitX') && isfield(GlassHys,'muX')
                    PFitX = GlassHys.PFitX;
                    muX = GlassHys.muY;

                    PFitY = GlassHys.PFitY;
                    muY = GlassHys.muY;
                else
                    % 13th order polynomial seems to be good fit without going overboard
                    [PFitX, ~, muX] = polyfit(GlassHys.Field,GlassHys.X,13);
                    %13th order not needed below image effect field, and is detrimental when glasshys datapoints are few
                    %[PFitX, ~, muX] = polyfit(GlassHys.Field,GlassHys.X,1);
                    GlassHys.PFitX = PFitX;
                    GlassHys.muX = muX;

                    % Do also for Y bg

                    % using line for now
                    [PFitY, ~, muY] = polyfit(GlassHys.Field,GlassHys.Y,1);
                    GlassHys.PFitY = PFitY;
                    GlassHys.muY = muY;

                    % assign to base workspace so fit need not be repeated
                    assignin('base', GlassFile, GlassHys);
                end

                %if Fangle == 90
                % poly fit subtraction
                GlassX = polyval(PFitX,Field,[],muX);
                GlassY = polyval(PFitY,Field,[],muY);
                %else
                %GlassX = GlassHys.X;  % Direct subtraction
                %GlassY = GlassHys.Y;  % Must use same recipe
                %end

                 if length(GlassHys.Field) == length(Field) && all((GlassHys.Field - Field) < 5)
                     disp(['Using direct subtraction of glass background for ' DMSvars{z}])
                     GlassX = GlassHys.X;  
                     GlassY = GlassHys.Y;  
                 end

                YBackground = GlassY;
                % Y bg is not observed to change with glass thickness variation
                RawCorrY = Y - GlassY;

                Background = GlassX;
                RawCorrX = X - Background;

                % Flag used to indicate that the same glass was used as the previous DMSvar, so 
                % same bg slope (2) should be subtracted
                if(options.sameglass)
                    Background = Background + BGSlope2*Field;

                    % center independently 
                    %FitLine1(2) = polyfit(Field(FitMask1),RawCorrX(FitMask1)-BGSlope2*Field(FitMask1),0);                      
                    %FitLine2(2) = polyfit(Field(FitMask2),RawCorrX(FitMask2)-BGSlope2*Field(FitMask2),0);
                    %BGIntercept2 = (FitLine1(2) + FitLine2(2))/2;
                    %Background = Background + BGIntercept2;
                else
                    % small straight line correction to account for glass thickness variation
                    FitLine1 = polyfit(Field(FitMask1),RawCorrX(FitMask1),1);
                    FitLine2 = polyfit(Field(FitMask2),RawCorrX(FitMask2),1);
                    BGSlope2 = (FitLine1(1) + FitLine2(1))/2;
                    %temporarily forcing BGSlope2 to constant value
                    %BGSlope2 = -0.0023292;
                    BGSlopeDiff2 = (FitLine1(1) - FitLine2(1))/2;
                    BGIntercept2 = (FitLine1(2) + FitLine2(2))/2; 
                    Background = Background + BGSlope2*Field + BGIntercept2 + BGSlopeDiff2*mean(options.fitrange);
                end
            else
                disp('Glass hysteresis data was found, but at a different field angle. Using straight line correction')
                Background = BGSlope*Field + BGIntercept + BGSlopeDiff*mean(options.fitrange);
                RawCorrY = Y;
            end
        elseif options.subtract
            disp('Glass hysteresis data was not found, or was specified to be ignored. Using straight line correction')
            % Call the BG a line
            Background = BGSlope*Field + BGIntercept + BGSlopeDiff*mean(options.fitrange);
            DataOut.BGSlope = BGSlope;
            RawCorrY = Y;
        else
            Background = 0; 
        end
    end
    
    %% Optional visual subtraction of line from X
    if options.visualsubtraction
        % nominalline = Field.*20./20000; % 10 uemu / 20,000 Oe
        nominalline = Field.*100./20000;
        if options.sul
            factor = 5*0.005/0.08;
            nominalline = Field.*0.08;
            [multiplier, inputs] = visualsubtraction(Field, X - Background, Field, nominalline, [1-factor:factor/20:1+factor]);
        else
            [multiplier, inputs] = visualsubtraction(Field, X - Background, Field, nominalline, [-1:0.05:1]);
        end
        
        % if the following errors, then there is not exactly one 'g' as an
        % input.
        if length(multiplier) == 1
            visualline = nominalline .* multiplier;
            DataOut.Analysis.VisualMultiplier = multiplier;
            Background = Background + visualline;
        else
            warning('There was not exactly 1 occurance of ''g'' in user inputs.  No visual correction applied.')
        end
    end

    %% Find Hs, Hn, and Hc or Hk by the extrapolation method
        RawCorrX = X - Background;
    % Need to refit after all corrections ...
     FitLine1 = polyfit(Field(FitMask1),RawCorrX(FitMask1),0);
     FitLine2 = polyfit(Field(FitMask2),RawCorrX(FitMask2),0);
     RawMs = (FitLine1(1) - FitLine2(1))/2;
    
    if abs(Fangle) < 30
        % Do fit for Hk
        % TODO: detect SUL
        if options.sul
            % Fit positive fields,    (For SUL)
            % negative fields
            FitMask1 = Field > options.hkfitrange(1) & Field < options.hkfitrange(2);
            FitMask2 = Field < -options.hkfitrange(1) & Field > -options.hkfitrange(2);
        else
            % Fit positive fields,    (For no SUL)
            % negative fields
            FitMask1 = RawCorrX > RawMs*.3 & RawCorrX < RawMs*.7;
            FitMask2 = RawCorrX < -RawMs*.3 & RawCorrX > -RawMs*.7;
        end
        
        % Fit lines for use in Hk calc.
        FitLine1 = polyfit(Field(FitMask1),RawCorrX(FitMask1),1);
        FitLine2 = polyfit(Field(FitMask2),RawCorrX(FitMask2),1);
        
        Hk1 = (RawMs - FitLine1(2)) / FitLine1(1);
        Hk2 = (-RawMs - FitLine2(2))/ FitLine2(1);
        Hk = (Hk1 - Hk2)/2;
        
        if options.sul
            SULMs1 = FitLine1(2);
            SULMs2 = FitLine2(2);
            SULMs = (SULMs1 - SULMs2)/2;
            SULSignal = SULMs*((Field > 0) - (Field < 0));
            
            DataOut.SULMs = SULMs;
            DataOut.SULSignal = SULSignal;
            DataOut.RawCorrXMinusSUL = RawCorrX - SULSignal;
        end
        
        DataOut.Analysis.HkLine1 = FitLine1(1)*Field + FitLine1(2);
        DataOut.Analysis.HkLine2 = FitLine2(1)*Field + FitLine2(2);
        
        DataOut.Hk = Hk;
        
        % This is here for a very stupid reason
        RawMsFull = RawMs;
    elseif options.minor
        % Shift loops up
        if exist('RawMsFull','var')
            RawCorrX = RawCorrX + RawMsFull;
        end
        
        FitMask1 = Field < 0 & Field > -5000;
        RawStartingM = polyfit(Field(FitMask1),RawCorrX(FitMask1),0);
        DataOut.RawStartingM = RawStartingM;
    else
        % Perp loop
        % Fit positive fields,  -RawMs*.7 < M < RawMs*.7
        % negative fields 
        FitMask1 = abs(Field) < options.fitrange(2) & RawCorrX < RawMs*.75 & RawCorrX > -RawMs*.75 & ~BranchMask;
        FitMask2 = abs(Field) < options.fitrange(2) & RawCorrX < RawMs*.75 & RawCorrX > -RawMs*.75 & BranchMask;
        
        %FitMask1 = RawCorrX < RawMs*.75 & RawCorrX > -RawMs*.75 & ~BranchMask;
        %FitMask2 = RawCorrX < RawMs*.75 & RawCorrX > -RawMs*.75 & BranchMask;
        
        % Don't bother if there aren't enough points in fitmask to do
        % reasonable fit
        if sum(FitMask1) > 1 && sum(FitMask2) > 1
            % Fit lines for use in Hs, Hn, Hc calc.
            FitLine1 = polyfit(Field(FitMask1),RawCorrX(FitMask1),1);
            FitLine2 = polyfit(Field(FitMask2),RawCorrX(FitMask2),1);
        else
            disp(['Not enough points found for Hc line fit in ' DMSvars{z}])
            % Turns everything to NaN
            FitLine1 = [NaN, NaN];
            FitLine2 = [NaN, NaN];
        end

        Hc1 = -FitLine1(2)/FitLine1(1);
        Hc2 = -FitLine2(2)/FitLine2(1);
        Hc = (Hc1 - Hc2)/2;
        
        Hs1 = (RawMs - FitLine1(2)) / FitLine1(1);
        Hs2 = (-RawMs - FitLine2(2))/ FitLine2(1);
        Hs = (Hs1 - Hs2)/2;
        
        Hn1 = (-RawMs - FitLine1(2)) / FitLine1(1);
        Hn2 = (RawMs - FitLine2(2))/ FitLine2(1);
        Hn = (Hn1 - Hn2)/2;
        
        DataOut.Hc = Hc;
        DataOut.Hc1 = Hc1;
        DataOut.Hc2 = Hc2;
        DataOut.Hn = Hn;
        DataOut.Hs = Hs;
        DataOut.Analysis.RawHcLine1 = FitLine1(1)*Field + FitLine1(2);
        DataOut.Analysis.RawHcLine2 = FitLine2(1)*Field + FitLine2(2);
        
        % Used if there are subsequent minor loops
        RawMsFull = RawMs;
    end

    if FoundGlass == 1
        DataOut.Analysis.GlassData = GlassFile;
    else
        DataOut.Analysis.GlassData = 'None';
    end
   
    if isfield(DataIn,'Area') && isfield(DataIn,'Thickness')
        Ms = RawMs*1000/DataIn.Area/DataIn.Thickness;
        CorrX = RawCorrX*1000/DataIn.Area/DataIn.Thickness;
        %DataOut.CorrX = CorrX;
        %DataOut.Ms = Ms;
        
        DataOut.Analysis.MsLine1 = ones(length(Field),1);
        DataOut.Analysis.MsLine1(1:end) = Ms;
        DataOut.Analysis.MsLine2 = -DataOut.Analysis.MsLine1;
        
        if Fangle ~= 0
            DataOut.Analysis.HcLine1 = DataOut.Analysis.RawHcLine1*1000/DataIn.Area/DataIn.Thickness;
            DataOut.Analysis.HcLine2 = DataOut.Analysis.RawHcLine2*1000/DataIn.Area/DataIn.Thickness;
        end
    end
    
%% Assign output values to DataOut if they exist

    OutputVars = {
                'RawCorrX'
                'RawCorrY'
                'YBackground'
                'RawMs'
                'Background'
                'Ms'
                'CorrX'
                };

    for i=1:length(OutputVars)
        if exist(OutputVars{i},'var')
            eval(['DataOut.' OutputVars{i} '=' OutputVars{i} ';']);
        end
    end

    %HalfField = interp1(RawCorrX(Field > 0 & Field < options.fitrange(1)),Field(Field > 0 & Field < options.fitrange(1)),RawMsFull/2);
    %DataOut.RawCorrX = RawCorrX;
    %DataOut.RawCorrY = RawCorrY;
    %DataOut.YBackground = YBackground;  %error when ybackground not defined due to straight line correction
    %DataOut.HalfField = HalfField;
    %DataOut.RawMs = RawMs;
    DataOut.Analysis.Field = Field;
    DataOut.Analysis.RawMsLine1 = ones(length(Field),1);
    DataOut.Analysis.RawMsLine1(1:end) = RawMs;
    DataOut.Analysis.RawMsLine2 = -DataOut.Analysis.RawMsLine1;
    DataOut.Analysis.FitRange = options.fitrange;
    DataOut.Analysis.HkFitRange = options.hkfitrange;
    DataOut.Analysis.options = options;
    %DataOut.Analysis.BGSlope2 = BGSlope2;  % error when no glass file
        %DataOut.Background = Background;

    assignin('base',DMSvars{z},DataOut);
    if options.plot
        % Plot loop separately
        figure
        plot(DataOut.Field, DataOut.RawCorrX)
        xlabel('Field (Oe)')
        ylabel('Raw Corr X (uemu)')
        title(strrep(DMSvars{z}, '_', ' '))
    end
end

end
