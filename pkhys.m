function DataOut = pkhys(DataIn , varargin)
%PKHYS Analyze pk hysteresis (svd) data
% Not actually finished...

if ~isstruct(DataIn)
    error('Input is not a struct.')
end

DataOut = DataIn;

options.red = 1; % Choose red or blue laser data for analysis. 0 for red, 1 for blue
options.kink = 0;
options.slope = 0;
options.plot = 0;

options = THargparse(varargin, options);

% Determine fit range for fitting saturation levels
FitRange = FindFitRange(DataIn);

Field = DataIn.Field;

if options.red == 1
    Kerr = DataIn.Kerr_y;
elseif options.red == 0
    Kerr = DataIn.Kerr_x;
end

numloops = size(Field,2);
%% Center and normalize loops (independently)
for i=1:numloops
    [Kerr(:,i) , Ampl(i), Offs(i), ~, ~] = normalizehys(Field(:,i),Kerr(:,i),FitRange,'kink',options.kink,'slope',options.slope,'branch',1);
end


%% Find loop parameters (Hn, Hc, Hs) 
HnDef = 0.9;
HsDef = 0.9;

BranchMask = true(length(Field),1);
BranchMask(floor(length(Field)/2):end) = false;
if Field(1) < 0
    % In case loop starts from negative fields
    BranchMask = ~BranchMask;
end

% preallocate

Hn1 = zeros(1,numloops);
Hn2 = zeros(1,numloops);
Hn = zeros(1,numloops);
Hc1 = zeros(1,numloops);
Hc2 = zeros(1,numloops);
Hc = zeros(1,numloops);
Hs1 = zeros(1,numloops);
Hs2 = zeros(1,numloops);
Hs = zeros(1,numloops);
SFD1 = zeros(1,numloops);
SFD2 = zeros(1,numloops);
SFD = zeros(1,numloops);


for i=1:numloops
    % flipud because lininterp1 assumes monotonic increase
    Field1 = flipud(Field(BranchMask,i));
    Field2 = Field(~BranchMask,i);
    Loop1 = flipud(Kerr(BranchMask,i));
    Loop2 = Kerr(~BranchMask,i);
    
    Hn1(i) = lininterp1(Loop1,Field1,HnDef);
    Hn2(i) = lininterp1(Loop2,Field2,-HnDef);
    Hn(i) = (-Hn1(i) + Hn2(i)) /2;
    Hc1(i) = lininterp1(Loop1,Field1,0);
    Hc2(i) = lininterp1(Loop2,Field2,0);
    Hc(i) = (-Hc1(i) + Hc2(i))/2;
    Hs1(i) = lininterp1(Loop1,Field1,-HsDef);
    Hs2(i) = lininterp1(Loop2,Field2,HsDef);
    Hs(i) = (-Hs1(i) + Hs2(i))/2;
    SFD1(i) = lininterp1(Loop1,Field1,0.5) - lininterp1(Loop1,Field1,-0.5);
    SFD2(i) = lininterp1(Loop2,Field2,0.5) - lininterp1(Loop2,Field2,-0.5);
    SFD(i) = (SFD1(i) + SFD2(i))/2;
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
              'Ampl'
              'Offs'
              };
for i=1:length(OutputVars)
    if exist(OutputVars{i},'var')
        eval(['DataOut.' OutputVars{i} '=' OutputVars{i} ';']);
    end
end

DataOut.Analysis.FitRange = FitRange;

%% Plot
if options.plot
    figure;
    plot(DataOut.Field,DataOut.Kerr);
    xlabel('Field')
    ylabel('Norm Kerr')
end

end
