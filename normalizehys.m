function [NormHys, Amplitude, Offset, Sat1, Sat2, Slope] = normalizehys(Field,HysVals,FitRange,varargin)
% Centers and normalizes hysteresis loop
%
% Options for Slope correction, "kink" correction (pk SUL signal), branch masking
%
% Allows symmetric fit ranges
% FitRange = [10000, 15000])
% and Asymmetric fit ranges
% FitRange = [-15000, -10000, 5000, 10000])
%
% Sweep needs to be maxfield-minfield-maxfield or reverse
% (Cannot start from 0 field)
%

NormHys = [];
Amplitude = [];
Offset = [];
Sat1 = [];
Sat2 = [];
Slope = [];

%% Parse input options
P = inputParser;
P.addOptional('slope',false); % true, false, 'left', 'right'
P.addOptional('kink',false); % true, false, [2000, 5000] (linear region after kink)
P.addOptional('branch',false); % true, false
P.parse(varargin{:});
params = P.Results;

%% Construct fitting masks
% Mask according to fit ranges
if length(FitRange) == 2
    FitRange = [-FitRange(2) -FitRange(1) FitRange(1) FitRange(2)];
end
FitMask1 = Field >= FitRange(3) & Field <= FitRange(4);
FitMask2 = Field >= FitRange(1) & Field <= FitRange(2);
% Mask branches
BranchMask = true(size(Field));
BranchMask(floor(length(Field)/2):end) = false;
if Field(1) < 0
    % In case loop starts from negative fields
    BranchMask = ~BranchMask;
end
if params.branch
    FitMask1 = FitMask1 & BranchMask;
    FitMask2 = FitMask2 & ~BranchMask;
end


%% Try to fit lines to the saturation regions
fit1 = sum(FitMask1) > 1;
fit2 = sum(FitMask2) > 1;
if fit1
    SatLine1 = polyfit(Field(FitMask1),HysVals(FitMask1),1);
    Sat1 = polyfit(Field(FitMask1),HysVals(FitMask1),0);
end
if fit2
    SatLine2 = polyfit(Field(FitMask2),HysVals(FitMask2),1);
    Sat2 = polyfit(Field(FitMask2),HysVals(FitMask2),0);
end
if ~fit1 || ~fit2
    % normalizehys can be used to get saturation levels even if the whole
    % loop is not there.
    return;
end

%% Subtract Slope
if params.slope
    if strcmpi(params.slope,'left')
        Slope = SatLine2(1);
        SlopeBG = Slope.*Field;
    elseif strcmpi(params.slope,'right')
        Slope = SatLine1(1);
        SlopeBG = Slope.*Field;
    elseif strcmpi(params.slope,'asym')
        Slope = [SatLine1(1) SatLine2(1)];
        SlopeBG = (Field >= 0).*Field.*SatLine1(1) + (Field < 0).*Field.*SatLine2(1);
    else
        Slope = (SatLine1(1) + SatLine2(1)) / 2;
        SlopeBG = Slope.*Field;
    end
else
    Slope = 0;
    SlopeBG = Slope.*Field;
end

%% Normalize
NormHys = HysVals - SlopeBG;
FitLine1 = polyfit(Field(FitMask1),NormHys(FitMask1),0);
FitLine2 = polyfit(Field(FitMask2),NormHys(FitMask2),0);
Amplitude = abs((FitLine1(1) - FitLine2(1))/2);
Offset = (FitLine1(1) + FitLine2(1))/2;
NormHys = (NormHys - Offset) ./ Amplitude;


%% Remove kink
% Assumes full hysteresis is given
if params.kink
    if size(params.kink,1) == 1 && size(params.kink,2) == 2 && isfloat(params.kink)
        KinkRange = params.kink;
    else
        % Defines linear region after kink
        %KinkRange = [4000, 7000];
        KinkRange = [3000 4500];
    end
    
% "dual bg subtraction"
KinkMask1 = Field >= KinkRange(1) & Field <= KinkRange(2) & BranchMask;
KinkMask2 = Field >= -KinkRange(2) & Field <= -KinkRange(1) & ~BranchMask;
KinkLine1 = polyfit(Field(KinkMask1),NormHys(KinkMask1),1);
KinkLine2 = polyfit(Field(KinkMask2),NormHys(KinkMask2),1);
KinkSlope = (KinkLine1(1) + KinkLine2(1)) / 2;
KinkIntercept = (KinkLine1(2) - KinkLine2(2))/2;
KinkAmp = 1 - KinkIntercept;
KinkField = KinkAmp/KinkSlope;
KinkBG = ((Field >= KinkField) .* KinkAmp) + ((Field <= -KinkField) .* -KinkAmp) + (abs(Field) < KinkField).*KinkSlope.*Field;
NormHys = (NormHys - KinkBG) * 1/(1-KinkAmp);

% Tyler's kink BG subtraction


% renormalize?


end


end
