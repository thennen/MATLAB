function Report = pk(varargin)
%PK Automatically find and analyze PK data in the base workspace.
% Based on an input search string or cell array of
% strings, the function will find variables which contain the string,
% determine if they are PK variables and if so, call the correct PK
% analysis function. The correct analysis function should be specified by
% keywords in the variable names such as "PKRecoil" or "PKDynamic"
%
% Variable name          Analysis Function
% ----------------------------------------
% SomeSample_PKLoop      pkhys
% SomeSample_PKRecoil    pkrecoil
% SomeSample_PKDynamic   pkdyn
%
% If input is a cell array, additional columns may be used to set analysis 
% parameters as fields in the data structure (see setfields function)
% e.g. 
% C = {'Name',     'red',     'FitRange';
%      'sample1',      1, [10000, 15000];
%      'sample2',      0, [11000, 16000]}
%
% An attempt is made to match PKDynamic data with PKRecoil data to pass
% mean field information for use for "intrinsic" calculation (e.g. iKVkT)
% in the Dynamic data.  Matching is done by matching keywords in the
% variable names separated by underscores.
%
% e.g.
% Recoil var             Recognized as corresponding to
% ---------------------------------------------------
% Sample1_PKRecoil       Sample1_PKDynamic
% Sample1_PKRecoil_2     Sample1_PKDynamic_2
%
% Any additional arguments will be passed through to the analysis functions
% e.g.
% pk('Sample1_PKRecoil', 'fitrange', [15000, 20000])
% 
% More examples:
% Analyze all recognized pk data
% pk
%
% Analyze only variables containing the string 'Sample1'
% pk('Sample1')
% 
% See Also:
% PKHYS, PKRECOIL, PKDYN
%
% Tyler Hennen 2014

%% Create cell array of all the variable names to be processed
if nargin == 0
    % Structs for PK analysis indicated by being a struct with 'PK' in the
    % name.
    PKVars = findbasevars('PK','struct');
    remvarargin = {};
elseif nargin >= 1
    switch class(varargin{1})
        % If pk is passed a string, do a search of the workspace for 
        % structs containing the string. Wildcard * allowed
        case 'char'
            PKVars = findbasevars(varargin{1},'struct');
            if isempty(PKVars)
                error('No variables match your search string')
            end
        % If pk is passed a cell array, the first column should be filled
        % with sample names, other columns can set analysis parameters by
        % naming the parameter in the first row (e.g. 'FitRange') and the
        % value in the corresponding sample row (e.g. [15000, 20000])
        case 'cell'
            PKVars = {};
            CellArray = varargin{1};
            % Loop through first column of cell array and search for
            % variables containing that string.  ALL matching variables 
            for i = 1:size(CellArray,1) 
                VarMatches = findbestmatch(findbasevars(CellArray{i},'struct'),CellArray{i});
                if ~isempty(VarMatches)
                    PKVars = [PKVars; VarMatches]; 
                end
            end
            % If cell array has more than one column, assume that the other
            % columns are parameters to be added to sample variables, and 
            % that the parameter names are contained in the first row
            if size(CellArray,2) > 1
                setfields(CellArray);
            end
        otherwise
            error('Invalid input variable type')
    end
    
    if nargin > 1
        remvarargin = varargin(2:end);
    else
        remvarargin = {};
    end
end

DynVars = findbestmatch(PKVars,'PKDynamic');
RecoilVars = findbestmatch(PKVars,'PKRecoil');
HysVars = findbestmatch(PKVars,'PKLoop');

totalnumvars = length(DynVars) + length(RecoilVars) + length(HysVars);

for h=1:length(HysVars)
    evalin('base',[HysVars{h} '= pkhys(' HysVars{h} vararg2strlist(remvarargin) ');']);
    disp([HysVars{h} ' done (' num2str(h) '/' num2str(totalnumvars) ')']);
end
for i=1:length(RecoilVars)
        evalin('base',[RecoilVars{i} '= pkrecoil(' RecoilVars{i} vararg2strlist(remvarargin) ');']);
        disp([RecoilVars{i} ' done (' num2str(length(HysVars)+i) '/' num2str(totalnumvars) ')']);
end
for j=1:length(DynVars)
    % Look for Nd, Ms from matching recoil variable, and if found, copy them to the PKDynamic
    % structure for use in demag corrected Sharrock extrapolation
    RecoilMatch = findbestmatch(RecoilVars,DynVars{j});
    if ~isempty(RecoilMatch)
        if isbasefield(RecoilMatch{1},'Nd')
            evalin('base',[DynVars{j} '.Nd =' RecoilMatch{1} '.Nd;']);
        end
        if isbasefield(RecoilMatch{1},'Ms')
            evalin('base',[DynVars{j} '.Ms =' RecoilMatch{1} '.Ms;']);
        end
        if isbasefield(RecoilMatch{1},'SFD_ext')
            evalin('base',[DynVars{j} '.SFD_ext =' RecoilMatch{1} '.SFD_ext;']);
        end
        if isbasefield(RecoilMatch{1},'Analysis')
            try
                evalin('base',[DynVars{j} '.Analysis.SULfield =' RecoilMatch{1} '.Analysis.SULfield;']);
                evalin('base',[DynVars{j} '.Analysis.SULsignal =' RecoilMatch{1} '.Analysis.SULsignal;']);
            end
        end
        % Save to the structure the variable name that was determined to be
        % the recoil data for the corresponding sample
        evalin('base',[DynVars{j} '.Analysis.RecoilMatch = '' ' RecoilMatch{1} ''';']);
    end
    if length(RecoilMatch) > 1
        warning(['Ambiguous PKRecoil matches for ' DynVars{j} '.  Used ' RecoilMatch{1}])
    end
                
                
    disp([DynVars{j} ' fitting ... (' num2str(length(HysVars) + length(RecoilVars) + j) '/' num2str(totalnumvars) ')']);
    
    %findbestmatch(RecoilVars,DynVars{i})
        evalin('base',[DynVars{j} '= pkdyn(' DynVars{j} vararg2strlist(remvarargin) ');']);
        
end

end

function stringout = vararg2strlist(cellin)
    stringout = '';
    if ~isempty(cellin)
        for i=1:length(cellin)
            switch class(cellin{i})
                case 'char'
                    stringout = [stringout ',''' cellin{i} ''''];
                case 'double'
                        stringout = [stringout ', [' num2str(cellin{i}) ']'];
            end
        end
    else
        stringout = '';
    end
end
