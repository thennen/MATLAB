function dms(varargin)
% DMS  Automatically find and analyze VSM data from the base workspace.
%  dms('searchstring') Analyze data containing searchstring in the name
%  dms(CellArrayOfStrings) Analyze
%
% Based on an input search string or cell array of strings, the function will 
% find variables which contain the string, determine if they are DMS variables 
% and if so, call the corresponding DMS analysis function.  The analysis 
% function is determined by keywords in the variable names such as "PerpHys" or 
% "Torque"
%
% Variable name          Analysis Function
% ----------------------------------------
% SomeSample_PerpHys     DMSHys
% SomeSample_IPHys       DMSHys
% SomeSample_DCD         DMSDCD
% SomeSample_Torque      DMSTorque
%
% If input is a cell array, additional columns may be used to set analysis 
% parameters as fields in the data structure (see setfields function)
% e.g. 
% C = {'Name',   'Area', 'Thickness',     'FitRange';
%      'sample1',  50.1,           8, [10000, 15000];
%      'sample2',  48.9,           9, [11000, 16000]}
% dms(C) Will calculate Ms based on given area and thickness, and will subtract 
% the background using the given fit range
%
% An attempt is made to match PerpHys data with background measurements and
% subtract accordingly.  Torque data will also be matched with corresponding
% PerpHys data and background torque data.  Matching is done by matching keywords
% in the variable names separated by underscores.  Background data is denoted by
% containing the word "glass" in the variable name.
%
% e.g.
% Background data name   Recognized as background for
% ---------------------------------------------------
% Glass_PerpHys          sample1_PerpHys, Sample2_PerpHys
% Glass_PerpHys_2        sample1_PerpHys_2, 2_sample1_PerpHys
% Glass_IPHys            sample1_IPHys, SAMPLE5_IPHys
% Glass_Torque           sample1_Torque, sample2_Torque
% Glass_Torque_18koe     sample1_Torque_18koe
% 
% PerpHys data name      Recognized for pairing with torque vars
% --------------------------------------------------------------
% sample1_PerpHys        sample1_Torque, sample1_Torque_18  
%
% Any additional arguments will be passed through to the analysis functions
% e.g.
% dms('Sample1_PerpHys', 'fitrange', [15000, 20000])
% 
% More examples:
% Analyze all recognized dms data
% dms
%
% Analyze only variables containing the string 'Sample1'
% dms('Sample1')
%
% Analyze variables containing either PerpHys or Torque
% a = {'PerpHys';'Torque'}
% dms(a)
% 
% See Also:
% DMSHYS, DMSTORQUE, DMSDCD
%
% Tyler Hennen 2013

%% Create cell array of all the variable names to be processed
if nargin == 0
    DMSVars = evalin('base','who');
    for i=length(DMSVars):-1:1
        if(~isempty(regexpi(DMSVars{i},'Glass')))
            DMSVars(i) = [];
        elseif evalin('base',['~isstruct(' DMSVars{i} ');'])
            DMSVars(i) = [];
        end
    end
else
    switch class(varargin{1})
        case 'char'
            % Do a search of the workspace for variables containing the input character string
            DMSVars = findbasevars(varargin{1},'struct');
            if isempty(DMSVars)
                error('No variables match your search string')
            end
        case 'cell'
            % Cell array should be filled with sample names
            DMSVars = {};
            CellArray = varargin{1};     
            for i = 1:size(CellArray,1)
                % Loop through first column of cell array and search for variables containing that string
                VarMatches = findbestmatch(findbasevars(CellArray{i},'struct'),CellArray{i});
                if ~isempty(VarMatches)
                    DMSVars = [DMSVars; VarMatches];
                end
            end
            
            if size(CellArray,2) > 1
                % If cell array has more than one column, assume that the other columns are parameters
                % to be added to sample variables, and that the parameter names are contained in the first row
                setfields(CellArray);
            end
        otherwise
            error('Invalid input variable type')
    end
end

if nargin > 1
    remvarargin = varargin(2:end);
else
    remvarargin = {};
end

%% Find and categorize all the names of DMS variables in workspace
vars = evalin('base','who');

whostring = 'who (''-regexp'' ,''\w*(G|g)(lass|LASS)\w*'')'; 
GlassVars = evalin('base',whostring);
% Remove any variables containing the word glass
PerpHysVars = findbestmatch(DMSVars,'PerpHys');
IPHysVars = findbestmatch(DMSVars,'IPHys');
DCDVars = findbestmatch(DMSVars,'DCD');
TorqueVars = findbestmatch(DMSVars,'Torque');

%% Pass the DMS data to the corresponding analysis function
for i=1:length(PerpHysVars)
    if isempty(regexpi(PerpHysVars{i},'glass'))
        DMSHys(PerpHysVars{i}, varargin{2:end});
        disp([PerpHysVars{i} ' done.'])
    end
end

for i=1:length(IPHysVars)
    if isempty(regexpi(IPHysVars{i},'glass'))
        % Always fit all with IPHys, and pass remaining varargs to DMSHys
        DMSHys(IPHysVars{i}, 'fitbranch', false, varargin{2:end})
        disp([IPHysVars{i} ' done.'])
    end
end

for i=1:length(DCDVars)
    % TODO: Identify groups of DCD belonging to same sample
    % subtract same BG from all of these
    evalin('base',[DCDVars{i} '= DMSDCD(' DCDVars{i} vararg2strlist(remvarargin) ');']);
    disp([DCDVars{i} ' done.'])
end

for i=1:length(TorqueVars)
    PerpHys = {};
    if isempty(regexpi(TorqueVars{i},'Glass'))
        PerpHys = findbestmatch(PerpHysVars,TorqueVars{i});
        GlassTorque = findbestmatch(GlassVars,TorqueVars{i});
        
        if isempty(GlassTorque)
            warning(['Skipping analysis for ' TorqueVars{i} '. No Torque background data identified.']);
            continue;
        end
        
        if isempty(PerpHys)
            evalin('base',[TorqueVars{i} '=DMSTorque(' TorqueVars{i} ',' GlassTorque{1} vararg2strlist(remvarargin) ');']);  
            % Indicate no PerpHys data was used
            evalin('base',[TorqueVars{i} '.Analysis.PerpHys = ''Nominal'';']);
        else
            evalin('base',[TorqueVars{i} '=DMSTorque(' TorqueVars{i} ',' GlassTorque{1} ',' PerpHys{1} vararg2strlist(remvarargin) ');']);
            % Indicate which PerpHys data was used
            evalin('base',[TorqueVars{i} '.Analysis.PerpHys = ''' PerpHys{1} ''';']);
        end

        % Indicate which angular bg data was used
        evalin('base',[TorqueVars{i} '.Analysis.AngBG = ''' GlassTorque{1} ''';']);

        disp([TorqueVars{i} ' done.'])
    end
end

%generate a report
%if ~isempty(PerpHysVars)
%    assignin('base','Perp_Hys_Results',collectfields(PerpHysVars,'RawMs','Ms','Hc','Hn','Hs'));
%end
%if ~isempty(IPHysVars)
%    assignin('base','IP_Hys_Results',collectfields(IPHysVars,'RawMs','Ms','Hk'));
%end
%if ~isempty(TorqueVars)
%    assignin('base','Torq_Results',collectfields(TorqueVars,'N2amp','Ms','K1p','K1','Hk'));
%end

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
