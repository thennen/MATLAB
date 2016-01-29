function setfields(SearchString, varargin)
% Sets field values of multiple structures
% e.g. setfields('PKDynamic','fitrange',[15000, 20000])
% will create a field in all structures "matching" the word PKDynamic
% fitrange = [15000, 20000]
%
% Can be called with a cell array 'Cell', in which case all fields matching
% Cell{i,1} will be assigned a field named Cell{1,j} with value Cell{i,j}
%
% function is somewhat intelligent about matching variables to the search
% string.  'var_1' will match only 'var_1' if 'var_10' is present, 'var'
% will match both

switch(class(SearchString))
    case 'char'
        MatchingVars = findbestmatch(findbasevars(SearchString,'struct'),SearchString);
        if mod(nargin-1,2) ~= 0
            error('Invalid number of arguments');
        end
        FieldNames = varargin(mod([1:nargin-1],2) == 1);
        FieldVals = varargin(mod([1:nargin-1],2) == 0);
        for i = 1:length(MatchingVars)
            dummy = evalin('base',MatchingVars{i});
            for j = 1:(nargin-1)/2
                dummy.(FieldNames{j}) = FieldVals{j};
            end
            assignin('base',MatchingVars{i},dummy);
        end
        
    case 'cell'
        dimensions = size(SearchString);
        if nargin == 1
            % User didn't try to name the parameters, so assume that the
            % names are in the first row
            ParamNames = SearchString(1,2:end);     
        elseif dimensions(2)  == nargin
            % Can also name the parameters by varargin
            ParamNames = varargin;
        else
            error('Wrong number of arguments')
        end

        for i = 1:dimensions(1)
            MatchingVars = findbestmatch(findbasevars(SearchString{i,1},'struct'),SearchString{i,1});
            
            if ~isempty(MatchingVars)
                for l = 1:length(MatchingVars)
                    Matchl = evalin('base',MatchingVars{l});
                    for k = 1:length(ParamNames)
                        if ~isempty(ParamNames{k}) && ~isempty(SearchString{i,k+1})
                            % Add the new fields
                            Matchl.(genvarname(ParamNames{k})) = SearchString{i,k+1};
                        elseif ~isempty(ParamNames{k}) && isfield(Matchl,genvarname(ParamNames{k}))
                            % If the field exists and its value in SearchString is empty, remove that field
                            Matchl = rmfield(Matchl,genvarname(ParamNames{k}));
                        end
                    end
                    % Assign the variable back into the base workspace
                    assignin('base',MatchingVars{l},Matchl);
                end
            end
        end
    otherwise
        error('invalid SearchString variable class')
end
