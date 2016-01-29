function [MatrixOut] = extractexcel(SearchString,varargin)
%EXTRACTEXCEL create cell array from data structures to paste into excel
% Extract any number of vectors or 1d cell arrays inside any number of
% structures and put them all in a cell array which can be pasted into excel

MatchingVars = findbasevars(SearchString);
nfields = nargin - 1;
numvars = length(MatchingVars);

MatrixOut = {};

for i = numvars:-1:1
    MatrixOut{1,(i-1)*nfields+1} = MatchingVars{i};
    for j = nfields:-1:1
        ColumnNumber = (i-1)*nfields+j;
        MatrixOut{2,ColumnNumber} = varargin{j};
        if isbasefield(MatchingVars{i},varargin{j})
            Columni = evalin('base',[MatchingVars{i} '.' varargin{j}]);
            if isnumeric(Columni)
                MatrixOut(3:length(Columni)+2,ColumnNumber) = num2cell(Columni);
            elseif iscell(Columni)
                error('you have not written a script to handle cell arrays')
            end
        else
            %MatrixOut(2:end,ColumnNumber) = [];
        end
    end
    
end

end

