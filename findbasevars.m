function [varNames, firstVar] = findbasevars(searchString,varargin)
%FINDBASEVARS return a cell array containing the names of all variables in 
% the base workspace which contain it.  If an exact match to SearchString
% is found, it will be returned as the first cell in the array.

% The function also returns variable with the name that appears first in
% varNames

% varargin optionally specifies the class of variable to search for

if isempty(searchString)
    varNames = '';
    firstVar = '';
else
    varNames = evalin('base',['who(''*' searchString '*'')']);

    exactMatch = strcmp(varNames,searchString);

    if any(exactMatch)
        varNames = [{searchString}; varNames(~exactMatch)];
    end
    
    if ~isempty(varNames)
        firstVar = evalin('base',varNames{1});
    else
        firstVar = '';
    end
end

if(nargin == 2)
% Second argument should be a string naming the desired class of the
% returned variable names
    for i=length(varNames):-1:1
        if ~strcmp(evalin('base',['class(' varNames{i} ')']),varargin{1})
            varNames(i) = [];
        end
    end
elseif(nargin > 2)
    error('Too many arguments.')
end

end

