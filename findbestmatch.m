function BestMatch = findbestmatch(ArrayOfStrings,StringToMatch)
%BESTMATCH return the best matching string within a cell array of strings
%ArrayOfStrings by comparing them to substrings separated in StringToMatch 
%by underscores.  The string(s) in the ArrayOfStrings which contain the
%most substrings will be returned as a cell array of strings BestMatch.
%This function is case insensitive.
%If there is ambiguity about the best match, all match candidates will be
%returned, sorted from shortest to longest string
%If no substrings of StringToMatch are found in any string in
%ArrayOfStrings, then an empty cell is returned

if isempty(ArrayOfStrings) || isempty(StringToMatch)
    BestMatch = {};
    return
end

Substrings = regexp(StringToMatch, '_', 'split');

% For each string in ArrayOfStrings, count number of matches
matches = zeros(length(ArrayOfStrings),1);
for i = 1:length(ArrayOfStrings)
    Substringi = regexp(ArrayOfStrings{i}, '_', 'split');
    matches(i) = length(intersect(Substrings, Substringi));
    
    %for j = 1:length(Substrings)
    %    matches(i) = matches(i) + ~isempty(regexpi(ArrayOfStrings{i},Substrings{j}));
    %end
end

%Return string with highest number of matches
MaxMatches = max(matches);
if MaxMatches > 0
    BestMatch = ArrayOfStrings(matches == max(matches));
else
    % Return empty cell:
    % BestMatch = {};
    % OR return all in ArrayOfStrings which contain StringToMatch:
    % for some reason the following is the first thing I thought of to
    % accomplish this ..
    BestMatch = ArrayOfStrings(~cellfun(@isempty, regexpi(ArrayOfStrings,['\w*' StringToMatch '\w*'])));
end

%If ArrayOfStrings contains StringToMatch, return StringToMatch only

if any(strcmp(StringToMatch,ArrayOfStrings))
    BestMatch = {StringToMatch};
end

%Sort from shortest to longest strings

matchlengths = cellfun(@length, BestMatch);
[a,indices] = sort(matchlengths);
BestMatch = BestMatch(indices);


end