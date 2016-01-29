function [MatrixOut,Colnames] = collectvect(SearchString,Field)
%COLLECTVECT combine vectors contained in multiple structures into a matrix
% Useful if for some horrible reason you decide to paste into a spreadsheet

% You may have screwed this up.
MatchingVars = findbasevars(SearchString);

numvars = length(MatchingVars);
MatrixOut = zeros(1,numvars);

for i = numvars:-1:1
    Colnames{1,i} = MatchingVars{i};
    
    if evalin('base',['isfield(' MatchingVars{i} ',''' Field ''')'])
        Columni = evalin('base',[MatchingVars{i} '.' Field]);
        for j = 1:length(Columni)
            MatrixOut(j,i) = Columni(j);
        end
    else
        Colnames(i) = [];
        MatrixOut(:,i) = [];
    end
    
end

end

