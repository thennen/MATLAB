function renamevars(searchstring, renamecommand)
% Supposed to rename variables that contain searchstring according to
% renamecommand, which can use parts of the original variable name like
% this renamecommand = '(4:10)_appendthis'    (end keyword is allowed)
% the word 'ii' will be replaced by the variables alphabetical index among
% all the matches
% Hennen 2014

switch(class(searchstring))
    case 'char'
        oldvarnames = findbasevars(searchstring);
    case 'cell'
        % can pass variable names directly
        % e.g. renamevars(who('C*'),'something_(1:end)')
        oldvarnames = searchstring;
end

for i=1:length(oldvarnames)
    newvarname{i} = strrep(renamecommand, 'ii', num2str(i));
    
    %find ( ), don't screw this up because there's no check
    leftbracket = strfind(renamecommand,'(');
    rightbracket = strfind(renamecommand,')');
    
    for j = 1:length(leftbracket) 
        slicestr = renamecommand(leftbracket(j):rightbracket(j));
        oldsubs = eval(['oldvarnames{' num2str(i) '}' slicestr]);
        newvarname{i} = strrep(renamecommand, slicestr, oldsubs);
    end
    
    newvarname{i} = genvarname(newvarname{i});
end
% check if any newvarnames are the same, if so, error
if length(unique(newvarname)) ~= length(newvarname)
    error('Some of the new var names generated are duplicates')
end
for i=1:length(oldvarnames)
    evalin('base',[newvarname{i} '=' oldvarnames{i} ';']);
    evalin('base',['clear ' oldvarnames{i} ';']);
end

end

