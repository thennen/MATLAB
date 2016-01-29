function pareval(searchString,command)
%PAREVAL "Parallel Evaluate"
%Evaluates 'command' in the base workspace for all structure names
%containing searchString.  The string '()' in command will be replaced by
%the matching variable names. (i) will be replaced by the index of the
%variable

BaseVars = findbasevars(searchString,'struct');

for i = 1:length(BaseVars)
    newcommand = strrep(command,'()',BaseVars{i});
    newcommand = strrep(newcommand,'(i)',num2str(i));
    try
        evalin('base',newcommand);
    catch
        disp(['Command failed for var ' BaseVars{i}])
    end
end

end