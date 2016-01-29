function structspace(structname)
% Enter debug mode inside of a structure, so that you can modify its
% contents without typing the structure name a million times

dbstop at 19;
structure = evalin('base',structname);
fields = fieldnames(structure);

for i=1:length(fields)
    eval([fields{i} '= structure.' fields{i} ';']);
end

clear structure;
clear fields;
clear i;

% function is useless if this breakpoint doesn't work

cfields = who;

for i = 1:length(cfields)
    eval(['cstructure.' cfields{i} '=' cfields{i}]);
end

assignin('base',structname,rmfield(cstructure,'structname'));

end