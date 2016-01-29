function bool = isbasefield(basevarname,field)
% Checks if the base workspace contains a struct with name basevarname
% and field with name field

if evalin('base',['exist(''' basevarname ''',''var'')'])
    bool = evalin('base',['isfield(' basevarname ',''' field ''')']);
else
    bool = false;
end

end

