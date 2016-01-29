function renamemamr(searchstring)
%puts power and freq in mamr structure's name

mamrvars = findbasevars(searchstring,'struct');

for i = 1:length(mamrvars)
    
    if isbasefield(mamrvars{i},'Power') && isbasefield(mamrvars{i},'Freq')
        
        
        % decimals will be rounded to nearest whole number
        % bad things will happen if multiple mamr structs have power and freq
        % that round to the same whole number.
        
        
        freqi = ndig(evalin('base',[mamrvars{i} '.Freq']),2);
        poweri = evalin('base',[mamrvars{i} '.Power']);
        
        if poweri < 0
            poweri = ['n' ndig(poweri,2)];
        else
            poweri = ndig(poweri,2);
        end
            
        
        
        newvarname = genvarname([searchstring '_' poweri 'dbm_' freqi 'ghz']);
        
        
        if ~strcmp(newvarname,mamrvars{i})
            %rename var
            evalin('base',[newvarname '=' mamrvars{i}]);

            %delete old var
            evalin('base',['clear ' mamrvars{i}]);
        end
    end
    
    
end


end

