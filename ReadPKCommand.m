function cmd=ReadPKCommand(commandline)

%{
This function reads a command line string in the form: string;number;number;string;string;number;number

Example: 'X Field Sweep and Measure MOKE Signal;6;15000;;True;0.015;10;'

It returns a cell array contanining each string and number, in the same order they
appear in the original command line.
%}


% Determine index (in string) 
idx = [1 (strfind(commandline,';')+1)];

% Allocate memory
cmd = cell(1,length(idx)-1);
cmd(:) = {'NaN'};

% Customized vector determining whether each element of the cell array will be string or number.
% *** MUST BE CHANGED IF COMMAND LINE FORMAT VARIES ***
isnum = [0 1 1 0 0 1 1];

% Saves numbers and strings into cell array "cmd"
for i = 1:(length(idx)-1)
    if (idx(i+1)-idx(i)-1)
        if isnum(i)
            cmd{i} = str2double(commandline(idx(i):(idx(i+1)-2)));
        else
            cmd{i} = commandline(idx(i):(idx(i+1)-2));
        end    
    end    
end

end