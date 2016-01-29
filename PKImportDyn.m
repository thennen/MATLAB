function varname = PKImportDyn(filepath)
%PKIMPORTDYN Import data from a PK dynamic file
% Assuming that 'filepath' is the path to a data file of the correct format.
%
% T Hennen 2013


[a, filename, extension] = fileparts(filepath);
    
DELIMITER = '\t';
% Should be larger than total lines in file
HEADERLINES = 10000;
    
% Import entire file, line by line, into cell array "Imported"
Imported = importdata(filepath,DELIMITER,HEADERLINES);

numlines = length(Imported);

% Create "Info" struct to hold values entered in PK program
Info.MeasurementDate = Imported{3}(19:end);
Info.MacroFile = Imported{4}(18:end);
Info.DepToolID = Imported{8};
Info.LotID = Imported{10};
Info.DiskID = Imported{11};
Info.OperationID = Imported{12};

% Format for variable name
varname = [Info.DiskID '_PKDynamic'];
% Change - to _
varname = strrep(varname,'-','_');
% Ensure variable name is valid
varname = genvarname(varname);
j=1;
% If varname is already a base variable, append _# to varname
while evalin('base',['exist(''' varname ''',''var'')'])
    j = j+1;
    if j == 2
        varname = [varname '_' ndig(j,2)];
    else
        varname = [varname(1:end-2) ndig(j,2)];
    end
end

% These column numbers in the data file will be saved in the output struct
savecols = [1 3 4 7];
% with these names, respectively
colnames = {'Field','Kerr_x','Kerr_y','Time'};

numsections = 0;
% Loop through all lines to identify section starts
 for i = 1:numlines
     % Dyn Hc file happens to have a 'wait' command before each new loop
    if ~isempty(strfind(Imported{i},'Wait'));
    numsections = numsections + 1;
    sectionstart(numsections) = i+1;
    end
 end


for n = 1:numsections    
    % Import everything after section n start (this function will stop import at the next text line)    
    Imported = importdata(filepath,DELIMITER,sectionstart(n));
    for i = 1:length(savecols)
        % Save the specified column numbers of the data file to correspondingly named arrays in OutputStruct
        DataStruct.(colnames{i})(:,n) = Imported.data(1:end-1,savecols(i));
    end
end

DataStruct.Info = Info;
DataStruct.Field = -DataStruct.Field;
% PK recipe runs in strange direction
DataStruct.Kerr_x = -DataStruct.Kerr_x;
% Multiply by -1 to change sweep direction
DataStruct.Kerr_y = -DataStruct.Kerr_y;

% Run analysis routine
% DataStruct = PKDyn(DataStruct);

% Add the data structure to the workspace
assignin('base',varname,DataStruct);
