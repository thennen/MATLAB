function varname = PKImportSVD(filepath)
% This function will import data from a PK SVD file, assuming that 'filepath' is the
% path to a data file of the correct format.
%
% T Hennen 2013

% These column numbers in the data file will be saved in the output struct
savecols = [1 3 4];
% with these names, respectively
colnames = {'Field','Kerr_x','Kerr_y'};

[~, filename, extension] = fileparts(filepath);
    
DELIMITER = '\t';
% Should be larger than total lines in file
HEADERLINES = 10000;
    
% Import entire file into cell array "Imported"
Imported = importdata(filepath,DELIMITER,HEADERLINES);
numlines = length(Imported);

% Create "Info" struct to hold values entered in PK program
Info.MeasurementDate = Imported{3}(19:end);
Info.SetupFile = Imported{6};
Info.DepToolID = Imported{7};
Info.LotID = Imported{11};
Info.DiskID = Imported{12};
Info.OperationID = Imported{13};

Commas = strfind(Imported{14},',');
R = Imported{14}(1:Commas(1)-1);
Theta = Imported{14}(Commas(1)+1:Commas(2)-1);
%R_Theta = [str2double(R) str2double(Theta)].';

% Format for variable name
varname = [Info.DiskID '_PKLoop'];
% Change - to _
varname = strrep(varname,'a-','_');
% Change . to ''
varname = strrep(varname,'.','');
% Ensure variable name is valid
varname = genvarname(varname);

varnameexists = evalin('base',['exist(''' varname ''')']);
if varnameexists
    DataStruct = evalin('base',varname);
    % Find the number of loops already imported into the struct
    sizeofexistingstruct = size(DataStruct.Field,2);
    
    
    % Checking for duplicate points seemed like a good idea at the time...
    %for j=1:sizeofexistingstruct
        %if DataStruct.R_Theta(:,j) == R_Theta;
        %    error('Already imported this loop location')
        %end
    %end
else
    sizeofexistingstruct = 0;
end



numsections = 0;
% Loop through all lines to identify section starts
for i = 1:numlines
    % SVD file happens to have 'N/A' before each new section
    if ~isempty(strfind(Imported{i},'N/A'));
    numsections = numsections + 1;
    sectionstart(numsections) = i;
    end
end

% Need to count the number of data points already concatenated
Totaldatapoints = 0;
for n = 1:numsections    
    % Import everything after section n start (this function will stop import at the next text line)    
    Imported = importdata(filepath,DELIMITER,sectionstart(n));
    ImportedLength = size(Imported.data,1);
    for i = 1:length(savecols)
        % Save the specified column numbers of the data file to
        % correspondingly named arrays in OutputStruct, concatenating sections
        DataStruct.(colnames{i})(Totaldatapoints+1:Totaldatapoints+ImportedLength,sizeofexistingstruct+1) = Imported.data(1:end,savecols(i));
    end
    Totaldatapoints = Totaldatapoints + ImportedLength;
end


DataStruct.Info = Info;
DataStruct.R(sizeofexistingstruct+1) = R;
DataStruct.Theta(sizeofexistingstruct+1) = Theta;

% Add the data structure to the workspace
assignin('base',varname,DataStruct);
disp(['Added ' Info.DiskID ' loop location R,Theta = ' R ',' Theta])
