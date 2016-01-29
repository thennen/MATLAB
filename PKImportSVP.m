function varname = PKImportSVP(filepath)
% This function will import data from a PK SVP file, assuming that 'filepath' is the
% path to a data file of the correct format.
%
% T Hennen 2013

[pathstr, name, ext] = fileparts(filepath);
    
fid = fopen(filepath);

%% Read the header
for i = 1:13
    Header{i} = fgetl(fid);
end

% Create "Info" struct to hold values entered in PK program
Info.MeasurementDate = Header{3}(19:end);
Info.SetupFile = Header{6};
Info.DepToolID = Header{7};
Info.LotID = Header{11};
Info.DiskID = Header{12};
Info.OperationID = Header{13};

% Format for variable name
varname = [Info.DiskID '_PKParameters'];
% Change - to _
varname = strrep(varname,'-','_');
% Change . to ''
varname = strrep(varname,'.','');
% Ensure variable name is valid
varname = genvarname(varname);


% Get column names
tscan = textscan(fid,'%s',2);
colnames{1} = [tscan{1}{1} tscan{1}{2}];
while 1
    tscan = textscan(fid,'%s',1);
    if ~strcmp(tscan{1},';')
        colnames{end+1} = tscan{1}{1};
    else
        break
    end
end
fclose(fid);

% Import data table
DELIMITER = '\t';
Imported = importdata(filepath,DELIMITER,14);
NumDataCols = size(Imported.data,2);

if NumDataCols ~= length(colnames)
    error('catastrophic failure')
end

for i = 1:length(colnames)
    DataStruct.([genvarname(colnames{i}) '_blue']) = Imported.data(1:end/2,i);
    DataStruct.([genvarname(colnames{i}) '_red']) = Imported.data(end/2 + 1:end,i);
end

DataStruct.Info = Info;

% Add the data structure to the workspace
assignin('base',varname,DataStruct);
