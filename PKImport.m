function varname = PKImport(filepath)
%PKIMPORT General PK data file import
% 
% Saves all PK commands to cell array 'AllPKCommands'
% Saves PK commands that correspond to the imported data to cell array 'PKDataCommands'

[pathstr, name, ext] = fileparts(filepath);

fid = fopen(filepath);

%% Read the header
for i = 1:14
    Header{i} = fgetl(fid);
end

% Create "Info" struct to hold values entered in PK program
Info.MeasurementDate = Header{3}(19:end);
Info.MacroFile = Header{4}(18:end);
Info.DepToolID = Header{8};
Info.LotID = Header{10};
Info.DiskID = Header{11};
Info.OperationID = Header{12};


%% Generate a unique variable name
varname = [Info.DiskID '_PKData'];
varname = strrep(varname,'-','_');
varname = genvarname(varname);
j=1;
% If varname is already a base variable, append _## to varname
while evalin('base',['exist(''' varname ''',''var'')'])
    j = j+1;
    if j == 2
        varname = [varname '_' ndig(j,2)];
    else
        varname = [varname(1:end-2) ndig(j,2)];
    end
end

% Command number
cn = 0;

%% Read initial commands until data begins
while ~feof(fid);
% Note beginning of line position
bol = ftell(fid);
% Try to read 7 numbers
dscan = textscan(fid,'%f',7);
    % if 7 numbers are read
    if length(dscan{1}) == 7
        % Move back to bol, to prepare for data import
        fseek(fid,bol,'bof');
        break;
    else
        cn = cn + 1;
        % Read the rest of the line (should be the pk command)
        AllPKCommands{cn} = fgetl(fid);
    end
end

% These column numbers in the data file will be saved in the output struct
savecols = [1 3 4 7];                                        
colnames = {'Field','Kerr_x','Kerr_y','Time'};

%% Read datablocks
% datablock number
j = 0;
while ~feof(fid)
    % Data point number
    n = 1;
    % Flag indicating whether this loop iteration found data
    founddata = 0;
    % Try to read 7 numbers
    dscan = textscan(fid,'%f',7);
    if length(dscan{1}) == 7
        j = j + 1;
        founddata = 1;
    end
    % If there's a data line, save it to DataStruct
    while length(dscan{1}) == 7
        for i = 1:length(savecols)                        
            DataStruct.(colnames{i}){j}(n) = dscan{1}(savecols(i));
        end
        n = n+1;
        dscan = textscan(fid,'%f',7);
    end
    
    cn = cn + 1;
    % Read the rest of the text line (should be PK command)
    AllPKCommands{cn} = fgetl(fid);
    if founddata
        % If the command was read immediately after a data block, save it
        PKDataCommands{j} = AllPKCommands{cn};
    end
end

fclose(fid);

% Find measurement duration
Duration = DataStruct.Time{end}(end) - DataStruct.Time{1}(1);
Hours = floor(Duration/3600);
Minutes = floor((Duration - Hours*3600)/60);
Seconds = Duration - Minutes * 60 - Hours * 3600;
DurationStr = [ndig(Hours,2) ':' ndig(Minutes,2) ':' ndig(Seconds,2)];
Info.MeasDuration = DurationStr;

%% Assign in imported data with generated varname
DataStruct.Info = Info;
DataStruct.AllPKCommands = AllPKCommands';
DataStruct.PKDataCommands = PKDataCommands;

assignin('base',varname,DataStruct);
