function varname = PKImportRecoil(filepath)
%PKIMPORECOIL Import recoil data from PK data file
%
% T Hennen 2013

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
% Format for variable name
varname = [Info.DiskID '_PKRecoil'];
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


%% Read lines until data starts

while ~feof(fid);
    % note beginning of line position
    bol = ftell(fid);
    % try to read 7 numbers
    dscan = textscan(fid,'%f',7);
        if length(dscan{1}) == 7
            % Move back to bol
            fseek(fid,bol,'bof');
            break;
        else
            % Read the rest of the line (should be the pk command) into command
            command = fgetl(fid);
            % Can do something with the line here if desired
        end
end

%% Read the major loops (first two)

% These column numbers in the data file will be saved in the output struct
savecols = [1 3 4 7];
% with these names, respectively
colnames = {'Field','Kerr_x','Kerr_y','Time'};

%Read 7 columns at a time, stop at text
n=1;
dscan = textscan(fid,'%f',7);
% First half major loop
while length(dscan{1}) == 7
    % save specified columns
    for i = 1:length(savecols)
        DataStruct.(colnames{i}){1}(n) = dscan{1}(savecols(i)); 
    end
    n = n+1;
    dscan = textscan(fid,'%f',7);
end
%Read the rest of the text line to get to next data line
fgetl(fid);
n=1;
dscan = textscan(fid,'%f',7);
% Second half major loop
while length(dscan{1}) == 7
    for i = 1:length(savecols)                        
        DataStruct.(colnames{i}){2}(n) = dscan{1}(savecols(i));
    end
    n = n+1;
    dscan = textscan(fid,'%f',7);
end
%Read the rest of the text line
fgetl(fid);

%% Read the minor loops (data file format changes)

% First two loops are major loops, so no reverse field associated
ReverseFields(1:2) = NaN;

% Minor loop number
j = 3;
while ~feof(fid)
    
    % Read all commands between minor loops, in case there are more than
    % one for some reason
    while ~feof(fid);
        % note beginning of line position
        bol = ftell(fid);
        % try to read 7 numbers
        dscan = textscan(fid,'%f',7);
        % if succeeded
        if length(dscan{1}) == 7
            % Move back to bol, to prepare to import data
            fseek(fid,bol,'bof');
            break;
        else
            % Read the rest of the line (should be the pk command) into command
            command = fgetl(fid);
        end
    end
    % store the last read command to read the reverse field
    tline = command;
    
    % Data point number
    n = 1;
    % Read reverse field from text line
    colons = strfind(tline,';');
    ReverseFields(j) = str2double(tline(colons(end-1)+1:colons(end)-1));
    
    dscan = textscan(fid,'%f',7);
    % Second half major loop
    while length(dscan{1}) == 7
        for i = 1:length(savecols)                        
            DataStruct.(colnames{i}){j}(n) = dscan{1}(savecols(i));
        end
        n = n+1;
        dscan = textscan(fid,'%f',7);
    end
    %Read the rest of the text line
    fgetl(fid);
    j = j + 1;
end

fclose(fid);

%% Assign in imported data with generated varname
DataStruct.Info = Info;
%The last "reverse field" read did not correspond to a minor loop (blame the pk recipe)
DataStruct.ReverseFields = ReverseFields(1:end-1);
assignin('base',varname,DataStruct);
