function varname = mamrimport(filepath)

if strcmp(class(filepath),'double')
    directory = dir;
    filepath = [pwd '\' directory(filepath+2).name];
end

[~, filename, extension] = fileparts(filepath);
varname = filename;                             %Extract filename from path
varname = strrep(varname,'-','_');              %Change - to _
varname = genvarname(varname);                  %Ensure variable name is valid
    
DELIMITER = '\t';
HEADERLINES = 5000;                                          % Should be larger than total lines in file
    
ImportedData = importdata(filepath,DELIMITER,HEADERLINES);          % Import everything into dummy

numlines = length(ImportedData);
 
 numsections = 0;
 P = 0;
 F = 0;
 for i = 1:numlines                                  % Find sections
    if strfind(ImportedData{i},'Power [dBm]')
        P = P+1;
        Power{P} = str2double(ImportedData{i+1});
    end
    if strfind(ImportedData{i},'Frequency [GHz]')
        F = F+1;
        Freq{F} = str2double(ImportedData{i+1});        
    end
    if strfind(ImportedData{i},'V_Hall')
        numsections = numsections + 1;
        sectionstart(numsections) = i;
    end
 end
 
DataStruct.Filepath = filepath;
fileinfo = dir(filepath);
DataStruct.FileDate = fileinfo.date;

for n = 1:numsections
    
    ImportedData = importdata(filepath,DELIMITER,sectionstart(n));                 % Import everything after @@Data +1 (this function will stop import at the next text line)

    savecols = [1 2];
    colnames = {'Field','V_Hall'};               
    

    for i = 1:length(savecols)
        DataStruct.(colnames{i}) = ImportedData.data(:,savecols(i));
    end
    
    if P>=n && F>=n
        DataStruct.Power = Power{n};
        DataStruct.Freq = Freq{n};
    end
    
    if n==1
        assignin('base',varname,DataStruct);                       %Import file as a struct into workspace
    else
        assignin('base',[varname '_' num2str(n)],DataStruct);      % Append section number
    end
    
end