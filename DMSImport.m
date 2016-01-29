function varname = DMSImport(filepath)
%DMSIMPORT import DMS data file as a Matlab structure

[~, filename, extension] = fileparts(filepath);
%Extract filename from path
varname = filename;
%Change - to _
varname = strrep(varname,'-','_');
%Ensure variable name is valid
varname = genvarname(varname);
    
DELIMITER = ' ';
% Should be larger than total lines in file
HEADERLINES = 10000;
    
% Import everything into dummy
Dummy = importdata(filepath,DELIMITER,HEADERLINES);

numlines = length(Dummy);

FlagsToFind = {'@Date','@Time','@Time at end','@Test ID','@Number of averages','@Emu/v','@Y Coils Correction Factor','@Sample Shape Correction Factor','@Coil Angle Alpha','@Coil Angle Beta'};
FlagVariableNames = {'Date','Time','TimeAtEnd','TestID','NumAvgs','EmuPerVolt','YCC','SSC','Alpha','Beta'};

for i = 1:length(FlagsToFind)
    for j = 1:114
        if strfind(Dummy{j},FlagsToFind{i})
            colonposition = strfind(Dummy{j},':');
            if(i>4)
                DataStruct.MeasurementParameters.(FlagVariableNames{i}) = str2num(Dummy{j}(colonposition+2:end));
            else
                DataStruct.MeasurementParameters.(FlagVariableNames{i}) = Dummy{j}(colonposition+2:end);
            end
            break
        end
    end
end


 for startdata = 1:200
     % Find the first line that says "@@Data"
     if strcmp(Dummy{startdata},'@@Data')
        break
     end
 end
 for enddata = startdata:numlines
     if strcmp(Dummy{enddata},'@@END Data.')
        break
     end
 end
 
 numsections = 0;
 % Find sections
 for i = startdata:enddata
    if strncmp(Dummy{i},'New Section',10)
    numsections = numsections + 1;
    sectionstart(numsections) = i;
    end
 end

for n = 1:numsections
    
    % Import everything after @@Data +1 (this function will stop import at the
    % next text line)
    Dummy = importdata(filepath,DELIMITER,sectionstart(n));
    % determine points measured per measurement variable step, which are each
    % recorded in 18 columns
    % i.e. DCD curve may measure in field and zero field hysteresis values
    NMeasurements = size(Dummy.data,2) / 18;
    
    savecols = [1 3 7 6 12 13];
    colnames = {'Timestamp','Temp','Angle','Field','X','Y'};
    % if there are more data points, save them too
    for i = 2:NMeasurements
        savecols = [savecols 12+(i-1)*18 13+(i-1)*18];
        colnames = [colnames ['X' num2str(i)] ['Y' num2str(i)]];
    end
    

    for i = 1:length(savecols)
        DataStruct.(colnames{i}) = Dummy.data(:,savecols(i));
        if strncmp(colnames{i},'X',1) || strncmp(colnames{i},'Y',1)
            % Convert to microemu
            DataStruct.(colnames{i}) = DataStruct.(colnames{i}).* 1000000;
        end
    end
    
    if n==1
        %Import file as a struct into workspace
        assignin('base',varname,DataStruct);
    else
        % Append section number
        assignin('base',[varname '_' num2str(n)],DataStruct);
    end
    
end
