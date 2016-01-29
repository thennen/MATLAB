function PKImportRemDecay_v2(filepath,varargin)

%{
************************************************************************************************************
*   Version 2, started: 12/20/2012, last update: 03/05/2013                                                *                                 *
*                                                                                                          *
*   This function imports data from a PK Remanent Decay file.                                              *
*   Assumes 'filepath' is the path to a data file in the correct format:                                   *
*     1) One major loop (split into two sections, from positive to negative field and viceversa)           *
*     2) Several recoil partial loops (from negative field to 0); these come after setting the field to    *
*        a high positive value (fixed), then to the desired negative value (variable); length may vary     *
*     3) Remanent decay data: acquisition as a function of time, at H = 0, right after each minor loop     *
*     4) One final major loop (splitted in two sections)                                                   *
*                                                                                                          *
*   Function now allows for additional arguments: if second argument is 0, data are only stored,           *
*   not plotted                                                                                            *
*                                                                                                          *
************************************************************************************************************
%}


% CONSTANTS
DELIMITER = '\t';
HEADERLINES = 50000; % This value should be greater than number of lines in the file

%{
*** LOAD DATA - PART 1 ***

- Import whole file as text
- Extract basic text information
- Customize columns to be stored (THIS CAN BE MODIFIED ACCORDING TO ONE'S NEEDS)
- Count number of sections and evaluate max section length
%}

%Initially, import entire file as if it was only header (i.e., line by line, as text) into Imported structure
Imported = importdata(filepath,DELIMITER,HEADERLINES);
numlines = length(Imported);


% Create "Info" struct to hold values entered in PK program
Info.MeasurementDate = Imported{3}(19:end);                  
Info.MacroFile = Imported{4}(18:end);                        
Info.DepToolID = Imported{8};                                
Info.LotID = Imported{10};                                   
Info.DiskID = Imported{11};                                  
Info.OperationID = Imported{12};

% Specify the columns to be saved
savecols = [1 3 4 7];                                 % Column number              
colnames = {'Field','Kerr_x','Kerr_y','Time'};        % Column name

% Count the number of sections and store the maximum length of a section
MAXSECTIONS = 0;
MAXLINES = 0;                                         % max number of lines within a section

% Initially read files once to determine number of sections and max number of lines within a section
for i = 1:numlines
    
    if any(strfind(Imported{i},'X Field Sweep'))      % look for instances of measurement command
        
        data = ReadPKCommand(Imported{i});            % acquire data from command line
        if MAXLINES < floor(data{2}/data{6}-1E-10)    % data{2} = total measure time, data{6} = sampling time interval
                                                      % --> floor(data{2}/data{6}) = number of measurements (or lines)
            MAXLINES = floor(data{2}/data{6}-1E-10);
        end
        MAXSECTIONS = MAXSECTIONS + 1;   
    
    end
end

% Allocate memory for the cell array used to store the data once converted from strings 
CellArrayofVectors = cell(MAXLINES,MAXSECTIONS);
CellArrayofVectors(:)={NaN};


%{ 
*** LOAD DATA - PART 2 ***

Load data from all-text "Imported" structure, convert them into numbers, assign them to "DataStruct" structure.
DataStruct has one field (a matrix) for each column of the original data one wants to copy.
Each measurement is saved in a different column. Unused cells are filled with 'NaN'.
%}

% current section number (in the end, it represents the total number of sections)
sectionnum = 1;

for i = 1:numlines

    if any(strfind(Imported{i},'X Field Sweep'))      % when measurement command is found...
        data = ReadPKCommand(Imported{i});            % interpret command
        data_num = floor(data{2}/data{6}-1E-10);      % save the number of measurements the command triggered
    
        for j = 1 : data_num                   

                % import data from line into a cell array (textscan returns cell arrays) and convert the cell array into a vector
                % then save each vector of numbers into a cell array, one column per section
                % additional text commands, stored at the end of the last line of data, are ignored
                CellArrayofVectors{j,sectionnum} = cell2mat(textscan(Imported{i - data_num + j},'%f %f %f %f %f %f %f %*[\n]'));               
                
                % load data into structure; each structure field is a matrix
                for k = 1:length(savecols)
                    DataStruct.(colnames{k})(j,sectionnum) = CellArrayofVectors{j,sectionnum}(savecols(k));   
                end
                
        end
        
        % fills remaining rows with NaN (otherwise they would be filled with zeros by default)
        for j = data_num+1 : MAXLINES    
            for k = 1:length(savecols)
                DataStruct.(colnames{k})(j,sectionnum) = NaN;
            end
        end
        
        sectionnum = sectionnum + 1;
   end        
end

%{
*** PART 2 ***

Organize DataStruct columns into three different structures containing data from
(1) Major loops
(2) Remanent Decay Measurements
(3) Recoil partial loops

Each structure contains various fields (CUSTOMIZABLE), each of which is a cell array of vectors.
Recoil partial loops vectors length varies.
%}

% Define various counters
l = 1;  % number of Major Loops
m = 1;  % number of Remanent decay measurements
p = 1;  % number of Recoil partial loops

for n = 1:sectionnum-1
    j = 1;

    % Select major loops data (first two and last two sections)
    if n == 1 || n == 2 || n == sectionnum-2 || n == sectionnum-1
        
        % Scan the whole column. Stops at the end or when 'NaN' is found
        while j <= MAXLINES && ~isnan(DataStruct.(colnames{1})(j,n))
    
            % Save each field from column n of DataStruct to a temporary structure having the same fields but one column only            
            for k = 1:length(savecols)
                Tmp.(colnames{k})(j,1) = DataStruct.(colnames{k})(j,n);
            end            
            j = j+1;
            
        end        
        
        % Save each field of the temporary structure to the final recoil loops structure, whose fields are cell array
        for k = 1:length(savecols)
            MajorLoops.(colnames{k}){l} = Tmp.(colnames{k})(:,1);
        end
        
        clear Tmp;        
        l = l + 1;
    
    % odd sections identify remanent decay measurements    
    elseif mod(n,2)
        
        % Scan the whole column. Stops at the end or when 'NaN' is found
        while j <= MAXLINES && ~isnan(DataStruct.(colnames{1})(j,n))
    
            % Save each field from column n of DataStruct to a temporary structure having the same fields but one column only            
            for k = 1:length(savecols)
                Tmp.(colnames{k})(j,1) = DataStruct.(colnames{k})(j,n);
            end            
            j = j+1;
            
        end        
        
        % Save each field of the temporary structure to the final recoil loops structure, whose fields are cell array
        for k = 1:length(savecols)
            RemDecay.(colnames{k}){m} = Tmp.(colnames{k})(:,1);
        end
        
        clear Tmp;        
        m = m + 1;
    
    % even sections identify recoil partial loops            
    elseif ~mod(n,2)
        
        % Scan the whole column. Stops at the end or when 'NaN' is found
        while j <= MAXLINES && ~isnan(DataStruct.(colnames{1})(j,n))
    
            % Save each field from column n of DataStruct to a temporary structure having the same fields but one column only            
            for k = 1:length(savecols)
                Tmp.(colnames{k})(j,1) = DataStruct.(colnames{k})(j,n);
            end            
            j = j+1;
            
        end        
        
        % Save each field of the temporary structure to the final recoil loops structure, whose fields are cell array
        for k = 1:length(savecols)
            RecoilLoops.(colnames{k}){p} = Tmp.(colnames{k})(:,1);
        end
        
        clear Tmp;
        p = p + 1;
    
    end
    
end

% Create name of base workspace variable (a structure containing all the relevant data)
varname = [Info.DiskID '_PKDecay'];
varname = strrep(varname,'-','_');            % Change - to _
varname = genvarname(varname);                % Ensure variable name is valid

DataStruct.Info = Info;
DataStruct.MajorLoops = MajorLoops;
DataStruct.RemDecay = RemDecay;
DataStruct.RecoilLoops = RecoilLoops;

% Check for additional arguments when function was called
% A second argument equal to 0 means data will solely be imported, not processed or plotted

if (nargin > 1)
    
    for i=1:length(varargin)
    
        if ischar(class(varargin{i})) && strcmp(varargin{i},'plot')
                
                if i < length(varargin)
                
                    if (varargin{i+1} ~= 0) % note varargin{1} does not exist if function is called with 1 argument only (the file path)
            
                        % Run analysis routine
                        DataRemDec = RemDecayAnalysis_v2(MajorLoops,RemDecay,RecoilLoops,Info);
        
                        DataStruct.Ms_09 = DataRemDec.Ms_09;
                        DataStruct.RemDecay_Norm = DataRemDec.RemDecay_Norm;
        
                        % Plot recoil loops (inside first major loop)
                        PlotMinorLoops(MajorLoops,RecoilLoops,Info);
                
                    end
                    
                else
                        % Run analysis routine
                        DataRemDec = RemDecayAnalysis_v2(MajorLoops,RemDecay,RecoilLoops,Info);
        
                        DataStruct.Ms_09 = DataRemDec.Ms_09;
                        DataStruct.RemDecay_Norm = DataRemDec.RemDecay_Norm;
        
                        % Plot recoil loops (inside first major loop)
                        PlotMinorLoops(MajorLoops,RecoilLoops,Info);                        
                    
                end
        end   
        
    end
    
else
    
    % Run analysis routine
    % DataRemDec = RemDecayAnalysis_v2(MajorLoops,RemDecay,RecoilLoops,Info);                 
    % DataStruct.Ms_09 = DataRemDec.Ms_09;
    % DataStruct.RemDecay_Norm = DataRemDec.RemDecay_Norm;
    
    % Plot recoil loops (inside first major loop)
    % PlotMinorLoops(MajorLoops,RecoilLoops,Info);

end

% Store variable in base workspace
assignin('base',varname,DataStruct);

end